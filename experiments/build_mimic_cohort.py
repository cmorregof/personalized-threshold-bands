"""build_mimic_cohort.py -- builds the MIMIC-IV cohort of Section 7 of the paper.

MIMIC-IV is distributed by PhysioNet to credentialed users under a data use agreement, so the
data are not included in this repository; this script must be run locally by such a user.

Usage
    python3 build_mimic_cohort.py MIMIC_IV_DIR OUT_DIR                 # both steps
    python3 build_mimic_cohort.py MIMIC_IV_DIR OUT_DIR --labs FILE     # reuse an extraction

Step 1  extracts the twelve laboratory tests below from hosp/labevents.csv.gz into
        OUT_DIR/labs12.parquet (subject_id, itemid, charttime, valuenum; numeric values only).
Step 2  builds OUT_DIR/cohort_labs12.parquet, one row per patient:
        * the admission panel: for every admission and test, the first time at which the test
          was charted between 6 hours before and 24 hours after the admission time, and the
          value charted then (the median, if several values share that time);
        * the personal baseline of every test: the median of the patient's values charted from
          3 years before the admission up to, but excluding, 7 days before it, together with the
          number of such values;
        * eligible admissions: full panel and at least two earlier values of every test;
        * for every patient, the latest eligible admission.
Columns: hadm_id, subject_id, admittime, the twelve panel values, and b_<test> (baseline) and
n_<test> (number of earlier values) for each test.
"""
import os
import sys
import time
import tempfile
import duckdb

TESTS = {  # MIMIC-IV hosp/d_labitems itemid -> column name (all with fluid = Blood)
    50983: 'sodium', 50971: 'potassium', 50902: 'chloride', 50882: 'bicarbonate',
    51006: 'bun', 50912: 'creatinine', 50931: 'glucose', 50893: 'calcium',
    51222: 'hemoglobin', 51221: 'hematocrit', 51301: 'wbc', 51265: 'platelets',
}
LABEVENTS_COLUMNS = (
    "{'labevent_id': 'BIGINT', 'subject_id': 'BIGINT', 'hadm_id': 'BIGINT', 'specimen_id': 'BIGINT', "
    "'itemid': 'INTEGER', 'order_provider_id': 'VARCHAR', 'charttime': 'TIMESTAMP', "
    "'storetime': 'TIMESTAMP', 'value': 'VARCHAR', 'valuenum': 'DOUBLE', 'valueuom': 'VARCHAR', "
    "'ref_range_lower': 'DOUBLE', 'ref_range_upper': 'DOUBLE', 'flag': 'VARCHAR', "
    "'priority': 'VARCHAR', 'comments': 'VARCHAR'}")


def connect():
    """A disk-backed scratch database, so that the joins also run on machines with little memory."""
    tmp = tempfile.mkdtemp(prefix='duckdb_')
    con = duckdb.connect(os.path.join(tmp, 'scratch.duckdb'))
    con.execute(f"SET temp_directory = '{tmp}'")
    con.execute("SET preserve_insertion_order = false")
    return con


def step1(mimic, out):
    src = os.path.join(mimic, 'hosp', 'labevents.csv.gz')
    dst = os.path.join(out, 'labs12.parquet')
    ids = ','.join(str(i) for i in TESTS)
    connect().execute(f"""
        COPY (SELECT subject_id, itemid, charttime, valuenum
              FROM read_csv('{src}', header = true, columns = {LABEVENTS_COLUMNS})
              WHERE itemid IN ({ids}) AND valuenum IS NOT NULL)
        TO '{dst}' (FORMAT parquet)""")
    return dst


def step2(mimic, out, labs):
    adm = os.path.join(mimic, 'hosp', 'admissions.csv.gz')
    dst = os.path.join(out, 'cohort_labs12.parquet')
    con = connect()
    con.execute(f"CREATE TABLE adm AS SELECT hadm_id, subject_id, admittime FROM read_csv_auto('{adm}')")
    con.execute(f"CREATE TABLE labs AS SELECT * FROM '{labs}' WHERE valuenum IS NOT NULL")
    # admission panel: first charting time in [admittime - 6 h, admittime + 24 h], value at that time
    for itemid, name in TESTS.items():
        con.execute(f"""
            CREATE TABLE p_{name} AS
            WITH first_time AS (
                SELECT a.hadm_id, a.subject_id, l.charttime AS t1
                FROM (SELECT hadm_id, subject_id, admittime,
                             admittime - INTERVAL 6 HOUR AS t0 FROM adm) a
                ASOF JOIN (SELECT DISTINCT subject_id, charttime FROM labs WHERE itemid = {itemid}) l
                  ON a.subject_id = l.subject_id AND a.t0 <= l.charttime
                WHERE l.charttime <= a.admittime + INTERVAL 24 HOUR)
            SELECT f.hadm_id, median(l.valuenum) AS {name}
            FROM first_time f JOIN labs l
              ON l.itemid = {itemid} AND l.subject_id = f.subject_id AND l.charttime = f.t1
            GROUP BY f.hadm_id""")
    joins = ' '.join(f'JOIN p_{n} USING (hadm_id)' for n in TESTS.values())
    cols = ', '.join(f'p_{n}.{n}' for n in TESTS.values())
    con.execute(f"CREATE TABLE panel AS SELECT adm.*, {cols} FROM adm {joins}")
    # personal baselines: median of the values from 3 years to 7 days before admission
    con.execute("""
        CREATE TABLE base AS
        SELECT p.hadm_id, l.itemid, median(l.valuenum) AS b, count(*) AS n
        FROM panel p JOIN labs l
          ON p.subject_id = l.subject_id
         AND l.charttime >= p.admittime - INTERVAL 3 YEAR
         AND l.charttime <  p.admittime - INTERVAL 7 DAY
        GROUP BY p.hadm_id, l.itemid""")
    bcols = ', '.join(f"max(CASE WHEN itemid = {i} THEN b END) AS b_{n}, "
                      f"max(CASE WHEN itemid = {i} THEN n END) AS n_{n}" for i, n in TESTS.items())
    full = ' AND '.join(f'n_{n} >= 2' for n in TESTS.values())
    con.execute(f"""
        CREATE TABLE eligible AS
        SELECT * FROM panel JOIN (SELECT hadm_id, {bcols} FROM base GROUP BY hadm_id) USING (hadm_id)
        WHERE {full}""")
    con.execute(f"""
        COPY (SELECT * EXCLUDE (rk) FROM (
                SELECT *, row_number() OVER (PARTITION BY subject_id
                                             ORDER BY admittime DESC, hadm_id DESC) AS rk
                FROM eligible) WHERE rk = 1 ORDER BY hadm_id)
        TO '{dst}' (FORMAT parquet)""")
    n = con.execute(f"SELECT count(*) FROM '{dst}'").fetchone()[0]
    return dst, n


if __name__ == '__main__':
    mimic, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    if '--labs' in sys.argv:
        labs = sys.argv[sys.argv.index('--labs') + 1]
    else:
        labs = step1(mimic, out)
        print(f'step 1 -> {labs} ({time.time() - t0:.0f} s)')
    dst, n = step2(mimic, out, labs)
    print(f'step 2 -> {dst}: {n} admissions ({time.time() - t0:.0f} s)')
