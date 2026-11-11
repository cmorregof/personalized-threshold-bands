"""tsweep_common.py -- the cohort, thresholds and bands shared by remeasure_stats.py and run_tsweep.py.

COHORT is cohort_labs12_v2.parquet from build_mimic_cohort.py (step 3). Attribute m of admission g is
"lab m above the patient's own baseline", with threshold C = b + resolution/4 as in run_mimic.py.
The band of a cell is [C - d, C + d], so its certainty is s = (V - C) / d, where
  * band = 'personal': d is the patient's own dispersion d_<test> (1.4826 MAD of the earlier
    values), floored at a quarter of the reporting resolution;
  * band = 'per_test': d is one half-width per test, the cohort median of the floored d_<test>.
The band is centred on the threshold C rather than on the baseline b, so that S_t = {s >= t} and
Neg_t = {s <= -t} are always on the right side of C; the two differ by a quarter of the resolution.
"""
import numpy as np
import duckdb

LABS = ['sodium', 'potassium', 'chloride', 'bicarbonate', 'bun', 'creatinine',
        'glucose', 'calcium', 'hemoglobin', 'hematocrit', 'wbc', 'platelets']
RES = dict(sodium=1, potassium=0.1, chloride=1, bicarbonate=1, bun=1, creatinine=0.1,
           glucose=1, calcium=0.1, hemoglobin=0.1, hematocrit=0.1, wbc=0.1, platelets=1)
TS = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
BANDS = ['personal', 'per_test']


def load(cohort):
    df = duckdb.sql(f"SELECT * FROM '{cohort}'").df()
    col = lambda p: df[[p + l for l in LABS]].to_numpy(float)
    V, B, D, R = col(''), col('b_'), col('d_'), col('r_')
    floor = np.array([RES[l] for l in LABS]) / 4.0
    C = B + floor
    assert not np.any(np.isclose(V, C, rtol=0, atol=1e-9))
    Dp = np.maximum(D, floor)
    widths = {'personal': Dp, 'per_test': np.broadcast_to(np.median(Dp, axis=0), Dp.shape)}
    return dict(V=V, B=B, C=C, R=R, widths=widths, floored=D < floor)
