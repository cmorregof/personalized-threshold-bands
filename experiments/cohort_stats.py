"""cohort_stats.py -- aggregate description of the MIMIC-IV cohort quoted in Section 7.

Usage: python3 cohort_stats.py COHORT.parquet OUT.csv
Writes only aggregate numbers: the size of the cohort, the number of earlier values behind the
baselines, how often reading a test against the patient's own baseline instead of the cohort
median changes the reading, and how many cells lie within one half-width of their baseline.
"""
import sys
import numpy as np
import pandas as pd
import duckdb

TESTS = ['sodium', 'potassium', 'chloride', 'bicarbonate', 'bun', 'creatinine',
         'glucose', 'calcium', 'hemoglobin', 'hematocrit', 'wbc', 'platelets']
RES = dict(sodium=1, potassium=0.1, chloride=1, bicarbonate=1, bun=1, creatinine=0.1,
           glucose=1, calcium=0.1, hemoglobin=0.1, hematocrit=0.1, wbc=0.1, platelets=1)

if __name__ == '__main__':
    df = duckdb.sql(f"SELECT * FROM '{sys.argv[1]}'").df()
    V = df[TESTS].to_numpy(float)
    B = df[['b_' + t for t in TESTS]].to_numpy(float)
    n_prior = df[['n_' + t for t in TESTS]].to_numpy(float)
    r = np.array([RES[t] for t in TESTS])
    personal = B + r / 4                       # the thresholds used in run_mimic.py
    cohort = np.median(V, axis=0) + r / 4      # one threshold per test: the cohort median
    assert not np.isclose(V, personal, rtol=0, atol=1e-9).any()
    assert not np.isclose(V, cohort, rtol=0, atol=1e-9).any()
    flip = (V > personal) != (V > cohort)
    q75, q25 = np.percentile(V, [75, 25], axis=0)
    delta = 0.1 * (q75 - q25)
    rows = [('admissions', V.shape[0]),
            ('median_earlier_values_per_test', float(np.median(n_prior))),
            ('share_cells_read_differently', flip.mean()),
            ('min_share_by_test', flip.mean(0).min()),
            ('max_share_by_test', flip.mean(0).max()),
            ('share_admissions_with_at_least_one', flip.any(1).mean()),
            ('share_admissions_with_at_least_three', (flip.sum(1) >= 3).mean()),
            ('share_cells_above_own_baseline', (V > personal).mean()),
            ('share_cells_within_one_half_width', (np.abs(V - personal) < delta).mean())]
    rows += [(f'share_read_differently_{t}', flip[:, j].mean()) for j, t in enumerate(TESTS)]
    rows += [(f'half_width_{t}', delta[j]) for j, t in enumerate(TESTS)]
    out = pd.DataFrame(rows, columns=['quantity', 'value'])
    out.to_csv(sys.argv[2], index=False)
    print(out.to_string(index=False))
