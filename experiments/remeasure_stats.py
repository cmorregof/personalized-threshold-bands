"""remeasure_stats.py -- how far a second measurement within 24 hours moves the certainty table.

For every test and band variant, over the admissions of the cohort that have a re-measurement:
the fraction of admissions with one, the median gap in hours, the share of floored dispersions,
the quantiles of |s - s_re| (s from the panel value, s_re from the re-measurement, both in units
of the band half-width), and, for every tolerance t, the share of cells with |s - s_re| >= t and
the share of decided cells (|s| >= t) that the re-measurement puts on the other side of C, which
is the cell-level error rate of certification against re-measurement. Only aggregates are written.

Usage: python3 remeasure_stats.py COHORT_V2.parquet results/mimic_remeasure_stats.csv
"""
import sys, os
import numpy as np
import pandas as pd
import duckdb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tsweep_common import LABS, TS, BANDS, load

if __name__ == '__main__':
    cohort, out = sys.argv[1], sys.argv[2]
    X = load(cohort)
    H = duckdb.sql(f"SELECT * FROM '{cohort}'").df()[['h_' + l for l in LABS]].to_numpy(float)
    rows = []
    for band in BANDS:
        W = X['widths'][band]
        s, sr = (X['V'] - X['C']) / W, (X['R'] - X['C']) / W
        for j, lab in enumerate(LABS + ['all']):
            cols = slice(None) if lab == 'all' else slice(j, j + 1)
            has = ~np.isnan(X['R'][:, cols])
            a, b = s[:, cols][has], sr[:, cols][has]
            dif = np.abs(a - b)
            row = dict(band=band, test=lab, n_cells=int(has.size), frac_remeasured=float(has.mean()),
                       gap_h_median=float(np.nanmedian(H[:, cols])),
                       frac_floored=float(X['floored'][:, cols].mean()),
                       **{f'absdiff_q{q}': float(np.quantile(dif, q / 100)) for q in (50, 90, 99, 99.9)})
            for t in TS:
                dec = np.abs(a) >= t
                row[f'frac_absdiff_ge_{t}'] = float((dif >= t).mean())
                row[f'frac_decided_{t}'] = float(dec.mean())
                row[f'flip_given_decided_{t}'] = float((np.sign(a) != np.sign(b))[dec].mean()) if dec.any() else np.nan
            rows.append(row)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f'{len(rows)} rows -> {out}')
