"""run_public.py -- the experiments of Section 7 on the three public datasets.
Scenario A: missing values only.  Scenario B: missing values plus bounded measurement noise
(|noise| < delta = 0.1 IQR per attribute); there certification uses the band t = 1.
Run:  python3 run_public.py results/public_runs.csv      (900 runs; a few minutes)"""
import sys, os, time, itertools
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exp_core import (load_datasets, fixed_thresholds, half_widths, make_mask, impute,
                        concepts_of, evaluate_point, evaluate_certified)

RATES = [0.0, 0.05, 0.10, 0.20, 0.30]
MECHS = ['MCAR', 'MAR', 'MNAR']
SEEDS = list(range(10))
IMPUTERS = ['mean', 'median', 'knn', 'mice']

DATA = load_datasets()
PREP = {}
for name, V in DATA.items():
    c = fixed_thresholds(V)
    R_true = V > c
    T_con, T_int, T_sel = concepts_of(R_true)
    PREP[name] = dict(V=V, c=c, delta=half_widths(V, 0.1), truth=(T_int, T_sel, T_con),
                      s0=int(np.ceil(0.05 * V.shape[0])))


def run_one(args):
    name, mech, p, seed, scen = args
    d = PREP[name]
    V, c, delta, truth, s0 = d['V'], d['c'], d['delta'], d['truth'], d['s0']
    rng = np.random.default_rng(1000 * seed + 17)
    M = make_mask(V, mech, p, rng)                     # missing cells
    obs = ~M
    if scen == 'B':
        noise = rng.uniform(-1, 1, size=V.shape) * 0.999 * delta
        Vm = V + noise
    else:
        Vm = V
    Vobs = np.where(M, np.nan, Vm)
    rows = []
    base = dict(dataset=name, mechanism=mech, rate=p, seed=seed, scenario=scen,
                real_rate=float(M.mean()))
    for meth in IMPUTERS:
        Vi = impute(Vobs, meth, seed) if p > 0 else Vm
        rows.append(dict(base, method=meth, **evaluate_point(truth, Vi > c, s0)))
    rows.append(dict(base, method='absent', **evaluate_point(truth, obs & (Vm > c), s0)))
    if scen == 'A':
        S, Neg = obs & (Vm > c), obs & (Vm < c)
        r = evaluate_certified(truth, S, Neg, s0)
        assert r['covered'] and (np.isnan(r['prec_int']) or r['prec_int'] == 1.0), r
        assert np.isnan(r['prec_con']) or r['prec_con'] == 1.0, r
        rows.append(dict(base, method='certified', **r))
    else:
        S, Neg = obs & (Vm - c >= delta), obs & (Vm - c <= -delta)
        r = evaluate_certified(truth, S, Neg, s0)
        assert r['covered'] and (np.isnan(r['prec_int']) or r['prec_int'] == 1.0), r
        assert np.isnan(r['prec_con']) or r['prec_con'] == 1.0, r
        rows.append(dict(base, method='certified', **r))
        S0, Neg0 = obs & (Vm > c), obs & (Vm < c)
        rows.append(dict(base, method='certified_noband', **evaluate_certified(truth, S0, Neg0, s0)))
    return rows


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'results/public_runs.csv'
    jobs = [(n, m, p, s, sc) for sc in ('A', 'B') for n in DATA for m in MECHS for p in RATES for s in SEEDS]
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=2) as ex:
        for k, res in enumerate(ex.map(run_one, jobs, chunksize=4)):
            rows.extend(res)
            if (k + 1) % 100 == 0:
                print(f'{k + 1}/{len(jobs)} configs, {time.time() - t0:.0f} s', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    for name, d in PREP.items():
        print(name, 'objects', d['V'].shape[0], 'attributes', d['V'].shape[1],
              'true intents', len(d['truth'][0]), 'min support', d['s0'])
    print(f'done: {len(df)} rows in {time.time() - t0:.0f} s')
