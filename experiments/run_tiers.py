"""run_tiers.py -- where do the errors of imputation live?
Splits the frequent patterns reported by an imputation method (MICE or kNN) into those that
are also certified (tier 1) and those that are not (tier 2), and measures the precision of
each tier. Tier 1 is error-free by Proposition 7.1; the question is how large it is and how
unreliable tier 2 is (Section 7.3).
Run:  python3 run_tiers.py results/public_tiers.csv      (uses the same masks as run_public.py)"""
import sys, os, time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exp_core import closed_intents, certified_intents, derive_all, make_mask, impute
from run_public import PREP, RATES, MECHS, SEEDS


def run_one(args):
    name, mech, p, seed, scen = args
    if p == 0:
        return []
    d = PREP[name]
    V, c, delta, (T_int, T_sel, _), s0 = d['V'], d['c'], d['delta'], d['truth'], d['s0']
    rng = np.random.default_rng(1000 * seed + 17)          # same masks and noise as run_public.py
    M = make_mask(V, mech, p, rng)
    obs = ~M
    Vm = V + rng.uniform(-1, 1, size=V.shape) * 0.999 * delta if scen == 'B' else V
    Vobs = np.where(M, np.nan, Vm)
    if scen == 'A':
        S, Neg = obs & (Vm > c), obs & (Vm < c)
    else:
        S, Neg = obs & (Vm - c >= delta), obs & (Vm - c <= -delta)
    cert, _ = certified_intents(S, Neg)
    cert = set(np.flatnonzero(cert).tolist())
    rows = []
    for meth in ('mice', 'knn'):
        Rh = impute(Vobs, meth, seed) > c
        ints, sel = closed_intents(Rh)
        supp = sel.sum(1)
        rep = {i for i in ints if i != 0 and supp[i] >= s0}
        t1, t2 = rep & cert, rep - cert
        prec = lambda X: (len(X & T_int) / len(X)) if X else np.nan
        rows.append(dict(dataset=name, mechanism=mech, rate=p, seed=seed, scenario=scen, method=meth,
                         n_rep=len(rep), share_tier1=len(t1) / len(rep) if rep else np.nan,
                         prec_all=prec(rep), prec_tier1=prec(t1), prec_tier2=prec(t2),
                         errors_in_tier2=(len(t2 - T_int) / max(len(rep - T_int), 1)) if rep else np.nan))
    return rows


if __name__ == '__main__':
    jobs = [(n, m, p, s, sc) for sc in ('A', 'B') for n in PREP for m in MECHS for p in RATES for s in SEEDS]
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=2) as ex:
        for res in ex.map(run_one, jobs, chunksize=4):
            rows.extend(res)
    df = pd.DataFrame(rows)
    df.to_csv(sys.argv[1] if len(sys.argv) > 1 else 'results/public_tiers.csv', index=False)
    print(f'done: {len(df)} rows in {time.time() - t0:.0f} s')
