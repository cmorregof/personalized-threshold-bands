"""run_mimic.py -- the experiments of Section 7 on MIMIC-IV (run locally by a credentialed user).

Objects: hospital admissions (one per patient) with a complete admission panel of 12 labs.
Attribute m of admission g: "lab m above the patient's own baseline", where the baseline is the
median of the patient's values from 7 days to 3 years before the admission. The true context is
the complete panel binarized against those personal baselines. Only aggregate metrics are written.

Usage: python3 run_mimic.py COHORT.parquet OUT.csv SCENARIO MECH[,MECH...] [N]
The runs of the paper (COHORT from build_mimic_cohort.py; set OMP_NUM_THREADS=1):
    python3 run_mimic.py cohort_labs12.parquet results/mimic_A_MCAR_MAR.csv   A MCAR,MAR 300
    python3 run_mimic.py cohort_labs12.parquet results/mimic_A_MNAR_MNARc.csv A MNAR,MNARc 300
    python3 run_mimic.py cohort_labs12.parquet results/mimic_B_MCAR_MNARc.csv B MCAR,MNARc 300
    python3 run_mimic.py cohort_labs12.parquet results/mimic_A_N1000.csv      A MCAR,MAR,MNAR,MNARc 1000
Each writes OUT.csv (patterns), OUT_tiers.csv and OUT_queries.csv (patient-pattern questions).
"""
import sys, os, time
import numpy as np
import pandas as pd
import duckdb
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exp_core import (make_mask, impute, concepts_of, closed_intents, certified_intents,
                        evaluate_point, evaluate_certified, derive_all)
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import KNNImputer, IterativeImputer

LABS = ['sodium', 'potassium', 'chloride', 'bicarbonate', 'bun', 'creatinine',
        'glucose', 'calcium', 'hemoglobin', 'hematocrit', 'wbc', 'platelets']
RES = dict(sodium=1, potassium=0.1, chloride=1, bicarbonate=1, bun=1, creatinine=0.1,
           glucose=1, calcium=0.1, hemoglobin=0.1, hematocrit=0.1, wbc=0.1, platelets=1)
RATES = [0.0, 0.05, 0.10, 0.20, 0.30]
SEEDS = list(range(10))

COHORT = sys.argv[1] if len(sys.argv) > 1 else 'cohort_labs12.parquet'
_df = duckdb.sql(f"SELECT * FROM '{COHORT}'").df()
V_ALL = _df[LABS].to_numpy(float)
B_ALL = _df[['b_' + l for l in LABS]].to_numpy(float)
# "above baseline" means strictly above the median; the quarter-resolution offset avoids ties
C_ALL = B_ALL + np.array([RES[l] for l in LABS]) / 4.0
assert not np.any(np.isclose(V_ALL, C_ALL, rtol=0, atol=1e-9))
Q75, Q25 = np.percentile(V_ALL, [75, 25], axis=0)
DELTA = 0.1 * (Q75 - Q25)                       # bounded-noise level for scenario B
del _df


def _ranks(X):
    N = X.shape[0]
    return (np.argsort(np.argsort(X, axis=0, kind='stable'), axis=0) + 0.5) / N


def mask_near_baseline(V, C, p, rng):
    """Clinical MNAR: values close to the patient's own baseline are more likely to be missing."""
    if p == 0:
        return np.zeros(V.shape, bool)
    r = _ranks(np.abs(V - C))
    return rng.random(V.shape) < np.clip(2 * p * (1 - r), 0, 1)


def impute_aux(Vobs, Aux, method, seed):
    """kNN or MICE using the personal baselines as fully observed auxiliary predictors."""
    X = np.hstack([Vobs, Aux])
    mu, sd = np.nanmean(X, 0), np.nanstd(X, 0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    if method == 'knn':
        Zi = KNNImputer(n_neighbors=5).fit_transform(Z)
    else:
        Zi = IterativeImputer(max_iter=10, random_state=seed).fit_transform(Z)
    return (Zi * sd + mu)[:, :Vobs.shape[1]]


def run_one(args):
    mech, p, seed, scen, N = args
    rng = np.random.default_rng(10007 * seed + 3)
    idx = rng.choice(V_ALL.shape[0], N, replace=False)
    V, B, C = V_ALL[idx], B_ALL[idx], C_ALL[idx]
    T_con, T_int, T_sel = concepts_of(V > C)
    truth = (T_int, T_sel, T_con)
    s0 = int(np.ceil(0.05 * N))
    M = mask_near_baseline(V, C, p, rng) if mech == 'MNARc' else make_mask(V, mech, p, rng)
    obs = ~M
    Vm = V + rng.uniform(-1, 1, size=V.shape) * 0.999 * DELTA if scen == 'B' else V
    Vobs = np.where(M, np.nan, Vm)
    base = dict(dataset='mimic_labs', N=N, mechanism=mech, rate=p, seed=seed, scenario=scen,
                real_rate=float(M.mean()), n_true_intents=len(T_int))
    rows, tiers = [], []
    if scen == 'A':
        S, Neg = obs & (Vm > C), obs & (Vm < C)
    else:
        S, Neg = obs & (Vm - C >= DELTA), obs & (Vm - C <= -DELTA)
    cert_mask, selS = certified_intents(S, Neg)
    cert = set(np.flatnonzero(cert_mask).tolist())
    # patient-level view: who has each frequent true pattern of up to three labs?
    supp_true = T_sel.sum(1)
    P_eval = [i for i in T_int if i != 0 and supp_true[i] >= s0 and bin(i).count('1') <= 3]
    Tm = T_sel[P_eval]
    Sm = selS[P_eval]
    Pm = derive_all(~Neg)[P_eval]
    undecided = Pm & ~Sm
    assert not (Sm & ~Tm).any() and not (Tm & ~Pm).any()      # certified answers are never wrong
    memb = []
    width_small = float(np.mean((Pm.sum(1) - Sm.sum(1)) / Tm.sum(1))) if P_eval else np.nan
    for meth in ('mean', 'median', 'knn', 'mice'):
        if p == 0:
            Vi = Vm
        elif meth in ('mean', 'median'):
            Vi = impute(Vobs, meth, seed)
        else:
            Vi = impute_aux(Vobs, B, meth, seed)
        Rh = Vi > C
        rows.append(dict(base, method=meth, **evaluate_point(truth, Rh, s0)))
        if P_eval:
            Im = derive_all(Rh)[P_eval]
            wrong = Im != Tm
            memb.append(dict(base, method=meth, n_patterns=len(P_eval),
                             abstain_cert=float(undecided.mean()), err_all=float(wrong.mean()),
                             err_where_cert_abstains=float(wrong[undecided].mean()) if undecided.any() else np.nan,
                             err_where_cert_decides=float(wrong[~undecided].mean()),
                             supp_relerr=float(np.mean(np.abs(Im.sum(1) - Tm.sum(1)) / Tm.sum(1))),
                             cert_width=width_small))
        if meth in ('mice', 'knn') and p > 0:
            ints, sel = closed_intents(Rh)
            supp = sel.sum(1)
            rep = {i for i in ints if i != 0 and supp[i] >= s0}
            t1, t2 = rep & cert, rep - cert
            pr = lambda X: (len(X & T_int) / len(X)) if X else np.nan
            tiers.append(dict(base, method=meth, n_rep=len(rep),
                              share_tier1=len(t1) / len(rep) if rep else np.nan,
                              prec_all=pr(rep), prec_tier1=pr(t1), prec_tier2=pr(t2)))
    rows.append(dict(base, method='absent', **evaluate_point(truth, obs & (Vm > C), s0)))
    r = evaluate_certified(truth, S, Neg, s0)
    assert r['covered'] and (np.isnan(r['prec_int']) or r['prec_int'] == 1.0), r
    assert np.isnan(r['prec_con']) or r['prec_con'] == 1.0, r
    rows.append(dict(base, method='certified', **r))
    if scen == 'B':
        S0, Neg0 = obs & (Vm > C), obs & (Vm < C)
        rows.append(dict(base, method='certified_noband', **evaluate_certified(truth, S0, Neg0, s0)))
    return rows, tiers, memb


if __name__ == '__main__':
    out, scen, mechs = sys.argv[2], sys.argv[3], sys.argv[4].split(',')
    N = int(sys.argv[5]) if len(sys.argv) > 5 else 300
    jobs = [(m, p, s, scen, N) for m in mechs for p in RATES for s in SEEDS]
    t0 = time.time()
    rows, tiers, membs = [], [], []
    with ProcessPoolExecutor(max_workers=4) as ex:
        for r, t, mb in ex.map(run_one, jobs, chunksize=2):
            rows.extend(r)
            tiers.extend(t)
            membs.extend(mb)
    pd.DataFrame(rows).to_csv(out, index=False)
    pd.DataFrame(tiers).to_csv(out.replace('.csv', '_tiers.csv'), index=False)
    pd.DataFrame(membs).to_csv(out.replace('.csv', '_queries.csv'), index=False)
    print(f'{len(jobs)} configs, {len(rows)} rows, cohort size {V_ALL.shape[0]}, {time.time() - t0:.0f} s')
