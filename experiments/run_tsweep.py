"""run_tsweep.py -- experiment E1: personal bands and a sweep of the tolerance t on MIMIC-IV.

Same objects, attributes and imputation baselines as run_mimic.py, but every observed cell now has
a band of positive width (see tsweep_common.py), so the certainty table s is finite and the
tolerance t decides which cells are sure: S_t = obs & (s >= t), Neg_t = obs & (s <= -t); a missing
cell has s = 0 and is unknown at every t.

Two truths:
  * truth = 'obs': R* = V > C, the complete panel (as in run_mimic.py). The complete data are a
    completion at every t, so certification never errs (asserted); t only costs abstention.
  * truth = 'remeasure': R* = r > C, where r is the next value of the same test within 24 h of the
    panel value. The admissions are sampled from those with a re-measurement of all twelve tests,
    so R* is a complete context and no cell needs special treatment. Certification can now err,
    and its error is reported.
Both truths are evaluated on the same sampled admissions, which therefore come from the subset with
complete re-measurement, not from the whole cohort as in run_mimic.py.

Writes only aggregate metrics:
  OUT.csv          patterns and concepts (as run_mimic.py), one row per run, method, band, t, truth;
                   imputation rows do not depend on band or t and carry them empty;
  OUT_queries.csv  patient-pattern questions (true patterns of at most three labs, support >= 5%);
  OUT with 'sweep' replaced by 'rho_int': certified intents at t = 0.25 grouped by deciles of rho_int (equal-count groups
                   after sorting by rho_int), with the share that are intents of R*.

Usage: python3 run_tsweep.py COHORT_V2.parquet results/mimic_T_sweep.csv [--quick]
       (--quick: N = 200 and 2 repetitions; set OMP_NUM_THREADS=1)
"""
import sys, os, time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exp_core import (make_mask, impute, concepts_of, certified_intents,
                      evaluate_point, evaluate_certified, derive_all)
from run_mimic import mask_near_baseline, impute_aux
from tsweep_common import TS, BANDS, load

QUICK = '--quick' in sys.argv
N = 200 if QUICK else 1000
SEEDS = list(range(2 if QUICK else 10))
CONFIGS = [('MCAR', 0.0), ('MCAR', 0.10), ('MNARc', 0.10)]
RHO_T = 0.25

X = load(sys.argv[1])
FULL = np.flatnonzero(~np.isnan(X['R']).any(axis=1))      # admissions re-measured on all twelve tests


def question_errors(Rh, P_eval, Tm, undecided):
    Im = derive_all(Rh)[P_eval]
    wrong = Im != Tm
    return dict(err_all=float(wrong.mean()),
                err_where_cert_abstains=float(wrong[undecided].mean()) if undecided.any() else np.nan,
                err_where_cert_decides=float(wrong[~undecided].mean()) if (~undecided).any() else np.nan)


def run_one(args):
    mech, p, seed = args
    rng = np.random.default_rng(10007 * seed + 3)
    idx = FULL[rng.choice(len(FULL), N, replace=False)]
    V, B, C, Rm = X['V'][idx], X['B'][idx], X['C'][idx], X['R'][idx]
    s0 = int(np.ceil(0.05 * N))
    M = mask_near_baseline(V, C, p, rng) if mech == 'MNARc' else make_mask(V, mech, p, rng)
    obs = ~M
    Vobs = np.where(M, np.nan, V)
    base = dict(dataset='mimic_labs', N=N, mechanism=mech, rate=p, seed=seed, real_rate=float(M.mean()))
    # imputed contexts do not depend on the band or on t
    imputed = {}
    for meth in ('mean', 'median', 'knn', 'mice'):
        if p == 0:
            Vi = V
        elif meth in ('mean', 'median'):
            Vi = impute(Vobs, meth, seed)
        else:
            Vi = impute_aux(Vobs, B, meth, seed)
        imputed[meth] = Vi > C
    imputed['absent'] = obs & (V > C)
    S_band = {}
    for band in BANDS:
        s = np.where(obs, (V - C) / X['widths'][band][idx], 0.0)
        S_band[band] = s
    rows, queries, rhos = [], [], []
    for truth, T in (('obs', V > C), ('remeasure', Rm > C)):
        T_con, T_int, T_sel = concepts_of(T)
        info = (T_int, T_sel, T_con)
        tb = dict(base, truth=truth, n_true_intents=len(T_int))
        for meth, Rh in imputed.items():
            rows.append(dict(tb, method=meth, band='', t=np.nan, **evaluate_point(info, Rh, s0)))
        supp_true = T_sel.sum(1)
        P_eval = [i for i in T_int if i != 0 and supp_true[i] >= s0 and bin(i).count('1') <= 3]
        Tm = T_sel[P_eval]
        for band in BANDS:
            s = S_band[band]
            rho = None
            for t in TS:
                S, Neg = s >= t, s <= -t
                if t == RHO_T:
                    cert, _, rho = certified_intents(S, Neg, s)
                    assert np.array_equal(cert, rho >= t)
                r = evaluate_certified(info, S, Neg, s0)
                if truth == 'obs':
                    assert r['covered'] and (np.isnan(r['prec_int']) or r['prec_int'] == 1.0), r
                    assert np.isnan(r['prec_con']) or r['prec_con'] == 1.0, r
                rows.append(dict(tb, method='certified', band=band, t=t, **r))
                if not P_eval:
                    continue
                Sm, Pm = derive_all(S)[P_eval], derive_all(~Neg)[P_eval]
                undecided = Pm & ~Sm
                cert_wrong = (Sm & ~Tm) | (~Pm & Tm)
                if truth == 'obs':
                    assert not cert_wrong.any()            # certified answers are never wrong
                qb = dict(tb, band=band, t=t, n_patterns=len(P_eval),
                          abstain_cert=float(undecided.mean()),
                          cert_err_all=float(cert_wrong.mean()),
                          cert_err_decided=float(cert_wrong[~undecided].mean()) if (~undecided).any() else np.nan)
                for meth, Rh in imputed.items():
                    queries.append(dict(qb, method=meth, **question_errors(Rh, P_eval, Tm, undecided)))
            # rho_int as a predictor: certified intents at t = RHO_T, frequent in the upper context
            upper = derive_all(~(s <= -RHO_T)).sum(1)
            cand = np.flatnonzero((rho >= RHO_T) & (upper >= s0))
            cand = cand[cand != 0]
            if len(cand):
                order = cand[np.argsort(rho[cand], kind='stable')]
                for k, grp in enumerate(np.array_split(order, 10)):
                    if len(grp) == 0:
                        continue
                    rg = rho[grp]
                    rhos.append(dict(tb, band=band, decile=k + 1, n=len(grp),
                                     rho_min=float(rg.min()), rho_mean=float(np.mean(np.minimum(rg, 1e6))),
                                     rho_max=float(rg.max()),
                                     prec=float(np.mean([g in T_int for g in grp]))))
    return rows, queries, rhos


if __name__ == '__main__':
    out = sys.argv[2]
    jobs = [(m, p, s) for m, p in CONFIGS for s in SEEDS]
    t0 = time.time()
    rows, queries, rhos = [], [], []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for r, q, h in ex.map(run_one, jobs):
            rows.extend(r)
            queries.extend(q)
            rhos.extend(h)
    pd.DataFrame(rows).to_csv(out, index=False)
    pd.DataFrame(queries).to_csv(out.replace('.csv', '_queries.csv'), index=False)
    pd.DataFrame(rhos).to_csv(out.replace('sweep', 'rho_int'), index=False)
    print(f'{len(jobs)} runs, N = {N}, subset with full re-measurement {len(FULL)} '
          f'of {X["V"].shape[0]}, {time.time() - t0:.0f} s')
