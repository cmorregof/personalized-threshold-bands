"""run_abstain.py -- experiment E2: an imputation competitor that abstains, on MIMIC-IV.

Patient-pattern questions as in run_mimic.py (true patterns of at most three tests, support >= 5%).
Two abstaining competitors, each built from m = 20 completed tables:
  * 'mi':       multiple imputation, IterativeImputer(sample_posterior=True), the baselines as
                fully observed predictors, random_state = 1000 * seed + k for k = 0..19;
  * 'knn_boot': kNN imputation (k = 5) with donors drawn from a bootstrap resample of the admissions.
For each question the vote is the share of the m tables that answer yes. At level alpha the competitor
answers yes if vote >= 1 - alpha, no if vote <= alpha, and abstains in between; alpha = 0 abstains
unless the m tables agree, and alpha = 0.5 always answers (majority, ties counted as yes).
Each row reports coverage (share of questions answered) and risk (share of wrong answers among them),
also restricted to the questions that certification leaves open (in the sweep design, those left open
by personal bands at t = 1). Single imputations ('mice', and 'absent', the least completion) answer every
question and appear with alpha = 0.5.

Two designs, both writing only aggregates:
  * design = 'main':  the samples, masks and certification (bands of width zero) of run_mimic.py with
                      N = 1000, truth = the complete panel; certification is one point.
  * design = 'sweep': the samples of run_tsweep.py (admissions re-measured on all twelve tests),
                      truths 'obs' and 'remeasure'; certification with personal bands and bands per
                      test is one point per tolerance t.
Mechanisms MCAR and clinical MNAR at 10%, ten repetitions each.

Usage: python3 run_abstain.py COHORT_V2.parquet results/mimic_abstain.csv [--quick]
       (--quick: N = 200, 2 repetitions, m = 5; set OMP_NUM_THREADS=1)
"""
import sys, os, time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exp_core import make_mask, concepts_of, derive_all
from run_mimic import mask_near_baseline, impute_aux
from tsweep_common import TS, BANDS, load
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import KNNImputer, IterativeImputer

QUICK = '--quick' in sys.argv
N = 200 if QUICK else 1000
SEEDS = list(range(2 if QUICK else 10))
M_IMP = 5 if QUICK else 20
ALPHAS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]
MECHS = [('MCAR', 0.10), ('MNARc', 0.10)]

X = load(sys.argv[1])
FULL = np.flatnonzero(~np.isnan(X['R']).any(axis=1))


def completed_tables(Vobs, B, seed):
    """The m completed tables of each competitor, in original units."""
    Xa = np.hstack([Vobs, B])
    mu, sd = np.nanmean(Xa, 0), np.nanstd(Xa, 0)
    sd[sd == 0] = 1.0
    Z = (Xa - mu) / sd
    n = Vobs.shape[1]
    back = lambda Zi: (Zi * sd + mu)[:, :n]
    out = {'mi': [], 'knn_boot': []}
    rng = np.random.default_rng(7919 * seed + 11)
    for k in range(M_IMP):
        imp = IterativeImputer(max_iter=10, sample_posterior=True, random_state=1000 * seed + k)
        out['mi'].append(back(imp.fit_transform(Z)))
        boot = rng.integers(0, len(Z), len(Z))
        out['knn_boot'].append(back(KNNImputer(n_neighbors=5).fit(Z[boot]).transform(Z)))
    return out


def vote_rows(base, answers, Tm, open_q):
    """answers: list of boolean matrices (patterns x objects), one per completed table."""
    vote = np.mean(answers, axis=0)
    rows = []
    for a in ALPHAS:
        yes = vote >= 1 - a if a < 0.5 else vote >= 0.5
        no = vote <= a if a < 0.5 else ~yes
        answered = yes | no
        wrong = (yes & ~Tm) | (no & Tm)
        ao = answered & open_q
        rows.append(dict(base, alpha=a, coverage=float(answered.mean()),
                         risk=float(wrong[answered].mean()) if answered.any() else np.nan,
                         coverage_open=float(ao.sum() / open_q.sum()) if open_q.any() else np.nan,
                         risk_open=float(wrong[ao].mean()) if ao.any() else np.nan))
    return rows


def run_one(args):
    design, mech, p, seed = args
    rng = np.random.default_rng(10007 * seed + 3)
    if design == 'main':
        idx = rng.choice(X['V'].shape[0], N, replace=False)
    else:
        idx = FULL[rng.choice(len(FULL), N, replace=False)]
    V, B, C, Rm = X['V'][idx], X['B'][idx], X['C'][idx], X['R'][idx]
    s0 = int(np.ceil(0.05 * N))
    M = mask_near_baseline(V, C, p, rng) if mech == 'MNARc' else make_mask(V, mech, p, rng)
    obs = ~M
    Vobs = np.where(M, np.nan, V)
    tables = completed_tables(Vobs, B, seed)
    single = {'mice': impute_aux(Vobs, B, 'mice', seed) > C, 'absent': obs & (V > C)}
    truths = [('obs', V > C)] if design == 'main' else [('obs', V > C), ('remeasure', Rm > C)]
    if design == 'main':
        certs = [('zero', np.nan, obs & (V > C), obs & (V < C))]
    else:
        certs = []
        for band in BANDS:
            s = np.where(obs, (V - C) / X['widths'][band][idx], 0.0)
            certs += [(band, t, s >= t, s <= -t) for t in TS]
    rows = []
    for truth, T in truths:
        _, T_int, T_sel = concepts_of(T)
        supp = T_sel.sum(1)
        P_eval = [i for i in T_int if i != 0 and supp[i] >= s0 and bin(i).count('1') <= 3]
        Tm = T_sel[P_eval]
        base = dict(design=design, N=N, mechanism=mech, rate=p, seed=seed, truth=truth,
                    n_patterns=len(P_eval), m=M_IMP)
        # the questions left open by certification with the band of the main comparison
        band0, t0, S0, Neg0 = certs[0] if design == 'main' else next(c for c in certs if c[0] == 'personal' and c[1] == 1.0)
        open_q = derive_all(~Neg0)[P_eval] & ~derive_all(S0)[P_eval]
        for name, tabs in tables.items():
            ans = [derive_all(Vi > C)[P_eval] for Vi in tabs]
            rows += [dict(r, method=name) for r in vote_rows(base, ans, Tm, open_q)]
        for name, Rh in single.items():
            rows += [dict(r, method=name) for r in vote_rows(base, [derive_all(Rh)[P_eval]], Tm, open_q) if r['alpha'] == 0.5]
        for band, t, S, Neg in certs:
            Sm, Pm = derive_all(S)[P_eval], derive_all(~Neg)[P_eval]
            und = Pm & ~Sm
            wrong = (Sm & ~Tm) | (~Pm & Tm)
            if truth == 'obs':
                assert not wrong.any()
            rows.append(dict(base, method='certified', band=band, t=t, coverage=float((~und).mean()),
                             risk=float(wrong[~und].mean()) if (~und).any() else np.nan))
    return rows


if __name__ == '__main__':
    out = sys.argv[2]
    jobs = [(d, m, p, s) for d in ('main', 'sweep') for m, p in MECHS for s in SEEDS]
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for r in ex.map(run_one, jobs):
            rows.extend(r)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f'{len(jobs)} runs, N = {N}, m = {M_IMP}, {time.time() - t0:.0f} s')
