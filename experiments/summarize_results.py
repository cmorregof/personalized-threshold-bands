"""summarize_results.py -- every number quoted in Section 7 of the paper, recomputed from the
aggregate result tables in results/ (one row per configuration and method; no patient data).
Run:  python3 summarize_results.py"""
import os, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda f: pd.read_csv(os.path.join(HERE, 'results', f))
IMP = ['mean', 'median', 'knn', 'mice']
f3 = lambda x: f'{x:.3f}'

# ------------------------------------------------------------------ public datasets
P = R('public_runs.csv')
cfg = P[['dataset', 'mechanism', 'rate', 'seed', 'scenario']].drop_duplicates()
print(f'public: {len(cfg)} configurations')
C = P[P.method == 'certified']
assert (C.prec_int.fillna(1) == 1).all() and (C.prec_con.fillna(1) == 1).all() and C.covered.astype(bool).all()
print('  certified: precision 1 on intents and concepts, intervals cover, in every configuration')

def tab(D, rate, metric, methods):
    d = D[np.isclose(D.rate, rate)]
    return {m: d[d.method == m][metric].mean() for m in methods}

A = P[P.scenario == 'A']
for mech in ['MCAR', 'MAR', 'MNAR']:
    d = A[A.mechanism == mech]
    for rate in [0.05, 0.10, 0.20, 0.30]:
        pr = tab(d, rate, 'prec_int', IMP + ['certified'])
        rc = tab(d, rate, 'rec_int', IMP + ['certified'])
        f1 = {m: 2 * pr[m] * rc[m] / (pr[m] + rc[m]) for m in pr}
        print(f'  A {mech} {int(rate*100):2d}%  prec ' + ' '.join(f'{m}={f3(v)}' for m, v in pr.items()))
        print(f'               rec  ' + ' '.join(f'{m}={f3(v)}' for m, v in rc.items()))
        print(f'               F1   ' + ' '.join(f'{m}={f3(v)}' for m, v in f1.items()))
print('  exact concepts, mean recall at 10% MCAR:', {m: f3(v) for m, v in tab(A[A.mechanism == 'MCAR'], 0.10, 'rec_con', IMP + ['certified']).items()})
print('  exact concepts, recall at 0%:', {m: f3(v) for m, v in tab(A[A.mechanism == 'MCAR'], 0.0, 'rec_con', IMP + ['certified']).items()})
B = P[(P.scenario == 'B') & (P.mechanism == 'MCAR')]
for rate in [0.0, 0.05, 0.10, 0.20]:
    pr = tab(B, rate, 'prec_int', IMP + ['certified', 'certified_noband'])
    rc = tab(B, rate, 'rec_int', IMP + ['certified', 'certified_noband'])
    print(f'  B MCAR {int(rate*100):2d}% prec ' + ' '.join(f'{m}={f3(v)}' for m, v in pr.items()))
    print(f'                rec  ' + ' '.join(f'{m}={f3(v)}' for m, v in rc.items()))
T = R('public_tiers.csv')
t = T[(T.scenario == 'A') & (T.mechanism == 'MCAR') & np.isclose(T.rate, 0.10)]
for m in ['mice', 'knn']:
    d = t[t.method == m]
    assert (d.prec_tier1.fillna(1) == 1).all()
    print(f'  tiers A MCAR 10% {m}: certified share {f3(d.share_tier1.mean())}, precision of the rest {f3(d.prec_tier2.mean())}, overall {f3(d.prec_all.mean())} (means of per-run ratios)')
    n1 = d.n_rep * d.share_tier1; n2 = d.n_rep - n1; ok2 = d.prec_tier2.fillna(1) * n2
    print(f'                      pooled over runs (quoted in the paper): certified share {f3(n1.sum() / d.n_rep.sum())}, '
          f'precision of the rest {f3(ok2.sum() / n2.sum())}, overall {f3((n1.sum() + ok2.sum()) / d.n_rep.sum())}')

# ------------------------------------------------------------------ MIMIC-IV
st = R('mimic_cohort_stats.csv').set_index('quantity').value
print(f"\nMIMIC cohort: {int(st['admissions'])} admissions, median {st['median_earlier_values_per_test']:.0f} earlier values per test")
print(f"  read differently against the own baseline than against the cohort median: {st['share_cells_read_differently']:.3f} of the cells "
      f"({st['min_share_by_test']:.3f} to {st['max_share_by_test']:.3f} by test); admissions with at least one: {st['share_admissions_with_at_least_one']:.3f}, "
      f"with at least three: {st['share_admissions_with_at_least_three']:.3f}")
print(f"  cells within one half-width (0.1 IQR) of the baseline: {st['share_cells_within_one_half_width']:.3f}")
pat = pd.concat([R(f) for f in ['mimic_A_MCAR_MAR.csv', 'mimic_A_MNAR_MNARc.csv', 'mimic_A_N1000.csv', 'mimic_B_MCAR_MNARc.csv']])
q = pd.concat([R(f) for f in ['mimic_A_MCAR_MAR_queries.csv', 'mimic_A_MNAR_MNARc_queries.csv', 'mimic_A_N1000_queries.csv', 'mimic_B_MCAR_MNARc_queries.csv']])
tiers = pd.concat([R(f) for f in ['mimic_A_MCAR_MAR_tiers.csv', 'mimic_A_MNAR_MNARc_tiers.csv', 'mimic_A_N1000_tiers.csv', 'mimic_B_MCAR_MNARc_tiers.csv']])
cfg = pat[['N', 'mechanism', 'rate', 'seed', 'scenario']].drop_duplicates()
print(f'\nMIMIC: {len(cfg)} configurations')
c = pat[pat.method == 'certified']
assert (c.prec_int.fillna(1) == 1).all() and (c.prec_con.fillna(1) == 1).all() and c.covered.astype(bool).all()
assert (q.err_where_cert_decides == 0).all()
print('  certified: precision 1, intervals cover; imputation error where certification decides = 0 in every configuration')
for N in [300, 1000]:
    d = pat[(pat.N == N) & (pat.scenario == 'A')]
    print(f'  N={N}: true intents (mean over configurations) {d.n_true_intents.mean():.0f} of 4096 = {d.n_true_intents.mean()/4096:.3f}')
    for mech in sorted(d.mechanism.unique()):
        dd = d[d.mechanism == mech]
        for rate in [0.10, 0.20]:
            pr = tab(dd, rate, 'prec_int', IMP + ['certified']); rc = tab(dd, rate, 'rec_int', IMP + ['certified'])
            print(f'    {mech} {int(rate*100)}% prec ' + ' '.join(f'{m}={f3(v)}' for m, v in pr.items()))
            print(f'             rec  ' + ' '.join(f'{m}={f3(v)}' for m, v in rc.items()))
    dn = d[(d.method.isin(IMP)) & (d.rate > 0)]
    print(f'    N={N}: lowest average precision of an imputation method over all settings {f3(dn.groupby(["mechanism", "rate", "method"]).prec_int.mean().min())}'
          f' (lowest single run {f3(dn.prec_int.min())})')
print('  query level (N = 1000, scenario A):')
d = q[(q.N == 1000) & (q.scenario == 'A')]
for mech in ['MCAR', 'MAR', 'MNAR', 'MNARc']:
    for rate in [0.05, 0.10, 0.20, 0.30]:
        dd = d[(d.mechanism == mech) & np.isclose(d.rate, rate)]
        ab = dd[dd.method == 'mice'].abstain_cert.mean()
        e_all = {m: dd[dd.method == m].err_all.mean() for m in IMP}
        e_ab = {m: dd[dd.method == m].err_where_cert_abstains.mean() for m in IMP}
        print(f'    {mech:5s} {int(rate*100):2d}%: open {f3(ab)} | wrong overall ' + ' '.join(f'{m}={f3(v)}' for m, v in e_all.items())
              + ' | wrong among open ' + ' '.join(f'{m}={f3(v)}' for m, v in e_ab.items()))
dd = d[(d.mechanism == 'MCAR') & np.isclose(d.rate, 0.10)]
print('    MCAR 10%: support relative error of imputation', {m: f3(dd[dd.method == m].supp_relerr.mean()) for m in IMP},
      '| certified interval width / true support', f3(dd[dd.method == 'mice'].cert_width.mean()))
print('    patterns per query set (mean):', d.n_patterns.mean())
db = pat[(pat.scenario == 'B')]
for mech in ['MCAR', 'MNARc']:
    for rate in [0.0, 0.10]:
        dd = db[db.mechanism == mech]
        pr = tab(dd, rate, 'prec_int', IMP + ['certified', 'certified_noband']); rc = tab(dd, rate, 'rec_int', IMP + ['certified', 'certified_noband'])
        print(f'  B {mech} {int(rate*100)}% prec ' + ' '.join(f'{m}={f3(v)}' for m, v in pr.items()))
        print(f'              rec  ' + ' '.join(f'{m}={f3(v)}' for m, v in rc.items()))
tt = tiers[(tiers.N == 300) & (tiers.scenario == 'A') & (tiers.mechanism == 'MCAR') & np.isclose(tiers.rate, 0.10)]
for m in ['mice', 'knn']:
    x = tt[tt.method == m]
    print(f'  tiers N=300 MCAR 10% {m}: certified share {f3(x.share_tier1.mean())}, precision of the rest {f3(x.prec_tier2.mean())}')
ca = pat[(pat.method == 'certified') & (pat.scenario == 'A') & (pat.N == 1000) & (pat.mechanism == 'MCAR') & np.isclose(pat.rate, 0.10)]
print('  N=1000 MCAR 10% certified-intent interval width (frequent recovered intents):', f3(ca.interval_width.mean()))
ca = pat[(pat.method == 'certified') & (pat.scenario == 'A') & (pat.N == 300) & (pat.mechanism == 'MCAR') & np.isclose(pat.rate, 0.10)]
print('  N=300 MCAR 10% certified-intent interval width:', f3(ca.interval_width.mean()))
im = pat[(pat.method.isin(IMP)) & (pat.scenario == 'A') & (pat.N == 300) & (pat.mechanism == 'MCAR') & np.isclose(pat.rate, 0.10)]
print('  N=300 MCAR 10% Jaccard support error of imputation:', {m: f3(im[im.method == m].supp_err.mean()) for m in IMP})
