"""make_table_results.py -- the numbers of Table 3 (results at 10% MCAR), from results/."""
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda f: pd.read_csv(os.path.join(HERE, 'results', f))
M = ['mean', 'median', 'knn', 'mice', 'certified']
P = R('public_runs.csv'); P = P[(P.scenario == 'A') & (P.mechanism == 'MCAR') & np.isclose(P.rate, 0.10)]
Q = pd.concat([R('mimic_A_MCAR_MAR.csv')]); Q = Q[(Q.N == 300) & (Q.mechanism == 'MCAR') & np.isclose(Q.rate, 0.10)]
U = R('mimic_A_N1000_queries.csv'); U = U[(U.mechanism == 'MCAR') & np.isclose(U.rate, 0.10)]
f2 = lambda x: f'{x:.2f}'
pc = lambda x: f'{100 * x:.1f}\\%'
row = lambda name, vals: name + ' & ' + ' & '.join(vals) + r'\\'
print(row('PUB precision', [f2(P[P.method == m].prec_int.mean()) for m in M]))
print(row('PUB recall', [f2(P[P.method == m].rec_int.mean()) for m in M]))
print(row('MIMIC precision', [f2(Q[Q.method == m].prec_int.mean()) for m in M]))
print(row('MIMIC recall', [f2(Q[Q.method == m].rec_int.mean()) for m in M]))
print(row('wrong all', [pc(U[U.method == m].err_all.mean()) for m in M[:4]] + ['0']))
print(row('open', ['--'] * 4 + [pc(U[U.method == 'mice'].abstain_cert.mean())]))
print(row('wrong open', [pc(U[U.method == m].err_where_cert_abstains.mean()) for m in M[:4]] + ['--']))
print('exact values:', {m: Q[Q.method == m].rec_int.mean() for m in M})
