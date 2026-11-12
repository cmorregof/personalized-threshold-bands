"""make_table_results.py -- Table 3 (results at 10% MCAR) and its per-dataset supplement, from results/.

Every cell is a mean +- sample standard deviation over the ten repetitions. For the public data each
repetition is first averaged over the three datasets. The 'absent' baseline (every missing cell
declared absent, the least completion) answers the patient-pattern questions in mimic_abstain.csv
(design 'main', which reuses the samples and masks of run_mimic.py with N = 1000).
Writes results/table3.tex and results/table3_datasets.tex and prints the plain numbers.
Run:  python3 make_table_results.py
"""
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda f: pd.read_csv(os.path.join(HERE, 'results', f))
OUT = lambda f: os.path.join(HERE, 'results', f)
M = ['mean', 'median', 'knn', 'mice', 'absent', 'certified']
NAMES = dict(mean='Mean', median='Median', knn='kNN', mice='MICE', absent='Absent (least completion)',
             certified='Certification')
ten = lambda x: np.isclose(x, 0.10)

P = R('public_runs.csv'); P = P[(P.scenario == 'A') & (P.mechanism == 'MCAR') & ten(P.rate)]
Q = R('mimic_A_MCAR_MAR.csv'); Q = Q[(Q.N == 300) & (Q.mechanism == 'MCAR') & ten(Q.rate)]
U = R('mimic_A_N1000_queries.csv'); U = U[(U.mechanism == 'MCAR') & ten(U.rate)]
A = R('mimic_abstain.csv')
A = A[(A.design == 'main') & (A.truth == 'obs') & (A.mechanism == 'MCAR') & ten(A.rate) & (A.method == 'absent')]


def ms(x):
    x = np.asarray(x, float)
    return x.mean(), x.std(ddof=1)


def pub(m, col):          # per-seed average over the three datasets, then mean and SD over seeds
    return ms(P[P.method == m].groupby('seed')[col].mean())


f2 = lambda t: f'${t[0]:.2f}\\pm{t[1]:.2f}$'
pc = lambda t: f'${100 * t[0]:.1f}\\pm{100 * t[1]:.1f}\\,\\%$'
rows = []
for m in M:
    cells = [f2(pub(m, 'prec_int')), f2(pub(m, 'rec_int')),
             f2(ms(Q[Q.method == m].prec_int)), f2(ms(Q[Q.method == m].rec_int))]
    if m in ('mean', 'median', 'knn', 'mice'):
        u = U[U.method == m]
        cells += [pc(ms(u.err_all)), pc(ms(u.err_where_cert_abstains))]
    elif m == 'absent':
        cells += [pc(ms(A.risk)), pc(ms(A.risk_open))]
    else:
        cells += ['$0$', '--']
    rows.append((m, cells))
open_share = ms(U[U.method == 'mice'].abstain_cert)

body = '\n'.join(('\\addlinespace[2pt]\n' if m == 'absent' else '') + NAMES[m] + ' & ' + ' & '.join(c) + r'\\'
                 for m, c in rows)
tex = r"""\begin{table}[t]
\centering
\footnotesize
\setlength{\tabcolsep}{2.4pt}
\begin{tabular}{lcccccc}
\toprule
 & \multicolumn{2}{c}{Public data} & \multicolumn{2}{c}{MIMIC-IV, $N=300$} & \multicolumn{2}{c}{MIMIC-IV, $N=1\,000$}\\
 & \multicolumn{2}{c}{patterns} & \multicolumn{2}{c}{patterns} & \multicolumn{2}{c}{questions}\\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}
 & Precision & Recall & Precision & Recall & Wrong, all & Wrong, open\\
\midrule
""" + body + r"""
\bottomrule
\end{tabular}
\caption{Results with $10\%$ of the values missing completely at random, as mean $\pm$ standard deviation over the ten repetitions. For the public data each repetition is first averaged over the three datasets. ``Wrong, all'' is the share of all patient--pattern questions answered wrongly. ``Wrong, open'' is the share of wrong answers among the questions that certification leaves open, which are """ + f'${100 * open_share[0]:.1f}\\pm{100 * open_share[1]:.1f}\\%$' + r""" of the total. ``Absent'' declares every missing cell absent, which gives the least completion. Patterns are shown for samples of $300$ admissions, because with $1\,000$ the precision of every method exceeds $0.98$ on average. Certification never answers wrongly, and every wrong answer of imputation falls on an open question (Corollary~\ref{cor:imputation}).}
\label{tab:results}
\end{table}
"""
open(OUT('table3.tex'), 'w').write(tex)

D = [('wine', 'Wine'), ('breast_cancer', 'Breast cancer'), ('diabetes', 'Diabetes')]
lines = []
for m in M:
    cells = []
    for d, _ in D:
        x = P[(P.method == m) & (P.dataset == d)]
        cells += [f2(ms(x.prec_int)), f2(ms(x.rec_int))]
    lines.append(('\\addlinespace[2pt]\n' if m == 'absent' else '') + NAMES[m] + ' & ' + ' & '.join(cells) + r'\\')
tex2 = r"""\begin{table}[t]
\centering
\footnotesize
\setlength{\tabcolsep}{2.4pt}
\begin{tabular}{lcccccc}
\toprule
 & \multicolumn{2}{c}{Wine} & \multicolumn{2}{c}{Breast cancer} & \multicolumn{2}{c}{Diabetes}\\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}
 & Precision & Recall & Precision & Recall & Precision & Recall\\
\midrule
""" + '\n'.join(lines) + r"""
\bottomrule
\end{tabular}
\caption{Public data with $10\%$ of the values missing completely at random, by dataset, as mean $\pm$ standard deviation over the ten repetitions.}
\label{tab:results-datasets}
\end{table}
"""
open(OUT('table3_datasets.tex'), 'w').write(tex2)
for m, c in rows:
    print(f'{NAMES[m]:28s}', ' | '.join(c))
print('open questions:', open_share)
print('-> results/table3.tex, results/table3_datasets.tex')
