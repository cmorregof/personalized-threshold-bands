"""plot_tsweep.py -- the figure of experiment E1 (fig_tsweep), from the tables written by run_tsweep.py.
Run:  python3 plot_tsweep.py ../figures/fig_tsweep.pdf [RESULTS_DIR]     (FIG_LANG=es for Spanish)
Left: against re-measurement, 10% MCAR, personal bands: abstention and error of certification
where it decides, against the tolerance t, with the error of MICE on the questions certification
leaves open. Right: share of the intents certified at t = 0.25 that are intents of the
re-measured context, by decile of rho_int, for both band variants (mean and SD over repetitions)."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_results import STYLE, INK, MUTED, T, style_axes, line

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]
RES = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'results')
q = pd.read_csv(os.path.join(RES, 'mimic_T_sweep_queries.csv'))
h = pd.read_csv(os.path.join(RES, 'mimic_T_rho_int.csv'))
sel = lambda d: d[(d.truth == 'remeasure') & (d.mechanism == 'MCAR') & np.isclose(d.rate, 0.10)]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 2.9))
d = sel(q)
d = d[(d.band == 'personal') & (d.method == 'mice')].groupby('t').mean(numeric_only=True)
lab, col, mk = STYLE['certified']
a1.plot(d.index, 100 * d.abstain_cert, color=col, lw=1.8, ls='--', marker=mk, ms=5.5, mfc='white',
        mec=col, mew=1.2, label=T('Certification: abstains', 'Certificación: se abstiene'))
line(a1, d.index, 100 * d.cert_err_decided, 'certified', T('Certification: errs where it decides',
                                                         'Certificación: yerra donde decide'))
line(a1, d.index, 100 * d.err_where_cert_abstains, 'mice', T('MICE: errs on the open questions',
                                                            'MICE: yerra en las preguntas abiertas'))
a1.set_xlabel(T('tolerance $t$', 'tolerancia $t$'))
a1.set_ylabel('%')
a1.set_ylim(0, 100)
a1.set_title(T('Patient–pattern questions', 'Preguntas paciente–patrón'), fontsize=8.5, color=INK, loc='left')
a1.legend(frameon=False, fontsize=7, loc='upper left')
style_axes(a1)

e = sel(h)
for band, key, name in (('personal', 'certified', T('personal bands', 'bandas personales')),
                        ('per_test', 'mean', T('one band per test', 'una banda por examen'))):
    g = e[e.band == band].groupby('decile')
    m, s = g.prec.mean(), g.prec.std()
    _, col, mk = STYLE[key]
    a2.errorbar(m.index, 100 * m, yerr=100 * s, color=col, lw=1.8, marker=mk, ms=5.5, mec='white',
                mew=1.2, capsize=2, elinewidth=0.8, label=name)
a2.set_xlabel(T(r'decile of $\rho_{\mathrm{int}}$ among the intents certified at $t=0.25$',
                r'decil de $\rho_{\mathrm{int}}$ entre las intensiones certificadas en $t=0.25$'))
a2.set_ylabel('%')
a2.set_xticks(range(1, 11))
a2.set_ylim(top=101)
a2.set_title(T('Certified intents that are re-measured intents', 'Intensiones certificadas que lo son al re-medir'),
             fontsize=8.5, color=INK, loc='left')
a2.legend(frameon=False, fontsize=7, loc='lower right')
style_axes(a2)
fig.tight_layout()
fig.savefig(OUT)
fig.savefig(OUT.replace('.pdf', '.png'), dpi=220)
print('->', OUT)
