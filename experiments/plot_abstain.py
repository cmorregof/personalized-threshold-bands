"""plot_abstain.py -- risk-coverage figure of experiment E2 (fig_risk_coverage), from results/mimic_abstain.csv.
Run:  python3 plot_abstain.py ../figures/fig_risk_coverage.pdf [RESULTS_DIR]     (FIG_LANG=es for Spanish)
Patient-pattern questions on MIMIC-IV, 10% MCAR, means over ten repetitions. Coverage is the share of
questions answered, risk the share of wrong answers among them. Left: the complete panel as truth
(samples of run_mimic.py, bands of width zero), where certification is one point. Right: the re-measured
panel as truth (samples of run_tsweep.py), where certification with bands traces a curve as t varies.
The abstaining competitors trace curves as their level alpha varies."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_results import STYLE, INK, MUTED, T, style_axes

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]
RES = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'results')
d = pd.read_csv(os.path.join(RES, 'mimic_abstain.csv'))
d = d[(d.mechanism == 'MCAR') & np.isclose(d.rate, 0.10)]
SER = {  # key: (label, colour, marker, linestyle)
    'mi': (T('Multiple imputation, $m=20$', 'Imputación múltiple, $m=20$'), STYLE['mice'][1], 's', '-'),
    'knn_boot': (T('kNN with bootstrap, $m=20$', 'kNN con bootstrap, $m=20$'), STYLE['knn'][1], 'D', '-'),
    'mice': ('MICE', STYLE['mice'][1], 's', 'none'),
    'absent': (T('Absent', 'Ausente'), MUTED, 'x', 'none'),
}


def curve(ax, x, key):
    lab, col, mk, ls = SER[key]
    g = x[x.method == key].groupby('alpha')[['coverage', 'risk']].mean().sort_values('coverage')
    ax.plot(100 * g.coverage, 100 * g.risk, color=col, ls=ls, lw=1.6, marker=mk, ms=5,
            mfc='white' if ls == 'none' else col, mec=col, mew=1.2, label=lab)


fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
for ax, design, truth, title in ((axes[0], 'main', 'obs', T('Truth: complete panel', 'Verdad: panel completo')),
                                 (axes[1], 'sweep', 'remeasure', T('Truth: re-measured panel', 'Verdad: panel re-medido'))):
    x = d[(d.design == design) & (d.truth == truth)]
    for key in ('mi', 'knn_boot', 'mice', 'absent'):
        curve(ax, x, key)
    c = x[x.method == 'certified']
    _, col, mk = STYLE['certified']
    if design == 'main':
        g = c[['coverage', 'risk']].mean()
        ax.plot(100 * g.coverage, 100 * g.risk, marker=mk, ms=7, color=col, ls='none',
                label=T('Certification', 'Certificación'))
    else:
        for band, ls, lab in (('personal', '-', T('Certification, personal bands', 'Certificación, bandas personales')),
                              ('per_test', '--', T('Certification, one band per test', 'Certificación, banda por examen'))):
            g = c[c.band == band].groupby('t')[['coverage', 'risk']].mean()
            ax.plot(100 * g.coverage, 100 * g.risk, color=col, ls=ls, lw=1.6, marker=mk, ms=4.5, label=lab)
            if band == 'personal':
                for t, r in g.iterrows():
                    ax.annotate(f'$t={t:g}$', (100 * r.coverage, 100 * r.risk), textcoords='offset points',
                                xytext=(3, 4), fontsize=6.5, color=MUTED)
    ax.set_xlabel(T('coverage: questions answered (%)', 'cobertura: preguntas respondidas (%)'))
    ax.set_ylabel(T('risk: wrong among answered (%)', 'riesgo: erradas entre las respondidas (%)'))
    ax.set_title(title, fontsize=8.5, color=INK, loc='left')
    ax.set_ylim(bottom=0)
    style_axes(ax)
    ax.legend(frameon=False, fontsize=6.5, loc='upper left')
fig.tight_layout()
fig.savefig(OUT)
fig.savefig(OUT.replace('.pdf', '.png'), dpi=220)
print('->', OUT)
