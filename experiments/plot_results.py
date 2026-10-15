"""plot_results.py -- Figures 3 and 4 of the paper, from the aggregate tables in results/.
Run:  python3 plot_results.py ../figures            (English: fig_public.pdf, fig_mimic.pdf)
      FIG_LANG=es python3 plot_results.py ../figures (Spanish: fig_public_es.pdf, fig_mimic_es.pdf)
Categorical colours in fixed order; the palette passes the colour-vision checks, and marker
shapes repeat the identity of every series."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda f: pd.read_csv(os.path.join(HERE, 'results', f))
ES = os.environ.get('FIG_LANG') == 'es'
T = (lambda en, es: es if ES else en)
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
STYLE = {  # method: (label, colour, marker)
    'certified': (T('Certification', 'Certificación'), '#2a78d6', 'o'),
    'mice': ('MICE', '#eb6834', 's'),
    'knn': ('kNN', '#1baf7a', 'D'),
    'mean': (T('Mean', 'Media'), '#eda100', '^'),
    'median': (T('Median', 'Mediana'), '#e87ba4', 'v'),
}
plt.rcParams.update({'font.size': 8.5, 'font.family': 'serif', 'mathtext.fontset': 'cm',
                     'axes.edgecolor': AXIS, 'axes.linewidth': 0.8, 'axes.labelcolor': INK2,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK})


def style_axes(ax):
    ax.grid(axis='y', color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=0, pad=3)


def line(ax, x, y, key, label=None):
    lab, col, mk = STYLE[key]
    ax.plot(x, y, color=col, lw=1.8, marker=mk, ms=5.5, mec='white', mew=1.2,
            solid_capstyle='round', solid_joinstyle='round', label=label or lab,
            zorder=3 if key == 'certified' else 2)


def fig_public(out):
    P = R('public_runs.csv')
    A = P[P.scenario == 'A']
    methods = ['certified', 'mice', 'knn', 'mean', 'median']
    titles = {'MCAR': T('MCAR: completely at random', 'MCAR: completamente al azar'),
              'MAR': T('MAR: at random', 'MAR: al azar'),
              'MNAR': T('MNAR: not at random', 'MNAR: no al azar')}
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.1), sharex=True, sharey=True)
    for j, mech in enumerate(['MCAR', 'MAR', 'MNAR']):
        d = A[A.mechanism == mech]
        for i, (metric, ylab) in enumerate([('prec_int', T('Precision', 'Precisión')),
                                           ('rec_int', T('Recall', 'Recuperación'))]):
            ax = axes[i, j]
            for m in methods:
                s = d[d.method == m].groupby('rate')[metric].mean()
                line(ax, 100 * s.index, s.values, m)
            style_axes(ax)
            ax.set_ylim(0, 1.04); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            if ES:  # decimal comma in the Spanish figure
                ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda y, _: f'{y:.2f}'.replace('.', ',')))
            ax.set_xticks([0, 5, 10, 20, 30])
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(titles[mech], fontsize=8.5, color=INK, loc='left')
            else:
                ax.set_xlabel(T('Missing values (%)', 'Datos faltantes (%)'))
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.0),
               handlelength=2.0, columnspacing=1.4, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.93), h_pad=0.8, w_pad=0.6)
    fig.savefig(out)
    fig.savefig(out.replace('.pdf', '.png'), dpi=220)
    plt.close(fig)


def fig_mimic(out):
    Q = R('mimic_A_N1000_queries.csv')
    mechs = [('MCAR', T('MCAR: completely at random', 'MCAR: completamente al azar')),
             ('MNARc', T('Clinical MNAR', 'MNAR clínico'))]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.6), sharex=True)
    for i, (mech, title) in enumerate(mechs):
        g = Q[Q.mechanism == mech].groupby(['rate', 'method'])[
            ['abstain_cert', 'err_all', 'err_where_cert_abstains']].mean().reset_index()
        rates = sorted(r for r in g.rate.unique() if r > 0)
        x = [0] + [100 * r for r in rates]
        ax = axes[i, 0]
        ab = g[g.method == 'mice'].set_index('rate').abstain_cert
        line(ax, x, [0] + [100 * ab[r] for r in rates], 'certified',
             label=T('Certification: questions left open', 'Certificación: preguntas abiertas'))
        for m in ('mice', 'knn', 'mean', 'median'):
            e = g[g.method == m].set_index('rate').err_all
            line(ax, x, [0] + [100 * e[r] for r in rates], m,
                 label=STYLE[m][0] + T(': wrong answers', ': respuestas erradas'))
        ax.set_ylim(0, 27)
        ax.set_ylabel(T('% of all questions', '% de todas las preguntas'))
        ax.set_title(title, fontsize=8.5, color=INK, loc='left')
        style_axes(ax)
        ax = axes[i, 1]
        for m in ('mice', 'knn', 'mean', 'median'):
            e = g[g.method == m].set_index('rate').err_where_cert_abstains
            line(ax, [100 * r for r in rates], [100 * e[r] for r in rates], m)
        ax.axhline(50, color=AXIS, lw=0.9)
        ax.text(0.6, 51.5, T('chance (50%)', 'azar (50 %)'), fontsize=7.5, color=MUTED, va='bottom', ha='left')
        ax.set_ylim(0, 56)
        ax.set_ylabel(T('% of the open questions', '% de las preguntas abiertas'))
        ax.set_title(T('Wrong answers of imputation on the open questions',
                       'Errores de la imputación en las preguntas abiertas'), fontsize=8.5, color=INK, loc='left')
        style_axes(ax)
    for ax in axes[1]:
        ax.set_xlabel(T('Missing values (%)', 'Datos faltantes (%)'))
        ax.set_xticks([0, 5, 10, 20, 30])
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=7.8,
               columnspacing=1.2, handlelength=2.0)
    fig.tight_layout(rect=(0, 0, 1, 0.9), h_pad=1.0, w_pad=1.2)
    fig.savefig(out)
    fig.savefig(out.replace('.pdf', '.png'), dpi=220)
    plt.close(fig)


if __name__ == '__main__':
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    sfx = '_es' if ES else ''
    fig_public(os.path.join(outdir, f'fig_public{sfx}.pdf'))
    fig_mimic(os.path.join(outdir, f'fig_mimic{sfx}.pdf'))
    print('figures written to', outdir)
