"""make_figure.py -- draws Figure 2 of the paper (the running example of Section 6):
(a) sure profiles and robustness of the sure concepts, (b) the Dowker barcode.
Run:  python3 make_figure.py figures/fig_crown.pdf
      FIG_LANG=es python3 make_figure.py figures/fig_crown_es.pdf     (Spanish labels)
Needs ptb.py, gudhi and matplotlib; also writes a PNG next to the PDF."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from ptb import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

S = [[Fr(4), Fr(-1,2), Fr(3,2)], [Fr(1), Fr(7,2), Fr(-5,2)], [Fr(-2), Fr(3), Fr(9,2)]]
P = ['p_1','p_2','p_3']; A = ['a_1','a_2','a_3']
BLUE, ORANGE, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#d9d8d4'
TMAX = 5.2
ES = os.environ.get('FIG_LANG') == 'es'
T = (lambda en, es: es if ES else en)

def lab(X, B):
    x = (r'\{' + ','.join(P[i] for i in sorted(X)) + r'\}') if X else r'\varnothing'
    b = (r'\{' + ','.join(A[i] for i in sorted(B)) + r'\}') if B else r'\varnothing'
    return '$(' + x + ',' + b + ')$'

sure = existence_sets(S, 'sure')
order = sorted(sure, key=lambda k: (-(float(rho(S, *k)) if rho(S, *k) > 0 else -1), -len(k[0]) if rho(S,*k)>0 else 0, sorted(k[0]), sorted(k[1])))
crisp = [k for k in order if rho(S, *k) > 0]
late = [k for k in order if rho(S, *k) <= 0]
rows = crisp + late

plt.rcParams.update({'font.size': 8.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': INK, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={'width_ratios': [1.35, 1]})

for i, k in enumerate(rows):
    y = len(rows) - 1 - i
    for (lo, hi) in sure[k]:
        hi = TMAX if hi == float('inf') else float(hi)
        ax.plot([float(lo), hi], [y, y], color=BLUE, lw=2.2, solid_capstyle='round',
                label=T('sure concept','concepto seguro') if i == 0 else None)
        if hi == TMAX:
            ax.annotate('', xy=(TMAX + 0.12, y), xytext=(TMAX - 0.05, y),
                        arrowprops=dict(arrowstyle='->', color=BLUE, lw=1.2))
    r = rho(S, *k)
    if r > 0:
        ax.plot([0, float(r)], [y + 0.28, y + 0.28], color=ORANGE, lw=2.2, solid_capstyle='round',
                label=T(r'certified: $t\leq\rho$', r'certificado: $t\leq\rho$') if i == 0 else None)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([lab(*k) for k in reversed(rows)])
ax.axhline(len(late) - 0.5, color=GRID, lw=0.8)
ax.set_xlim(0, TMAX + 0.2); ax.set_ylim(-0.7, len(rows) - 0.2)
ax.set_xlabel(T('tolerance $t$','tolerancia $t$'))
ax.set_xticks([0, 0.5, 1, 1.5, 2, 3, 4, 4.5])
ax.set_xticklabels(['0', '1/2', '1', '3/2', '2', '3', '4', '9/2'])
ax.grid(axis='x', color=GRID, lw=0.5); ax.set_axisbelow(True)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.set_title(T('(a) sure profiles and robustness','(a) perfiles seguros y robustez'), fontsize=9, loc='left', color=INK)
ax.legend(loc='upper right', bbox_to_anchor=(1.0, 0.62), frameon=True, facecolor='white', edgecolor=GRID, fontsize=7.5)

bars = dowker_persistence(S)
bars = sorted(bars, key=lambda z: (z[0], -z[1]))
UMIN = -1.3
for j, (dim, ub, ud) in enumerate(bars):
    y = len(bars) - 1 - j
    lo = UMIN if ud == -float('inf') else ud
    c = BLUE if dim == 0 else ORANGE
    bx.plot([lo, ub], [y, y], color=c, lw=2.2, solid_capstyle='round')
    if ud == -float('inf'):
        bx.annotate('', xy=(UMIN - 0.15, y), xytext=(UMIN + 0.05, y), arrowprops=dict(arrowstyle='->', color=c, lw=1.2))
    bx.text(ub + 0.12, y, f'$H_{dim}$', va='center', fontsize=8, color=INK)
bx.axvspan(-0.5, 0.5, color=ORANGE, alpha=0.10, lw=0)
bx.axvline(0, color=MUTED, lw=0.6)
bx.text(0, len(bars) - 0.35, r'$u=\pm t,\ t=1/2$', ha='center', fontsize=7.5, color=MUTED)
bx.set_yticks([]); bx.set_ylim(-0.7, len(bars) - 0.1)
bx.set_xlim(UMIN - 0.25, 5.1)
bx.set_xlabel(T('level $u$  (sure side $u>0$, possible side $u<0$)','nivel $u$  (lado seguro $u>0$, lado posible $u<0$)'))
bx.set_xticks([-0.5, 0, 1, 1.5, 3, 3.5, 4, 4.5]); bx.set_xticklabels(['$-1/2$', '0', '1', '3/2', '3', '7/2', '4', '9/2'], fontsize=7.5)
for sp in ('top', 'right', 'left'): bx.spines[sp].set_visible(False)
bx.grid(axis='x', color=GRID, lw=0.5); bx.set_axisbelow(True)
bx.set_title(T('(b) Dowker barcode','(b) código de barras de Dowker'), fontsize=9, loc='left', color=INK)
fig.tight_layout()
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'fig_crown.pdf')
fig.savefig((sys.argv[1] if len(sys.argv) > 1 else 'fig_crown.pdf').replace('.pdf', '.png'), dpi=200)
print('rows:', [lab(*k) for k in rows]); print('dowker:', bars)
