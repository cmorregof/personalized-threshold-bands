"""make_figure_lattices.py -- draws Figure 1 of the paper (lattices and holes of the crown):
(a) the concept lattice of the crisp context, (b) its Dowker complex, (c)-(e) the certified
intents at three tolerances. Every node and edge is recomputed from the table with ptb.py.
Run:  python3 make_figure_lattices.py figures/fig_lattices.pdf
      FIG_LANG=es python3 make_figure_lattices.py figures/fig_lattices_es.pdf   (Spanish labels)
Also writes a PNG next to the PDF."""
import sys, os
from itertools import permutations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from ptb import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

S = [[Fr(4), Fr(-1, 2), Fr(3, 2)], [Fr(1), Fr(7, 2), Fr(-5, 2)], [Fr(-2), Fr(3), Fr(9, 2)]]
nG, nM = 3, 3
BLUE, ORANGE, INK, MUTED, GRID, SURF = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#c3c2b7', '#ffffff'
ES = os.environ.get('FIG_LANG') == 'es'
T = (lambda en, es: es if ES else en)
P = lambda g: f'$p_{g + 1}$'
A = lambda m: f'$a_{m + 1}$'


def intent_label(B):
    if not B:
        return r'$\varnothing$'
    return '$' + ''.join(f'a_{m + 1}' for m in sorted(B)) + '$'


def certified_intents(t):
    return [B for B in powerset(nM) if t <= rho_int(S, B)]


def covers(sets):
    """Covering pairs (smaller, larger) of a family of frozensets ordered by inclusion."""
    out = []
    for a in sets:
        for b in sets:
            if a < b and not any(a < c < b for c in sets):
                out.append((a, b))
    return out


plt.rcParams.update({'font.size': 8.5, 'font.family': 'serif', 'mathtext.fontset': 'cm',
                     'axes.edgecolor': MUTED, 'text.color': INK})
W, H = 24.0, 8.6                     # one canvas in data units, 0.3 inch per unit
fig = plt.figure(figsize=(0.3 * W, 0.3 * H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect('equal'); ax.set_axis_off()
RAD, TITLE_Y = 0.24, 8.15


def node(x, y, color):
    ax.add_patch(plt.Circle((x, y), RAD, facecolor=SURF, edgecolor=color, lw=1.4, zorder=3))


def edge(p, q, color, lw=1.1, z=1):
    ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, solid_capstyle='round', zorder=z)


def title(x, text):
    ax.text(x, TITLE_Y, text, ha='left', va='center', fontsize=8.5, color=INK)


# ---------------------------------------------------------------- (a) the concept lattice at t <= 1/2
R = sure_rel(S, Fr(1, 2))
assert R == poss_rel(S, Fr(1, 2))                    # complete context: one completion
C = sorted(concepts(R), key=lambda c: (-len(c[0]), sorted(c[0])))
top = [c for c in C if len(c[0]) == 3]; bot = [c for c in C if len(c[0]) == 0]
mid2 = [c for c in C if len(c[0]) == 2]; mid1 = [c for c in C if len(c[0]) == 1]
assert len(top) == len(bot) == 1 and len(mid2) == len(mid1) == 3


def ncross(o2, o1):
    pos2 = {c: i for i, c in enumerate(o2)}; pos1 = {c: i for i, c in enumerate(o1)}
    es = [(pos2[a], pos1[b]) for a in o2 for b in o1 if b[0] < a[0]]
    return sum(1 for i, (a, b) in enumerate(es) for (c, d) in es[i + 1:] if (a - c) * (b - d) < 0)


o2, o1 = min(((p2, p1) for p2 in permutations(mid2) for p1 in permutations(mid1)),
             key=lambda z: (ncross(*z), [sorted(c[0]) for c in z[0]]))
X0 = 3.4
pos = {top[0]: (X0, 7.0), bot[0]: (X0, 1.0)}
for i, c in enumerate(o2):
    pos[c] = (X0 - 2.0 + 2.0 * i, 5.0)
for i, c in enumerate(o1):
    pos[c] = (X0 - 2.0 + 2.0 * i, 3.0)
for c in C:
    for d in C:
        if c[0] < d[0] and not any(c[0] < e[0] < d[0] for e in C):
            crown = len(c[0]) == 1 and len(d[0]) == 2
            edge(pos[c], pos[d], BLUE if crown else GRID, lw=1.6 if crown else 1.0, z=2 if crown else 1)
for c in C:
    node(*pos[c], BLUE if len(c[0]) in (1, 2) else MUTED)
# reduced labelling: each attribute above its attribute concept, each object below its object concept
for m in range(nM):
    ext = down(R, frozenset([m]), nG)
    x, y = pos[(ext, up(R, ext, nM))]
    ax.text(x, y + 0.48, A(m), ha='center', va='bottom', fontsize=8.5)
for g in range(nG):
    itt = up(R, frozenset([g]), nM)
    x, y = pos[(down(R, itt, nG), itt)]
    ax.text(x, y - 0.48, P(g), ha='center', va='top', fontsize=8.5)
title(0.4, T(r'(a) concept lattice, $t\leq1/2$', r'(a) retículo de conceptos, $t\leq1/2$'))

# ---------------------------------------------------------------- (b) the Dowker complex at t <= 1/2
CX, CY, RT = 10.3, 3.9, 2.1
V = {0: (CX, CY + RT), 1: (CX - 0.866 * RT, CY - 0.5 * RT), 2: (CX + 0.866 * RT, CY - 0.5 * RT)}
cx = dowker_complex(R)
assert frozenset([0, 1, 2]) not in cx and all(frozenset(e) in cx for e in [(0, 1), (0, 2), (1, 2)])
for (g, h) in [(0, 1), (0, 2), (1, 2)]:
    shared = [m for m in range(nM) if R[g][m] and R[h][m]]
    assert len(shared) == 1
    edge(V[g], V[h], BLUE, lw=1.6, z=2)
    mx, my = (V[g][0] + V[h][0]) / 2, (V[g][1] + V[h][1]) / 2
    dx, dy = mx - CX, my - CY
    nrm = (dx * dx + dy * dy) ** 0.5
    ax.text(mx + 0.55 * dx / nrm, my + 0.55 * dy / nrm, A(shared[0]), ha='center', va='center', fontsize=8.5)
for g, (x, y) in V.items():
    node(x, y, BLUE)
    if g == 0:
        ax.text(x, y + 0.48, P(g), ha='center', va='bottom', fontsize=8.5)
    else:
        ax.text(x + (-0.5 if g == 1 else 0.5), y - 0.35, P(g), ha='center', va='top', fontsize=8.5)
ax.text(CX, CY - 0.05, T('hollow', 'hueco'), ha='center', va='center', fontsize=7.5, color=MUTED, style='italic')
title(7.6, T('(b) Dowker complex', '(b) complejo de Dowker'))

# ---------------------------------------------------------------- (c)-(e) certified intents
panels = [(15.6, Fr(1), r'(c) $t\in(1/2,1]$', 14.0), (19.3, Fr(2), r'(d) $t\in(1,2]$', 18.3),
          (21.95, Fr(3), r'(e) $t>2$', 21.3)]
expected = {Fr(1): 4, Fr(2): 3, Fr(3): 1}
for xc, t, ttl, xt in panels:
    I = certified_intents(t)
    assert len(I) == expected[t]
    levels = sorted({len(B) for B in I})
    ys = [7.0, 4.0, 1.0] if len(levels) == 3 else [4.0]
    p = {}
    for k, y in zip(levels, ys):
        row = sorted([B for B in I if len(B) == k], key=sorted)
        for j, B in enumerate(row):
            p[B] = (xc - 0.9 * (len(row) - 1) + 1.8 * j, y)
    for a, b in covers(I):
        edge(p[a], p[b], ORANGE, lw=1.4, z=1)
    for B in I:
        x, y = p[B]
        node(x, y, ORANGE)
        left = x < xc - 0.01
        ax.text(x + (-0.45 if left else 0.45), y, intent_label(B), ha='right' if left else 'left',
                va='center', fontsize=8)
    title(xt, ttl)
ax.text(18.65, 0.12, T('certified intents as the tolerance grows', 'intensiones certificadas al crecer la tolerancia'),
        fontsize=7.5, color=MUTED, ha='center', va='bottom', style='italic')
out = sys.argv[1] if len(sys.argv) > 1 else 'figures/fig_lattices.pdf'
fig.savefig(out)
fig.savefig(out.replace('.pdf', '.png'), dpi=220)
print('written', out)
