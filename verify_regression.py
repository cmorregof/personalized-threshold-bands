"""
verify_regression.py -- a broader randomized regression suite, on tie-heavy tables
with zeros, for the statements of Sections 2 and 5 of the paper (numbers as in the paper):
Prop. 2.3, Thm. 5.2, Prop. 5.4, Thm. 5.5, Thm. 5.7, Prop. 5.11, the profile formula
behind Prop. 5.12, Ex. 5.13, Prop. 5.15, Thm. 5.17 and Cor. 5.18.
Run:  python3 verify_regression.py      (needs ptb.py, numpy and gudhi)
"""
from fractions import Fraction as Fr
from itertools import product
import random, math
import numpy as np
import gudhi
from ptb import *

random.seed(7)
np.random.seed(7)

# ---------------------------------------------------------------- a fixed 5x5 test table with ties
T = Fr(1, 3)
S_TEST = [
    [Fr(1),      Fr(2),  -2 * T,  -4 * T,  Fr(-2)],   # P1
    [Fr(5, 2),   Fr(-2), Fr(1),    4 * T,  -T],       # P2
    [T,          T,      Fr(1),   -4 * T,   5 * T],   # P3
    [-T,         Fr(-2), Fr(1),    4 * T,  -5 * T],   # P4
    [4 * T,      Fr(2),  Fr(2),   -5 * T,  -2 * T],   # P5
]


def rand_matrix(nG, nM, vals):
    return [[Fr(random.choice(vals)) for _ in range(nM)] for _ in range(nG)]


def all_sigs(nG, nM):
    return [(X, B) for X in powerset(nG) for B in powerset(nM)]


# ---------------------------------------------------------------- Thm 5.2: closed-form profiles
def check_T1(s):
    nG, nM = len(s), len(s[0])
    sure = existence_sets(s, 'sure')
    poss = existence_sets(s, 'possible')
    for (X, B) in all_sigs(nG, nM):
        b, d = birth(s, X, B), death(s, X, B)
        # sure: (max(b,0), d]
        lo, hi = max(b, 0), d
        exp = [(lo, hi)] if hi > lo else []
        got = sure.get((X, B), [])
        assert got == exp, ('sure', X, B, got, exp)
        # possible: (max(-d,0), -b]
        lo, hi = max(-d, 0), -b
        exp = [(lo, hi)] if hi > lo else []
        got = poss.get((X, B), [])
        assert got == exp, ('poss', X, B, got, exp)
    return len(sure), len(poss)


print('Thm. 5.2  closed-form sure and possible profiles')
n1, n2 = check_T1(S_TEST)
print('    5x5 test table OK   #sure signatures =', n1, ' #possible signatures =', n2)
for trial in range(300):
    nG, nM = random.choice([(2, 2), (3, 3), (4, 3), (3, 4), (4, 4)])
    s = rand_matrix(nG, nM, [-3, -2, -1, 0, 1, 2, 3])        # many ties and zeros
    check_T1(s)
print('    300 random tie-heavy matrices (zeros included) OK')


# ---------------------------------------------------------------- Thm 5.7: certified concepts
def completions(Rs, Rp):
    nG, nM = len(Rs), len(Rs[0])
    free = [(g, m) for g in range(nG) for m in range(nM) if Rp[g][m] and not Rs[g][m]]
    for bits in product([False, True], repeat=len(free)):
        R = [list(r) for r in Rs]
        for (g, m), v in zip(free, bits):
            R[g][m] = v
        yield R


def check_T2(s, ts):
    nG, nM = len(s), len(s[0])
    for t in ts:
        Rs, Rp = sure_rel(s, t), poss_rel(s, t)
        nfree = sum(Rp[g][m] and not Rs[g][m] for g in range(nG) for m in range(nM))
        if nfree > 12:
            continue
        inter = None
        for R in completions(Rs, Rp):
            c = concepts(R)
            inter = c if inter is None else inter & c
        expected = {(X, B) for (X, B) in all_sigs(nG, nM) if 0 < t <= rho(s, X, B)}
        assert inter == expected, (t, inter ^ expected)


print('Thm. 5.7  concept of every completion  <=>  0 < t <= rho(X,B)')
for trial in range(120):
    nG, nM = random.choice([(2, 2), (3, 3), (3, 4), (4, 3)])
    s = rand_matrix(nG, nM, [-3, -2, -1, 0, 1, 2, 3])
    check_T2(s, [Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(5, 2), Fr(3), Fr(7, 2)])
print('    120 random matrices x 7 tolerances, all completions enumerated OK')


# ---------------------------------------------------------------- Prop 5.4 / Thm 5.5: stability without genericity
def perturb(s, eps):
    return [[float(x) + random.uniform(-eps, eps) for x in row] for row in s]


print('Prop. 5.4 / Thm. 5.5  inner and outer levels and robustness are 1-Lipschitz (no genericity)')
sf = [[float(x) for x in row] for row in S_TEST]
sigs = all_sigs(5, 5)
for eps in [0.05, 0.10, 0.1666, 0.20, 0.30, 0.50, 1.00]:
    worst, worst_gone = 0.0, 0.0
    for rep in range(200):
        s2 = perturb(S_TEST, eps)
        for (X, B) in sigs:
            b1, d1, b2, d2 = birth(sf, X, B), death(sf, X, B), birth(s2, X, B), death(s2, X, B)
            for a, c in ((b1, b2), (d1, d2)):
                if math.isfinite(a):
                    worst = max(worst, abs(a - c))
            r1, r2 = rho(sf, X, B), rho(s2, X, B)
            if math.isfinite(r1):
                worst = max(worst, abs(r1 - r2))
            l1, l2 = d1 - max(b1, 0), d2 - max(b2, 0)
            if (l1 > 0) != (l2 > 0):                     # created or destroyed sure concept
                worst_gone = max(worst_gone, max(l1, l2))
    assert worst <= eps + 1e-12 and worst_gone <= 2 * eps + 1e-12
    print(f'    eps={eps:<6}  max endpoint shift={worst:.4f} <= eps ;  '
          f'max lifespan of created/destroyed sure concept={worst_gone:.4f} <= 2eps')


# ---------------------------------------------------------------- Prop 2.3: SE-ISI = apposition
print('Prop. 2.3  SE-ISI extents = extents of the apposition  K^sure | K^poss')
for (lo, hi, t) in grid_pieces(S_TEST):
    Rs, Rp = sure_rel(S_TEST, t), poss_rel(S_TEST, t)
    e1 = {X for (X, B, C) in seisi_concepts(Rs, Rp)}
    e2 = {X for (X, B) in concepts(apposition(Rs, Rp))}
    assert e1 == e2
print('    OK on every piece of the critical grid of the test table')
bars = existence_sets(S_TEST, 'seisi')
nbars = sum(len(v) for v in bars.values())
multi = {k: v for k, v in bars.items() if len(v) > 1}
print(f'    test table: distinct SE-ISI signatures = {len(bars)}, maximal intervals = {nbars}, '
      f'signatures with more than one interval = {len(multi)}')
per_piece = [(float(lo), float(hi) if hi != INF else 'inf',
              len(seisi_concepts(sure_rel(S_TEST, t), poss_rel(S_TEST, t))))
             for (lo, hi, t) in grid_pieces(S_TEST)]
print('    concepts per piece (lo, hi, #):', per_piece)


# ---------------------------------------------------------------- Thm 5.17 / Cor 5.18: Dowker persistence
def dgm_by_dim(pers, dim, cap=None):
    out = []
    for (k, ub, ud) in pers:
        if k != dim:
            continue
        # work in (-level) coordinates so that births < deaths, as gudhi expects
        out.append((-ub, (INF if ud == -INF else -ud)))
    return out


print('Thm. 5.17 / Cor. 5.18  Dowker persistence: duality, stability, agreement with the concept poset')
pers_obj = dowker_persistence(sf, 'objects')
pers_att = dowker_persistence(sf, 'attributes')
for k in (0, 1, 2):
    dB = gudhi.bottleneck_distance(dgm_by_dim(pers_obj, k), dgm_by_dim(pers_att, k))
    assert dB < 1e-9, (k, dB)
print('    object-side and attribute-side barcodes coincide (Dowker duality) OK')
print('    bars with positive length (dim, u_birth, u_death):',
      sorted([(k, round(ub, 3), (round(ud, 3) if ud != -INF else '-inf'))
              for (k, ub, ud) in pers_obj if ub - (ud if ud != -INF else -1e9) > 1e-9]))
for eps in [0.05, 0.1666, 0.30, 0.50, 1.00]:
    worst = 0.0
    for rep in range(200):
        p2 = dowker_persistence(perturb(S_TEST, eps), 'objects')
        for k in (0, 1, 2):
            worst = max(worst, gudhi.bottleneck_distance(dgm_by_dim(pers_obj, k), dgm_by_dim(p2, k)))
    assert worst <= eps + 1e-9
    print(f'    eps={eps:<6}  max bottleneck distance over H0,H1,H2 = {worst:.4f} <= eps')

levels = sorted({x for row in S_TEST for x in row})
for u in levels:
    R = level_rel(S_TEST, u)
    D = dowker_complex(R)
    Q, chains = concept_poset_order_complex(R)
    bD = betti_of_complex(D) if D else []
    bQ = betti_of_complex(chains) if chains else []
    strip = lambda b: [x for x in b] if any(b) else []
    while bD and bD[-1] == 0: bD.pop()
    while bQ and bQ[-1] == 0: bQ.pop()
    assert bD == bQ, (u, bD, bQ)
for trial in range(150):
    nG, nM = random.choice([(3, 3), (4, 4), (5, 4), (4, 5)])
    R = [[random.random() < 0.55 for _ in range(nM)] for _ in range(nG)]
    D = dowker_complex(R)
    Q, chains = concept_poset_order_complex(R)
    bD = betti_of_complex(D) if D else []
    bQ = betti_of_complex(chains) if chains else []
    while bD and bD[-1] == 0: bD.pop()
    while bQ and bQ[-1] == 0: bQ.pop()
    assert bD == bQ, (R, bD, bQ)
print('    Betti numbers of D(R) and of the order complex of Q(R) agree '
      '(all levels of the test table + 150 random relations) OK')


# ---------------------------------------------------------------- Ex 5.13 / Prop 5.15: the interval mode
print('Ex. 5.13 / Prop. 5.15  the interval (SE-ISI) mode at a tie and near it')


def full_barcode(s, mode, cap):
    return [(float(lo), float(min(hi, cap))) for v in existence_sets(s, mode).values() for (lo, hi) in v]


x_sig = (frozenset({0}), frozenset({0}), frozenset({0, 1}))
A = [[Fr(10), Fr(-1)], [Fr(5), Fr(-5)]]
print('    Example 5.13 (tie at 5): profile of ({x},{m1},{m1,m2}) =', existence_sets(A, 'seisi')[x_sig],
      ' crit =', critical_values(A))
for eps in [Fr(1, 100), Fr(1, 10), Fr(2, 5)]:
    A2 = [[Fr(10), Fr(-1)], [Fr(5) + eps, Fr(-5) + eps]]
    bars2 = existence_sets(A2, 'seisi')[x_sig]
    dB = gudhi.bottleneck_distance(full_barcode(A, 'seisi', 20), full_barcode(A2, 'seisi', 20))
    dBs = gudhi.bottleneck_distance(full_barcode(A, 'sure', 20), full_barcode(A2, 'sure', 20))
    print(f'      eps={float(eps):<5} profile -> {[(float(a), float(b)) for a, b in bars2]}   '
          f'd_B(SE-ISI barcodes)={dB:.3f}   d_B(sure barcodes)={dBs:.3f}')

Bm = [[Fr(10), Fr(-1)], [Fr(9, 2), Fr(-11, 2)]]
print('    Proposition 5.15 with eta = 1 (no tie): crit =', [float(c) for c in critical_values(Bm)],
      ' profile of ({x},{m1},{m1,m2}) =', existence_sets(Bm, 'seisi')[x_sig])
mid = (frozenset({0, 1}), frozenset(), frozenset({0}))
print('      profile of ({x,y},{},{m1}) =', existence_sets(Bm, 'seisi').get(mid))
for eps in [Fr(2, 5), Fr(1, 2), Fr(3, 5)]:
    B2 = [[Fr(10), Fr(-1)], [Fr(9, 2) + eps, Fr(-11, 2) + eps]]
    bars2 = existence_sets(B2, 'seisi').get(x_sig)
    mid2 = existence_sets(B2, 'seisi').get(mid)
    dB = gudhi.bottleneck_distance(full_barcode(Bm, 'seisi', 20), full_barcode(B2, 'seisi', 20))
    print(f'      eps={float(eps):<4} profile of the x-signature -> {[(float(a), float(b)) for a, b in bars2]}; '
          f'({{x,y}},{{}},{{m1}}) -> {mid2};  d_B={dB:.3f}')

pass


# ---------------------------------------------------------------- Prop 5.11: certified extents (and intents)
def check_T8(s, ts):
    nG, nM = len(s), len(s[0])
    n = 0
    for t in ts:
        Rs, Rp = sure_rel(s, t), poss_rel(s, t)
        ext_all, int_all = None, None
        for R in completions(Rs, Rp):
            cs = concepts(R)
            E, I = {X for X, B in cs}, {B for X, B in cs}
            ext_all = E if ext_all is None else ext_all & E
            int_all = I if int_all is None else int_all & I
        assert ext_all == {X for X in powerset(nG) if t <= rho_ext(s, X)}
        assert int_all == {B for B in powerset(nM) if t <= rho_int(s, B)}
        n += 1
    for X in powerset(nG):
        for B in powerset(nM):
            assert rho_ext(s, X) >= rho(s, X, B) and rho_int(s, B) >= rho(s, X, B)
    return n


print('Prop. 5.11  extent (dually, intent) of every completion  <=>  t <= rho_ext (rho_int);  rho_ext, rho_int >= rho')
random.seed(3)
cnt = 0
for trial in range(150):
    nG, nM = random.choice([(2, 2), (3, 3), (3, 4), (4, 3)])
    s = rand_matrix(nG, nM, [-3, -2, -1, 0, 1, 2, 3])
    cnt += check_T8(s, [Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(5, 2), Fr(3), Fr(7, 2)])
print(f'    {cnt} (matrix, tolerance) instances, all completions enumerated OK')

# ---------------------------------------------------------------- ties in the sure mode
s1 = [[Fr(1)], [Fr(1)]]
e = Fr(1, 10)
s2 = [[Fr(1) + e], [Fr(1) - e]]
new = set(existence_sets(s2, 'sure')) - set(existence_sets(s1, 'sure'))
print('Ties in the sure mode: breaking a tie creates the sure concept', new, 'with profile', existence_sets(s2, 'sure')[next(iter(new))])

print('\n')


# ---------------------------------------------------------------- Prop 5.12 (general form) and Prop 5.11
print('Prop. 5.12 / Prop. 5.11  SE-ISI profile formula (general form) and certified-extent formula')
random.seed(11)
mats = 0
for trial in range(400):
    nG, nM = random.choice([(2, 2), (3, 2), (3, 3), (4, 3)])
    s = [[Fr(random.choice([-4, -3, -2, -1, 1, 2, 3, 4])) + Fr(random.randint(-3, 3), 10) for _ in range(nM)] for _ in range(nG)]
    bars = existence_sets(s, 'seisi')
    mats += any(len(v) > 1 for v in bars.values())
    pts = sorted({abs(x) for row in s for x in row})
    test = []
    for k in range(len(pts) + 1):
        lo = pts[k - 1] if k > 0 else 0
        hi = pts[k] if k < len(pts) else pts[-1] + 1
        test += [(lo + hi) / 2, hi]
    for (X, B, C), v in bars.items():
        ps = lambda m: psi(s, X, m)
        b1 = max((ps(m) for m in range(nM) if m not in B), default=-INF)
        d1 = min((ps(m) for m in B), default=INF)
        b2 = -min((ps(m) for m in C), default=INF)
        d2 = -max((ps(m) for m in range(nM) if m not in C), default=-INF)
        d3 = min((phi(s, B, g) for g in X), default=INF)
        b3 = -min((phi(s, C, g) for g in X), default=INF)
        beta, delta = max(b1, b2, b3, 0), min(d1, d2, d3)
        wins = [(-phi(s, C, g), phi(s, B, g)) for g in range(nG) if g not in X]
        for t in test:
            inE = beta < t <= delta and not any(a < t <= b for a, b in wins)
            assert inE == any(lo < t <= hi for lo, hi in v), (s, X, B, C, t)
print(f'    SE-ISI formula OK on 400 random tables; {mats} of them have a signature that disappears and reappears')
random.seed(3)
cnt = 0
for trial in range(150):
    nG, nM = random.choice([(2, 2), (3, 3), (3, 4), (4, 3)])
    s = rand_matrix(nG, nM, [-3, -2, -1, 0, 1, 2, 3])
    for t in [Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(5, 2), Fr(3), Fr(7, 2)]:
        Rs, Rp = sure_rel(s, t), poss_rel(s, t)
        if sum(Rp[g][m] and not Rs[g][m] for g in range(nG) for m in range(nM)) > 12:
            continue
        ext_all = None
        for R in completions(Rs, Rp):
            E = {X for X, B in concepts(R)}
            ext_all = E if ext_all is None else ext_all & E
        assert ext_all == {X for X in powerset(nG) if t <= rho_ext(s, X)}, (s, t)
        cnt += 1
print(f'    certified-extent formula OK on {cnt} (matrix, tolerance) pairs, all completions enumerated')
print('\nALL CHECKS PASSED')
