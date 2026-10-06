"""verify_examples.py -- recomputes every worked example of the paper, in exact arithmetic.
Statement numbers are those of the paper.
Run:  python3 verify_examples.py      (needs ptb.py; a few seconds)"""
from fractions import Fraction as Fr
from itertools import product, permutations
import math
import random
from ptb import *

INF = math.inf
F = frozenset
ok = lambda msg: print("  OK ", msg)


def match_dist(A, B):
    """Matching distance between two finite families of intervals (a, b] (Section 2)."""
    A, B = list(A), list(B)
    n = len(A) + len(B)
    Ad, Bd = A + [None] * len(B), B + [None] * len(A)
    best = INF
    for perm in permutations(range(n)):
        c = 0
        for i, j in enumerate(perm):
            a, b = Ad[i], Bd[j]
            if a is None and b is None:
                continue
            if a is None:
                c = max(c, (b[1] - b[0]) / 2)
            elif b is None:
                c = max(c, (a[1] - a[0]) / 2)
            else:
                c = max(c, abs(a[0] - b[0]), abs(a[1] - b[1]))
        best = min(best, c)
    return best


# ---------------------------------------------------------------- Examples 2.1 and 2.2: animals
print("Examples 2.1 and 2.2  animals")
G_ = ['duck', 'dog', 'fish']; M_ = ['flies', 'swims', 'walks']
R = [[1, 1, 1], [0, 0, 1], [0, 1, 0]]
C = concepts([[bool(x) for x in r] for r in R])
assert len(C) == 4
assert {(F({0, 1, 2}), F()), (F({0, 1}), F({2})), (F({0, 2}), F({1})), (F({0}), F({0, 1, 2}))} == C
assert (F({0}), F({0})) not in C                       # ({duck},{flies}) is not a concept
ok("four concepts; ({duck},{flies}) is not one")
unk = [(2, 0), (1, 1)]                                  # does the fish fly? does the dog swim?
survives = []
for bits in product([0, 1], repeat=2):
    Rc = [r[:] for r in R]
    for (g, m), b in zip(unk, bits):
        Rc[g][m] = b
    Cc = concepts([[bool(x) for x in r] for r in Rc])
    survives.append((F({0, 2}), F({1})) in Cc)
    if bits == (1, 0):                                  # fish flies, dog does not swim
        assert (F({0, 2}), F({0, 1})) in Cc
assert survives == [True, False, False, False]
ok("({duck,fish},{swims}) is a concept in exactly one of the four completions")

# ---------------------------------------------------------------- Section 2: matching distance
assert match_dist([(1, 10)], [(1, 9)]) == 1
assert match_dist([(Fr(1), Fr(10))], [(Fr(1), Fr(3, 2))]) == Fr(9, 2)
ok("matching distance examples: 1 and 9/2")

# ---------------------------------------------------------------- Example 3.3 and Table 1
print("Example 3.3 and Table 1  the running example")
S = [[Fr(4), Fr(-1, 2), Fr(3, 2)], [Fr(1), Fr(7, 2), Fr(-5, 2)], [Fr(-2), Fr(3), Fr(9, 2)]]
raw = {0: ([80, 90, 70], [64, 86, 78], Fr(4)),
       1: ([15, 21, 24], [16, 14, 18], Fr(2)),
       2: ([Fr('37.2'), Fr('36.0'), Fr('38.3')], [Fr('36.6'), Fr('37.0'), Fr('36.5')], Fr('0.4'))}
for m, (v, c, d) in raw.items():
    assert [(Fr(v[g]) - c[g]) / d for g in range(3)] == [S[g][m] for g in range(3)]
assert sorted(abs(x) for r in S for x in r) == [Fr(k, 2) for k in range(1, 10)]
ok("readings, baselines and half-widths give the budgets of Table 1; |s| = 1/2, 1, ..., 9/2")
row_p1 = [incomplete_context(S, t)[0][1] for t in [Fr(k, 8) for k in range(1, 60)]]
assert all(x == '-' for x, t in zip(row_p1, [Fr(k, 8) for k in range(1, 60)]) if t <= Fr(1, 2))
assert all(x == '?' for x, t in zip(row_p1, [Fr(k, 8) for k in range(1, 60)]) if t > Fr(1, 2))
assert all(incomplete_context(S, t)[0][0] == '+' for t in (Fr(1, 2), Fr(4))) and incomplete_context(S, Fr(41, 10))[0][0] == '?'
ok("respiratory cell of p1: decided negatively for t <= 1/2, doubtful after, never positive; heart-rate cell of p1 positive up to t = 4")

# ---------------------------------------------------------------- Remark 3.7: unrestricted bands produce any incomplete context
random.seed(1)
for _ in range(50):
    nG, nM = random.randint(1, 4), random.randint(1, 4)
    I = [[random.choice('+-?') for _ in range(nM)] for _ in range(nG)]
    v = [[Fr(random.randint(-9, 9)) for _ in range(nM)] for _ in range(nG)]
    s_ = [[{'+': Fr(1), '-': Fr(-1), '?': Fr(0)}[I[g][m]] for m in range(nM)] for g in range(nG)]
    c = [[v[g][m] - s_[g][m] for m in range(nM)] for g in range(nG)]
    sTheta = [[(v[g][m] - c[g][m]) / 1 for m in range(nM)] for g in range(nG)]
    for t in (Fr(1, 3), Fr(1, 2), Fr(1)):
        assert incomplete_context(sTheta, t) == I
ok("Remark 3.7: delta = 1 and c = v - s reproduce any prescribed incomplete context for t in (0,1]")

# ---------------------------------------------------------------- Example 3.9: the expressiveness criterion
v, c = [10, 12, 9], [8, 11, 7]
s1 = [Fr(v[i] - c[i], 1) for i in range(3)]
assert s1 == [2, 1, 2]
assert s1[0] > max(s1[1], 0) and v[0] <= v[1]
for t in (Fr(11, 10), Fr(3, 2), Fr(2)):
    assert {g for g in range(3) if s1[g] >= t} == {0, 2}
assert all({g for g in range(3) if v[g] >= th} != {0, 2} for th in range(0, 15))
ok("Example 3.9: {x,z} is the sure extent for t in (1,2] and no global threshold produces it")

# ---------------------------------------------------------------- Proposition 4.1: the critical grid of the running example
crit = critical_values(S)
assert len(crit) == 9
ctxs = [tuple(map(tuple, incomplete_context(S, t))) for t in [Fr(1, 4)] + [c_ + Fr(1, 4) for c_ in crit]]
assert len(set(ctxs)) == 10
assert all(x != '?' for r in incomplete_context(S, Fr(1, 2)) for x in r)
assert [r.count('?') for r in incomplete_context(S, Fr(3, 4))] == [1, 0, 0]
assert incomplete_context(S, Fr(5, 4))[1][0] == '?'
ok("nine critical values, ten distinct contexts, K(t) = K(0+) for t <= 1/2")

# ---------------------------------------------------------------- Example 5.3 and the discussion after Theorem 5.7
print("Example 5.3  the 3 x 2 table")
s = [[Fr(3), Fr(1)], [Fr(2), Fr(-1)], [Fr(-2), Fr(4)]]
X, B = F({0, 1}), F({0})
assert inn(s, X, B) == 2 and out(s, X, B) == -1 and rho(s, X, B) == 1
assert existence_sets(s, 'sure')[(X, B)] == [(0, Fr(2))]
assert existence_sets(s, 'possible')[(X, B)] == [(0, Fr(1))]
assert (X, B) in concepts(level_rel(s, 1)) and (X, B) not in concepts(level_rel(s, Fr(5, 2)))
assert (X, B) not in concepts(level_rel(s, -Fr(3, 2), strict=True))
t = Fr(3, 2)
assert incomplete_context(s, t)[0][1] == '?' and incomplete_context(s, t)[1][1] == '?'
Rc = level_rel(s, t); Rc[1][1] = True
assert (X, B) in concepts(Rc)                        # (g2,b) alone does not suffice: g1 still lacks b
Rc[0][1] = True
assert (X, B) not in concepts(Rc) and (X, F({0, 1})) in concepts(Rc)
X3, B3 = F({2}), F({1})
assert inn(s, X3, B3) == 4 and out(s, X3, B3) == 1 and rho(s, X3, B3) == -1
assert existence_sets(s, 'sure')[(X3, B3)] == [(Fr(1), Fr(4))]
assert (X3, B3) not in existence_sets(s, 'possible')
ok("({g1,g2},{a}): sure (0,2], possible (0,1], rho = 1; at t = 3/2 both doubtful cells of b must be +")
ok("({g3},{b}): sure (1,4], empty possible profile, rho = -1, never certified")

# ---------------------------------------------------------------- Proposition 5.9(ii): the converse fails
s2 = [[Fr(2)], [Fr(0)]]
assert (F({0}), F({0}), F({0})) in seisi_concepts(level_rel(s2, 1), level_rel(s2, -1, True))
assert rho(s2, F({0}), F({0})) < 1
ok("Proposition 5.9(ii): ({g},{m},{m}) is SE-ISI at t = 1 but ({g},{m}) is not certified")

# ---------------------------------------------------------------- Example 5.13 and Proposition 5.15
print("Example 5.13 and Proposition 5.15  ties in the interval mode")
sig = (F({0}), F({0}), F({0, 1}))
A = [[Fr(10), Fr(-1)], [Fr(5), Fr(-5)]]
assert existence_sets(A, 'seisi')[sig] == [(Fr(1), Fr(10))]
for eps in (Fr(1, 10), Fr(1, 100), Fr(1, 1000)):
    Ae = [A[0][:], [x + eps for x in A[1]]]
    E = existence_sets(Ae, 'seisi')[sig]
    assert E == [(Fr(1), 5 - eps), (5 + eps, Fr(10))]
    assert match_dist([(Fr(1), Fr(10))], E) > 4
ok("tie at 5: the profile (1,10] splits for every eps > 0, at matching distance > 4")
s_1 = [[Fr(10), Fr(-1)], [Fr(9, 2), Fr(-11, 2)]]
E0 = existence_sets(s_1, 'seisi')[sig]
E1 = existence_sets([s_1[0][:], [x + Fr(3, 5) for x in s_1[1]]], 'seisi')[sig]
assert E0 == [(Fr(1), Fr(10))] and match_dist(E0, E1) == Fr(41, 10)
ok("eta = 1, eps = 0.6: the profile splits and the matching distance is 4.1")

# ---------------------------------------------------------------- Section 6.1: the crown, Table 2 and Figure 1
print("Section 6.1  the crown")
crisp = concepts(level_rel(S, 0, strict=True))
table2 = {
    (F({0, 1, 2}), F()):        (INF, Fr(-1, 2), (0, INF), INF),
    (F({0, 1}), F({0})):        (Fr(1), Fr(-1, 2), (0, Fr(1)), Fr(1)),
    (F({0, 2}), F({2})):        (Fr(3, 2), Fr(-1, 2), (0, Fr(3, 2)), Fr(3, 2)),
    (F({1, 2}), F({1})):        (Fr(3), Fr(-1, 2), (0, Fr(3)), Fr(1, 2)),
    (F({0}), F({0, 2})):        (Fr(3, 2), Fr(-1, 2), (0, Fr(3, 2)), Fr(3, 2)),
    (F({1}), F({0, 1})):        (Fr(1), Fr(-1, 2), (0, Fr(1)), Fr(1, 2)),
    (F({2}), F({1, 2})):        (Fr(3), Fr(-1, 2), (0, Fr(3)), Fr(1, 2)),
    (F(), F({0, 1, 2})):        (INF, Fr(-1, 2), (0, INF), Fr(1, 2)),
}
assert set(table2) == crisp
sure, poss = existence_sets(S, 'sure'), existence_sets(S, 'possible')
for (X, B), (i, o, prof, rext) in table2.items():
    assert inn(S, X, B) == i and out(S, X, B) == o and rho(S, X, B) == Fr(1, 2)
    assert sure[(X, B)] == [prof] and poss[(X, B)] == [(0, Fr(1, 2))]
    assert rho_ext(S, X) == rext
ok("Table 2: inner and outer levels, sure and possible profiles, rho = 1/2 and rho_ext for the eight crisp concepts")
assert {g for g in range(3) if S[g][0] >= 1} == {0, 1} and {g for g in range(3) if S[g][0] >= Fr(3, 2)} == {0}
assert S[0][0] > max(S[1][0], 0) and raw[0][0][0] <= raw[0][0][1]
assert S[1][1] > max(S[2][1], 0) and raw[1][0][1] <= raw[1][0][2]
assert not any(S[g][2] > max(S[h][2], 0) and raw[2][0][g] <= raw[2][0][h] for g in range(3) for h in range(3))
ok("Step 1: columns a1 and a2 satisfy the criterion of Proposition 3.8, column a3 does not")
for t in (Fr(3, 4), Fr(1)):
    assert sum(r.count('?') for r in incomplete_context(S, t)) == 1   # only (p1,a2) is doubtful
    Rc = level_rel(S, t); Rc[0][1] = True                           # the completion in which p1 has a2
    assert all((X, B) not in concepts(Rc) for (X, B) in crisp)
ok("Step 4: for t in (1/2,1] the completion in which p1 has a2 destroys each of the eight concepts")
Rc = level_rel(S, Fr(3, 4)); Rc[0][1] = True
assert (F({0, 1}), F({0, 1})) in concepts(Rc)
ok("Step 5: rho_ext({p1,p2}) = 1 > rho = 1/2; in the completion with (p1,a2) = + the intent grows to {a1,a2}")
assert sure[(F({0}), F({0}))] == [(Fr(3, 2), Fr(4))]
assert sure[(F({1}), F({1}))] == [(Fr(3), Fr(7, 2))]
assert sure[(F({2}), F({2}))] == [(Fr(3), Fr(9, 2))]
assert inn(S, F({0, 1, 2}), F({1})) == Fr(-1, 2) and out(S, F({0, 1, 2}), F({1})) == -2
assert poss[(F({0, 1, 2}), F({1}))] == [(Fr(1, 2), Fr(2))]
se = existence_sets(S, 'seisi')
assert len(se) == 24 and all(len(v_) == 1 for v_ in se.values())
ok("Step 6: late sure concepts, the possible concept (G,{a2}) on (1/2,2], 24 SE-ISI signatures each on one interval")
Fv = {sig_: dowker_F(S, sig_) for sig_ in [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]}
assert Fv == {(0,): 4, (1,): Fr(7, 2), (2,): Fr(9, 2), (0, 1): 1, (0, 2): Fr(3, 2), (1, 2): 3, (0, 1, 2): Fr(-1, 2)}
ok("Step 7: the values of F on vertices, segments and the triangle")

# ---------------------------------------------------------------- Section 6.2: a tie
print("Section 6.2  a tie")
T = [r[:] for r in S]; T[1][0] = Fr(96 - 86, 4)
assert T[1][0] == Fr(5, 2) == -T[1][2]
sigp = (F({0}), F({0}), F({0, 1, 2}))
assert existence_sets(T, 'seisi')[sigp] == [(Fr(3, 2), Fr(4))]
eps = Fr(1, 100)
Te = [r[:] for r in T]; Te[1] = [x + eps for x in Te[1]]
Ee = existence_sets(Te, 'seisi')[sigp]
assert Ee == [(Fr(3, 2), Fr(5, 2) - eps), (Fr(5, 2) + eps, Fr(4))]
assert match_dist([(Fr(3, 2), Fr(4))], Ee) == 1 + eps
ok("the SE-ISI profile (3/2,4] of ({p1},{a1},M) splits at 5/2, a jump of 1 + eps in matching distance")

# ---------------------------------------------------------------- Section 7: Galois closures are not additive
print("Section 7  closures of the modal-style example")
Rm = [[True, False], [False, True], [True, True], [False, False]]  # objects a, b, c, d; attributes m1, m2
def clos(X):
    return down(Rm, up(Rm, X, 2), 4)
assert clos(F({0})) == F({0, 2}) and clos(F({1})) == F({1, 2})
assert up(Rm, F({0, 1}), 2) == F() and clos(F({0, 1})) == F({0, 1, 2, 3})
assert clos(F()) == F({2})
ok("{a}'' = {a,c}, {b}'' = {b,c}, {a,b}'' = G and the closure of the empty set is {c}")

print("\nALL EXAMPLES VERIFIED")
