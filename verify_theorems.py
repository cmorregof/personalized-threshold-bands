"""verify_theorems.py -- exhaustive and randomized checks, in exact rational arithmetic,
of the statements of Sections 3 to 7 of the paper (statement numbers as in the paper):
Prop. 3.8, Thm. 5.2, Thm. 5.7, Props. 5.8, 5.9, 5.11, 5.12, Thm. 5.14, Prop. 5.15, Thm. 5.19,
Props. 7.1, 7.2, Cor. 7.3 and the running example of Section 6. The counts printed for
Thm. 5.2 and Thm. 5.19 are those quoted in Remark 5.20.
Run:  python3 verify_theorems.py      (needs ptb.py and gudhi; about a minute)"""
import random, math
from fractions import Fraction as Fr
from itertools import product, permutations
from ptb import *
import gudhi
random.seed(2026)
INF = math.inf
ok = lambda msg: print("  OK ", msg)

# ---------------------------------------------------------------- helpers
def rank_map(s, t, k):
    Rs = level_rel(s, t); Rp = level_rel(s, -t, strict=True)
    Ds = set(dowker_complex(Rs)); Dp = dowker_complex(Rp)
    if not Dp: return 0
    st = gudhi.SimplexTree()
    for sig in Dp: st.insert(list(sig), filtration=0.0 if sig in Ds else 1.0)
    st.make_filtration_non_decreasing()
    st.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
    pb = st.persistent_betti_numbers(0.0, 1.0)
    return pb[k] if k < len(pb) else 0

def quadrant(s, t, k):
    return sum(1 for d, ub, ud in dowker_persistence(s) if d == k and ub >= t and ud <= -t)

def betti(R, k):
    cx = dowker_complex(R)
    if not cx: return 0
    b = betti_of_complex(cx); return b[k] if k < len(b) else 0

def bneck(A, B):
    A = list(A); B = list(B); n = len(A) + len(B)
    Ad = A + [None]*len(B); Bd = B + [None]*len(A); best = INF
    for perm in permutations(range(n)):
        c = 0
        for i, j in enumerate(perm):
            a, b = Ad[i], Bd[j]
            if a is None and b is None: continue
            if a is None: c = max(c, (b[1]-b[0])/2)
            elif b is None: c = max(c, (a[1]-a[0])/2)
            else: c = max(c, abs(a[0]-b[0]), abs(a[1]-b[1]))
        best = min(best, c)
    return best

def completions(s, t):
    nG, nM = len(s), len(s[0])
    unk = [(g, m) for g in range(nG) for m in range(nM) if -t < s[g][m] < t]
    base = level_rel(s, t)
    for bits in product([0, 1], repeat=len(unk)):
        R = [row[:] for row in base]
        for (g, m), b in zip(unk, bits): R[g][m] = bool(b)
        yield R

vals = [Fr(x, 2) for x in range(-8, 9)]
def rand_s(nG, nM, extra=()):
    return [[random.choice(vals + list(extra)) for _ in range(nM)] for _ in range(nG)]

# ---------------------------------------------------------------- Prop 3.8 (exact expressiveness criterion)
print("Prop. 3.8  exact criterion for the expressiveness gap")
cnt = 0
for _ in range(3000):
    n = random.randint(2, 5)
    v = [random.randint(0, 6) for _ in range(n)]
    s = [random.choice(vals) for _ in range(n)]
    scal = {frozenset(g for g in range(n) if v[g] >= th) for th in range(-1, 9)}
    ts = sorted({x for x in s if x > 0} | {Fr(1, 4)})
    pers = {frozenset(g for g in range(n) if s[g] >= t) for t in ts}
    gap = not pers <= scal
    crit = any(s[g] > max(s[h], 0) and v[g] <= v[h] for g in range(n) for h in range(n))
    assert gap == crit, (v, s); cnt += 1
ok(f"{cnt} random columns: gap <=> exists g,h with s(g) > max(s(h),0) and v(g) <= v(h)")

# ---------------------------------------------------------------- Thm 5.2 with extended-real tables
print("Thm. 5.2  closed-form profiles on tables with +-inf entries")
cnt = 0
for _ in range(300):
    nG, nM = random.choice([(2, 2), (3, 3), (3, 4), (4, 3)])
    s = rand_s(nG, nM, extra=(INF, -INF))
    ts = sorted({abs(x) for r in s for x in r if 0 < abs(x) < INF} | {Fr(1, 4)})
    ts = sorted(set(ts) | {t + Fr(1, 4) for t in ts})
    for t in ts:
        Cs = concepts(level_rel(s, t)); Cp = concepts(level_rel(s, -t, strict=True))
        for X in powerset(nG):
            for B in powerset(nM):
                i, o = death(s, X, B), birth(s, X, B)
                assert ((X, B) in Cs) == (max(o, 0) < t <= i)
                assert ((X, B) in Cp) == (max(-i, 0) < t <= -o)
                cnt += 1
ok(f"{cnt} (table, t, signature) cases, entries include +-inf")

# ---------------------------------------------------------------- Thm 5.7 + Prop 5.8 (three-way reading) + Prop 5.9 (Li-Mei-Lv)
print("Thm. 5.7 / Prop. 5.8 / Prop. 5.9  certification, three-way reading, Li-Mei-Lv bridge")
cnt = n_exact_not_cert = 0
for _ in range(200):
    nG, nM = random.choice([(2, 3), (3, 3), (3, 2)])
    s = rand_s(nG, nM)
    for t in [Fr(1, 4), Fr(1, 2), Fr(1), Fr(3, 2), Fr(2)]:
        comps = [concepts(R) for R in completions(s, t)]
        Rs = level_rel(s, t); Rneg = [[x <= -t for x in r] for r in s]
        se = seisi_concepts(Rs, level_rel(s, -t, strict=True))
        for X in powerset(nG):
            for B in powerset(nM):
                cert = all((X, B) in C for C in comps)
                assert cert == (t <= rho(s, X, B))
                tw = (all(Rs[g][m] for g in X for m in B)
                      and all(any(Rneg[g][m] for m in B) for g in range(nG) if g not in X)
                      and all(any(Rneg[g][m] for g in X) for m in range(nM) if m not in B))
                assert cert == tw
                if cert: assert (X, B, B) in se
                elif (X, B, B) in se: n_exact_not_cert += 1
                cnt += 1
ok(f"{cnt} cases: certified <=> t<=rho <=> three-way condition; certified => exact SE-ISI")
ok(f"exact SE-ISI concepts that are not certified: {n_exact_not_cert} (converse fails, as stated)")
s = [[Fr(2)], [Fr(0)]]; t = Fr(1)
assert (frozenset({0}), frozenset({0}), frozenset({0})) in seisi_concepts(level_rel(s, t), level_rel(s, -t, True))
assert rho(s, frozenset({0}), frozenset({0})) < t
ok("counterexample of the converse (G={g,h}, M={m}, s=(2,0), t=1)")

# ---------------------------------------------------------------- Prop 5.15 (sharp radius)
print("Prop. 5.15  sharpness family for the interval mode")
sig = (frozenset({0}), frozenset({0}), frozenset({0, 1}))
for eta in [Fr(1), Fr(1, 2), Fr(1, 10)]:
    s = [[Fr(10), Fr(-1)], [5 - eta/2, -5 - eta/2]]
    E = existence_sets(s, 'seisi')[sig]; assert E == [(Fr(1), Fr(10))]
    for eps in [eta/2 + Fr(1, 1000), eta, Fr(1), Fr(39, 10)]:
        sp = [s[0][:], [x + eps for x in s[1]]]
        Ep = existence_sets(sp, 'seisi')[sig]
        assert Ep == [(Fr(1), 5 + eta/2 - eps), (5 - eta/2 + eps, Fr(10))]
        assert bneck([(float(a), float(b)) for a, b in E], [(float(a), float(b)) for a, b in Ep]) > 4
ok("eta in {1,1/2,1/10}, eps in (eta/2,4): the profile splits and d_sig > 4")

# ---------------------------------------------------------------- Thm 5.19 (certified Betti numbers)
print("Thm. 5.19  certified Betti numbers")
n_rank = n_lb = n_stab = 0
for _ in range(150):
    nG, nM = random.choice([(3, 3), (3, 4), (4, 3), (4, 4)])
    s = rand_s(nG, nM)
    crit = sorted({abs(x) for r in s for x in r if x != 0})
    ts = sorted(set(crit) | {c + Fr(1, 4) for c in crit} | {Fr(1, 4)})
    for t in ts:
        for k in (0, 1, 2):
            assert rank_map(s, t, k) == quadrant(s, t, k); n_rank += 1
        unk = sum(1 for r in s for x in r if -t < x < t)
        if unk <= 8:
            for R in completions(s, t):
                for k in (0, 1): assert betti(R, k) >= rank_map(s, t, k)
                n_lb += 1
    eps = Fr(1, 4)
    sp = [[x + random.choice([-eps, 0, eps]) for x in r] for r in s]
    for t in ts:
        for k in (0, 1):
            assert rank_map(sp, t, k) >= rank_map(s, t + eps, k)
            assert rank_map(s, t, k) >= rank_map(s, t + Fr(1, 2), k)
            n_stab += 1
ok(f"rank = quadrant count on {n_rank} cases")
ok(f"beta_k(D(R)) >= certified Betti for all {n_lb} completions enumerated")
ok(f"monotone in t and beta'_k(t) >= beta_k(t+eps) on {n_stab} cases")

# ---------------------------------------------------------------- Section 6: the crown
print("Section 6  the running example (the crown)")
S = [[Fr(4), Fr(-1, 2), Fr(3, 2)], [Fr(1), Fr(7, 2), Fr(-5, 2)], [Fr(-2), Fr(3), Fr(9, 2)]]
raw = {0: ([80, 90, 70], [64, 86, 78], 4), 1: ([15, 21, 24], [16, 14, 18], 2),
       2: ([Fr('37.2'), Fr('36.0'), Fr('38.3')], [Fr('36.6'), Fr('37.0'), Fr('36.5')], Fr('0.4'))}
for m, (v, c, d) in raw.items():
    assert [Fr(v[g] - c[g]) / d for g in range(3)] == [S[g][m] for g in range(3)]
ok("raw readings, baselines and half-widths reproduce the certainty table")
absv = sorted(abs(x) for r in S for x in r)
assert absv == [Fr(k, 2) for k in range(1, 10)]
ok("|s| = {1/2, 1, ..., 9/2}: generic with eta = 1/2")
crisp = concepts(level_rel(S, 0, strict=True)); assert len(crisp) == 8
inn = {}
for X, B in crisp:
    assert birth(S, X, B) == Fr(-1, 2) and rho(S, X, B) == Fr(1, 2); inn[(X, B)] = death(S, X, B)
F = frozenset
assert inn[(F({0, 1}), F({0}))] == 1 and inn[(F({0, 2}), F({2}))] == Fr(3, 2) and inn[(F({1, 2}), F({1}))] == 3
assert inn[(F({0}), F({0, 2}))] == Fr(3, 2) and inn[(F({1}), F({0, 1}))] == 1 and inn[(F({2}), F({1, 2}))] == 3
ok("8 crisp concepts; out = -1/2 and rho = 1/2 for all; inner levels 1, 3/2, 3 (pairs) and 3/2, 1, 3 (singletons)")
assert rho_ext(S, F({0, 1})) == 1 and rho_ext(S, F({0, 2})) == Fr(3, 2) and rho_ext(S, F({1, 2})) == Fr(1, 2)
ok("rho_ext: {p1,p2} -> 1, {p1,p3} -> 3/2, {p2,p3} -> 1/2")
ri = {B: rho_int(S, B) for B in powerset(3)}
assert ri == {F(): Fr(1, 2), F({0}): Fr(1, 2), F({2}): Fr(1, 2), F({1}): 2, F({0, 2}): Fr(1, 2),
              F({0, 1}): 1, F({1, 2}): 2, F({0, 1, 2}): INF}
ok("Table 2, last column: rho_int = 1/2, 1/2, 1/2, 2, 1/2, 1, 2, +inf")
def cert_by_enum(t, kind):
    fams = []
    for R in completions(S, t):
        if kind == 'int':
            fams.append({up(R, down(R, B, 3), 3) for B in powerset(3)})
        else:
            fams.append({down(R, up(R, X, 3), 3) for X in powerset(3)})
    return set.intersection(*fams)
assert cert_by_enum(Fr(1), 'int') == {F({1}), F({0, 1}), F({1, 2}), F({0, 1, 2})}
assert cert_by_enum(Fr(2), 'int') == {F({1}), F({1, 2}), F({0, 1, 2})}
assert cert_by_enum(Fr(3), 'int') == {F({0, 1, 2})}
assert cert_by_enum(Fr(1), 'ext') == {F({0, 1, 2}), F({0, 1}), F({0, 2}), F({0})}
assert cert_by_enum(Fr(3, 2), 'ext') == {F({0, 1, 2}), F({0, 2}), F({0})}
assert cert_by_enum(Fr(2), 'ext') == {F({0, 1, 2})}
assert not set.intersection(*[concepts(R) for R in completions(S, Fr(1))])
ok("Step 6: certified intents and extents at t = 1, 3/2, 2, 3 by enumeration; no certified concept at t = 1")
sure = existence_sets(S, 'sure')
assert sure[(F({0}), F({0}))] == [(Fr(3, 2), Fr(4))] and sure[(F({1}), F({1}))] == [(Fr(3), Fr(7, 2))] and sure[(F({2}), F({2}))] == [(Fr(3), Fr(9, 2))]
poss = existence_sets(S, 'possible')
assert poss[(F({0, 1, 2}), F({1}))] == [(Fr(1, 2), Fr(2))]
se = existence_sets(S, 'seisi'); assert len(se) == 24 and all(len(v) == 1 for v in se.values())
ok("late sure concepts, possible concept (G,{a2}) on (1/2,2], 24 SE-ISI signatures each on one interval")
bars = sorted(dowker_persistence(S))
assert bars == sorted([(0, 4.5, -INF), (0, 4.0, 1.5), (0, 3.5, 3.0), (1, 1.0, -0.5)])
assert [rank_map(S, t, 1) for t in (Fr(1, 2), Fr(51, 100))] == [1, 0]
assert all(rank_map(S, t, 0) == 1 for t in (Fr(1, 4), Fr(9, 2))) and rank_map(S, Fr(5), 0) == 0
ok("Dowker bars; certified beta_1 = 1 exactly for t <= 1/2; certified beta_0 = 1 for t <= 9/2")
# column-wise personalization
for m, (g, h) in [(0, (0, 1)), (1, (1, 2))]:
    v = raw[m][0]
    assert S[g][m] > max(S[h][m], 0) and v[g] <= v[h]
ok("columns a1 and a2 satisfy the criterion of Prop. 3.8; a3 does not")
assert not any(S[g][2] > max(S[h][2], 0) and raw[2][0][g] <= raw[2][0][h] for g in range(3) for h in range(3))

# ---------------------------------------------------------------- Section 6.2: the tie variant
print("Section 6.2  tie variant (p2 heart rate 96 instead of 90)")
T = [r[:] for r in S]; T[1][0] = Fr(96 - 86, 4); assert T[1][0] == Fr(5, 2) == -T[1][2]
sigp = (F({0}), F({0}), F({0, 1, 2}))
assert existence_sets(T, 'seisi')[sigp] == [(Fr(3, 2), Fr(4))]
for eps in [Fr(1, 100), Fr(1, 1000)]:
    Tp = [r[:] for r in T]; Tp[1] = [x + eps for x in Tp[1]]
    Ep = existence_sets(Tp, 'seisi')[sigp]
    assert Ep == [(Fr(3, 2), Fr(5, 2) - eps), (Fr(5, 2) + eps, Fr(4))]
    d = bneck([(1.5, 4.0)], [(float(a), float(b)) for a, b in Ep]); assert abs(d - 1) < 0.02
    for mode in ('sure', 'possible'):
        A, B = existence_sets(T, mode), existence_sets(Tp, mode)
        for k in set(A) | set(B):
            a = A.get(k, []); b = B.get(k, [])
            fa = [(float(x), float(y)) for x, y in a]; fb = [(float(x), float(y)) for x, y in b]
            if fa and fb and fa[0][1] == INF: assert fb[0][1] == INF; continue
            assert bneck(fa, fb) <= float(eps) + 1e-12
ok("the SE-ISI profile of ({p1},{a1},M) splits at 5/2 (jump ~1) while sure/possible profiles move <= eps")
print()

# ---------------------------------------------------------------- Prop 5.12 (profiles of SE-ISI signatures) and Thm 5.14 on random generic tables
print("Prop. 5.12 / Thm. 5.14  SE-ISI profile formula; conditional stability on generic tables")
def seisi_formula(s, X, B, C, t):
    nG, nM = len(s), len(s[0])
    tm = max([0] + [psi(s, X, m) for m in range(nM) if m not in B] + [-death(s, X, C)])
    tp = min([death(s, X, B)] + [-psi(s, X, m) for m in range(nM) if m not in C])
    if not (tm < t <= tp): return False
    return all(not (-phi(s, C, g) < t <= phi(s, B, g)) for g in range(nG) if g not in X)
n1 = 0
for _ in range(200):
    nG, nM = random.choice([(2, 2), (3, 3), (3, 2), (2, 3)])
    s = rand_s(nG, nM)
    crit = sorted({abs(x) for r in s for x in r if x != 0})
    for t in sorted(set(crit) | {c + Fr(1, 4) for c in crit} | {Fr(1, 4)}):
        se = seisi_concepts(level_rel(s, t), level_rel(s, -t, strict=True))
        for X in powerset(nG):
            for B in powerset(nM):
                for C in powerset(nM):
                    assert ((X, B, C) in se) == seisi_formula(s, X, B, C, t); n1 += 1
ok(f"tau^-, tau^+ and windows Omega_g on {n1} (table, t, triple) cases")
n2 = 0
for _ in range(150):
    nG, nM = random.choice([(2, 2), (3, 3), (2, 3)])
    mags = random.sample(range(1, 40), nG * nM)
    s = [[Fr(mags[g * nM + m], 2) * random.choice([-1, 1]) for m in range(nM)] for g in range(nG)]
    absv = sorted(abs(x) for r in s for x in r)
    eta = min([absv[0]] + [b - a for a, b in zip(absv, absv[1:])])
    eps = eta / 2 - Fr(1, 100)
    if eps <= 0: continue
    sp = [[x + random.choice([-eps, eps, eps / 3]) for x in r] for r in s]
    E, Ep = existence_sets(s, 'seisi'), existence_sets(sp, 'seisi')
    assert set(E) == set(Ep)
    for k in E:
        assert len(E[k]) == len(Ep[k])
        for (a, b), (a2, b2) in zip(E[k], Ep[k]):
            assert abs(a - a2) <= eps and (b == b2 == INF or abs(b - b2) <= eps)
    n2 += 1
ok(f"Theorem 5.14 on {n2} random generic tables perturbed by eps < eta/2")

# ---------------------------------------------------------------- GUDHI input: float conversion is injective on the values used
print("GUDHI input  exact-to-float conversion is injective")
def injective(s):
    vals = {dowker_F(s, sig) for r in range(1, len(s) + 1) for sig in combinations(range(len(s)), r)}
    return len({float(x) for x in vals}) == len(vals)
from itertools import combinations
assert injective(S) and injective(T)
assert all(injective(rand_s(4, 4)) for _ in range(200))
ok("on the crown, its tie variant and 200 random tables")
# ---------------------------------------------------------------- Prop 5.11 (ii)-(iv): certified intents and closure systems
print("Prop. 5.11  certified extents and intents: formulas, B^(sure poss) = B, closure systems, monotonicity")
rg = random.Random(511)
vals2 = [Fr(x, 2) for x in range(-6, 7)]
def rnd_table(nG, nM):
    out = []
    for _ in range(nG):
        row = []
        for _ in range(nM):
            r = rg.random()
            row.append(Fr(0) if r < 0.15 else (rg.choice([INF, -INF]) if r < 0.2 else rg.choice(vals2)))
        out.append(row)
    return out
n_int = n_cl = n_ineq = 0
for _ in range(500):
    nG, nM = rg.randint(1, 4), rg.randint(1, 4)
    s_ = rnd_table(nG, nM)
    ts = sorted({abs(x) for r in s_ for x in r if x not in (0, INF, -INF)} | {Fr(1, 3)})
    ts = ts + [ts[-1] + 1]
    for t in ts:
        comps = list(completions(s_, t))
        if len(comps) > 4096: continue
        Rs, Rp = level_rel(s_, t), level_rel(s_, -t, strict=True)
        ci = set.intersection(*[{up(R, down(R, B, nG), nM) for B in powerset(nM)} for R in comps])
        ce = set.intersection(*[{down(R, up(R, X, nM), nG) for X in powerset(nG)} for R in comps])
        for B in powerset(nM):
            assert (B in ci) == (up(Rp, down(Rs, B, nG), nM) == B) == (t <= rho_int(s_, B)); n_int += 1
        for X in powerset(nG):
            assert (X in ce) == (t <= rho_ext(s_, X))
        assert frozenset(range(nM)) in ci and frozenset(range(nG)) in ce
        assert all(A & C in ci for A in ci for C in ci) and all(A & C in ce for A in ce for C in ce); n_cl += 1
    for t1, t2 in zip(ts, ts[1:]):
        assert {B for B in powerset(nM) if t2 <= rho_int(s_, B)} <= {B for B in powerset(nM) if t1 <= rho_int(s_, B)}
        assert {X for X in powerset(nG) if t2 <= rho_ext(s_, X)} <= {X for X in powerset(nG) if t1 <= rho_ext(s_, X)}
    for X in powerset(nG):
        for B in powerset(nM):
            r = rho(s_, X, B)
            assert rho_ext(s_, X) >= r and rho_int(s_, B) >= r; n_ineq += 1
ok(f"(ii) on {n_int} (table, t, B) cases; (iv) closure systems on {n_cl} (table, t) cases; monotone in t")
ok(f"(iii) rho_ext(X) >= rho(X,B) and rho_int(B) >= rho(X,B) for all {n_ineq} signatures")

# ---------------------------------------------------------------- Section 7: Props 7.1, 7.2 and Cor 7.3
print("Props. 7.1, 7.2 and Cor. 7.3  missing values as zero budgets")
n71 = n72 = 0
for _ in range(800):
    nG, nM = rg.randint(1, 4), rg.randint(1, 4)
    sstar = [[rg.choice([v for v in vals2 if v != 0] + [INF, -INF]) for _ in range(nM)] for _ in range(nG)]
    t = rg.choice([Fr(1, 4), Fr(1, 2), Fr(1), Fr(3, 2), Fr(5, 2)])
    shat, miss = [], []
    for g in range(nG):
        row, mrow = [], []
        for m in range(nM):
            if rg.random() < 0.3:
                row.append(Fr(0)); mrow.append(True)
            else:
                x = sstar[g][m]
                row.append(x if x in (INF, -INF) else x + Fr(rg.randint(-99, 99), 100) * t); mrow.append(False)
        shat.append(row); miss.append(mrow)
    Rstar = [[x > 0 for x in r] for r in sstar]
    Rs, Rp = level_rel(shat, t), level_rel(shat, -t, strict=True)
    assert all(Rs[g][m] <= Rstar[g][m] <= Rp[g][m] for g in range(nG) for m in range(nM))   # Prop 7.1
    comps = list(completions(shat, t))
    if len(comps) > 4096: continue
    for B in powerset(nM):                                                                    # Prop 7.2
        supp = [len(down(R, B, nG)) for R in comps]
        assert min(supp) == len(down(Rs, B, nG)) and max(supp) == len(down(Rp, B, nG))
        assert len(down(Rs, B, nG)) <= len(down(Rstar, B, nG)) <= len(down(Rp, B, nG))
        for g in range(nG):
            has = [all(R[g][m] for m in B) for R in comps]
            yes, no = all(Rs[g][m] for m in B), not all(Rp[g][m] for m in B)
            assert not (yes and no)
            if yes: assert all(has)
            elif no: assert not any(has)
            else: assert any(has) and not all(has)
            n72 += 1
    for _ in range(3):                                                                        # Cor 7.3
        Rimp = [[(rg.random() < 0.5) if miss[g][m] else (shat[g][m] > 0) for m in range(nM)] for g in range(nG)]
        assert all(Rs[g][m] <= Rimp[g][m] <= Rp[g][m] for g in range(nG) for m in range(nM))
        cert_int = {B for B in powerset(nM) if t <= rho_int(shat, B)}
        true_int = {up(Rstar, down(Rstar, B, nG), nM) for B in powerset(nM)}
        for B in {up(Rimp, down(Rimp, B, nG), nM) for B in powerset(nM)}:
            if B not in true_int: assert B not in cert_int
        for B in powerset(nM):
            for g in range(nG):
                if all(Rimp[g][m] for m in B) != all(Rstar[g][m] for m in B):
                    assert not all(Rs[g][m] for m in B) and all(Rp[g][m] for m in B)
    n71 += 1
ok(f"Prop. 7.1 and Cor. 7.3 on {n71} random tables with missing cells, errors below t and infinite entries")
ok(f"Prop. 7.2: exact three-way answers and attained support bounds on {n72} (table, B, g) cases")

print("\nALL STATEMENTS VERIFIED")
