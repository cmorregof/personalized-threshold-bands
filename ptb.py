"""
ptb.py -- reference implementation of personalized threshold bands for
incomplete formal contexts, in exact rational arithmetic.

Companion code of: M. Castaneda, J. C. Riano-Rojas, C. M. Orrego Franco,
"Personalized threshold bands for incomplete formal contexts: tolerance,
certified concepts and stability".

Conventions (as in the paper)
-----------------------------
s                 : |G| x |M| table of certainty budgets s = (v - c) / delta
level relations   : R_{>=u} = {(g,m) : s(g,m) >= u},  R_{>u} = {(g,m) : s(g,m) > u}
tolerance t > 0   : sure relation {s >= t}, possible relation {s > -t},
                    negative relation {s <= -t}
inn(s, X, B)      : inner level  = min of s over the rectangle X x B
out(s, X, B)      : outer level  = strongest reason for X x B not to be maximal
rho(s, X, B)      : robustness   = min(inn, -out)
Objects and attributes are indexed 0..n-1; sets are frozensets.
"""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations, product
import math

INF = math.inf


# ----------------------------------------------------------------------------
# basic set helpers (objects / attributes are indexed 0..n-1, sets are frozensets)
# ----------------------------------------------------------------------------
def powerset(n):
    for r in range(n + 1):
        for c in combinations(range(n), r):
            yield frozenset(c)


def up(R, X, nM):
    """X' : attributes shared by all objects of X in the boolean relation R."""
    return frozenset(m for m in range(nM) if all(R[g][m] for g in X))


def down(R, B, nG):
    """B' : objects having all attributes of B in the boolean relation R."""
    return frozenset(g for g in range(nG) if all(R[g][m] for m in B))


def concepts(R):
    """All formal concepts (X,B) of a boolean matrix R (list of lists)."""
    nG, nM = len(R), len(R[0])
    out = set()
    for B0 in powerset(nM):
        X = down(R, B0, nG)
        B = up(R, X, nM)
        out.add((X, B))
    return out


# ----------------------------------------------------------------------------
# relations induced by a budget matrix
# ----------------------------------------------------------------------------
def level_rel(s, u, strict=False):
    if strict:
        return [[x > u for x in row] for row in s]
    return [[x >= u for x in row] for row in s]


def sure_rel(s, t):
    assert t > 0
    return level_rel(s, t)


def poss_rel(s, t):
    assert t > 0
    return level_rel(s, -t, strict=True)


def incomplete_context(s, t):
    """three-valued table at tolerance t>0: '+', '-', '?'."""
    return [['+' if x >= t else ('-' if x <= -t else '?') for x in row] for row in s]


def critical_values(s):
    vals = sorted({abs(x) for row in s for x in row if x != 0})
    return vals


# ----------------------------------------------------------------------------
# SE-ISI concepts (Yao) of an incomplete context
# ----------------------------------------------------------------------------
def seisi_concepts(Rs, Rp):
    """Triples (X,B,C): B = X^sure, C = X^poss, X = B^sure  ∩  C^poss."""
    nG, nM = len(Rs), len(Rs[0])
    out = set()
    for X in powerset(nG):
        B = up(Rs, X, nM)
        C = up(Rp, X, nM)
        clo = down(Rs, B, nG) & down(Rp, C, nG)
        if clo == X:
            out.add((X, B, C))
    return out


def apposition(Rs, Rp):
    return [list(Rs[g]) + list(Rp[g]) for g in range(len(Rs))]


# ----------------------------------------------------------------------------
# inner level, outer level and robustness of a rectangle (X,B)
# ----------------------------------------------------------------------------
def phi(s, B, g):
    return min((s[g][m] for m in B), default=INF)


def psi(s, X, m):
    return min((s[g][m] for g in X), default=INF)


def inn(s, X, B):
    """Inner level: the weakest cell inside the rectangle X x B (+inf if empty)."""
    return min((s[g][m] for g in X for m in B), default=INF)


def out(s, X, B):
    """Outer level: the strongest outside object or attribute (-inf if none)."""
    nG, nM = len(s), len(s[0])
    a = max((phi(s, B, g) for g in range(nG) if g not in X), default=-INF)
    b = max((psi(s, X, m) for m in range(nM) if m not in B), default=-INF)
    return max(a, b)


def rho(s, X, B):
    """Robustness: (X,B) is a concept of every completion of K(t) iff t <= rho."""
    return min(inn(s, X, B), -out(s, X, B))


# names used by the verification scripts: the sure profile is (max(out,0), inn]
death = inn
birth = out


# ----------------------------------------------------------------------------
# brute-force existence sets over the critical grid
# ----------------------------------------------------------------------------
def grid_pieces(s):
    """Pieces of the critical grid on which K(t) is constant: (lo, hi], hi possibly inf.
    Returns list of (lo, hi, sample_t)."""
    lam = [0] + critical_values(s)
    pieces = []
    for k in range(1, len(lam)):
        pieces.append((lam[k - 1], lam[k], lam[k]))          # sample at right end (included)
    pieces.append((lam[-1], INF, lam[-1] + 1))
    return pieces


def existence_sets(s, mode):
    """Profiles by brute force: dict signature -> list of maximal intervals (lo, hi]
    of tolerances, for mode in {'sure', 'possible', 'seisi'}."""
    pieces = grid_pieces(s)
    alive = {}
    for idx, (lo, hi, t) in enumerate(pieces):
        Rs, Rp = sure_rel(s, t), poss_rel(s, t)
        if mode == 'sure':
            sigs = concepts(Rs)
        elif mode == 'possible':
            sigs = concepts(Rp)
        elif mode == 'seisi':
            sigs = seisi_concepts(Rs, Rp)
        else:
            raise ValueError(mode)
        for sg in sigs:
            alive.setdefault(sg, []).append(idx)
    bars = {}
    for sg, idxs in alive.items():
        runs, start, prev = [], idxs[0], idxs[0]
        for i in idxs[1:]:
            if i != prev + 1:
                runs.append((start, prev))
                start = i
            prev = i
        runs.append((start, prev))
        bars[sg] = [(pieces[a][0], pieces[b][1]) for a, b in runs]
    return bars


# ----------------------------------------------------------------------------
# Dowker filtration
# ----------------------------------------------------------------------------
def dowker_F(s, sigma):
    """F(sigma) = max over attributes of the weakest budget within sigma."""
    nM = len(s[0])
    return max(min(s[g][m] for g in sigma) for m in range(nM))


def dowker_persistence(s, side='objects', max_dim=None):
    """Persistence of the super-level filtration of F on the full simplex.
    Returns list of (dim, u_birth, u_death) with u_birth >= u_death (levels),
    u_death = -inf for essential classes."""
    import gudhi
    if side == 'attributes':
        s = [list(col) for col in zip(*s)]
    n = len(s)
    st = gudhi.SimplexTree()
    for r in range(1, n + 1):
        for sig in combinations(range(n), r):
            st.insert(list(sig), filtration=float(-dowker_F(s, sig)))
    st.make_filtration_non_decreasing()
    st.compute_persistence(homology_coeff_field=2, min_persistence=-1.0)
    out = []
    for dim, (bv, dv) in st.persistence(homology_coeff_field=2, min_persistence=-1.0):
        out.append((dim, -bv, (-dv if dv != INF else -INF)))
    return out


def betti_of_complex(simplices):
    import gudhi
    st = gudhi.SimplexTree()
    for sx in simplices:
        st.insert(list(sx), filtration=0.0)
    st.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
    return st.betti_numbers()


def dowker_complex(R):
    nG, nM = len(R), len(R[0])
    return [sig for sig in powerset(nG) if sig and any(all(R[g][m] for g in sig) for m in range(nM))]


def concept_poset_order_complex(R):
    """Order complex of Q(R) = concepts with non-empty extent and non-empty intent."""
    Q = [c for c in concepts(R) if c[0] and c[1]]
    Q = sorted(Q, key=lambda c: (len(c[0]), sorted(c[0])))
    idx = {c: i for i, c in enumerate(Q)}
    chains = []

    def extend(chain):
        chains.append(tuple(idx[c] for c in chain))
        last = chain[-1]
        for c in Q:
            if last[0] < c[0]:
                extend(chain + [c])
    for c in Q:
        extend([c])
    return Q, chains


def rho_ext(s, X):
    """Largest tolerance at which X is an extent of every completion (Prop. 5.11)."""
    nG, nM = len(s), len(s[0])
    return min((max(min(psi(s, X, m), -s[g][m]) for m in range(nM))
                for g in range(nG) if g not in X), default=INF)


def rho_int(s, B):
    """Dual of rho_ext: largest tolerance at which B is an intent of every completion."""
    st = [list(col) for col in zip(*s)]
    return rho_ext(st, B)
