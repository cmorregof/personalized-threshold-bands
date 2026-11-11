"""test_rho_int.py -- brute-force check of rho_int_all and certified_intents (exp_core.py).

On small random certainty tables with ties and zeros, enumerates every completion of the
incomplete context K(t) at every tolerance t of the critical grid (and between grid points)
and checks that B is an intent of every completion  <=>  t <= rho_int(B)  <=>  certified_intents.
Also compares rho_int_all with the reference implementation ptb.rho_int.
Run:  python3 test_rho_int.py
"""
import itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from exp_core import rho_int_all, certified_intents
import ptb


def intents(R):
    G, n = R.shape
    out = set()
    for code in range(1 << n):
        B = [m for m in range(n) if code >> m & 1]
        ext = [g for g in range(G) if all(R[g, m] for m in B)]
        clo = sum(1 << m for m in range(n) if all(R[g, m] for g in ext))
        out.add(clo)
    return out


rng = np.random.default_rng(0)
checked = 0
for trial in range(300):
    G, n = rng.integers(1, 5), rng.integers(1, 4)
    s = rng.choice([-2, -1, -0.5, 0, 0.5, 1, 2], size=(G, n)).astype(float)
    rho = rho_int_all(s)
    sl = s.tolist()
    for code in range(1 << n):
        B = [m for m in range(n) if code >> m & 1]
        assert rho[code] == ptb.rho_int(sl, B), (s, B)
    grid = sorted({abs(x) for x in s.ravel() if x != 0} | {0.25, 0.75, 1.5, 3.0})
    for t in grid:
        if t <= 0:
            continue
        S, Neg = s >= t, s <= -t
        unk = [(g, m) for g in range(G) for m in range(n) if not S[g, m] and not Neg[g, m]]
        common = None
        for fill in itertools.product([False, True], repeat=len(unk)):
            R = S.copy()
            for (g, m), v in zip(unk, fill):
                R[g, m] = v
            I = intents(R)
            common = I if common is None else common & I
        cert, _, rho2 = certified_intents(S, Neg, s)
        for code in range(1 << n):
            assert (code in common) == (t <= rho[code]) == bool(cert[code]), (s, t, code)
            checked += 1
print(f'rho_int_all and certified_intents agree with brute force on {checked} (table, t, B) cases: OK')
