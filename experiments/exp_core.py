"""exp_core.py -- shared code for the experiments of Section 7: imputation versus certification.

A missing value is treated as a cell of total doubt ('?' at every tolerance). The true
context (complete, clean data binarized with the same fixed thresholds) is then always one
of the completions, so every certified intent or concept is a true one, whatever the
missingness mechanism. Imputation methods fill the gaps and binarize; they report more
patterns but some of them are wrong. The experiments measure both sides of that trade-off.
"""
import numpy as np
from sklearn.datasets import load_wine, load_breast_cancer, load_diabetes
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import SimpleImputer, KNNImputer, IterativeImputer


# ---------------------------------------------------------------- data and thresholds
def load_datasets():
    bc = load_breast_cancer()
    return {
        'wine': load_wine().data,                 # 178 x 13
        'breast_cancer': bc.data[:, :10],         # 569 x 10, the ten "mean" features
        'diabetes': load_diabetes().data,         # 442 x 10
    }


def fixed_thresholds(V):
    """One threshold per attribute, fixed in advance: the midpoint between the two distinct
    values that bracket the median, so that no observation lies exactly on it."""
    c = np.empty(V.shape[1])
    for m in range(V.shape[1]):
        u = np.unique(V[:, m])
        k = np.searchsorted(u, np.median(V[:, m]), side='right')
        lo = u[k - 1]
        hi = u[k] if k < len(u) else u[k - 1] + 1.0
        c[m] = 0.5 * (lo + hi)
        assert not np.any(V[:, m] == c[m])
    return c


def half_widths(V, a=0.1):
    q75, q25 = np.percentile(V, [75, 25], axis=0)
    return a * (q75 - q25)


# ---------------------------------------------------------------- missingness mechanisms
def _ranks(V):
    N = V.shape[0]
    return (np.argsort(np.argsort(V, axis=0, kind='stable'), axis=0) + 0.5) / N


def make_mask(V, mech, p, rng):
    """Boolean mask of missing cells with overall rate close to p."""
    N, n = V.shape
    if p == 0:
        return np.zeros((N, n), bool)
    if mech == 'MCAR':
        return rng.random((N, n)) < p
    r = _ranks(V)
    if mech == 'MNAR':                      # high values are more likely to be missing
        return rng.random((N, n)) < np.clip(2 * p * r, 0, 1)
    if mech == 'MAR':                       # missingness driven by fully observed attributes
        n_drv = max(1, int(round(0.3 * n)))
        drivers = rng.choice(n, n_drv, replace=False)
        targets = [m for m in range(n) if m not in set(drivers)]
        peff = min(p * n / len(targets), 0.49)
        M = np.zeros((N, n), bool)
        for m in targets:
            d = rng.choice(drivers)
            M[:, m] = rng.random(N) < np.clip(2 * peff * r[:, d], 0, 1)
        return M
    raise ValueError(mech)


# ---------------------------------------------------------------- imputation baselines
def impute(Vobs, method, seed):
    mu = np.nanmean(Vobs, axis=0)
    sd = np.nanstd(Vobs, axis=0)
    sd[sd == 0] = 1.0
    Z = (Vobs - mu) / sd
    if method == 'mean':
        Zi = SimpleImputer(strategy='mean').fit_transform(Z)
    elif method == 'median':
        Zi = SimpleImputer(strategy='median').fit_transform(Z)
    elif method == 'knn':
        Zi = KNNImputer(n_neighbors=5).fit_transform(Z)
    elif method == 'mice':
        Zi = IterativeImputer(max_iter=10, random_state=seed).fit_transform(Z)
    else:
        raise ValueError(method)
    return Zi * sd + mu


# ---------------------------------------------------------------- lattices with bit tricks
_ALLB = {}


def all_subsets(n):
    if n not in _ALLB:
        codes = np.arange(1 << n)
        bits = ((codes[:, None] >> np.arange(n)) & 1).astype(np.float32)
        _ALLB[n] = (bits, (1 << np.arange(n)).astype(np.int64))
    return _ALLB[n]


def derive_all(R):
    """For every attribute subset B (indexed by its bit code): sel[B] = B' (bool, objects)."""
    bits, _ = all_subsets(R.shape[1])
    return (bits @ (~R).astype(np.float32).T) == 0


def closed_intents(R):
    """Intents (closed attribute sets) of the context R, with extents and supports."""
    bits, pw = all_subsets(R.shape[1])
    sel = derive_all(R)
    lack = sel.astype(np.float32) @ (~R).astype(np.float32)
    closure = (lack == 0).astype(np.int64) @ pw          # code of B''
    intents = np.unique(closure)
    return set(intents.tolist()), sel


def concepts_of(R):
    intents, sel = closed_intents(R)
    return {(sel[i].tobytes(), i) for i in intents}, intents, sel


def certified_intents(S, Neg):
    """B is an intent of every completion iff every m outside B has a witness g that surely
    has all of B and surely lacks m (Proposition 5.11(ii))."""
    bits, _ = all_subsets(S.shape[1])
    selS = derive_all(S)
    wit = (selS.astype(np.float32) @ Neg.astype(np.float32)) > 0
    return np.all(wit | (bits > 0), axis=1), selS


def _size(code):
    return bin(code).count('1')


def recall_by_size(Tf, found):
    out = {}
    for lab, cond in (('1', lambda k: k == 1), ('2', lambda k: k == 2), ('3', lambda k: k == 3), ('4p', lambda k: k >= 4)):
        tot = [i for i in Tf if cond(_size(i))]
        out['rec_k' + lab] = (sum(1 for i in tot if i in found) / len(tot)) if tot else np.nan
        out['n_k' + lab] = len(tot)
    return out


# ---------------------------------------------------------------- evaluation
def evaluate_point(R_true_info, R_hat, s0):
    """Metrics for a method that returns an ordinary (complete) context R_hat."""
    T_int, T_sel, T_con = R_true_info
    con_hat, int_hat, sel_hat = concepts_of(R_hat)
    supp_hat = sel_hat.sum(1)
    supp_true = T_sel.sum(1)
    Tf = {i for i in T_int if i != 0 and supp_true[i] >= s0}
    Rf = {i for i in int_hat if i != 0 and supp_hat[i] >= s0}
    prec = len(Rf & T_int) / len(Rf) if Rf else np.nan
    found = Tf & int_hat
    rec = len(found) / len(Tf) if Tf else np.nan
    # extent error on the recovered frequent intents: Jaccard distance of supports
    errs = []
    for i in found:
        a, b = sel_hat[i], T_sel[i]
        errs.append(1 - (a & b).sum() / max((a | b).sum(), 1))
    supp_err = float(np.mean(errs)) if errs else np.nan
    # exact concepts (extent and intent), frequent and non-trivial
    Tc = {c for c in T_con if c[1] != 0 and supp_true[c[1]] >= s0}
    Rc = {c for c in con_hat if c[1] != 0 and supp_hat[c[1]] >= s0}
    cprec = len(Rc & T_con) / len(Rc) if Rc else np.nan
    crec = len(Tc & con_hat) / len(Tc) if Tc else np.nan
    return dict(prec_int=prec, rec_int=rec, supp_err=supp_err, n_rep=len(Rf), n_true=len(Tf),
                prec_con=cprec, rec_con=crec, interval_width=np.nan, covered=np.nan,
                **recall_by_size(Tf, found))


def evaluate_certified(R_true_info, S, Neg, s0):
    """Metrics for certification: S = sure relation, Neg = sure-negative relation."""
    T_int, T_sel, T_con = R_true_info
    P = ~Neg
    cert, selS = certified_intents(S, Neg)
    selP = derive_all(P)
    lower, upper = selS.sum(1), selP.sum(1)
    supp_true = T_sel.sum(1)
    covered = bool(np.all((lower <= supp_true) & (supp_true <= upper)))
    cert_set = set(np.flatnonzero(cert).tolist())
    Tf = {i for i in T_int if i != 0 and supp_true[i] >= s0}
    Cf = {i for i in cert_set if i != 0 and upper[i] >= s0}
    prec = len(Cf & T_int) / len(Cf) if Cf else np.nan
    found = Tf & cert_set
    rec = len(found) / len(Tf) if Tf else np.nan
    widths = [(upper[i] - lower[i]) / supp_true[i] for i in found]
    width = float(np.mean(widths)) if widths else np.nan
    # certified concepts = concepts of both extreme completions (Theorem 5.7)
    conS, _, _ = concepts_of(S)
    conP, _, _ = concepts_of(P)
    certc = conS & conP
    Tc = {c for c in T_con if c[1] != 0 and supp_true[c[1]] >= s0}
    Cc = {c for c in certc if c[1] != 0 and upper[c[1]] >= s0}
    cprec = len(Cc & T_con) / len(Cc) if Cc else np.nan
    crec = len(Tc & certc) / len(Tc) if Tc else np.nan
    return dict(prec_int=prec, rec_int=rec, supp_err=np.nan, n_rep=len(Cf), n_true=len(Tf),
                prec_con=cprec, rec_con=crec, interval_width=width, covered=covered,
                **recall_by_size(Tf, found))
