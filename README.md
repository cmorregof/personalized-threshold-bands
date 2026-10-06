# Personalized threshold bands for incomplete formal contexts

Companion code for the paper

> Mateo Castañeda, Juan Carlos Riaño-Rojas and Carlos Manuel Orrego Franco.
> *Personalized threshold bands for incomplete formal contexts: tolerance, certified concepts and stability.*
> Universidad Nacional de Colombia, Sede Manizales, 2026.

Each object–attribute pair receives its own threshold band. A reading above the band counts as present, one below it as absent, and one inside it as unknown. A single tolerance `t > 0` widens or narrows every band at once, and the whole construction depends on the data only through the certainty table `s = (v - c) / δ`. At tolerance `t` the sure, possible and negative relations are `{s ≥ t}`, `{s > -t}` and `{s ≤ -t}`.

All computations use exact rational arithmetic (`fractions.Fraction`), so equalities and boundary cases are decided exactly. Persistent homology is computed with [GUDHI](https://gudhi.inria.fr), which takes floating-point filtration values. `verify_theorems.py` checks that the conversion from the exact values is injective on the values that occur, so the order of every filtration, and with it every persistence diagram, is exact.

## Contents

| File | Purpose | Running time |
|---|---|---|
| `ptb.py` | Reference implementation: tolerance contexts, sure, possible and SE-ISI concepts, inner and outer levels, robustness `rho` and `rho_ext`, profiles over the critical grid, Dowker filtration | |
| `verify_examples.py` | Recomputes every worked example of the paper, including Tables 1 and 2 | about 1 s |
| `verify_theorems.py` | Exhaustive and randomized checks of Sections 3 to 6; prints the counts quoted in Remark 5.20 | about 30 s |
| `verify_regression.py` | A broader randomized suite on tie-heavy tables with zeros | about 2 min |
| `make_figure.py` | Draws Figure 1, with English or Spanish labels | a few seconds |
| `figures/` | Figure 1 in English and Spanish, as PDF and PNG | |

## Quick start

```bash
pip install -r requirements.txt
python3 verify_examples.py
python3 verify_theorems.py
python3 verify_regression.py
python3 make_figure.py figures/fig_crown.pdf
FIG_LANG=es python3 make_figure.py figures/fig_crown_es.pdf
```

Each script stops with an `AssertionError` at the first statement that fails, and ends with `ALL ... VERIFIED` or `ALL CHECKS PASSED` when everything holds.

## Where each statement is checked

Statement numbers are those of the paper.

| Statement | Script | What is checked |
|---|---|---|
| Examples 2.1 and 2.2 | `verify_examples.py` | The four concepts of the animals context; the pair `({duck, fish}, {swims})` survives in exactly one of four completions |
| Proposition 2.3 | `verify_regression.py` | SE-ISI extents coincide with the extents of the apposition, on every piece of the critical grid |
| Example 3.3, Table 1 | `verify_examples.py`, `verify_theorems.py` | The budgets follow from the readings, baselines and half-widths |
| Remark 3.7 | `verify_examples.py` | Unrestricted bands reproduce any prescribed incomplete context |
| Proposition 3.8, Example 3.9 | `verify_theorems.py`, `verify_examples.py` | The expressiveness criterion on 3 000 random columns |
| Proposition 4.1 | `verify_examples.py` | Nine critical values and ten distinct contexts in the running example |
| Theorem 5.2, Example 5.3 | `verify_theorems.py`, `verify_regression.py`, `verify_examples.py` | Closed-form sure and possible profiles on 314 000 cases, including entries `±∞` |
| Proposition 5.4, Theorem 5.5 | `verify_regression.py` | Inner and outer levels and robustness move by at most `ε`, without genericity |
| Theorem 5.7, Propositions 5.8 and 5.9 | `verify_theorems.py`, `verify_regression.py` | Certification, its three-way reading and the bridge to Li–Mei–Lv approximate concepts, enumerating every completion |
| Proposition 5.11 | `verify_regression.py`, `verify_examples.py` | Certified extents and the formula for `rho_ext`, enumerating every completion |
| Proposition 5.12, Theorem 5.14 | `verify_theorems.py`, `verify_regression.py` | The SE-ISI profile formula on 517 568 cases; conditional stability on 150 generic tables |
| Example 5.13, Proposition 5.15 | `verify_examples.py`, `verify_theorems.py`, `verify_regression.py` | The split at a tie and the sharpness of the radius `η/2` |
| Theorem 5.17, Corollary 5.18 | `verify_regression.py` | Betti numbers of the Dowker complex and of the order complex of the concept poset agree; Dowker duality; bottleneck stability |
| Theorem 5.19 | `verify_theorems.py` | Certified Betti numbers on 5 928 cases; the lower bound on 87 638 enumerated completions |
| Section 6, Table 2, Figure 1 | `verify_examples.py`, `verify_theorems.py`, `make_figure.py` | Every number of the running example and of its tie variant |
| Section 7 | `verify_examples.py` | The closure computations of the modal-style example |

## Conventions in `ptb.py`

Objects and attributes are indexed from `0`, sets are `frozenset`s, and a table `s` is a list of rows. `inn(s, X, B)` is the inner level, the weakest cell of the rectangle `X × B`; `out(s, X, B)` is the outer level, the strongest outside object or attribute; `rho(s, X, B) = min(inn, -out)` is the robustness. The sure profile of `(X, B)` is `(max(out, 0), inn]` and the possible profile is `(max(-inn, 0), -out]`. `existence_sets(s, mode)` computes the same profiles by brute force over the critical grid, for `mode` in `'sure'`, `'possible'` and `'seisi'`, and is what the scripts compare the formulas against.

## License

MIT; see `LICENSE`.
