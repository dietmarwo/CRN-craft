# Review check: Direct SSA against the exact chemical master equation

Declared on 1 November 2026 during the iteration 4 review, before any of these
simulations were run. `validation/check_cme.py` implements it.

The delivered analytic campaign uses three zero/first-order fixtures. This check
covers the rate paths those cannot: second-order mass action with the binomial
homodimer divisor, a heterodimer, Michaelis–Menten and Hill propensities, a null
reaction, and transient as well as near-stationary times.

## Fixtures

All propensities use the stochastic-count convention c · ∏ C(x, ν).

| Fixture | Reactions | Initial state | Tested marginal |
|---|---|---|---|
| F1 dimerization | 2A → B (c 0.1), B → 2A (c 1.0), null B → B (c 2.0) | A 20, B 0 | B |
| F2 heterodimer | A + B → C (c 0.05), C → A + B (c 0.5) | A 15, B 10, C 0 | C |
| F3 Hill/MM loop | S → P, Michaelis–Menten on S (vmax 5, km 20); P → S, Hill activation by P (basal 0.2, vmax 3, k 10, n 2.5) | S 30, P 0 | P |
| F4 open | ∅ → X (c 10), 2X → ∅ (c 0.2) | X 0 | X |

Save times 0, 0.5, 2 and 10 s. 100,000 replicates each, 8 threads, root seeds
2026110451 to 2026110454.

## Exact law

The reachable state space is enumerated from the initial state, with F4
truncated at X ≤ 80. The probability vector at each save time is computed by
uniformization of the CME generator, with propensities written independently
in Python. In review, the oracle was compared with `scipy.linalg.expm` (declared
agreement 1e-12) and with a 40-digit mpmath matrix exponential.

## Test

One marginal per fixture at the three positive save times: 12 tests. Each
uses the Dvoretzky–Kiefer–Wolfowitz bound with a Bonferroni family
α = 0.001, so ε = sqrt(ln(2 · 12 / 0.001) / (2 · 100000)) ≈ 0.00710. Every
maximum CDF error must be ≤ ε. The DKW bound is conservative for discrete laws.

Power check: the F1 and F4 samples are also compared with laws built with the
falling-factorial convention (no 1/2 divisor for 2A or 2X). Each of those two
families must fail, that is, at least one of its three tests must exceed ε.

No threshold changes after seeing results; failed rows are retained.
