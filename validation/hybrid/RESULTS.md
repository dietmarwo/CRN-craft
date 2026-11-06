# Continuous variables and the `hybrid` method: validation

Script: `validation/check_hybrid.py` (installed wheel; standard library and
NumPy only). Reports: `validation/hybrid/hybrid-gates.json` for the default
sizes and `validation/hybrid/hybrid-gates-large.json` for `--large`, each
identical under `python -O`. The script prints run times; the reports hold
fixed-work evidence only.

## Protocol (declared before the first run)

Every `hybrid` run uses the Python defaults `rtol = 1e-6`, `atol = 1e-9`,
two threads and rules-run/0.5. No target, size, seed, test or alpha is
changed after a run; every row is kept, including failures.

### A. Exact-law gates

Each gate runs 20,000 replicates. The distribution in the last column is
compared with the empirical distribution by the one-sample Kolmogorov
distance D, the largest |F_n − F| over the support (for a time censored at
the horizon T, over [0, T]). One DKW family of K = 5 comparisons at family
alpha 0.001: a gate passes when D ≤ ε = sqrt(ln(2K/alpha) / (2n)) = 0.01517.

| Gate | Model (seed) | Exact law |
|---|---|---|
| A1 decaying variable drives a death | x' = −k x, k = 0.5, x(0) = 2; one molecule dies at rate c·x, c = 1; T = 10 (6101) | death time: F(t) = 1 − exp(−c x(0) (1 − e^−kt) / k), defective (no death with probability e^−4) |
| A2 constant derivative | x' = b, b = 2, x(0) = 0; immigration with propensity `"x"`; T = 3 (6102) | count at T: Poisson(b T² / 2 = 9) |
| A3 bound stop | x' = −1, x(0) = 3, lower bound 0 with `stop`; 20 molecules die at rate 0.2 each; T = 5 (6103) | count at the stop: Binomial(20, e^−0.6) |
| A4 bound clamp | x' = −1, x(0) = 2, lower bound 0 with `clamp`; immigration with propensity `"x"`; T = 5 (6104) | count at T: Poisson(∫ x dt = 2) |
| A5 no continuous variable | immigration 10, death 1 per molecule, from 0; T = 3 (6105) | count at T: Poisson(10 (1 − e^−3)) |

Deterministic checks, each a pass or fail row: in A3 every replicate stops
at |t − 3| ≤ 1e-9 and reports x = 1.5 at t = 1.5 (within 1e-9) and x = 0 at
t = 4 and 5; in A4 every replicate completes and reports x = 0 at t = 3 and
5; with no event (A6: x' = −x from 1, z' = z (1 − z) from 0.1), x and z at
t = 0, 1, ..., 10 lie within 100 (atol + rtol |exact|) of e^−t and
1 / (1 + 9 e^−t).

Power (reported, not gating): the same samples against wrong laws, c = 1.1
in A1, b = 2.2 in A2 and immigration 10.5 in A5; each should give D > ε.

### B. ATP equivalence with exact thinning

`thinning` on `models/rules/mtdna-three-class-atp.rules.json` (rules-model
0.2, the legacy `atp` block) against `hybrid` on
`models/rules/mtdna-three-class-atp-continuous.rules.json` (rules-model 0.3,
continuous variable `ATP`). Horizon 25,920 hours, report times 0, 8,640,
17,280 and 25,920, one tracker on `mutants` with first passage
`N > 0 && mutants >= 0.9 * N`, budgets of 10^8 events and 10^9 sampling
steps per replicate.

| Size | Replicates per arm | Seeds (thinning, hybrid) |
|---|---|---|
| default (CI) | 400 | 7101, 7102 |
| `--large` (optional) | 4,000 | 7201, 7202 |

| Test | Statistic | Test |
|---|---|---|
| B1 takeover time | first passage of the tracker, +∞ when it never holds | two-sample Kolmogorov–Smirnov, asymptotic p-value (Stephens' correction) |
| B2 final ATP | `ATP` at 25,920 hours (0 after an ATP stop) | same |
| B3 final population | `N` at 25,920 hours | same |
| B4 ATP exhaustion | share of replicates `stopped` | Fisher's exact test, two-sided |

One family of the four tests at family alpha 0.01, Bonferroni: each passes
when p ≥ 0.0025. Every replicate must complete or stop (no budget limit);
otherwise the row fails.

## Results

Pending the first run.
