# Continuous variables and the `hybrid` method: validation

Script: `validation/check_hybrid.py` (installed wheel; standard library and
NumPy only). Reports: `validation/hybrid/hybrid-gates.json` for the default
sizes and `validation/hybrid/hybrid-gates-large.json` for `--large`. The
script prints run times; the reports hold fixed-work evidence only.

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

First runs, 2026-11-06, on the 2-vCPU development container with the wheel
of the engine stream (commit `9f59917`). **Every gate, check and test
passed at both sizes; no failed row to retain.** The default report was
also produced under `python -O` and is identical. After the following
commit (a floor on the step cap and clamp ties, edge cases only), both
reports were reproduced byte for byte.

### A. Exact-law gates (20,000 replicates each, ε = 0.01517)

| Gate | Largest distance D | Outcome |
|---|---|---|
| A1 decaying variable drives a death | 0.00492 | pass |
| A2 constant derivative | 0.00399 | pass |
| A3 bound stop | 0.00265 | pass |
| A4 bound clamp | 0.00155 | pass |
| A5 no continuous variable | 0.00315 | pass |

The largest distance, 0.00492, is about a third of ε.

| Deterministic check | Outcome | Worst value |
|---|---|---|
| `A1_complete` | pass | |
| `A1_observed` | pass | |
| `A2_complete` | pass | |
| `A3_all_stopped` | pass | |
| `A3_stop_time` | pass | 1.78e-15 |
| `A3_x_before_stop` | pass | |
| `A3_x_after_stop` | pass | |
| `A4_all_completed` | pass | |
| `A4_x_clamped` | pass | |
| `A5_complete` | pass | |
| `A6_ode_accuracy` | pass | 0.0216 |

`A3_stop_time` gives the largest |stop time − 3| in hours;
`A6_ode_accuracy` the largest error as a fraction of 100 (atol + rtol |x|),
so the integrator stays within 2.2 (atol + rtol |x|) of the exact
solutions over ten time units.

Power (same samples against wrong laws, not gating):

| Gate | Wrong law | D | Against ε |
|---|---|---|---|
| A1 | c = 1.1 | 0.0372 | rejected |
| A2 | b = 2.2 | 0.1155 | rejected |
| A5 | immigration 10.5 | 0.0580 | rejected |

### B. ATP equivalence with exact thinning

Default size, 400 replicates per arm (seeds 7101 and 7102):

| Arm | Events | Stopped | Taken over | Median takeover (h) | Mean final ATP | Mean final N |
|---|---|---|---|---|---|---|
| thinning | 60,721,498 | 41 | 41 | 14030 | 404.1 | 5011 |
| hybrid | 60,862,156 | 42 | 42 | 15287 | 402.8 | 5176 |

| Test | Statistic | p | Outcome (p ≥ 0.0025) |
|---|---|---|---|
| B1_takeover_time | D = 0.0150 | 1.000 | pass |
| B2_final_atp | D = 0.0475 | 0.748 | pass |
| B3_final_population | D = 0.0550 | 0.569 | pass |
| B4_atp_exhaustion | 41 vs 42 stopped | 1.000 | pass |

Large size (`--large`), 4,000 replicates per arm (seeds 7201 and 7202):

| Arm | Events | Stopped | Taken over | Median takeover (h) | Mean final ATP | Mean final N |
|---|---|---|---|---|---|---|
| thinning | 606,956,342 | 416 | 416 | 15450 | 402.2 | 5106 |
| hybrid | 608,765,560 | 476 | 476 | 15705 | 396.1 | 5656 |

| Test | Statistic | p | Outcome (p ≥ 0.0025) |
|---|---|---|---|
| B1_takeover_time | D = 0.0150 | 0.756 | pass |
| B2_final_atp | D = 0.0150 | 0.756 | pass |
| B3_final_population | D = 0.0192 | 0.445 | pass |
| B4_atp_exhaustion | 416 vs 476 stopped | 0.036 | pass |

In every arm the number of replicates stopped by ATP exhaustion equals the
number whose mutants reached 90% (41 and 42, 416 and 476); a likely reason
is that once the wild type is gone, production f·wt no longer covers
maintenance v2·N. B1 and B4 then measure closely related events (the
per-replicate sets were not compared). The smallest p-value of the eight ATP tests is
0.036 (B4 at the large size: 10.4% against 11.9% stopped, a difference of
2.1 standard errors). It passes the predeclared level 0.0025; it is the
first place to look if a larger independent run is made. Means of final N
differ by about 11% at the large size, which the KS test of the whole
distribution does not flag (p = 0.45); the population distribution is
heavy-tailed because takeover lets N grow.

Run time for 400 replicates on two threads: thinning about 6 s, hybrid about
62 s (11 times as long: every event needs a Dormand–Prince step of seven
right-hand-side evaluations). The large gate took about 11 minutes.

Limitations: the gates cover one variable at a time with closed-form laws
and the ATP model at one set of parameters and the default tolerances; they
do not measure the error as a function of rtol and atol, and the two-sample
tests detect only differences of a few percent in the ATP model at these
sizes.
