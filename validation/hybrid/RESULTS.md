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

### C. Independent variables (rules-model/0.4): exact-law gates

Declared 2026-11-08, before the first run of the rules-model/0.4 engine.
Every model is rules-model/0.4, runs `hybrid` at the Python defaults
(`rtol = 1e-6`, `atol = 1e-9`), two threads and rules-run/0.5, with 20,000
replicates per gate. As in A, the one-sample Kolmogorov distance D against
the exact law is gated in one DKW family, here of K = 4 comparisons at
family alpha 0.001: a gate passes when D ≤ ε = sqrt(ln(2K/alpha) / (2n)) =
0.014989.

| Gate | Model (seed) | Exact law |
|---|---|---|
| C1 Gompertz mortality | m' = g m, g = 0.5, m(0) = 0.1, `independent: true`; 20 molecules die at per-molecule rate m (`rate: "m"`), no births; T = 4 (6201) | count at T: Binomial(20, S) with S = exp(−m(0)(e^{gT} − 1)/g) = 0.278647 |
| C2 dependent variable reads an independent one | as C1, plus a dependent x' = m from 0; the molecules die at rate x instead (6202) | count at T: Binomial(20, exp(−(m(0)/g)((e^{gT} − 1)/g − T))) = Binomial(20, 0.172800) |
| C3 independent stop bound | x' = −1 from 3, lower bound 0 with `stop`, `independent: true`; 20 molecules die at rate 0.2 each; T = 5 (6203) | count at the stop: Binomial(20, e^−0.6) |
| C4 eight-compartment pulse | u1' = −k u1 from 1, u_i' = k (u_{i−1} − u_i) from 0 for i = 2..8, all independent, k = 2; arrivals with propensity c max(0, u8), c = 50; deaths at rate d = 0.1 each, from 0; T = 12 (6204) | count at T: Poisson(μ), μ = c e^{−dT} (k/(k−d))^7 P(8, (k−d)T)/(k−d) = 11.34878 (P the regularized lower incomplete gamma function) |

Amendment before the first gate run: C4 was declared with the propensity
`c * u8`. A smoke run of three replicates (not a gate run, no statistic
computed) failed every replicate of both arms with `CRN-R030` at
`/rules/0/propensity`: near t = 0 the integrated u8, exactly 0 there and
positive after, comes out a little below 0, and a negative propensity fails
the replicate. The propensity is `c * max(0, u8)`, which has the same exact
law. No size, seed, test or alpha changed.

Deterministic checks, each a pass or fail row: C1, C2 and C4 complete;
in C1 the reported m is bit-identical across replicates at every report
time and within 100 (atol + rtol |m|) of m(0) e^{gt} at t = 0, 1, ..., 4;
in C3 every replicate stops at |t − 3| ≤ 1e-9, reports x = 1.5 at t = 1.5
(within 1e-9) and x = 0 at t = 4 and 5; the C4 model with 200 replicates
gives the same result with one and with two threads (apart from
`run_sha256`).

### D. Independent against dependent integration

The same four models with `independent` set to `false` on every variable,
run with their own seeds, against the C samples (the `true` arms). Each
compares the counts of C1–C4 (at T, or at the stop in D3) by the
two-sample Kolmogorov–Smirnov test with the asymptotic p-value of B. One
family of four tests at family alpha 0.01, Bonferroni: each passes when
p ≥ 0.0025. Every replicate must complete (D1, D2, D4) or stop (D3) in both
arms; otherwise the row fails.

| Test | Model | `false` arm (seed) |
|---|---|---|
| D1 | C1 | 20,000 replicates (6211) |
| D2 | C2 | 20,000 replicates (6212) |
| D3 | C3 | 20,000 replicates (6213) |
| D4 | C4 | 20,000 replicates (6214) |

### E. Speed of independent integration

Script `validation/benchmark_independent_ode.py` (timing, so not part of
the fixed-work report). Model C4 with 4,000 replicates, seed 6301, two
threads, the Python default tolerances: `rules.ensemble(...)` plus
`simulate()` timed for `independent: true` and for `false`, three runs
each in alternating order. Declared target: the median with `true` is at
least 2× faster than the median with `false`. The two arms' results are
reported (events, mean count at T), not compared, because they differ
within the tolerances.

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

### C–E. Independent variables (rules-model/0.4)

First run, 2026-11-08, on the 2-vCPU development container with the wheel
of the engine stream (I6), after the C4 amendment above. **Every gate,
check and test passed; no failed row to retain.** The report was also
produced under `python -O` and is identical, and its A and B sections are
identical to the earlier report: the rules-model/0.3 path did not change.

| Gate (20,000 replicates, ε = 0.014989) | Largest distance D | Outcome |
|---|---|---|
| C1 Gompertz mortality | 0.00229 | pass |
| C2 dependent variable reads an independent one | 0.00372 | pass |
| C3 independent stop bound | 0.00288 | pass |
| C4 eight-compartment pulse | 0.00459 | pass |

| Deterministic check | Outcome | Worst value |
|---|---|---|
| `C1_complete`, `C2_complete`, `C4_complete` | pass | |
| `C1_m_shared` | pass | |
| `C1_m_accuracy` | pass | 0.0043 |
| `C3_all_stopped` | pass | |
| `C3_stop_time` | pass | 8.9e-16 |
| `C3_x_before_stop`, `C3_x_after_stop` | pass | |
| `C4_threads` | pass | |

`C1_m_accuracy` is the largest error of the shared m as a fraction of
100 (atol + rtol |m|).

| Test (20,000 per arm) | Mean count (`true`, `false`) | D | p | Outcome (p ≥ 0.0025) |
|---|---|---|---|---|
| D1 C1 | 5.577, 5.551 | 0.00765 | 0.600 | pass |
| D2 C2 | 3.472, 3.440 | 0.00935 | 0.345 | pass |
| D3 C3 | 10.986, 10.961 | 0.00790 | 0.559 | pass |
| D4 C4 | 11.322, 11.337 | 0.00335 | 1.000 | pass |

The exact means are 5.573 (C1), 3.456 (C2), 10.976 (C3) and 11.349 (C4).

Speed (E, `validation/hybrid/independent-timing.json`): on model C4 with
4,000 replicates and two threads, the median of three runs took 0.195 s
with `independent: true` and 0.541 s with `false` (2.77×; declared target
2×, met). Both arms had 155,061 events and the same mean count at T,
11.414. On `examples/mtdna/14_independent_ode.py`, a 24-hour oscillator of
two variables that modulates the three-class model's degradation, the
times were about equal (1.2 s for 32 replicates of one year): there the
two extra components cost little next to the per-event work. The shared
solution was more accurate there, because it takes long steps over the
whole year, while per replicate the oscillator restarts at every event: the
largest error of the mean rhythm was 1.1e-4 against 1.8e-3.

Limitations: the gates use models with one to eight independent variables
and closed-form laws at the default tolerances. The two-sample tests detect
only differences of about 1.5% of the count distribution at these sizes.
The speed gain depends on the model: it is large when the independent
system is large or forces small steps, and negligible when events dominate.
