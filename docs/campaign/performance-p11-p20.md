# Performance iterations P11–P20: declaration

Declared on 1 November 2026, before any P11–P20 code or measurement. Base: the
reviewed campaign source (`complete-campaign` branch, P10 accepted). Each
iteration has a hypothesis, an implementation, a test that fails without it, and
a paired before/after measurement against the previous accepted binary. A
rejected candidate is reverted; its diff, test outcome and timings are kept.

## Exactness

P11–P18 and P20 must not change any numerical result. Their gates are:

- identical numerical digests on every workload below (SSA result without
  `run_sha256`; the complete FSP payload, including `work` and error fields;
  archive `result.json` bytes);
- the release test suite, the SSA and CME checkers byte-identical to their
  previous reports, and `simulate` equal to `simulate-cached`;
- the same random words in the same order, the same ordered sums, and the same
  diagnostic codes and messages on failure.

P19 changes an approximate method (tau leaping) and is judged statistically,
see below.

## Workloads

| ID | Workload | Used by |
|---|---|---|
| W1–W4 | The four frozen P01–P10 workloads (two channels, sparse 66, three-gene circuit, dense 66), seed 4242, 64 replicates, horizon 100, saves at 0 and 100, one worker, `simulate` and `simulate-cached` | all SSA iterations |
| W5 | Three-gene circuit, horizon 100, 2001 save points (every 0.05), statistics only | P15 |
| W6 | As W5, keeping all 64 trajectories | P15, P20 |
| W7 | Two channels, horizon 0.1, 64 replicates (the P06/P07 short-path control, 21 pairs) | P16 |
| W8a | FSP birth–death, 401 states, horizon 40 (the P10 control) | P17, P18 |
| W8b | FSP, two independent birth–death species with bounds 0..99 (10,000 states), from (0,0), horizon 5 | P17, P18 |
| W9 | Tau leaping, immigration–death with stationary mean 10,000, epsilon 0.05, horizon 5, 256 paths | P19 |

## Candidates

| Iteration | Hypothesis |
|---|---|
| P11 | Record the prefix sums during the strict total pass and select the channel by binary search over them; identical channel, O(log R) selection |
| P12 | `simulate-cached`: recompute the strict prefix sums, checks and positive minimum only from the lowest invalidated channel |
| P13 | `simulate-cached`: a recomputed rate that is bitwise unchanged does not invalidate the prefix |
| P14 | One scheduler poll per event instead of up to three |
| P15 | One flat save buffer per path instead of one vector per save point |
| P16 | Reuse per-worker rate, prefix and state buffers across replicates |
| P17 | FSP: skip rows whose current series term is exactly zero in the sparse product |
| P18 | FSP: column-oriented gather with the series accumulation fused into the same pass, same ascending summation order |
| P19 | Tau: transformed-rejection Poisson sampling (Hörmann's PTRS) for leap means above 16; the product method stays below |
| P20 | Archive: serialize the result once and reuse the bytes checked against the output budget |

## Timing protocol and acceptance

- Release builds on this host (2 vCPUs, shared), one worker, serial runs with no
  other measurement or build running. Seven alternating before/after pairs (21
  for W7, five for W8b and W9); core time from provenance, process time recorded
  as well.
- **Accept** if the median paired core speedup on the iteration's target
  workload(s) is at least 1.03 and the lower end of its paired bootstrap 95%
  interval (20,000 resamples, seed 20261101) is above 1.0, and no W1–W4
  workload regresses with an interval upper end below 0.97.
- **Otherwise reject**, revert the code and retain the evidence. A rejected
  iteration is not counted as an improvement.
- P19 additionally requires: (a) a sampler test of 200,000 draws at each of
  the means 16.5, 30, 100, 1000 and 100,000 from fixed streams, with a two-sided
  DKW band at family α = 0.001 (Bonferroni over the five means) and every sample
  mean within five standard errors; (b) `validation/check_tau.py` passing its
  unchanged predeclared envelopes; (c) a median W9 speedup of at least 1.5.
  Because the drawn values change, the tau algorithm identifier gets a new
  version.

## Amendments before measurement

- **P20 timer scope (recorded after implementing P20, before any P20
  measurement).** The core timer in `provenance.json` stops after the output
  budget check. Before P20 that check was a counting serialization and the
  archive serialized the result a second time after the timer stopped; P20
  keeps the checking serialization and reuses its bytes, so the work it
  removes lies outside the core timer. Core time cannot observe P20, so P20
  is judged on paired process time with the same thresholds, and its W1–W4
  regression guard uses process time as well. Core time is still recorded.
- **P14–P16 base (recorded before the confirmation run).** P14–P16 were
  measured on a base that still contained P13, which was rejected and
  reverted afterwards. A confirmation run compares the P16 binary (with P13)
  with the reverted source on W1–W4 under `simulate-cached` and is reported
  with the results; it is evidence for the revert, not a new candidate.

## Results

Measured on 1 November 2026 on the declared host (2 vCPUs, shared), release
builds, one worker, serial runs; the reviewer agent working on the optional
integrations paused while a `QUIET` file existed. Every exact iteration
reproduced the numerical digest of its predecessor on every workload (the
harness refuses a timing otherwise). Raw rows, digests and decisions are in
`validation/campaign/p11-p20/`; the timing logs are
`timing-p11-p16.log` and `timing-r13-p17-p20.log`.

Paired median speedup (after / before, >1 is faster) with the bootstrap 95%
interval. Targets are in bold.

| Iteration | Command | Targets | W1 | W2 | W3 | W4 | Decision |
|---|---|---|---|---|---|---|---|
| P11 binary search over prefix sums | cached | **W2** 1.073 [1.070, 1.107], **W4** | 1.048 [1.046, 1.053] | **1.073** | 0.988 [0.971, 1.036] | **1.003** [0.987, 1.018] | accepted |
| P11 (same binaries) | simulate | **W2** 1.011 [0.994, 1.018] | 1.052 [1.034, 1.058] | **1.011** | 0.997 [0.992, 1.001] | **0.997** [0.963, 1.014] | (targets not gained under `simulate`) |
| P12 resume running sums | cached | **W2** 1.589 [1.584, 1.678] | 0.992 [0.922, 1.020] | **1.589** | 1.039 [1.025, 1.042] | **1.021** [1.010, 1.042] | accepted |
| P13 skip unchanged rates | cached | **W3**, **W4** | 0.927 [0.916, 0.945] | 0.968 [0.952, 0.989] | **0.990** [0.981, 1.001] | **1.008** [0.984, 1.023] | **rejected**, reverted (W1 regression) |
| P14 one poll per event | cached | **W1** 1.095 [1.082, 1.132], **W3** | **1.095** | 1.042 [1.036, 1.051] | **1.030** [1.017, 1.069] | 1.022 [0.992, 1.034] | accepted |
| P15 flat save buffer | cached | **W5** 1.005 [0.998, 1.006], **W6** 1.003 [0.991, 1.006] | 1.026 [0.960, 1.313] | 0.999 | 0.999 | 1.013 [1.003, 1.022] | **rejected** as a speedup; kept as the memory-accounting fix |
| P16 reuse worker buffers | cached | **W7** 1.041 [1.003, 1.075] | 1.011 [1.002, 1.018] | 1.002 [0.974, 1.108] | 1.000 [0.996, 1.006] | 0.971 [0.957, 1.059] | accepted |
| R13 confirmation (P16 → P13 reverted) | cached | — | 1.008 [0.998, 1.017] | 1.010 [1.009, 1.014] | 1.002 [1.001, 1.003] | 1.031 [1.028, 1.039] | revert is neutral to positive |
| P17 FSP zero-row skip | fsp | **W8a** 1.245 [1.231, 1.254], **W8b** 1.140 [1.129, 1.190] | 0.987 [0.980, 1.001] | 0.988 [0.985, 0.989] | 1.002 [1.001, 1.008] | 0.995 [0.983, 0.999] | accepted |
| P18 FSP column gather, fused sum | fsp | **W8a** 1.218 [1.160, 1.229], **W8b** 1.723 [1.582, 1.741] | 0.996 [0.978, 1.001] | 1.007 [1.004, 1.013] | 0.998 [0.992, 1.006] | 1.008 [1.007, 1.017] | accepted |
| P19 tau PTRS | tau | **W9** 312 [296, 330] | 1.016 [1.010, 1.018] | 1.003 [0.997, 1.004] | 1.001 [0.985, 1.028] | 1.006 [1.004, 1.019] | accepted |
| P20 single serialization (process time) | cached | **W6** 0.989 [0.987, 1.001] | 0.994 [0.982, 1.004] | 0.996 [0.992, 1.001] | 0.990 [0.981, 1.000] | 1.007 [0.988, 1.014] | **rejected**, reverted |

Seven of ten candidates were accepted (P11, P12, P14, P16, P17, P18, P19);
P13 and P20 were rejected and reverted, and P15 gave no speedup but stays
because it fixes the memory reservation (see below). No iteration regressed a
frozen workload beyond the declared guard.

Notes:

- **P11** is accepted on W2 under `simulate-cached`. Under `simulate` the
  same change helps W1 (1.052) but not its targets: the full recomputation
  dominates there.
- **P12** gives the largest SSA gain of the series: on the sparse 66-reaction
  network a reaction changes few propensities, and the strict sums now resume
  at the lowest changed channel instead of starting over.
- **P13** compared each recomputed rate bitwise with the old one to resume
  later still. The extra comparison costs more than the rare skipped work
  (W1 0.927). P14–P16 were measured on a base that still contained P13; the
  R13 run measures the revert on top of P16 and finds it neutral to positive.
- **P15** replaced one vector per save point by one flat buffer per path. It
  is timing-neutral on W5/W6, so it is not counted as a speedup, but it fixes
  the review finding that the old layout allocated about 375 MiB against a
  241 MiB reservation for many save points; the reservation formula is now
  tested (`workspace_reservations_follow_the_published_formula`).
- **P17 and P18** are bitwise exact: unit tests compare the new products with
  the previous row product bit for bit, including exact zeros and subnormal
  terms. The product of their paired medians on W8b is 1.96; measured
  directly, the final source runs W8b 2.34 times and W8a 1.42 times faster
  than P10 (see the cumulative comparison below).
- **P19** passes all three declared gates: (a) the sampler test (200,000
  draws at each of the five means, DKW band at family alpha 0.001, sample
  means within five standard errors); (b) `check_tau.py` with its unchanged
  envelopes (its leap means stay at or below 16, so its draws are unchanged);
  (c) W9 speedup 312, far above 1.5: before P19 a leap with mean of order
  10^4 drew hundreds of chunks of the multiplication method, each needing
  about 17 uniforms. As supplementary
  evidence added after the declaration, `validation/compare_tau_samplers.py`
  compares 4,000 W9 paths under each sampler: two-sample Kolmogorov–Smirnov
  distance 0.0145 against a critical value of 0.0436 at alpha 0.001. The same
  comparison shows that tau leaping at epsilon 0.05 inflates the W9 variance by
  a factor of about 1.7 under both samplers, a property of the leap rule when
  epsilon exceeds the relative fluctuation size (here 0.01), not of the
  sampler; see the tutorial's tau section.
- **P20** was judged on process time, as amended before measurement. Keeping
  the bytes of the checking serialization did not save time on W6 (process
  0.989); the core timer, which now includes building the buffer, rose by
  about 2%. The candidate diff is kept as `P20-rejected-candidate.diff`.

## Cumulative comparison

The final source (P11, P12, P14–P19 and the review fixes) against the P10
binary and against the original baseline 5a07334 (M2, `simulate` only), with
the same harness and protocol. Paired median speedup [95% interval]; ns per
event from median core time.

| Workload | final vs P10, `simulate` | final vs P10, `simulate-cached` | final vs 5a07334, `simulate` | final `simulate-cached` vs final `simulate` |
|---|---|---|---|---|
| W1 two channels | 1.146 [0.952, 1.387] | 1.123 [1.085, 1.154] | 1.647 [1.429, 1.663] | 1.014 [0.999, 1.029] |
| W2 sparse 66 | 1.058 [1.054, 1.078] | 1.714 [1.679, 1.735] | 1.974 [1.956, 1.989] | 3.256 [3.241, 3.264] |
| W3 circuit 12 | 1.023 [0.994, 1.033] | 1.043 [1.043, 1.050] | 1.418 [1.404, 1.426] | 2.141 [2.047, 2.197] |
| W4 dense 66 | 1.063 [1.054, 1.136] | 1.037 [1.027, 1.057] | 1.369 [1.320, 1.392] | 1.176 [1.164, 1.275] |
| W5 circuit, 2001 saves | 1.030 [1.018, 1.063] | 1.052 [1.046, 1.053] | 1.424 [1.400, 1.438] | 2.136 [2.121, 2.206] |
| W6 W5 keeping 64 paths | 1.025 [1.013, 1.047] | 1.059 [1.051, 1.061] | 1.401 [1.396, 1.414] | 2.127 [2.115, 2.223] |
| W7 short paths | 1.038 [0.998, 1.083] | 0.988 [0.955, 1.024] | 19.86 [18.57, 20.82] | 0.974 [0.949, 1.022] |
| W8a FSP 401 states | 1.424 [1.413, 1.444] (fsp) | | | |
| W8b FSP 10,000 states | 2.340 [2.220, 2.396] (fsp) | | | |
| W9 tau, mean 10^4 | 335 [314, 343] (tau) | | | |

Final cost per event (median core time / events, one worker):

| Workload | 5a07334 `simulate` | final `simulate` | final `simulate-cached` |
|---|---|---|---|
| W1 two channels | 80.6 ns | 49.1 ns | 48.3 ns |
| W2 sparse 66 | 991.6 ns | 502.3 ns | 151.7 ns |
| W3 circuit 12 | 458.2 ns | 323.2 ns | 151.3 ns |
| W4 dense 66 | 902.7 ns | 659.5 ns | 553.3 ns |

All SSA comparisons reproduce the numerical digests exactly, including
against 5a07334. The FSP comparison is not digest-exact: the probabilities
are identical, but the review fixes changed the `work` count and the
omitted-weight and roundoff fields. Tau is approximate (`--approximate`). W7
paths last about 0.1 ms, so its intervals are wide; its 19.9 against
5a07334 is P06's single-worker path without a thread. Raw reports are in
`validation/campaign/p11-p20/cumulative/`.
