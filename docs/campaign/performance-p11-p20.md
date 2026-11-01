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
