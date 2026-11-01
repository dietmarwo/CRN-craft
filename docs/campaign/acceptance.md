# Remaining-milestone and performance campaign

Baseline: GitHub main `5a0733444aa46c8cea2a5b61bd8601daf5890d79`, fetched 2026-11-01.

This document is committed before new numerical experiments. Existing strict inputs,
Direct random-word streams, result/archive versions and incomplete-ensemble rules
remain compatibility gates. New capabilities use separate versioned requests/results;
no existing contract silently acquires a new meaning.

## Functional iterations, in order

1. M3: reachable sparse FSP with an escaping sink, time-subdivided uniformization,
   explicit series-tail and roundoff accounting, work/state/memory/output/cancellation
   limits; analytic per-state errors <= 1e-12 on small linear/conserved fixtures;
   independent matrix-exponential comparison <= 2e-12. Retain failures.
2. M4a: freeze Direct workload and benchmark against pinned GillesPy2 and rebop
   where installable. Same host, same dynamic model, seed policy, work and outputs;
   five paired timings with raw data, medians and ranges. Record unavailable Julia
   and physical-core scaling measurements explicitly. No speed claim from unlike work.
3. M5: isolated PyO3 binding and installable Python package; CLI/binding equality;
   observable autocorrelation/spectrum and censored first-passage summaries; labelled
   synthetic three-gene study with disjoint development/validation seeds. A small
   study does not establish the original universal 95%/99% behavior targets.
4. M4b: dependency-cached Direct selected from the review's measured full propensity
   bottleneck. Identical summation and selection order and pinned trajectory identity.
5. M6: inspect then-current optimization-jobs/fcmaes APIs, use their existing protocol,
   fixed seed-set fit-moments/circuit objectives, explicit validation calls, bounded
   threads and failed/incomplete objective rejection. Run an end-to-end fitting study.
6. M7: isolated pinned optional diffsol path, deterministic count-kinetic convention,
   forward sensitivities and global-to-local fitting if the actual dependency/toolchain
   permits. Record measured dependency limitations rather than inventing a BDF solver.
7. M8: provenance-preserving explicit SBML and expanded BioNetGen .net subsets;
   strict rejection of unsupported semantics. Round-trip fixtures and rejection tests.
8. M9: demand-driven bounded non-negative tau leaping only after an exact throughput
   bottleneck is measured; declared approximation and bias study at epsilon 0.1,
   0.05, 0.025, 10000 independent paths each, family CDF DKW alpha=0.001, plus
   analytic mean/bias checks. A failure is retained, not reclassified as exact.
9. M10: build wheel/sdist, install in a fresh environment, run package acceptance,
   preserve owner-controlled publishing. No push, tag or upload to package indexes.

## Ten performance iterations after functional work

Each has a hypothesis, code or rejected candidate, correctness check and paired
before/after measurement. A failed candidate counts as an experiment, not a speedup.
Freeze model generation, event cap, seed=4242, end=100, 64 replicates, and sparse
save grid [0,100]. Models: 2-channel immigration/death, 66-channel independent
birth/death, 12-channel three-gene circuit, dense-coupling control. Use release
builds, 1 worker, 7 paired repetitions (3 for exploratory runs). Exact event totals
and numerical digests must match before timing comparisons. Report process and
core elapsed separately. This shared host may only support limited core scaling.

P01 specialize order-zero/one mass-action; P02 order-two fast path with range fallback;
P03 avoid full state revalidation on trusted scheduler firings; P04 no redundant
re-evaluation of selected propensity; P05 reuse path buffers/allocation; P06 bypass
thread/channel startup for one worker; P07 batch bounded worker queues; P08 faster
strict ordered sum/min check; P09 precompute Hill constants for safe common cases;
P10 FSP sparse-product/output reuse. Adopt only numerically valid candidates, retain
negative timing rows. Broader work may change a candidate's implementation, not its
predeclared correctness threshold.

Existing regression/mutation suites, fmt, clippy -D warnings, release --locked,
normal/-O checker equivalence and dependency policy are final gates. Optional
integration dependencies stay outside the default workspace and have separate pins.
