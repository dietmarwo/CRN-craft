# Exact SSA acceleration implementation log

GitHub base: baed0e4590f93e4c1dcc77cd7f039a2604a1594c (fetched 2 November 2026).
The unmerged comparison at dbb6809 is also included. Protocol: ba76fe2,
`docs/decisions/ssa-acceleration-2026-11-02.md`, committed before experiments.

## Step 1 — profiling and reusable Python preparation

Added a diagnostic Rust component profiler and a Python orchestration script.
No hardware profiler was installed; component microbenchmarks are labelled and
are not claimed to be an additive breakdown of event-loop time.

Exposed the existing immutable Rust `Prepared` object through `crn.prepare` and
`PreparedSimulation.simulate`. Each execution creates its own worker/RNG state,
so repeated calls replay the same request and concurrent reuse is independent.
Both properties passed on all four frozen models. A new request is needed to
change a seed or other input; no silent mutation or seed advancement is implied.

The whole-batch prepared gains are small: W2 0.200→0.197 s, W3 0.669→0.667 s,
W4 0.867→0.856 s. Preparation is not the principal explanation of the previous
Julia gap. Component results and all five Python timings are retained in
`validation/ssa-acceleration/step1-profile.json`. The engine itself is unchanged.

Release build/tests and wheel build passed. Python lint initially found a
loop-variable closure warning in the concurrent-reuse checker; explicitly
binding the shared handle fixed it before execution. Ruff then passed.

## Subsequent steps

Results are appended as each implementation/test/measurement loop completes.

## Step 2 — DirectCR

Implemented an original composition–rejection variant: power-of-two envelope
mass chooses a bucket, unbiased integer sampling chooses a member, and a global
rejection step gives accepted channel probabilities proportional to actual
rates. Waiting times use the actual total. Dependency updates and a balanced
sum tree avoid repeated prefix scans and accumulated subtraction drift.

Added separately versioned `ssa-run/0.1` and `ssa-result/0.1` contracts, native
and Python `ssa` commands, and immutable prepared Python SSA handles. The common
scheduler/reducer is shared with legacy execution. New results expose random
word and propensity-work counts; exhaustion is `sampling_limit` and withholds
statistics. Every random-word rejection loop polls cancellation and its budget.

Release tests, clippy with warnings denied and the Python wheel build passed.
A million categorical samples passed the predeclared DKW bound after bucket
activation/removal/migration. Direct/cached inverse paths match legacy outputs;
new DirectCR paths replay across 1/3 workers and repeated prepared solves.
All four workload smoke tests completed and prepared/one-shot outputs matched.
Two initial full-suite runs failed only because discovery contract assertions
still expected the old command list/milestone; those assertions were updated and
both failed logs retained. Full statistical and performance gates follow after
all five methods are present; smoke timings are not comparative conclusions.

## Step 3 — SortingDirect and Ziggurat

Added the published Sorting Direct move-one-position heuristic, reset for each
replicate and coupled to dependency updates and checked balanced totals. Added
a separately identified exponential Ziggurat sampler using Julia 1.13.3's MIT
tables, with the full upstream license in `THIRD_PARTY_NOTICES.md`. Its bounded
iterative rejection/tail loops use CRN-craft's existing per-path ChaCha streams;
it neither changes legacy sampling nor promises Julia-identical random paths.

Release tests and clippy passed. All eight available algorithm/wait combinations
replay across threads. One million exponential draws passed all 90 declared CDF
checks (mean 0.998291, second moment 1.993573; descriptive only). Forced zero,
tail and repeated rejection inputs respect the work limit and cancellation.
Prepared and one-shot results agree on all four development workloads. These
single-observation smoke timings show that tree maintenance can cost more than
prefix sums on the dense custom-rate model; no universal speedup is claimed.
The final five-round comparison will retain these slower variants.

## Step 4 — bounded RSSACR

Implemented exact thinning with composition/rejection over upper propensities,
cheap acceptance below certified lower propensities, and reverse-reader updates
only when a changed population exits its bracket. Every rejected candidate
advances the clock; only accepted reactions increment the event counter.

Bounds cover integer boxes of the actual gated evaluator. Elementary normal
mass action uses proved monotone arithmetic; custom one-variable rates enumerate
at most 129 feasible integer values. Higher-order/range-risk mass action uses
point bounds. A neighboring numerical error can shrink and pin read coordinates
when the current propensity is valid. Gated repression, distinct availability
species, zero activation and point fallback/refresh passed exhaustive checks.
Retained bucket capacity and combined RNG/evaluation accounting are bounded.

All 210 predeclared marginal comparisons passed at N=100,000 per fixture/variant;
maximum CDF error 0.00491534 < 0.00804613. Every wrong homodimer convention was
rejected. The complete Python -O repeat produced a byte-identical report after
the final capacity/accounting hardening. No thresholds were adjusted. Release
tests, clippy, 69 schema/CLI checks, ten native/Python/prepared/concurrent variant
checks, legacy Python checks and the dependency gate passed. Normal/-O contract
and dependency reports agree. Old contract snapshots remain unchanged.

A new contract test initially expected a zero event cap to execute; the existing
positive-cap validator correctly rejected it. The test now uses one event to
exercise incomplete execution. Ruff caught unused benchmark loop names and a
concurrent-test closure binding; these were fixed. Failed logs are retained.
Development smoke timings are not the final performance comparison.

### Step 4 review follow-up: absorption and packaging

Review found that positive RSSA envelopes can surround a state where every
actual propensity is zero. Added an explicit absorption check after state
changes when no lower bound guarantees a positive channel. It stops at the
first positive rate, preserves the fast path when a positive lower bound exists,
and completes truly absorbing states without random words. A zero-budget gated
pair-reaction regression covers this case for both wait samplers. Earlier law
reports remain in place; the full fixed gate is rerun on the final engine.

The initially started seven-workload validation was intentionally interrupted
after W3 and is retained under `superseded-validation/`. The final comparison
will require hashes of the actual Rust/Python source as well as the harness.
The Python binding workspace was formatted and passed its own clippy gate.
The upstream Ziggurat notice is included both in the core crate and in the built
wheel's license directory; this was checked in the wheel ZIP.

### Final correctness gates

The final absorption-aware engine passed all 210 fixed marginal checks again,
and its normal/-O reports are byte-identical (7,000,000 paths per pass). Added a
cooperative external-stop regression while RSSACR is actively rejecting candidate
events; both wait samplers cancel with work already performed and no statistics.
Its first request mistakenly asked for sample variance from one path; the strict
validator rejected it. The corrected request asks only for the mean. No engine
validation rule changed. The failed test log remains available.

Final release suites pass **181 tests on Rust 1.88 and Rust 1.85.0**. Core and
separate Python binding clippy pass with warnings denied. Final schema checks
pass 69 cases; all ten variant contract/replay/concurrent cases pass in normal
and optimized Python. Legacy SSA (17 checks), nonlinear CME (14 checks), Direct/
cached equality, Python interfaces, installed-package checks and dependency
checks pass. The installed-package checker needed its expected discovery version
updated from 0.5 to 0.6; the first failure is retained. CI now runs the new
contract and 210-comparison gates, including optimized-Python equality.

## Step 5 — repeated Julia comparison and larger connected networks

Frozen W1–W4 were retained, and W5–W7 add a 256-species/510-reaction diffusion
chain, a 32-species/992-reaction conversion network, and a 96-species/126-reaction
compartment network. All 17 engines passed the declared 512-path-per-model
comparison: 60,928 paths, 6,946 marginal checks, unchanged family alpha 0.01.
All stoichiometric/propensity probes and conservation checks passed. Raw final
samples are compressed losslessly with hashes. SciPy emitted a warning when its
unused exact KS p-value calculation fell back to an asymptotic calculation;
acceptance uses only the unchanged KS distance and the predeclared DKW bound.
The warning is retained. Both completed Julia stderr logs are empty.

After all builds and heavy checks finished, a fresh Julia process performed
five shuffled/rotated timing rounds on all seven workloads. The report retains
**1,015 observations** across 12 Rust variants at two API boundaries and five
Julia algorithms. Every batch reproduces its own numerical digest. No outliers
were removed. Source and sample-file hashes match the delivered files.

| Workload | Fastest CRN prepared | Speedup over legacy cached | Fastest CRN / fastest Julia |
|---|---|---:|---:|
| W1 | Legacy cached | 1.00× | 1.54× |
| W2 | SortingDirect/Ziggurat | 1.19× | 1.89× |
| W3 | Legacy cached | 1.00× | 2.39× |
| W4 | RSSACR/Ziggurat | 4.04× | 2.24× |
| W5 | RSSACR/Ziggurat | 2.39× | 2.20× |
| W6 | RSSACR/Ziggurat | 6.71× | 2.21× |
| W7 | RSSACR/Ziggurat | 1.63× | 1.79× |

These ratios compare medians on this host, not universal algorithm guarantees.
Julia remains faster on every measured network. RSSACR greatly reduces rate
work on W4/W6, while tree maintenance/rejection/accounting can outweigh savings
on smaller networks. The original cached implementation remains valuable.
Preparation mattered little for W1–W4 but W6's best new method improves from
0.1141 s one-shot to 0.0874 s prepared. The complete comparison retains every
slower variant, timing ranges, startup phases, work counts and scope limitations.

The comparison Markdown and its scientific figure are generated directly from
the retained measurements. The plot was visually checked. Earlier comparison
text and raw evidence are preserved. Reproduction instructions include pinned
Julia environments, Python dependencies, wheel installation and source-hash
checks; CI covers the new statistical and contract gates.

## Packaging

The deliverable contains an aggregate patch from GitHub base baed0e4 and an
alternative from the previous comparison base dbb6809, plus changed source files,
raw evidence, the comparison, this log, attribution and application instructions.
Only one patch should be applied. The packaging script verifies both patches in
clean detached worktrees, requires their resulting Git trees to equal the final
committed tree, computes per-file SHA256 hashes and checks ZIP CRC integrity.
The generated patch-verification report records the actual outcomes. No remote
was pushed and no package was published.
