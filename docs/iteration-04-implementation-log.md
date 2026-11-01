# CRN-craft iteration 4 implementation log

Date: 1 November 2026. Milestone: **M2 — Direct SSA, ensembles and archives**.
Version: `0.1.0-dev.4`.

Repository: https://github.com/dietmarwo/CRN-craft

Baseline: **`bb2b9aed6129734bbe5856a709e0f2d1e32ca198`**, main after PR #3.

## Outcome and GitHub review

Fetched the actual GitHub main and read the attached iteration-3 review,
repository instructions, current plan, source and CI. Main includes M1b and
the review's primitive-ray, initial-ray-cap and converging-linkage regressions.
Its 17-mutant gate and the five-targeted-model structure oracle are preserved.
The untouched baseline passed 100 release tests and 58 contract checks here.
A final remote-ref check still reported `bb2b9ae`.

Implemented the next planned milestone: executable Direct SSA with the existing
mass-action, Hill and Michaelis–Menten rates; pinned random streams and owned
samplers; deterministic parallel ensembles and mean/sample variance/mean SE;
event, population, deadline, serialized-output and workspace admission limits;
cooperative cancellation; exact input snapshots and verifiable result archives.
The delivery is a Git mail patch ZIP plus this standalone Markdown log. No
push, tag, package publication or deployment was performed.

## Changes

| Area | Implementation and behavior |
|---|---|
| Sampling | `sampling.rs`: pinned ChaCha8, SHA-256 root-seed key, replicate stream IDs, explicit uniform endpoints, libm wait inversion and strict ordered selection |
| Scheduler | `simulation.rs`: Direct recomputation, transactional count updates, absorption, event/save equality and pending-event reuse across the grid |
| Ensembles | Bounded worker queues, ascending-ID collection and Welford reduction; all IDs accounted for, including unstarted paths |
| Limits | Preflight workspace/evidence reservation, live event/population/deadline/cancel checks, actual serialized-byte counting and bounded archive writes |
| Results | Result 0.2, per-replicate stop/work records, first-ID retained path prefixes, statistics only for complete ensembles |
| Archives | `archive.rs`: exact input bytes, four fixed payloads, SHA-256 manifest, marker written last, overwrite refusal and independent verification entry point |
| CLI | `simulate MODEL RUN NEW_DIRECTORY [--cancel-file PATH]` and `verify DIRECTORY`; coded JSON results and honest exit semantics |
| Contracts | Discovery 0.4; archive/archive-check 0.1; result-v1 preserves the reserved 0.1 schema; eleven generated snapshots |
| Validation | New RNG, trajectory, analytic and archive checks; nine new source mutants; expanded CI and updated plan/agent documentation |

Short Rust paths are under `crates/crn-craft/src/`. The public network APIs
retain their state checks. The scheduler uses a private known-valid-state rate
path to avoid repeating O(species) validation for each channel. Rate formulas,
conservation/cone/linkage algorithms, exporter and reviewed regressions are
otherwise unchanged.

RNG dependencies are pinned `rand_chacha 0.3.1` and `rand_core 0.6.4`, with
default features disabled. The lockfile also adds ppv-lite86 0.2.21 and
zerocopy/zerocopy-derive 0.8.63. Their permissive licenses and MSRV metadata are
recorded in `validation/iteration-04/rng-dependencies.json`; all three Rust
toolchains passed. No native solver/build dependency or Python runtime was
added to CRN-craft. Project code continues to forbid unsafe Rust; this does not
claim that dependencies contain no unsafe code.

I inspected current `optimization-jobs` storage and artifact handling at
`cc6b03d594f73d1f3b38f088e7e3bd8663814d75`. M2 uses the same principles of
checksummed evidence and synchronized writes, without importing the job layer
or copying its implementation. The shared runner integration remains M6.

## Numerical and execution decisions

The sampler contract is `crn-craft.chacha8-direct/0.1`. Each replicate gets the
same stream regardless of thread allowance, output density or model metadata.
Wait uniforms reject zero and exclude one; channel uniforms include zero and
exclude one. Waits use `-log(U)/total`. An event exactly at a save/horizon time
is included, and pending events are not redrawn at intermediate save points.

Positive rate underflow remains CRN-E042, as the review requested. Total-rate
overflow, lost positive increments, probabilities below the 53-bit selection
resolution and rounded-to-total thresholds fail with CRN-E061. Nonpositive,
nonfinite or stagnant event times fail with CRN-E062. No silent zeroing,
clipping or last-channel fallback was introduced. Null reactions count as
events; population safety stops do not reflect or apply a rejected event.

Every requested replicate has a record. Incomplete ensembles have no
unconditional statistics, because completion may depend on state. Fully
completed results use sample variance with denominator n-1 and mean SE
sqrt(variance/n); a single replicate has no estimated SE. Variance uncertainty
is not invented. Exact i64 counts can still lose precision when converted to
f64 summaries; JSON consumers must preserve 64-bit integers to retain counts.

The engine admits at most 1,000,000 requested replicates, 64 workers and a
256 MiB conservative workspace reservation; feasible combinations may be much
smaller. Encoded archive contents are also capped at 256 MiB and by the run's
output budget. These are content/allocation checks, not an OS RSS or real-time
guarantee. Deadline/cancellation checks cover trajectory/rate/save work;
parsing, compilation and final archive I/O are outside the wall-time budget.

Output counting occurs without building an unbounded serialized buffer. When
the numerical payload does not fit, trajectories/statistics are removed,
CRN-E064 and `output_limit` are recorded, and complete replicate accounting is
retained if it fits. Archive writing enforces aggregate content bytes again.
Forced process death may leave only input snapshots or partial files and no
valid marker. M2 has no resumable or per-event recovery journal.

## Executed validation

Targets were recorded in `docs/decisions/iteration-04-validation.md` before
the first simulations. No statistical threshold was relaxed.

| Gate | Result |
|---|---|
| Rust 1.85.0, 1.88.0 and stable 1.99.1 | **126 release tests on each**, fmt and clippy with warnings denied |
| Contract/schema/CLI checker | **61 passed**, normal and `-O` identical |
| Preserved M1a rate/network oracle | **218 passed**, normal and `-O` identical |
| Reviewed M1b cone/graph oracle | **309 passed**, normal and `-O` identical |
| GillesPy2 1.8.3 exporter checks | **6 passed**, 15,717 propensity comparisons; normal and `-O` identical |
| New SSA acceptance | **17 grouped checks**, normal and `-O` identical |
| RNG oracle | **360 words across nine streams**, including block and high stream-ID boundaries |
| Source mutation gate | Passing unmutated baseline; **26/26 compiled mutants caught** |
| Dependency policy, Ruff | Passed; dependency and RNG reports identical under `-O` |

Rust breakdown: 100 reviewed baseline tests, three sampler/boundary unit tests,
14 SSA integration tests and nine archive tests. New cases cover sampler
endpoints, exact event/save equality, absorption, null events, state/overflow
failures, all requested IDs, mixed completion without conditional statistics,
budget boundaries, corruption, duplicate keys/paths, missing markers, symlinks,
overwrite refusal, false accounting/statistics and oversized-file rejection
before payload reads. The last oversized-file regression was added after the
mutation campaign; source implementation and all mutation cases were unchanged.

### Predeclared analytic campaign

100,000 independent replicates each; times 0, 0.25, 1 and 2 seconds. Root seeds
2026110401, 2026110402 and 2026110403. The nine positive-time marginal CDFs use
DKW with Bonferroni family alpha 0.001: **epsilon = 0.006999331052635781**.

| Fixture | Analytic law | Executed events | Largest CDF error |
|---|---|---:|---:|
| Pure decay, X(0)=20, death 0.7 | Binomial(20, exp(-0.7t)) | 1,506,415 | 0.003334014071384228 |
| Immigration/death, X(0)=0, birth 3, death 0.8 | Poisson((3/0.8)(1-exp(-0.8t))) | 898,423 | 0.002944659262227145 |
| Same rates, X(0)=12 | Surviving-initial binomial + independent immigration Poisson | 1,858,609 | 0.0022524040621131336 |

All **300,000 paths** completed; all nine tests passed. The family executed
4,263,447 events. Independent integer reductions reproduced mean, sample
variance and mean SE within the declared numerical tolerance. The checker
rejects an all-zero solver for all three models. Histograms, fixture hashes,
result hashes and every CDF row are retained in `ssa.json`.

A separate Python Direct-path reference uses scalar ChaCha8, integer binomial
factors and its own event-history lookup. It matched 64 two-species paths and
event counts at 1, 2 and 8 threads, including homodimer conversion and non-unit
Omega. Numerical JSON was identical across those thread allowances after
excluding the raw run hash. Enzyme paths preserve E+ES=10 and S+ES+P=50.

Actual CLI cancellation interrupted active event generation. Final normal and
`-O` measurements were about **5.3 ms and 5.9 ms**, below the predeclared local
one-second target. Their event counts differ, as expected for timed stops.
Timings are stored separately from reproducible numerical reports. Live
wall-time stopping and valid incomplete archives also passed.

### Mutations and retained development issues

The original 17 mutations remain. Nine additions alter wait scaling,
selection equality, stream separation, save-before-equal-event behavior,
event-cap comparison, variance denominator, conditional statistics, payload
digest checking and completion digest checking. Each compiled and caused
named test failures after a passing baseline; compile errors never count.

The initial Python lint pass found an implicit zip strictness argument and a
loop-variable style issue in the independent oracle; both were corrected.
No numerical/structural fixture failed. Later multi-toolchain builds in the
workspace encountered damaged output files: a zero-byte/non-executable test
binary and a zero-length object mapping error. Those logs are retained under
`development-environment-*`. Clean temporary target directories passed all
gates without changing solver behavior or numerical thresholds. The exact
delivery patch is rebuilt independently again after freezing it.

## Compatibility, limitations and next iterations

Eight legacy schema snapshots (including reserved result as result-v1) and
all original model/run fixture bytes remain unchanged. Current result changes
from reserved 0.1 to emitted 0.2, with migration and new schemas documented.
Historical check/analysis operation flags remain unchanged; discovery 0.4 is
the current capability catalog. Archive write completion is explicitly
separate from simulation completion.

Validation ran on Linux x86_64, Python 3.12.14, nine visible/affinitized logical
CPUs and an eight-CPU quota. Same-platform reproducibility was tested; other
platforms, physical-core scaling and hard real-time cancellation were not.
Archive checks establish integrity and internal consistency, not signatures
or an adversarial filesystem security boundary. With only some paths retained,
full summary recomputation is unavailable. Analytic agreement and mutant
detection are evidence of correctness, not a proof for arbitrary models.

**Next is M3:** reachable-state FSP, explicit escaping-mass sink, stable sparse
uniformization, separate leakage/truncation/numerical accounting and selected
DSMTS/independent validation. **M4a** then supplies pinned same-host quality/cost
comparisons. **M5** adds bindings and the circuit application. M4b advanced
aggregators remain bottleneck-driven; M6–M10 cover integration, optional ODE/
sensitivity/import/approximation work and release acceptance.

There is now a working Direct baseline, but **no defensible competitive
performance-gap estimate yet**. The campaign above is correctness evidence;
it is not a controlled throughput comparison with GillesPy2, rebop or Julia.

## Delivery verification

The ZIP includes one mail patch against the exact main baseline, application
instructions, this log, checksums and fresh-checkout verification evidence.
The post-commit appendix records application, resulting-tree identity and
rebuild/acceptance results for the frozen patch. Those hashes live outside the
source commit to avoid self-referential patch metadata.
