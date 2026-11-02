# Ten further Julia-gap iterations — implementation log

Date: 2 November 2026, UTC. Fresh GitHub base:
`add6cca` (merge PR #7). Protocol committed before candidate measurements:
[ten-iteration protocol](../decisions/julia-gap10-2026-11-02.md).

The attached log's prior optimizations are already merged. This campaign
investigates the installed JumpProcesses 10.0.1 source and adds C01–C10. The
benchmark host is AMD EPYC 9V74, Linux x86-64, one numerical worker. No hardware
instruction profiler is installed; source inspection and work counts support
mechanism explanations, while timings establish observed performance only.
Julia was rerun directly on this host, alongside independently built baseline
and final wheels. No cross-host projection is used.

Each entry below records a complete implementation/test/timing/decision loop.
The original sampler v2, precision policy and legacy stream identities are kept.

## Completed scope and result

All ten further iterations were implemented, tested, measured and decided in
sequence. Three changes remain: C05 lower-bound reuse, C06 exact custom-range
sharing and C09 compiled elementary kernels. Seven failed the rule declared
before measurement and were reverted. Their candidate patches and failed
development logs remain as evidence. The prior A01–A10 work is already in the
GitHub base and is not counted again.

| Iteration | Implemented experiment | Final decision |
|---|---|---|
| C01 | Inline the fast stop-poll path and move the slow path out of line | Revert |
| C02 | Poll every 64 random words, retaining per-word work budgets and other checkpoints | Revert |
| C03 | Batch changed tree ancestors in canonical summation order | Revert |
| C04 | Cache invalid-leaf metadata | Revert |
| C05 | Reuse an enabled lower mass-action corner as the minimum positive bound | Keep |
| C06 | Share exact custom-law integer ranges by bitwise law, regulator and feasible interval | Keep |
| C07 | Precompute power-of-two bucket ceilings | Revert |
| C08 | Scan SortingDirect rates contiguously with inverse positions | Revert |
| C09 | Compile common mass-action shapes into specialized enum payloads | Keep |
| C10 | Thin LTO and one codegen unit in native and binding release profiles | Revert |

The direct same-host geometric-mean CRN/fastest-Julia gap fell from **1.540×
to 1.441×** across seven frozen workloads. The aspirational 1.25× target was
missed. Julia remained faster on all seven workloads. W4's fastest prepared
method gained 1.350×; the fastest W5 and W7 choices were slightly slower.
One slower W4 DirectCR variant also regressed in a follow-up, documented below.
The retained implementation is a measured tradeoff, not a universal speedup.

The fixed candidate rule required target geometric gain ≥1.02, improvement in
at least four of five paired rounds, all-row geometric gain ≥0.995, and the
declared one-time check for any individual guard below 0.95. Each candidate
passed the complete Rust release suite and all 35 numerical identity rows
before timing. Every primary timing observation is retained.

## C01

Decision: **reverted**. Target geometric speedup 0.9880×; all-row geometric speedup 0.9880×. Target aggregate rounds: 0.9934×, 0.9748×, 0.9740×, 1.0051×, 0.9999×. All 35 numerical identities and the full release suite passed. Build/test time 16.9 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C01.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | yes | 0.9933× |
| W2/direct_cr | yes | 0.9668× |
| W2/sorting_direct | yes | 0.9666× |
| W3/legacy_cached | yes | 0.9858× |
| W3/sorting_direct | yes | 0.9944× |
| W4/cached_direct | yes | 0.9600× |
| W4/rssa_cr | yes | 1.0233× |
| W5/direct_cr | yes | 1.0003× |
| W5/rssa_cr | yes | 0.9572× |
| W6/rssa_cr | yes | 1.0700× |
| W7/rssa_cr | yes | 0.9566× |

C02 development check: the first attempt changed cancellation precedence at an
already exhausted sampling budget. The existing forced-tail test caught it.
Polling at budget exhaustion preserves that precedence; a retry-only zero-word
case now also checks the declared 64-word stop bound. Failed logs are retained.

## C02

Decision: **reverted**. Target geometric speedup 1.0060×; all-row geometric speedup 1.0207×. Target aggregate rounds: 0.9752×, 0.9681×, 0.9091×, 0.9907×, 1.0479×. All 35 numerical identities and the full release suite passed. Build/test time 15.9 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C02.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | guard | 1.1242× |
| W2/direct_cr | yes | 0.9963× |
| W2/sorting_direct | yes | 1.0245× |
| W3/legacy_cached | guard | 1.0557× |
| W3/sorting_direct | yes | 1.0042× |
| W4/cached_direct | yes | 0.9736× |
| W4/rssa_cr | yes | 1.0237× |
| W5/direct_cr | yes | 1.0066× |
| W5/rssa_cr | yes | 0.9785× |
| W6/rssa_cr | yes | 0.9967× |
| W7/rssa_cr | yes | 1.0525× |

## C03

Decision: **reverted**. Target geometric speedup 0.9898×; all-row geometric speedup 0.9891×. Target aggregate rounds: 1.0112×, 1.0040×, 0.9730×, 0.9960×, 1.0478×. All 35 numerical identities and the full release suite passed. Build/test time 17.8 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C03.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | guard | 0.9863× |
| W2/direct_cr | yes | 0.9750× |
| W2/sorting_direct | yes | 1.0001× |
| W3/legacy_cached | guard | 0.9948× |
| W3/sorting_direct | yes | 0.9512× |
| W4/cached_direct | guard | 0.9801× |
| W4/rssa_cr | yes | 1.0229× |
| W5/direct_cr | yes | 0.9561× |
| W5/rssa_cr | yes | 0.9574× |
| W6/rssa_cr | yes | 1.0257× |
| W7/rssa_cr | yes | 1.0346× |

## C04

Decision: **reverted**. Target geometric speedup 1.0138×; all-row geometric speedup 1.0167×. Target aggregate rounds: 1.0073×, 0.9961×, 1.0172×, 1.0162×, 0.9975×. All 35 numerical identities and the full release suite passed. Build/test time 17.5 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C04.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | guard | 1.0664× |
| W2/direct_cr | yes | 1.0204× |
| W2/sorting_direct | yes | 0.9916× |
| W3/legacy_cached | guard | 1.0245× |
| W3/sorting_direct | yes | 1.0063× |
| W4/cached_direct | guard | 0.9848× |
| W4/rssa_cr | yes | 1.0345× |
| W5/direct_cr | yes | 1.0592× |
| W5/rssa_cr | yes | 1.0192× |
| W6/rssa_cr | yes | 0.9917× |
| W7/rssa_cr | yes | 0.9892× |

## C05

Decision: **kept**. Target geometric speedup 1.0448×; all-row geometric speedup 1.0147×. Target aggregate rounds: 1.1002×, 1.1731×, 1.0205×, 0.9853×, 1.1098×. All 35 numerical identities and the full release suite passed. Build/test time 17.3 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C05.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | guard | 0.9818× |
| W2/direct_cr | guard | 1.0342× |
| W2/sorting_direct | guard | 1.0053× |
| W3/legacy_cached | guard | 0.9583× |
| W3/sorting_direct | guard | 1.0088× |
| W4/cached_direct | guard | 1.0757× |
| W4/rssa_cr | guard | 0.9552× |
| W5/direct_cr | guard | 1.0148× |
| W5/rssa_cr | yes | 0.9860× |
| W6/rssa_cr | yes | 1.0947× |
| W7/rssa_cr | yes | 1.0569× |

## C06

Decision: **kept**. Target geometric speedup 1.3568×; all-row geometric speedup 1.0312×. Target aggregate rounds: 1.3378×, 1.3568×, 1.3416×, 1.4127×, 1.4402×. All 35 numerical identities and the full release suite passed. Build/test time 17.2 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C06.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | guard | 1.0455× |
| W2/direct_cr | guard | 1.0231× |
| W2/sorting_direct | guard | 0.9890× |
| W3/legacy_cached | guard | 1.0251× |
| W3/sorting_direct | guard | 0.9906× |
| W4/cached_direct | guard | 0.9993× |
| W4/rssa_cr | yes | 1.3568× |
| W5/direct_cr | guard | 0.9674× |
| W5/rssa_cr | guard | 1.0647× |
| W6/rssa_cr | guard | 0.9597× |
| W7/rssa_cr | guard | 0.9741× |

## C07

Decision: **reverted**. Target geometric speedup 0.9927×; all-row geometric speedup 0.9924×. Target aggregate rounds: 0.9603×, 1.0080×, 1.0084×, 1.0119×, 1.0021×. All 35 numerical identities and the full release suite passed. Build/test time 16.4 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C07.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | guard | 0.9475× |
| W2/direct_cr | yes | 0.9872× |
| W2/sorting_direct | guard | 0.9960× |
| W3/legacy_cached | guard | 0.9885× |
| W3/sorting_direct | guard | 1.0517× |
| W4/cached_direct | guard | 0.9796× |
| W4/rssa_cr | yes | 1.0528× |
| W5/direct_cr | yes | 0.9364× |
| W5/rssa_cr | yes | 1.0205× |
| W6/rssa_cr | yes | 0.9712× |
| W7/rssa_cr | yes | 0.9923× |

## C08

Decision: **reverted**. Target geometric speedup 1.0156×; all-row geometric speedup 1.0085×. Target aggregate rounds: 1.0686×, 0.9779×, 1.0932×, 0.9646×, 0.9866×. All 35 numerical identities and the full release suite passed. Build/test time 16.9 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C08.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | guard | 1.0260× |
| W2/direct_cr | guard | 1.0849× |
| W2/sorting_direct | yes | 1.0593× |
| W3/legacy_cached | guard | 1.0035× |
| W3/sorting_direct | yes | 0.9737× |
| W4/cached_direct | guard | 1.0031× |
| W4/rssa_cr | guard | 0.9839× |
| W5/direct_cr | guard | 1.0177× |
| W5/rssa_cr | guard | 1.0016× |
| W6/rssa_cr | guard | 0.9883× |
| W7/rssa_cr | guard | 0.9586× |

## C09

Decision: **kept**. Target geometric speedup 1.0254×; all-row geometric speedup 1.0135×. Target aggregate rounds: 1.0360×, 1.1006×, 1.0210×, 0.9614×, 1.0348×. All 35 numerical identities and the full release suite passed. Build/test time 17.4 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C09.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | yes | 1.0577× |
| W2/direct_cr | yes | 1.0174× |
| W2/sorting_direct | guard | 0.9771× |
| W3/legacy_cached | guard | 0.9974× |
| W3/sorting_direct | guard | 1.0146× |
| W4/cached_direct | guard | 1.0193× |
| W4/rssa_cr | guard | 0.9893× |
| W5/direct_cr | yes | 1.1126× |
| W5/rssa_cr | yes | 0.9871× |
| W6/rssa_cr | yes | 1.0235× |
| W7/rssa_cr | yes | 0.9607× |

## C10

Decision: **reverted**. Target geometric speedup 0.9817×; all-row geometric speedup 0.9817×. Target aggregate rounds: 0.9639×, 0.9680×, 1.0246×, 0.9900×, 0.9680×. All 35 numerical identities and the full release suite passed. Build/test time 87.0 s, outside timing. Candidate diff and every observation are retained under `validation/julia-gap10/C10.*`.

| Row | Target | Paired median speedup |
|---|---|---:|
| W1/legacy_cached | yes | 1.0323× |
| W2/direct_cr | yes | 0.9507× |
| W2/sorting_direct | yes | 1.0190× |
| W3/legacy_cached | yes | 1.0224× |
| W3/sorting_direct | yes | 0.9714× |
| W4/cached_direct | yes | 0.9697× |
| W4/rssa_cr | yes | 1.0578× |
| W5/direct_cr | yes | 0.8865× |
| W5/rssa_cr | yes | 0.9623× |
| W6/rssa_cr | yes | 0.9764× |
| W7/rssa_cr | yes | 0.9623× |

## Final correctness gates

The retained C05/C06/C09 tree passed **188 release tests** on both Rust 1.88
and the 1.85 minimum. Core and separate binding clippy passed with warnings
denied. Formatting and Python lint passed. All 69 schema/CLI checks, ten variant
contract/concurrent/replay checks, Python and installed-package checks passed.
Legacy SSA (17 checks) and nonlinear CME (14 checks) passed in normal/optimized
Python and with cached execution; their byte-identical reports preserve legacy
numerical identities. The dependency gate passed with no new runtime dependency.

The full ten-variant exact-law campaign passed **210/210** fixed marginal checks
at 100,000 paths per fixture (7 million paths per pass), including wrong-law
power checks. Normal and Python -O reports are byte-identical. No threshold was
changed. The original per-random-word stop checks remain, since C02 was reverted.
C06 tests cover separate availability gates, consumed-regulator interval changes
and per-path cache reset. C09 compares elementary kernels with the checked
general evaluator at normal/subnormal coefficients, extreme integer counts and
overflow diagnostics. Neither optimization introduces a new sampler identity.

The C10 Thin-LTO build was rejected and release binaries/wheel were rebuilt with
the retained source and original release profile before final validation. The
GitHub-main baseline wheel is loaded in a separate Python process; a smoke solve
confirmed it uses sampler v2 and the intended baseline extension. The final
comparison includes all twelve prepared baseline methods, all final one-shot and
prepared variants, and the five pinned Julia algorithms.

## Direct baseline / final / Julia comparison

Actual versions: Julia 1.13.3, Catalyst 17.0.0 and JumpProcesses 10.0.1, from the
retained Manifest. Every row below uses the fastest of twelve prepared CRN
options or five Julia algorithms. Times are warm medians from five observations
per engine on the same host, one numerical worker, seed 4242. W1–W4 use 64 paths
to t=100; W5–W7 use 16 paths to t=20. Initial and final states are saved.

| Model | Baseline seconds | Final seconds | Baseline / final | Julia seconds | Final / Julia |
|---|---:|---:|---:|---:|---:|
| W1 | 0.07541 | 0.07084 | 1.065× | 0.05061 | 1.400× |
| W2 | 0.13661 | 0.13414 | 1.018× | 0.09802 | 1.369× |
| W3 | 0.36009 | 0.35007 | 1.029× | 0.29200 | 1.199× |
| W4 | 0.19371 | 0.14345 | 1.350× | 0.10719 | 1.338× |
| W5 | 0.08091 | 0.08233 | 0.983× | 0.04840 | 1.701× |
| W6 | 0.06830 | 0.06263 | 1.091× | 0.03883 | 1.613× |
| W7 | 0.02083 | 0.02112 | 0.986× | 0.01380 | 1.531× |

Ratios use unrounded medians; values above 1 in the last column favor Julia.
The full comparison retains **1,435 timing observations**, all 41 engine/API
rows per workload, ranges, replay hashes and preparation costs. No measurement
was removed. Small differences on this shared host should not be treated as
precise rankings. Builds and all heavy correctness checks finished before warm
timing. Julia's timing process took 67.68 seconds to reach readiness for all
seven models and five aggregators; this is not a per-model startup charge.

Independent-seed validation ran **60,928 paths and 6,946 comparisons**, all
passing, before timing. W1/W2 use exact Poisson references; other marginals use
the legacy cached reference with the unchanged conservative two-sample DKW
rule. Baseline/final numerical results match for every timed method. Finite
marginal checks do not establish arbitrary-model or rare-tail accuracy; the
known steep-Hill RSSACR precision limits remain documented.

The Rust prepared boundary includes simulation, result serialization and Python
decoding; baseline IPC is outside its timer. Julia includes serial solves,
saved-state checks and final reduction, with pristine problem copies outside
the timer. Julia retains learned SortingDirect order within each batch; CRN
resets per replicate for worker-independent replay. These are application
boundaries with different allocation/check pipelines. No multithread, GPU,
other-host or memory-scaling result is claimed.

## Investigated regressions

Six prepared rows had baseline/final ratios below 0.95 in the original complete
comparison. Five new alternating diagnostic pairs, two repeats per side, were
run separately. They did not overwrite primary results or alter acceptance.

| Prepared row | Original ratio | Diagnostic paired ratio |
|---|---:|---:|
| W3 RSSACR/inverse | 0.944× | 0.987× |
| W4 cached Direct/Ziggurat | 0.943× | 1.014× |
| W4 DirectCR/inverse | 0.907× | 0.990× |
| W4 DirectCR/Ziggurat | 0.940× | 0.878× |
| W5 cached Direct/inverse | 0.919× | 0.996× |
| W7 legacy cached | 0.916× | 0.977× |

Five drops did not reproduce at 5%; W4 DirectCR/Ziggurat did. Its diagnostic
0.878× throughput means approximately **13.9% more elapsed time**. RSSACR
remains the tested W4 choice and improves 1.350×, but a user specifically choosing
DirectCR/Ziggurat can incur a regression. This variant was outside the eleven
declared candidate rows; its result is retained rather than excluded.

Native W4 DirectCR was also measured at the base and each kept commit, three
samples per binary and wait sampler. Base/final medians were 1.36985/1.35862
seconds for inverse and 1.51546/1.37652 for Ziggurat. Native numerical digests
and the 1,226,555 events match exactly. These sequential samples did not
reproduce a final/base native slowdown. Generated-code layout, allocation and
host variation are possible explanations, not established causes. The Python
boundary regression remains unattributed without an instruction profiler.
All original and diagnostic observations remain in the package.

## Julia source findings and remaining gap

The pinned Julia source uses compact incrementally maintained group sums,
typed inlined mass-action loops, suitable endpoint brackets, Xoshiro/Ziggurat,
and learned SortingDirect order retained across solves. CRN already had the
previous campaign's joint member/acceptance sampler and envelope cache.
Source inspection therefore guided ten concrete new experiments rather than
assuming that changing language or adding inlining would guarantee a gain.

C05 and C06 remove repeated bound work while retaining exact integer ranges,
availability gates, per-path cache reset and errors. C09 compiles elementary
shapes once, retaining the original rounded operations and checked fallback.
CRN still pays for canonical sums, positive-channel resolution checks,
transactional population arithmetic, cooperative stops, work budgets and stable
worker replay. Copying Julia's unchecked integer updates or sum shortcuts would
change those guarantees; they were not part of this campaign.

No hardware profiler is installed on this host. Julia source-file hashes and
specific functions are recorded in `source-investigation.json` and the source
investigation document. Prior-host callgrind percentages are not current
attribution. No dependency, archive/request/result schema, RNG or sampler
identity changed. No remote was pushed and no package was published.

## Package contents and reproduction

The ZIP root contains one aggregate patch against
`add6cca0fb146e8a40418432e963219514086024`, this implementation log, an apply
README, original license/notices and a per-file SHA256 manifest. The `repo/`
folder contains changed files, all raw evidence, rejected candidate diffs,
full comparison, Julia source analysis and reproduction instructions. It is a
changed-files snapshot, not a full checkout. Apply only the root aggregate
patch; nested CXX patches are rejected-experiment evidence.

In the patched repository, start with `validation/julia-gap10/README.md`.
`validation/validate_gap10.py` records the full executed gate list;
`validation/compare_julia_gap10.py` performs independent validation followed by
the same-host timing run; `validation/report_gap10.py` renders without retiming;
`validation/diagnose_gap10_regressions.py` retains diagnostic follow-up.
Campaign-specific local tool and baseline paths should be adjusted on another
host. Rebuild and force reinstall the wheel because the development version
number is unchanged. Original MIT attribution and third-party notices remain.
