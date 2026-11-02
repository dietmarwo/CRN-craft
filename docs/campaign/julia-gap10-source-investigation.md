# What the pinned Julia code does differently

Inspected locally on 2 November 2026: Julia 1.13.3, Catalyst 17.0.0 and
JumpProcesses 10.0.1, retained by the Julia Manifest. These are source-level
observations. No hardware profiler is installed on this host; do not interpret
them as percentages of current runtime.

| Source | Mechanism | Relevance to this campaign |
|---|---|---|
| `aggregators/prioritytable.jl`, `update!` / `sample` | In-place group and overall sums; descending group scan; one floating draw supplies member and acceptance via integer/fractional parts | CRN already caches envelope prefixes and has sampler v2 joint draws. Remaining per-leaf checked tree work motivates C03/C04; C07 removes repeated ceiling construction |
| `aggregators/rssacr.jl`, `generate_jumps!` | Rejection candidates accumulate unit exponentials, divide by the envelope sum once; no per-random-word cooperative stop callback | CRN must still observe save times and bounded stops between candidates. C01/C02 reduce stop overhead while keeping work limits |
| `aggregators/bracketing.jl:82–83` | Mass-action brackets evaluate lower/upper endpoints | C05 reuses an already computed enabled lower rate. CRN custom bounds enumerate integer states to cover gated/nonmonotone floating evaluators; C06 shares identical successful ranges instead of weakening their proof |
| `massaction_rates.jl:6–18` | Inlined, typed mass-action loop; integer reactant product then scaled floating coefficient | C09 specializes common elementary shapes while preserving CRN's existing rounded multiplication order, overflow/underflow diagnostics and combinatorial divisor |
| `massaction_rates.jl`, `executerx!` | `@inbounds` direct integer additions to sparse state | CRN retains transactional checked count updates and configured safety limits; no unchecked integer update is proposed |
| `aggregators/sortingdirect.jl:initialize!` / `execute_jumps!` | Learned order persists across solves; fired jump moves one place forward | CRN resets per replicate for thread-independent replay. C08 keeps that reset but makes rate scans contiguous |
| Julia `Random/src/normal.jl`, `randexp` | Ziggurat fast path; Xoshiro RNG in the benchmark | Already represented by CRN's versioned Ziggurat option; ChaCha8 and its stream contract are retained |

Julia's `PriorityTable` also has absolute-epsilon shortcuts and does not rebuild
our balanced canonical sum or enforce our positive-channel resolution policy.
Copying it wholesale would change the numerical/execution contract. The same
applies to removing stop budgets, unchecked population arithmetic, or replacing
the stable Hill evaluator with an unguarded power formula. Those are not the
experiments in this campaign.

The source suggests efficient specialized loops and compact mutable state, not
an intrinsic language guarantee. C10 tests cross-crate Thin LTO so the release
compiler can see more of the calling context. It records compilation overhead
separately. Individual source-file hashes are captured with the final evidence.

Reference documentation: [JumpProcesses aggregators and bounds](https://docs.sciml.ai/JumpProcesses/stable/jump_types/)
and [problem reuse / SortingDirect replay](https://docs.sciml.ai/JumpProcesses/stable/tutorials/ensembles_and_problem_reuse/).
The original code is MIT licensed; this campaign implements its own structures
and keeps the existing Julia-derived exponential-table notice.

## What the ten experiments established

C05 (enabled lower-corner reuse), C06 (shared exact custom ranges) and C09
(elementary kernels) met the fixed paired timing rule. The other seven changes
were reverted, including polling changes and Thin LTO. Thus this experiment does
not support treating more inlining, batching, contiguous scans or LTO as automatic
speedups. Compact data and less repeated work are useful hypotheses that still
need measurement on the actual workloads. Runtime stop polling, checked population
arithmetic and the canonical precision policy are retained.

The final same-host rerun reduced the geometric-mean gap from 1.540× to 1.441×;
Julia remained faster on all seven workloads. W4's best RSSACR choice gained
1.350×, but DirectCR/Ziggurat regressed in the Python prepared boundary
(diagnostic ratio 0.878×). Native localization samples did not reproduce a
final/base slowdown, despite identical paths and event counts. Layout,
allocation and host variability are hypotheses; without a profiler, attribution
is unresolved. See the full comparison and retained diagnostic observations.
