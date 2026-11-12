# Reusable parameters, generated kernels and rebop: full optimization

> **Measured with the 0.6.0 network engine (note of 12 November 2026).**
> The CRN-craft routes on this page (`ParameterTemplate`, generated
> propensity kernels, the `fast` policy) belong to the network engine of
> releases 0.2.0 to 0.6.0, which release 0.7.0 replaced with one engine. The
> full optimization run was not repeated, and the tables below are left as
> measured. In 0.7.0 a parameter search binds values on a prepared run
> (`bind_parameters`, see [parameter sweeps](reusable-parameters.md)), and
> the benchmark worker has the routes `fresh` and `reuse`. Their worker
> calls were timed against the 0.6.0 workers in one session
> (`validation/merge/parameter-workers-m4d.json`, declared and reported in
> `validation/merge/RESULTS.md`, "M4d (law gate of the run/1.0 default,
> worker speed rerun)"), with a declared target of 1.25 times the 0.6.0
> time. New `reuse` against old `reuse` took 1.10 (birth_death) and 1.05
> (sparse33) times as long, and **1.55 times on binding12, which failed the
> target**: there `bind_parameters` prepares the whole model again on each
> call. Both routes of each engine returned identical moments.

> **Measured before the 0.2.0 release.** These measurements come from the
> import commit `85a085b` of 2 November 2026 and have not been repeated on
> release 0.2.0 (4 November 2026). Review fixes since the import
> (`fdfd3d8`) changed the fallback that generated kernels use for rates they
> do not specialize (other laws than mass action, reactions above second
> order) and for values that are not normal floats; the three mass-action
> models here have no reaction above second order. `fast` also became the
> default of `ssa` requests (`crn-craft.ssa-run/0.3`), and the Ziggurat
> waiting-time sampler was refactored so the rules engine can share it.
> Results are unchanged; the timings on 0.2.0 may differ.

Measured 2 November 2026, America/La_Paz. Base GitHub main `c16fcc5`, including
PR #8's optional `method.fast`. Actual engine versions: CRN-craft 0.2.0.dev1,
rebop 0.10.0, fcmaes-rust 0.2.1.
This is a same-host experiment; it does not project the previous Julia timings.

All **54 complete optimization runs** finished: three fixed
models, six execution routes and three optimizer seeds. There were **6,912
actual objective evaluations**. All nine three-route CRN optimizer trace
comparisons were identical, point for point and loss for loss. Compiled code
reused each topology across all candidates; coefficients remained runtime data.

The reusable interface reduced total optimizer time by about 7% on sparse33
and 1.73× on binding12 (per seed 1.69–1.77×); birth/death was essentially flat.
Three seeds with one repetition each cannot resolve smaller effects: on
birth/death, where fresh and reusable CRN simulate identical paths, their ratio
still ranged 0.85–1.16× across seeds. Generated propensities gave a median
1.18× on birth/death (per seed 0.99–1.27×), 0.92–0.97× on sparse33 (median
about 8.5% slower) and 0.73–1.01× on binding12; none of these is resolved.
Generation did not repay its bundle compilation in these 128-evaluation runs.
Rebop macros were faster than reusable CRN by 3.44× on birth/death (per seed
2.99–4.35×) and 1.72× on binding12 (1.64–1.77×), and approximately tied on
sparse33 (0.94–1.00×), using optimizer time. Across all nine runs per route,
complete cost including one build was 13.21s for reusable CRN, 17.53s for
generated CRN and 12.88s for rebop macros; the 2.5% difference between reusable
CRN and rebop macros is within the run-to-run variation. These aggregated
results also include different cross-engine search paths and final points.

The new parameter interface and optional generator are runnable. Numerical
speed and optimization quality are reported separately below; the generated
backend specializes propensities and retains the core solver.

## birth_death

| Route | Optimizer s | Held-out s | First-run total s | Exact loss | Held-out loss |
|---|---|---|---|---|---|
| CRN fresh, fast | 0.1246 | 0.0831 | 2.0171 | 0.000776 | 0.000878 |
| CRN reusable, fast | 0.1272 | 0.0829 | 2.0069 | 0.000776 | 0.000878 |
| CRN generated, fast | 0.1078 | 0.0739 | 6.0127 | 0.000776 | 0.000878 |
| rebop expression | 1.8624 | 1.8965 | 3.7032 | 0.002768 | 0.003229 |
| rebop rebuilt LMA | 2.1787 | 2.0818 | 4.1188 | 0.002768 | 0.003229 |
| rebop macro | 0.0370 | 0.0195 | 2.2970 | 0.002768 | 0.003229 |

## sparse33

| Route | Optimizer s | Held-out s | First-run total s | Exact loss | Held-out loss |
|---|---|---|---|---|---|
| CRN fresh, fast | 1.9365 | 1.6922 | 5.4477 | 0.001761 | 0.001885 |
| CRN reusable, fast | 1.8080 | 1.6751 | 5.2956 | 0.001761 | 0.001885 |
| CRN generated, fast | 1.9626 | 1.6901 | 9.3571 | 0.001761 | 0.001885 |
| rebop expression | 44.6813 | 45.4836 | 90.1658 | 0.003450 | 0.003648 |
| rebop rebuilt LMA | 23.6343 | 24.4120 | 49.9980 | 0.003450 | 0.003648 |
| rebop macro | 1.8018 | 1.6212 | 5.6664 | 0.003450 | 0.003648 |

## binding12

| Route | Optimizer s | Held-out s | First-run total s | Exact loss | Held-out loss |
|---|---|---|---|---|---|
| CRN fresh, fast | 0.0453 | 0.0121 | 1.8642 | 0.001000 | 0.001132 |
| CRN reusable, fast | 0.0262 | 0.0126 | 1.8465 | 0.001000 | 0.001132 |
| CRN generated, fast | 0.0266 | 0.0113 | 5.8652 | 0.001000 | 0.001132 |
| rebop expression | 3.3750 | 3.3818 | 6.5518 | 0.000559 | 0.000765 |
| rebop rebuilt LMA | 3.3178 | 3.2864 | 6.6262 | 0.000559 | 0.000765 |
| rebop macro | 0.0153 | 0.0034 | 2.2617 | 0.000559 | 0.000765 |

## Binding cost

The following per-candidate native preparation costs are medians of each run's mean; the optimizer tables also include simulation and application overhead.

| Case | Fresh preparation µs | Reusable binding µs | Preparation gain |
|---|---|---|---|
| birth_death | 90.2 | 12.2 | 7.42× |
| sparse33 | 1003.6 | 150.8 | 6.65× |
| binding12 | 105.4 | 10.4 | 10.18× |

Times/losses are medians over the three declared optimizer seeds. Lower is
better. “Optimizer s” includes all objective calls, transport, reduction and
optimizer overhead. “First-run total” includes process setup/teardown, source
generation, the entire relevant bundle's dependency-cache-warm compilation,
optimizer time, exact-oracle assessment and 4096-path held-out validation.
It excludes dependency installation/compilation. Every seed and timing
observation remains in [the raw report](../validation/parameter-optimization/evidence/optimization.json).

## Repeated-use cost

The following totals cover all nine runs per route, charging the bundle's build
once. They illustrate reuse across models/seeds; they are not per-model compiler
cost estimates. Existing applications can also amortize the generic Rust worker
build over arbitrary later models.

| Route | Generation + build s | Nine optimizations s | Nine complete runs + one build s |
|---|---|---|---|
| CRN fresh, fast | 1.7608 | 6.3143 | 13.6324 |
| CRN reusable, fast | 1.7608 | 5.8928 | 13.2068 |
| CRN generated, fast | 5.7806 | 6.1865 | 17.5279 |
| rebop expression | 0.0000 | 146.1770 | 289.9288 |
| rebop rebuilt LMA | 0.0000 | 89.6286 | 180.5700 |
| rebop macro | 2.1944 | 5.7772 | 12.8773 |

## Quality and evaluation accounting

Targets use independent Poisson moments for the two linear systems and a
13-state generator/matrix exponential for reversible binding. Each objective
uses 32 trajectories with a fixed training seed; each final point receives 4096
new trajectories on an independent held-out seed. Exact loss evaluates the
fitted point's mathematical moments, without simulation noise. The objective
normalizes squared mean/variance residuals and excludes t=0. Parameter errors
below are the largest relative error over parameters and optimizer seeds.
A small moment loss does not guarantee parameter identifiability.

| Case | Route | Exact loss range | Largest parameter error | Actual objective calls |
|---|---|---|---|---|
| birth_death | CRN fresh, fast | 0.000471–0.008478 | 24.8% | 384 |
| birth_death | CRN reusable, fast | 0.000471–0.008478 | 24.8% | 384 |
| birth_death | CRN generated, fast | 0.000471–0.008478 | 24.8% | 384 |
| birth_death | rebop expression | 0.001464–0.013255 | 91.2% | 384 |
| birth_death | rebop rebuilt LMA | 0.001464–0.013255 | 91.2% | 384 |
| birth_death | rebop macro | 0.001464–0.013255 | 91.2% | 384 |
| sparse33 | CRN fresh, fast | 0.000857–0.002384 | 23.2% | 384 |
| sparse33 | CRN reusable, fast | 0.000857–0.002384 | 23.2% | 384 |
| sparse33 | CRN generated, fast | 0.000857–0.002384 | 23.2% | 384 |
| sparse33 | rebop expression | 0.002739–0.006245 | 34.8% | 384 |
| sparse33 | rebop rebuilt LMA | 0.002739–0.006245 | 34.8% | 384 |
| sparse33 | rebop macro | 0.002739–0.006245 | 34.8% | 384 |
| binding12 | CRN fresh, fast | 0.000859–0.001132 | 51.9% | 384 |
| binding12 | CRN reusable, fast | 0.000859–0.001132 | 51.9% | 384 |
| binding12 | CRN generated, fast | 0.000859–0.001132 | 51.9% | 384 |
| binding12 | rebop expression | 0.000376–0.000939 | 56.9% | 384 |
| binding12 | rebop rebuilt LMA | 0.000376–0.000939 | 56.9% | 384 |
| binding12 | rebop macro | 0.000376–0.000939 | 56.9% | 384 |

All route/model median exact losses met the predeclared <=.05 aspiration.
Parameter recovery was less precise: worst errors reached 51.9% for CRN binding
and 91.2% for a rebop birth/death run. These are finite-budget fits to noisy
32-path moments; a low moment loss is insufficient evidence of accurate rate
identification. Held-out loss and parameter error remain separate results.

The exact-loss <=.05 aim was aspirational, not a reason to repeat or drop runs.
Cross-engine streams differ, so different final points need not imply a solver
bias or a superior optimizer. Three search seeds and these three synthetic
models support only the reported scope.

## Compilation and boundaries

| Worker bundle | Measured rebuild s |
|---|---|
| dynamic | 1.7608 |
| generated | 5.7786 |
| rebop_macro | 2.1933 |

Model-independent dependencies were primed consistently. Each worker package
was cleaned and rebuilt separately, with the same standard Cargo release
profile. The generated bundle contains the three model kernels and an arithmetic
guard. A single-model production build could have a different cost. No
clean-machine installation, LTO tuning or other-platform claim is made.

CRN fresh/reusable/generated all use cached Direct, Ziggurat and **fast=True**,
one numerical worker. Fast's state checks, work budgets, stop polling and
accounting remain. Rebop macros use generated model-specific state/reaction
code; rebop dynamic uses runtime expressions or rebuilt mass-action definitions.
Rebop macros specialize the full reaction code and state representation; our
AOT route specializes propensity evaluation within the existing checked solver.
Native startup also includes the common worker framework/reference seed-list
preparation; it is not the minimum overhead of an in-process Rust application.
Dynamic rebop explicitly uses sparse storage on the 33-species model. Both
Python routes construct per-trajectory xarray outputs; native routes communicate
with a persistent worker and return reduced moments. These are complete
application costs with different allocation/check pipelines. Native bind and
simulation times are retained for CRN; event counts are not, so the engines'
different search paths cannot be normalized by work. The CRN worker also
serializes and hashes each result. Ratios of optimization times are not
isolated event-kernel throughput measurements.

The host is AMD EPYC 9V74 80-Core Processor, affinity allowance
9, CPU quota `800000 100000`.
All builds and heavy gates finished before optimizer timing. No multithread,
GPU, biological-network generalization or Julia performance claim is made.

## Validation and reproduction

The final engine passed **196 release tests** on Rust 1.88 and 1.85,
core/binding/worker clippy with warnings denied, strict schema/CLI, installed
Python, legacy SSA/CME and dependency gates. Both the checked-default and fast
families passed **210 exact-law comparisons each**, 100,000 paths per
fixture/variant, with the existing .001 family alpha and .008046120180018303
DKW epsilon unchanged. Normal/-O outputs match. Independent generated checks
cover 180 bindings/executions, 3150 arithmetic evaluations, cache refresh,
invalid inputs, cancellation, concurrency, worker replay and wrong topology.
The final gate logs record `check_parameter_interface.py` with SHA-256
`9f08cca6…`; the committed file (`8ee7cc7f…`) was edited afterwards and only
`parameter-validation.json` was refreshed with it. The evidence describes the
import commit `85a085b`; later commits change hashed sources.

[API guide](reusable-parameters.md),
[fixed protocol](decisions/reusable-parameters-2026-11-02.md),
[reproduction instructions](../validation/parameter-optimization/README.md) and
[implementation log](campaign/reusable-parameters-implementation-log-2026-11-02.md)
accompany the patch. Failed development outputs remain separate. At import,
the checked mode was the default and nothing had been pushed or published.
Since `crn-craft.ssa-run/0.3`, `fast` is the default of `ssa` requests, and
the interface was released in 0.2.0.
