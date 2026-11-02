# Reusable parameters and generated kernels: declared protocol

2 November 2026, America/La_Paz. Base: GitHub main `c16fcc5` (merged PR #8).
Read the attached fast-switch log; its projected Julia figures are historical,
not measurements for this campaign. No push, publication or branch deletion.

## Implementation

Add Rust/Python parameter templates for fixed topology. A binding supplies all
declared parameter values in a stable named order, validates their domains,
shares immutable stoichiometry/dependencies, resolves fresh rates and Hill
memos, and regenerates custom-bound groups. It produces complete model/request
snapshots and their real hashes. Bindings/executions are independent, including
concurrent calls, seeds and worker counts. Initial populations remain fixed.

Add an optional Rust source generator for model-specific mass-action propensity
kernels. Generated kernels use the existing SSA schedulers, fast/default RNG,
checked updates, stops, budgets and reduction. Guard the topology identity and
retain the checked evaluator for exceptional arithmetic and unsupported orders.
This specializes propensities, not the entire solver or every rate-law family.
The default engine gains no build-time compiler, unsafe code or dependency.

## Correctness gates

Before timing: fmt, core/binding clippy -D warnings, release tests --locked on
Rust 1.88 and 1.85, schema/CLI and dependency gates normally and under Python -O,
Python installed/prepared/concurrency tests, legacy SSA/CME identities and both
unchanged exact-law gates: default and --fast, 210 marginal comparisons each,
100,000 paths, family alpha .001, epsilon .008046120180018303. Normal/-O reports
must match. Keep all failures; never loosen these rules.

Test binding/fresh preparation and generated/dynamic equivalence across every
algorithm/wait/fast combination; parameter changes including zero/subnormal
rates and custom-law domain boundaries; cache-group splits/joins and Hill memo
reset; malformed values; seed/thread replay; wrong topology, cancellation and
work/output limits. Numerical identity is separate from model/run hashes and
timings. Generated mass-action arithmetic must preserve the existing order.

## Complete optimization experiments

Use fcmaes-rust 0.2.1 BiteOpt, 128 requested evaluations, optimizer seeds
20261102/20261103/20261104, fixed starting point and log-rate bounds. Record
actual objective calls and every evaluated point/loss; never drop a failed run.
Fit two positive rates on three frozen models: immigration/death, 33 independent
immigration/death species (66 reactions), and reversible A+B <-> C (finite
13-state chain). Exact time-course mean/variance targets come from the Poisson
laws or a separately constructed finite-state generator/matrix exponential.
Model definitions, targets and hashes are frozen before the first experiment.

Each candidate uses 32 paths, a fixed common training seed within its engine,
and the same uniform six-point save grid and normalized mean/variance loss.
**Amendment (before any optimization experiment):** uniform grids match rebop
Python's public output API: horizons 10, 20 and 2 respectively (binding12 moved
from horizon 10 with seven points to horizon 2 with six, so its transient
identifies both rates). Initial nonuniform draft targets are retained as
development evidence. Cross-engine random
streams differ. Three native CRN fast/Ziggurat routes: fresh preparation,
reusable parameter binding, generated propensities. Use cached Direct for all
three so specialization is separated from algorithm choice. Include rebop
0.10.0 dynamic parameter-expression and dynamic rebuilt-mass-action APIs and
generated rebop macro kernels. Rebop dynamic uses sparse=True for the
33-species model and False for the two small models, fixed before timing.
One numerical worker, no overlapping heavy jobs.

For each final point report training loss, exact-moment loss, parameter error,
and independent-seed 4096-path held-out loss. These quality measures are
descriptive; median exact-moment loss <=.05 is an aspirational fit target, not
permission to rerun or discard a seed. CRN routes must evaluate identical points
bit for bit and produce the same optimizer trace. A mismatch blocks timing.

Record generation, compilation, process startup, all objective evaluations,
optimizer overhead, held-out evaluation, and total time. Build the dynamic,
generated-CRN and rebop-macro workers separately with model-independent
dependencies primed and each worker package rebuilt once. Label these as
dependency-cache-warm compilation costs, retain all build logs, and charge the
whole three-model bundle build in each single-run total. Also show repeated-run
amortization. Do not claim a clean-machine installation benchmark or extrapolate
the previous Julia gap. Persist raw observations and report slower kernels.

Deliver one aggregate patch ZIP applying to the fetched base, an implementation
log, the comparison, code generation sources and reproduction/evidence files.
Verify clean-base patch application, exact Git tree reconstruction, ZIP CRC and
per-member hashes, and provide direct downloadable files.
