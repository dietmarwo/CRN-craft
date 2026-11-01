# Agent context: CRN-craft M2

Version 0.1.0-dev.4. Start with `crn-craft describe` (discovery 0.4).
Implemented: Direct SSA; pinned ChaCha8 sampling; ordered parallel ensembles;
mean/sample variance/mean SE; event, population, wall-time, output and workspace
limits; cooperative cancellation; checksummed result archives. Reviewed network,
rate, linear/cone/linkage analysis and GillesPy2 elementary export remain.
FSP, ODE, fitting, bindings and SBML import are later milestones.

Read README.md, docs/direct-ssa.md, docs/archive-contract.md and
validation/iteration-04/RESULTS.md. `simulate MODEL RUN NEW_DIRECTORY` executes
Direct SSA. `--cancel-file PATH` enables cooperative external cancellation.
`verify DIRECTORY` checks integrity and relational consistency, not biological
validity or stochastic accuracy. A complete archive may record an incomplete
simulation. `simulate` exits 0 only for full completion; `verify` exits 0 for
valid incomplete evidence as well. Never overwrite an existing output folder.

Result is 0.2, archive/archive-check 0.1. Reserved result 0.1 remains result-v1.
Model/run/check/state and current/legacy analysis/export contracts remain
unchanged; historical operation flags do not replace the discovery catalog.
Raw hashes identify consumed bytes. Preserve LF and exact model/run snapshots.

Mass action uses binomial count factors; custom Hill/MM laws are full gated
propensities. Omega never silently rescales coefficients. Positive underflow
stays CRN-E042. Unresolved selection ranges or non-advancing time fail explicitly.
No last-channel fallback, reflective population boundary or hidden clipping.

Every replicate is accounted for. Statistics are withheld for any incomplete
ensemble, because completion may be state-dependent. Mean/variance reduction is
ordered by replicate ID, with one pinned stream per root-seed/replicate pair.
Compare numerical payloads separately from run hashes and timing/build/thread
provenance. Same-platform 1/2/8-thread identity does not prove cross-platform identity.

Next prompt:

> Implement M3 FSP with reachable-state enumeration, an explicit escaping-mass
> sink and stable sparse uniformization. Keep physical leakage, Poisson-series
> truncation and numerical error separate. Predeclare analytic/DSMTS gates;
> independently check them and retain failures. Preserve all M2 statistical,
> sampler, archive, mutation and review regressions. Do not claim a rigorous
> floating-point error bound or a speed comparison without the relevant evidence.

M4a then measures pinned simulator quality/cost; M5 supplies bindings and the
circuit application. Advanced SSA aggregators follow measured bottlenecks.
