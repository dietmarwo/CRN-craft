# Agent context: CRN-craft campaign

Start with README.md, docs/campaign/README.md and `crn-craft describe` (0.6).
Default engine MSRV is 1.85; Python/jobs/diffsol are isolated optional workspaces.
Legacy exact SSA commands are simulate and simulate-cached, with unchanged original
model/run/result/archive contracts. The separate ssa command implements Direct,
cached Direct, DirectCR, SortingDirect and RSSACR with inverse or Ziggurat waits.
Read docs/alternative-ssa.md and the five-step acceleration protocol before changes. FSP executes via fsp; tau is explicitly
approximate and has separate request/result schemas. Optional crn-diffsol uses
BDF and forward sensitivities with deterministic count kinetics, not exact
nonlinear stochastic means. Importers accept explicit subsets/conventions only.

Completed exact fixed-work results must preserve stream/reduction identity.
Incomplete ensembles withhold statistics; wall-time/cancellation work is outside
fixed-work reproducibility. Keep sampler, archive, conservation, CME and mutation
review regressions. Check discovery to distinguish executable capabilities from
optional/source-only interfaces and experimental approximations.

Next prompts should name a concrete fixture, gate and budget. Read the original
implementation plan and campaign acceptance before numerical changes. Retain failed
rows; do not loosen tolerances after a failure. The Vilar and circuit diagnoses are
exploratory; broad behavioral claims, other-platform wheels and
owner-controlled publishing remain additional validation work.

Fixed-topology parameter searches can use the Rust/Python parameter template API;
read docs/reusable-parameters.md and docs/comparison-parameter-optimization.md.
Optional generated mass-action propensities retain the existing checked solver;
fast is still explicit and the dynamic numerical engine has no compiler dependency.
