# Agent context: CRN-craft campaign

Start with README.md, docs/architecture.md, docs/campaign/README.md and `crn-craft describe` (0.10; 0.13 from the native CLI).
Default engine MSRV is 1.85; Python/jobs/diffsol are isolated optional workspaces.
Legacy exact SSA commands are simulate and simulate-cached, with unchanged original
model/run/result/archive contracts. The separate ssa command implements Direct,
cached Direct, DirectCR, SortingDirect and RSSACR with inverse or Ziggurat waits,
and next reaction (ssa-run/0.4, inverse waits, checked policy only). Its
`method.fast` policy is on by default in ssa-run/0.3 and 0.4 (`crn.ssa_spec`
writes 0.3, or 0.4 with `fast=False` for `next_reaction`); set `fast: false`
for the checked policy. ssa-run/0.1 and 0.2 keep their checked defaults and
result bytes.
Dynamic attributed species use the experimental rules engine (crate
crn-craft-rules, CLI rules-check and rules, Python crn_craft.rules): read
docs/rules-language.md and docs/mtdna-rules.md, and docs/specialized-kernels.md
before adding a model-specific rules kernel. docs/mtdna-user-manual.md,
docs/mtdna-optimization-tutorial.md and docs/mtdna-performance-analysis.md,
with examples/mtdna/, show the mtDNA workflows end to end.
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
exploratory; broad behavioral claims remain additional validation work. From
0.2.0 the Python package is on PyPI (`pip install crn-craft`) and the crates
`crn-craft`, `crn-craft-rules` and `crn-craft-cli` (`cargo install
crn-craft-cli --locked`) are on crates.io; the diffsol integration and the job
worker build from the repository.

Fixed-topology parameter searches can use the Rust/Python parameter template API;
read docs/reusable-parameters.md and docs/comparison-parameter-optimization.md.
Optional generated mass-action propensities retain the existing solver and its
policy choice; the dynamic numerical engine has no compiler dependency.

## Optimization lessons

Use [C1–C4](../tutorials/optimization/README.md) for runnable calibration, optional BDF/SQP fitting, native MODE/MAP-Elites oscillator searches, and bounded rules inference. These use library or command interfaces; no new stock worker profiles are implied. Their SSA runs use the default `fast` policy; C1's `--include-checked` adds a separate checked-policy arm. [Validation](../tutorials/optimization/VALIDATION.md) distinguishes completed checks from fitted-model quality.
