# Reusable parameters and generated-kernel optimization — implementation log

2 November 2026, America/La_Paz. Fresh GitHub base: `c16fcc5` (PR #8).
The imported C01–C10 work and optional `method.fast` switch are already in this
base. No remote was pushed, no package published and no branch deleted.

## Declared scope

Implement a reusable parameter interface, optional generated propensity kernels,
and complete optimization comparisons with rebop. The protocol was committed
before numerical/optimization experiments in
`docs/decisions/reusable-parameters-2026-11-02.md`. Default checked semantics and
fast request/sampler identities remain unchanged. The attached fast-switch log
was read as context; its Julia projection is not a new measurement.

## Implementation

1. **Reusable bindings.** Added Rust `parameters::ParameterTemplate` and Python
   `prepare_parameters`, with named or ordered values, read-only names/defaults,
   optional seed overrides, independent bound executions and full model/request
   snapshots. Compiled species, stoichiometry, dependencies and reverse readers
   use shared immutable storage. Each binding validates the candidate, resolves
   fresh rates/Hill memos, recomputes custom-bound groups and reruns admission.
   It still clones/validates the model and serializes snapshots per binding.
2. **Generated kernels.** Added a mass-action Rust source generator and CLI
   example. Elementary zero/first/dimer/pair propensities are specialized;
   exceptional arithmetic and higher orders retain the checked fallback.
   `ssa::execute_with_kernel` guards fixed topology. The core schedulers,
   bounds, counters, checked population updates, stops and reduction are shared.
   Custom Hill/MM laws use the dynamic interface. This is propensity
   specialization, not a generated whole solver or a runtime compiler.
3. **Optimizer comparison.** Added an isolated locked worker workspace with
   dynamic CRN, generated CRN and rebop macro build variants, plus two Python
   rebop routes. Frozen models cover two-rate small/sparse linear systems and
   reversible binding. Independent targets use Poisson moments and a separately
   constructed 13-state generator/matrix exponential. All native CRN timing
   routes use cached Direct/Ziggurat with fast=True.
4. **Quality and cost.** Complete fcmaes-rust BiteOpt runs retain every point,
   actual evaluation count, training loss, exact-moment loss, parameter error,
   independent held-out loss and timing components. Single-run totals include
   source generation, dependency-cache-warm compilation of the complete bundle,
   worker startup, optimizer/objectives and held-out evaluation. Repeated-run
   results amortize compilation once.

Existing request/result/archive schemas and dependency gates are unchanged.
The default engine gains no optimizer, rebop, compiler or unsafe code. Existing
M6 job adapters retain their legacy route; the new public interfaces can be
adopted without replacing the job protocol. Initial populations are fixed within
one template; structural/model/run changes require a new template.

## Development findings retained

- The initial API grids were nonuniform. Before the first optimization
  experiment they were replaced with uniform six-point grids to match rebop's
  Python output API; the draft target file is retained. No statistical rule was
  changed. Sparse rebop storage for the 33-species model was also fixed before
  optimization timing.
- Initial worker builds caught missing Diagnostic formatting/borrow handling
  and unnecessary macro parentheses. Corrected source and failed logs remain.
- The first parameter checker passed all 180 candidate/mode identities, then
  found a missing dispatch arm for the extra generated arithmetic guard. The
  checker branch was fixed and the full check rerun; the partial report remains.
- The pinned Ruff version required union syntax for a new isinstance check;
  the failed gate output is preserved. The final gate runner marks a report
  complete only after every check and normal/-O identity has finished.

## Validation and measured results

Final executed gates, optimization results and package verification are appended
below after the declared runs. Raw evidence and reproduction instructions live
under `validation/parameter-optimization/`.

## Final executed validation

- **196 release tests** passed on Rust 1.88 and the Rust 1.85 minimum.
- Core, Python binding and optional worker clippy passed with warnings denied;
  formatting and pinned Ruff checks passed. The rebuilt wheel was installed and
  checked through isolated Python, CLI and prepared APIs.
- Schema/CLI, dependency, legacy SSA/CME and cached identities passed. Default
  and fast exact-law families each passed **210/210** fixed comparisons at
  100,000 paths per fixture/variant (7 million paths per family pass). Both
  families were repeated under Python -O with byte-identical reports. The
  family alpha .001, epsilon .008046120180018303 and wrong-law power checks were
  retained unchanged. No acceptance rule was loosened.
- New checks passed **180** binding/fresh/generated execution identities over
  all 20 algorithm/wait/fast combinations, **3150** generated arithmetic checks,
  **13** malformed/domain/unsupported cases, concurrent candidate isolation,
  cancellation, budget stops, seed/worker replay and wrong topology rejection.
  Normal/-O reports match. The final optional worker request parser also rejects
  unknown fields and malformed optional fields; those nonnumerical harness
  changes were rebuilt/rechecked before timing.

## Complete optimization results

All **54 runs** completed: three frozen models × six routes × three optimizer
seeds, **6912 actual objective calls**. All nine comparisons of the three CRN
routes have identical full point/loss traces and final solutions. Source and
binary hashes accompany every measured build. No heavy check or build ran
alongside optimizer timing. Nothing was removed from the measured report.

Median optimizer time, first-run total and exact-moment loss are shown below.
The total includes full-bundle dependency-cache-warm compilation, generation,
startup, optimization, exact assessment and 4096-path independent held-out work.
Dependencies/installations are outside this timing boundary.

| Model | Route | Optimizer seconds | First-run total seconds | Exact loss | Held-out loss |
|---|---|---:|---:|---:|---:|
| birth_death | crn_fresh_fast | 0.124566 | 2.017115 | 0.000776 | 0.000878 |
| birth_death | crn_reuse_fast | 0.127162 | 2.006920 | 0.000776 | 0.000878 |
| birth_death | crn_generated_fast | 0.107825 | 6.012666 | 0.000776 | 0.000878 |
| birth_death | rebop_expression | 1.862447 | 3.703214 | 0.002768 | 0.003229 |
| birth_death | rebop_rebuild_lma | 2.178711 | 4.118827 | 0.002768 | 0.003229 |
| birth_death | rebop_macro | 0.036952 | 2.296990 | 0.002768 | 0.003229 |
| sparse33 | crn_fresh_fast | 1.936490 | 5.447709 | 0.001761 | 0.001885 |
| sparse33 | crn_reuse_fast | 1.807974 | 5.295595 | 0.001761 | 0.001885 |
| sparse33 | crn_generated_fast | 1.962553 | 9.357135 | 0.001761 | 0.001885 |
| sparse33 | rebop_expression | 44.681275 | 90.165805 | 0.003450 | 0.003648 |
| sparse33 | rebop_rebuild_lma | 23.634331 | 49.997977 | 0.003450 | 0.003648 |
| sparse33 | rebop_macro | 1.801764 | 5.666436 | 0.003450 | 0.003648 |
| binding12 | crn_fresh_fast | 0.045295 | 1.864211 | 0.001000 | 0.001132 |
| binding12 | crn_reuse_fast | 0.026207 | 1.846490 | 0.001000 | 0.001132 |
| binding12 | crn_generated_fast | 0.026581 | 5.865205 | 0.001000 | 0.001132 |
| binding12 | rebop_expression | 3.375010 | 6.551790 | 0.000559 | 0.000765 |
| binding12 | rebop_rebuild_lma | 3.317760 | 6.626167 | 0.000559 | 0.000765 |
| binding12 | rebop_macro | 0.015255 | 2.261667 | 0.000559 | 0.000765 |

Reusable binding improves optimizer time about 7% on sparse33 and 1.73× on
binding12; birth/death is flat within this small experiment. Generated
propensities gain 1.18× on birth/death but regress about 8.5% on sparse33 and
are flat on binding12. Their compilation cost was not recovered in these runs.
They remain an explicit optional backend. The existing dynamic engine already
has elementary rate kernels; this experiment does not establish an advantage
for specializing just propensity dispatch on every model.

Native preparation itself improved substantially: fresh/reusable per-candidate
times were 90.2/12.2 microseconds on birth/death (7.42×), 1003.6/150.8 on
sparse33 (6.65×), and 105.4/10.4 on binding12 (10.18×). Simulation, accounting
and transport limit how much this changes full optimizer time.

Rebop's macro route is much stronger than its Python API results suggest:
optimizer time is 3.44× faster than reusable CRN on birth/death, 1.72× on
binding12 and roughly equal on sparse33. Complete costs across nine runs per
route, with compilation charged once, are 13.21s for reusable CRN, 17.53s for
generated CRN and 12.88s for rebop macros. Python rebop routes take 180.57s
(rebuilt mass action) and 289.93s (runtime expressions), including held-out work.
Their per-trajectory xarray construction is part of that application cost.
These ratios describe full fitting workflows, not isolated kernel throughput.

The bundle rebuilds took 1.7608s (dynamic), 5.7786s (generated) and 2.1933s
(rebop macro), with model-independent dependencies primed and each worker
package rebuilt separately in the same release profile. A single-model build
or a different release profile can change this tradeoff. The three-model bundle
and arithmetic guard were charged in full; repeated-use totals charge once.
Reference seed preparation and worker transport are included in native startup
and objective times; an in-process Rust caller can have lower overhead.

## Final solution quality and limitations

All route/model median exact losses met the declared <=.05 aspiration. This
supports moment fitting on these fixtures, rather than precise rate recovery:
worst relative parameter errors reached **51.9%** on CRN binding and **91.2%**
on one rebop birth/death run. The fixed 128-evaluation budget and noisy 32-path
training objective leave substantial uncertainty. The report retains every
parameter error, exact loss and 4096-path held-out loss. Cross-engine random
streams differ, so different points/losses are not an optimizer ranking.

Generated CRN propensities retain the existing schedulers, budgets, accounting,
checked state changes and reduction. Rebop macros specialize their full model
code. No instruction-profile attribution is claimed. Parameter bindings still
validate/clone/serialize the candidate model; topology compilation is shared.
Initial populations remain fixed per template. AOT covers mass action and keeps
checked fallbacks; custom-law bindings are supported dynamically. Existing M6
job adapters were not migrated, and there is no automatic compiler/dynamic
loader in the Python runtime. No multithread, other-platform, biological-model
or new Julia-gap claim is made.

The full comparison is `docs/comparison-parameter-optimization.md`; instructions
and raw evidence are under `validation/parameter-optimization/`. The shipped
log contains all outcomes, including slower generated variants and development
failures. The default remains checked and fast is explicit. No remote was
pushed or package published.

## Review notes (2 November 2026, after import)

Two independent reviews (core code; evidence and documentation) found no
blocking correctness problem: bound runs, generated-kernel runs and fresh
preparation were byte-identical in every tested combination, and every number
in the comparison recomputes from the raw evidence. Fixed afterwards in commits
on `reusable-parameters`:

- Generated kernels fell back to the state-scanning `propensity` for
  higher-order reactions (about 17× slower than the default path on a
  1,500-reaction review model with `direct`); they now call
  `propensity_in_scheduler_state` (1.7–1.8× slower there; still measure before
  using generated kernels on very large networks).
- Ordinary preparation had started compiling the dependency graph before the
  admission checks (`CRN-E040` instead of `CRN-E060` for an oversized request);
  the old order is restored with a regression test.
- Binding mismatches use `CRN-E014` and list missing and unknown names; Python
  seed errors use `CRN-E020`; Python rejects sets and maps float overflow to
  `CrnError`.
- The comparison now gives per-seed ranges. Fresh and reusable CRN simulate
  identical paths, yet their optimizer-time ratio varied 0.85–1.16× across
  seeds on birth/death, so the generated-kernel effects (median 1.18× on
  birth/death, about 8.5% slower on sparse33) and the 2.5% gap between
  reusable CRN and rebop macros are not resolved. The claim that event counts
  were retained was wrong and is corrected.
- Provenance: the import is one squashed commit, so "protocol committed
  before experiments" cannot be checked from history; the grid amendment is now
  marked in the protocol. The final gate logs record
  `check_parameter_interface.py` with SHA-256 `9f08cca6…`, the committed file
  is `8ee7cc7f…`. Reproduction steps gained `cargo fetch`, the 1.85.0
  toolchain, the run and report steps; the run no longer requires
  `/sys/fs/cgroup/cpu.max`; the packager no longer leaves a success claim when
  ZIP verification fails.
- Not changed: the kernel fingerprint covers the model, not the CRN-craft
  version (documented: regenerate kernels after upgrading); `describe` does not
  list the library-level parameter API.

