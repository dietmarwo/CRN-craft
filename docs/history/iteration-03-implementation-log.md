# CRN-craft iteration 3 implementation log

Date: 1 November 2026

Milestone: **M1b — conservation cones, linkage analysis and reference export**

Version: `0.1.0-dev.3`

Repository: https://github.com/dietmarwo/CRN-craft

Baseline: **`c167add312a271671bdc31c95162b5b5c1227ee1`**, GitHub main after PR #2.

## Outcome

Implemented and tested the next planned iteration from current GitHub main. The analyzer now computes a complete nonnegative conservation cone when explicit resource/arithmetic limits permit, derives the missing conservation bounds highlighted in the review, and reports reaction linkage/strong/terminal classes and structural deficiency. A pinned GillesPy2 1.8.3 exporter produces executable Python for elementary stochastic mass-action models with explicit coefficient conversion and unsupported-feature errors.

The delivery is a patch ZIP against the exact baseline above and this Markdown log. No remote push, tag, package publication or deployment is part of the work. CRN-craft simulation remains unavailable; **M2 Direct SSA is next**.

## GitHub and review findings

Fetched actual GitHub main and read the repository instructions, implementation plan, source, CI and attached iteration-2 review. Main contains the delivered M1a implementation plus all three review regressions: non-unit steep-Hill `log1p`, weighted conservation bounds and high-order binomial symmetry. The mutation harness now requires a passing unmutated baseline and contains seven mutants. The attached review records successful hosted CI for the merge; this iteration does not claim to have run its new workflow remotely.

The untouched baseline passed 76 release tests, 48 contract checks and 218 independent network checks locally. Its existing rates, compiled-network operations, linear conservation algorithm and reviewed network tests were preserved. A later remote-ref check still reported `c167add`.

Accepted the review's recommendation to calculate the nonnegative conservation cone. Added explicit intermediate-ray and deterministic-work limits because a 64-species cap does not by itself make extreme-ray enumeration tractable. Kept the existing basis calculation as an independent result and sound fallback. The first Direct SSA should retain the current CRN-E042 positive-rate-underflow error; any future zero-with-counter approximation needs its own explicit option and assessment.

## Implemented changes

| Area | Delivered behavior | Main implementation |
|---|---|---|
| Nonnegative conservation | Orthant/hyperplane double description, support-based adjacency, primitive exact rays, checked residuals and initial totals | `crates/crn-craft/src/cone.rs` |
| Effective bounds | Minimum bound over complete cone rays; retain sound basis bounds when the cone is unavailable; stronger W001 suppression | `src/analysis.rs` |
| Reaction graph | Complete complexes including catalysts/empty complex, directed reaction edges, linkage/SCC/terminal classes, weak reversibility and deficiency | `src/linkage.rs` |
| Reference export | Generated Python for GillesPy2 1.8.3 elementary mass action, safe aliases, resolved constants, unit target volume, homodimer factor conversion | `src/export.rs` |
| Compatibility | Analysis 0.2; explicit `analyze-v1` for unchanged 0.1 output; discovery 0.3; new export 0.1 report | `src/analysis_v1.rs`, library/CLI/schema files |
| Validation | Independent support/reachability oracle, actual target-package comparisons, 14 source mutants with passing-baseline requirement, expanded CI | `validation/`, tests, workflow |
| Documentation | M1b analysis guide, exporter guide, migration/diagnostics, updated plan and agent context, evidence | `docs/`, README, changelog |

Short `src/` paths in the table are under `crates/crn-craft/`. No new external Rust dependencies were added; Cargo.lock changes only the workspace package versions. Python reference packages are optional test dependencies, not runtime dependencies. No GPL simulator/enumerator source was copied or linked.

## Conservation gap resolved

Run `crn-craft analyze models/analysis/cone-gap.model.json`. Its reactions are X2 → X0 + 2X1 and X3 → 2X2 + X0, starting at (0,0,0,1). The old mixed-sign basis proves no individual bounds. Exact rays (0,1,2,4) and (1,0,1,3) give effective bounds **(3,4,2,1)**. Current analysis suppresses the applicable growth warnings; legacy analysis retains its old warning policy.

Across the 304-model independent corpus (300 seeded networks plus four targeted models), bounds increase from **507 to 521 species**, including all four in the review counterexample. This is a newly pinned corpus, so those counts are separate from the review's 65-miss experiment. All 372 returned extreme rays matched independent exact support enumeration.

A successful cone is complete for nonnegative linear conservation laws. It is not a reachability decision, a nonlinear invariant solver or a proof that an unbounded species will grow without limit. Bounds refer to the initial conservation class even when another state is supplied for rate evaluation.

## Explicit numerical and resource policy

Cone enumeration uses checked i128 and supports at most 64 species, 256 reactions, 4,096 intermediate rays and 2,000,000 deterministic work units. It reports `size_limit`, `ray_limit`, `work_limit`, `arithmetic_overflow` or `verification_failed` without partial rays/bounds. Rust callers may lower the ray/work limits; higher values are clamped. Structural/rate analysis and independently computed basis bounds remain available on cone failure.

The old linear-basis report retains its original meaning. New `analysis.initial_class_upper_bounds` selects exact cone bounds when available, otherwise sound basis bounds. New `nonnegative_conservation` exposes completeness, limits, consumed work and rays. New `linkage` is structural, includes zero-rate channels and does not apply a deficiency theorem to arbitrary kinetics. Graph traversal is iterative and tested on a 4,096-complex chain.

The rate evaluator, availability gates, exact i64 updates, transactional failures and strict underflow policy are unchanged.

## Reference exporter boundary

```bash
cargo build --release --workspace --locked
./target/release/crn-craft export models/examples/enzyme-mechanism.model.json gillespy2
```

The command emits JSON with raw model hash and generated Python source; it does not execute the target. `build_model()` requires GillesPy2 1.8.3. Target volume is fixed at one because source coefficients already use stochastic counts; source Omega does not rescale them. A homodimer coefficient is divided by two to match CRN-craft's binomial convention. Original species/reactions map to deterministic collision-free aliases.

Supported reactant order is zero, one or two. Custom Hill/MM laws, higher order, declared population safety limits and initial counts above 2^53 are rejected with CRN-E050. Inexact/underflowed homodimer halving is CRN-E051. Failure emits no source artifact. The target uses floating-point counts and different extreme-rate arithmetic; target ODE semantics, live budgets and CRN-craft range/underflow diagnostics are not preserved.

Actual installed target models passed **15,717 propensity comparisons**. All small dyadic count-grid cases were exact; the maximum relative error over the remaining examples was **1.1328806373726087e-16**. A two-trajectory seeded NumPySSASolver smoke run also passed. This verifies the exporter path, not an implemented CRN-craft simulator or a performance comparison.

## Contract migration

- Model/run/check/result/state 0.1 schemas are unchanged, with existing model/run byte hashes preserved.
- `analyze` now emits **analysis 0.2**, adding cone, graph and effective-bound fields.
- `analyze-v1 MODEL [STATE]` and `analysis_v1::analyze_bytes` retain the original 0.1 output and warnings. Its generated schema is byte-identical to the old snapshot.
- Discovery is **0.3** and lists the implemented commands, limits and reference-export subset.
- Export starts at **0.1**. All data output remains JSON with success/error exit codes 0/1.

Clients depending on the old analysis shape can select `analyze-v1`; clients using current analysis should inspect cone status before treating its ray set as complete. All eight schema snapshots are checked against generated output.

## Tests and evidence

| Check | Outcome |
|---|---|
| Baseline before changes | 76 Rust tests; 48 contract and 218 network checks passed |
| Rust 1.85.0 / 1.88.0 / stable 1.99.1 | 97 release tests on each; fmt and clippy with warnings denied |
| Contract/CLI checker | 58 passed, normal and `-O` reports identical |
| Preserved M1a network oracle | 218 passed, normal and `-O` reports identical |
| New exact cone/graph oracle | 308 passed, normal and `-O` reports identical |
| Real target exporter | Six grouped checks, 15,717 comparisons and reference smoke run; normal and `-O` identical |
| Mutation gate | 14 compiled mutants caught after a passing unmutated baseline |
| Dependency gate / Ruff | Passed; dependency reports normal and `-O` identical |

New Rust tests cover exact primitive rays, the review gap, weighted bounds, adjacency, duplicate/reversed equations, open/empty/null networks, both sides of size/ray/work caps, arithmetic/total overflow, graph identity/direction, legacy output, exporter restrictions and CLI errors. An additional regression forces cone overflow while linear analysis still succeeds, proving the aggregate bound fallback stays intact.

The structural oracle uses exhaustive support enumeration and Fraction nullspaces, rather than copying the Rust double-description algorithm. Its graph checker uses transitive closure. Four corrupted-output probes exercise the checker itself. Seven new source mutants cover cone combination/adjacency/bound division, graph direction/complex identity, homodimer conversion and loss of fallback bounds; the seven reviewed M1a mutants remain.

Structural fixture seed: 20261103; SHA-256 `704c72253352754584904bb700702326a992bc6f78fbbe4ff3b803da0dec94ee`. Numerical tolerances remain the predeclared M1a 2e-12 relative plus eight minimum subnormals. Targets are recorded in `docs/decisions/iteration-03-validation.md`. Raw reports and environment/package provenance are under `validation/iteration-03/`; earlier evidence is preserved.

## Development issues retained

Rust 1.85 clippy requested a match for a comparison chain; stable clippy requested a reference slice rather than a clone in one test. Both were fixed without changing numerical semantics, and failed lint logs are retained. One subsequent checker launch found its local built executable at mode 0644; restoring execute permission resolved it without source changes. Its traceback is retained, and final delivery verification uses an independent fresh build.

The new fallback mutant initially matched both the cone and fallback expressions, so the harness rejected its ambiguous anchor before scoring it. The anchor now includes the line prefix to select only the fallback. The failed harness log is retained; the full gate was rerun with a passing baseline and all 14 mutants caught.

No numerical/structural fixture failed and no tolerance was relaxed. Intentional mutation failures are recorded separately and count only after the unmutated suite passes. Logs retain diagnostic content; surplus blank lines at EOF are normalized for clean patch whitespace.

## Remaining work and performance position

**Next: M2**, the first CRN-craft Direct SSA scheduler, pinned samplers, deterministic ensembles, live budgets/cancellation and verifiable result archives. M3 adds FSP and broader stochastic/statistical validation. M4a provides the first meaningful pinned simulator quality/cost comparisons. M5 adds bindings and the circuit application; M4b advanced aggregators are demand-driven. Ecosystem workers, ODE/sensitivity paths, imports and optional tau-leaping remain later milestones.

The analyzer and one elementary exporter now cover the M1b work needed to proceed. Wider target/custom-rate exporters and expressions should follow concrete fixtures. There is still **no CRN-craft throughput measurement or performance-gap estimate**, because it has no executable simulator. Bound coverage and target propensity agreement are correctness evidence, not simulation speed claims.

## Delivery verification

The ZIP contains one Git mail patch, `APPLY.md`, this log, checksums and post-commit verification evidence. Apply to a clean branch at the exact baseline above. The exact final patch is checked with `git apply --check`, applied with `git am` to a fresh baseline checkout using `core.autocrlf=true`, compared by resulting tree and rebuilt/tested there. The final delivery appendix records those results after the patch is frozen, avoiding self-referential commit and patch hashes.
