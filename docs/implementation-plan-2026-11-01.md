# CRN-craft: reviewed implementation plan, revision 2

**Date:** 1 November 2026  
**Basis:** the attached 31 October CRN-craft plan, GillesPy2 alternatives note, and 29 October ecosystem review.  
**Iteration 1:** M0, contracts and preflight validation. Numerical simulation is not part of this iteration.  
**Original iteration-1 delivery:** a new local project suitable for a future `dietmarwo/CRN-craft` repository. No repository creation, push or release is required to test it.

## Iteration 2 progress — 1 November 2026

M1a is implemented in `0.1.0-dev.2`, based on GitHub main `8134f51`: compiled indices, sparse stoichiometry, availability gates, dependency graphs, all four rate laws, checked reaction updates and `analyze`. Bounded exact conservation also ships now (64 species / 256 reactions, checked i128); unavailable outcomes are explicit. The four existing 0.1 schemas remain unchanged, new state/analysis contracts are 0.1, and discovery advances to 0.2. See the iteration 2 log and acceptance evidence for executed checks and limitations.

M1b remains linkage classes, reference exporters and only fixture-driven expressions; arbitrary-precision or complete nonnegative conservation analysis is not claimed. M2 remains the first executable Direct SSA. The original milestone text below is retained as the design baseline.

## 1. Review and decisions

Keep the central design: a small native Rust engine for well-mixed reaction networks, explicit models and budgets, exact SSA before approximations, Python for the circuit study, and existing optimization-jobs/experiment-craft interfaces for execution and campaigns. Keep Julia, rebop and GillesPy2 as optional, pinned references. Do not reproduce Catalyst's symbolic framework or add an orchestration server.

The original plan is technically promising but too broad in several gates. The following changes make it implementable and keep its evidence claims defensible.

| Original issue | Revised decision | Consequence |
|---|---|---|
| The illustrative repressilator includes only one of its three gene modules | Ship a complete six-species, twelve-reaction example | Examples can become actual fixtures in M1/M2 |
| Population bounds serve both SSA and FSP | Separate `population_limit` safety stops from FSP `bounds` | Crossing a safety limit stops execution; FSP leakage goes to a sink, never a reflected state |
| FSP declared impossible for unbounded networks | Unbounded networks can use a finite truncation, with measured escaping mass | Reject missing truncation bounds, not unbounded physical state spaces |
| ODE scaling is promised but not defined | Specify combinatorial factors, Ω mapping and units below | Prevent homodimer factor-of-two errors and incorrect mean claims |
| Custom rate laws might be mistaken for mass-action coefficients | Hill/MM/expression rates are complete propensities, with reactant-availability gating | Dependency analysis includes both rate reads and reactant availability |
| “Cyclic expression” fixture for a JSON tree | No references in the initial AST; test excessive depth and unsupported reference syntax | Do not invent a graph language merely to test cycles |
| Whole summaries must be bit-identical across thread counts | Compare numerical payloads for completed fixed-work runs; separate timing/provenance | Deadlines, cancellations and metadata are outside byte-identity guarantees |
| Conservation bases must equal an oracle basis | Compare null spaces and rank, plus exact residuals | Equivalent bases are not falsely rejected |
| A binomial count of DSMTS pointwise failures is proposed | Do not assume independence across times/species | Predeclare suite-level rules that allow dependent observations |
| All observables promise standard errors and all scores fixed success rates | Specify an estimator and uncertainty method per observable; predeclare held-out operating conditions | No unsupported SE for variance, quantiles, correlations or censored passage times |
| RSSA-CR/NRM described as universally sublinear | Benchmark time per event versus size, dependency degree and rate spread | Dense coupling can defeat sparse-update advantages |
| A single throughput gate both blocks and runs alongside M5 | M4a baseline measurements after M3; M5 can use verified Direct while M4b is demand-driven | The useful circuit application does not wait for three new algorithms |
| Tau-leap bias must decrease monotonically in noisy estimates | Compare uncertainty-aware bias bounds over tighter tolerances | Sampling noise is not mistaken for numerical regression |
| Raw model hashes and normalized archives are conflated | Preserve consumed bytes; hash normalized representations separately, if used | The bytes checked are the bytes identified in evidence |
| Current package versions and “only maintained tool” statements are too strong | Pin actual reference environments at the milestone; avoid exclusivity or blanket speed claims | The alternatives note is context, not benchmark evidence |

The work is feasible without the home lab. Baseline and correctness runs can be executed in this environment or CI. Seven-node commissioning and home-hardware performance remain separate work after access is available.

## 2. Scope and mathematical contract

The initial process is a time-homogeneous continuous-time Markov chain on non-negative integer species counts. One well-mixed compartment, constant system size, and time in seconds are supported. No spatial model, delays, events, variable volume, rule expansion, GPU execution or biological validity claim is implied.

### Mass action and Ω

For reactant multiplicities νᵢ and counts xᵢ, the SSA propensity is

\[
a_j(x)=c_j\prod_i {x_i\choose\nu_{ij}}.
\]

Thus `2A -> ∅` has propensity `c*A*(A-1)/2`, while `A+B -> C` has propensity `c*A*B`. Source reactions have propensity `c`. Stored coefficients are stochastic count coefficients; changing Ω alone does not silently change the SSA.

Let q be the total reaction order and z=x/Ω. If the macroscopic concentration law is `v(z)=k*product(z_i^ν_i)`, then

\[
c=k\,\Omega^{1-q}\prod_i \nu_i!,\qquad
k=\frac{c\,\Omega^{q-1}}{\prod_i\nu_i!}.
\]

The reaction-rate ODE uses powers, not falling factorials. At a fixed model Ω, its count-space rate is `c*product(x_i^ν_i/ν_i!)`, and the concentration derivative divides the stoichiometric count flux by Ω. The ODE is a macroscopic approximation for nonlinear networks; it is not generally the exact mean equation. Affine propensities provide an exact mean cross-check when moments exist and safety stops do not intervene. Catalyst documents the same distinction between combinatorial jump and power-law ODE rates [S1].

M0 validates Ω and coefficients but does not calculate rates or convert units. M1 implements these formulas with exact small-count tests. M7 must reject unsupported ODE interpretations rather than guess them.

### Other rates

Hill activation: `basal + vmax*x^n/(k^n+x^n)`. Repression: `basal + vmax*k^n/(k^n+x^n)`. Michaelis–Menten: `vmax*x/(km+x)`. Here x, k and km are in counts; vmax and basal are events/second; n is dimensionless. These are complete effective propensities for a declared coarse-grained model, not proof that a mechanistic reduction is valid.

All reactions require sufficient reactants, including catalysts listed on both sides. If unavailable, the channel is disabled before evaluating a custom rate. Negative or non-finite evaluated rates, count overflow and inability to advance floating-point time produce errors, never clipping. Analytic derivatives for Hill laws need domain handling at zero, especially fractional n; M7 must not claim differentiability there unconditionally.

The expression tree belongs in M1b only if a concrete fixture needs it. It contains owned nodes, no strings executed as code and no graph references. Limits on depth and size are required. `min`/`max` introduce nonsmooth points; gradients must carry that qualification.

## 3. Contracts and identity

The Rust types and `describe` output are canonical. Checked-in schemas are regenerated and compared in CI; Markdown examples do not form a second schema source.

M0 includes model, run, reserved result and check-report contracts. Rates have a `kind` tag. Scalars are finite numeric literals or parameter names. Unknown fields are rejected recursively, including inside rate/method variants. Duplicate JSON keys are rejected before conversion to maps; duplicate species/reaction identifiers and cross-namespace collisions are separate errors. Diagnostics contain severity, stable code, RFC 6901 field pointer, message and repair hint.

M0 has no arbitrary observable registry or placeholder simulator. `describe` lists rate contracts and methods with explicit `contract_only`/unavailable execution status. A successful `check` means structural/domain validity only; it neither evaluates propensities nor proves boundedness, non-explosion, stability or biological correctness.

Run references use a portable relative model path, resolved from the run directory, and a lowercase SHA-256 of the exact source bytes. The CLI checks that the reference resolves to the supplied model below that directory. The in-memory API checks bytes and hash, but cannot check filesystem identity. Original byte snapshots are preserved in M2; any canonical form gets its own separately labelled hash.

Use explicit increasing save times, including 0 and the declared end time. This removes cumulative `save_every` rounding from the contract. At a save time, retain the state after events at or before that time. M2 tests the equality boundary and absorbing states.

Model population limits are safety stops. FSP bounds are inclusive truncation coordinates; omitted bounds are errors in the initial explicit-bounds contract. FSP's reachable-state cap does not need to exceed the full Cartesian box. A large box may have a small reachable subset due to conservation laws.

The reserved M0 result envelope is documentation and a generated type only. M2 will define its full evidence/archive version before emitting actual results. Format changes after this delivery get a new version and migration notes, even when adding variants.

## 4. Reproducibility, budgets and evidence

M2 uses a pinned ChaCha8 implementation with one stream per root-seed/replicate-index pair. Pin known-answer vectors and dependency versions. Rand's documentation distinguishes PRNG portability from floating-point distribution portability [S2]. A pure Rust math implementation and fixed reduction order are useful, but cross-platform identity remains a measured target.

Use a uniform strictly inside (0,1) for exponential waits and a separately specified [0,1) uniform for channel selection. Reject zero when mapping a 53-bit integer for the exponential sampler; allowing U=1 would permit a zero waiting time. Check total propensity, finite positive waits and `t+dt>t`. Do not hide dynamic-range loss by silently selecting the last channel.

Completed runs with fixed event/replicate work must have identical numerical trajectories and statistics at 1, 2 and 8 threads on the same binary/platform. Wall times, thread counts, paths and build identities are compared separately. Deadlines and cancellation change the completed subset; they are inherently outside this guarantee. Incomplete ensembles must not inherit an unqualified `exact_ensemble` assessment: completion can be state-dependent, so selecting only completed trajectories may bias statistics.

Apply event, state/population, wall-time, retained-output and memory/workspace limits. M0 checks configurations and a minimum numeric-output payload; it does not promise exact serialized size or execution memory. M2 must enforce live byte/work limits, preserve partial evidence, account for every requested replicate and provide bounded cooperative cancellation latency measured in time as well as events.

Common random numbers for optimization mean a declared fixed replicate seed set for all candidates. The same RNG stream is not automatically a strong coupling when topologies or reaction ordering change; measure the variance of paired objective differences. Freeze reaction order within each candidate representation and use disjoint validation seeds. Count validation evaluations separately.

Archives contain original inputs, numerical payloads, termination, assessment, checks, actual work, build/environment identity and hashes. Write a completion marker last. Verify path uniqueness, file hashes and completion state. Interrupted output must remain distinguishable from a complete archive. Reuse the sibling projects' mechanisms after inspecting their then-current implementation; do not copy an unverified protocol sketch.

## 5. Verification and performance policy

### Exact and structural tests

M0: strict parser, references, numeric domains, deterministic diagnostics, schemas, CLI behavior, file errors, and at least eight negative fixtures with expected codes and pointers. Include `1e400`, duplicate keys, wrong versions, unsupported rates/methods, oversized/deep documents, invalid counts, path/hash mismatches and incomplete FSP bounds. JSON cannot encode a cyclic tree; no cycle test is claimed.

M1: compare stoichiometry and dependencies against a simple independent oracle on 200 seeded small networks. Verify a full-rank conservation basis by exact multiplication and rational-span/rank comparison. Normalize only for display. Checked i128 overflow must yield a diagnostic or use a justified arbitrary-precision fallback; it cannot silently truncate the answer. Dependencies include reads by custom rates and reactant gates. Exact homodimer/heterodimer tests precede exporters.

### Statistical tests

Predeclare hypotheses, time points, marginal species, replicate counts, seeds and failure rules in `docs/decisions/statistics.md` before stochastic output exists. Fixed seeds make CI repeatable, not proof of correct coverage.

For analytic marginal CDF checks use DKW with Bonferroni family allocation, initially `n=100000` and family α=`0.001`. A union bound does not require independence between test points. Pure decay is binomial; immigration–death from zero is Poisson; with deterministic nonzero initial population it is a convolution of surviving initial particles and immigration. State the initial law in each fixture.

Do not apply a binomial exceedance-count model to correlated DSMTS observations from the same paths. Preserve upstream Z/Y diagnostics, then choose a separately justified suite rule (for example simultaneous Bonferroni thresholds with justified marginal approximations, or independently generated ensembles with a calibrated test). Handle zero expected variance explicitly. Rotating-seed checks report false-alarm rates; they do not retroactively loosen a failed gate.

Mean SE uses replicate variation. Variance, quantile, spectral and correlation uncertainty require their own stated estimators and assumptions. A single replicate does not produce an estimated ensemble SE. First-passage times that exceed the horizon are censored; report their probability and use censoring-aware summaries rather than discarding them. Statistical agreement with another simulator is evidence of consistency, not a certificate of correctness.

### FSP

Use reachable-state enumeration, an explicit absorbing sink for transitions out of the truncation and sparse uniformization. The original FSP result supplies a missing-mass bound under its assumptions [S3]. Document the non-explosion/model assumptions and matrix orientation. Avoid naive Poisson recursion from `exp(-Λt)` when Λt is large; subdivide time or use a stable centered-weight method.

Separate physical truncation leakage, discarded Poisson-series weight and numerical error. If missing mass is measured from an already truncated positive series, it includes that discarded series contribution; do not label it independently and add it twice. Ordinary f64 calculations are numerical error estimates, not a certified interval bound without roundoff analysis. Small distribution error does not by itself bound unbounded moments.

Use analytic fixtures, an independently implemented small-state checker, and later FiniteStateProjection.jl. Shared propensity code can hide the same error in SSA and FSP, so mutation tests must alter combinatorial factors and stoichiometric signs and be detected by an independent path. Compare analytic per-state values at an initial absolute target `1e-12`, recording the stated error accounting.

### Performance and application evidence

M4a measures Direct before adding algorithms. Freeze models, seed policy, duration and output density; report events, trajectories/second, completed work, startup and end-to-end cost separately. Use the same host, affinity/thread allowance, compiler settings and reference versions, paired repeated runs, and uncertainty. Avoid comparing rebop's specialized macro path with our dynamic path without labelling that difference.

Retain the original factor-of-two dynamic-API throughput goal as a workload-specific engineering target, not a forecast. Eight-thread efficiency of 0.8 is a target only where eight dedicated physical cores and sufficiently long work are available. CPU quotas, shared hosts and memory bandwidth must be recorded. Missing hardware is a missing measurement, not a failure or pass.

Select each extra aggregator from a measured bottleneck. Cost includes propensity evaluation, updates, selection and output. NRM/RSSA-CR cannot guarantee sublinear overall cost on densely coupled or unsuitable networks.

Oscillation criteria combine trajectory-wise autocorrelation and spectral structure, but finite noisy switches can mimic peaks. Split development and held-out topologies/parameters/seeds. Report false positives, false negatives, uncertainty, horizon/burn-in sensitivity and ambiguous classifications. The original 95%/99% goals are provisional benchmark targets until the fixture distribution is frozen, not universal properties.

Tau-leaping starts only after a workload establishes a need. Measure bias with confidence intervals against exact/FSP references, across decreasing ε. Require non-negative populations, failure accounting and an acceptable declared bias envelope; do not require sample estimates to move monotonically.

## 6. Milestones and exit gates

The week estimates in the original plan are replaced by concrete gates. Workload, numerical and integration discoveries determine iteration size; one numerical feature per reviewable change is preferred.

| Milestone | Deliverable | Exit evidence |
|---|---|---|
| **M0 / iteration 1** | Two-crate Rust workspace, strict JSON, typed model/run/reserved-result contracts, generated schemas, `describe`, `check`, complete examples, negative fixtures, CI and agent docs | Local fmt/clippy/release tests on installed stable and Rust 1.85/1.88; Python schema/CLI checker both normally and under `-O`; dependency gate; precise limitations |
| **M1a / iteration 2** | Compiled network indices, checked stoichiometry, reactant gates, dependencies, mass-action/Hill/MM rate evaluation, `analyze` | Exact count tests, 200 random networks, known enzyme conservation laws, independent oracle and mutation catches |
| **M1b** | Conservation/null-space completion and linkage classes if not in M1a; only needed expression nodes and reference exporters | Rank/residual checks, overflow cases, explicit rate convention exporter checks; no unnecessary symbolic differentiation framework |
| **M2** | Direct SSA, own pinned samplers, fixed-order ensembles, basic mean/variance, budgets/cancellation and verifiable archives | Analytic marginal tests; save-time/absorption fixtures; 1/2/8-thread numerical identity; interrupted/failed/oversized archive tests |
| **M3** | FSP, independently checked traces/statistics, selected DSMTS fixtures | Analytic FSP accuracy, leakage/error-accounting tests, predeclared statistical suite, checker mutation tests |
| **M4a** | Pinned rebop/Julia/GillesPy2 comparisons and Direct baseline | Reproducible same-host quality/cost table with failed and unavailable rows |
| **M5** | Minimal PyO3 binding, trajectory observables, three-gene circuit tutorial | Same JSON semantics; numerical CLI/binding equivalence; held-out behavioral tests; recorded end-to-end study |
| **M4b, as needed** | One of SortingDirect, NRM or RSSA-CR per iteration | Full M3 suite and demonstrable workload benefit before another method |
| **M6** | Existing command-job protocol, fcmaes worker, fitting and topology campaigns | Actual current SDK pins, `evaluate` recheck, fixed seed-set objective, held-out validation accounting, failure/cancellation evidence |
| **M7** | Optional diffsol ODE path, forward sensitivities, SQP fitting, cheap circuit screen | Analytic linear checks, nonlinear reference solutions, domain-aware derivatives, numerical gradient checks, measured feature MSRV |
| **M8** | Explicit SBML subset and expanded BioNetGen `.net` import | Supported cases pass; unsupported semantics receive coded errors; provenance retained |
| **M9, on demand** | Tau-leaping | Measured SSA bottleneck, non-negativity and declared bias/uncertainty study |
| **M10** | Full package, wheels and release documentation | Installed-package acceptance and clear supported-capability table; owner-controlled publishing |

M4a follows M3. M5 need not wait for M4b. If Direct misses the performance goal, profile and record the gap; a correct Python migration and observable study may still proceed. Source availability for the existing circuit study is an explicit M5 input; a synthetic substitute must be labelled as such.

## 7. Dependency and integration boundaries

Default library/CLI: edition 2024, MSRV 1.85, unsafe forbidden in project code. M0 uses only serde, serde_json, serde_path_to_error, schemars and sha2 plus transitive Rust dependencies. Lockfiles are shipped. The dependency gate rejects known native solver/build packages, link scripts with a Cargo `links` declaration, and unexpected additions pending review; it is not a proof that all dependencies are free of unsafe code or all assembly instructions.

RNG, math, ensemble and ODE dependencies enter only with their feature. Reconfirm diffsol's actual version/features/MSRV in M7; do not promise the earlier 1.89 claim as verified here. A BDF rewrite is a substantial numerical project, not a routine two-week fallback. Prefer a pinned working version or defer ODE work if the dependency tradeoff is unfavorable.

Integration crates live outside the core workspace with their own exact pins. Reuse the then-current job API and result schema, and keep job lifecycle distinct from numerical assessment. Nested parallelism must honor the supplied thread allowance. No current package name or remote repository availability is assumed from an earlier search.

Retain the two proposed M6 applications: `app.crn.fit-moments` for log-scaled rate fitting to mean/variance time courses, and `app.crn.circuit-score` for parameter tuning within a fixed topology. An incomplete simulation is a failed objective evaluation, not a deceptively good partial score. Integrate the worker's cancellation/deadline checkpoint in the event loop, and prove the completed objective identical for thread allowances 1 and 8. Planned campaigns compare seed-set sizes, aggregators and equal-budget BiteOpt retry versus surrogate search; preserve failures and count disjoint validation calls separately.

Retain Tutorial 2 in M7: global fcmaes exploration followed by SQP refinement from explicit best-point starts, with independently checked sensitivities. A deterministic circuit screen may prioritize topologies but must not reject them by itself: noise-driven oscillation can occur where an ODE converges to a fixed point. Predeclare a reduced SSA allocation for those topologies and include a Vilar-type noise-driven fixture to verify that it still reaches stochastic scoring. This screening claim needs its own measured false-negative analysis.

## 8. Iteration 1 acceptance and delivery

1. Create the new local project without modifying the attached source documents.
2. Implement contracts, generated schemas and structured discovery/checking.
3. Include immigration–death, an explicit enzyme mechanism, a reduced Michaelis–Menten rate model and the complete repressilator, with hash-bound runs.
4. Preserve expected diagnostics for broken models and runs. Add parser/resource/path/CLI regressions and an independent JSON Schema/CLI acceptance script.
5. Add MSRV/stable CI, deterministic schema snapshots, a dependency policy, licensing, agent guidance, model-card templates, statistical policy and release guidance.
6. Run and retain local evidence; distinguish locally tested checks from hosted CI configuration.
7. Deliver this revised plan, an implementation log with failures/fixes and limits, and a source ZIP including lockfile, docs, tests, fixtures and evidence. No build directories or cached dependencies in the ZIP.

M0 completion does **not** close a simulator performance gap: it has no solver to time. M1a is now implemented; remaining M1b analysis/exporter work precedes M2's executable exact SSA baseline.

## 9. Sources and provenance

The attached documents are inputs, not independently reproduced performance results. The review preserves their ecosystem roles and selected project name. No code is copied from GillesPy2 or other simulator implementations.

- **[S1]** Catalyst, official API documentation: [ODE and jump rate laws](https://docs.sciml.ai/Catalyst/stable/api/core_api/). Consulted 1 November 2026; documents factorial scaling versus binomial jump laws.
- **[S2]** Rust Rand Book: [Reproducibility](https://rust-random.github.io/book/crate-reprod.html). Consulted 1 November 2026; distinguishes portable generators, value-breaking changes and floating-point limitations.
- **[S3]** Munsky and Khammash (2006), [The finite state projection algorithm for the solution of the chemical master equation](https://doi.org/10.1063/1.2145882); [author-hosted paper](https://www.engr.colostate.edu/~munsky/Papers/JChemPhys_124_044104.pdf).
- **[S4]** [DSMTS source repository](https://github.com/darrenjw/dsmts). Pin the selected fixture revision, licenses and hashes in M3; no external data are vendored in M0.
- **[S5]** [Schemars documentation](https://docs.rs/schemars/0.8.22/schemars/). The implementation pins 0.8.22 and verifies emitted schemas, rather than assuming latest-version examples match its API.
- Algorithm references retained for later implementation: Gillespie (1977), DOI [10.1021/j100540a008](https://doi.org/10.1021/j100540a008); Gibson–Bruck (2000), DOI [10.1021/jp993732q](https://doi.org/10.1021/jp993732q); Cao–Gillespie–Petzold (2006), DOI [10.1063/1.2159468](https://doi.org/10.1063/1.2159468).
