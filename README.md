# CRN-craft

![CRN-craft: chemical reaction networks in native Rust with stochastic simulation. Models, budgets, evidence: define, simulate, verify.](docs/images/crn-craft-banner.jpg)

> **One of 18 related repositories.** The [repository overview](https://github.com/dietmarwo/fcmaes-rust/blob/main/docs/repositories.md) shows what each one does and how they connect.

Native Rust simulation of well-mixed chemical reaction networks, designed for reproducible simulation and optimization workflows.

**Version 0.7.0.** Install the Python package with `pip install crn-craft` and the command line with `cargo install crn-craft-cli --locked`; the library is the crate `crn-craft` (`crn-craft-rules` remains as a deprecated re-export). One engine runs every model: fixed reaction networks and rule-based models whose species appear and disappear (mtDNA clones, clusters, polymers). Exact methods: Static Direct with four event selectors and `auto`, next reaction, the tree method, thinning for the ATP block and a certified mtDNA fast path. Approximate methods: tau and binomial leaping, and `hybrid` for continuous variables. Also included: finite state projection, conservation and structure analysis, ODE right-hand sides and Jacobians, quantile statistics, sealed result archives and compiled event loops. Optional integrations add Python and R workflows, diffsol BDF sensitivities, fcmaes job applications and SBML/BioNetGen import subsets. The default Rust engine has no Python, Julia or C/C++ dependency. A successful numerical run is not biological validation.

Start with the [tutorial](docs/tutorial.md), the [architecture overview](docs/developer/architecture.md) and `crn-craft describe`. Users of the mtDNA models start with the [mtDNA user manual](docs/mtdna-user-manual.md). Release 0.7.0 merged the former network and rules engines; the [changelog](CHANGELOG.md#migration) lists what changed and how to migrate.

## One engine

A model is a `crn-craft.model/1.0` document: types of molecules, optionally with numeric attributes; parameters; observables; and rules. Each rule names reactant patterns, a rate (or a whole-rule propensity) and products. A fixed reaction network is the special case of types without attributes. `crn-craft convert` turns a 0.6.0 network model (`crn-craft.model/0.1`) into this form, and the documents of the rules language (`crn-craft.rules-model/0.1`–`0.4`, `rules-run/0.1`–`0.6`) still run and replay their results byte for byte.

```json
{"name": "mutation",
 "reactants": [{"label": "x", "type": "mtDNA", "where": "x.wt == 1", "factor": "x.s * mu"}],
 "rate": "max(0, K - N) / D",
 "products": [{"type": "mtDNA", "copy": "x"},
              {"type": "mtDNA", "attributes": {"clone": "new_id()", "s": "choice(fracAdv, sm, sw)", "wt": 0}}]}
```

A request (`crn-craft.run/1.0`) names the method, the report times, the replicates and the budgets:

| Method | Kind | For |
| --- | --- | --- |
| `static_direct` | exact | closed models (a fixed set of species); selectors `linear`, `composition_rejection`, `sorting`, `rssa_cr` and `auto` (the default in `run/1.0`) |
| `next_reaction` | exact | closed models |
| tree (no `method`) | exact | every model, including species that appear and disappear |
| `mtdna_direct` | exact | the certified three-class mtDNA structure |
| `tau`, `binomial` | approximate | closed models |
| `hybrid` | approximate, declared tolerances | models with continuous variables |
| `thinning` | exact for its model class | models with an ATP block |
| `fsp` | distribution | closed models within declared bounds |

Every run draws from ChaCha8 streams keyed by the root seed and the replicate, so results do not depend on the thread count, and a stored model and request replay byte for byte. Results report mean, sample variance, standard errors and, on request, quantiles (`outputs.quantiles`) for every report time and observable; Python adds pandas and xarray views and `Result.save` to markdown or JSON.

## Optimization tutorials

The [four CRN optimization lessons](tutorials/optimization/README.md) cover reusable stochastic calibration, audited BDF sensitivities, oscillator Pareto/behavior searches, and bounded aggregation inference. Each has executable smoke/study presets, independent checks and held-out validation; the [shared catalog](https://github.com/dietmarwo/optimization-jobs/blob/main/tutorials/simulators/README.md) connects them to CFD and campaign lessons.

## Documentation

**Start here**

- [Tutorial](docs/tutorial.md): from a first model to exact ensembles, FSP, tau leaping and the Python interface, with screenshots and figures
- [Architecture overview](docs/developer/architecture.md): the crates, the one engine, contracts, parallelism, budgets, archives and generated loops
- [Technical report](docs/evidence/technical-report.md): numerical methods, validation evidence and measured performance
- [Comparison with GillesPy2](docs/evidence/comparison-gillespy2.md): same models, timed against GillesPy2's C++ and NumPy solvers
- [Event selectors](docs/alternative-ssa.md): the Static Direct selectors (`linear`, composition–rejection, sorting, RSSA-CR, `auto`) and how they were validated
- [Reusable parameters](docs/reusable-parameters.md): binding new parameter values to a prepared run, for optimization loops
- [Parameter-optimization comparison](docs/evidence/comparison-parameter-optimization.md): measured with 0.6.0 against rebop, with a pointer to the 0.7.0 worker rerun
- [Comparison with Catalyst and JumpProcesses](docs/evidence/comparison-catalyst-jumpprocesses.md) (measured with 0.6.0): seven workloads and five Julia SSA algorithms, with pinned environments and distribution checks; the [earlier comparison](docs/evidence/comparison-catalyst-jumpprocesses-before-acceleration.md) is preserved
- [Ten further Julia-gap iterations](docs/evidence/comparison-julia-gap10.md) (measured with 0.6.0): same-host rerun of all methods, exposed regressions and [inspected Julia mechanisms](docs/campaign/julia-gap10-source-investigation.md)
- [CRN-craft versus Julia at one and eight threads](docs/evidence/comparison-julia-8-threads.md) (measured with 0.6.0's `fast` policy): five JumpProcesses methods, 2,240 timed observations, scaling and replay checks
- [Agent context](docs/developer/ai-context.md): what an AI agent needs to use and extend CRN-craft
- [Model language](docs/reference/model-language.md): the `crn-craft.model/1.0` reference: types with attributes, patterns, rates and propensities, observables, continuous variables and higher-order rules, with [K1 law gates](validation/rules/RESULTS.md), [builder checks](validation/builder/RESULTS.md), [hybrid gates](validation/hybrid/RESULTS.md) and the [merge gates](validation/merge/RESULTS.md)
- [Using CRN-craft from R](docs/using-from-r.md): the Python package through reticulate, with a tested example in `examples/r/`
- [mtDNA user manual](docs/mtdna-user-manual.md), [optimization tutorial](docs/mtdna-optimization-tutorial.md) and [performance analysis](docs/evidence/mtdna-performance-analysis.md): simulating mtDNA mutant accumulation from Python, the command line and Rust, with runnable examples in `examples/mtdna/`, fitting with fcmaes, and what each method and option costs
- [Three-class mtDNA ensembles](docs/mtdna-rules.md): parallel Python ensembles, event trackers and callbacks, Static Direct, a specialized fixed-population sampler, next reaction, approximate tau and binomial steps, ATP-coupled thinning and the `hybrid` method for continuous variables ([review log](docs/campaign/mtdna-a0-a6-review-log-2026-11-03.md))
- [Specialized kernels and generated loops](docs/developer/specialized-kernels.md): when a model-specific kernel pays off, and compiled event loops

**Using CRN-craft**

- [Exact simulation](docs/reference/methods.md): Static Direct, termination states and budgets
- [Archives](docs/reference/archives.md): sealed `archive/1.0` result directories and the verifier
- [Contracts](docs/reference/contracts.md): models, requests, results and analysis reports, with the versions still read
- [Results and claims](docs/reference/results.md): what a result does and does not establish
- [Diagnostics](docs/reference/diagnostics.md): every error and warning code
- [Rust API](docs/reference/analysis.md): using the engine from Rust
- [Conservation cones and reaction graphs](docs/history/m1b-analysis.md): exact analysis and its limits
- [GillesPy2 reference export](docs/reference/gillespy2-export.md): exporting models for cross-checks
- [Extending the project](docs/developer/extending.md): how to add rate laws, methods and evidence

**Model and algorithm cards**

- [Card template](docs/evidence/model-cards/README.md)
- [Direct SSA card](docs/evidence/model-cards/direct.md)
- [Rate-law contracts](docs/evidence/model-cards/rate-laws.md)

**Campaign M3–M10 and performance iterations**

- [Campaign guide](docs/campaign/README.md) and [acceptance plan](docs/campaign/acceptance.md)
- [M3: finite state projection](docs/campaign/M3.md)
- [M4a: same-host Direct baseline](docs/campaign/M4a.md) and [M4b: dependency-cached Direct](docs/campaign/M4b.md)
- [M5: Python and circuit workflows](docs/campaign/M5.md)
- [M6: job SDK and optimization studies](docs/campaign/M6.md)
- [M7 targets](docs/campaign/M7-targets.md) and [M7: BDF, sensitivities, fitting and screening](docs/campaign/M7.md)
- [M8: import subsets](docs/campaign/M8.md)
- [M9 targets](docs/campaign/M9-targets.md) and [M9: approximate Poisson leaping](docs/campaign/M9.md)
- [M10: installation and release preparation](docs/campaign/M10.md)
- [Performance experiments P01–P10](docs/campaign/performance.md) and [P11–P20](docs/campaign/performance-p11-p20.md)
- [Campaign implementation log](docs/campaign/implementation-log-2026-11-01.md) and [review of 1 November 2026](docs/campaign/review-2026-11-01.md)
- [Exact-SSA acceleration log](docs/campaign/ssa-acceleration-implementation-log-2026-11-02.md) and [Julia-gap campaign log](docs/campaign/julia-gap-implementation-log-2026-11-02.md)
- [Ten further iterations implementation log](docs/campaign/julia-gap10-implementation-log-2026-11-02.md), [fixed protocol](docs/decisions/julia-gap10-2026-11-02.md) and [raw evidence](validation/julia-gap10/README.md)
- [The `fast` switch: implementation log](docs/campaign/fast-implementation-log-2026-11-02.md), with evidence in `validation/fast/`
- [Reusable parameters: implementation log](docs/campaign/reusable-parameters-implementation-log-2026-11-02.md) and [protocol](docs/decisions/reusable-parameters-2026-11-02.md)
- [Rules engine K0/K1: implementation log](docs/campaign/rules-k1-implementation-log-2026-11-02.md), [protocol](docs/decisions/rules-k0-k1-2026-11-02.md) and [results](validation/rules/RESULTS.md)
- [mtDNA A0–A6: implementation log](docs/campaign/mtdna-a0-a6-implementation-log-2026-11-03.md), [protocol](docs/decisions/mtdna-a0-a6-2026-11-03.md) and [review](docs/campaign/mtdna-a0-a6-review-log-2026-11-03.md)
- [Faster generic rules (S1–S6): implementation log](docs/campaign/fast-generic-rules-implementation-log-2026-11-03.md) and [review](docs/campaign/fast-generic-rules-review-log-2026-11-03.md)

**Plans and decisions**

- [Implementation plan, revision 2](docs/history/implementation-plan-2026-11-01.md)
- [M0 design decisions](docs/decisions/plan.md) and [statistical verification policy](docs/decisions/statistics.md)
- Acceptance declared before implementation: [iteration 2](docs/decisions/iteration-02-validation.md), [iteration 3](docs/decisions/iteration-03-validation.md), [iteration 4](docs/decisions/iteration-04-validation.md), [iteration 4 CME review check](docs/decisions/iteration-04-review-cme.md)
- Julia comparison and SSA acceleration protocols: [first Julia comparison](docs/decisions/julia-comparison-2026-11-02.md), [exact-SSA acceleration](docs/decisions/ssa-acceleration-2026-11-02.md), [closing the Julia gap](docs/decisions/julia-gap-2026-11-02.md)

**Implementation logs and validation results**

- Implementation logs: [iteration 1](docs/history/iteration-01-implementation-log.md), [iteration 2](docs/history/iteration-02-implementation-log.md), [iteration 3](docs/history/iteration-03-implementation-log.md), [iteration 4](docs/history/iteration-04-implementation-log.md)
- Validation results: [summary](validation/RESULTS.md), [iteration 2](validation/iteration-02/RESULTS.md), [iteration 3](validation/iteration-03/RESULTS.md), [iteration 4](validation/iteration-04/RESULTS.md)

**Project**

- [Changelog](CHANGELOG.md), [releasing](RELEASING.md) and [agent working rules](AGENTS.md)

## Quick start

From Python (`pip install crn-craft`):

```python
import crn_craft as crn

model = crn.load_example("immigration-death")
result = crn.ensemble(
    model,
    end=30.0,
    report_times=[0.0, 10.0, 20.0, 30.0],
    replicates=400,
    seed=42,
    method="static_direct",
).simulate()
print(result.statistics["mean"])
print(result.statistics["median"])
result.save("immigration-death.md")
```

A model can also be built in code. A numeric rate is a mass-action
constant; a string is the whole propensity:

```python
import crn_craft as crn

model = crn.Model("dimerization")
model.add_species({"A": 100, "A2": 0})
model.add_parameter({"k_off": 0.1})
model.add_reaction(0.001, ["A", "A"], ["A2"])
model.add_reaction("k_off * A2", ["A2"], ["A", "A"])
result = crn.ensemble(
    model, end=50.0, report_times=[0.0, 25.0, 50.0], replicates=200, seed=1
).simulate()
print(result.statistics["mean"])
```

The command line (`cargo install crn-craft-cli --locked`) reads the same documents:

```bash
crn-craft describe
crn-craft check models/examples/immigration-death.model.json
crn-craft analyze models/examples/enzyme-mechanism.model.json
crn-craft convert models/legacy/repressilator.model.json > repressilator.model.json
crn-craft export models/examples/enzyme-mechanism.model.json gillespy2
```

`crn-craft run MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]` runs a request (`crn-craft.run/1.0`), and `crn-craft verify DIRECTORY` checks an archive; the [tutorial](docs/tutorial.md#3-checking-inputs-and-reading-diagnostics) writes a request step by step. Every command prints JSON and exits 0 on success and 1 for invalid input. Diagnostics name a code, a JSON path and a hint:

```json
{
  "severity": "error",
  "code": "CRN-R012",
  "path": "/rules/0/rate",
  "message": "unknown name `k_missing`; expected a parameter or observable in `k_missing`",
  "hint": "See docs/rules-language.md for the names and functions each field may use."
}
```

## What is included

| Component | Implemented now |
|---|---|
| Rust library and CLI | Strict JSON parsing, typed contracts, semantic checks, generated schemas, one engine in the crate `crn-craft` |
| Model language | Types with numeric attributes, patterns with conditions and local factors, rates and whole-rule propensities, observables, continuous variables, higher-order rules, stop conditions; network models through `convert` |
| Exact methods | Static Direct with four selectors and `auto`, next reaction, the tree method for open species sets, the certified mtDNA fast path, thinning |
| Approximate methods | Tau and binomial leaping with critical-channel fallback; `hybrid` with declared tolerances |
| Distributions | FSP: reachable sparse projection, escaping sink, stable split uniformization, separate error accounting |
| Analysis | Exact conservation laws and cones within explicit limits, linkage classes and deficiency, ODE right-hand sides and Jacobians |
| Results | Mean, sample variance, standard errors and quantiles; kept paths, event logs and trackers; sealed `archive/1.0` directories with a verifier |
| Reproducibility | ChaCha8 streams keyed per replicate, results independent of the thread count, byte-identical replay of stored documents |
| Speed | Compiled event loops (opt-in, local Rust toolchain), byte-identical to the interpreter |
| Python and R | Builder, pandas and xarray results, `Result.save`, importers for SBML and BioNetGen subsets, R through reticulate |
| Optional integrations | Native fcmaes worker and job applications; diffsol BDF and sensitivities; global-to-SQP fitting |
| Evidence | Law gates (4,814 comparisons over twelve models for every selector), exact master-equation and analytic campaigns, a 47-mutant source mutation gate, DSMTS FSP references, the frozen 0.6.0 manual run in CI, and comparisons with GillesPy2 |

Unknown fields and duplicate JSON keys are errors. Names are unique ASCII identifiers. Counts are checked integers. A request names its model by path and by the SHA-256 of the model's exact bytes, so reformatting a model changes its hash; the repository's `.gitattributes` keeps LF line endings, so the hashes also hold with Git for Windows' `core.autocrlf=true`.

## Release history

Release 0.3.0 answers the first user feedback on the Python interface
([plan](docs/history/plan-0.3.0-feedback.md)). rules-model/0.3 lets a rule give its
whole `propensity` instead of a per-tuple `rate`, and declares `continuous`
variables: ODEs between events that rates read, simulated by the
approximate `hybrid` method at declared tolerances (rules-run/0.5,
rules-result/0.3; Dormand–Prince 5(4) with integrated hazards). The ATP
model is now an ordinary model file
(`models/rules/mtdna-three-class-atp-continuous.rules.json`); the legacy
`atp` block and `thinning` keep working. In Python, `crn_craft.rules.Model`
builds models from `add_reaction(rate, reactants, products)` calls, results
are pandas DataFrames (`keep_paths="all"`, `to_xarray()`), and event logs
reserve output from the model's own sizes, so the reported 50-entry
mutation log fits the default budget. See the
[mtDNA user manual](docs/mtdna-user-manual.md) and the
[changelog](CHANGELOG.md) for migration notes.

Release 0.4.0 makes Ziggurat waiting times and sparse count trackers the
Static Direct defaults in rules-run/0.6, 1.17× faster per event on the
three-class model, and Python's `method="static_direct"` selects it. The
new default uses the random stream differently;
`{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False}`
reproduces the numbers of 0.3.0, and stored requests replay unchanged
([Static Direct options](docs/mtdna-rules.md#generic-static-direct-options)).
Opt-in, `rules.compile_loop` compiles a model's Static Direct event loop
with a local Rust toolchain into a cached program whose results are
byte-identical to the interpreter's, about 2.2× faster on that model
([manual](docs/mtdna-user-manual.md#optional-compiled-event-loop)).

Release 0.5.0 compiles string rates too: the generated loop now covers
every model Static Direct runs (expression rates, the builder's whole-rule
propensities, `distinct`, `max` and `min` observables), about 2× faster on
the builder's SIR model with byte-identical results. rules-model/0.4 marks
continuous variables `independent` when they read only parameters and
other independent variables; `hybrid` integrates them once per run and
shares the solution with every replicate (`Model.add_ode(...,
independent=True)`, 2.8× on an eight-compartment forcing). In Python,
`add_species` and `add_parameter` take mappings, and kept paths and event
logs are xarray DataArrays (`paths_frame()` and `event_frame()` return the
old DataFrames); see the [changelog](CHANGELOG.md) for migration.

Release 0.6.0 orders kept paths like event logs, replicate first
(`result.paths` has the dimensions replicate × time × variable), adds
`Model.set_stop` for stop conditions in builder models, and documents
[using CRN-craft from R](docs/using-from-r.md) through reticulate.

Release 0.7.0 merges the network engine and the rules engine into one
engine in the crate `crn-craft` ([plan](docs/history/plan-0.7.0-engine-merge.md),
[implementation log](docs/campaign/merge-0.7/README.md)). Every rules
document of 0.6.0 replays byte for byte; network models convert to
`crn-craft.model/1.0`; quantiles, `Result.save`, selectors, higher-order
rules and larger limits are new. See the [changelog](CHANGELOG.md) for the
migration.


## Developer checks

CI (`.github/workflows/ci.yml`) is the complete list. The core of it:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
cargo build --release --locked -p crn-craft-cli
python3 -m pip install -r validation/requirements.txt
python3 -m ruff format --check validation python
python3 -m ruff check validation python
python3 validation/check_contracts.py
python3 validation/check_ssa.py --output ssa.json
python3 -O validation/check_ssa.py --output ssa-optimized.json
diff ssa.json ssa-optimized.json
python3 validation/check_cme.py --output cme.json
python3 validation/check_rules.py --output rules.json
python3 validation/check_mtdna.py --binary target/release/crn-craft --output mtdna.json
python3 validation/check_links.py
# After installing the wheel (maturin build --release --locked --out dist; pip install dist/*.whl):
python3 validation/check_python.py --binary target/release/crn-craft --output python.json
python3 validation/check_rules_python.py --output rules-python.json
python3 validation/check_manual.py --binary target/release/crn-craft --output manual.json
python3 validation/check_manual.py --binary target/release/crn-craft --manual validation/compat/manual-0.6.0.md --output manual-0.6.0.json
python3 validation/check_tutorial.py --binary target/release/crn-craft --output tutorial.json
(cd /tmp && python3 -I "$OLDPWD/validation/check_installed.py" --output installed.json)
python3 validation/mutation_check.py
```

Each Python check also runs under `python -O`, and CI compares the two reports. The engine-merge gates (`validation/merge/check_laws.py`, the speed gate and the worker comparison) are run by hand and recorded in [validation/merge/RESULTS.md](validation/merge/RESULTS.md). GitHub Actions runs fmt, clippy and release tests on Rust 1.85.0, 1.88.0 and stable, and the installed-wheel acceptance on Linux, Windows and macOS. The default build has no C/C++ simulator, Julia or Python runtime dependency; it pins pure-Rust libm for rate evaluation.

## Campaign and remaining validation

The functional M3–M10 work and ten performance iterations are documented in
[docs/campaign](docs/campaign/README.md), with predeclared gates and retained failures
under validation/campaign. Optional integrations have independent lockfiles/MSRVs.
Release wheels pass the installed-package acceptance on Linux, Windows and macOS
([RELEASING.md](RELEASING.md)); a Julia reference run and broad biological
behavioral validation remain owner/platform-specific work. The circuit demonstrations
are labelled synthetic; detector diagnostics do not establish universal sensitivity.

The ecosystem boundary remains: CRN-craft supplies numerical models; fcmaes searches
parameters; optimization-jobs executes jobs; experiment-craft coordinates campaigns.

Read [the agent context](docs/developer/ai-context.md), [contract details](docs/reference/contracts.md), [diagnostics](docs/reference/diagnostics.md), [archive contract](docs/reference/archives.md) and [release instructions](RELEASING.md) before extending the project. License: MIT; dependencies retain their respective licenses.
