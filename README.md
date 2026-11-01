# CRN-craft

Native Rust foundations for well-mixed chemical reaction networks, designed for reproducible simulation and optimization workflows.

**Current status: remaining-milestone campaign (`0.2.0-dev.1`).** Native Direct and exact dependency-cached SSA, sparse FSP and explicitly approximate Poisson leaps execute. Optional Python bindings, trajectory diagnostics, pinned fcmaes job applications, diffsol BDF/forward sensitivities, global-to-SQP fitting and explicit SBML/BioNetGen import subsets are included. The default Rust engine remains independent of Python, jobs and ODE dependencies. A successful numerical run is not biological validation.

Start with [the campaign guide](docs/campaign/README.md), [measured milestone results](docs/campaign/) and `crn-craft describe`. Legacy input/result/archive contracts remain available; new approximate/FSP results use distinct schemas.

## Quick start

Requires Rust 1.85 or newer. From this directory:

```bash
cargo build --release --workspace --locked
./target/release/crn-craft describe
./target/release/crn-craft describe model
./target/release/crn-craft check models/examples/immigration-death.model.json models/examples/immigration-death.run.json
./target/release/crn-craft analyze models/examples/enzyme-mechanism.model.json
./target/release/crn-craft analyze models/examples/immigration-death.model.json models/examples/immigration-death.state.json
./target/release/crn-craft analyze models/analysis/cone-gap.model.json
./target/release/crn-craft export models/examples/enzyme-mechanism.model.json gillespy2
./target/release/crn-craft simulate models/examples/immigration-death.model.json models/examples/immigration-death.run.json /tmp/crn-ssa-example
./target/release/crn-craft verify /tmp/crn-ssa-example
cargo test --release --workspace --locked
```

`check` emits JSON with `valid`, `execution_available: false`, hashes and diagnostics. `analyze MODEL [STATE]` emits `crn-craft.analysis/0.2` with indexed stoichiometry, evaluated rates, conservation cones, linkage classes and structural deficiency. `analyze-v1` retains the 0.1 report. `export MODEL gillespy2` returns generated Python in a JSON report; it does not execute that source. Without STATE it uses initial counts; otherwise STATE must name every species. These commands exit **0** on success and **1** for invalid inputs, usage, I/O or numerical errors. Warnings do not change exit 0. Conservation size/overflow limits are explicit partial-analysis warnings. `simulate MODEL RUN NEW_DIRECTORY [--cancel-file PATH]` writes result 0.2 and archive 0.1. It exits 0 only for a completed ensemble; incomplete work exits 1 with evidence when finalization succeeds. `verify DIRECTORY` exits 0 for a valid archive even if the simulation is incomplete. The output directory must be new. Read the [Direct SSA guide](docs/direct-ssa.md) and [archive contract](docs/archive-contract.md). Historical preflight/analysis flags remain unchanged; use `describe` for executable capabilities.

```json
{
  "severity": "error",
  "code": "CRN-E014",
  "path": "/reactions/0/rate/coefficient",
  "message": "unknown parameter \"missing\"",
  "hint": "Declare the parameter or provide a numeric literal."
}
```

## What is included

| Component | Implemented now |
|---|---|
| Rust library and CLI | Strict JSON parsing, typed contracts, semantic preflight checks, generated schemas |
| Network core | Immutable indices, sparse stoichiometry, availability gates, dependencies, transactional checked updates |
| Rate evaluation | Stochastic mass action, Hill activation/repression and Michaelis–Menten; finite-range diagnostics |
| Analyzer | Exact linear basis and nonnegative cone within explicit limits; complete conservation bounds when the cone succeeds; linkage/SCC/terminal classes and deficiency |
| Reference exporter | GillesPy2 1.8.3 elementary mass action, explicit homodimer conversion, fixed count convention and safe aliases |
| Models | Immigration–death; explicit enzyme mechanism; reduced Michaelis–Menten; complete three-gene repressilator |
| Simulator | Direct SSA; pinned ChaCha8 streams; deterministic ordered ensembles; mean/sample variance/mean SE; event, population, wall-time and output limits; cooperative cancellation |
| Results | Exact input snapshots, complete replicate accounting, checksummed archives and a marker written last; independent verifier |
| FSP | Reachable sparse projection, escaping sink, stable split uniformization, separate leakage/series/roundoff accounting |
| Experimental approximation | Checked Poisson leaps, critical Direct fallback, own approximate result schema |
| Optional applications | Python diagnostics/importers; native fcmaes worker; diffsol BDF and sensitivities; global-to-SQP fitting |
| Evidence | Rust regression tests, 13 invalid fixtures, independent exact/Decimal oracle for 200 seeded networks, 300-network exact cone/graph oracle with five targeted models, 43 source mutations against a passing baseline, three 100,000-path analytic campaigns, four 100,000-path exact master-equation checks of second-order, Hill and Michaelis–Menten channels, thread identity and real GillesPy2 exporter checks, schema and dependency gates |
| Documentation | Revised plan, diagnostics, rate conventions, statistical policy, extension/model-card templates and agent guide |

Unknown fields and duplicate JSON keys are errors. Names are unique ASCII identifiers. Counts are checked i64 values. Rate parameters have finite domains. Custom rates are complete propensities, with separate reactant-availability gates. No extra mass-action factor is applied to a custom law.

`CompiledNetwork` exposes pure rate evaluation and an atomic `fire_reaction` operation without a clock or random generator. Conservation arithmetic uses checked i128; unavailable results never imply absence of laws. Derived bounds describe the **initial conservation class**, even when rates are evaluated at another state. Read [the network API](docs/network-api.md), [M1b analysis limits](docs/m1b-analysis.md) and [reference export guide](docs/reference-export.md).

`run.model.sha256` identifies the model's exact bytes. Model paths are relative to the **run file directory**, must stay below it, and must resolve to the same model supplied to `check`. Reformatting a model changes its hash. Recompute after reviewing an intentional edit:

```bash
sha256sum models/examples/immigration-death.model.json
```

Copy the returned digest into the corresponding run file. Do not assume that semantically equivalent JSON has the same byte identity. The repository's `.gitattributes` keeps LF line endings in every checkout, so the example hashes also hold with Git for Windows' `core.autocrlf=true`.

## Developer checks

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
python3 -m pip install -r validation/requirements.txt
python3 -m ruff format --check validation
python3 -m ruff check validation
python3 validation/check_dependencies.py
python3 validation/check_contracts.py
python3 -O validation/check_contracts.py
python3 validation/check_network.py --output network.json
python3 -O validation/check_network.py --output network-optimized.json
diff network.json network-optimized.json
python3 validation/check_structure.py --output structure.json
python3 -O validation/check_structure.py --output structure-optimized.json
diff structure.json structure-optimized.json
python3 validation/reference_rng.py
python3 validation/check_ssa.py --output ssa.json
python3 -O validation/check_ssa.py --output ssa-optimized.json
diff ssa.json ssa-optimized.json
python3 validation/check_cme.py --output cme.json
python3 -O validation/check_cme.py --output cme-optimized.json
diff cme.json cme-optimized.json
python3 validation/mutation_check.py
```

The optional real-reference checks use a separate Python environment:

```bash
python3 -m venv .venv-reference
.venv-reference/bin/pip install -r validation/reference-requirements.txt
.venv-reference/bin/python validation/check_export.py
.venv-reference/bin/python -O validation/check_export.py
```

The eleven original schema snapshots are preserved, and three new native schemas are compared with the CLI. Model/run/check/state, analysis 0.2, legacy analysis 0.1 and export 0.1 remain unchanged. Reserved result 0.1 is retained as `result-v1`; real simulation results use 0.2. Archive and archive-check are 0.1; discovery is 0.5. Strict parsers and semantic checks remain mandatory. The optional exporter checker runs a short GillesPy2 smoke test; the independent SSA campaign validates CRN-craft against declared analytic distributions, and `check_cme.py` against exact master-equation laws for nonlinear and custom rates.

GitHub Actions is configured for fmt, clippy and release tests on Rust 1.85.0, 1.88.0 and stable; dependency, schema/CLI and independent network/structure/SSA/CME checks normally and under Python `-O` with identical reports; source mutation checks; and contract/network/structure/SSA/archive tests with `core.autocrlf=true`. See [iteration 4 evidence](validation/iteration-04/RESULTS.md) for executed local checks. The default build has no C/C++ simulator, Julia or Python runtime dependency. It pins pure-Rust libm for rate evaluation. Python is only needed for independent acceptance checks.

## Campaign and remaining validation

The functional M3–M10 work and ten performance iterations are documented in
[docs/campaign](docs/campaign/README.md), with predeclared gates and retained failures
under validation/campaign. Optional integrations have independent lockfiles/MSRVs.
Publishing, multi-platform wheel acceptance, a Julia reference run and broad biological
behavioral validation remain owner/platform-specific work. The circuit demonstrations
are labelled synthetic; detector diagnostics do not establish universal sensitivity.

The ecosystem boundary remains: CRN-craft supplies numerical models; fcmaes searches
parameters; optimization-jobs executes jobs; experiment-craft coordinates campaigns.

Read [the agent context](docs/ai-context.md), [contract details](docs/contracts.md), [diagnostics](docs/diagnostics.md), [archive contract](docs/archive-contract.md) and [release instructions](RELEASING.md) before extending the project. License: MIT; dependencies retain their respective licenses.
