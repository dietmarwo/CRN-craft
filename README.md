# CRN-craft

Native Rust foundations for well-mixed chemical reaction networks, designed for reproducible simulation and optimization workflows.

**Current status: iteration 3 / M1b (`0.1.0-dev.3`).** This version validates model/run documents, compiles indexed networks, evaluates propensities and analyzes conservation cones and reaction graphs. It also exports elementary mass-action models to pinned GillesPy2 1.8.3. It does not simulate, optimize or import SBML. A successful analysis is not a biological validation result.

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
cargo test --release --workspace --locked
```

`check` emits JSON with `valid`, `execution_available: false`, hashes and diagnostics. `analyze MODEL [STATE]` emits `crn-craft.analysis/0.2` with indexed stoichiometry, evaluated rates, conservation cones, linkage classes and structural deficiency. `analyze-v1` retains the 0.1 report. `export MODEL gillespy2` returns generated Python in a JSON report; it does not execute that source. Without STATE it uses initial counts; otherwise STATE must name every species. These commands exit **0** on success and **1** for invalid inputs, usage, I/O or numerical errors. Warnings do not change exit 0. Conservation size/overflow limits are explicit partial-analysis warnings. No command emits a simulation result or exit code 2.

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
| Run contracts | Direct SSA requests and bounded FSP requests; execution unavailable |
| Evidence | Rust regression tests, 13 invalid fixtures, independent exact/Decimal oracle for 200 seeded networks, 300-network exact cone/graph oracle with five targeted models, 17 source mutations against a passing baseline and real GillesPy2 exporter checks, schema and dependency gates |
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
python3 validation/mutation_check.py
```

The optional real-reference checks use a separate Python environment:

```bash
python3 -m venv .venv-reference
.venv-reference/bin/pip install -r validation/reference-requirements.txt
.venv-reference/bin/python validation/check_export.py
.venv-reference/bin/python -O validation/check_export.py
```

All eight schema snapshots are compared with the CLI. Model/run/check/result/state 0.1 remain unchanged. The original analysis 0.1 schema is retained as `analysis-v1`; current analysis is 0.2, export is 0.1 and discovery is 0.3. JSON Schema validates shape and selected local bounds; Rust checks cross-references, positivity, cardinalities, hash/path identity and other semantic constraints. Core checks evaluate rates and structure. The optional exporter checker also runs a short GillesPy2 SSA smoke test; this does not implement or benchmark a CRN-craft simulator.

GitHub Actions is configured for fmt, clippy and release tests on Rust 1.85.0, 1.88.0 and stable; dependency, schema/CLI and independent network/structure checks normally and under Python `-O` with identical reports; source mutation checks; and contract/network/structure tests with `core.autocrlf=true`. See [iteration 3 evidence](validation/iteration-03/RESULTS.md) for executed local checks. The default build has no C/C++ simulator, Julia or Python runtime dependency. It pins pure-Rust libm for rate evaluation. Python is only needed for independent acceptance checks.

## Next work

[The revised plan](docs/implementation-plan-2026-11-01.md) tracks completed M1a/M1b analysis work. Broader exporters and fixture-driven expressions remain optional extensions. **Next is M2**, which adds the first executable Direct SSA baseline. Performance comparisons begin once a verified simulator exists.

The planned ecosystem boundary remains: CRN-craft supplies numerical models; fcmaes searches parameters; optimization-jobs executes jobs; experiment-craft coordinates campaigns. Their integrations are future milestones and are not dependencies of the core project.

Read [the agent context](docs/ai-context.md), [contract details](docs/contracts.md), [diagnostics](docs/diagnostics.md), [result contract](docs/result-contract.md) and [release instructions](RELEASING.md) before extending the project. License: MIT; dependency notices remain under their respective licenses.
