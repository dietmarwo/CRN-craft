# Parameter optimization comparison: reproduction

Base: GitHub main `c16fcc5`, merged PR #8. This adds reusable parameter bindings
and optional model-specific propensity source generation. The fast policy is
inherited from the base; the checked default and old schemas remain available.

Read `docs/decisions/reusable-parameters-2026-11-02.md` for the predeclared
protocol and `docs/reusable-parameters.md` for the public Rust/Python API.
This directory is a separate optional benchmark workspace; rebop and optimizers
are not dependencies of the default numerical engine.

## Build

Use Rust 1.88 (the core is also tested on 1.85) and Python 3.12. Install the
benchmark dependencies in an isolated environment, build and install the wheel:

```sh
python -m pip install -r validation/requirements.txt -r validation/parameter-optimization/requirements.txt
maturin build --release --locked --out dist
python -m pip install --force-reinstall --no-deps dist/*.whl
python validation/generate_parameter_workers.py
cargo fmt --manifest-path validation/parameter-optimization/worker/Cargo.toml
python validation/parameter_optimization.py --build --target /tmp/crn-parameter-target
python validation/check_parameter_interface.py
python validation/validate_parameters.py --target /tmp/crn-parameter-target
```

The model/target definitions are already frozen. `parameter_optimization_models.py`
can recompute their independent Poisson/finite-state targets and rejects a
changed freeze rather than silently replacing the oracle. The initial
nonuniform-grid draft is retained as development evidence; before any
optimization measurement the six-point grids were made uniform to match
rebop's public Python API. Binding's horizon is 2, so its transient identifies
both rates. The other horizons are 10 and 20.

The `--build` command primes model-independent dependencies, then rebuilds only
the worker package separately for dynamic CRN, generated CRN and rebop macros.
It copies the three binaries into its target directory and records their hashes
and compilation times. These are dependency-cache-warm package builds using the
same standard Cargo release profile, not clean-machine dependency installation
or each repository's separately tuned release profile. Each source bundle
contains all three benchmark models; the generated CRN bundle also compiles an
arithmetic guard. Single-run totals charge the entire relevant bundle build.
Repeated-run totals amortize it once. Setup costs therefore depend on the
application packaging decision.

## Run complete optimization experiments

Finish every build and heavy check first, then run without concurrent workloads:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python validation/parameter_optimization.py --run
```

The run requires successful final gates and matching source hashes. It uses
fcmaes-rust 0.2.1 BiteOpt, three fixed optimizer seeds and 128 requested objective
calls per run. Actual counts and every point/loss are retained. All native CRN
routes use cached Direct, Ziggurat and `fast=True`, one numerical worker. The
comparison isolates fresh preparation, reusable binding and generated
propensities. Repeated preparations still use identical streams at each point.

Rebop 0.10.0 has three routes: a persistent expression model with runtime
parameters, mass-action reconstruction per candidate, and compiled macro code
with mutable coefficients. The two Python routes use sparse=True for the
33-species model and False for small models. Macro seeds use the same NumPy
integer-seed conversion as the Python API. All routes reduce the same six saved
states to unbiased sample variance and mean; API allocation/serialization and
IPC costs differ and remain included in complete optimizer time. Rebop and CRN
use different RNG streams; exact path identity is required only between the
three CRN routes.

Each candidate uses 32 paths with a fixed training seed. Final points are
reported against exact independent moments and against 4096 new paths on a
held-out seed. The normalized mean/variance objective excludes t=0. Compilation,
source generation, worker startup, all objectives, optimizer overhead, held-out
validation and total elapsed time are separate fields. Finite sample losses and
three search seeds do not establish an optimizer ranking or a universal speedup.

## Evidence

- `fixtures.json`, `models/`: frozen targets, bounds, guesses, seeds and hashes.
- `worker/`: optional locked Rust workspace and generated source modules.
- `evidence/generation.json`, `builds.json`: generation/compilation observations.
- `evidence/parameter-validation*.json`: normal/-O binding/AOT equality and
  extreme arithmetic, negative cases, concurrency and replay checks.
- `evidence/final/`: complete default/fast law, contract, legacy, MSRV and lint gates.
- `evidence/optimization.json`: full optimization traces, quality and costs.
- Development failures remain under explicitly named logs/JSON; no failed row
  is removed or converted into a performance pass.

No experiment in this directory reruns Julia. The attached prior fast-switch
log's projected Julia gap is not treated as a same-host result here.
