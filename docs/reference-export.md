# GillesPy2 reference export

M1b exports model source for **GillesPy2 1.8.3**, with a runtime version guard. It is an optional reference tool; the Rust core does not import, link or execute GillesPy2.

```bash
cargo build --release --workspace --locked
./target/release/crn-craft export models/examples/enzyme-mechanism.model.json gillespy2 > export.json
python3 -c 'import json; from pathlib import Path; r=json.loads(Path("export.json").read_text()); Path("reference.py").write_text(r["artifact"]["source"])'
python3 -m venv .venv-reference
.venv-reference/bin/pip install -r validation/reference-requirements.txt
.venv-reference/bin/python -c 'from reference import build_model; print(build_model())'
```

The JSON envelope contains success, raw model SHA-256, diagnostics and an optional artifact. The artifact has generated Python source, target/version, original-to-target species and reaction aliases, and numerical limitations. A failure never emits a partial source file. `build_model()` constructs the model; it does not choose a solver, seed or time grid. Use GillesPy2's `NumPySSASolver` for the tested reference path.

## Supported semantics

The input must pass normal CRN-craft validation. Supported channels are mass action with total reactant order 0, 1 or 2, including catalysts and zero coefficients. Product stoichiometry retains its original integer values. Parameter references resolve to one frozen target rate parameter per channel. Species become S0, S1, …; parameters become K0, K1, …; reactions become R0, R1, …, avoiding target/Python reserved-name collisions.

| Source propensity | Target rate constant at volume 1 |
|---|---|
| c | c |
| c × A | c |
| c × A × B | c |
| c × A × (A−1) / 2 | c / 2 |

CRN-craft's coefficients already use stochastic counts. Target volume is deliberately fixed at one even when source Omega differs; Omega does not automatically rescale source coefficients. GillesPy2 1.8.3's elementary homodimer expression omits the factorial divisor, so its rate constant must be halved. This convention was checked in the installed package's `gillespy2/core/reaction.py` and against its actually generated propensity expressions. Package provenance is recorded with the validation evidence. The [GillesPy2 project](https://github.com/StochSS/GillesPy2) and [basic usage documentation](https://gillespy2.readthedocs.io/en/latest/getting_started/basic_usage/basic_usage.html) describe the target API.

Custom Hill/MM rates and reactant order above two are rejected with CRN-E050. A declared population safety limit is also rejected: the target model does not preserve CRN-craft's safety-stop behavior. Initial counts above 2^53 are rejected rather than silently rounded. Homodimer coefficients whose halving would underflow or round inexactly are rejected with CRN-E051. No model feature is silently removed.

## Numerical boundary

This is a stochastic-propensity translation, not an ODE exporter or a complete execution-contract translation. Target ODE coefficients, live event/output/time budgets, checked i64 updates, underflow diagnostics and cross-runtime bit identity are not preserved. The target solver uses floating-point populations; subsequent counts can leave its exact range even if initial counts are supported. Extreme target expressions can overflow or underflow differently from CRN-craft's scaled evaluator. Use representable reference fixtures and inspect target results; do not infer universal equivalence from successful source generation.

The exporter intentionally creates no run document or simulation result. Additional targets and custom-rate exports need separate reference fixtures and explicit convention/limit checks in M4a.

## Validation

`validation/check_export.py` compiles and imports the generated source, constructs actual pinned GillesPy2 models and checks initial counts, stoichiometry, alias mappings and target volume. It evaluates the package's own generated propensities for all 51×51 small-count pairs on a fixture covering every supported order, a catalyst, zero rates, parameter substitution, symbol collisions and non-unit Omega. It compares them with exact integer combinations and CRN-craft CLI evaluation. The two shipped elementary models are also checked.

A two-trajectory seeded immigration–death run is a reference SSA smoke test only. It establishes that generated source works with the tested target solver; it does not validate or benchmark CRN-craft simulation. Normal and Python `-O` reports must match byte-for-byte. See the iteration evidence for comparison counts and errors.
