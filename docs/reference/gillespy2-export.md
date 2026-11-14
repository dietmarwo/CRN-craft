# GillesPy2 reference export

M1b exports model source for **GillesPy2 1.8.3**, with a runtime version guard. It is an optional reference tool; the Rust core does not import, link or execute GillesPy2.

```text
cargo build --release --workspace --locked
./target/release/crn-craft export models/examples/enzyme-mechanism.model.json gillespy2 > export.json
python3 -c 'import json; from pathlib import Path; r=json.loads(Path("export.json").read_text()); Path("reference.py").write_text(r["artifact"]["source"])'
python3 -m venv .venv-reference
.venv-reference/bin/pip install -r validation/reference-requirements.txt
.venv-reference/bin/python -c 'from reference import build_model; print(build_model())'
```

The JSON envelope contains success, raw model SHA-256, diagnostics and an optional artifact. The artifact has generated Python source, target/version, original-to-target species and reaction aliases, and numerical limitations. A failure never emits a partial source file. `build_model()` constructs the model; it does not choose a solver, seed or time grid. Use GillesPy2's `NumPySSASolver` for the tested reference path.

## Supported semantics

The input is any model version that has a flat view (a closed model; `CRN-R017` otherwise): `crn-craft.model/1.0`, a `rules-model`, or a network `crn-craft.model/0.1` through the converter. Each species label becomes one target species and each channel one reaction. A channel whose rule has a constant `rate` (reading no observable) and takes at most two molecules is GillesPy2 mass action, including catalysts and zero coefficients; the constant, times the pattern `factor`, is one frozen target rate parameter per channel. Any other channel (a rate that reads an observable, a whole-rule `propensity`, more than two molecules) is a custom propensity string: the rule's expression with parameter values inlined, each observable as the weighted sum of its species and `^` as `**`, times the falling factorial n(n − 1)…/k! per reactant for a `rate` rule, or times one availability factor `x/(x + 1e-300)` per molecule taken for a whole-rule propensity. Product stoichiometry retains its integer values. Species become S0, S1, …; parameters become K0, K1, …; reactions become R0, R1, …, avoiding target/Python reserved-name collisions; the aliases name the rules (`rule[j]` for the channels of a rule with several).

| Source propensity | Target rate constant at volume 1 |
|---|---|
| c | c |
| c × A | c |
| c × A × B | c |
| c × A × (A−1) / 2 | c / 2 |

CRN-craft's coefficients already use stochastic counts. Target volume is deliberately fixed at one even when source Omega differs; Omega does not automatically rescale source coefficients. GillesPy2 1.8.3's elementary homodimer expression omits the factorial divisor, so its rate constant must be halved. This convention was checked in the installed package's `gillespy2/core/reaction.py` and against its actually generated propensity expressions. Package provenance is recorded with the validation evidence. The [GillesPy2 project](https://github.com/StochSS/GillesPy2) and [basic usage documentation](https://gillespy2.readthedocs.io/en/latest/getting_started/basic_usage/basic_usage.html) describe the target API.

Rejected with CRN-E050, without a partial source: comparisons, logic and functions in an expression (GillesPy2's Python solvers have no `exp`, its C++ solver no `min`, `max` or `abs`, and its expression printer drops comparisons), observables reduced other than by `sum`, more than 64 molecules of one species in a custom propensity, a negative constant rate, and a `stop` condition, which GillesPy2 cannot apply. Population limits of a network model are refused by the converter (`CRN-R016`). Homodimer coefficients whose halving would underflow or round inexactly are rejected with CRN-E051. No model feature is silently removed. Until release 0.6.0 the exporter refused the Hill and Michaelis–Menten laws; since 0.7.0 they are custom propensities.

## Numerical boundary

This is a stochastic-propensity translation, not an ODE exporter or a complete execution-contract translation. Target ODE coefficients, live event/output/time budgets, checked i64 updates, underflow diagnostics and cross-runtime bit identity are not preserved. The target solver uses floating-point populations; subsequent counts can leave its exact range even if initial counts are supported. Extreme target expressions can overflow or underflow differently from CRN-craft's scaled evaluator. Use representable reference fixtures and inspect target results; do not infer universal equivalence from successful source generation.

The exporter intentionally creates no run document or simulation result. Additional targets need separate reference fixtures and explicit convention and limit checks.

## Validation

`validation/check_export.py` compiles and imports the generated source, constructs actual pinned GillesPy2 models and checks initial counts, stoichiometry, alias mappings and target volume. It evaluates the package's own generated propensities for all 51×51 small-count pairs on a fixture covering every supported order, a catalyst, zero rates, parameter substitution, symbol collisions and non-unit Omega. It compares them with exact integer combinations and CRN-craft CLI evaluation. The shipped example models are also checked, including the custom propensities of the Michaelis–Menten and repressilator examples against the exact rate laws and the engine's rates, and a rate with a function and a stop condition must be refused with CRN-E050.

A two-trajectory seeded immigration–death run is a reference SSA smoke test only. It establishes that generated source works with the tested target solver; it does not validate or benchmark CRN-craft simulation. Normal and Python `-O` reports must match byte-for-byte. See the iteration evidence for comparison counts and errors.
