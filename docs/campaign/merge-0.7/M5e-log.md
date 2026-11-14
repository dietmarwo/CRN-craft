# M5e log: ODE right-hand side and Jacobians on the flat view

Task M5e of [the plan](../../history/plan-0.7.0-engine-merge.md): `kinetics.rs`
takes a `FlatView`; a new `rules/dual.rs` evaluates rule expressions with
dual numbers; the right-hand side and both Jacobians come from them;
`integrations/diffsol` and `validation/check_ode.py` keep working.

## What changed

- **`crates/crn-craft/src/kinetics.rs`** rewritten on the flat view.
  `Kinetics::new(&FlatView)` reads the view's species, channels (rule,
  reactants, net changes, pattern factor, whole-rule flag), rule
  expressions, observables (reduction and weighted terms) and parameters.
  `evaluate(y, p)` computes, per channel, the propensity at real counts:
  - a rule with a `rate`: rate × factor × Π yᵢ^kᵢ / kᵢ!, the macroscopic
    form documented in 0.6.0, computed as the same running products of
    y / j (so mass action of order 0 and 1 is bit for bit the 0.6.0 value);
  - a rule with a whole `propensity`: the expression × factor, without the
    availability gate (see the decisions below);
  and accumulates net change × propensity into the right-hand side, and
  net change × derivatives into the state Jacobian and the parameter
  Jacobian (parameters in sorted order, as before). Derivatives come from
  the dual parts: with respect to a parameter directly, with respect to a
  species through the observables that read it (sum: its weight; distinct,
  max, min: 0, since they are piecewise constant in the counts) and through
  the combination term. `Derivatives` is unchanged.
- **`crates/crn-craft/src/rules/dual.rs`** (new, crate-internal):
  `Program::new(&Expr)` lists the parameters and observables an expression
  reads; `Program::evaluate(params, observables)` returns the value
  (computed with `expr.rs`'s own `eval_binary` and `deterministic_call`, so
  bit for bit `Expr::eval`), the gradient over those inputs and a flag for
  derivatives that are singular at 0; `Program::vanishes` decides from the
  expression's structure whether it is 0 whenever given observables are 0.
  Every operator and function of `expr.rs` is covered: `+ - * / ^`, unary
  `-` and `!`, comparisons, `&&`/`||` (short-circuit as `eval`), `if`,
  `min`, `max`, `exp`, `log`, `sqrt`, `abs`, `floor`, `ceil`; random
  functions cannot occur in rates and give NaN.
- **`rules/expr.rs`**: visibility only (`Node`, `Op`, `Func`, `Expr::root`,
  `truth`, `flag`, `eval_binary`, `deterministic_call` are `pub(super)`).
  **`rules/mod.rs`**: `pub(crate) mod dual;`. **`rules/flat_view.rs`**:
  `rate_expression` is used now (its `dead_code` allowance removed).
  **`lib.rs`**: the API table lists the model versions the kinetics read.
- **`integrations/diffsol`**: builds the kinetics with
  `Kinetics::from_bytes`, so it solves network models (with their
  population limits, read from the network model as before), rules models
  and `crn-craft.model/1.0`. Species names in diagnostics are the network
  model's names for a network model, else the view's labels; initial counts
  and parameter values come from the kinetics. Two new CLI tests.
- **`validation/mutation_check.py`**: the `kinetics_parameter_derivative_zeroed`
  anchor follows the new code; two new mutants
  (`kinetics_combination_derivative_dropped`, `dual_power_base_factor_dropped`).
- **`validation/merge/`**: the M5e declaration and results in `RESULTS.md`;
  the gate's report `kinetics/references.json`.

## Decisions (small and reversible, taken here)

- **Old network-type entry points.** `Kinetics::compile(&contracts::Model)`
  stays until M8, as a thin wrapper: it converts the network model (with
  its population limits removed, since the converter refuses them and the
  kinetics never read them) and builds the kinetics from the flat view.
  The old *implementation* goes now: unlike M5a's FSP command, the kinetics
  have no output that must stay byte-identical (D2), and one implementation
  serves every model version. The public fields `Kinetics::model` and
  `Kinetics::network` are removed now. New: `Kinetics::new(&FlatView)`,
  `Kinetics::from_bytes(&[u8])` (any model version; a network model goes
  through `compile`), `species()`, `initial()`, `parameter_values()`.
- **Whole-rule propensities and their gate.** The stochastic engine
  applies a whole-rule propensity only while every taken count is present.
  That gate is a step with no smooth extension, so the ODE drops it, but
  only where dropping it is safe: a whole-rule propensity that takes
  molecules is accepted when its expression is 0 whenever a taken species
  is 0, judged from the structure (the species' observable as a factor of a
  product, the numerator of a quotient, …), as the converted
  Michaelis–Menten form `vmax * (S / (km + S))` is. Otherwise `CRN-E080` at
  `/rules/i/propensity`. This keeps the 0.6.0 acceptance on converted
  models (Hill laws with reactants and Michaelis–Menten laws consuming more
  than their substrate stay refused; W4 stays refused) and keeps the
  diffsol domain argument valid: every accepted channel that consumes a
  species is 0 when that species is 0.
- **Branch rules.** Where an expression has no derivative, the selected
  branch's derivative is used: `if` and `&&`/`||` follow the evaluation;
  `min(x, y)` takes x when x ≤ y and `max(x, y)` when x ≥ y (ties: the first
  argument); `abs(x)` takes +x at 0 (as `if(x < 0, -x, x)`); `floor`,
  `ceil`, comparisons and logic have derivative 0. A chain-rule factor of 0
  gives 0 even where the other factor is infinite (a literal exponent at a
  base of 0, where ln 0 = −∞).
- **Singular points and infinite intermediates.** At a base of 0, x^y has
  the one-sided derivatives (0 for y > 1, 1 for y = 1, ∂/∂y = 0 for y > 0);
  for 0 < y < 1, and for `sqrt` and `log` at 0, the derivative is infinite
  and the evaluation fails with "singular at zero". The converted Hill
  activation `basal + vmax / (1 + (k / R)^n)` has k / R = ∞ at R = 0; an
  infinite intermediate carries the derivatives of its reciprocal, so the
  result gets the limits of the quotient rule: the 0.6.0 one-sided
  derivatives (0 for n > 1, vmax / k for n = 1, singular for n < 1). Without
  this the derivative would be NaN at every regulator count of 0.
- **Caps.** As 0.6.0: at most 64 species, 32 parameters and 256 channels
  (reactions), and at most eight molecules of one species in a channel with
  a `rate` (the macroscopic form is a running product of y / j). Whole-rule
  propensities have no order cap (no combination term).
- **Inputs.** `y` must be finite and non-negative, `p` finite; negative
  parameters are allowed now (rules models may have them), and a negative
  or non-finite propensity is `CRN-E080` naming the rule.

## Tests written first

The implementation was drafted before the tests; the tests were then run
first against stubs with the new API (`Kinetics::new` and `evaluate`
returning an error; `Program::evaluate` returning `Expr::eval`'s value with
zero gradients and `vanishes` returning false). Log:
`/home/claude/m07/m5e/first-fail.log`.

- `crates/crn-craft/tests/kinetics.rs`: 7 of 7 failed against the stubs.
  The 0.6.0 test is kept with its body unchanged (`Kinetics::compile` of a
  parsed network model, central differences); new: network, parsed, converted and
  view-built kinetics are bit for bit equal for five examples, with the
  network's species, initial counts and parameters; finite-difference
  Jacobians of a `model/1.0` with `^`, `/`, `exp` and `max` (both branches
  of `max`), a count-2 pattern and a consuming whole-rule propensity, and
  its right-hand side by hand; higher-order (`Multi`) channels against
  closed forms (3A, 2A + B with an observable rate, A + B + C), zero counts,
  the converted M4b models, the order cap; whole-rule propensities
  (Michaelis–Menten bit for bit, Hill activation at R = 0 for n = 1, 2,
  2.5 and the singular n = 0.5, refusals of Hill with reactants and of
  Michaelis–Menten consuming another species, the dropped gate of a count-2
  propensity); first-order channels equal Σ changes × `FlatView::propensity`
  at integer states for four examples and an attributed model with pattern
  factors and weighted observables; invalid inputs, negative and non-finite
  propensities, observable overflow, the size cap, `CRN-R017`, `CRN-E001`,
  and ignored population limits.
- `rules/dual.rs` unit tests: 5 of 6 failed against the stubs (values equal
  to `Expr::eval` bit for bit on 2,048 random expressions passed, as the
  stub keeps the value). The others: central differences of 13 smooth
  expressions, the branch rules, powers at a base of 0, infinite
  intermediates (Hill at R = 0, overflow of (x / k)^n, exp of −∞, exp(log x)
  at 0 undefined), and `vanishes`.
- `integrations/diffsol/tests/cli.rs`: the 2 new tests failed against the
  previous `main.rs` (it parsed only network models).
- `crates/crn-craft/tests/kinetics_references.rs`: the gate, written after
  the declaration was committed and run once.

## Measured

| Check | Result |
| --- | --- |
| Reference gate (declared in `RESULTS.md` before the run) | 6,252 values and 12 acceptance checks, **0 failures**; 6,192 bit for bit equal; largest relative difference 2.6e-16 (Hill), 1.7e-16 (Michaelis–Menten); all 5,422 zero references exact; W4–W7 refused with `CRN-E080` as in 0.6.0 |
| `integrations/diffsol` (`--locked`, temporary target) | builds; 10 of 10 tests (8 CLI, 2 new; 2 kinetics) |
| `check_ode.py` | 50 of 50; identical under `-O`; 31 rows equal to the last stored report, 19 differ in trailing digits |
| Mutants (`kinetics_parameter_derivative_zeroed`, `kinetics_combination_derivative_dropped`, `dual_power_base_factor_dropped`), applied one at a time, `cargo test -p crn-craft --lib --test kinetics` | 3 of 3 caught (6, 5 and 6 failing tests) |

Evaluation time of `Kinetics::evaluate` at the probe state (µs per call,
median of five batches of 20,000, two alternating runs each, under the CPU
lock; the same program built against the 0.6.0 crate and this one; no
target declared, information only):

| Model | 0.6.0 | M5e |
| --- | ---: | ---: |
| immigration-death | 0.32 | 0.26 |
| michaelis-menten | 0.27 | 0.38 |
| enzyme-mechanism | 0.72–0.77 | 0.52 |
| repressilator | 2.02–2.08 | 1.93 |
| vilar-noise | 3.37–3.66 | 1.92–1.95 |
| W2 (33 species, 66 reactions) | 11.2 | 5.0–5.4 |

## Checks

On the M5e commits, release build, under the CPU lock; Python: `venv-m5e`
with the wheel built from the root.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass (also in `integrations/diffsol`) |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass after rewriting two negated float comparisons in the new gate test; `integrations/diffsol` clippy with `-D warnings` also passes |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 127 (6 new), `kinetics` 7 (6 new), `kinetics_references` 1 (new), all others unchanged; I1 (`stored_requests` 2, `static_direct_defaults` 9), I3 (`generated_loop` 26, `generated_kernel` 3), doctests 36 + 1 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (the known output-filename collision warning only) |
| `ruff format --check --exclude '*.md' .`, `ruff check .` (0.11.13) | pass |
| `cargo build --release --workspace --locked`; wheel from the root | pass |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` (`CRN_CRAFT_SOURCE`, `CARGO_TARGET_DIR`) | 22 blocks pass; identical under `-O` |
| `check_ode.py` | see above |

Not run: the full `mutation_check.py` (it rebuilds a copy of the tree; the
three kinetics and dual mutants were run as above, and only their anchors
changed), and the Python scripts of paths this task does not touch (no
schema, contract, Python or CLI code changed; the kinetics have no Python
or CLI entry point).

**CI.** CI runs neither `check_ode.py` nor the diffsol build and tests
(`.github/workflows/ci.yml` has no step for them; `docs/technical-report.md`
says the integrations are checked locally because they need newer
toolchains). The new gate test runs in CI as part of `cargo test`.

## Limitations

- Hill laws are evaluated as the converted fraction, not the 0.6.0
  logistic form; values agree to 2.6e-16 relative in the references.
- `vanishes` is structural: an expression that is 0 at a zero count only
  by its values (for example `k * (S > 0)`) is refused.
- Presence-based observables (`distinct`, `max`, `min`) are piecewise
  constant in the counts; their derivative is 0, including at the jump.
- A rule's rate is evaluated also where its channel's combination term is 0
  (its derivative is needed there); a rate that is undefined at a zero
  count fails the evaluation, while the stochastic engine skips it.
- Models with ATP, continuous variables or products that create attribute
  tuples have no flat view and no kinetics (`CRN-R017`).
- The kinetics of a network model go through the converter and bring its
  refusals: a network without reactions or with an initial count above
  2^53 is refused (`CRN-R016`), and reserved names are renamed (`exp` →
  `exp_2`), so `Kinetics::parameters` holds the converted names.

## Deviations from the plan

- The gate is a Rust integration test, not a Python script: the kinetics
  have no CLI or Python entry point.
- The plan does not say how whole-rule propensities and their gate extend
  to real counts, nor the size and order caps; the decisions above fill
  that in. Infinite intermediates (reciprocal derivatives) are not in the
  plan; the converted Hill activation needs them at R = 0.
- `diffsol` now reads rules models and `model/1.0` (the plan only asks that
  it builds and passes against the new interface).
- Two mutants were added to `mutation_check.py`.

## Migration notes

- `crn_craft::kinetics::Kinetics` reads any closed model:
  `Kinetics::from_bytes(bytes)` (network `model/0.1`, rules models,
  `model/1.0`) or `Kinetics::new(&FlatView)`. `Kinetics::compile(&Model)`
  still takes a parsed network model (until the network types are removed);
  `evaluate` and `Derivatives` are unchanged. Removed: the public fields
  `Kinetics::model` and `Kinetics::network`; use `species()`, `initial()`,
  `parameters` and `parameter_values()`.
- Behavior: channels replace reactions in the size cap; parameters may be
  negative; a negative propensity is `CRN-E080`; Hill laws with reactants
  and Michaelis–Menten laws consuming other species are refused at
  `/rules/i/propensity` (0.6.0: `/ode`) with a new message; a Hill law
  with k ≤ 0 is no longer refused by name (its expression's value decides);
  the evaluation messages name the rule.
- `crn-diffsol MODEL REQUEST` accepts rules models and `crn-craft.model/1.0`
  (no population limits there); output unchanged in shape.

## For later tasks

- M8: `Kinetics::compile` and the `model/0.1` branch of
  `Kinetics::from_bytes` use `crate::parse_model` and `contracts::Model`
  (to drop population limits before converting). When `contracts` moves
  into `convert.rs`, give the converter a way to ignore population limits
  for the kinetics, and remove `compile`. `integrations/diffsol` reads
  population limits with `crn_craft::parse_model`; it needs its own reader
  then (M9).
- `validation/merge/kinetics-reference` still compiles against this API but
  would print this engine's values; the M0c references were made with the
  kept 0.6.0 binary `/home/claude/m07/bin060/kinetics-reference`.
- M10: no user document describes the kinetics' model versions yet
  (`docs/architecture.md` and `docs/tutorial.md` mention the diffsol
  integration only).
- Run evidence (logs of every check, the failing-first runs, the mutant
  and timing runs, the ODE reports) is in `/home/claude/m07/m5e/`.
