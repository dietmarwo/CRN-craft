# M5a log: the flat view and finite state projection

Task M5a of [the plan](../../history/plan-0.7.0-engine-merge.md): the common basis of
M5 (`FlatView`, a read-only view of a closed model lowered to a flat
network) and the `fsp` method of `crn-craft.run/1.0` with its result
`crn-craft.fsp-result/1.0`.

## What changed

- **`crates/crn-craft/src/rules/flat_view.rs`** (new, `pub mod flat_view`,
  `crn_craft::rules::FlatView`). Built from model bytes
  (`FlatView::from_bytes`) or a parsed `rules::Model` (`FlatView::new`):
  it compiles the model and lowers it with the same `flat::Network::compile`
  the static methods run. It reads the lowered network only; `flat.rs` is
  unchanged. Species are in *view order*: by type in declaration order,
  then in lowering order (initial tuples, then the tuples products create).
  Labels are the type name (`A`) or the type with every attribute
  (`mtDNA{wt=1,mut=0}`; integers below 2^53 without a fraction, other values
  in Rust's shortest round-trip form). For a converted network model the
  view order is the network's species order. `propensity` loads the state,
  refreshes every observable with `Network::refresh` and calls
  `RateCache::channel`, the function Static Direct calls, so the value is
  the engine's, including `Multi` (higher-order) channels and whole-rule
  `propensity` expressions. The full API is listed below.
- **New diagnostic `CRN-R017`**: a model that has no flat view. At `/rules`
  for a lowering failure (products that create attribute tuples during a
  run, random product attributes, a whole-rule propensity over a type with
  several tuples, the static tuple and channel limits), at `/atp` for an ATP
  model, at `/continuous` for continuous variables. `docs/diagnostics.md`
  lists it and the new `CRN-R020` and `CRN-E070` cases.
- **`crates/crn-craft/src/fsp.rs`** rewritten around one numerical core,
  `project`, generic over a private `Generator` trait (initial state, net
  changes, the rates of a state, an admission check). The enumeration,
  uniformization (split so λ·dt ≤ 16), error accounting, limits and
  messages are the old code, moved; the only generalization is that the
  step count runs over the intervals from time 0 through the report times,
  which equals the old count for a grid that starts at 0. Two generators:
  - the network `CompiledNetwork` for the old `fsp::solve` (`run/0.1` →
    `fsp-result/0.1`, the `crn-craft fsp` command, kept until M8);
  - the `FlatView` for `fsp::v1::solve` (`run/1.0` with method `fsp` →
    `fsp-result/1.0`). A state where the model's `stop` condition holds is
    absorbing (no exits), as a stopped replicate holds its state; the
    network models had no stop.
- **`crn_craft::fsp::v1`** (new module in `fsp.rs`): `FspRequest`
  (`schema`, `model`, `time` with `end` and `report_times`, `method`
  `{"kind": "fsp", "bounds": {LABEL: max}, "max_states": N, "tolerance":
  t}`, `budget` with `max_output_bytes` and optional `wall_time_s`; no
  ensemble, outputs, trackers or event log), `FspResult`, `FspSnapshot`,
  `FspVersion` (the fields and meanings of `fsp-result/0.1`; `species` are
  the flat labels, `states` unsigned, one snapshot per report time),
  `solve`, `solve_files`, `is_fsp_request`, `is_fsp_request_file`. Request
  errors are `CRN-R020` (version at `/schema`, binding, time grid with the
  simulation request's messages, `/method/max_states`, `/method/tolerance`,
  budget, and at `/method/bounds/LABEL` an unknown label, a missing bound, a
  bound below the initial count, or a type that never holds molecules);
  numerical failures keep `CRN-E070` at `/fsp`; rate and stop failures are
  `CRN-R030`. The output limit is `max_output_bytes`, capped at 1 GiB. The
  run is polled through `rules::Control` (flag every poll, cancel file and
  wall clock at the first and every 1,024th poll).
- **Contracts and schemas**: `schemas/fsp-result-1.0.schema.json` (new,
  `crn-craft describe fsp-result-1.0`, from `crn_craft::schema`) and
  `schemas/run-1.0.schema.json` regenerated: its top level is now
  `oneOf` of `definitions/RunRequest` (the simulation request, unchanged)
  and `definitions/FspRequest`, with a description of both. The other 23
  schema files are byte-identical to `describe` output, the `rules-*` and
  `fsp-result` (0.1) schemas included. `rules::describe()` gains
  `contracts_1_0.fsp` (`method`, `result`, `schema`, `python`).
- **CLI**: `crn-craft rules MODEL REQUEST [--cancel-file PATH]` reads the
  request once leniently; a request whose `method.kind` is `fsp` goes to
  `fsp::v1::solve_files` and prints the compact `fsp-result/1.0` (exit 0),
  or `{"success": false, "diagnostics": [...]}` (exit 1). The help text
  lists `fsp-result-1.0` and the dispatch.
- **`rules::prepare`** (and so `rules.ensemble(..., method="fsp")`) refuses
  an FSP request with `CRN-R020` at `/method` naming the FSP entry points,
  instead of the strict parser's unknown-variant `CRN-E003`. The check runs
  only when the strict parse fails.
- **Python**: native `rules_fsp_json(model, request, cancel_file=None)`
  (Ctrl-C cancels); `crn_craft.rules.fsp(model, end, report_times, bounds,
  max_states=None, tolerance=None)` returns the parsed `fsp-result/1.0`
  document. `model` is bytes, text, a dictionary or a `Model` of any rules
  version; defaults `FSP_MAX_STATES = 100_000`, `FSP_TOLERANCE = 1e-12`,
  output limit 1 GiB. `Model.fsp` and DataFrames are M6a.
- **`rules/mod.rs`** edits, kept small for the merge with M4c: `pub mod
  flat_view; pub use flat_view::FlatView;` between `mod flat;` and `mod
  generated;`; `fn read` made `pub(crate)` and the model-path check of
  `solve_files` moved into `pub(crate) fn references_model`; the FSP
  refusal in `prepare`'s request parse; the `run-1.0` hook at the end of
  `schema()` and its doc paragraph; one `contracts_1_0.fsp` line in
  `describe()`.
- **Validation**: `validation/merge/check_fsp_references.py` (the gate,
  report `validation/merge/fsp/references.json`),
  `validation/check_fsp_rules.py` (every case of `check_fsp.py` on the new
  command; CI runs it next to `check_fsp.py` in the `reference-export`
  job), an FSP case and its refusal in `check_contracts.py`, `rules.fsp`
  checks in `check_rules_python.py`. `mutation_check.py`: the anchor of
  `fsp_reflects_sink` follows the code's new indentation.

## Inputs of the network features (read for M5b, M5e, M5f)

| File | Reads from `CompiledNetwork` / the network model | From `FlatView` |
| --- | --- | --- |
| `analysis.rs` | `species()`, `initial_counts()`, `reactions()` (name, `reactants`, `products`, `net_changes` as `Term`s, `reads`), `dependencies()`, `propensities(state)`, `named_state`, `warnings()` (`CRN-W001` net producer without net consumer, `CRN-W002` zero mass-action coefficient, `CRN-W003` no net effect) | `species`, `initial`, `channels`, `rule`/`rule_name`, `reactants`, `products`, `changes`, `reads`, `dependencies`, `propensities`, `index`; W001/W003 from `changes`, W002 from `rate_constant` and `whole` |
| `conservation.rs` | species count, `reactions()[r].net_changes`, `initial_counts()` | `changes`, `initial` |
| `cone.rs` | species count, reaction count, `net_changes`, `initial_counts()` | `changes`, `initial` |
| `linkage.rs` | `reactants` and `products` (complexes as sorted (species, coefficient)) | `reactants`, `products` (sorted by species, merged counts) |
| `kinetics.rs` | the model's rate laws and `Scalar` parameters, reactant orders, `net_changes`, sorted parameter names | `reactants` (orders for y^k/k!), `whole`, `factor`, `rate_expression` (crate-internal `Expr`, slots follow `parameter_names` and `observable_names`), `parameter_names`/`parameter_values`, `observable(o)` (reduction and (species, weight) terms), `rate_observables`, `rate_constant`, `changes` |
| `export.rs` (GillesPy2) | species and initial counts, reactions with reactants/products, mass-action coefficients | `species`, `initial`, `reactants`, `products`, `rate_constant` with `factor` for mass action, `rate_source` and `observable` for custom propensities |

## FlatView API summary (for M5b, M5e, M5f)

```rust
use crn_craft::rules::FlatView;           // also crn_craft::rules::flat_view::NO_FLAT_VIEW ("CRN-R017")
FlatView::from_bytes(&[u8]) -> Result<FlatView, Vec<Diagnostic>>
FlatView::new(&rules::Model) -> Result<FlatView, Vec<Diagnostic>>
fn model(&self) -> &Model
fn species(&self) -> &[String]            // labels, view order
fn index(&self, label: &str) -> Option<usize>
fn species_type(&self, s) -> usize        // index into model.types
fn attributes(&self, s) -> &[f64]
fn initial(&self) -> &[u64]
fn channels(&self) -> usize
fn rule(&self, c) -> usize; fn rule_name(&self, rule) -> &str; fn rules(&self) -> usize
fn changes(&self, c) -> &[(usize, i64)]   // net change, sorted, zeros left out
fn reactants(&self, c) -> &[(usize, u64)] // (species, molecules taken), sorted, merged
fn products(&self, c) -> &[(usize, u64)]
fn whole(&self, c) -> bool                // whole-rule propensity vs rate × Π C(n, k)
fn factor(&self, c) -> f64                // product of pattern factors
fn reads(&self, c) -> &[usize]            // species the propensity reads (reactant gates included)
fn dependencies(&self, c) -> &[usize]     // channels Static Direct updates after c fires
fn rate_source(&self, rule) -> &str       // `propensity` text, else `rate` text
fn rate_constant(&self, rule) -> Option<f64>  // value when parameters only
fn rate_observables(&self, rule) -> Vec<usize>
pub(crate) fn rate_expression(&self, rule) -> &Expr
fn parameter_names(&self) -> Vec<&str>; fn parameter_values(&self) -> &[f64]   // sorted
fn observable_names(&self) -> Vec<&str>   // sorted, the order expressions index
fn observable(&self, o) -> (Reduce, Vec<(usize, f64)>)
fn has_stop(&self) -> bool; fn stops(&self, counts) -> Result<bool, Diagnostic>
fn propensity(&self, c, counts: &[u64]) -> Result<f64, Diagnostic>
fn propensities(&self, counts: &[u64]) -> Result<Vec<f64>, Diagnostic>
pub(crate) fn evaluator(&self) -> Evaluator  // load(counts), channel(c), stops(): scratch reuse
```

Counts are in view order. Propensity errors are `CRN-R030` (at `/counts`
for a wrong length or a count above 2^53, `/observables` for an overflow,
`/rules/i/rate` or `/rules/i/propensity` for a negative or non-finite value).

## Tests written first

The implementation was drafted before the tests were written; the tests
were then run first against stubs (`FlatView::new` and `fsp::v1::solve`
returning an error, schemas not regenerated, no CLI dispatch) and failed:
`tests/flat_view.rs` 6 of 6, the 12 new `v1_*` tests of `tests/fsp.rs` 11 of
12 (the `prepare` refusal was already implemented), the schema snapshot test
of `contracts_v1.rs`, and 2 CLI tests (the new FSP test and the run-1.0
snapshot). The 6 old `fsp.rs` tests passed throughout.

- `tests/flat_view.rs` (new, 6 tests):
  `converted_models_have_the_network_structure_and_propensities` (every
  example, W1–W7 and both higher-order models: labels, initial counts,
  channel per reaction, changes, reactants, products, reads, dependencies,
  factor and `whole` equal the network engine's; propensities at the
  initial, a probe, the empty and 12 pseudo-random states within the
  tolerances M2 and M4b measured: bit-identical for mass action of order 0,
  1 and ≥ 3 and Michaelis–Menten, 4 ε relative over pairs, 1e-14 · max(1, n)
  relative for Hill with the vmax · 2^-1020 floor); labels and API of an
  attributed model; higher-order and whole-rule propensities and `stops`;
  the `CRN-R017` refusals; ATP models refused; coded propensity errors.
- `tests/fsp.rs`: the six network tests kept unchanged (they guard the old
  command until M8) and ported to converted models through `run/1.0`
  (`v1_binomial_decay_and_large_lambda_t`, `v1_leakage_…`,
  `v1_capacity_cancellation_bounds_and_null_reaction`, `v1_omitted_weight_…`,
  `v1_the_tail_bound_…`, `v1_many_snapshots_…`), plus
  `v1_agrees_with_the_network_fsp` (immigration-death, michaelis-menten,
  enzyme-mechanism and both higher-order models: same states and edges,
  probabilities within 1e-12 of the old command), `v1_stop_states_are_absorbing`
  (birth with stop `A >= 3` against the Poisson law),
  `v1_attributed_models_use_flat_labels` (Binomial law over `M{wt=0}`,
  `M{wt=1}`), `v1_requests_are_checked_with_coded_diagnostics`,
  `v1_fsp_requests_are_recognized_and_refused_by_prepare`,
  `v1_solve_files_resolves_the_model_and_honors_a_cancel_file`.
- `crates/crn-craft-cli/tests/cli.rs`:
  `rules_solves_fsp_requests_and_describes_their_contracts` (output equals
  the library's result as compact sorted JSON, refusal exit 1, cancel file,
  `describe fsp-result-1.0` and `run-1.0` equal the files).
- `tests/contracts_v1.rs`: the snapshot test follows `oneOf` for `run-1.0`
  (both request shapes reference the filtered version definition) and
  checks `fsp-result-1.0`.

## Measured

| Check | Result |
| --- | --- |
| FlatView propensities against the network engine | 27,150 compared, 27,013 bit-identical, all within the declared tolerances |
| Reference gate `check_fsp_references.py` (declared in RESULTS.md) | 11,278 checks, **0 failures**; largest probability difference 1.4e-17 (enzyme-mechanism), the other six models bit-identical; λ and edge counts equal; `-O` identical |
| `check_fsp_rules.py` (all `check_fsp.py` cases, its tolerances, DSMTS `e79fd58`) | 16 of 16; largest distribution error 4.49e-14 (as the old script); `-O` identical |
| `check_fsp.py` on the network command | 16 of 16; report byte-identical to the 0.6.0 binary's (it holds every result's SHA-256) |
| FSP mutants of `mutation_check.py` (applied one at a time, `cargo test -p crn-craft --lib --test fsp`) | 3 of 3 caught, each by old and ported tests |

The declaration counted 38 sink comparisons; the references have 39 report
times, so N is 11,278, not 11,277 (RESULTS.md says so).

## Checks

On the M5a commits, release build, under the CPU lock. Python: venv-m5a with
the wheel built from the root (`maturin build --release --locked --out
dist`).

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass (one loop in a new test fixed after the first run) |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 121, `flat_view` 6 (new), `fsp` 18 (12 new), `contracts_v1` 7, CLI 15 (1 new), `convert` 10, `higher_order` 8, `quantiles` 13, I1 (`stored_requests` 2, `static_direct_defaults` 9), I3 (`generated_loop` 26, `generated_kernel` 3), all others, doctests 36 + 1 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (the known output-filename collision warning only) |
| Schemas against `crn-craft describe KIND` | 23 byte-identical, `run-1.0` regenerated, `fsp-result-1.0` new |
| `ruff format --check --exclude '*.md' .`, `ruff check .` (0.11.13) | pass |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` (`CRN_CRAFT_SOURCE`, `CARGO_TARGET_DIR`) | 22 blocks pass; `-O` identical |
| `check_contracts.py` | 93 of 93 (3 new); `-O` identical |
| `check_rules_python.py` | pass (new `rules.fsp` check); `-O` identical |
| `check_rules.py` | 47 of 47; `-O` identical |
| `check_mtdna.py` | 151 checks pass; `-O` identical |
| `check_installed.py` (`python -I`, outside the tree) | 17 checks pass; `-O` identical |
| `check_dependencies.py` | pass; `-O` identical |

Not run: the full `mutation_check.py` (it rebuilds a copy of the tree; only
`fsp.rs` anchors moved and those three mutants were run), and the scripts of
paths this task does not touch (`check_compiled_loop.py`: the generator is
unchanged and I3 holds by the fixtures; `check_python.py`,
`check_builder.py`, `check_hybrid.py`, `check_tutorial.py`).

## Limitations

- The FSP core is single-threaded, as before; its cost is the old core's
  plus one observable refresh per state.
- A type that has no initial instance and that no rule produces has no
  species in the view, and rules that need it have no channel. The network
  engine listed such species and reactions; none of the 14 network models
  of the repository has one. M5b's analysis of such a model would differ
  from the 0.6.0 report.
- Models with ATP or continuous variables have no flat view (`CRN-R017`).
- Bounds of attributed types are typed as labels (`M{wt=0}`), so the
  attribute formatting is part of the contract.
- `rate_expression` is crate-internal and unused until M5e
  (`#[allow(dead_code)]`).

## Deviations from the plan

- **Old command on the old types, one core** (the task's fallback, chosen
  over converting): `fsp::solve` keeps the network types behind the
  `Generator` trait. Converting the model first would change the old
  command's behavior before M8: the converter refuses population limits and
  models without reactions (`CRN-R016`), rate failures would become
  `CRN-R030` instead of `CRN-E042`, a population limit could no longer fail
  a reachable state with `CRN-E041`, and Hill propensities differ in the
  last bits. With the trait the old output is byte-identical to 0.6.0.
- **The FSP request has its own shape inside `run/1.0`** instead of the
  simulation request plus a method kind: an FSP request has no ensemble,
  outputs, trackers or event log, and its budget only `max_output_bytes`
  and `wall_time_s` (the old FSP contract also refused ensembles and event
  caps). Hence the `oneOf` top level of `run-1.0.schema.json` and the
  adapted M3 snapshot test. `Method` keeps `Copy` and no `fsp` variant, so
  M4c's `static_direct.selector` is untouched. Reversible until 0.7.0.
- `FlatView::initial` returns `&[u64]` and `reactants` `&[(usize, u64)]`
  (the plan: `Vec<u64>`, `Vec<(usize, u32)>`); counts of patterns over one
  type add, as in M4b. The view has more methods than the plan lists (above).
- Report times need not start at 0 (the rules grid); the distribution starts
  at 0 and the step count covers [0, first report time].
- Stop conditions make states absorbing (not in the plan; the network
  models had none).
- `fsp-result-1.0` is generated by `crn_craft::schema` next to
  `fsp-result`, not by `rules::schema`.
- `check_fsp.py` stays on the network command until M8 (CI runs it), and
  the migrated checks are the new `validation/check_fsp_rules.py`, also in
  CI.
- `prepare` refuses FSP requests with `CRN-R020` at `/method` (not in the
  plan).

## Migration notes

- New: `crn-craft.run/1.0` requests with method `fsp` (`bounds`,
  `max_states`, `tolerance`; `time`, `budget.max_output_bytes`, optional
  `budget.wall_time_s`), result `crn-craft.fsp-result/1.0`, run by
  `crn-craft rules MODEL REQUEST`, `crn_craft::fsp::v1::solve`/`solve_files`
  and `crn_craft.rules.fsp(model, end, report_times, bounds,
  max_states=None, tolerance=None)`. `crn-craft describe fsp-result-1.0`.
- `schemas/run-1.0.schema.json` now lists the simulation request and the
  FSP request under `oneOf`; the simulation request moved to
  `definitions/RunRequest`.
- New: `crn_craft::rules::FlatView` and diagnostic `CRN-R017`.
- `rules.ensemble(..., method="fsp")` and `rules::prepare` with an FSP
  request raise `CRN-R020` at `/method`.
- Unchanged: `crn-craft fsp`, `crn_craft.fsp`, `fsp::solve` and
  `fsp-result/0.1` (output byte-identical).

## For the coordinator and later tasks

- Merge with M4c: both change `schemas/run-1.0.schema.json`; regenerate it
  with `crn-craft describe run-1.0` after the merge. The `rules/mod.rs`
  edits are at the places listed above.
- M5b/M5e/M5f: build the view with `FlatView::from_bytes`; channel indices
  equal network reaction indices for converted models (one rule per
  reaction, one channel per rule), and view species follow the network
  order.
- Run evidence (reports of every check, the failing-first log) is in
  `/home/claude/m07/m5a/`.
