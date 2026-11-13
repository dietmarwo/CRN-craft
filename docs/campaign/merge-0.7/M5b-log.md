# M5b log: network analysis on the flat view

Task M5b of [the plan](../../history/plan-0.7.0-engine-merge.md): `analysis.rs`,
`conservation.rs`, `cone.rs` and `linkage.rs` take a `FlatView`;
`crn-craft analyze MODEL` reads `crn-craft.model/1.0`, the rules models and
(through the converter) `crn-craft.model/0.1`; `analysis_v1.rs` and the
`analyze-v1` command are removed.

## What changed

- **`conservation.rs`, `cone.rs`, `linkage.rs`**: `analyze_conservation`,
  `analyze_cone`, `analyze_cone_with_limits` and `analyze_linkage` take
  `&FlatView` instead of `&CompiledNetwork`. They read the species count,
  the channels' net changes (`changes`), the initial counts and, for the
  complex graph, the channels' reactants and products. The arithmetic, caps,
  statuses and diagnostics are unchanged (the mutation anchors of
  `mutation_check.py` still match once each).
- **`analysis.rs`** rewritten around the view. `analyze_bytes(model,
  state)` builds the view with `FlatView::of_any_model`, reads an optional
  `crn-craft.state/0.1` document by species label, evaluates every channel's
  propensity with `FlatView::propensities` (the Static Direct value), and
  reports `crn-craft.analysis/1.0`. The report has the fields of
  `analysis/0.2`: species are the view's labels, reactions its channels
  (named by their rule), dependencies the channels Static Direct updates.
  `Term` and `ReactionStructure` are now defined in `analysis.rs` (the
  network engine keeps its own copies until M8), with the same names, so the
  generated schema has the same shape as 0.2; only the version string and
  descriptions differ.
- **Warnings** come from the view, in 0.6.0's order: per reaction
  `CRN-W003` (no net change) and `CRN-W002` (a rate that is the constant 0),
  then per species `CRN-W001` (net produced, never net consumed), dropped
  for species with a proven bound as before; then the conservation and cone
  warnings.
- **`rules/flat_view.rs`**: `FlatView::with_every_type(model)` lowers the
  model with a zero-count initial entry for every attribute-free type that
  has none, so a type that never holds molecules keeps its species and the
  rules over it their channels (M5a's view drops them). Analysis and the
  GillesPy2 export (M5f) use it; FSP keeps `FlatView::new`, so M5a's FSP
  requests are unchanged. `FlatView::of_any_model(bytes)` builds that view
  from any model version, converting a `model/0.1` document first.
- **Contracts**: `schemas/analysis-1.0.schema.json` (new, `crn-craft
  describe analysis-1.0`); `schemas/state.schema.json` regenerated
  (descriptions only: labels, counts at most 2^53). The `describe` kinds
  `analysis` and `analysis-v1` are gone; their files stay until M8 removes
  them. Core and Python discovery move to `crn-craft.describe/0.11`
  (commands without `analyze-v1`, schema map key `analysis-1.0` instead of
  `analysis` and `analysis-v1`), the native CLI to `describe/0.15`; the
  rules entry's `contracts_1_0` gains `analysis`.
- **CLI and Python**: `crn-craft analyze` prints `analysis/1.0`;
  `analyze-v1` is an unknown command (`CRN-E031`); the help text says which
  commands read which models. `crn_craft.analyze(model)` accepts every model
  version (bytes, text or a dictionary). `python -m crn_craft analyze-v1`
  is removed; the native `analyze_state_json(model, state=None)` lost its
  `legacy` argument.
- **Validation**: `validation/merge/check_analysis_references.py` (the
  gate) and its report; `check_network.py`, `check_structure.py`,
  `check_contracts.py` and `check_installed.py` migrated in place (the
  `analyze` command itself changed, so there is no old path to keep testing);
  `mutation_check.py`: the anchor of `analysis_drops_fallback_bounds`
  follows the new indentation and variable name.

## Tests written first

Run against a stub `analyze_bytes` that returned a failed report, with
conservation, cone and linkage already ported (logs in
`/home/claude/m07/m5bf/first-fail-*.log`):

- `crates/crn-craft/tests/analysis.rs` (new, 10 tests): 9 of 10 failed. It
  checks for the five examples, W1–W7 and both higher-order models that the
  report of the network file equals the report of the converted model apart
  from the hash, and that species, initial counts, reactions (names,
  reactants, products, net changes, reads) and dependencies equal the
  network engine's, with propensities within the M2/M4b tolerances;
  attributed models with labels, states by label, a rule with several
  channels; higher-order and whole-rule coefficients; types that never hold
  molecules; warnings and their paths; refusals (`CRN-R016`, `CRN-E012`,
  `CRN-E003`, `CRN-E001`, `CRN-R017` at `/rules` and `/atp`, `CRN-R012`);
  state errors (`CRN-E041` by label, `CRN-E002`, `CRN-E003`); `CRN-R030` on
  overflow and 0 on underflow; files, the schema snapshot and discovery. The
  test that passed (`the_analysis_view_of_a_converted_model_equals_the_converter_output`)
  checks `FlatView::of_any_model`, written before the tests.
- `crates/crn-craft/tests/structure.rs` rewritten on `model/1.0` models and
  `FlatView::with_every_type`, with the conservation tests of
  `tests/network.rs` moved in: 6 of 27 failed (those that call
  `analyze_bytes`); the 21 ported structure tests passed against the ported
  modules, as ports should.
- `crates/crn-craft-cli/tests/cli.rs`: 3 failed (`analysis/1.0`, state
  errors through the new path, the new
  `analyze_reads_network_rules_and_1_0_models_and_analyze_v1_is_gone`).
- `tests/contracts_v1.rs`: the network-entry-point test now expects analysis
  to accept 1.0 and rules models (failed first).
- The four GillesPy2 export tests of the old `tests/structure.rs` moved
  unchanged to `tests/export.rs` (rewritten in M5f).

Porting notes: a rules model needs a rule, so the 0.6.0 fixtures "without
reactions" get one rule without reactants or products (a zero column, which
leaves conservation and the cone unchanged; the linkage test of an empty
network keeps its other assertions); a rule has at most 64 products, so the
overflow chains put the coefficient 2147483647 on the reactant side
(`count`); counts are at most 2^53 (0.6.0 used `i64::MAX` in one overflow
fixture).

## Measured

Declared in `validation/merge/RESULTS.md` before the runs; rows there.

| Check | Result |
| --- | --- |
| Reference gate (twelve models, network file and converted model) | 194,388 leaves, **0 failures**, no propensity difference; `-O` identical |
| `check_network.py` | 219 of 219; same fixture stream as 0.6.0; 1,030 rates, largest relative error 4.36e-15; 21 refusals without reactions, 1 overflow (`CRN-R030`), 2 underflows to 0, 1 Hill floor case (see below); `-O` identical |
| `check_structure.py` | 309 of 309; 37 refusals without reactions; 211 rays, 14 additional cone bounds; `-O` identical |
| `check_contracts.py` | 101 of 101 (the `analysis-v1` schema row gone, 9 rules-model analysis rows new); `-O` identical |
| Mutants of `mutation_check.py` in `analysis.rs`, `conservation.rs`, `cone.rs`, `linkage.rs`, `export.rs` (applied one at a time, `cargo test -p crn-craft --test structure --test analysis --test export`) | 11 of 11 caught; `cone_initial_ray_cap_ignored` survived the first run (the ported test gave the empty network a zero column), caught after a test of a view without channels was added |

## Checks

On the M5b and M5f commits together (one tree), release build, under the
CPU lock; Python: `venv-m5bf` with the wheel built from the root and
installed. Each Python script ran normally and under `python -O` with
identical reports.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 121, `analysis` 10 (new), `structure` 27, `export` 8, `contracts` 35, `contracts_v1` 7, `network` 25, `flat_view` 6, `fsp` 18, `ssa_alternatives` 15, CLI 15, I1 (`stored_requests` 2, `static_direct_defaults` 9), I3 (`generated_loop` 26, `generated_kernel` 3), all others, doctests 36 + 1 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (the known output-filename collision warning only) |
| `cargo build --release --workspace --locked`; wheel from the root | pass |
| Schemas against `crn-craft describe KIND` | 24 byte-identical (`analysis-1.0` new; `state`, `export` regenerated, descriptions only) |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` | 22 blocks pass |
| `check_network.py`, `check_structure.py`, `check_contracts.py` | see above |
| `check_export.py`, `check_importers.py`, `check_installed.py` (`python -I`, outside the tree) | see the M5f log; 9, 49 and 17 checks |
| `check_python.py` | 12 checks pass |
| `check_rules_python.py` | 26 checks pass |
| `check_rules.py` | 47 of 47 |
| `check_tutorial.py` | 17 blocks pass |
| `check_dependencies.py` | pass |

Not run: the full `mutation_check.py` (it rebuilds a copy of the tree; the
11 mutants of the changed files were run as above) and the scripts of paths
these tasks do not touch (`check_ssa.py`, `check_cme.py`, `check_fsp*.py`,
`check_tau.py`, `check_mtdna.py`, `check_compiled_loop.py`,
`check_builder.py`, `check_hybrid.py`, `check_alternative_*.py`).

## Deviations from the plan

- `FlatView::with_every_type` and `FlatView::of_any_model` are new; without
  them a converted network with a species that starts at 0 and is never
  produced would lose that species and its reactions (M5a's note).
- The analysis keeps the optional `crn-craft.state/0.1` input (by label).
  The plan's command line lists `analyze MODEL` and M8b deletes the `state`
  schema; `check_export.py` and `compare_julia.py` need rates at a given
  state, so M8 must keep that schema or replace it.
- Warning paths point into the report's lists: `CRN-W001` at
  `/species/i` and `CRN-W003` at `/reactions/i` as before, `CRN-W002` now at
  `/reactions/i` (0.6.0: `/reactions/i/rate/coefficient`) with the message
  "reaction rate is zero".
- Discovery moved to `describe/0.11` (core and Python) and `0.15` (native)
  because a command and two schema kinds were removed; M7 rewrites
  discovery anyway.
- `schemas/analysis.schema.json` is no longer generated; M8b's list of
  schemas to delete names `analysis-v1` but not `analysis`, so M8 should
  delete both.
- `check_network.py`: one Hill case (`hill_activation_1_1e+300`) is accepted
  within M2's declared floor, which the M5b declaration did not list; see
  RESULTS.md.

## Limitations

- Analysis of a network model goes through the converter: models without
  reactions and models with population limits are refused (`CRN-R016`);
  reserved names are renamed as the converter does (`exp` → `exp_2`), and
  the report and state documents use the converted names.
- A propensity that overflows fails with `CRN-R030` (0.6.0: `CRN-E042`); a
  positive propensity that underflows is reported as 0 (0.6.0 refused it).
- Hill laws with an induced term below vmax · 2^-1020 give the basal rate
  (M2).
- Models with ATP, continuous variables or products that create attribute
  tuples during a run have no analysis (`CRN-R017`).

## Migration notes

- `crn-craft analyze MODEL [STATE]` prints `crn-craft.analysis/1.0` and
  reads `model/0.1`, `rules-model/0.1`–`0.4` and `model/1.0`; the state
  names species by label (`A`, or `M{wt=1}` for an attributed type).
- Removed: `crn-craft analyze-v1`, `python -m crn_craft analyze-v1`,
  `crn_craft::analysis_v1`, `crn_craft::analyze_legacy_files`, the
  `describe` kinds `analysis` (0.2) and `analysis-v1` (0.1). New kind:
  `analysis-1.0`.
- Rust: `analyze_conservation`, `analyze_cone`, `analyze_cone_with_limits`
  and `analyze_linkage` take `&crn_craft::rules::FlatView` (build it with
  `FlatView::of_any_model(bytes)` or `FlatView::with_every_type(&model)`);
  `crn_craft::analysis::{Term, ReactionStructure}` replace the network
  types in the report.
- Python `crn_craft.analyze(model)` returns `analysis/1.0` for any model
  version.
- Discovery: `crn-craft.describe/0.11` (core, Python), `0.15` (native CLI).

## For later tasks

- M5e: `FlatView::with_every_type` gives the network's species for
  converted models, which kinetics needs for the same reason as analysis.
- M6a: `Model.analyze()` can call `crn_craft.analyze(self.to_dict())`.
- M8: delete `schemas/analysis.schema.json` and
  `schemas/analysis-v1.schema.json`; keep or replace `state`.
