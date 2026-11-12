# M4e log: quantile statistics (FR1a)

Task M4e of [the plan](../../plan-0.7.0-engine-merge.md), the first part of
the feature request FR1: the median and other percentiles of every
observable at every report time, computed in the engine. Done on the branch
`merge/M4e` (worktree `/home/claude/m07/wt-m4e`) from `merge-0.7` after M0c,
M1, M2, M3 and M6b. The work was interrupted once by a server restart; the
section "Run history" says what that affected.

## What changed

### Contracts and engine (`crates/crn-craft/src/rules/`)

- `crn-craft.run/1.0`: `outputs.quantiles` (`Outputs::quantiles`,
  `Option<Vec<f64>>`, with `Outputs::quantiles()` returning the list or an
  empty slice). At most 16 (`MAX_QUANTILES`) strictly increasing
  probabilities in (0, 1); absent or empty means none. Refusals, all
  `CRN-R020`:
  - the field in any `rules-run` version, also when empty: at `/schema`
    ("outputs.quantiles requires crn-craft.run/1.0"), the convention of the
    other version gates;
  - more than 16 entries: at `/outputs/quantiles`;
  - an entry not in (0, 1) or not above its predecessor: at
    `/outputs/quantiles/<i>`, first offending index;
  - `null`: `CRN-E003`, as for every present optional field.
- `crn-craft.result/1.0`: `statistics.quantiles`, a list of
  `{"p": ..., "values": [[...], ...]}` (`Quantile`) in request order,
  `values` indexed like `statistics.mean` (report time, then observable,
  continuous columns and `ATP` included; an empty list without observables,
  as `mean` is).
- Withholding: `statistics` is reported when every replicate completed or
  stopped and the request asks for any statistic or any quantile; it is
  `null` otherwise and after an output-limit trim, as before. So the
  quantiles are withheld by exactly the rule of `mean`. A request with only
  quantiles gets `statistics` with `mean`, `sample_variance` and
  `mean_standard_error` `null`.
- Definition (`rules/quantile.rs`, documented on `Outputs::quantiles` and
  in the run-1.0 schema): NumPy's default `linear` quantile (type 7). With
  the n values sorted by `f64::total_cmp`, h = (n − 1)·p, j = ⌊h⌋,
  g = h − j: x[j] + g·(x[j+1] − x[j]) for g < 0.5,
  x[j+1] − (1 − g)·(x[j+1] − x[j]) for g ≥ 0.5, x[n−1] for h ≥ n − 1, NaN
  if any value is NaN. This is how `numpy.quantile` evaluates the plan's
  formula (`_lerp` in NumPy 2.3.5), so the values equal NumPy's bit for bit.
  The formula as written differs from NumPy in the last bit for about 1.2%
  of random integer samples (for example `numpy.quantile([16, 25], 0.925)`
  is 24.325, the formula as written 24.325000000000003); infinities follow
  NumPy too (`[1, inf]` at 0.5 is NaN).
- Buffer: `Reducer.quantiles`, replicates × report times × observables
  `f64`, allocated zeroed when quantiles are requested. `Reducer::fold`
  copies each completed replicate's values to offset `id × cells`; folds
  run in replicate order, so the buffer does not depend on the thread count.
  After the ensemble each cell's column is gathered into one sort vector,
  sorted, and evaluated for every probability.
- Memory: `admission` adds n × report times × observables × 8 bytes plus
  n × 8 for the sort column to the execution reservation when quantiles are
  requested. A request whose reservation fits without them but not with
  them is refused at `/outputs/quantiles` with the hint "Request fewer
  report times, fewer replicates or no quantiles." A request already over
  the limit without quantiles keeps its old `/budget` refusal.
- Every method that computes statistics gets quantiles from the same
  reducer: `tree`, `static_direct` (both waiting-time policies, sparse and
  per-event trackers), `next_reaction`, `mtdna_direct` (and its fallback),
  `thinning`, `tau`, `binomial`, `hybrid`, generated rules kernels and
  generated event loops (`execute_with_rules_loop`), and
  `Prepared::bind_parameters` candidates (the candidate's request is the
  template's, re-serialized, so it keeps `outputs.quantiles`).
- Schemas: `rules::schema` removes the 1.0-only parts from the `rules-*`
  schemas (`ONLY_1_0`: `Outputs.quantiles`, `Statistics.quantiles`;
  `ONLY_1_0_DEFINITIONS`: `Quantile`), so `rules-run` and `rules-result` stay
  byte-identical. The run-1.0 and result-1.0 descriptions name the new
  fields. `schemas/run-1.0.schema.json` and `result-1.0.schema.json` were
  regenerated with `crn-craft describe`; the other 22 schema files compare
  byte-identical with `describe` output.
- Discovery: `describe()["contracts_1_0"]["quantiles"]` =
  `{"request": "outputs.quantiles", "result": "statistics.quantiles",
  "max": 16, "definition": "numpy-linear-type-7"}`. Rustdoc of
  `RunVersion::V10`, `ResultVersion::V10`, `RunRequest`, `RulesResult` and
  the contract table name the fields (all `concat!` docs, outside the 0.x
  schemas).
- `docs/diagnostics.md`: the new `CRN-R020` cases.

### Python (`python/crn_craft/rules.py`)

- `ensemble(..., quantiles=None)` and `prepare(..., quantiles=None)`: a
  sequence of probabilities (a NumPy array works; a scalar, a string, a
  mapping, booleans or non-numbers raise `TypeError`).
- `None` or an empty sequence: no quantiles, and the documents are exactly
  those written without the argument (tested byte for byte for `None`, `()`,
  `[]` and an empty array). This is the 0.6.0 behavior (I2).
- With quantiles, both documents are 1.0:
  - the request is written as `crn-craft.run/1.0` with
    `outputs.quantiles` last in `outputs` (the method and every other field
    as before; `ensemble` writes `static_direct` as rules-run/0.6 anyway, so
    the paths do not change);
  - a rules-model is written as `crn-craft.model/1.0`, rewriting only the
    schema string: a `Model` or dictionary is serialized as before with the
    `schema` value replaced in place, which equals a byte rewrite of its
    usual serialization; model bytes or text (`rules.load`) get the one
    occurrence of the quoted schema string replaced, every other byte kept
    (more than one occurrence raises `ValueError` asking for a
    dictionary); a model/1.0 document or a document that is no rules model
    is left alone (the engine judges it);
  - `prepare` with a user request: the request (bytes, text or dictionary;
    duplicate keys in bytes or text raise `CrnError` `CRN-E002`) gets
    `outputs.quantiles` (`ValueError` if it already has the field), schema
    `run/1.0`, and the hash of the rewritten model when its `model.sha256`
    was the hash of the supplied bytes. A Static Direct method of
    rules-run/0.2 to 0.5 gets that version's defaults written out
    (`"waiting_time": "inverse"`, `"sparse_trackers": false`, only where
    absent), because run/1.0 has the 0.6 defaults; the test runs a
    rules-run/0.4 request both ways and compares sampler, frames and paths.
- `Result.statistics["quantiles"]`: a dict from each probability (the
  float the result carries) to a DataFrame with the index and columns of
  `statistics["mean"]`; `statistics["median"]` is the same DataFrame object
  as `quantiles[0.5]` when 0.5 was requested.
- Compiled loops serve exact model bytes. A loop compiled for the model/1.0
  document (`d = m.to_dict(); d["schema"] = "crn-craft.model/1.0"`) runs
  quantile requests byte-identically, also with `ensemble(m, ...)` on the
  builder model, since the rewritten bytes are the same. A loop compiled
  for the rules-model bytes is refused with the existing `CRN-R020` at
  `/model`; when the loop's `model_sha256` is that of the bytes the
  quantiles rewrote, the hint says so and names the model/1.0 document
  (`PreparedEnsemble._unified_from`).
- `Result.save` (M6b) needed no change: run against engine quantiles it
  writes the subsections after `sample_variance`, ascending, 0.5 as
  `### median`, and the JSON `statistics.quantiles` list.

### Validation

- `validation/check_rules_python.py`: `check_quantiles` (below); the report
  checks now use the engine's quantiles instead of the hand-made ones of
  M6b (`with_quantiles` removed); `result_bytes` covers the quantile frames.
- `validation/check_installed.py` (`rules_result_report`): the SIR run asks
  for the quantiles 0.025, 0.5 and 0.975; the JSON report's model is the
  model/1.0 document and its request `run/1.0`.
- `validation/fixtures/sir-report.md`: regenerated (see "Fixture diff").
- `validation/check_compiled_loop.py`: a loop compiled for the model/1.0
  document runs quantile requests with 1 and 2 threads byte-identically to
  the interpreter, and `verify()` returns them; the rules-model loop with
  quantiles is refused at `/model`. 53 rows (50 before).
- `validation/merge/check_quantiles.py` (new): the declared exact-law gate
  through the CLI. `validation/merge/measure_quantiles.py` (new): the
  declared measurement.

## Tests written first, and how they failed first

- `crates/crn-craft/tests/quantiles.rs` (12 tests at first) and
  `crates/crn-craft-rules/tests/quantiles_loop.rs` (1) were written against
  the unchanged engine and use only the existing API, inspecting results as
  JSON. All 13 failed: every request with the field was refused with
  `CRN-E003` (unknown field `quantiles`), and the schema test found no
  `quantiles` property.
- `check_rules_python.py`, run against the unchanged `rules.py` (the base
  commit's file in the venv, the new wheel otherwise): `TypeError:
  ensemble() got an unexpected keyword argument 'quantiles'`.
- `check_installed.py` with the quantile request against the old fixture:
  "the SIR report differs".

The tests:

- (a) `quantiles_equal_the_linear_definition_on_the_kept_paths_bit_for_bit`:
  12 probabilities (0.001 to 0.999, including 1/3), 200 replicates with
  every path kept, two observables, `static_direct` with both waiting-time
  policies, `tree`, `next_reaction`, `tau` and `binomial`; the reported
  values equal the definition applied to the kept paths bit for bit, and the
  result without its quantiles equals the result of the same request
  without the field (apart from `run_sha256`).
  `every_stored_request_as_run_1_0_with_quantiles_matches_its_paths`: the 17
  stored requests (never edited) rewritten as run/1.0 with all paths kept
  and five probabilities, which covers `tree`, every `static_direct` option,
  `next_reaction`, `tau`, `mtdna_direct` and `hybrid` on their models;
  `thinning_with_atp_reports_quantiles_of_every_column` (with `ATP`);
  `a_single_replicate_and_a_stop_give_the_held_values`;
  `widely_spread_values_take_numpys_upper_branch` (added after the mutation
  check below): ten SIR epidemics, 16 probabilities, where the formula as
  written differs from NumPy's evaluation in some cells, and the test
  asserts that it does.
- (b) `quantiles_do_not_depend_on_the_thread_count` (1, 2, 4 threads, 500
  replicates, whole result equal apart from `run_sha256`) and
  `a_generated_loop_reports_the_interpreters_quantiles` (generated loop
  against interpreter, 1, 2 and 4 threads, both waiting-time policies).
- (c) `old_request_versions_refuse_the_field` (rules-run/0.1 to 0.6, with a
  list and with an empty list; `null`), `invalid_probability_lists_are_refused`
  (16 accepted; 17, 0, 1, −0.5, 1.5, repeated, decreasing, a subnormal after
  0.1), `the_quantile_buffer_counts_against_the_workspace_limit`
  (1,000,000 replicates × 150 report times × 2 observables refused at
  `/outputs/quantiles` with the three remedies in the hint; the same
  without quantiles admitted; 100 report times admitted).
- Also: `quantiles_without_other_statistics_and_without_quantiles`,
  `quantiles_are_withheld_with_the_other_statistics` (event limit,
  output-limit trim calibrated from the result size),
  `bound_candidates_keep_the_quantiles`,
  `the_1_0_schemas_describe_the_quantiles_and_the_rules_schemas_do_not`,
  and four unit tests in `rules/quantile.rs` (NumPy values for small
  samples, the upper branch, NaN and infinities, column layout).
- (d) `check_rules_python.py` `check_quantiles`: for the SIR builder model,
  the ten-epidemic case (asserting that it tells the two evaluations
  apart), `tree`, `static_direct` (default and inverse), `next_reaction`,
  `tau`, `mtdna_direct`, `binomial`, `thinning`, `hybrid`, `prepare` and a
  bound candidate, `statistics["quantiles"][p]` equals
  `np.quantile(result.paths.values, p, axis=0)` byte for byte, with the
  index and columns of `"mean"`, and `"median"` is `quantiles[0.5]`; plus
  documents, thread independence (1, 2, 4), withholding, the refusals and
  the compiled-loop hint (stand-in loops, no cargo needed).
- (e) `validation/merge/check_quantiles.py`, declared first (below).

### Mutation check

Five mutations of the committed engine, each tested against `quantiles.rs`
and `contracts_v1.rs`, then reverted:

| Mutation | Caught by |
| --- | --- |
| formula as written instead of NumPy's evaluation | at first **nothing** in the integration tests (only the unit tests); after adding the ten-epidemic test: `widely_spread_values_take_numpys_upper_branch` |
| old withholding rule (no statistics without `statistics`) | `quantiles_without_other_statistics_and_without_quantiles` |
| buffer not counted in the admission | `the_quantile_buffer_counts_against_the_workspace_limit` |
| no version gate | `old_request_versions_refuse_the_field` |
| no schema filter | `schemas_match_their_snapshots_and_list_only_their_versions` (M3) and the new schema test |

Why the first mutation slipped through: when adjacent order statistics
differ by at most 1, as they do for 200 replicates of small counts,
1 − g is exact for g ≥ 0.5 and both evaluations round the same real
number. They differ only when x[j+1] − x[j] ≥ 2 and the products round
differently, which needs few replicates with widely spread values.

## Declared check and measurement

Declared in `validation/merge/RESULTS.md` and committed (`8dca966`) before
any gate run. The declared failure margin was wrong in the first commit
(0.024 and 48 standard errors, from one tail of one probability); it was
corrected in `5f34dbd`, still before any run, to 0.0178 and 36 standard
errors over both tails of all three probabilities. The gate and its
tolerance did not change.

- **Exact law (e)**: immigration–death (k_in = 10, k_out = 1, X(0) = 0) as
  `model/1.0`, 100,000 replicates, seed 20261111, t = 10, quantiles 0.025,
  0.5, 0.975 under `static_direct`, `tree`, `next_reaction`, through the CLI.
  Declared: |q̂ − Q(p)| ≤ 1 against the Poisson quantiles 4, 10, 17;
  0 failures expected. **Result: 9 of 9 pass, all exactly equal to Q(p)**;
  quantiles at t = 0 are 0; requests and results valid against the 1.0
  schemas; report byte-identical under `python -O`.
- **Measurement** (no target): 4,000 cells × 37 report times × 6 observables
  (three-class CME, `static_direct`, 2 threads, 1.30e9 events per ensemble):
  median `simulate()` 66.21 s without and 64.31 s with quantiles (spread
  about ±3 s; the difference is noise), peak RSS +7,180 KiB against a buffer
  of 7,104,000 bytes (6,938 KiB) plus a 32,000-byte sort column. Full table
  in RESULTS.md, raw samples in `validation/merge/quantile-cost.json`.

## Run history (interrupted session)

- A server restart stopped the task during the first run of
  `check_quantiles.py`; that run wrote nothing. The uncommitted engine,
  Python and test work survived in the worktree; the engine sources were
  unchanged after their last test run, and the CLI binary was newer than
  them.
- The first complete gate run passed the 9 law rows and failed one contract
  check: the script read the CLI's refusal as a bare list, while the CLI
  prints `{"success": false, "diagnostics": [...]}`. The engine's refusal
  was correct (exit 1, `CRN-R020` at `/schema`). The parser was fixed and
  the gate rerun with the same seed: the 9 rows are byte-identical to the
  first complete run, and every check passes. Both runs are recorded in
  RESULTS.md.
- Before the final gate the wheel was rebuilt and reinstalled from the
  committed tree. `check_compiled_loop.py`, `check_quantiles.py` and the
  measurement ran before that rebuild, on a wheel and CLI built from the
  same engine sources; the Rust tests, I2 and the other Python checks ran
  after it.

## Fixture diff (`validation/fixtures/sir-report.md`)

`git diff` of the regenerated fixture against M6b's:

- `model_sha256` changed: the model is the builder's document with its
  schema string rewritten to `crn-craft.model/1.0` (verified: the new hash
  is the SHA-256 of `model_bytes({**sir.to_dict(), "schema":
  "crn-craft.model/1.0"})`).
- `run_sha256` changed: the request is `crn-craft.run/1.0` (was
  `rules-run/0.6`) and has `outputs.quantiles`.
- Three new subsections after `### sample_variance`: `### 2.5% percentile`,
  `### median`, `### 97.5% percentile`, each a fenced block with one line
  per observable in the order S, I, R. Their values were checked
  independently: `numpy.quantile` of the trajectory lines printed in the
  same file gives exactly the printed lines.

No other line changed: the sampler line is the same (`run/1.0` has the
Ziggurat default of `rules-run/0.6`), as are the method line, the
trajectories, the terminations, the rule events and the three moment
subsections.

## Checks

Venv `/home/claude/m07/venv-m4e` (wheel built from the worktree root with
`maturin build --release --locked --out dist`), CLI
`target/release/crn-craft` of the worktree, everything under the CPU lock.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 115 (4 new), `quantiles` 13 (new), `contracts_v1` 7, CLI 14, facade 3, generated kernels 3, generated loops 25 (I3), `quantiles_loop` 1 (new), `static_direct_defaults` 9 and `stored_requests` 2 (I1), all other tests, doctests 34 + 1 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass |
| Every schema file against `crn-craft describe KIND` | 24 of 24 equal (22 unchanged, `run-1.0` and `result-1.0` regenerated) |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| `check_contracts.py` | 90 of 90 (the 1.0 snapshots include the regenerated schemas); `-O` identical |
| `check_rules.py` | 47 of 47; `-O` identical |
| `check_rules_python.py` | pass, 25 checks (new: "quantiles equal numpy.quantile of the paths"); `-O` identical |
| `check_installed.py` (`python -I`, from `/tmp`, outside the tree) | pass, 17 checks (`rules_result_report` with quantiles and the new fixture); `-O` identical |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` (`CRN_CRAFT_SOURCE` = worktree) | pass, 22 of 22 blocks, including both compiled-loop blocks; `-O` identical |
| `check_compiled_loop.py --require-cargo` (`CRN_CRAFT_SOURCE`, `CARGO_TARGET_DIR` as in CI) | pass, 53 rows; `-O` identical |
| `validation/merge/check_quantiles.py` | pass, 9 comparisons, 0 failures; `-O` identical |

## Limitations

- The buffer holds every replicate's series, so memory grows with
  replicates × report times × observables (8 bytes each); the 2 GiB
  workspace limit bounds it. A streaming quantile estimator would need less
  memory but would not be exact; not done.
- Sorting runs on one thread after the ensemble. For the measured
  ensemble it is negligible; for 1,000,000 replicates × 100 cells it would
  be about 100 sorts of 10^6 values (seconds).
- NaN values are reported as NaN (NumPy's rule) and serialize as `null`,
  like a NaN mean.
- A compiled loop must be compiled for the model/1.0 document to run
  quantiles; the hint says so, but `compile_loop` has no option to do the
  rewrite itself.
- `prepare(..., quantiles=)` re-serializes a request given as bytes or text
  (compact JSON); the model keeps its bytes apart from the schema string.

## Deviations from the plan

- `statistics.quantiles` is **absent** (not `null`) in a result/1.0 whose
  request asked for none. Reason: a run/1.0 request without quantiles then
  gives exactly the result of rules-run/0.6 apart from schema strings and
  hashes, which M3's equivalence test (`contracts_v1.rs`) checks. When
  requested, the field is withheld together with the other statistics (the
  whole `statistics` object is `null`), which is the plan's rule. The plan's
  phrase "null whenever statistics.mean is null" is read as that rule: a
  request with quantiles but without `mean` gets its quantiles and a `null`
  `mean`.
- The definition is NumPy's evaluation of the plan's formula (two branches),
  not the single expression; this is what the plan's test (d) requires
  ("equals np.quantile exactly"), and the field documentation states both.
- The old-version refusal points at `/schema`, as the other version gates
  do, not at `/outputs/quantiles`.
- The workspace reservation also counts one sort column (n × 8 bytes).
- With `quantiles=` the Python layer writes the model as model/1.0 as well
  as the request as run/1.0 (coordinator instruction; the plan requires
  only run/1.0, which accepts every rules-model). A run/1.0 request written
  by hand may still carry quantiles with a rules-model.
- The top-level default `(0.025, 0.5, 0.975)` of `crn.ensemble` is not
  implemented: the top-level `ensemble` does not exist before M6a. M6a sets
  that default; `crn_craft.rules` keeps "none" (I2).
- The manual section is not written (I2 freezes the manual until M10b);
  a draft follows.

## Migration notes

- New: `outputs.quantiles` in `crn-craft.run/1.0` and
  `statistics.quantiles` in `crn-craft.result/1.0` (`crn-craft describe
  run-1.0`, `result-1.0`); `crn_craft::rules::{Quantile, MAX_QUANTILES}`,
  `Outputs::quantiles` (field and method) and `Statistics::quantiles`.
  Rust code that builds `Outputs` or `Statistics` with a struct literal
  must add `quantiles: None`.
- A `rules-run` request with `outputs.quantiles` is refused (`CRN-R020` at
  `/schema`).
- Python: `rules.ensemble(..., quantiles=...)`,
  `rules.prepare(..., quantiles=...)`, `rules.MODEL_1_0`, `rules.RUN_1_0`,
  `result.statistics["quantiles"]` and `["median"]`. Without `quantiles`
  nothing changes. With it, the documents, the result schema and the hashes
  are the 1.0 ones; the numbers are unchanged.

## For later tasks

- M4a/M4b edit `rules/mod.rs` too: the new code is in `Outputs`,
  `Statistics`, `Quantile`, `validate_request` (one block before the
  duplicate-statistic check), `admission` (an `else if` after the workspace
  check), `Reducer`, `statistics()`, `execute_inner` (reducer construction
  and the statistics condition), `schema()` (the `ONLY_1_0` filter, which
  M4b can reuse for the pattern `count`), the 1.0 descriptions and
  `describe()["contracts_1_0"]`.
- M6a: the top-level `crn.ensemble` default `(0.025, 0.5, 0.975)`,
  `quantiles=()` to switch them off; the top-level path writes model/1.0
  anyway, so compiled loops built through the top level match.
- M10b manual draft ("Quantiles", after the statistics section):

  > `quantiles=(0.025, 0.5, 0.975)` adds the median and percentiles of
  > every observable at every report time:
  > `result.statistics["quantiles"][0.975]` is a DataFrame like
  > `statistics["mean"]`, and `statistics["median"]` is
  > `quantiles[0.5]`. They are NumPy's default quantiles
  > (`numpy.quantile`), computed in the engine from every replicate, not
  > only the kept paths. The run then uses the `crn-craft.run/1.0` and
  > `crn-craft.model/1.0` contracts, so its hashes differ from a run
  > without quantiles while its numbers do not. `result.save("sir.md")`
  > writes them as `### 2.5% percentile`, `### median` and
  > `### 97.5% percentile`.
