# M8b log: delete the network engine

Task M8b of [the plan](../../history/plan-0.7.0-engine-merge.md) (M8, "Delete"), on
the integration checkout `/home/claude/fb2/engine`, branch `merge-0.7`, from
`342c322` (the merge of M8a, M8c and M9). Evidence, reports and logs of every
run below: `/home/claude/m07/m8b/` (`step0/`, `dev/`, `final/`, `gate/`,
`mutation/`, `integrations/`).

## Commits

| Commit | Title |
| --- | --- |
| `11b4d45` | M8b: Run the tutorial's blocks on the unified commands and API |
| `44d3189` | M8b: Delete the network engine, its contracts and the old API |
| `db06f87` | M8b: Format the changed validation scripts with ruff |
| `7da5060` | M8b: Pin absorbing states and thread-count replay for every method |
| `e970617` | M8b: Pin tau leaping at a critical homodimer |
| `0ed9b24` | M8b: Law gate after the deletion |
| (this log) | M8b: Add the task log |

The host restarted once during the task (after the deletion was in the
working tree, before the first deletion commit); work resumed from the
working tree, and the checks below ran on the committed states named.

## Step 0: the merge at 342c322

Built the release CLI and the wheel from `342c322` (wheel installed into the
engine venv), then, under the CPU lock:

| Check | Result |
| --- | --- |
| `cargo test --release --workspace --locked` | pass: 35 result lines, 538 tests, 0 failures (I1 `stored_requests` 2, `static_direct_defaults` 9; I3 `generated_loop` 28, `generated_kernel` 3) |
| `check_installed.py` (`python -I`, from `/tmp`) | 18 checks pass |
| `check_compiled_loop.py --require-cargo` | 58 rows, pass |
| `check_python.py` | 17 rows, pass |
| `check_rules_python.py` | 27 checks, pass |
| `check_contracts.py` | 126 of 126 |
| `check_manual.py --require-cargo` (I2) | 22 of 22 blocks |
| `check_tutorial.py` | **failed** at the bash block of line 152 |

The automatic merge of M8a's and M8c's `check_installed.py` and
`check_compiled_loop.py` is fine. `check_tutorial.py` failed for another
reason: M9 made `examples/tutorial/workdir.py` write `crn-craft.model/1.0`
models and `crn-craft.run/1.0` requests, while `docs/tutorial.md` still ran
the hidden network commands on them (`crn-craft simulate` refused the 1.0
model with `CRN-E032`) and imported `crn_craft._legacy`. Fixed in a separate
commit (`11b4d45`) before any deletion, see "Tutorial" below; afterwards 17
blocks pass, normally and under `-O` with identical reports.

## What was deleted

**Rust (`crates/crn-craft/src`).** `simulation.rs`, `ssa.rs` (except the
Ziggurat), `ssa/selectors.rs`, `ssa/bounds.rs`, `next_reaction.rs`,
`tau.rs`, `rates.rs`, `network.rs`, `parameters.rs`, `validate.rs`,
`contracts.rs`, `codegen.rs`; in `archive.rs` the `crn-craft.archive/0.1`
writer and verifier (`verify_any` and `AnyReport` too, since one version is
left); in `fsp.rs` the network entry point and `fsp-result/0.1`; in
`sampling.rs` the network engine's Direct sampler pieces, its streams, the
xoshiro generator and the `fast` prefix sums; in `lib.rs` `CheckReport`,
`parse_model`, `parse_run`, `check_bytes`, `check_files`,
`describe_network`, the `CompiledNetwork` re-export and the old `schema`
kinds; in `unified.rs` the network `check/0.1` branch and the `AnyCheck`
wrapper (`check_any` and `check_any_files` return `RulesCheck`);
`Kinetics::compile`. The examples `bench_ssa.rs`, `generate_kernel.rs` and
`profile_ssa.rs` (old engine only). `grep -rn "xoshiro" crates` finds
nothing (also case-insensitively).

**Moved or kept.**

- `next_reaction::Clocks` (the flat path's next reaction heap) moved to
  `rules/clocks.rs` (crate-private).
- `ssa.rs` is a 12-line module around the kept `ssa/ziggurat*.rs`;
  `crn_craft::ssa::exponential_from_words` keeps its path (generated loops
  reach it through `LoopHost`). The ziggurat's two unit tests now drive
  `exponential_from_words` with closures (same words, same bounds; the
  budget and stop poll are the caller's).
- `contracts::{ModelRef, Ensemble, Statistic}` (parts of every rules and
  1.0 request) moved unchanged into `crn_craft::rules`.
- The network model reader is now private to the converter:
  `rules/convert/network.rs` holds `Model`, `Species`, `Reaction`, `Scalar`,
  `Rate`, `ModelVersion` and `Convention` (with their doc comments, so the
  `model` schema is unchanged), the strict parse with `CRN-E032`, and the
  model validation of 0.6.0 without its warnings (the converter dropped
  them; `analyze` reports `CRN-W001`–`W003`). The plan says "move into
  `convert.rs`"; a private submodule of `convert` keeps `convert.rs`
  readable.
- `validate::valid_relative_path` moved into `rules` (crate-private).
- `rules::convert::from_network_without_limits` (crate-private) replaces
  `Kinetics::compile`: `Kinetics::from_bytes` of a network model drops the
  population limits before validation, as 0.6.0 did.
- `sampling.rs` keeps `UNIFORM_STEP`, `uniform_closed_open` and
  `poisson_with`; the network stream stays only as a `#[cfg(test)]` random
  source of the selector and sampler unit tests (same words as before).
- `archive.rs` keeps the file helpers of `archive::v1`.
- Every schema is unchanged: the 15 files under `schemas/` equal `crn-craft
  describe KIND` (compared as JSON).

**Every reader of `model/0.1` still works** (`check`, `analyze`, `export`,
`convert`, `Kinetics::from_bytes` for `crn-diffsol` and `check_ode.py`) and
reads through the converter.

**CLI.** The hidden commands `simulate`, `simulate-cached`, `ssa`, `tau` and
`fsp` and `legacy()` are gone; they are usage errors (`CRN-E031`). The usage
error report is now `crn-craft.check/1.0` (was the deleted `check/0.1`).
`verify` reads `archive/1.0` only; `--help` says so. The CLI depends on
`crn-craft` only.

**Python.** `python/crn_craft/_legacy.py`, its packaged network examples
(`python/crn_craft/models/legacy/`) and the PyO3 functions and classes only it
used are gone. `__main__.py` lost its hidden commands (its usage report is
`check/1.0` like the native one) and `_common._call`; `importers.py`
validates through `crn_craft.convert` instead of `_legacy.check` (same
diagnostics for an invalid model: the converter returns the network
validation errors; warnings are no longer in the import failure message).
The extension depends on `crn-craft` only.

**Schemas** (13 files): `run`, `result`, `result-v1`, `ssa-run`,
`ssa-result`, `tau-request`, `tau-result`, `archive`, `archive-check`,
`check`, `fsp-result`, `analysis`, `analysis-v1`.

**Models.** The five original examples and their network requests moved
with `git mv` to `models/legacy/` (converter test inputs; the `run.json`
files bind the original bytes by hash). `models/examples/` holds the
converted `crn-craft.model/1.0` files (byte-identical to `crn-craft convert`
of the originals and to the packaged `python/crn_craft/models/*`) and the
state file. The four network `run/0.1` fixtures of `models/invalid/` and
their rows in `expected.json` are gone (their contract is); the model rows
stay. Readers updated in the same commit: the Rust tests, the `flat.rs`
corpus, `integrations/diffsol/tests/{cli,kinetics}.rs` (network inputs, as
M9 asked), `validation/check_{contracts,export,network,ode,python,
rules_python,ssa,structure}.py`, `benchmark_campaign.py` and
`validation/merge/models.json` (paths only; the SHA-256 values are those of
the unchanged network bytes). `kinetics_references.rs` maps the stored
references' old example paths to `models/legacy`; the stored references are
not edited.

**Facade.** The CLI, the extension and the root `Cargo.toml` (workspace
dependency `crn-craft-rules` removed) no longer use `crn-craft-rules`;
`Cargo.lock` and `bindings/python/Cargo.lock` follow. The facade crate stays
a workspace member, a deprecated re-export with its smoke test (3 tests and
the quickstart doctest pass).

**Validation and CI.** `validation/legacy/` was already gone (M8a). The
line-endings job ran the deleted tests `contracts`, `network`, `ssa`,
`archives`; it now runs `--lib --test network_models --test convert --test
structure --test stored_requests --test archives_v1` (all read committed
models and hash them). `validation/rng-vectors.json` is deleted with the
test that included it; `reference_rng.py` checks only the rules stream.
`check_contracts.py` (`REMOVED`, `REMOVED_SCHEMAS`) and `check_python.py`
(`NATIVE_REMOVED`, `REMOVED_COMMANDS`, `crn_craft._legacy` not importable)
assert that the old names are gone. `check_rules_python.py` asserts that an
`archive/0.1` directory is refused.

**Tutorial.** `docs/tutorial.md`: only the runnable blocks changed (prose
left to M10a): `simulate` became `run --archive`, `simulate-cached` a second
identical run, `fsp` and `tau` became `run` with the workdir's 1.0 requests,
and every Python block uses `crn.ensemble`, `Model.fsp`, `Model.analyze` and
`result.paths`/`result.statistics` instead of `crn_craft._legacy`'s
`run_spec`, `simulate`, `fsp`, `tau` and `analyze`.

## Tests: deleted, ported and added

Deleted with the engine: `tests/ssa.rs`, `archives.rs`, `network.rs`,
`ssa_alternatives.rs`, `tau.rs`; the run/0.1 half of `tests/contracts.rs`;
the unit tests of the deleted files. Ported:

- `tests/contracts.rs` → `tests/network_models.rs` (21 tests): every model
  case of 0.6.0 (strict parsing, depth, size, duplicate keys, versions,
  identifiers, stoichiometry, limits, rate-law domains, the invalid
  fixtures) through `convert::from_network` and `check_any`, plus request
  identity (`run/1.0`) and decimal literals through the converter.
- `flat_view.rs`, `analysis.rs`, `rules/convert.rs` and
  `rules/combination.rs` compared with the network engine; they now compare
  with `rules/law_oracle.rs`, the 0.6.0 structure (stoichiometry, reads,
  dependencies of `network.rs`) and rate laws (`rates.rs`, fast branches and
  multiplication order included), kept for tests only. The unit tests include
  it as a module, the integration tests by `#[path]`.
- `kinetics.rs`: `Kinetics::compile` → `Kinetics::from_bytes` (same
  assertions); `contracts_v1.rs`: the network entry-point test is now the
  converter's refusal of 1.0 and rules documents and `CRN-E033` for any
  request with a network model; `fsp.rs`: the six network cases are
  dropped (their `v1_` copies stay) and so is `v1_agrees_with_the_network_fsp`
  (see the gaps); `archives_v1.rs`: `verify` refuses a 0.1 file set;
  `crn-craft-cli/tests/cli.rs`: path tests on `run/1.0` (`CRN-R020` at
  `/model/path` instead of `CRN-E022`), `tau` with a cancel file through
  `run`, removed commands and schemas are usage errors.

Added to `tests/boundaries_v1.rs` (properties of the deleted tests that
apply to the unified engine and were not pinned; they document existing
behavior, so they passed at once on the unchanged engine, except
`unstarted_replicates_are_accounted_for`, see below): the replicate stream
against `validation/rng-vectors-rules.json` (the independent ChaCha8),
report-grid independence, rules without reactants or products as events
under the cap, kept paths are the first replicates and carry the
statistics, partly completed ensembles report no statistics, unstarted
replicates are accounted for, absorbing states draw no word, every method
replays across thread counts, tau at a critical homodimer.

`unstarted_replicates_are_accounted_for` failed first as written from
`tests/ssa.rs` (every replicate unstarted after a pre-cancel): the unified
engine marks the first replicate of each worker as `started` (it stops at
its first poll, with 0 events and 0 words). That is the rules engine's
0.6.0 behavior and stays (M8b must not change results); the test pins it
(at most `threads` started, no event, no word, every replicate listed).

### Gap analysis of the deleted tests

"Pinned" names where the property is held now; "n/a" marks properties of
removed features.

| Deleted test (file) | Now |
| --- | --- |
| ChaCha8 known answers (ssa) | `boundaries_v1::the_replicate_stream_matches_independent_known_answers`; `check_ssa.py` oracle reproduces trajectories (CI) |
| sampler endpoints, waits, selection, precision (ssa) | `sampling` doctest (endpoints); `tests_review::a_zero_uniform_does_not_fail_the_static_wait`, `static_resolution_failures_use_the_resolution_code`; `boundaries_v1` resolution; `clock.rs` (event times) |
| identical payload at 1, 2, 8 threads (ssa); replay across threads (ssa_alternatives) | `rules::tests::replay_is_independent_of_threads`, `selectors_do_not_depend_on_the_thread_count`, `boundaries_v1::every_method_replays_across_thread_counts`, `check_ssa.py`, `check_python.py` |
| statistics equal the retained paths (ssa) | `boundaries_v1::kept_paths_are_the_first_replicates_and_carry_the_statistics`, `check_ssa.py`, `quantiles.rs` |
| save density (ssa) | `boundaries_v1::a_denser_report_grid_changes_no_event_or_common_state` |
| absorption, zero rates, one-replicate SE (ssa) | `selectors::absorbing_states_end_every_selector`, `boundaries_v1` (kept paths test, one replicate) |
| event cap, mixed completion (ssa) | `boundaries_v1` event-cap tests and `a_partly_completed_ensemble_reports_no_statistics`; `rules::tests::budgets_and_cancellation_withhold_statistics` |
| cancellation and wall time, unstarted replicates (ssa) | `boundaries_v1::unstarted_replicates_are_accounted_for` (with the difference above); `check_ssa.py` wall-time stop (CI) |
| population stop, count overflow (ssa); transactional updates (network) | n/a for population limits (the converter refuses them); `flat::a_refused_update_leaves_the_counts_unchanged`, `higher_order::overflowing_combination_counts_fail_the_replicate_with_a_coded_diagnostic` |
| underflow and dynamic range (ssa) | dynamic range: `R031` tests; a positive propensity that underflows is 0 in the unified engine (M2, M5b), not refused: n/a |
| null reactions count as events (ssa) | `boundaries_v1::a_rule_without_reactants_or_products_counts_events_and_meets_the_cap` |
| memory, work and FSP requests refused before running (ssa) | `rules::tests::diagnostics_paths_and_admission`, `fsp::v1_fsp_requests_are_recognized_and_refused_by_prepare` |
| admission before dependency compilation, dependency edge limit (ssa, network) | n/a (no dependency-edge limit in the unified engine) |
| live output cap (ssa); output budget of the compact JSON (ssa_alternatives) | `rules::tests::budgets_and_cancellation_withhold_statistics`, `quantiles::quantiles_are_withheld_with_the_other_statistics` |
| existing and live cancel files (ssa) | `cli::tau_accepts_a_cancel_file_and_stops_before_any_replicate`, `fsp::v1_solve_files_resolves_the_model_and_honors_a_cancel_file`, `boundaries_v1::a_cancel_during_a_long_run_stops_it_within_the_poll_interval` (every method, `rssa_cr` included), `check_ssa.py` live cancel file (CI) |
| lost channel, underflowed wait (ssa) | `boundaries_v1` resolution, `clock::run_1_0_still_fails_a_runaway_rate` |
| retained ids (ssa) | `boundaries_v1` kept paths test |
| cached Direct (4 tests, ssa) | n/a (cached Direct removed in M4c/M8a) |
| saved rows keep species order (ssa) | observables are named and sorted; `check_ssa.py` enzyme conservation on every path (CI) |
| workspace reservations (ssa) | `quantiles::the_quantile_buffer_counts_against_the_workspace_limit`, `tests_feedback::reservation_covers_a_full_log` |
| homodimer, heterodimer, high order, subnormal, overflow, binomial symmetry (network) | `combination.rs` tests, `rules/convert.rs::mass_action_matches_the_network_kernels` (against the oracle) |
| Hill and Michaelis–Menten laws, near threshold, extremes (network) | `rules/convert.rs::hill_propensities_match_the_network_law`, `michaelis_menten_is_bit_identical` (against the oracle; the log1p branch has no counterpart, M8a) |
| reactant gate before evaluation (network) | `tests_feedback::whole_rule_propensities_need_enough_molecules`, `rules::tests::extinction_does_not_evaluate_undefined_rates` |
| columns, signs, reads, dependencies (network) | `flat_view::converted_models_have_the_network_structure_and_propensities` (14 models, against the oracle) |
| a no-change reaction invalidates no propensity (network) | **gap**: no corpus model has such a reaction; only `flat::index_based_lowering_equals_the_quadratic_construction` covers the lowering |
| state domains, named states (network) | `analysis.rs` and `cli.rs` state errors (`CRN-E041`, `E002`, `E003`) |
| compiled parameters, typed models, indices, counts above 2^53 (network) | n/a (no public compiled network; counts are at most 2^53 since M4a) |
| alternative-SSA replay, strict parsing, sampling limit (ssa_alternatives) | `selectors.rs`, `rules::tests::strict_request_input`, `budgets_and_cancellation_withhold_statistics`, `check_alternative_contracts.py` (CI) |
| absorbing network uses no sampling budget, RSSA-CR absorption (ssa_alternatives) | `boundaries_v1::absorbing_states_draw_no_random_words` (every method) |
| `fast`, legacy inverse paths, `ssa-run` defaults, discovery 0.11 (ssa_alternatives) | n/a (removed, D2 and D4) |
| archive tests (archives) | `archives_v1.rs` (the same ten names on `archive/1.0`); the 0.1 memory-admission test: n/a, its 1.0 analog `large_report_grids_with_many_threads_remain_verifiable_when_cancelled` |
| tau: approximate assessment, non-negativity, repeatability, critical homodimer (tau) | `boundaries_v1::the_sampler_identity_says_which_methods_are_approximate`, `tau_never_passes_the_last_reachable_count_of_a_critical_homodimer`, `tests_review::tau_leaps_respect_the_vacancy_capacity`, `check_tau.py` (CI) |
| tau: direct steps at critical counts, halved rejected leaps, capped leap means, distinct keys (tau) | n/a: rules of the removed `tau-request/0.1` algorithm (the unified CGP leap has one competing critical event and one stream, D4) |
| network FSP equals the 1.0 path (fsp `v1_agrees_with_the_network_fsp`) | **weaker**: `validation/merge/check_fsp_references.py` against the stored 0.6.0 results (not in CI; 0 failures here) |
| run/0.1 contract tests (contracts) | n/a; the `run/1.0` strictness is in `rules::tests::strict_request_input`, `check_contracts.py` and `fsp::v1_requests_are_checked_with_coded_diagnostics` |
| unit tests of `simulation.rs`, `parameters.rs`, `rates.rs`, `network.rs`, `tau.rs`, deleted `sampling.rs` tests | open uniform and cancel-file rate limit: `tests_review`; poll per event: `boundaries_v1` cancel test; workspace reuse and parameters: `check_alternative_contracts.py`, `check_python.py`, `quantiles::bound_candidates_keep_the_quantiles`; memo, xoshiro, fast prefix, tau keys: n/a; selector and bounds tests: ported earlier to `rules/select.rs` and `rules/rssa.rs` |

## Checks

Rust on `0ed9b24` (the engine code is that of `44d3189`; later commits add
tests, the gate report and ruff formatting), under the CPU lock; Python on
the CLI and wheel built from `44d3189`, venv `/home/claude/fb2/engine-venv`
(Python 3.12, NumPy 2.3.5, pandas 3.1.1, SciPy 1.17.0, GillesPy2 1.8.3),
`CRN_CRAFT_SOURCE` and `CARGO_TARGET_DIR` set to the checkout, a fresh
`CRN_CRAFT_CACHE_DIR`, DSMTS at `e79fd58`. Every Python script ran normally
and under `python -O`; all 22 report pairs are byte-identical. Reports:
`/home/claude/m07/m8b/final/`.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass; also `bindings/python` with `-D warnings` |
| `cargo test --release --workspace --locked` | pass: 30 result lines, 397 tests, 0 failures (538 at step 0; the old engine's tests left). I1 `stored_requests` 2, `static_direct_defaults` 9; I3 `generated_loop` 28, `generated_kernel` 3 (fixtures and stored requests untouched); facade 3 + doctest |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (only the known output-filename collision warning) |
| `cargo build --release --workspace --locked`; `maturin build --release --locked --out dist`, wheel installed | pass |
| Schemas against `crn-craft describe KIND` | all 15 files equal; the 13 removed kinds are usage errors |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass, 123 files |
| `check_dependencies.py` | pass |
| `check_contracts.py` | 139 of 139 (126 at step 0: removed-name rows added, legacy rows replaced) |
| `check_network.py`, `check_structure.py` | 219 of 219, 309 of 309 |
| `reference_rng.py` | pass, 9 streams, 360 words |
| `check_ssa.py`, `check_cme.py` | 17 of 17, 14 of 14 |
| `check_export.py` | 9 of 9 |
| `check_fsp.py` (DSMTS) | 16 rows, largest distribution error 4.1e-14 |
| `check_tau.py` | pass |
| `check_rules.py`, `check_mtdna.py` | 47 of 47; 151 checks, 32,768 paths |
| `check_python.py` | 17 rows (removed names, native functions, `_legacy` and commands asserted gone) |
| `check_rules_python.py` | 27 checks |
| `check_compiled_loop.py --require-cargo` | 58 rows |
| I2 `check_manual.py --require-cargo` | 22 of 22 blocks |
| `check_builder.py`, `check_hybrid.py` | 18 tests; 11 checks |
| `check_alternative_contracts.py`, `check_alternative_ssa.py` | 9 checks; 168 of 168 |
| `check_importers.py` | 49 rows |
| `check_tutorial.py` | 17 blocks |
| `check_installed.py` (`python -I`, from `/tmp`) | 18 checks |
| `mutation_check.py --lock` (once, on the tree of `7da5060`) | **47 of 47 caught**, complete |
| `integrations/diffsol`, `integrations/jobs`, `validation/parameter-optimization/worker` (`--locked`, temporary target, removed afterwards) | build; tests 8 + 2, 3, 0 pass |
| `check_ode.py` (the rebuilt `crn-diffsol`, not in CI) | 50 of 50; `-O` identical |
| `validation/merge/check_fsp_references.py` (stored 0.6.0 FSP results) | 11,278 checks, 0 failures, largest probability difference 1.4e-17 |
| `validation/merge/check_analysis_references.py` | 194,388 leaves, 0 failures |
| `validation/merge/check_export_references.py` (0.6.0 binary as reference) | 9 artifacts, 160 propensities, 0 failures |

## Law gate after the deletion

The plan's risk table requires it after M8. `validation/merge/check_laws.py`
with the release CLI backend (a copy of the `44d3189` binary,
`/home/claude/m07/m8b/gate/crn-craft`), `--method '{"kind":
"static_direct"}'`, all twelve models, 2 threads, the gate seeds; one
resumable piece per model (`--models NAME`), each under the CPU lock
(`/home/claude/m07/m8b/gate/`; the CLI backend has no chunk state, so a
model is the unit; vilar-noise, the longest, took about six minutes
including lock waits). The pieces
are combined in `validation/merge/laws/m8b-cli.json`.

| Result | |
| --- | --- |
| Comparisons | 4,814, **0 failures**, largest \|z\| 3.88 (W5), 2,318 cells, 940 at t = 0 |
| Equal to `laws/m4d-auto.json` | yes, for every model: converted model hash, request schema, sampler, the 20 chunks (seeds and events), events, power sums S1–S4, statistics, summary and every row; only `path` differs for the five examples (`models/legacy`), and `models_json_sha256` and `engine` in the header |

The deletion changed no result.

## Removed public items (for the CHANGELOG's Removed section)

**Rust, crate `crn-craft`.**

- Modules: `simulation`, `tau`, `network`, `parameters`, `codegen`,
  `validate`, `contracts`, `next_reaction` (hidden).
- `crn_craft::ssa` except `exponential_from_words`: `Request`,
  `RequestVersion`, `ResultVersion`, `Algorithm`, `WaitingTime`, `Method`,
  `Budget`, `Work`, `Termination`, `Record`, `SsaResult`, `Prepared`,
  `prepare`, `execute`, `execute_with_kernel`, `solve`, `solve_files`,
  `PropensityKernel`.
- Crate root: `CheckReport`, `CompiledNetwork` (re-export), `parse_model`,
  `parse_run`, `check_bytes`, `check_files`, `describe_network`; `schema`
  kinds `run`, `result`, `result-v1`, `archive`, `archive-check`, `check`,
  `tau-request`, `tau-result`, `fsp-result`, `ssa-run`, `ssa-result`.
- `crn_craft::archive`: `ArchiveVersion`, `ArchiveFile`, `Manifest`,
  `ArchiveReport` (0.1), `simulate_files`, `simulate_cached_files`,
  `verify_archive`, `validate_result`, and the unreleased `verify_any` and
  `AnyReport` (use `archive::v1`).
- `crn_craft::fsp`: `FspVersion`, `FspSnapshot`, `FspResult` (0.1), `solve`
  (use `fsp::v1`).
- `crn_craft::sampling`: `SAMPLER`, `stream`, `tau_stream`, `uniform_open`,
  `total_propensity`, `waiting_time`, `event_time`, `select_channel`
  (kept: `UNIFORM_STEP`, `uniform_closed_open`, `poisson_with`).
- `crn_craft::kinetics::Kinetics::compile` (use `Kinetics::from_bytes`).
- Unreleased since 0.6.0: `crn_craft::AnyCheck` (`check_any` and
  `check_any_files` return `rules::RulesCheck`).
- Moved: `contracts::{ModelRef, Ensemble, Statistic}` →
  `crn_craft::rules::{ModelRef, Ensemble, Statistic}` (same schema
  definitions); the network model types are private to the converter.

**Python.** The private module `crn_craft._legacy` (the 0.6.0 top level:
`run_spec`, `simulate`, `ssa_spec`, `simulate_ssa`, `prepare_ssa`,
`prepare_parameters`, `ParameterTemplate`, `mass_action_kernel_source`,
`tau`, `fsp`, `analyze`, `check`, `describe`, `load_example` of network
bytes, ...) and its packaged network examples; the native functions
`legacy_describe_json`, `check_json`, `check_files_json`, `simulate_json`,
`simulate_cached_json`, `prepare_json`, `prepare_ssa_json`, `ssa_json`,
`ssa_files_json`, `prepare_parameters_json`, `mass_action_kernel_source`,
`fsp_json`, `tau_json`, `tau_files_json`, `archive_json` and the classes
`PreparedRun`, `PreparedSsaRun`, `ParameterTemplateRun`,
`BoundParameterRun`; `python -m crn_craft simulate|simulate-cached|ssa|tau|fsp`.

**Command line.** `simulate`, `simulate-cached`, `ssa`, `tau`, `fsp` (hidden
since M7) are usage errors; `check MODEL RUN` with a network `run/0.1`
request gives `check/1.0` with `CRN-E033` instead of `check/0.1`; `verify`
refuses `archive/0.1` directories (`CRN-E065`); the usage error report is
`crn-craft.check/1.0`; `describe/1.0` drops `legacy.archives` and the
`archive`, `archive-check` kinds.

**Contracts and schemas.** No longer read: `crn-craft.run/0.1`,
`crn-craft.ssa-run/0.1`–`0.4`, `crn-craft.tau-request/0.1`,
`crn-craft.archive/0.1`. No longer written: `crn-craft.result/0.2` (and the
reserved `result/0.1`), `crn-craft.ssa-result/0.1`–`0.4`,
`crn-craft.tau-result/0.1`, `crn-craft.fsp-result/0.1`,
`crn-craft.check/0.1`, `crn-craft.archive-check/0.1`,
`crn-craft.describe/0.11`. Schema files removed: `run`, `result`,
`result-v1`, `ssa-run`, `ssa-result`, `tau-request`, `tau-result`,
`archive`, `archive-check`, `check`, `fsp-result`, `analysis`,
`analysis-v1`. Still read: `crn-craft.model/0.1` through the converter
(`check`, `analyze`, `export`, `convert`, the kinetics).

## Migration notes

- Network requests: convert the model (`crn-craft convert MODEL`) and write
  a `crn-craft.run/1.0` request; `simulate` → `run --archive DIR`, `ssa` →
  `run` with `static_direct` and a `selector` or `next_reaction`, `tau` →
  `run` with `{"kind": "tau", "epsilon": ...}`, `fsp` → `run` with an `fsp`
  method.
- Rust: `crn_craft::contracts::{ModelRef, Ensemble, Statistic}` are
  `crn_craft::rules::{...}`; `Kinetics::compile(&model)` →
  `Kinetics::from_bytes(bytes)`; `archive::verify_archive` and `verify_any`
  → `archive::v1::verify`.
- Python: `crn_craft._legacy` is gone; the M6a table maps every name to the
  unified API.
- Repository paths: the original network examples are in `models/legacy/`;
  `models/examples/` holds `crn-craft.model/1.0` files.

## Line counts

`crates/*/src`, every `.rs` file:

| | `crn-craft` | `crn-craft-cli` | `crn-craft-rules` | Total |
| --- | ---: | ---: | ---: | ---: |
| Before (`342c322`) | 36,056 (56 files) | 897 | 204 | 37,157 |
| After | 28,490 (48 files) | 788 | 204 | 29,482 |

The tests under `crates/*/tests` went from 16,011 to 13,752 lines. The
whole task: 113 files, +3,082 / −17,677 lines (`git diff --stat 342c322`
at `0ed9b24`).

## Deviations from the plan and the prompt

- `state` schema kept (the plan lists it for deletion): since M5b it is the
  label-based state of `crn-craft analyze MODEL STATE`, listed by
  `describe/1.0`, and `check_export.py` and `compare_julia.py` need it.
- `model` schema kept: it is the network `model/0.1` the converter still
  reads (generated from the private types).
- `validation/dependency-policy.json` unchanged: it lists workspace
  packages, and the facade stays a workspace member, so `crn-craft-rules`
  must stay allowed (`check_dependencies.py` passes).
- The network types live in `rules/convert/network.rs`, a private
  submodule of `convert`, not in `convert.rs` itself.
- The tutorial needed more than its one `_legacy` Python block: the shared
  `crn` namespace of nine Python blocks and five bash blocks ran removed
  commands (and already failed at 342c322, see step 0).
- The `models/invalid` network run fixtures were deleted with their
  contract (the plan does not mention them).
- `rng-vectors.json` was read by `reference_rng.py` too, only to cross-check
  the deleted Rust test; both go, and a Rust test now pins the rules stream
  against `rng-vectors-rules.json`.
- `crn_craft::archive::verify_any` and `AnyReport` (added in M5d, never
  released) are removed rather than kept as a one-version dispatcher.

## Limitations and open issues

- `validation/merge/kinetics-reference/` builds only against the 0.6.0 crate
  (it says so; it records 0.6.0 values for the M0c references) and no
  longer builds against this checkout.
- Docs still describe the network engine and the removed commands (M10):
  `README.md`, `RELEASING.md`, `docs/*.md` (`contracts.md`,
  `direct-ssa.md`, `alternative-ssa.md`, `network-api.md`,
  `reusable-parameters.md`, `archive-contract.md`, `architecture.md`, ...),
  the tutorial prose, `python/README.md`, `crates/crn-craft/README.md`,
  `crates/crn-craft-cli/README.md` and the facade's crate docs.
- `CRN-E032`'s hint still names `rules-check` and `rules` (aliases of
  `check` and `run`); unchanged to keep the diagnostic text.
- `mutation_check.py` still lists the deleted test file `tau` in
  `TEST_FILES` (skipped automatically); `network_models` is not in the
  list.
- One gap stays (table above): no test has a reaction without net change,
  so "it invalidates no propensity" is pinned only through the lowering
  construction. The FSP comparison with the network engine is a manual
  reference gate now, not a CI test.

## For later tasks

- M10: the docs listed above; the tutorial prose around the ported blocks;
  `docs/contracts.md` and `archive-contract.md` describe removed contracts.
- M11: the facade's dependencies and the `crn-craft-rules` entries of the
  lockfiles can be trimmed (M8c note); `publish-crates.yml` still publishes
  the facade, as planned.
