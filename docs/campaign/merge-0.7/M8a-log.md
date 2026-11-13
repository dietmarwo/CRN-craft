# M8a log: validation scripts on the unified engine

Task M8a of [the plan](../../history/plan-0.7.0-engine-merge.md) (M8, "validation
scripts"), on the integration checkout `/home/claude/fb2/engine`, branch
`merge-0.7`, from `21f71b8`. Evidence, reports and logs of every run below:
`/home/claude/m07/m8a/` (`final/` holds the cargo runs).

## Result

No script under `validation/` that CI or the docs use imports
`crn_craft._legacy` or calls a hidden old command (`simulate`,
`simulate-cached`, `ssa`, `tau`, `fsp` with a network run, `check` with a
network run, `describe_network`). Every script now uses `import crn_craft as
crn` or `crn_craft.rules` and the unified CLI (`run`, `check`, `convert`,
`describe`, `verify`) with converted `crn-craft.model/1.0` models and
`crn-craft.run/1.0` requests. `check_python.py` and `check_installed.py` only
name the removed functions to assert that they are gone, and
`check_contracts.py` lists the hidden commands to assert that they are
hidden.

## Commits

| Commit | Title |
| --- | --- |
| `60e4c19` | validation: exact laws on the unified engine (ssa, cme, alternatives, rules, mtdna, reference stream) |
| `cb9b075` | validation: features on the unified API (fsp, tau, contracts, importers, jobs, parameters) |
| `8abee79` | M8a: Python and installed checks without crn_craft._legacy; CI steps |
| `c68a398` | validation: port benchmarks to the unified run command |
| `3f62cf1` | validation: port comparison scripts to the unified API |
| `8117839` | validation: port Julia comparisons and selector timing to the unified API |
| `ef899e5` | validation: delete retired gap-10 campaign drivers |
| `86a0e87` | tests: pin the unified engine boundaries the old engine tests covered |
| `02c5b27` | validation: re-point the mutation gate to the unified engine |
| (this log) | M8a log |

## Decisions that change what a report contains

**D4, the stream.** The unified engine draws from ChaCha8 keyed by
SHA-256(`crn-craft.rules/0.1` and the root seed as 8 little-endian bytes),
the replicate index as stream id, little-endian 64-bit words. The old
network domain was `crn-craft.chacha8/0.1` and one NUL byte, and the old
`fast` policy used xoshiro256++. Every seed in the scripts is unchanged, so
every trajectory is a new one; no test was changed to make it pass. Where a
script compared against a stored sample or digest, the reference was
regenerated from the unified engine. The independent oracle in `check_ssa.py`
(Python, no engine code) reproduces the engine's trajectories exactly for 1, 2
and 8 threads: the wait is `-log(1 - u) / total` with `u = (word >> 11) /
2^53`, redrawn when 0, and the channel is selected by subtracting rates from
`u * total`.

**Reference vectors.** `validation/reference_rng.py` now holds the ChaCha8
reference for the rules stream (`DOMAIN = b"crn-craft.rules/0.1"`) and its
known vectors in the new `validation/rng-vectors-rules.json`. The old
`validation/rng-vectors.json` stays, and `reference_rng.py` still checks it,
because the Rust test `crates/crn-craft/tests/ssa.rs` includes it until M8b
deletes the old engine and that test.

**Law gate size.** `check_alternative_ssa.py` ran 210 positive-time tests
(21 per variant, 10 variants). The variants `cached_direct` and `fast` (the
removed algorithm and the removed xoshiro policy) are gone, so there are 8
variants (4 selectors times inverse and ziggurat) and `TESTS = 168`. The
sample size (N = 100,000), the family-wise error rate (0.001) and every
check stay; the DKW bound `sqrt(ln(2 TESTS / ALPHA) / (2 N))` is slightly
tighter because `TESTS` is smaller. No tolerance was loosened. The power
rows (a 1-in-1000 shift is detected) pass.

**Timing key.** An `archive/1.0` has no `provenance.json`, so the benchmarks
that read `core_seconds` from it now report `process_seconds`, the wall time
of the `crn-craft run` process (start, run, print). Report schemas that
changed: `performance/0.2`, `ssa-validation/0.2`, `fsp-validation/0.2`,
`tau-validation/0.2`, `gillespy2-comparison/0.3`. Comparisons of process
seconds with older core seconds are not like for like.

## Scripts

Status: ported (P), deleted (D), unchanged (U). Check: result of the run under
the CPU lock, normal and `python -O` with byte-identical reports unless noted.

| Script | Status | What happened | Check |
| --- | --- | --- | --- |
| `reference_rng.py` | P | ChaCha8 rules-stream reference; `__main__` checks both vector files | passes (CI step) |
| `check_ssa.py` | P | Converted fixtures, `run --archive`, independent `archive/1.0` checker, oracle with `-log(1-u)`, cancellation through `run --archive --cancel-file`, wall-time stop exits 1 with `wall_time_limit`; no `--command` | 17 passed, 0 failed; `-O` identical; cancel latency 1.6 ms |
| `check_cme.py` | P | `network_model`, `convert_network`, `run_request`; same exact-law tests | 14 passed, 0 failed; `-O` identical |
| `check_alternative_ssa.py` | P | Python API, 4 selectors times 2 waits, `TESTS = 168`, `--fast` removed | 168 tests, 0 failures, power rows pass; `-O` identical |
| `check_alternative_contracts.py` | P | run/1.0 and result/1.0 schemas, CLI and Python identity, budgets, pre-cancel, unknown, null and duplicate fields, default is linear with ziggurat | all 9 cases pass; `-O` identical |
| `check_rules.py` | P | `network_dimerization` through convert and `run` | 47 passed, 0 failed; `-O` identical |
| `check_mtdna.py` | P | `network_sample` through convert and `run` (next reaction) | passes; `-O` identical |
| `check_fsp.py` | P | Replaces `check_fsp_rules.py` (renamed); `run` with fsp requests; DSMTS through `python -I .../isolated.py` | passes, max error 4e-14; `-O` identical |
| `check_tau.py` | P | Replaces `check_tau_rules.py` (renamed); `run` with tau requests | passes; `-O` identical |
| `check_contracts.py` | P | Four run-based invalid cases now generated run/1.0 requests, each expecting CRN-R020 at the path | 126 passed, 0 failed; `-O` identical |
| `check_importers.py` | P | `crn_craft.Model.from_dict(model).analyze()` | passes; `-O` identical |
| `check_jobs.py` | P | `rules.model_bytes`, `rules.CrnError` | not run: needs the external runner binary (see open issues) |
| `check_parameter_binding.py` | P | `rules.model_bytes`, `rules.CrnError`, `crn_craft.convert`; the only parameter-binding script | 18 bindings, 36 worker evaluations, 12 invalid inputs, 12 concurrent; `-O` identical (worker built from `validation/parameter-optimization/worker`) |
| `validate_parameters.py` | P | The parameter-validation gate, `--fast` and `--command simulate-cached` runs removed | import and ruff only (a campaign driver that rebuilds wheels) |
| `check_installed.py` | P | The legacy section replaced by the unified API and CLI: `c.ensemble`, `rules.fsp`, a tau run through the CLI, `Model.from_dict().analyze()`, CRN-E003 on a bad schema | 18 checks pass under `python -I`; `-O` identical |
| `check_python.py` | P | Docstring only: it already covers the unified API and asserts that the network functions are gone | passes; `-O` identical |
| `legacy/check_python_network.py` | D | Tested only the removed network functions of 0.6.0; `check_python.py` covers the unified API, and its CI step is removed | none |
| `check_fsp_rules.py`, `check_tau_rules.py` | D | Renamed to `check_fsp.py` and `check_tau.py`; their CI steps are the new names | see above |
| `check_parameter_interface.py` | D | Tested the network engine's `ParameterTemplate`, removed with that API; `check_parameter_binding.py` covers binding on the unified engine | none |
| `check_network.py`, `check_structure.py`, `check_export.py` | U | Already on the unified `analyze` and `export` through the converter; they read the old network files under `models/` and `examples/`, which M8b moves | 219, 309 and 9 passed; `-O` identical |
| `check_rules_python.py`, `check_hybrid.py`, `check_builder.py` | U | Unified API only | pass; `-O` identical (27 checks; 11 checks; 3 checks) |
| `check_manual.py` (I2), `check_compiled_loop.py` | U | Unified API only | I2 passes; compiled loop passes; `-O` identical |
| `check_tutorial.py` | U | Runs `docs/tutorial.md`, which still imports `_legacy` (M10a) | 17 blocks pass; `-O` identical |
| `check_dependencies.py`, `check_generated_loop.py`, `check_generated_rules.py`, `check_ode.py` | U | No engine API of the old engine; `check_ode.py` is not a CI step and needs the diffsol binary | `check_dependencies.py` passes, `-O` identical; the others not run (I3 covers the generated kernels in `cargo test`) |
| `benchmark_campaign.py` | P | `request_for`, `convert`, `run_case` with `run` and process seconds; schema `performance/0.2` | smoke passes |
| `benchmark_order_two.py`, `benchmark_short_paths.py`, `performance_step.py` | P | `run` replaces `simulate-cached` | smoke passes |
| `benchmark_iterations.py`, `benchmark_fsp.py` | P | run/1.0 ssa, fsp and tau requests; core seconds are process seconds | smoke passes |
| `compare_gillespy2.py` | P | `crn.ensemble` selectors, `--crn-only` smoke form, rows `cli_run_process`, `python_default`, `python_checked`, `python_<selector>`, `gillespy2_ssa_cpp` | full run (below); figure re-plotted |
| `compare_references.py` | P | `--command` default `run`, process seconds | full run: two_channel and sparse_66 (rebop over CRN 1.27 and 21.95) |
| `compare_tau_samplers.py` | P | run/1.0 tau | smoke with one binary twice: KS 0.0228 below the bound 0.0436 |
| `compare_julia.py` | P | `crn_job`, `final_counts`, `crn_payload`, `run_crn`; engines `CRN-linear` (reference) and `CRN-default`; the cached-versus-direct identity check is dropped (cached Direct no longer exists) | crn-craft side only |
| `compare_julia_accelerated.py`, `compare_julia_gap10.py` | P | Engines are the 8 selector/wait variants, `linear/inverse` is the reference of the two-sample gates, `TESTS = 5218` (was 6946, same alpha), seeds by engine index change with the list; the retired 0.6.0 baseline comparison is removed from the gap-10 script | crn-craft side only |
| `measure_alternative_ssa.py` | P | Selector timing and prepared equality through `crn.ensemble` | runs (16 paths, all 4 workloads, 8 variants) |
| `benchmark_mtdna.py`, `benchmark_mtdna_gap.py`, `benchmark_compiled_loop.py`, `benchmark_compiled_string_rates.py`, `benchmark_independent_ode.py`, `benchmark_generic_rules.py` | U | Rules API only | `--help` runs |
| `benchmark_static_direct_defaults.py` | U | Rules API only | one repeat: 1.108 speedup, target 1.1 met |
| `acceleration_workloads.py`, `probe_structural_rules.py` | U | Model generation; the `rules` alias keeps old binaries usable | import |
| `report_gap10.py`, `render_acceleration_report.py`, `plot_julia_comparison.py`, `package_*.py` | U | Renderers and packaging of retained evidence | import |
| `parameter_optimization.py`, `parameter_optimization_models.py`, `report_parameter_optimization.py`, `generate_parameter_workers.py` | U | The worker already uses `crn_craft::rules` | worker builds |
| `mutation_check.py` | P | Anchors re-pointed (below) | 47 of 47 caught |
| `profile_ssa_acceleration.py` | D | Drove the retired old-engine profiler (`crates/crn-craft/examples/bench_ssa.rs`) with cached Direct; prepared equality is in `check_alternative_contracts.py` | none |
| `measure_julia_gap.py`, `run_gap10.py`, `gap10_iteration.py`, `validate_gap10.py`, `diagnose_gap10_regressions.py` | D | Drivers of the closed gap-10 campaign: they ran the old engine through `bench_ssa`, cached Direct and `fast`, and compared with the 0.6.0 package | none |
| `gap10_python_reference.py` | D | Only served the removed 0.6.0 baseline comparison | none |
| `merge/baseline.py`, `merge/references.py` | U | Run on the 0.6.0 package, as the task requires | not touched |
| `merge/*` (the others) | U | Coordinator's suite | not touched |

I read "port the gap-10 scripts" as: port `compare_julia_gap10.py`, keep the
evidence and its renderers (`report_gap10.py`, `render_acceleration_report.py`,
`plot_julia_comparison.py`), and delete the campaign drivers whose subject
(the old engine under test) is gone. The evidence directories under
`validation/julia*`, `validation/ssa-acceleration` and
`validation/parameter-optimization/evidence` stay as they were.

## CI (`.github/workflows/ci.yml`)

- Removed: the cached `check_ssa` and `check_cme` steps, the `check_fsp_rules`
  and `check_tau_rules` steps, the legacy Python network steps, and the `--fast`
  `check_alternative_ssa` steps.
- Added or renamed: `check_fsp.py` and `check_tau.py` (the renamed scripts,
  same flags), the mutation-gate comment now names the unified code.
- Every job keeps its real checks, and each replacement step checks at least
  what the removed one did. A script that imports `crn_craft` runs only in the
  `python` job, which installs the wheel, NumPy, pandas, xarray and
  jsonschema; `check_ssa.py` and `check_cme.py` need only jsonschema, and the
  reference-export job installs SciPy for `check_export.py`.
- Every `validation/*.py` named in the workflow exists.

## Mutation gate

`validation/mutation_check.py` keeps 47 mutants. 30 mutate the unified code
now (21 keep their name, 9 are renamed because the old idea has another
form); 17 mutate code that stays (cone, linkage, conservation, export,
analysis, FSP, kinetics, dual). The runtime test list is the library tests
plus the integration tests that exist from a fixed list, so it survives the
removal of the old engine's tests. The run is resumable (`--resume`), writes
the report after every mutant, can take an advisory lock (`--lock`) around
each cargo run, and builds incrementally.

| Old anchor | Unified anchor |
| --- | --- |
| binomial divisor and symmetry (`rates.rs`) | `rules/combination.rs` |
| product sign, custom rate read, dependency updates (`network.rs`, `simulation.rs`) | `rules/flat.rs` (`changes_of`, `readers`, `affected`) |
| reactant gate (`rates.rs`) | `rules/flat.rs` (`RateCache::channel`) |
| Hill near the threshold (`rates.rs`) | `rules/expr.rs` (Hill denominator); the old log1p branch does not exist any more |
| wait scale, event cap, save rows, selection (`sampling.rs`, `simulation.rs`) | `rules/flat_path.rs` (`simulate_linear`, `scan`, `check_events`, `save_before`) |
| replicate stream, kept paths, variance, conditional statistics | `rules/mod.rs` |
| zero wait | `rules/engine.rs` (`Clock::advance`) |
| lost positive increment (3 mutants), minimum resolution | `flat.rs` `direct_with`, `flat_path.rs` linear loop, `select.rs` tree, `engine.rs` rule check (named `direct_scan_`, `linear_`, `tree_accepts_unresolved_rate` and `rule_resolution_check_removed`) |
| cancel file interval, poller (`simulation.rs`) | `rules/mod.rs` cancel-file test, `rules/engine.rs` generic poll (the unified engine has no cancel interval) |
| archive digest, completion marker, byte budget, unstarted records, initial state, stop time | `archive/v1.rs` |
| tau exactness and event cap | `rules/mod.rs` sampler identity, `rules/leap.rs` |

Two old mutants have no counterpart that a test can reach: a selection
boundary `target < a` versus `<=` and a save at `<` versus `<=` an event
time. They differ only when a draw equals a rate sum or a report time exactly,
which no stream reaches through the public API (the old tests injected
words). They are replaced by `ssa_scan_does_not_subtract` and
`ssa_save_row_not_written`.

**Finding of the first run.** With the anchors moved, 9 of the 47 mutants
survived: the unified tests did not pin the event-cap boundary of the flat
path, the tau event cap, the sampler identity of the approximate methods, the
polling of a long run, the 2^-53 refusals of the linear and tree selectors, or
a stop time past the horizon in an archive. The old engine's tests covered
these. They are pinned now by `crates/crn-craft/tests/boundaries_v1.rs` (5
tests) and one case in `archives_v1.rs` (tamper of `stop_time`); no engine
code changed. The first complete pass over the re-pointed list was 36 caught
and 11 survivors (two of the 11 were weak analogues I replaced before the next
pass); after the tests and the replacement anchors the pass is 47 of 47 caught
(`mutations.json`, `complete: true`).

## Comparisons run with a real peer

`compare_gillespy2.py` ran in full (GillesPy2 1.8.3, 64 paths, one thread,
median of 3, one NumPy-solver repeat), `crn-craft` built from this checkout.
Both validations pass on every workload: the exact-Poisson tests (W1, W2) for
the checked and the ziggurat configuration, and the two-sample tests against
GillesPy2 (smallest KS p-value 0.209). Median seconds:

| Workload | CLI process | default | linear, inverse | composition-rejection | sorting | RSSA-CR | GillesPy2 C++ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| W1 two channels | 0.048 | 0.047 | 0.049 | 0.082 | 0.063 | 0.063 | 0.136 |
| W2 sparse 66 | 0.103 | 0.108 | 0.196 | 0.127 | 0.154 | 0.109 | 1.077 |
| W3 repressilator | 0.290 | 0.276 | 0.273 | 0.553 | 0.371 | 0.432 | 0.970 |
| W4 dense 66 | 0.118 | 0.124 | 0.853 | 1.309 | 1.204 | 0.121 | 1.070 |

The run was on a shared 2-core host, so read the ratios, not the seconds. The
default (`auto` selector, ziggurat waits) is at or near the fastest column on
every workload.

The Julia scripts (`compare_julia*.py`) were not re-measured (no Julia here).
Their crn-craft side was exercised once: `run_crn` for both engines on W1 and
W2, and every selector variant of both gap scripts through `rust_job` on W1
to W3 with replay equality of the one-shot and the prepared form, 8 paths.

## Final standard set

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: 35 test binaries, 537 tests, 0 failures (I1 `stored_requests`, `static_direct_defaults`; I3 `generated_loop`, `generated_kernel`) |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass |
| `ruff format --check validation --exclude '*.md'`, `ruff check validation` | pass (79 files) |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` | pass; `-O` identical |
| Every touched script, normal and `python -O` | reports identical (table above) |

## Migration notes

- `--command` is gone from `check_ssa.py` and `check_cme.py`: the only run
  command is `run`. A script that passed `simulate-cached` now fails at the
  command line.
- A script that read `result["trajectories"]` or `counts` reads
  `result["paths"][i]["values"]`; species are observables of the same name,
  sorted by name; `save_times` is `report_times`, `keep_trajectories` is
  `keep_paths`.
- The request's `model.path` must name the model file next to the request.

## Open issues and hand-off to M8b

- Still using `_legacy` (not scripts under `validation/`): `docs/tutorial.md`
  (so `check_tutorial.py` runs it), `tutorials/optimization/{calibration,
  oscillators,sensitivity}.py`, `examples/circuit/screen.py`,
  `examples/tutorial/make_assets.py`, `python/crn_craft/__main__.py`,
  `importers.py`, `_common.py` and the module itself.
- Hidden commands: only the assertions in `check_contracts.py`
  (`HIDDEN`, `HIDDEN_SCHEMAS`) and `check_python.py` (names asserted absent)
  mention them; they must change with M8b.
- Scripts that read the old network files in `models/` and `examples/`
  (`check_network.py`, `check_structure.py`, `check_export.py`,
  `check_ode.py`, `benchmark_campaign.models`, `acceleration_workloads.py`)
  need their paths updated when M8b moves the files to `models/legacy`.
- `crates/crn-craft/tests/ssa.rs` reads `validation/rng-vectors.json`; delete
  both together. The integration tests `archives`, `network`, `ssa`,
  `ssa_alternatives` test the old engine; `mutation_check.py` does not list
  them.
- `probe_structural_rules.py` and `benchmark_generic_rules.py` keep the
  `rules` alias so that they run against older binaries.
- `check_jobs.py` needs the external fcmaes runner binary; it was ported but
  not run here. `check_ode.py` needs the diffsol binary and is not in CI.
- `validation/parameter-optimization` evidence text still names the deleted
  `check_parameter_interface.py` in places (one sentence of the README was
  updated).
- A selection boundary and a save at an equal time cannot be pinned through the
  public API; a test with an injected random source (the `Words` type is
  generic) could cover them if the mutation gate should include them.
