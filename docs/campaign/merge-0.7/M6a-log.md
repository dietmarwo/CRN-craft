# M6a log: one Python API

Task M6a of [the plan](../../plan-0.7.0-engine-merge.md) (section 3,
"Unified Python API", and M6 items 1 to 5), done together with M7 on the
branch `merge/M6` (worktree `/home/claude/m07/wt-m6`) from `merge-0.7` at
`c98dcfe` with `merge/M5cdg` merged in at `f3fcb6b`.

## Step 0: the coordinator's merge

The coordinator resolved the conflicts of `merge/M5cdg` (CLI `rules` command,
`rules/mod.rs` module lines, `RESULTS.md`) without building. On `f3fcb6b`,
under the CPU lock, before any change of this task:

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: 35 test binaries, 526 tests, 0 failures (I1 `stored_requests`, `static_direct_defaults`; I3 `generated_loop`, `generated_kernel`) |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (the known output-filename collision warning) |
| Every schema file against `crn-craft describe KIND` | 26 byte-identical, all eight 1.0 schemas among them; `analysis` and `analysis-v1` are no `describe` kind since M5b (their files stay until M8b) |
| `check_contracts.py` | 103 of 103; `-O` identical |
| `check_rules_python.py` | 27 checks pass; `-O` identical |
| `check_installed.py` (`python -I`, from `/tmp`) | 17 checks pass; `-O` identical |
| I2 `check_manual.py --require-cargo` (`CRN_CRAFT_SOURCE` = an export of `f3fcb6b`) | 22 of 22 blocks; `-O` identical |

No merge problem was found, so there is no fix commit. The manual check ran
against a `git archive` export of `f3fcb6b` and a copy of its CLI binary,
because this task's edits had started in the worktree meanwhile. Evidence:
`/home/claude/m07/m6a/step0/`.

## What changed

### Top level (`python/crn_craft/__init__.py`)

- `crn_craft.__all__` is exactly: `Model`, `ensemble`, `prepare`, `load`,
  `compile_loop`, `loop_source`, `describe`, `check`, `convert`,
  `load_example`, `read_sbml`, `read_bngl`, `load_report`, `archive`,
  `verify_archive`, `CrnError`, `CompiledLoopError`. Apart from submodules
  (`rules`, `builder`, `importers`, `observables`, `report`) and
  `__version__`, nothing else is public.
- Through these names the documents are `crn-craft.model/1.0` and
  `crn-craft.run/1.0`, `quantiles` defaults to `(0.025, 0.5, 0.975)` and
  `selector=` chooses the Static Direct selector (absent means `auto`).
  `quantiles=()` or `None` asks for none.
- `crn.Model` is a subclass of the builder whose class flag `_unified` is
  true: `to_dict()` writes `crn-craft.model/1.0`, otherwise the same
  document as `crn_craft.rules.Model` (only the schema value differs, in
  place, so the key order and every other byte are those of 0.6.0).
- `crn.ensemble` and `crn.prepare` call the shared implementations in
  `rules.py` (`_ensemble`, `_prepare`) with `unified=True`; the result is a
  `PreparedEnsemble` with `_unified` set. A rules-model gets the schema
  string `model/1.0`; a network `crn-craft.model/0.1` is converted
  (`crn_craft.convert`). `crn.prepare` rewrites a rules-run request as
  `run/1.0` (below), adds the default quantiles unless the request states
  `outputs.quantiles`, and keeps the exact bytes of documents that need no
  change (a `run/1.0` request with quantiles for `model/1.0` bytes) and of a
  request that does not parse (the engine reports it).
- `crn.load(path)` reads a model file of any version, checks it with
  `crn.check` and returns its exact bytes.
- `crn.describe()` is `crn-craft.describe/1.0` (M7), `crn.describe(kind)`
  the schema `crn-craft describe KIND` prints.
- `crn.check(model, request=None)` returns the report of
  `crn-craft check MODEL [REQUEST]` (M7, native `check_any_json`).
- `crn.load_example(name)` returns a `crn.Model` of the converted example
  packaged in `python/crn_craft/models/<name>.model.json` (the bytes
  `crn-craft convert models/examples/<name>.model.json` prints). The network
  files moved to `python/crn_craft/models/legacy/` for
  `crn_craft._legacy.load_example`. The wheel packages both directories
  (the whole `python/crn_craft` tree; `pyproject.toml` needed no change;
  checked with `unzip -l`).
- `crn.read_sbml(path, *, convention)` and `crn.read_bngl(path, *,
  convention)` return `crn.Model` objects of the importers' converted
  `crn-craft.model/1.0`; `model.import_report` keeps the import report
  without its `model` (provenance). The conventions are the importers'
  (`"stochastic_count_propensity"` for SBML, `"binomial"` or
  `"falling_factorial"` for BioNetGen).
- `archive`, `verify_archive`, `compile_loop`, `loop_source` and
  `load_report` are the objects `crn_craft.rules` exports.

### `crn_craft.rules` (I2)

- Every public name `rules.py` had in 0.6.0 is still there (checked by
  `check_python.py` against the list taken from `43105b2`), with the same
  defaults and documents: `rules.Model.to_dict()` writes
  `rules-model/0.3`/`0.4`, `rules.ensemble` the lowest `rules-run` version
  that fits, no quantiles by default.
- New there: `selector=` on `ensemble` and `prepare` (it needs `run/1.0`,
  so passing it writes the 1.0 documents, as `quantiles=` does),
  `FspResult`, `DEFAULT_QUANTILES`, `SELECTORS`, and the new `Model`
  methods (the builder is shared).
- Two fixes in the request rewrite of `quantiles=` (M4e's
  `_request_with_quantiles`, now `_request_1_0`): a rules-run Static Direct
  request rewritten as `run/1.0` gets `"selector": "linear"` written out,
  because `run/1.0` defaults to `auto` since M4c (composition–rejection from
  512 channels) and rules-run requests always scan linearly; and a request
  whose version is no rules-run version (for example `crn-craft.run/9`) is
  no longer rewritten to `run/1.0`, so the engine refuses it. The M4e check
  in `check_rules_python.py` now expects the selector in the rewritten
  method.

### Builder (`python/crn_craft/builder.py`)

- `add_reaction` with two reactant molecules writes the 0.6.0 patterns
  (`x`, `y`), also for `["A", "A"]` (bytes checked against a 0.6.0
  document). Three or more molecules give one pattern per species in
  first-appearance order with `count` (left out when 1), labeled `x`, `y`
  for up to two species and `x1`, `x2`, ... for more, as the converter
  labels network reactions. Such a model writes `crn-craft.model/1.0` also
  through `crn_craft.rules.Model` (a 1.0 feature).
- `Model.from_dict(document)` (dictionary, text or bytes; a network
  `model/0.1` is converted; the engine checks it, `CrnError` otherwise).
  Representation (the task's open choice): the model keeps the document
  and `to_dict()` returns it as given, key order included, until an edit;
  `crn.Model.from_dict` changes only the schema string of a rules-model to
  `model/1.0`, and `from_dict` of a `model/1.0` document sets the instance
  flag. When the builder can represent the document (attribute-free types,
  one count observable per type named after it, one numeric initial count
  per type; explicit defaults such as `"attributes": []`, `"reduce":
  "sum"`, integral float counts and a `description` accepted), its content
  is also loaded into the builder state, so `add_*`, `set_initial` and
  `set_stop` work and the first edit switches `to_dict()` to the builder's
  document. Every converted network model and every packaged example is
  editable this way. Other documents (attributes, filtered or weighted
  observables, ATP, ...) are kept whole and refuse edits with `ValueError`.
- `Model.fsp(end, report_times, bounds, *, max_states=None,
  tolerance=None)` returns `rules.FspResult`: `states` (DataFrame, index
  `state`, one int64 column per species label), `probabilities` (DataArray
  `("time", "state")`), `mean` (the projected first moment, DataFrame by
  time), `truncation` (DataFrame of `physical_leakage_estimate`,
  `series_tail_upper_bound`, `discarded_weight_estimate`, `mass_defect`,
  `roundoff_mass_residual`), `species` and `metadata` (the
  `fsp-result/1.0` document). `rules.fsp` keeps returning the document.
- `Model.analyze()` returns the `crn-craft.analysis/1.0` report.
- `Model.to_gillespy2()` runs the engine's GillesPy2 export and returns
  `build_model()` of the generated module (a `gillespy2.Model`); an export
  refusal raises `CrnError` (`CRN-E050` for a stop condition), a missing
  GillesPy2 `ImportError`.

### The legacy module (`python/crn_craft/_legacy.py`)

The 0.6.0 `__init__.py` moved there (`git mv`), unchanged apart from
importing the shared helpers from the new `python/crn_craft/_common.py`
(`CrnError`, `model_bytes`, `convert`, `_call`, `_crn_error`, `_path`,
`__version__`), `describe` reading the native `legacy_describe_json`
(`crn-craft.describe/0.11`), and `load_example` reading
`models/legacy/`. `import crn_craft._legacy as crn` gives every name of the
0.6.0 top level with its behavior; the PyO3 functions behind them are kept.
`crn_craft.importers` validates with `_legacy.check`, as before. M8b
deletes the module.

### Bindings (`bindings/python/src/lib.rs`)

New native functions: `describe_json` now returns `describe/1.0`;
`legacy_describe_json` (`describe/0.11`), `check_any_json`,
`check_any_files_json` and `run_files_json` (the text `crn-craft run`
prints, for `python -m crn_craft run`). No function was removed.

### Validation

- `validation/check_python.py` rewritten for the unified API (below). The
  0.6.0 script is `validation/legacy/check_python_network.py` on
  `crn_craft._legacy`, run in CI next to it until M8a.
- `validation/check_installed.py`: `check_unified` (below), the SIR report
  through `crn.ensemble` (its documents equal `rules.ensemble(...,
  quantiles=(0.025, 0.5, 0.975))`, so `validation/fixtures/sir-report.md` is
  unchanged), the Python `--version` with its milestone, and the network
  calls switched to `crn_craft._legacy`.
- `.github/workflows/ci.yml`: the legacy script runs in the `python` job,
  normally and under `-O`.

## Tests

The implementation was drafted before the tests. The tests were then run
against the step-0 wheel and binary (`f3fcb6b`, venv
`/home/claude/m07/m6a/venv-step0`) and failed:
`check_python.py` at its first check (`AttributeError: module 'crn_craft'
has no attribute '__all__'`), `check_installed.py` at its imports (`No
module named 'crn_craft._legacy'`). Logs:
`/home/claude/m07/m6a/*-first-fail.log`.

`check_python.py` keeps every check of the old script, on the unified API
and the same rules:

| 0.6.0 check | Unified form |
| --- | --- |
| `CLI_binding_threads_{1,2,8}`: binding results identical for 1, 2 and 8 threads apart from `run_sha256`; the archived CLI result equals the binding's | `crn.ensemble` (Static Direct, 64 replicates, seed 5001) for 1, 2 and 8 threads: metadata (without `run_sha256`), records and rule events identical; `crn-craft run --archive` of the same documents gives the same result, record by record |
| `reusable_parameters_and_codegen`: a bound candidate replays from its documents, repeated binding is identical, booleans are refused, the kernel source is generated | `bind_parameters` with every parameter × 1.2: preparing the candidate's reported documents again (request bound to their hash, M6b) gives the same numbers, records, rule events and paths; binding twice is identical; a boolean is refused (`CrnError` at `/parameters/birth`, through the top-level flag only); `crn.loop_source` returns the generated loop |
| `strict_inputs` | `crn.check` of a document with a duplicate key and of `crn-craft.model/99` is invalid |
| `binding_errors_and_cancel_file` | `crn.prepare` of a `crn-craft.run/9` request raises `CrnError` `CRN-E003` (a `ValueError`); a cancel file gives termination `cancelled` |
| `keyboard_interrupt` | Ctrl-C stops `simulate()` of a top-level ensemble within 2 s |
| synthetic oscillation and first-passage checks, circuit study | unchanged; the study runs the converted repressilator with Static Direct and classifies the first network species (`m1`) |

New checks: `names` (top-level names, removed names absent, every 0.6.0
name of `crn_craft.rules` present, shared objects), `documents_1_0` (0.6.0
documents through `crn_craft.rules`, byte for byte for the SIR builder
model; 1.0 documents and default quantiles through the top level with the
same numbers; `quantiles=None`/`()`; `selector=` through both names with
the composition–rejection sampler; `crn.prepare` of rules-run requests as
dictionary, text and bytes; exact bytes kept; network models converted by
`crn.ensemble` and refused with `CRN-E033` by `crn_craft.rules`), `models`
(pair and `count` patterns, `load_example` of all five examples against
`crn.convert`, editing a converted model, a fixed document refusing edits,
`from_dict` checks and conversion, `read_bngl`), `model_methods`
(`Model.fsp` against the Poisson mean, `Model.analyze`, `to_gillespy2` or
its `ImportError` without GillesPy2), `command_line` (M7).

`check_installed.py` gained `check_unified`: the top-level names and
`describe/1.0`; `Model.fsp` of immigration–death at t = 0, 2, 5, 10 against
the Poisson law (every state within 1e-10 absolute, the mean within 1e-9,
leakage below 1e-12); `Model.analyze()` of the enzyme model (exact
conservation of E + ES = 10 and S + ES + P = 50); a third-order reaction
3A → ∅ (k = 0.01, A(0) = 10): one pattern with `count` 3, its FSP against
the exact pure-death law 10 → 7 → 4 → 1 at t = 0, 1, 5, 20 (within
1e-10), and a Static Direct ensemble (2,000 replicates, seed 11, 2
threads) whose kept paths stay in {10, 7, 4, 1} and decrease, and whose
means at t = 1, 5, 20 lie within 5 standard errors of the FSP means.

Declaration, honestly placed: the ensemble comparison (2,000 replicates,
seed 11, 3 comparisons, 5 standard errors, 0 failures expected) was written
into the script before it first ran, but this declaration was written
after that first run, which passed; nothing was changed after it.

## Checks

One gate for M6a and M7. The validation scripts ran on the CLI and the
wheel built from `383cf1f`, under the CPU lock, venv
`/home/claude/m07/venv-m6` (Python 3.12, NumPy 2.3.5, pandas 3.1.1, rebop
0.10.0, GillesPy2 1.8.3 added for `to_gillespy2` and `check_export.py`),
`CRN_CRAFT_SOURCE` and `CARGO_TARGET_DIR` set to the worktree, DSMTS at
`e79fd58`. The last commit (`83369a0`, one discovery string and one
comment) was then built and checked again: fmt, clippy, all tests,
rustdoc, the wheel, the schemas and the seven Python scripts its change can
reach (last block of the table). Every Python script ran normally and under
`python -O` with byte-identical reports. Reports:
`/home/claude/m07/m6a/gate/` and `/home/claude/m07/m6a/final-gate/`.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass (`missing_docs` is a warning in `crn-craft` and the facade, so this also covers documentation) |
| `cargo test --release --workspace --locked` | pass: 35 test binaries, 532 tests, 0 failures (526 at step 0; new: 4 CLI tests, 2 doctests); I1 `stored_requests` 2, `static_direct_defaults` 9; I3 `generated_loop` 26, `generated_kernel` 3 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (only the known output-filename collision warning) |
| Schemas against `crn-craft describe KIND` | 26 byte-identical; `check-1.0` regenerated in `d15dfa2`; `analysis`, `analysis-v1` are no kind since M5b |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass, 134 files |
| `check_dependencies.py` | pass |
| `check_contracts.py` (rewritten, M7) | 126 of 126 |
| `check_network.py`, `check_structure.py` | 219 of 219, 309 of 309 |
| `reference_rng.py` | pass |
| `check_ssa.py`, `check_cme.py` (and `--command simulate-cached`, equal reports) | 17 of 17, 14 of 14 |
| `check_export.py` | 9 of 9 |
| `check_fsp.py`, `check_fsp_rules.py` (DSMTS) | 16 and 16 rows, largest distribution error 4.49e-14 |
| `check_tau.py`, `check_tau_rules.py` | pass, 30,000 paths each |
| `check_rules.py` | 47 of 47 |
| `check_mtdna.py` | 151 checks, 32,768 paths |
| `check_python.py` (rewritten) | 17 rows, `to_gillespy2` built |
| `legacy/check_python_network.py` (the 0.6.0 script on `_legacy`, new in CI) | 12 rows |
| `check_rules_python.py` | 27 checks |
| `check_compiled_loop.py --require-cargo` | 57 rows, every one identical to the interpreter |
| I2 `check_manual.py --require-cargo` | 22 of 22 blocks |
| `check_builder.py` | 18 tests, 0 failed |
| `check_hybrid.py` | 11 checks, 0 failed |
| `check_alternative_contracts.py` (Python CLI `ssa` hidden, still working) | 20 checks |
| `check_alternative_ssa.py`, `--fast` | 210 of 210, 210 of 210 |
| `check_importers.py` | 49 rows |
| `check_tutorial.py` | 17 blocks |
| `check_installed.py` (`python -I`, from `/tmp`) | 18 checks (new: `unified_api`) |
| On `83369a0`: build, tests, clippy, rustdoc, wheel, schemas | as above: 532 tests, 26 schemas identical |
| On `83369a0`: `check_contracts.py`, `check_python.py`, `legacy/check_python_network.py`, `check_rules_python.py`, `check_alternative_contracts.py`, `check_tutorial.py`, `check_installed.py` | 126 of 126, 17 rows, 12 rows, 27, 20, 17 blocks, 18; `-O` identical |

Not run: `mutation_check.py` (it rebuilds a copy of the tree); no mutated
file changed, and all 47 anchors were checked to match exactly once. The
optional workspaces (`integrations/*`, the parameter worker) were not
rebuilt: they use neither `crn_craft::describe` nor `RulesCheck`
(`grep`); `validation/merge/kinetics-reference` reads
`crn_craft::describe()["version"]`, a field `describe/1.0` keeps.

## Limitations

- Open decision: `rules.ensemble(..., quantiles=...)` (M4e) and
  `crn.ensemble` build a new `run/1.0` Static Direct request without a
  selector, so it runs `auto`; on a model of 512 or more lowered channels
  that is composition–rejection, and the run's numbers then differ from the
  same `rules.ensemble` call without quantiles (the SIR and every validated
  model are below 512 channels, so nothing measured changed). M4e's
  docstring promise "the paths are those of the request without quantiles"
  holds only below 512 channels; the `ensemble` docstring now says so.
  `prepare(..., quantiles=)` of a rules-run request writes
  `"selector": "linear"` (this task), so it keeps its paths. Writing
  `"linear"` in `rules.ensemble` too would keep the promise for every model
  but change the SIR fixture's `run_sha256` and the equality of
  `crn.ensemble` and `rules.ensemble(quantiles=...)` documents; it was not
  done.
- A `bind_parameters` candidate still reports documents rebuilt in Python
  (M6b): its request is bound to the engine's re-serialized bytes, so the
  reported pair replays only after rebinding the hash. A native accessor
  to the engine's bytes was not added.
- The top-level `ensemble` keeps the tree method as its default; a
  higher-order model needs `method="static_direct"` (or `next_reaction`,
  `tau`, `hybrid`).
- `from_dict` of an editable document writes the builder's document after
  the first edit, without the explicit defaults of the original (same
  meaning, other bytes and hash).
- `to_gillespy2` was tested with GillesPy2 1.8.3 installed in the task
  venv (model construction); CI's `python` job has no GillesPy2 and tests
  the `ImportError` path. Running the returned model needs GillesPy2's own
  toolchain.
- `read_sbml`/`read_bngl` keep the importers' explicit convention names;
  the plan's example `convention="binomial"` is a BioNetGen convention.
- The docs still show the removed names; M10 rewrites them:
  `README.md`, `python/README.md`, `docs/tutorial.md` (text; its code block
  now imports `crn_craft._legacy`), `docs/alternative-ssa.md`,
  `docs/architecture.md`, `docs/direct-ssa.md`, `docs/reusable-parameters.md`,
  `docs/comparison-gillespy2.md`,
  `docs/comparison-catalyst-jumpprocesses-before-acceleration.md`,
  `validation/parameter-optimization/README.md`. The manual is frozen (I2).

## Deviations from the plan

- Item 5 ("remove the functions ... and the PyO3 functions behind them") is
  replaced, as instructed, by moving them unchanged to the private
  `crn_craft._legacy`; the PyO3 functions stay until M8b.
- `crn.Model` is a subclass of `crn_craft.rules.Model` (class flag), not
  the same object: `crn_craft.rules.Model` must keep writing the 0.6.0
  versions. The "one flag" is `_unified` on the model class or instance
  and on `PreparedEnsemble`.
- Extra public names: none at the top level; in `crn_craft.rules`:
  `FspResult`, `DEFAULT_QUANTILES`, `SELECTORS`.
- Through the top-level names a bound parameter value must be a real
  number (booleans refused); `crn_craft.rules` keeps converting them as
  0.6.0 did.
- `rules.describe()["crate"]` is now `"crn-craft"` (M7), which changes that
  one field of a `crn_craft.rules` result.
- A network model through `crn_craft.rules` is refused with `CRN-E033`
  (M7) instead of `CRN-E003` at its first network field.

## Migration notes (Python)

| 0.6.0 top level | 0.7.0 |
| --- | --- |
| `crn.load_example(name)` (network bytes) | `crn.load_example(name)` returns a `crn.Model` of the converted model; `crn.convert(bytes)` converts any network model |
| `crn.run_spec(...)` + `crn.simulate(model, run)` | `crn.ensemble(model, end=..., report_times=..., replicates=..., method="static_direct").simulate()` |
| `crn.prepare(model, run, cached=...)` (network) | `crn.prepare(model, request)` takes a rules or `run/1.0` request |
| `crn.ssa_spec(...)`, `crn.simulate_ssa`, `crn.prepare_ssa` | `crn.ensemble(..., method="static_direct", selector="composition_rejection" / "sorting" / "rssa_cr")`, or `method="next_reaction"` |
| `crn.prepare_parameters(model, request)`, `ParameterTemplate.bind(values)` | `crn.ensemble(...)` (or `crn.prepare`) then `run.bind_parameters({name: value})` |
| `crn.mass_action_kernel_source(model)` | `crn.loop_source(model)` and `crn.compile_loop(model)` |
| `crn.fsp(model, run)` | `model.fsp(end=..., report_times=..., bounds={...})` (an `FspResult`); `crn_craft.rules.fsp` returns the document |
| `crn.tau(model, request)` | `crn.ensemble(..., method={"kind": "tau", "epsilon": ...})` |
| `crn.analyze(model)` | `model.analyze()` |
| `crn.check(model, run)` (check/0.1) | `crn.check(model, request)` (check/1.0, or rules-check/0.1 for a rules model alone) |
| `crn.describe()` (describe/0.11) | `crn.describe()` (describe/1.0) and `crn.describe(kind)` |
| `crn.model_bytes(document)` | `crn_craft.rules.model_bytes` (a 0.6.0 attribute of `crn_craft.rules`) |
| `crn_craft.importers.import_sbml/import_bng_net` | unchanged; `crn.read_sbml(path, convention=...)` and `crn.read_bngl(path, convention=...)` return `Model` objects |
| everything above | still available, unchanged, as `crn_craft._legacy.<name>` until M8b (private) |

`crn_craft.rules` code of 0.6.0 runs unchanged with the same documents and
results. New there: `selector=`, `Model.from_dict`, `Model.fsp`,
`Model.analyze`, `Model.to_gillespy2`, three or more reactant molecules in
`add_reaction`.

## For M8a: files switched to `crn_craft._legacy`

Mechanically (`import crn_craft as crn` became `import crn_craft._legacy as
crn`, nothing else), in commit `60d1eba`:

- `examples/circuit/screen.py`, `examples/tutorial/make_assets.py`
- `tutorials/optimization/calibration.py`, `oscillators.py`, `sensitivity.py`
- `docs/tutorial.md` (the Python block's import; run by `check_tutorial.py`)
- `validation/check_alternative_contracts.py`, `check_alternative_ssa.py`,
  `check_jobs.py`, `check_parameter_binding.py`,
  `check_parameter_interface.py`, `compare_gillespy2.py`, `compare_julia.py`,
  `compare_julia_accelerated.py`, `compare_julia_gap10.py`,
  `measure_alternative_ssa.py`, `profile_ssa_acceleration.py`
- `validation/check_importers.py` (`crn_craft._legacy.analyze`)
- `validation/check_installed.py` (`legacy.<name>` for `load_example`,
  `run_spec`, `simulate`, `describe`, `model_bytes`, `analyze`, `fsp`,
  `tau`)
- `validation/legacy/check_python_network.py` (the old `check_python.py`)

`validation/check_rules_python.py` uses `rules.model_bytes` instead of
`crn_craft.model_bytes`. Not switched, because they run on the 0.6.0
package: `validation/merge/baseline.py`, `validation/merge/references.py`.
`validation/gap10_python_reference.py`, `validation/merge/check_laws.py`
and `check_higher_order.py` use only names that stayed.
