# Changelog

All notable changes are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) with breaking contract changes
called out. Every contract is versioned on its own (`crn-craft.<kind>/<n>`).

## [Unreleased]

### Changed

- README: links the [overview of the 18 related repositories](https://github.com/dietmarwo/fcmaes-rust/blob/main/docs/repositories.md).

### Fixed

- Documentation that described 0.3 or pre-0.6 behaviour as current:
  - The crate READMEs' `Cargo.toml` snippets now use `"0.4"`; `"0.3"` resolved to 0.3.x, which lacks rules-run/0.6 and `compile_loop`.
  - The performance analysis no longer calls `rules.compile_loop` unreleased.
  - The kernel checklist labels its 69/135 ns and 1.12× figures as measured with the pre-0.6 default (inverse waits, per-event trackers).
  - The technical report gives the native discovery version per release.
  - The user manual shows `add_continuous`'s bound arguments as keyword-only.
  - The tutorial notes that next reaction defaults to `fast=False`.

## [0.4.0] - 2026-11-06

Faster generic Static Direct for Python: new defaults in rules-run/0.6
(Ziggurat waits and sparse trackers, 1.17× per event on the three-class
mtDNA model) and an opt-in compiled event loop (`rules.compile_loop`,
byte-identical to the interpreter, 2.2× on that model). Together 7.0 s →
2.6 s on the 1000-cell case of the performance analysis. Python
`static_direct` results use a different random stream by default; see
Migration.

Static Direct gets faster defaults in a new request version. Stored
rules-model 0.1–0.3 files and rules-run 0.1–0.5 requests keep their meaning
and replay byte-identically; Python's `method="static_direct"` now samples
with a different use of the random stream (see **Migration**).

### Added

- `crn-craft.rules-run/0.6`: the fields of 0.5, with new `static_direct`
  defaults. An absent `waiting_time` selects Ziggurat waits, and an absent
  `sparse_trackers` integrates every eligible tracker (an unweighted
  single-class count without a first-passage predicate) over its constant
  spans; without an eligible tracker the run integrates per event and does
  not warn. Explicit `waiting_time: "inverse"` and `sparse_trackers: false`
  restore the 0.5 default; an explicit `sparse_trackers: true` without an
  eligible tracker still warns with `CRN-R051`. The result is
  rules-result/0.3, whose shape is unchanged; its `method` echoes the
  request, and `sampler` names what ran. Every other method, including
  `hybrid` and the Static Direct fallback of `mtdna_direct` (inverse waits,
  per-event trackers), runs as under 0.5. The rules-run schema gains 0.6.
- On the 1000-cell three-class case of `examples/mtdna/performance.py`
  (two threads, two count trackers), the 0.6 default took 41.6 ns per event
  against 48.8 for the 0.5 default, 1.17× faster against a declared 1.10×
  target (`validation/mtdna/static-direct-defaults-timing.json`). It passed
  the K = 20 master-equation gate of the Ziggurat sampler with and without
  eligible trackers, and the worker-replay and event-log integral checks
  (`validation/mtdna/gates-static-direct-defaults.json`,
  [log](docs/campaign/static-direct-defaults-log-2026-11-06.md)).
- Discovery: the rules entry of the native CLI's `describe` gains
  `accepted_requests` (rules-run/0.1 to 0.6) and
  `static_direct_defaults_request`, so the native CLI reports
  `crn-craft.describe/0.13`; core and Python `describe/0.10` are unchanged.
- Tests: stored requests of every version from 0.1 to 0.5 replay to recorded
  result hashes (`crates/crn-craft-rules/tests/stored`); the generated event
  loop and kernel replay the 0.6 default byte for byte.
- `validation/check_mtdna.py --section generic` runs only the Ziggurat and
  sparse-tracker gate; `validation/benchmark_static_direct_defaults.py`
  times the two defaults.
- Opt-in compiled Static Direct event loops from Python, off by default.
  `crn_craft.rules.compile_loop(model, *, cache_dir=None, cargo="cargo",
  source=None, progress=True)` writes a small Cargo project around the
  generated rules loop, builds it with the user's Rust toolchain (1.85 or
  newer) against `crn-craft-rules` at exactly the installed version, from
  crates.io or from a checkout (`source=` or `CRN_CRAFT_SOURCE`), and
  returns a `CompiledLoop` (`fingerprint`, `path`, `version`, `cached`,
  `build_seconds`, `log`). `rules.ensemble(..., compiled=loop)`,
  `rules.prepare(..., compiled=loop)` and `simulate(compiled=loop)` run
  Static Direct requests in that program as a subprocess, with results
  byte-identical to the interpreter's; `CompiledLoop.verify(run)` checks a
  request both ways. Builds are cached in the user cache directory
  (`cache_dir=`, `CRN_CRAFT_CACHE_DIR`), keyed by the loop fingerprint, the
  installed version, `rustc -vV`, the target and Rust flags, and share one
  Cargo target directory; every call checks the program against the
  interpreter before returning it. `rules.loop_source(model)` returns the
  generated source without cargo. Refusals: `CRN-R020` at `/method` for
  other methods and for models the generator refuses, `CRN-R020` at
  `/model` for other model bytes and bound candidates, `ValueError` for
  `callback` and `batch_callback`, and `rules.CompiledLoopError` (reason
  `cargo_missing`, `build_failed` with the log path, or `mismatch`). It
  serves the rules-run/0.6 defaults and every explicit option. On the
  1000-cell three-class workload with two threads, end to end from Python,
  it was 2.21× faster than the interpreter under the 0.6 default (Ziggurat
  waits, sparse trackers) and 1.86× with inverse waits and per-event
  trackers (declared target 1.5×); the first build took 42 s on the two-CPU
  host. See
  `docs/decisions/compiled-rules-loop-python-2026-11-06.md` and
  `validation/mtdna/compiled-loop-python.json`.
- `validation/check_compiled_loop.py`, run in the `python` CI job against
  the checkout: byte-identical results over 1 and 2 threads, both waiting
  times, both tracker policies, event logs, the takeover stop, budgets and
  a cancel file, and every refusal. `validation/benchmark_compiled_loop.py`
  times the route. `validation/check_manual.py` runs manual blocks that
  call `compile_loop` against the checkout (`CRN_CRAFT_SOURCE`) and skips
  them without cargo unless `--require-cargo`.

### Changed

- The Python extension encodes rules results with
  `bindings/python/src/rules_payload.rs`, which compiled loops share; the
  bytes are unchanged.
- CI: `retry-runner-shortage.yml` reruns the failed jobs of a release or CI
  run when every job that did not succeed was cancelled before a runner picked
  it up (hosted macOS arm64 shortage), up to three attempts. A release is no
  longer held back until someone reruns it by hand.
- Python: `rules.ensemble` writes rules-run/0.6 for `static_direct`, given
  as a name or as a dictionary; options a dictionary leaves out take the 0.6
  defaults. Results carry `metadata["method"]` (rules-result/0.3).
  `rules.prepare` keeps the version of the request it is given.
- Rust: `RunVersion` gains `V06`; code that matches it exhaustively needs a
  new arm. A `Method::StaticDirect { waiting_time: None, sparse_trackers:
  None }` means the defaults of its request's version.
- The mtDNA user manual, the performance analysis, the rules guides and the
  READMEs describe the new defaults; the manual examples that use
  `static_direct` show their new output.

### Migration

- **Python `static_direct` results change.** The same seed gives other
  cell-by-cell numbers, with the same distribution: sampler
  `chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1` with an
  eligible tracker, `chacha8-u53-ziggurat/rules-static-direct/1` without.
  Tracker integrals can also differ in their last bits. To reproduce the
  numbers of 0.3.0, pass
  `method={"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False}`;
  the sampler is then `chacha8-u53/rules-static-direct/1` again. The result
  document is rules-result/0.3 instead of 0.2 and adds `method`; frames,
  paths and observations equal those of 0.3.0.
- **JSON and Rust requests** of rules-run/0.1 to 0.5 are unchanged. A
  request moved to 0.6 without options gets the new defaults; add
  `"waiting_time": "inverse", "sparse_trackers": false` to keep the old
  numbers.
- **`mtdna_direct` fallback.** A model `mtdna_direct` does not certify still
  runs Static Direct with inverse waits and per-event trackers. Its numbers
  equal those of the explicit options above, no longer those of Python's
  `method="static_direct"`.

## [0.3.0] - 2026-11-06

The answer to the first user feedback on the 0.2.0 Python interface, from a
researcher simulating mtDNA mutant models: six findings, fixed in the
increments I1–I7 of [the plan](docs/plan-0.3.0-feedback.md). Python results
become pandas DataFrames, which breaks code written for the NumPy results of
0.2.0; see **Migration** below. Stored rules-model 0.1/0.2 files and
rules-run 0.1–0.4 requests keep their meaning and replay byte-identically.

### Highlights

- **Continuous variables** (rules-model/0.3, I5): ODEs between events that
  rule rates read, with optional lower bounds that stop or clamp, simulated
  by the new approximate `hybrid` method at declared tolerances
  (rules-run/0.5, rules-result/0.3). The ATP model is now an ordinary model
  file, `models/rules/mtdna-three-class-atp-continuous.rules.json`; the
  legacy `atp` block and exact `thinning` keep working (finding 5).
- **`crn_craft.rules.Model`** (I4): build models from Python with
  `add_reaction(rate, reactants, products)`, as in rebop; a number is a
  mass-action constant, a string the whole propensity, backed by the new
  whole-rule `propensity` of rules-model/0.3 (finding 6).
- **pandas results** (I3): records, rule counts, statistics, paths, trackers
  and event logs are DataFrames; `keep_paths="all"` keeps every replicate and
  `Result.to_xarray()` returns the paths as a Dataset (finding 2).
- **Event logs fit the default budget** (I2): the reported 50-entry mutation
  log for 400 cells of the three-class model is admitted under the Python
  default of 64 MiB, and the `CRN-R020` message names the bytes needed and
  the fix (finding 3).
- **A runnable manual** (I1, I7): every Python block of the mtDNA user manual
  runs on its own and is checked against its printed output in CI; rule
  names are on prepared runs and results (finding 4).
- **Method guidance** (I6): the docs lead with the generic methods;
  `mtdna_direct` is documented as an opt-in fast path with its requirements
  and fallback (finding 1).

### Added

- rules-model/0.3 (I4, I5): a rule gives exactly one of `rate` (per reactant
  tuple) and `propensity` (its whole rate over parameters and observables;
  `CRN-R013` for both or neither, `CRN-R015` before 0.3). A propensity
  applies while each reactant pattern can take a distinct molecule, picks
  molecules uniformly, and fails the replicate with `CRN-R030` at
  `/rules/i/propensity` when negative or not finite. Every method samples
  propensities; the static methods need each reactant type to lower to one
  tuple (`CRN-R020` at `/method` otherwise), `mtdna_direct` falls back to
  Static Direct with `CRN-R051`, and the generated Rust kernels and loops
  refuse them.
- rules-model/0.3 `continuous` variables `{name, initial, derivative,
  lower_bound?, on_bound?}` with `on_bound` `"stop"` or `"clamp"` (I5).
  Rates and propensities may read them; the stop, observables, patterns,
  products and trackers may not (`CRN-R012`, `CRN-R020`). Results list them
  as observable columns after the model's observables.
- Method `hybrid` (rules-run/0.5, `{"kind": "hybrid", "rtol", "atol"}`, I5):
  integrated propensities (Davis 1984; Alfonsi et al. 2005) with the
  Dormand–Prince 5(4) pair, Hairer–Nørsett–Wanner step control and dense
  output; event times, reports and bound crossings located by the Illinois
  method. Sampler `chacha8-u53/rules-approximate-hybrid-dopri54/1`, listed
  as approximate. Its result is `crn-craft.rules-result/0.3`: rules-result/0.2
  plus `method` (with the tolerances) and `continuous` (the column names).
  Only `hybrid` runs continuous variables (`CRN-R020` at `/method` for every
  other method); it refuses ATP models and models that do not lower.
- `models/rules/mtdna-three-class-atp-continuous.rules.json`, the ATP model
  with the continuous variable `ATP` (I5).
- Python `crn_craft.rules.Model` with `add_parameter`, `add_species`,
  `set_initial`, `add_reaction(rate, reactants, products, name=None)`,
  `add_continuous`, `to_dict` and `to_json`; `rules.ensemble` and
  `rules.prepare` accept a `Model` (I4, I5). `method="hybrid"` writes
  `rtol=1e-6` and `atol=1e-9` unless given and selects rules-run/0.5.
- `PreparedEnsemble.rule_names` and `Result.rule_names` (I1);
  `Result.rule_events`, `Result.paths`, `Result.to_xarray()` and
  `event_columns()` without a replicate for every logged event (I3).
- Validation: `validation/check_manual.py` (every manual block, I1),
  `validation/check_builder.py` (18 predeclared z tests against exact means
  and rebop 0.10.0, all passed; [results](validation/builder/RESULTS.md), I4)
  and `validation/check_hybrid.py` (five DKW gates against exact laws and four
  predeclared two-sample tests against exact thinning on the ATP model, all
  passed at 400 and 4,000 replicates per arm;
  [results](validation/hybrid/RESULTS.md), I5). CI runs all three in the
  Python job, normally and under `-O`, and installs `rebop==0.10.0` there.
  `validation/check_installed.py` checks the builder, the pandas results,
  the reported event-log request at the default budget and a small `hybrid`
  run from the installed wheel (I7).
- Examples `examples/mtdna/04_callback.py`, `05_structure.py`,
  `12_builder.py` and `13_continuous_atp.py`; manual chapters on
  continuous variables (10) and on building models in Python (12) (I1, I4,
  I5, I6, I7).

### Changed

- Python results are pandas objects (I3); pandas (`>=2.0,<4`) is a required
  dependency and xarray an optional extra (`crn-craft[xarray]`). See
  **Migration**.
- The native CLI's discovery is `crn-craft.describe/0.12` (I7): its rules
  entry (`dynamic_species`) gains `hybrid_request` and `hybrid_result`, lists
  `hybrid` among the methods and approximate methods and
  `crn-craft.rules-model/0.3` among the accepted models. The core and Python
  `describe()` stay at `describe/0.10`.
- The `CRN-R020` reservation message itemizes replicates × bytes per
  replicate (log entries × bytes per entry, plus trackers), the total and
  `max_output_bytes`; its hint names the smallest `max_output_bytes` that
  admits the request (I2). The reservation and admission checks run after
  lowering, so only the order of diagnostics for requests with several
  errors can change.
- The manual, `docs/mtdna-rules.md`, `describe` and the `rules.ensemble`
  docstring lead with the generic methods; `mtdna_direct` and the generated
  Rust loops are described as opt-in, with what each requires and gains,
  and the performance analysis says which numbers depend on the fast path
  (I6). The manual's chapters 12–18 are now 13–19.
- The rules-model, rules-run and rules-result schemas gain the 0.3, 0.5 and
  0.3 shapes; the descriptions of other fields change in places (I2, I4, I5).

### Fixed

- An event log of 50 mutation events for the three-class model was refused
  with `CRN-R020: observation reservation exceeds budget` above 263 cells
  under the Python default budget (I2, finding 3). Each log entry is now
  bounded by the model's own sizes (the three-class model reserves 721
  bytes per entry instead of 5,088), never above the 0.2.0 bound, so no
  previously admitted request is refused. The 64 MiB default stays.
- The manual's chapter 7 callback used `names` without defining it, and the
  chapter 5 snippet used `result` from an earlier block (I1, finding 4).
  Every block now runs on its own.

### Migration

Python code written for the NumPy results of 0.2.0:

| 0.2.0 | 0.3.0 |
| --- | --- |
| `result.records` was a dict of NumPy arrays | a DataFrame, one row per replicate, indexed by replicate id (`records.index`, name `replicate`); `records["id"]` is gone |
| `records["rule_events"]` (a replicates × rules matrix) | `result.rule_events`, a DataFrame with one column per rule name |
| `statistics["mean"][:, names.index("N")]` | statistics are DataFrames indexed by report time with one column per observable: select by name, `statistics["mean"]["N"]` |
| `event_columns(id)` returned arrays `time`, `rule`, `observables` | a DataFrame indexed by `entry` with `time`, `rule`, `rule_name` and one column per observable; `event_columns()` gives every log, indexed by (`replicate`, `entry`); the old matrix is `event_columns(id)[list(metadata["observables"])].to_numpy()` |
| `tracker_columns(i)` returned arrays with `id` | a DataFrame indexed by replicate id |
| `keep_paths=True` kept one replicate | `keep_paths=True` (or `"all"`) keeps every replicate; write `keep_paths=1` for one |
| `metadata["paths"]` rows | `result.paths`, indexed by (`replicate`, `time`) |
| NumPy arrays | `.to_numpy()` on any frame or column; integer termination codes are `records["termination"].cat.codes` |

- **Silent break:** `termination` is categorical with names, so
  `records["termination"] == rules.TERMINATIONS.index("stopped")` is `False`
  for every replicate under pandas, with no error. Compare with the name,
  `records["termination"] == "stopped"`.
- **New dependency:** pandas is installed with the package. Environments
  that install the wheel with `--no-deps` must install pandas themselves.
- **Contracts:** `crn-craft.rules-model/0.3`, `crn-craft.rules-run/0.5` and
  `crn-craft.rules-result/0.3` are new. Stored rules-model 0.1/0.2 files and
  rules-run 0.1–0.4 requests keep their meaning and replay byte-identically;
  0.1–0.4 requests still return rules-result/0.1 or 0.2. A model cannot
  have both `atp` and `continuous` (`CRN-R015`).
- **Schema:** `rate` is no longer a required property of a rule in the
  rules-model schema, since a 0.3 rule may give `propensity` instead. The
  engine still rejects a 0.1 or 0.2 rule without `rate` (`CRN-E003`), but a
  tool that validates such files only against the schema no longer catches
  it.
- **Discovery:** clients of the native CLI's `describe` must accept
  `crn-craft.describe/0.12` and the new rules entries; core and Python
  discovery are unchanged at 0.10.
- **Rust (`crn-craft-rules`):** `Method` gains `Hybrid { rtol, atol }`,
  `RunVersion` gains `V05` and `ResultVersion` `V03`, so exhaustive matches
  need the new arms. `RulesResult` gains `method` and `continuous`, the
  model gains `continuous`, and a rule's `rate` is now `Option<String>`
  next to the new `propensity: Option<String>`.

## [0.2.0] - 2026-11-04

The first release: the crates `crn-craft`, `crn-craft-rules` and
`crn-craft-cli` on crates.io and the Python package `crn-craft` on PyPI, with
abi3 wheels for CPython 3.10+ on Linux x86_64, Windows x64 and macOS
x86_64/arm64 and a source distribution. The Python bindings, the diffsol
integration and the job worker stay source-only workspaces. The development
iterations that make up this release follow.

### Highlights

- Exact SSA: legacy checked Direct and dependency-cached Direct archives, and
  `direct`, `cached_direct`, `direct_cr`, `sorting_direct` and `rssa_cr` with
  inverse or Ziggurat waits (`crn-craft.ssa-run/0.3`), where `fast`, emulating
  JumpProcesses.jl, is the default and `fast: false` the checked policy, plus
  Gibson–Bruck next reaction (`crn-craft.ssa-run/0.4`, checked policy).
- Sparse finite state projection, explicitly approximate Poisson leaping,
  network analysis with conservation laws and cone enumeration.
- Experimental rules engine for dynamic, attributed species (native CLI and
  Python `crn_craft.rules`), with exact-law gates for mtDNA and aggregation
  models: parallel ensembles, event trackers and callbacks, Static Direct,
  next reaction, a specialized fixed-population mtDNA sampler, approximate tau
  and binomial steps, ATP-coupled thinning and generated Rust kernels.
- Python bindings with prepared simulations, reusable parameter templates,
  trajectory diagnostics and strict SBML/BioNetGen import subsets; optional
  diffsol BDF with forward sensitivities; fcmaes job applications.
- Optimization tutorials C1–C4 and measured comparisons with GillesPy2 and
  Catalyst/JumpProcesses.
- Packaging: `python-release.yml` builds and tests every wheel and the source
  distribution with the installed-package acceptance and publishes a `v*` tag
  to PyPI through Trusted Publishing; `publish-crates.yml` verifies the three
  crates from their own archives and publishes them in dependency order.
  `scripts/package_version.py` and `scripts/check_release.py` check versions,
  crate metadata and lockfiles. The manual `wheels.yml` workflow is replaced.
- Generated Rust event loops for Static Direct rules models (Rust only), and
  user guides for mtDNA mutant accumulation: a user manual, an optimization
  tutorial with fcmaes and a performance analysis, with runnable examples.
- Archives record `cargo_lock_sha256` from a build-script hash of the lockfile
  (the packaged crate's own, or the checkout's root one) instead of
  `include_bytes!` of a path outside the crate, which no packaged crate could
  build. Checkout builds record the same hash as before; see
  `docs/archive-contract.md`.

### Generated rules loops, API documentation and user guides

#### API documentation

- Every public item of `crn-craft` and `crn-craft-rules` is documented, with
  crate landing pages, module docs and doc examples for the main entry
  points; both crates warn on missing docs. Generated JSON schemas gain the
  same text as `description` entries; their shapes are unchanged.

#### Release fixes

- The PyPI description says correctly that the rules engine ships in the
  wheel as `crn_craft.rules`; its command-line tools are in the native CLI.
- The mtDNA manual explains that the wheel and `crn-craft-cli` both install a
  `crn-craft` command, and that only the native one has the rules commands.
- `validation/check_installed.py` checks the installed rules engine against
  the immigration-death law and across thread counts.
- The wheel and crate packaging workflows also run on pull requests that
  change crate sources.

#### Added

- mtDNA documentation:
  - `docs/mtdna-user-manual.md`: every application scenario, from takeover
    statistics and trackers to event logs, callbacks, the clone and ATP
    models, sweeps, the CLI and Rust.
  - `docs/mtdna-optimization-tutorial.md`: fitting `mu` and `sm` to
    heteroplasmy distributions with fcmaes, covering noise, common random
    numbers, staged refinement, held-out checks and parallel evaluation.
  - `docs/mtdna-performance-analysis.md`: methods, threads, copy number,
    observation, interfaces, other models, comparisons and study sizing.
  - Runnable scripts in `examples/mtdna/`, with timings in
    `validation/mtdna/user-scenarios-timing.json`.

- Generated rules event loops (Rust only): `generate_rules_loop` emits the
  whole Static Direct event for a lowered model with constant or vacancy
  rates and sum projections, and `execute_with_rules_loop` runs it through
  the new `RulesLoop` trait and `LoopHost`. Results are byte-identical to
  `execute` for every Static Direct request on the same model bytes. On the
  128-cell mtDNA probe it was 2.3× faster than the interpreter under
  Ziggurat waits and sparse trackers and within 1.36× of `mtdna_direct`.
  New example `generate_rules_loop`, tests with ten checked-in loops and
  `validation/check_generated_loop.py`. See
  `docs/decisions/generated-rules-loop-2026-11-03.md`.

#### Changed

- Comparison documents state that they were measured on development builds
  of 2 November, before the 0.2.0 release, and name the engine changes since.
  Wording that predates the release is corrected: `fast` is the default of
  `ssa-run/0.3`, the legacy `simulate` commands are labeled checked rather
  than default, and the status lines no longer say nothing was published.
  The GillesPy2 figure's legend is regenerated from the same report.
- Generated rules loops emit every stop predicate, not only the guarded
  fraction form, with the bytecode's semantics. An undefined value still
  reaches the interpreter's diagnostic. With the stop
  `N > 0 && mutants / N >= 0.9`, the mtDNA probe's loop went from 90 to
  53 ns per event. Existing loops must be regenerated, because the emitted
  source, and so the fingerprint, changes.
- Static Direct finds the smallest positive rate while summing in networks
  of at least 32 channels instead of testing every rate in a second pass.
  Results are byte-identical; the 96- and 512-channel probe networks ran
  1.25× and 1.34× faster. `validation/probe_structural_rules.py` gains the
  512-channel workload.

#### Validation

- A compact, table-driven loop layout for densely coupled networks above
  1 MiB replayed byte for byte but reached only 1.15× over the interpreter
  on a 256-channel network (criterion 1.25×), and was not adopted. Fusing the
  minimum into unrolled loops slowed them down. Both are recorded in
  `validation/mtdna/generated-loop-followups.json`, with the compact patch
  under `validation/mtdna/followup-attempts/`.
- Recorded a loop of structural changes meant to close the gap between
  generic Static Direct and `mtdna_direct`: inline rate factors, a fused
  total and resolution check, exact incremental sum projections and a
  grouped evaluator. None met its declared criterion, so samplers,
  identities and contracts are unchanged. The probe is
  `validation/probe_structural_rules.py`; the results and attempt patches
  are under `validation/mtdna/`. See
  `docs/decisions/generic-rules-structural-2026-11-03.md`.

### mtDNA ensembles and faster generic rules

#### Changed: documentation for both engines

- `docs/architecture.md` rewritten for the fixed-network and rules engines; discovery versions, method lists, diagnostics (`CRN-E067`, `CRN-R031`, `CRN-R051`) and status lines updated across README, AGENTS.md and docs/.
- Fixed: `rules.ensemble` with trackers but no method raised `KeyError`.

#### Changed: faster generic rules simulation (S1–S6)

- Static Direct, next reaction, ATP thinning and the leap methods cache rule rates and shared vacancy bases, recompute only the channels and observables an event affects, defer report-only observables, keep the live-instance count incrementally and evaluate the stop and predicates from compiled bytecode. Existing requests return byte-identical results.
- `rules-run/0.4` adds two opt-ins on `static_direct`: `waiting_time: "ziggurat"` and `sparse_trackers: true`. Each has its own sampler identity; the result stays `rules-result/0.2`. If no tracker is eligible, `sparse_trackers` reports the default identity and `CRN-R051`.
- `generate_rules_kernel` emits a Rust kernel for one closed model with constant or vacancy rates. `execute_with_rules_kernel` runs it under the Static Direct scheduler with byte-identical results. It refuses kernels generated for other model bytes or by another generator version, because the fingerprint covers the emitted source. Rust only; it is not available from Python.
- [docs/specialized-kernels.md](docs/specialized-kernels.md): when a model-specific kernel is worth it, and the review checklist. Review log: [docs/campaign/fast-generic-rules-review-log-2026-11-03.md](docs/campaign/fast-generic-rules-review-log-2026-11-03.md).

#### Added: mtDNA ensembles, event observation and alternative rules methods (A0–A6)

- Three-class mtDNA fixtures (`models/rules/mtdna-three-class*.rules.json`) with independent CME gates, and a usage guide in [docs/mtdna-rules.md](docs/mtdna-rules.md).
- `crn_craft.rules`: prepared parallel rules ensembles from Python with read-only NumPy record columns, cooperative cancellation, budgets, `bind_parameters` candidates, live and batch event callbacks, trackers (crossings, episodes, extrema, integrals, first passage) and a bounded event log. The Rust `EventObserver` interface exposes the same committed events.
- Rules methods: `static_direct`, `mtdna_direct` (exact grouped Direct for the fixed three-class model, with an explained Static Direct fallback, `CRN-R051`), `next_reaction` (Gibson–Bruck), approximate `tau` and `binomial`, and `thinning` for the ATP-coupled model (`rules-model/0.2` `atp` block).
- Contracts: `rules-run/0.2` (method, trackers, event log) and `0.3` (adds `mtdna_direct`) answer with `rules-result/0.2`; `rules-model/0.2` adds `atp`. `rules-run/0.1`, `rules-model/0.1` and the tree method keep their meaning and result bytes, and reject the new fields even when empty or null. Network `ssa-run/0.4` is 0.3 plus the `next_reaction` aggregator (checked policy and inverse waits only) and answers with `ssa-result/0.4`.
- Discovery: core and Python `describe/0.10`, native CLI `describe/0.11`, with the 0.4 request default and the rules methods.
- Migration: stored requests of earlier versions are unchanged. `crn.ssa_spec(aggregator="next_reaction")` writes 0.4 with `fast=False`; other aggregators still get 0.3. A present `method`, `event_log`, `trackers` or `atp` must hold a value; `null` is rejected.
- Review log: [docs/campaign/mtdna-a0-a6-review-log-2026-11-03.md](docs/campaign/mtdna-a0-a6-review-log-2026-11-03.md).

### Alternative exact SSA methods and the Julia-gap campaign (`0.2.0-dev.1`)

#### Changed: `fast` is the default (`crn-craft.ssa-run/0.3`)

- New request version `crn-craft.ssa-run/0.3` with the same fields as 0.2; `method.fast` defaults to `true`. It answers with `crn-craft.ssa-result/0.3`, which always states `method.fast`. Set `"fast": false` for the checked exact policy (ChaCha8 words, lost-increment and 53-bit resolution checks, stop polls before every random word). `fast` must be a boolean; `null` is rejected.
- `crn.ssa_spec` writes 0.3 requests and now defaults to `fast=True`; pass `fast=False` for the checked policy. Non-boolean values raise `TypeError`.
- Migration: stored `ssa-run/0.1` and `0.2` requests keep their meaning and return the same bytes as before (0.1 rejects `fast`, 0.2 keeps it off unless set). Python callers that relied on the old default get different (equally exact) paths unless they pass `fast=False`; the checked replay evidence scripts now pass it explicitly. With the same seed, a 0.3 request with `fast` on reproduces a 0.2 request with `fast: true`, and a 0.3 request with `fast: false` reproduces a 0.1 request; only the schema strings, `method.fast` and `run_sha256` differ.
- Discovery: core and Python `describe/0.8`, native CLI `describe/0.9`; `alternative_exact_ssa.fast` lists the default per request version.

#### Added: simulator optimization tutorials (contributed by Astra)

- Four lessons in [tutorials/optimization](tutorials/optimization/README.md) with native bindings, independent moment, BDF, generator and finite-state checks, disjoint validation and bounded smoke and study presets. SSA runs use the default `fast` policy; C1's `--include-checked` adds a checked arm and every row records policy, sampler and request schema.
- Review fixes: C3's cross-backend control uses its own seeds instead of the training seed that selected the design, a horizon whose spectral grid contains the 20-unit target period, and per-path oscillatory fractions; C1 records the training loss of the true rates and whether a fit stayed at its truth-derived start; C2 reports local condition numbers and a terminal-product twin instead of a dimension-bound rank, and treats failed native solves as undefined points; C4 reports each kernel's exact-law optimum. Options a chapter does not use are rejected, summaries carry a run state, and `results/verified.json` is written by `record.py` from complete runs. The C3 reference (fcmaes-rust `crn_reference`) takes the initial state from the request and rejects unknown fields.

#### Added: dynamic, attributed species (experimental rules engine, K0/K1)

- New crate `crn-craft-rules` and CLI commands `rules-check` and `rules`: species types with numeric attributes, instances created and removed on the fly, rule schemata over zero to two reactant patterns with product-form propensities, observables (`sum`, `distinct`, `max`, `min`), `new_id()` clones, random product attributes and stop conditions. Contracts `crn-craft.rules-model/0.1`, `rules-check/0.1`, `rules-run/0.1`, `rules-result/0.1`; the native CLI's discovery moves to `describe/0.8` (core and Python stay at 0.7). Reference: `docs/rules-language.md`.
- K1 gates (`validation/check_rules.py`, CI job `rules`): 30 exact-law comparisons over 100,000 replicates per fixture pass (largest distance 0.0034 against a bound of 0.0074), including Kowald's fixed-population mtDNA model and constant and multiplicative aggregation; all six wrong-law power checks reject. After an independent review, identical-pattern pairs use exact pair trees (sampler `chacha8-u53/rules-direct/2`), statistics stream in replicate order with admission limits, and overflows and undefined conditions fail instead of passing silently; the gates were rerun unchanged.

#### Added: reusable parameter bindings and generated kernels (imported, by Astra)

- `parameters::ParameterTemplate` (Rust) and `crn.prepare_parameters` (Python) compile a model's topology once and bind parameter vectors or exact name maps per candidate; each binding keeps byte snapshots for replay. Optional `codegen::mass_action_source` / `crn.mass_action_kernel_source` emit a Rust propensity kernel for `ssa::execute_with_kernel`, guarded by `CompiledNetwork::topology_sha256`. Results are byte-identical to fresh preparation. See `docs/reusable-parameters.md` and the optimization comparison with rebop.
- Review fixes: generated kernels fall back to the new `CompiledNetwork::propensity_in_scheduler_state` instead of the state-scanning `propensity` (a 1,500-reaction review model went from about 40 s to 4.2 s with `direct`); ordinary preparation again checks admission before compiling the dependency graph (`CRN-E060` before `CRN-E040`, as before the import); binding mismatches are `CRN-E014` with missing and unknown names, invalid Python seeds `CRN-E020`; Python `bind` rejects sets and turns float overflow into `CrnError`.

#### Measured: fast mode versus Julia at one and eight threads

- `docs/comparison-julia-8-threads.md`: same-host measurement of the ten `fast` configurations against five JumpProcesses 10.0.1 methods on W1–W7, original and 256-path batches, at one and eight workers (2,240 timed observations; 6,514 terminal-marginal checks passed). Geometric mean of best fast / best Julia: 0.69x (original) and 1.04x (256 paths) at eight workers; 1.15x and 1.12x at one worker (derived on import from the same observations), replacing the earlier serial projection of 1.02x. The harness and raw evidence were not part of the import; all summary statistics were recomputed from the reported observations.

#### Added: the `fast` switch (`crn-craft.ssa-run/0.2`)

- `method.fast` (default `false`) emulates JumpProcesses.jl: Xoshiro256++ words per replicate, running totals without lost-increment or 53-bit resolution checks, and stop polls every 256 event-loop iterations. New methods run 1.09–4.3x faster with it; it passes its own 210-comparison exact-law gate, now also in CI.
- `ssa-run/0.2` and `ssa-result/0.2`; 0.1 requests are still accepted, reject `fast`, and return byte-identical results. Discovery moves to `describe/0.7`. Python: `crn.ssa_spec(..., fast=True)`.
- Benchmarks time CRN-craft with `fast` and report the checked default alongside: the GillesPy2 comparison was rerun (fast cached Direct 1.5–8.2x and the best fast method 2.9–9.2x faster than GillesPy2's C++ solver on W1–W4); the Julia comparison gained a projected `fast` table (geometric mean 1.02x).

#### Imported: ten further Julia-gap iterations (C01–C10, by Astra)

- Kept C05 (lower-bound reuse), C06 (shared exact custom-law ranges) and C09 (compiled elementary mass-action kernels); same-host geometric-mean gap to Julia 1.54x to 1.44x. Evidence in `validation/julia-gap10/`.

#### Added (imported acceleration patch)

- `crn-craft.ssa-run/0.1` and `crn-craft.ssa-result/0.1` with `direct`, `cached_direct`, `direct_cr`, `sorting_direct` and `rssa_cr` aggregators, inverse or Ziggurat (Julia 1.13.3 tables, MIT) waits, and the `ssa` command in both CLIs; discovery `describe/0.6` lists them under `alternative_exact_ssa`.
- Reusable prepared handles in Python: `crn.prepare`, `crn.prepare_ssa`, `crn.ssa_spec`, `crn.simulate_ssa`.
- `docs/comparison-catalyst-jumpprocesses.md` on seven workloads, with the evidence in `validation/ssa-acceleration/`.

#### Performance (declared in docs/decisions/julia-gap-2026-11-02.md)

- Kept: ChaCha8 AVX2 dispatch (A02), exact Hill-rate memo by regulator count (A03, W3 about 2x on every method, legacy included), cached composition envelope (A05), sums-only lazily checked rate tree (A06), and one word for composition member and acceptance (A08; DirectCR and RSSACR sampler version 2, 210/210 exact-law checks).
- Rejected and reverted with evidence in `validation/julia-gap/`: A01 (poll per loop), A04 (redundant RSSACR evaluation), A07 (additions-only ordered sums), A09 (batched tree rebuild), A10 (bucket buffer policy).
- DirectCR, SortingDirect and RSSACR 1.07–2.4x faster than the import, cached_direct 1.01–2.05x; projected gap to JumpProcesses 1.11–1.83x (geometric mean 1.46x, from 2.02x measured). Legacy Direct and cached results are byte-identical.

#### Fixed

- SortingDirect no longer fails a replicate when its ordered scan sum rounds one ulp below the tree total; the threshold is redrawn.
- The alternative-SSA output budget is applied once, to the `ssa-result` JSON, instead of after a legacy archive allowance.
- `python -m crn_craft ssa` exits 1 for rejected or incomplete runs and reports diagnostics like the native CLI; `prepare`/`prepare_ssa` convert non-JSON native errors to `CrnError`.
- The native CLI `--version` milestone matches `describe`; the notices no longer point to a file that does not exist here.

### Campaign review and performance iterations P11–P20 (`0.2.0-dev.1`)

#### Performance (declared in docs/campaign/performance-p11-p20.md)

- P11: Direct selects the channel by binary search over the prefix sums recorded during the strict sum.
- P12: `simulate-cached` resumes the strict sums at the lowest changed channel (1.59x on the sparse 66-reaction network).
- P14: one scheduler poll per event; P16: per-worker buffers reused across replicates.
- P17, P18: FSP skips exact-zero rows, then stores the generator by column and fuses the series accumulation (W8b 2.34x against P10, bitwise identical probabilities).
- P19: tau leaping draws leap means above 16 with Hörmann's PTRS; the algorithm identifier becomes `/0.3` (W9 312x).
- Rejected and reverted with evidence: P13 (skip bitwise-unchanged rates) and P20 (single archive serialization). P15 (flat save buffers) gave no speedup and is kept as the memory-accounting fix.
- Final source: `simulate` 1.4–2.0x faster than 5a07334 on the frozen workloads, `simulate-cached` 2.1–3.3x faster than `simulate` on the 12- and sparse 66-reaction networks; all SSA digests unchanged.

#### Fixed

- A negative count from firing is reported as `CRN-E043`, not as a population-limit stop.
- Saved states use one flat buffer per path, so the published memory reservation is an upper bound (it was exceeded for many save points).
- FSP enforces `max_output_bytes` on the exact serialized size (`CRN-E070`) and sums the omitted series weight explicitly.
- Tau leaping caps leap means at 10^7, uses its own random-stream domain, accepts `--cancel-file`, and its CGP factors are tested.
- `sampled_first_passage` reports only grid-supported bounds; `oscillation` needs 16 samples after burn-in and handles constant series.
- Importers reject BioNetGen fixed species, missing SBML stoichiometry or units, `fast` reactions, identifier collisions and pathological literals within one second (`CRN-E090`).

#### Added

- `docs/comparison-gillespy2.md` and `validation/compare_gillespy2.py`: CRN-craft against GillesPy2 1.8.3 on W1–W4 (`simulate-cached` 1.1–3.3x faster than its C++ solver, 43–126x faster than its NumPy solver), with equivalence checks against each other and exact laws.
- README banner and an index linking every Markdown document.
- Python: Ctrl-C cancels native calls; `cancel_file=` for `simulate`, `fsp` and `tau`; `CrnError` with parsed diagnostics; `__version__`; Python CLI `--version`, `--cancel-file` for `fsp`/`tau` and native-CLI output formatting.
- `docs/tutorial.md` with screenshots and figures generated by `examples/tutorial/make_assets.py` and checked by `validation/check_tutorial.py`; `docs/architecture.md`; `docs/technical-report.md`.
- CI job building the wheel and running the Python, importer, tutorial and installed-package checkers normally and under `-O`.
- Two mutation-gate mutants for the lost-increment checks of the prefix-sum and resumed-sum passes (45 in total).
- `validation/benchmark_iterations.py` (paired iteration timings) and `validation/compare_tau_samplers.py`.

### Iteration 4 / M2 (`0.1.0-dev.4`)

#### Added

- Implement Direct SSA for existing count propensities, exact save-time boundaries and absorbing states.
- Pin ChaCha8 streams and owned uniform/wait/channel mappings; reject unresolved dynamic ranges and clock stagnation.
- Add bounded ordered parallel ensembles, reproducible mean/sample variance and mean SE.
- Enforce event/population/deadline/output/workspace limits and cooperative API/file cancellation; retain every replicate's accounting and withhold partial-ensemble statistics.
- Emit result 0.2 and checksummed archive 0.1; add `simulate`, `verify`, generated schemas and migration notes. Preserve result-v1 and existing inputs/analysis/export.
- Add the predeclared 300,000-path analytic campaign, independent scalar RNG and trajectory oracles, archive-corruption tests and nine new source mutants.
- Preserve all three iteration-3 review regressions and the 17-mutant baseline. No FSP or comparative performance claim.
- Review tests for lost positive increments and underflowed waits at the sampler boundary, retained-trajectory IDs, existing and live cancel files, tampered partial-evidence archives and the archive byte budget, and CLI `simulate`/`verify` exit codes. Eight probe mutants had passed every delivered test.
- `validation/check_cme.py`: four 100,000-path fixtures with second-order, Hill and Michaelis–Menten channels against exact master-equation laws (12 DKW tests, family alpha 0.001, declared in `docs/decisions/iteration-04-review-cme.md`). It runs in CI normally and under `-O`.
- Nine review mutants in the mutation gate, now 35.

#### Fixed

- `--cancel-file` no longer slows simulation about tenfold, and a wall-time limit no longer doubles the cost per event. Every poll used to stat the file, which costs about a microsecond against roughly 100 ns per event for small models, and read the clock. A per-worker poller now reads the cancellation flag at every poll and the wall clock and cancel file at every 64th poll, the file at most once per millisecond. A found file cancels all workers, and a worker keeps an observed stop. Numerical results are unchanged.

### Iteration 3 / M1b (`0.1.0-dev.3`)

#### Added

- Exact nonnegative conservation-cone rays with deterministic work/ray caps, checked i128 arithmetic and explicit unavailable statuses.
- Effective initial-class bounds that resolve the reviewed mixed-sign-basis gap when cone enumeration completes.
- Reaction complexes, linkage/strong/terminal classes, weak reversibility and structural deficiency, with iterative graph walks.
- GillesPy2 1.8.3 elementary mass-action Python export, explicit homodimer conversion, collision-free aliases and unsupported-feature diagnostics.
- Independent support-enumeration/reachability oracle for 300 seeded networks, real target propensity/SSA-smoke checks and seven additional source mutants.
- Tests and mutants for gcd reduction of combined cone rays, a caller ray cap below the number of initial orthant rays, and linkage classes joined only through converging edges, which three probe mutants had passed; the structure oracle gains the `primitive_combination` model.

#### Changed

- Analysis advances to 0.2; `analyze-v1` preserves 0.1. Discovery advances to 0.3; export starts at 0.1. Model/run/check/result/state schemas remain unchanged.
- W001 suppression uses exact cone bounds when available and retains sound basis bounds on unavailable outcomes.
- CI covers structure/oracle checks under normal Python and `-O`, 17 mutants after a passing baseline, and optional pinned reference validation.

#### Limits

- Cone enumeration can grow exponentially; size alone does not ensure completion. No partial cone is labelled complete.
- Export supports mass-action reactant order at most two, no population safety limits and exactly representable initial counts up to 2^53. It is for reference SSA, not target ODE semantics or CRN-craft budget/range equivalence.
- Direct SSA and all CRN-craft simulation remain M2 or later.

### Iteration 2 / M1a (`0.1.0-dev.2`)

#### Added

- Immutable compiled networks, sparse stoichiometry, reactant gates and propensity dependency graphs.
- Checked mass-action, Hill activation/repression and Michaelis–Menten evaluation using pinned pure-Rust libm 0.2.15; atomic i64 reaction updates.
- `analyze MODEL [STATE]`, new state/analysis 0.1 schemas, exact checked small-network conservation with explicit cap/overflow status and initial-class bounds.
- Independent 200-network Fraction/Decimal oracle, source mutation tests, numerical edge/atomicity regressions, documentation and retained evidence.
- Tests for the `log1p` Hill path near x = k, for weight division in conservation bounds and for binomial symmetry at high reactant orders, which three probe mutants had passed.

#### Changed

- Discovery is `crn-craft.describe/0.2` and advertises implemented rate evaluation; model/run/check/result 0.1 contracts are unchanged.
- Analyzer growth warnings are suppressed only by proven initial-class bounds. Structural `check` semantics are unchanged.
- CI includes the independent network oracle normally and under `-O`, source mutation checks and network tests in the LF checkout job.
- The mutation check runs the unmutated suite first and requires it to pass, then runs seven mutants.

#### Limitations

- No scheduler, RNG, SSA/FSP/ODE, linkage classes, expression AST, exporters or ecosystem integrations yet.
- Conservation uses checked i128 within 64 species / 256 reactions; its nonnegative basis heuristic is not a complete bound solver.
- f64 propensity evaluation is tested against declared tolerances, not certified correctly rounded or cross-platform identical.

### Iteration 1 / M0 (`0.1.0-dev.1`)

#### Added

- M0 Rust workspace with strict model/run contracts and reserved result envelope.
- JSON duplicate-key, unknown-field, domain, identity, budget and truncation checks.
- Structured `describe` and `check` CLI commands with explicit unavailable simulation status.
- Generated schemas, complete examples, negative fixtures and regression/acceptance suites.
- Local validation evidence, CI configuration, revised plan and agent/release documentation.
- Tests for unsupported run versions, methods and aggregators, and for both sides of the JSON nesting limit.
- Tests for both sides of the population-limit, FSP-bound and minimum-output-budget rules, which three probe mutants had passed.
- A CI job that checks out the repository with `core.autocrlf=true` and runs the contract suite on it.

#### Changed

- CI cancels superseded runs and runs the dependency gate once, normally and under Python `-O`, comparing the two reports.
- `SOURCE_MANIFEST.json` from the delivery archive is not tracked; git history identifies the files.

#### Fixed

- Decimal literals parse to the correctly rounded double (serde_json `float_roundtrip`). The default parser read about a quarter of random literals one ulp away, including shortest round-trip forms such as `3.983236103458967e-14`.
- Text files keep LF line endings in every checkout (`.gitattributes`), so exact-byte model hashes hold under Git for Windows' `core.autocrlf=true`.
- The nesting-limit hint states the rule the parser applies: depth at most 48, with the root at depth 0.

#### Limitations

- No propensity evaluator, network analyzer, SSA, FSP, ODE, bindings, importers or ecosystem worker yet.
- No simulator correctness, throughput or biological validation is claimed.

### 0.2.0-dev.1 campaign

- Sparse reachable FSP with escaping sink, stable split uniformization and separate error accounting.
- Exact dependency-cached Direct with original ordered streams/results; external cancellation checkpoint.
- Isolated Python ABI3 package, trajectory diagnostics and strict SBML/BioNetGen import subsets.
- Native fcmaes job applications, fixed/disjoint seed accounting and actual shared SDK/runner acceptance.
- Optional diffsol BDF, analytic Jacobians/forward sensitivities and global-to-SQP fitting.
- Experimental non-negative Poisson leaping with explicit approximate request/result schemas.
- Ten performance experiments and same-host reference/quality evidence; no package publishing.

[Unreleased]: https://github.com/dietmarwo/CRN-craft/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/dietmarwo/CRN-craft/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/dietmarwo/CRN-craft/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/dietmarwo/CRN-craft/releases/tag/v0.2.0
