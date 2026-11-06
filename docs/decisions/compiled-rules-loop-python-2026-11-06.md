# Compiled rules event loop from Python

6 November 2026. The target and the gates were declared before the
1000-cell measurement. A 64-cell smoke run during development, with two
threads, took 0.24 s compiled against 0.46 s interpreted.

## Goal

Python users with a Rust toolchain can compile the generated Static Direct
event loop (`generate_rules_loop`, `execute_with_rules_loop`) for their
model and run with it. The route is opt-in, off by default, and its results
are byte-identical to the interpreter's.

## Design

- `rules.compile_loop(model)` asks the extension for the loop source
  (`PreparedRules.loop_source`, which calls `generate_rules_loop`) and
  writes a small Cargo project: the generated loop, the exact model bytes,
  a driver (`bindings/python/src/rules_loop_main.rs`) and the result
  encoder the extension itself compiles
  (`bindings/python/src/rules_payload.rs`). The project depends on
  `crn-craft-rules` at exactly the installed version, from crates.io or
  from a checkout (`source=` or `CRN_CRAFT_SOURCE`).
- The build is a standalone program. Python sends it a rules-run request
  on standard input and reads the metadata and record columns in the form
  the extension returns them, so `Result` parses both alike. No Rust value
  crosses a library boundary, and the program's toolchain need not match
  the one that built the wheel.
- Builds are cached in the user cache directory, keyed by the loop
  fingerprint, the installed version, `rustc -vV`, the target and Rust
  flags, the driver sources and, for a checkout, a digest of its engine
  sources. One shared Cargo target directory serves every model.
- Every `compile_loop` call, cached or not, runs three short requests
  (default, inverse and Ziggurat waits, dense and sparse trackers, an event
  log, kept paths and instances, two threads) through the program and the
  interpreter and returns the handle only if the metadata and every record
  column are byte-identical. `CompiledLoop.verify(run)` does the same for a
  user's request.
- Refusals: models the generator refuses (`CRN-R020` at `/method`, before
  cargo runs), other methods (`CRN-R020` at `/method`), another model or a
  bound candidate (`CRN-R020` at `/model`), callbacks (`ValueError`), and
  `CompiledLoopError` for a missing cargo, a failed build (with the log
  path) or a differing result.

## Gates

- **Identity.** `validation/check_compiled_loop.py` builds the loop of the
  three-class model with the takeover stop against the checkout and
  requires byte-identical results over 1 and 2 threads, inverse and
  Ziggurat waits, dense and sparse trackers with an event log, the takeover
  stop (at least one cell stopped), the default method settings, event,
  sampling, instance, output and wall-time budgets and a cancel file. It
  runs in the `python` CI job, normally and under `-O`.
- **Speed.** Declared target: on the 1000-cell three-class workload of the
  [performance analysis](../mtdna-performance-analysis.md#2-methods)
  (takeover stop, three years, monthly reports, seed 17, two threads, two
  count trackers), `run.simulate(compiled=loop)` is at least **1.5×**
  faster than interpreted `run.simulate()` under `method="static_direct"`.
  Measured end to end from Python, including the subprocess and parsing,
  excluding the build, under the host's CPU lock, as the median of three
  runs in alternating order (`validation/benchmark_compiled_loop.py`).
  Every timed pair must agree, and `verify` must pass. Build times are
  reported, first and cached.
