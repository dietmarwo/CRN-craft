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

## Results

Recorded 6 November 2026 on the two-CPU host, under the CPU lock, with the
loop built against this checkout
(`validation/mtdna/compiled-loop-python.json`).

| Method | Interpreted, median s | Compiled, median s | ns per event | Speedup |
| --- | ---: | ---: | ---: | ---: |
| `static_direct` | 7.06 | 3.71 | 49.2 → 25.8 | **1.90×** |
| Ziggurat waits + sparse trackers | 5.84 | 2.71 | 40.6 → 18.8 | 2.16× |

The declared target (1.5×) is met. All three runs of each route were
within 3.3% of their median, and every pair agreed; `verify` passed for
both methods. 143.4 million events (default) and 143.8 million (Ziggurat)
per run.

Build times, with the crates already downloaded:

- first build into an empty cache and target directory: 39.9 s (cargo
  39.8 s) from the checkout, 38.5 s from crates.io;
- a second model with the shared target directory warm: 2.7 s;
- a cached call, including the self-check: 0.06 s.

The shared target directory took about 150 MB; each cached program about
2.5 MB.

## Implementation log

- `bindings/python/src/rules_payload.rs` is the extension's result encoder,
  moved out of `rules.rs` unchanged; the driver compiles the same file.
- `PreparedRules.loop_source()` and `rules_loop_driver()` in the extension;
  `python/crn_craft/_rules_loop.py` builds, caches, runs and checks the
  program; `rules.py` adds `compiled=` to `prepare`, `ensemble` and
  `simulate`, and re-exports `compile_loop`, `loop_source`, `CompiledLoop`
  and `CompiledLoopError`.
- `validation/check_compiled_loop.py` (38 rows: 15 identical runs, the
  reached wall time and 22 refusals) passed normally and under `-O` with
  equal reports. With `CARGO_TARGET_DIR` set to the checkout's release
  build, the third-party crates were reused; `crn-craft` and
  `crn-craft-rules` were compiled again, beside the workspace's artifacts
  (29 s instead of 40 s).
- `validation/check_manual.py` sets `CRN_CRAFT_SOURCE` to the checkout,
  skips `compile_loop` blocks without cargo and fails them with
  `--require-cargo`, which CI passes. `validation/check_installed.py`
  checks the source, the driver and two refusals without cargo.

## Limitations

- The fingerprint, the version pin and the self-check catch a program that
  disagrees with the installed extension, but a wheel built from an
  unreleased tree with the same version number as a crates.io release can
  only be checked by the self-check; such builds should use `source=`.
- Mid-run cancellation and reached wall-time limits end the run in both
  routes, but which replicates complete then depends on timing, so only a
  cancel file present before the run and an unreached wall-time limit are
  compared byte for byte.
- Ctrl-C ends the program; partial results are discarded, as in the
  interpreter.
- No release profile tuning (LTO, one codegen unit) was tried; the target
  was met with Cargo's defaults, which keep the target directory shareable.
- Measured on one host and one model; other models, in particular networks
  with many constant-rate channels, can gain less or lose
  ([checklist](../specialized-kernels.md#generated-rules-event-loops)).
