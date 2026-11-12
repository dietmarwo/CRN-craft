# M9 log: downstream

Task M9 of [the plan](../../plan-0.7.0-engine-merge.md), branch `merge/M9` (worktree
from `merge-0.7` at `ef899e5`). Nothing under `tutorials/`, `examples/` or
`integrations/` uses `crn_craft._legacy`, a hidden old command or the network Rust
types any more.

## What changed

- **`integrations/diffsol`**: no longer calls `crn_craft::parse_model` or the network
  types (`contracts::Model`, which M8b removes). The network `model/0.1` population
  limits and species names are read from the JSON itself (`network_species`);
  everything else goes through `Kinetics::from_bytes`, which reads any model version.
  `tests/kinetics.rs` uses `Kinetics::from_bytes` instead of `parse_model` and
  `Kinetics::compile`. Behavior and output unchanged (10 of 10 tests).
- **`integrations/jobs`, `validation/parameter-optimization/worker`**: already on
  `crn_craft::rules` and `rules::convert::from_network` (M5g); no change needed.
- **`tutorials/optimization`**:
  - `calibration.py` (C1): `crn.ensemble` plus `bind_parameters`. The 0.6.0
    `fast`/`checked` policy axis does not exist any more (one generator, D4); the
    option `--include-checked` became `--include-next-reaction`, a second exact method
    with its own stream, rows keyed `direct` / `next_reaction`. The "ordinary replay"
    check compares a fresh preparation of the fitted model with the binding: equal
    statistics, events, sampler, schema; only `model_sha256` and `run_sha256` may differ
    (a binding serializes the bound model itself).
  - `oscillators.py` (C3): the three-species Hill model is built as a `model/1.0`
    document (propensities as the converter writes them); generator audit through
    `Model.analyze()`; candidates through `bind_parameters`; paths from `result.paths`;
    the Vilar control uses `load_example`. The reference backend `crn_reference` is
    unchanged.
  - `sensitivity.py` (C2): converted enzyme example; SSA transfer through
    `bind_parameters`; the result's observables are sorted by name, so they are
    reordered to the ODE's species order (S, E, ES, P, which is also asserted).
  - `rules.py` (C4): `crn-craft run` with `crn-craft.run/1.0` and the model's schema
    string `crn-craft.model/1.0` (was the `rules` alias and `rules-run/0.1`).
  - `common.py` `checked` reads `Result.metadata`; `record.py`, `test_science.py`
    follow the C1 change; lesson texts C1-C4, README and VALIDATION updated;
    `results/verified.json` regenerated from complete smoke and study runs.
- **`examples/circuit`**: `screen.py` on `crn.ensemble` (converted models);
  `fit_ode.py` feeds the converted model to `crn-diffsol`.
- **`examples/tutorial`**: `workdir.py` writes converted models and `run/1.0`
  requests; `make_assets.py` sessions use `run`, `run --archive`, `verify`, `run` with
  method `fsp` and `tau`, and figures use `crn.ensemble` and `Model.fsp`. Session 06
  (`simulate-cached`) is replaced by `06-replay` (same request, same bytes).
  Figures `immigration-death`, `fsp-vs-ssa`, `repressilator`, `tau-accuracy`,
  `first-passage` and `numbers.json` regenerated and committed (content checked by eye:
  same curves and annotations; the sampled values differ, below).

## Tests first

`test_science.py` was adapted first and failed on the old code paths
(`binding_model()["rules"]`, `checked` on a result, the recorded samplers); the
verified-record test failed until `verified.json` was regenerated. No threshold was
changed.

## Checks

| Check | Result |
| --- | --- |
| `cargo build/fmt/clippy -D warnings/test --release --locked` in `integrations/diffsol`, `integrations/jobs`, `validation/parameter-optimization/worker` (shared temporary target) | all pass; tests 10 (diffsol), 3 (jobs), worker 0 |
| C1-C4 smoke (`--threads 2`; C1 with `--include-next-reaction`) | pass |
| C1-C4 study presets | pass (C1 7.6 s, C2 2.7 s, C3 9.5 s, C4 3.3 s here; they ran) |
| `test_science.py` (15 tests), normally and `-O` | pass |
| `validation/check_jobs.py` (worker `crn-job-worker`, pinned runner), normal and `-O` | pass, 15 rows, 17 jobs, 18 manifests; reports identical |
| `examples/circuit/screen.py`, `fit_ode.py` | pass (gate true) |
| `make_assets.py` | runs; figures regenerated |
| `cargo fmt --all -- --check`, clippy `-D warnings`, `cargo test --release --workspace --locked` (I1, I3 included), rustdoc `-D warnings`, `cargo build --release --workspace --locked` | pass |
| I2 `check_manual.py --binary ... --require-cargo`, normal and `-O` | pass, reports identical |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |

## Numbers that changed (D4: new stream)

Thresholds are unchanged. Old (0.6.0 record) -> new:

| Quantity | Old | New |
| --- | --- | --- |
| C1 largest held-out mean error, smoke / study (gate 6) | 2.0 / 2.2 | 2.3 / 2.8 |
| C3 control, largest standard errors, smoke / study (gate 6) | 1.4 / 1.7 | 1.1 / 1.8 |
| C3 occupied archive cells (of 16) | 3 / 3 | 2 / 2 |
| C3 quality-optimal design, held-out oscillatory share | 31% / 78% | 62% / 66% |
| tutorial `immigration_death` events, final mean | 87508, 19.705 | 88037, 20.09 |
| tutorial tau variance ratio, eps 0.003 / 0.01 / 0.03 / 0.1 | 0.957 / 0.991 / 1.183 / 1.678 | 1.113 / 1.010 / 1.093 / 1.834 |
| tutorial first-passage detected, grid 0.05 / 0.2 / 1.0 | 460 / 459 / 441 | 464 / 460 / 440 |

Unchanged: C2 (audits 9.1e-10, 6.5e-11, 3.2e-12; twin; fits), C4 (rules engine), the
FSP numbers of the tutorial, and the thread-independence check (equal digests for
1, 2, 8 threads; the digest value changed with the stream).
The tau `leaps` field of `numbers.json` is gone (the unified result has no leap count).

## Limitations and notes for later tasks

- **Terminal screenshots** (`docs/images/tutorial/01..11-*.png`) are not re-rendered:
  Playwright/Chromium is not installed here. The `*.terminal.html` and `*.txt`
  transcripts of the new sessions were generated; `render_terminals.py` must run
  where a browser exists. Until then `04-simulate` .. `11-export` PNGs show old
  commands, and `06-cached.png` should be replaced by `06-replay.png`. M10 must update
  `docs/tutorial.md` (commands, `06-cached`, the tau and first-passage numbers, the
  missing-parameter diagnostic, now `CRN-R012` at `/rules/0/rate`).
- `integrations/diffsol/tests/cli.rs` and `tests/kinetics.rs` include the network
  files of `models/examples/` (and inline network JSON for the population-limit
  tests). If M8b converts `models/examples/` to `model/1.0`, the limit tests need
  `models/legacy/` inputs; the kinetics test works with either.
- `validation/check_ode.py` (M8a) still feeds network models to `crn-diffsol`; this
  works because diffsol keeps reading `model/0.1` through `Kinetics::from_bytes`.
  M8b must keep that branch (via the private converter) or M8a must port the script.
- C3 needs `crn_reference`, which the stale `fcmaes-rust` checkout here lacks; it was
  built from the upstream `main` file `examples/crn_reference.rs` (109 lines) in a
  temporary copy.
- No script under `validation/` needed a change for M9.

## Migration notes

`tutorials/optimization/run.py C1 --include-checked` is now `--include-next-reaction`;
the recorded policies are `direct` and `next_reaction`.
