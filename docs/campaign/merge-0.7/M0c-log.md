# M0c log: baseline and references of the 0.6.0 engine

Task M0c of [the merge plan](../../plan-0.7.0-engine-merge.md). The network
engine of 0.6.0 is deleted in M8, so everything later tasks compare against
was recorded now, from the old engine, and committed under
`validation/merge/`.

## What changed

| File | Content |
| --- | --- |
| `validation/merge/RESULTS.md` | Declarations (law gate, speed gate, selector rule, grid, seeds), then the results |
| `validation/merge/models.json` | Machine-readable grid: path, hash, `end`, report times, replicates, chunk, seeds, FSP bounds, timing sizes |
| `validation/merge/merge_common.py` | Standard-library helpers: request builder, host facts, exact moment formulas, law-gate comparisons |
| `validation/merge/references.py` | Writes `validation/merge/reference/` |
| `validation/merge/baseline.py` | Writes `validation/merge/baseline-0.6.0.json` |
| `validation/merge/kinetics-reference/` | Standalone Cargo package printing `Kinetics::evaluate` as JSON (no Python or CLI path exposes `kinetics`) |
| `validation/merge/reference/` | 52 files, 513 KB (616 KB on disk) |
| `validation/merge/baseline-0.6.0.json` | 252 timing rows, 537 KB |

Order of commits: the declarations (`RESULTS.md`, `models.json`) first, before
any reference or baseline run; then the scripts; then the outputs.

## Old-engine build kept for later tasks

From the worktree root, with the CI commands:
`cargo build --release --workspace --locked` and
`maturin build --release --locked --out dist` (maturin of the engine venv).

- `/home/claude/m07/bin060/crn-craft`: the 0.6.0 CLI.
- `/home/claude/m07/bin060/crn_craft-0.6.0-cp310-abi3-manylinux_2_34_x86_64.whl`: the wheel.
- `/home/claude/m07/bin060/kinetics-reference`: the kinetics printer, built with
  `cargo build --release --locked --offline --target-dir <worktree>/target`
  inside `validation/merge/kinetics-reference`.
- `/home/claude/m07/venv060`: Python 3.12 with the wheel, numpy 2.3.5,
  pandas 3.1.1, xarray 2026.9.0, jsonschema 4.26.0, rebop 0.10.0 and ruff
  0.11.13 (the CI pins).

Both scripts run under any Python 3.10+ with the standard library and start
a worker under `--python` for the parts that import `crn_craft`; `-O` is passed
on to the worker. So M4d can rerun the 0.6.0 baseline on this host after the
network engine is gone.

## Sizing probes before the declarations

Two probes ran before `RESULTS.md` was committed, only to choose sizes; their
output was not stored.

1. Per-path cost of `direct` and `cached_direct` (checked, inverse, one
   thread) at the timing horizons. Milliseconds per path, `direct` /
   `cached_direct`: W1 1.44/1.33, W2 12.9/3.9, W3 8.4/7.0, W4 16.6/13.5,
   W5 97.7/21.8, W6 135.5/42.9, W7 8.4/4.4, immigration-death 0.015/0.013,
   michaelis-menten 0.005/0.006, enzyme-mechanism 0.008/0.007,
   repressilator (end 10) 0.64/0.51, vilar-noise (end 50) 17.0/11.0. Event
   counts of the two were equal for every model.
2. FSP feasibility of the bounds (states, edges, work, leakage) and the
   dynamics of vilar-noise at 10 h, W3 at 30 and W4 at 20 (2,000 replicates).

## Runs

All runs used the kept 0.6.0 build, under `flock /home/claude/craft5/cpu.lock`.

```sh
python3 validation/merge/references.py --python /home/claude/m07/venv060/bin/python \
  --binary /home/claude/m07/bin060/crn-craft \
  --kinetics-binary /home/claude/m07/bin060/kinetics-reference
python3 validation/merge/baseline.py --python /home/claude/m07/venv060/bin/python \
  --binary /home/claude/m07/bin060/crn-craft
```

- References: 39 min 57 s wall, 74.5 CPU minutes (user 74 min 2 s, sys
  29 s), 17.4 billion events. The estimate in `RESULTS.md` was about 3,300 CPU
  seconds of simulation plus overhead; vilar-noise alone took about 26 CPU
  minutes instead of the estimated 7, because its event rate is highest in
  the first hours (55,000 events per path in 10 h; scaling the probe at 50 h
  linearly had suggested about 18,500). The run stayed within the budget of
  about 90 CPU minutes.
- Baseline: 18 min 21 s wall.

## Results

### References

- Direct identity: for all 12 models the first 1,000 replicates of chunk 0
  rerun with `direct` gave the same saved states and per-replicate event
  counts as `cached_direct`.
- FSP: 5 models (W1, W4, enzyme-mechanism, immigration-death,
  michaelis-menten); largest leakage + mass defect 6.1e-13 (limit 1e-9).
- Self-check of the SSA reference against the FSP law: 178 comparisons,
  0 failures, largest |z| 2.08 (enzyme-mechanism, mean of S at t = 10).
- Kinetics: accepted for W1, W2, W3 and the five examples; refused, as
  declared, for W4 (`CRN-E080`, Michaelis–Menten law consuming more than its
  substrate) and W5, W6, W7 (`CRN-E080`, size cap).
- Analysis: 12 reports, all exit 0.

### Baseline (gated CLI boundary, median of three, one thread)

| Model | Best 0.6.0 method | Seconds | Events | Best checked method | Seconds | Best `direct_cr` | Best at the Python boundary | Seconds |
| --- | --- | ---: | ---: | --- | ---: | ---: | --- | ---: |
| W1 | `cached_direct/ziggurat/fast` | 0.0478 | 1,272,658 | `cached_direct/ziggurat/checked` | 0.0665 | 0.0733 | `direct/ziggurat/fast` | 0.0510 |
| W2 | `direct_cr/ziggurat/fast` | 0.0942 | 1,235,619 | `rssa_cr/ziggurat/checked` | 0.1355 | 0.0942 | `direct_cr/ziggurat/fast` | 0.0888 |
| W3 | `cached_direct/ziggurat/fast` | 0.2224 | 3,656,543 | `sorting_direct/ziggurat/checked` | 0.3607 | 0.3385 | `cached_direct/ziggurat/fast` | 0.2133 |
| W4 | `rssa_cr/ziggurat/fast` | 0.1062 | 1,227,674 | `rssa_cr/ziggurat/checked` | 0.1309 | 0.8529 | `rssa_cr/ziggurat/fast` | 0.1132 |
| W5 | `direct_cr/ziggurat/fast` | 0.0677 | 325,267 | `next_reaction/inverse/checked` | 0.0909 | 0.0677 | `direct_cr/ziggurat/fast` | 0.0504 |
| W6 | `rssa_cr/ziggurat/fast` | 0.0630 | 229,211 | `rssa_cr/ziggurat/checked` | 0.0722 | 0.2112 | `rssa_cr/ziggurat/fast` | 0.0437 |
| W7 | `rssa_cr/ziggurat/fast` | 0.0209 | 108,058 | `rssa_cr/ziggurat/checked` | 0.0258 | 0.0224 | `cached_direct/ziggurat/fast` | 0.0154 |
| enzyme-mechanism | `cached_direct/ziggurat/fast` | 0.0765 | 460,645 | `direct/ziggurat/checked` | 0.0854 | 0.1209 | `sorting_direct/ziggurat/fast` | 0.0424 |
| immigration-death | `cached_direct/ziggurat/fast` | 0.0724 | 626,862 | `cached_direct/ziggurat/checked` | 0.0824 | 0.1062 | `cached_direct/ziggurat/fast` | 0.0361 |
| michaelis-menten | `direct/ziggurat/fast` | 0.0612 | 307,345 | `next_reaction/inverse/checked` | 0.0658 | 0.0738 | `direct/ziggurat/fast` | 0.0247 |
| repressilator | `sorting_direct/ziggurat/fast` | 0.2867 | 4,194,928 | `sorting_direct/ziggurat/checked` | 0.4241 | 0.4035 | `cached_direct/ziggurat/fast` | 0.2693 |
| vilar-noise | `cached_direct/ziggurat/fast` | 1.0809 | 14,224,888 | `cached_direct/ziggurat/checked` | 1.5616 | 2.1004 | `cached_direct/ziggurat/fast` | 1.0963 |

The 0.6.0 best is always a `fast` variant with Ziggurat waits. If the unified
engine has no unchecked policy, the speed gate compares it against these
`fast` times, not against the checked column.

**Selector rule (M4c):** both kept.

- Sorting: t_cr / t_sort = 1.315 on W3 (0.3385 s against 0.2574 s), the only
  workload at or above 1.3 (W1 1.158, W4 1.090, W2 0.746, W7 0.890, W6 0.598,
  W5 0.367). The margin is small; at the informative Python boundary the
  ratio on W3 is 1.384, and all three repetitions of both methods on W3 lie
  apart (direct_cr 0.3314–0.3455 s, sorting 0.2506–0.2662 s), so the result
  is not a single outlier.
- RSSA-CR: t_cr / t_rssa = 8.03 on W4 and 3.35 on W6.

Timing noise: the relative spread (max − min)/median of the three CLI
repetitions has a median of 10% over the 252 rows; 33 rows exceed 30%, with
single outliers up to 84% (the machine is shared, and processes that do not
take the lock can still run). The median of three absorbs a single outlier per row. The
determinism reruns (W1, W7, immigration-death, about 30 minutes later)
reproduced the best times within 7% (W1 0.0477 and 0.0499 s against 0.0478 s;
W7 0.0201 and 0.0200 s against 0.0209 s; immigration-death 0.0719 and
0.0772 s against 0.0724 s); the best method changed between near-equal
variants (W1 `direct` against `cached_direct`, W7 `cached_direct` against
`rssa_cr`, all Ziggurat and fast).

### M0c gate

| Check | Result |
| --- | --- |
| Outputs exist and are committed, size | `reference/` 513 KB, baseline 537 KB, 1.05 MB together (limit about 20 MB) |
| `references.py` rerun on immigration-death, michaelis-menten, enzyme-mechanism, W1 | 24 files byte-identical to the committed ones |
| Same under `python -O` | 24 files byte-identical |
| `baseline.py` rerun on W1, W7, immigration-death, all non-timing fields | 63 of 63 rows identical, host and tools identical |
| Same under `python -O` | 63 of 63 rows identical |
| Reference self-check against FSP | 178 comparisons, 0 failures |
| `ruff format --check --exclude '*.md'` and `ruff check` (0.11.13) on `validation/merge` | pass |
| `cargo fmt --check`, `cargo clippy -D warnings` on `kinetics-reference` | pass |

No Rust code of the workspace changed. The brief's checks, run once at the
end on the worktree under the CPU lock: `cargo fmt --all -- --check` pass;
`cargo clippy --workspace --all-targets --locked -- -D warnings` pass;
`cargo test --release --workspace --locked` pass (I1 and I3 included);
`RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` pass;
`ruff format --check --exclude '*.md' .` and `ruff check .` pass;
`check_manual.py --binary target/release/crn-craft --require-cargo` (I2)
passed 22 of 22 blocks, normally and under `python -O`, identical reports.

## Deviations from the plan and choices made

1. **Reference law computed with `cached_direct`.** The plan says "200,000
   replicates of `direct`". `cached_direct` with the checked policy and inverse
   waits produces bit-identical paths (documented in `alternative-ssa.md`) at
   up to 4.5 times less cost on W5 and W6, which made 200,000 replicates
   affordable for every model. The identity is checked per model (1,000
   replicates) and recorded in each `*.moments.json` (`direct_identity`).
2. **Law horizons shorter than the timing horizons** for W1–W7 (W1 10, W2 20,
   W3 30, W4 20, W5 2, W6 2, W7 10 instead of 100 or 20), so that all twelve
   models keep 200,000 replicates within the budget. The timing sizes are those
   of `acceleration_workloads.py`.
3. **Ziggurat waits included in the baseline**: 21 methods instead of the
   11 named (5 × `fast` on/off + `next_reaction`). The 0.7.0 default uses
   Ziggurat waits, and without them the "0.6.0 best" would be slower than the
   0.6.0 engine really is, which would loosen the speed gate.
4. **Two timing boundaries.** The speed gate is declared on the CLI
   boundary (one `crn-craft ssa` process per run), which both engines share
   without depending on the shape of the Python API; the prepared Python call
   is recorded as information. M4d is to rerun `baseline.py` in its own
   session and compare against that rerun, because timings on this shared
   machine drift between sessions; the rerun's non-timing fields must equal the
   committed file.
5. **FSP comparisons use the exact law's standard errors** (σ²/n and the
   law's fourth moment), not the sample's, so a rare state absent from one
   sample cannot fail a comparison with a zero standard error.
6. **Kinetics also at a probe state** (count of species i = initial + i + 1),
   because at the initial state many terms vanish (zero counts).
7. **Seeds per chunk.** Each model's 200,000 replicates are 20 (W5: 40)
   chunks with root seeds 1,000,000 + k, so every replicate's states can be
   kept for the fourth moments. The engine under test uses 2,000,000 + k with
   the same chunks; this is declared, so seeds cannot be chosen after a result.
8. Example-model timing sizes are not in the plan: 10,000 replicates
   (enzyme-mechanism, immigration-death, michaelis-menten), 1,000
   (repressilator) and 256 (vilar-noise) at `end = 10`, so the fastest
   method takes 0.06–1.1 s.

## Limitations

- The law gate compares marginal means and variances per species and time;
  it does not compare covariances or whole distributions (FSP results are
  stored, so M5a compares probabilities directly).
- FSP references exist for 5 of the 12 models; the others exceed the
  100,000-state limit.
- The CLI timings of the example models are dominated by writing the
  per-replicate records of 10,000 replicates (immigration-death: 0.072 s at
  the CLI, 0.036 s at the Python boundary).
- `kinetics-reference` compiles only against the 0.6.0 crate API
  (`parse_model`, `kinetics::Kinetics` on a network model). After M5e and M8 it
  will no longer build; its output is committed, and the 0.6.0 binary is kept
  at `/home/claude/m07/bin060/kinetics-reference`.
- Timings come from one shared two-CPU virtual machine (Intel Xeon @ 2.10 GHz).

## Notes for later tasks

- `check_laws.py` (M2) should import `merge_common`: `load_manifest`,
  `cell_statistics` / `statistics_from_power_sums`, `compare_to_reference`,
  `compare_to_law`, `summarize`. References: `reference/<name>.moments.json`
  (fields `mean`, `variance`, `mean_se`, `variance_se`, indexed
  [report time][species] in the order of `species`) and
  `reference/<name>.fsp-moments.json` (`mean`, `variance`,
  `fourth_central_moment`).
- Keeping every path of 200,000 replicates in one run is too large for W5
  (200,000 × 5 × 256 counts); run the declared chunks with seeds
  `gate_seed + k` and accumulate integer power sums as `references.py` does.
- The species order of the references is the 0.6.0 model order; converted
  models must map observables by name.
- M5a: `reference/<name>.fsp.json.gz` holds the full `fsp-result/0.1`
  (states, probabilities, error terms); compare probabilities to 1e-10.
- M5b: `reference/<name>.analysis.json.gz` holds the exact `crn-craft analyze`
  output.
- M5e: `reference/<name>.kinetics.json` holds `rhs`, `state_jacobian`
  (rows species, columns species) and `parameter_jacobian` (columns in sorted
  parameter order, listed in `parameters`) at the `initial` and `probe` points.
- M4d: baseline rows carry the exact request (`request`), so the unified
  engine can be run on the same sizes; `timing_summary.cli` holds the 0.6.0
  best per model.

## Migration notes

None: M0c adds validation files only.
