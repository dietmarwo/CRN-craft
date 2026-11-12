# M4d log: speed gate

Task M4d of [the plan](../../plan-0.7.0-engine-merge.md): the unified
engine's best exact interpreted method takes at most 1.25 × the 0.6.0 best on
each of the twelve speed-gate models, timed in one session against a rerun of
`baseline.py` with the kept 0.6.0 build. The per-event work that preceded it is
in the [M4f log](M4f-log.md). After the gate, the coordinator asked for the
full check suite, a law gate of the run/1.0 default (`auto`) and a rerun of the
M5g worker speed comparison with a 1.25 × target; those are recorded here
too.

## Summary

| Item | Result |
| --- | --- |
| Run 1 (CLI of `9c4b5f9`) | **failed**: W5 1.28; the other eleven models 0.66–1.18 |
| Fix round (one) | `57cd574`: group ceilings kept with the composition masses; byte-identical, about 1 % of W5's instructions |
| Run 2, the rerun (CLI of `57cd574`) | **passed**: all twelve within 1.25, largest W3 1.20 |
| Law gate of the run/1.0 default (`auto`) | **passed**: 4,814 comparisons, 0 failures; results equal to the `rssa_cr` and `linear` gates of the selector it picks |
| Worker speed rerun (target 1.25 ×) | **failed** on one of six ratios: `binding12` new `reuse` / old `reuse` 1.55; the other five 1.04–1.15 |

Commits: `122c4e7` (gate ratio in `compare_speed.py`), `57cd574` (fix
round), `a21c358` (both runs' data), `da8a71f` (declaration of the law gate
and the worker rerun), and the commit that adds this log with the results
of both.

## Declared target

M0c's "Speed gate (M4d)" in `validation/merge/RESULTS.md`: twelve models
(W1–W7 at the workload sizes, the five examples), gated boundary = wall time
of one CLI process, median of three after a warm-up, one thread; 0.6.0 best =
the smallest median over its 21 methods; unified best = the smallest median
over the unified engine's methods and selectors (`linear`,
`composition_rejection`, `sorting`, `rssa_cr`, `auto`, `next_reaction`);
gate: unified best ≤ 1.25 × 0.6.0 best per model, against a same-session
rerun of `baseline.py` whose non-timing fields equal `baseline-0.6.0.json`.
The task prompt allows one round of profiling and fixes and one rerun after a
miss; the owner decided during the task that a remaining miss would be
recorded with its hot spots and the release would go ahead with it.

`compare_speed.py` gained `--baseline` (one or more `baseline.py` reports,
merged by model), `next_reaction` among the unified methods, and a `gate`
object per model (reference method and seconds, best method and seconds,
ratio) computed by the declared rule (`122c4e7`, before run 1).

## Runs

Each run: `baseline.py --python /home/claude/m07/venv060/bin/python --binary
/home/claude/m07/bin060/crn-craft --check-against
validation/merge/baseline-0.6.0.json` in four parts of three models (each
under ten minutes, resumable by part), then `compare_speed.py` against the
four parts; everything under the CPU lock (`/home/claude/m07/m4d/run.sh`).
Every part compared 63 rows with `baseline-0.6.0.json`: no difference, in
both runs. The parts are merged by model into
`validation/merge/baseline-0.6.0-m4d-run1.json` and
`baseline-0.6.0-m4d.json` (`/home/claude/m07/m4d/merge.py`); the speed
reports are `speed-m4d-run1.json` and `speed-m4d.json` (their `reference`
field names the part files, which for run 1 now lie in
`/home/claude/m07/m4d/run1/`).

Ratios (unified best / same-session 0.6.0 best):

| Model | Run 1: 0.6.0 best | Run 1: unified best | Run 1 | Run 2: 0.6.0 best | Run 2: unified best | Run 2 |
| --- | --- | --- | ---: | --- | --- | ---: |
| W1 | `direct/ziggurat/fast` 0.0654 | `linear` 0.0504 | 0.77 | `direct/ziggurat/fast` 0.0482 | `auto` (`linear`) 0.0511 | 1.06 |
| W2 | `direct_cr/ziggurat/fast` 0.0894 | `auto` (`rssa_cr`) 0.1032 | 1.15 | `direct_cr/ziggurat/fast` 0.0906 | `rssa_cr` 0.1031 | 1.14 |
| W3 | `cached_direct/ziggurat/fast` 0.2254 | `linear` 0.2651 | 1.18 | `cached_direct/ziggurat/fast` 0.2248 | `auto` (`linear`) 0.2707 | 1.20 |
| W4 | `rssa_cr/ziggurat/fast` 0.1046 | `rssa_cr` 0.1124 | 1.08 | `rssa_cr/ziggurat/fast` 0.1119 | `auto` (`rssa_cr`) 0.1084 | 0.97 |
| W5 | `direct_cr/ziggurat/fast` 0.0602 | `auto` (`rssa_cr`) 0.0770 | **1.28** | `direct_cr/ziggurat/fast` 0.0614 | `auto` (`rssa_cr`) 0.0685 | 1.12 |
| W6 | `rssa_cr/inverse/fast` 0.0641 | `auto` (`rssa_cr`) 0.0610 | 0.95 | `rssa_cr/inverse/fast` 0.0639 | `auto` (`rssa_cr`) 0.0599 | 0.94 |
| W7 | `direct_cr/ziggurat/fast` 0.0202 | `rssa_cr` 0.0195 | 0.97 | `cached_direct/ziggurat/fast` 0.0185 | `rssa_cr` 0.0216 | 1.17 |
| enzyme-mechanism | `direct/ziggurat/fast` 0.0761 | `auto` (`linear`) 0.0505 | 0.66 | `cached_direct/ziggurat/fast` 0.0731 | `sorting` 0.0622 | 0.85 |
| immigration-death | `cached_direct/ziggurat/fast` 0.0695 | `linear` 0.0488 | 0.70 | `cached_direct/ziggurat/fast` 0.0687 | `auto` (`linear`) 0.0474 | 0.69 |
| michaelis-menten | `cached_direct/ziggurat/fast` 0.0575 | `linear` 0.0381 | 0.66 | `direct/ziggurat/fast` 0.0561 | `auto` (`linear`) 0.0368 | 0.66 |
| repressilator | `cached_direct/ziggurat/fast` 0.2676 | `linear` 0.3152 | 1.18 | `cached_direct/ziggurat/fast` 0.2809 | `linear` 0.3178 | 1.13 |
| vilar-noise | `cached_direct/ziggurat/fast` 1.0723 | `linear` 0.9841 | 0.92 | `cached_direct/ziggurat/fast` 1.0514 | `auto` (`linear`) 0.9451 | 0.90 |

Every selector's time is in RESULTS.md ("M4d (speed gate): results").

**Fix round.** W5 under callgrind (RSSA-CR, 510 channels, 4 molecules per
species): 2.1 trials per event, about half of the composition draws rejected
(the group upper bounds, 8, are exact powers of two), every candidate
evaluated exactly, about 10 random words per event. Within the rules of the
task (same arithmetic and words, or run/1.0-only with law gates), the only
change found was in the composition sampler: the ceilings of the active
groups are computed when the cumulative masses are rebuilt instead of per
draw (`57cd574`; byte-identical on the 156 harness requests of
`/home/claude/m07/m4f/harness.py`, the I1 stored requests and the
generated-loop replays; about 1 % fewer instructions on W5). Fewer trials or
words per event would change RSSA-CR's results and were not attempted.

**Noise.** `auto` runs the same computation as the selector it picks, so
each model has two timings of one computation; in run 2 they differ by up to
1.55 × (W6 0.0599 against 0.0930 s, W4 0.1084 against 0.1620 s). Taking the
slower of the two, run 2 would give W1 1.32, W2 1.26, W3 1.21, W4 1.45, W5
1.26, W6 1.46 and W7 1.24. Pooled over both runs (each side's smaller median
of the two runs) the ratios are W1 1.05, W2 1.15, W3 1.18, W4 1.04, W5 1.14,
W6 0.94, W7 1.05, repressilator 1.18 and 0.66–0.90 for the other examples.
The 0.6.0 side, a minimum over 21 methods with close variants, profits from
the same effect. So the gate passes by its declared rule, but the margins of
W2, W3, W5, W7 and the repressilator (1.12–1.20) lie within this host's
run-to-run noise, and W5's move from 1.28 to 1.12 is mostly noise, not the
1 % fix.

**Compiled event loops (informative).** Python boundary, wheel of
`57cd574`, `selector: linear` (what the generated loops run), median of three
`simulate()` calls after a warm-up (`/home/claude/m07/m4d/compiled.py`,
`compiled.json`):

| Model | Interpreted (`linear`) s | Compiled s | Compiled / interpreted | 0.6.0 best (Python) s | Compiled / 0.6.0 best |
| --- | ---: | ---: | ---: | ---: | ---: |
| W1 | 0.0477 | 0.0353 | 0.74 | 0.0440 | 0.80 |
| W2 | 0.1915 | 0.1727 | 0.90 | 0.0882 | 1.96 |
| W3 | 0.2627 | 0.1808 | 0.69 | 0.2502 | 0.72 |
| W4 | 0.7850 | refused (source over 1 MiB) | | 0.1014 | |
| W5 | 0.2916 | refused | | 0.0527 | |
| W6 | 0.5022 | refused | | 0.0468 | |
| W7 | 0.0306 | 0.0350 | 1.15 | 0.0154 | 2.27 |
| enzyme-mechanism | 0.0454 | 0.0370 | 0.81 | 0.0411 | 0.90 |
| immigration-death | 0.0453 | 0.0410 | 0.90 | 0.0403 | 1.02 |
| michaelis-menten | 0.0389 | 0.0309 | 0.79 | 0.0266 | 1.16 |
| repressilator | 0.3319 | 0.2080 | 0.63 | 0.2604 | 0.80 |
| vilar-noise | 1.0012 | 0.6126 | 0.61 | 1.3030 | 0.47 |

The compiled loops run the linear scan, so on W2 and W7 (RSSA-CR models) they
are slower than the 0.6.0 best; W4–W6 exceed the generated-source limit.

## Remaining hot spots (for a follow-up release)

- W1 (`linear`, 1.06): two ChaCha8 words per event (about 3.7 ns per word
  against 1.5 ns for 0.6.0's xoshiro256++) and the Ziggurat wait.
- W3 and the repressilator (12 channels, Hill rates, `linear`, 1.20 and
  1.13): one mispredicted scan exit per event (scans of more than four rates
  keep the early exit), the running-sum resume, the Hill programs of the
  affected channels, two ChaCha8 words per event.
- W2, W5, W7 (`rssa_cr`, 1.14, 1.12, 1.17): trials per event, composition
  rejection (W5: half the draws, bounds at exact powers of two), exact
  evaluation of every candidate, about 10 words per event on W5.
- W4 (`rssa_cr`, 0.97): eight-species count updates and interval tests per
  event; range refreshes on 7 % of the events.
- Levers that change results or the build: RSSA-CR interval widths or group
  ceilings off the powers of two, prefix-sum selection as a new selector, a
  cheaper generator under a new sampler identity (D4 fixes ChaCha8), LTO or
  one codegen unit (which also changes the 0.6.0 comparison build).

## Law gate of the run/1.0 default

Declared in RESULTS.md (`da8a71f`) before the run. `check_laws.py` with the
Python backend on the wheel of the engine sources of `57cd574`, request
method `{"kind": "static_direct"}` without a selector, all twelve models,
the gate seeds, 2 threads, in resumable pieces under the CPU lock
(`/home/claude/m07/m4d/gate/`); report `validation/merge/laws/m4d-auto.json`.

**passed**: 4,814 comparisons, 0 failures, largest |z| 3.88 (W5). `auto` picked `rssa_cr` for W2 and W4–W7 and `linear` for the other seven models, as declared. Sampler, chunk layout, event counts and power sums S1–S4 equal those of `laws/m4f-rssa_cr.json` (the `rssa_cr` picks) and `laws/m4c-linear.json` (the `linear` picks) for all twelve models. The `linear` selector keeps the sampler identity `chacha8-u53-ziggurat/rules-static-direct/1` of rules-run Static Direct.

The run was interrupted once (the task's agent was stopped) and resumed from
its stored chunks (W1–W5 finished before, W6–vilar-noise after); a resumed
chunk is the same computation, and the equality with the earlier gates
confirms it.


## Worker speed rerun

Declared in RESULTS.md (`da8a71f`) before the run. The script, the
migrated worker and the fixtures come from `merge/M5cdg` at `97a60c3`
(`git archive` into `/home/claude/m07/m4d/workers/root/`); the worker's
`crn-craft` path dependency in that scratch copy points to
`/home/claude/fb2/engine/crates/crn-craft` (a symlinked `crates/` failed:
the crate inherits its edition from the workspace root), built with
`cargo build --release --locked` against the worker's own lockfile. Old
workers: the M5g binaries in `/home/claude/m07/m5cdg/old-worker/`. Run
2026-11-12 09:26–09:28 under the CPU lock; report
`validation/merge/parameter-workers-m4d.json`. Wall time of the 129 calls,
median of 3, seconds (bind and simulation: sums over the calls):

| Case | old fresh | old reuse | old generated | new fresh | new reuse | new fresh / old fresh | new reuse / old reuse | new reuse / old generated |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| birth_death | 0.169 | 0.168 | 0.149 | 0.194 | 0.185 | 1.15 | 1.10 | 1.24 |
| sparse33 | 2.909 | 2.894 | 2.979 | 3.146 | 3.043 | 1.08 | 1.05 | 1.02 |
| binding12 | 0.052 | 0.034 | 0.033 | 0.054 | 0.053 | 1.04 | **1.55, failed** | 1.59 |

Against the M5g row (before M4c/M4f): sparse33 new `reuse` 5.326 → 3.043 s
(`auto` now runs RSSA-CR on its 66 channels; simulation 2.85 s against
0.6.0's 2.82 s), birth_death 0.296 → 0.185 s. The target fails on
`binding12` `reuse`: of the 0.019 s gap, binding is 0.011 s (new 0.0137 s,
about 106 µs per call, against 0.0024 s: `Prepared::bind_parameters`
re-serializes and prepares the whole model) and simulation 0.008 s (0.0248
against 0.0170 s for 32 short paths per call: per-call preparation of the
run and per-path setup outweigh the events). Within each engine the routes
returned identical moments at every point. Not optimized further (owner
decision).

## Checks

| Check | Result |
| --- | --- |
| Speed gate, run 2 (`compare_speed.py`, declared rule) | passed, all twelve ≤ 1.25 (largest W3 1.20) |
| `baseline.py --check-against baseline-0.6.0.json`, every part of both runs | 63 of 63 rows equal per part |
| Byte identity of the fix round: 156 harness requests, I1 stored requests, generated-loop replays | equal |
| Law gate of the run/1.0 default (`laws/m4d-auto.json`) | 4,814 comparisons, 0 failures; equal to the gates of the picked selectors |
| Worker speed rerun (`parameter-workers-m4d.json`) | failed on `binding12` `reuse` (1.55); other five ratios ≤ 1.15 |
| Rust on the integration of `merge-0.7` (`da8a71f`) with `merge/M6` (`5c30f0a`): fmt, clippy `-D warnings`, rustdoc `-D warnings`, `cargo test --release --workspace --locked` | pass; 532 tests, 0 failures |
| Python validation suite | run on the integrated branch; recorded in the [log index](README.md) |

## Run history

- A rate limit and two server restarts stopped the task between steps; the
  run-1 `baseline.py` parts survived a restart and were used (the run script
  skips finished parts).
- Another agent's builds and tests held the CPU lock at times; the timed
  runs ran only under the lock.
- Working files: `/home/claude/m07/m4d/` (run scripts, parts, logs, gate
  state, worker scratch tree, compiled-loop timings).

## Deviations

- The worker comparison script is not on `merge-0.7` yet (it comes with
  M5g); it ran from a scratch copy of `merge/M5cdg` against this checkout's
  engine, with the path edit above. Its report is committed here as
  `parameter-workers-m4d.json`.
- The law gate of the default ran with the Python backend (resumable) on
  the wheel of the same sources rather than the CLI backend (see the
  declaration).

## Migration notes

None beyond M4f's (the fix round keeps every output).
