# Static Direct defaults (rules-run/0.6): implementation log

6 November 2026. rules-run/0.6 makes the S3 opt-ins of rules-run/0.4,
Ziggurat waits and sparse trackers, the defaults of `static_direct`, and
Python's `ensemble(method="static_direct")` writes 0.6 requests. This log
declares the gates and the speed target before any run, then records the
outcomes.

## Declaration (before any gate or timing run)

**Gates.** The samplers the 0.6 default selects are the S3 samplers,
`chacha8-u53-ziggurat/rules-static-direct/1` and
`chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1`. These gates
cover Ziggurat waits and sparse trackers:

1. `validation/check_mtdna.py`, `check_generic_sampler`: the K = 20
   master-equation gate on `models/rules/mtdna-three-class-cme.rules.json`
   (S3; rows `cme/static-ziggurat/*`), 32,768 paths, the per-row DKW bound
   `dkw(36)` at family alpha 0.001, and the worker-replay and event-log
   integral checks (1 and 8 workers identical, an event log changes no
   trajectory or tracker output, integrals reconstructed from the full log
   within 1e-12 relative; row `static-ziggurat/worker-replay-and-log-integrals`).
   Its first-passage tracker is not eligible for sparse integration, so the
   law rows ran the Ziggurat identity; the replay rows run the sparse one.
2. The Rust release tests: the generated loop and kernel replay every
   Static Direct policy byte for byte against the interpreter; the S3 and
   review tests of the sparse identity and of Ziggurat infinite waits.
3. `validation/check_rules_python.py`: the Ziggurat sparse-tracker callback
   replay.

For the 0.6 default the same gate runs three more times, unchanged in
bound, size and alpha, with fresh seeds:

- `cme/static-default-0.6/*`: the 0.6 default with the S3 first-passage
  tracker, seed `SEED + 5`; expected identity
  `chacha8-u53-ziggurat/rules-static-direct/1` and no `CRN-R051`;
- `cme/static-default-0.6-sparse/*`: the 0.6 default with the two count
  trackers, seed `SEED + 6`; expected identity
  `chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1`, so the law
  is also checked on the path with sparse projection refresh;
- `static-default-0.6/worker-replay-and-log-integrals`: the replay and
  integral checks for the 0.6 request.

The existing rows and their bounds stay as they are. If any row fails, the
work stops and the failure is reported; no bound is changed.

**Speed target.** On the 1000-cell three-class case of the methods section
of `examples/mtdna/performance.py` (takeover stop, three years, monthly
reports, seed 17, two threads, two count trackers), the 0.6 default
`static_direct` is at least **1.10×** faster than the 0.5 default, in
median wall nanoseconds per event over three runs each in alternating
order (0.5 first), under the CPU lock. Wall seconds are reported as well.
Script: `validation/benchmark_static_direct_defaults.py`. The background
measurement in the brief was 78 against 67 ns per event (1.16×).

## Outcomes

### Gates

`python validation/check_mtdna.py --binary target/release/crn-craft
--section generic --output validation/mtdna/gates-static-direct-defaults.json`,
normally and under `-O` (identical reports), with the CLI built from
`099c08d`. `--section generic` runs only `check_generic_sampler`; CI runs
the whole script. All 29 rows passed; nothing was rerun.

| Case | Request | Identity | Law rows | Largest distance | Bound | Replay row |
| --- | --- | --- | ---: | ---: | ---: | --- |
| `static-ziggurat` (S3) | 0.4, both opt-ins, first-passage tracker | Ziggurat, `CRN-R051` | 9 | 0.00609 (mutants, t = 3) | 0.01306 | passed |
| `static-default-0.6` | 0.6 default, first-passage tracker | Ziggurat, no warning | 9 | 0.00527 (mutants, t = 3) | 0.01306 | passed |
| `static-default-0.6-sparse` | 0.6 default, two count trackers | Ziggurat, sparse trackers | 9 | 0.00673 (mutants, t = 1) | 0.01306 | — |

The S3 rows reproduce the recorded S3 result (largest distance 0.00609).
The replay rows ran the sparse identity for both requests. The Rust replay
tests and `check_rules_python.py` (including the Ziggurat sparse-tracker
callback replay and the new checks of the Python default) passed.

### Speed

`python validation/benchmark_static_direct_defaults.py --output
validation/mtdna/static-direct-defaults-timing.json` under the CPU lock, with
a wheel built from this tree: Linux, two CPUs, Python 3.12.3, commit
`099c08d`. 1000 cells, two threads, two count trackers, takeover stop.

| Run | 0.5 default s | ns per event | 0.6 default s | ns per event |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 7.00 | 48.8 | 5.79 | 40.3 |
| 2 | 6.81 | 47.5 | 5.99 | 41.6 |
| 3 | 7.04 | 49.1 | 6.02 | 41.9 |
| median | 7.00 | 48.8 | 5.99 | 41.6 |

Events: 143,436,887 (0.5, 99 takeovers) and 143,795,432 (0.6, 97
takeovers); the two requests use the random stream differently. Speedup
**1.17×** in median ns per event (1.17× in wall seconds), above the declared
1.10× target. The absolute times are lower than the 78 and 67 ns per event
in the performance analysis, which were measured on another two-CPU host
with Python 3.13; the ratio agrees (1.16× there).

## Validation limitations

- The law gate is the K = 20 three-class chain; the sparse path's speed and
  identity were measured on the K = 1000 model only.
- The timing ran on a shared two-CPU host under its CPU lock; other hosts
  will give other absolute times.
