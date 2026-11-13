# Implementation log: one engine, one API (release 0.7.0)

This log records how the plan in
[`docs/plan-0.7.0-engine-merge.md`](../../history/plan-0.7.0-engine-merge.md) was
carried out. Each task has its own log in this directory; this page gives the
order, the outcome of each task and every deviation from the plan.

## Execution notes

- All tasks were committed on one integration branch, `merge-0.7`, and merged
  into `main` with a single pull request, followed by a release-preparation
  pull request. This replaces the plan's rule of one pull request per task.
- The stopgap branch `docs/two-engines` (task M0a) was not merged: the merged
  release replaces the explanation of two engines it added.

## Tasks

| Task | Log | Outcome |
| --- | --- | --- |
| M0c | [log](M0c-log.md) | 0.6.0 baselines and references recorded (law references, speed baseline, FSP references) before anything changed |
| M1 | [log](M1-log.md) | The rules engine moved into `crn_craft::rules`; `crn-craft-rules` became a facade; no behavior change (I1–I3) |
| M2 | [log](M2-log.md) | Converter `model/0.1` → rules model; law gate on every convertible model passed |
| M3 | [log](M3-log.md) | `model/1.0`, `run/1.0`, `result/1.0`, `check/1.0` as new versions of the rules contracts |
| M4a | [log](M4a-log.md) | Limits of `model/1.0` raised to 4,096 types and 65,536 rules |
| M4b | [log](M4b-log.md) | Higher-order rules (`count`, the `Multi` shape) |
| M4c | [log](M4c-log.md) | Static Direct selectors; the first RSSA-CR law gate failed on vilar-noise (clock), kept as a failed row |
| M4e | [log](M4e-log.md) | Quantile statistics (FR1a), NumPy type 7 bit for bit |
| M4f | [log](M4f-log.md) | Per-event cost cut; RSSA-CR clock rule for `run/1.0`; `auto` refit; RSSA-CR law gate passed |
| M4d | [log](M4d-log.md) | Speed gate: run 1 failed (W5 1.28), the rerun passed (all ≤ 1.25); law gate of `auto` passed; worker rerun failed on `binding12` `reuse` (1.55), recorded as a known limitation |
| M5a | [log](M5a-log.md) | `FlatView` and FSP on rules models; DSMTS and 0.6.0 FSP references reproduced |
| M5b | [log](M5b-log.md) | Analysis on the flat view (`analysis/1.0`) |
| M5c | [log](M5c-log.md) | Tau gates on converted models through the rules `tau` method |
| M5d | [log](M5d-log.md) | Archives of rules runs (`archive/1.0`) |
| M5e | [log](M5e-log.md) | ODE right-hand sides and Jacobians on the flat view (dual numbers) |
| M5f | [log](M5f-log.md) | Importers return `model/1.0`; GillesPy2 export on the flat view |
| M5g | [log](M5g-log.md) | Workers and jobs on `bind_parameters`; jobs at problem version 3 |
| M6a | [log](M6a-log.md) | One Python API: the rules API at the top level, `crn_craft.rules` kept |
| M6b | [log](M6b-log.md) | `Result.save` to markdown and JSON (FR1b) |
| M7 | [log](M7-log.md) | One command line (`run`, `check`, `convert`, `describe/1.0`) |
| M8a | [log](M8a-log.md) | Validation scripts on the unified engine; law gates and the mutation gate (47 of 47) pass |
| M8b | [log](M8b-log.md) | The network engine, its contracts and the old API deleted; law gate after the deletion identical to the M4d gate |
| M8c | [log](M8c-log.md) | Generated loops on the single crate; fixtures differ only in paths and fingerprints |
| M9 | [log](M9-log.md) | Tutorials, examples and integrations on the unified API |
| M10a | [log](M10a-log.md) | Tutorial rewritten; every block checked |
| M10b | [log](M10b-log.md) | Manual and references updated; the 0.6.0 manual frozen and run in CI (I2b); link checker |
| M10c | [log](M10c-log.md) | Architecture, agent and report documents; GillesPy2 re-measured; mtDNA workload 0.6.0 vs 0.7.0 |
| M10d | this page | README and CHANGELOG (written by the coordinator) |
| Review | [log](review-log.md) | Final review of the whole merge before the pull request |

## Deviations from the plan

- **M4f** was added before the speed gate, because the unified engine was
  1.6–5.4× slower than 0.6.0 after M4c.
- **M5h** (tau leaping speed on large networks) was found in M5c and
  deferred to the follow-up optimization release; the CHANGELOG lists it as a
  known limitation. So is the `binding12` `reuse` miss of M4d.
- **`crn_craft.rules.ensemble`** writes `"selector": "linear"` when it
  upgrades a request to `run/1.0`, so adding quantiles never changes the
  paths (coordinator, after M6a).
- **`crn.compile_loop`** at the top level builds for the `model/1.0`
  document that `crn.ensemble` runs (coordinator, found in M10b).
- **Agents and restarts.** The work ran as background agents, at most two
  at a time. Five server restarts and two API usage limits interrupted
  them; every long run was resumable, and no result was lost.
