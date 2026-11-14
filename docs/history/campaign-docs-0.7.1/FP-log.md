# FP log: the mtDNA case study runs the fast path

Date: 13 November 2026. Branch `docs/FP` (from `docs-0.7.1` at `c436e1a`),
worktree `/home/claude/docs071/wt-FP`. Nothing pushed. Task: since E1,
`mtdna_direct` certifies the builder-built three-class model under
`run/1.0`; the owner decided that the case study (pages 12 to 12e) uses
the fast path as its default method wherever the model and the feature
allow.

## Setup

- Venv `/home/claude/m07/venv-m10a`: `crn_craft.__version__` is `0.7.0`
  (the version bump comes with release prep), and `crn.ensemble(three_class(),
  ..., method="mtdna_direct")` reports
  `chacha8-u53-ziggurat/rules-mtdna-direct/1` with no diagnostics. The
  0.7.0 engine refuses the builder form (`CRN-R051`), so this is the E1
  engine; no wheel was built.
- The `bash` blocks ran the console script of that wheel (the worktree has
  no `target/release`). Compiled loops built against this worktree
  (`CRN_CRAFT_SOURCE`).
- Every run under the CPU lock.

## 1. Which blocks accept the fast path

Each feature tested with `method="mtdna_direct"` on the block's own model
and arguments (scripts in `/home/claude/docs071/fp-work/scripts/`, not
committed). "Yes" means the result's sampler is
`chacha8-u53-ziggurat/rules-mtdna-direct/1` and `metadata["diagnostics"]`
is empty; trackers, logs and callbacks were also checked to produce data
(appearances, first passages, log entries, callback events equal to the
event count, callback stops).

| Page, block | Uses | Fast path | Reason |
| --- | --- | --- | --- |
| 12, A hundred cells | plain model, two report times | yes | |
| 12a, A first ensemble | takeover stop 0.9, 13 report times, records | yes | |
| 12a, The mean hides the typical cell | quantiles, `keep_paths="all"` | yes | |
| 12a, Takeover fraction | takeover stop 0.9 | yes | |
| 12a, Variants | `K=500`, `w0=400`, `m0=100`, half-life 20 days, stop at 0.6 | yes | |
| 12a, Variants (misspelled name) | error before any run | n/a | no `method`; the same `CRN-R012` also with `mtdna_direct` |
| 12b, Class appearances | three trackers, `keep_episodes`, `first_passage` predicate, stop 0.9 | yes | |
| 12b, Mutation events in a log | `event_log` with two rules | yes | |
| 12b, Stop at the first advantaged mutant | callback that stops a cell, `threads=1`, `rule_names` | yes | |
| 12c, ATP as a continuous variable | `add_ode`, `hybrid` | no | not the three-class structure: `mtdna_direct` raises `CRN-R020: continuous variables require the hybrid method` |
| 12d, The same cells under every method | eight methods | unchanged | the page's subject; it already has an `mtdna_direct` row |
| 12d, The fast path | Static Direct against `mtdna_direct` | unchanged | the comparison is the block |
| 12d, A compiled model | `compile_loop`, trackers | no | `run.simulate(compiled=loop)` with `mtdna_direct` raises `CRN-R020: generated rules loops require Static Direct`; `compile_loop` has no method argument |
| 12e, Parameter sweeps | `bind_parameters`, stop 0.9 | yes | also checked at `sm = 1.0` and `mu = 1e-4, sm = 2.0` |
| 12e, From the command line | request file `run/1.0`, `"method": {"kind": "mtdna_direct"}`, tracker on `A`, quantile 0.5, `--archive`, `verify` | yes | `check` valid with no diagnostics; archived result names the fast-path sampler |

## 2. Changes, one commit per page

| Commit | Page | Change |
| --- | --- | --- |
| `9d3ae8b` | `12-mtdna.md` | block on `mtdna_direct`; Conventions states the rule once (fast path; features it runs with; where the pages use another method and why; 0.97 and 3.6 of 13 November; `CRN-R051` fallback outside the structure). Output: sampler, mean mutant fraction 0.1108 to 0.1111 |
| `22d2fa8` | `12a-mtdna-takeover.md` | four blocks on `mtdna_direct`; outputs rerun (23 of 200 taken over, was 15; 28,536,814 events, the number E1 reports for arm A and G). Prose has no number from the outputs; the claim "the same seed gives the same cells" holds (0.115 = 23/200, first day 273 = 272.7) |
| `cb74cf2` | `12b-mtdna-events.md` | three blocks on `mtdna_direct`; outputs rerun; the closing paragraph now says "the fast path" and day 739.2 (the log's advantaged mutation of cell 0, equal to the callback's first stop day, as before) |
| `951116c` | `12c-mtdna-atp.md` | one sentence: the fast path does not apply, the model is not the three-class structure; output unchanged |
| `dcf451d` | `12d-mtdna-methods.md` | blocks and outputs unchanged; says the other pages of the case study run the fast path; the compiled block says compiled loops run Static Direct; "Matching an hourly scheme" runs the exact model on the fast path |
| `9c8c664` | `12e-mtdna-large-studies.md` | sweep and request file on `mtdna_direct`; outputs rerun (sweep table; 12 candidates in 6.1 s, was 63.1 s under the 0.7.0 engine; archive: 20 of 200 stopped, fast-path sampler); prose: every candidate runs the fast path, Static Direct about 3.6 times as long, a compiled loop serves one parameter set |

Headings and "Applies:" lines are unchanged, so every anchor stays.
Block counts are unchanged, so `validation/guide-blocks.json` is unchanged.

## 3. Checks

See the report below.

## Report

Runner `/home/claude/docs071/fp-work/checks/run.sh` (logs beside it), each
`check_guide.py` run under the CPU lock, wheel of `venv-m10a` (E1 engine),
compiled loops built against this worktree:

| Check | Result |
| --- | --- |
| `check_guide.py --require-cargo` on the six pages | passed: 6 files, 22 runnable blocks |
| the same under `python -O` | passed; report byte-identical (`cmp`) |
| `check_guide.py --require-cargo`, all pages | passed: 26 files, 109 runnable blocks |
| `check_links.py` | 949 links, 36 stubs, 53 SUMMARY entries, 0 broken |

Deviations:

1. No `--binary`: the worktree has no release CLI, so `bash` blocks ran
   the wheel's console script (same engine, E1).
2. The 12e sweep time fell from 63.1 s to 6.1 s. Part of it is the fast
   path and part is E1's faster Static Direct; the page quotes only the
   13 November ratio 3.6.

Open questions:

1. Pages outside this task were not changed. `13-fitting.md` and
   `14a-mtdna-model-files.md` already run `mtdna_direct`. No other page
   says which method the case study uses.
