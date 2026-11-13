# D6 log: the reference pages

Date: 13 November 2026. Branch `docs/D6`, from `docs-0.7.1`. Task: the
reference pages of the documentation restructure (the table of contents,
sections 1.4, 4 and 6).

## What was produced

| Path | What |
| --- | --- |
| `validation/make_python_reference.py` | Generates `docs/reference/python-api.md` from the installed package. `--write` writes, `--check` fails on a difference and never writes, no flag prints. `--lint` prints docstring problems. |
| `docs/reference/python-api.md` | Generated page: 6 task groups, one entry per public name (signature, first paragraph, argument table, rest of the docstring). |
| `.github/workflows/ci.yml` | `python validation/make_python_reference.py --check` after the link check. |
| `docs/reference/methods.md` | Base `direct-ssa.md` plus the method material of `mtdna-rules.md`, the hybrid method and the numerical policy of the model reference, and the manual's certificate of `mtdna_direct` and compiled-loop requirements. |
| `docs/reference/contracts.md` | Adds the version tables of the request and the model (manual section 15, `mtdna-rules.md` 87-110) and the request, its budgets and admission (model reference 377-416). |
| `docs/reference/results.md` | Adds the result metadata table (manual section 18), the accumulation notes (model reference 418-432), event logs, tracker columns, callbacks. |
| `docs/reference/model-language.md` | Hybrid, request and result, numerical policy and the diagnostics table moved out; the K1 validation section cut to a pointer; the history section "Prepared Python ensembles" dropped. |
| `docs/reference/diagnostics.md` | Gets the `CRN-R` table of the model reference. |
| `docs/developer/rust-api.md` | New: manual section 16, the Rust part of the 0.7.0 migration, two complete programs. |
| `validation/check_rust_docs.py` | Compiles and runs the `rust` blocks of `rust-api.md` in a scratch crate. |
| `docs/mtdna-rules.md`, `docs/reusable-parameters.md`, `docs/alternative-ssa.md` | Stubs. |
| `docs/evidence/mtdna-performance-analysis.md` | "The Java comparison" appended to Reproducing. |

## Where each part went

`mtdna-rules.md` (373 lines), by the table in section 4.2:

| Lines | Destination |
| --- | --- |
| 1-16 | no new content against the manual; the model description is in guide pages 12 and 14a (another task) |
| 20-38 | dropped (an unchecked block; 12b has the checked version) |
| 40-62 | dropped (duplicate of chapter 4 and Migrating) |
| 64-83 | trackers, logs, callbacks: `reference/results.md` (and chapter 6) |
| 87-110 | `reference/contracts.md`, "Versions of the request and the model" |
| 112-147 | `reference/methods.md`: methods table, higher-order rules, whole-rule propensities, static lowering |
| 151-179 | `reference/methods.md`, "The fast path" |
| 181-204 | dropped (duplicate of chapter 5; unchecked block) |
| 206-235 | `reference/methods.md`: tau, binomial, next reaction |
| 237-258, 292-295 | `reference/methods.md`, "The ATP block and thinning" |
| 260-277, 297-303 | `evidence/mtdna-performance-analysis.md`, "The Java comparison" |
| 279-290 | `reference/results.md`, "Event logs and tracker columns" |
| 305-349 | `reference/methods.md`, "Waiting times and the options of Static Direct" |
| 351-373 | dropped (the performance analysis has the measurements) |

`reusable-parameters.md`: 16-48 are in guide page 5 (parameter sweeps, common
random numbers, limits); 1-14 (the 0.6.0 table) is in `guide/migrating.md`.
`alternative-ssa.md`: the table is in `guide/migrating.md`, which another task
writes; I checked the copy in that task's checkout.

From the manual: certificate 987-1033 and compiled loops 1285-1347 to
`methods.md`; the version tables to `contracts.md`; the metadata and
terminations to `results.md`; the argument and result tables of section 18
to `python-api.md` (the argument descriptions are the table `ARGUMENTS` in the
script, the result attributes are the docstring of `Result`); section 16 to
`rust-api.md`. Frequent diagnostics of the manual are in `diagnostics.md`.

From `rules-language.md` (now `model-language.md`): hybrid 306-375 and numerical
policy 434-458 to `methods.md`; request 377-416 to `contracts.md`; result
418-432 to `results.md`; diagnostics 460-478 to `diagnostics.md`; K1 480-486
replaced by a pointer; 527-542 dropped.

## Checks

See the report at the end of this file.

## Report

Checks, run on 13 November 2026 with the 0.7.0 wheel of
`/home/claude/m07/venv-m10a`:

- `python validation/make_python_reference.py --check`: up to date. The
  script fails when a public name is in no group, when a group lists a name
  that does not exist, and when an argument has no description.
- `ruff check validation` and `ruff format --check validation`: clean.
- `python validation/check_links.py`: 625 links, 34 stubs, 24 broken. All 24
  point to pages that other tasks write: guide pages `12*`, `13-fitting`,
  `14a`, `15-from-r`, `migrating`, `../howto.md` and `../comparison.md`. This
  includes the stub `alternative-ssa.md`, which points to
  `guide/migrating.md`, and the links of `reference/methods.md` to page 14a.
  Nothing else is broken.
- `validation/check_rust_docs.py` under the CPU lock: 2 rust blocks, 0 failed
  (both compile and run against the workspace crate `crn-craft` 0.7.0).
- `json` blocks of the reference pages are strict JSON without duplicate keys.
  The one model (the first model of the model reference) passes
  `crn-craft check` (0.7.0 binary of the other checkout, read only). The other
  blocks are fragments (method objects, the `continuous` field), not models.
- No reference page has a `python` or `bash` block. The shell commands of
  the GillesPy2 export page and the first-model page are `text` blocks.

Docstring problems (package sources not edited):

- `PreparedEnsemble` has no docstring. The page uses its own sentence.
- `crn_craft.importers.ImportFailure` has no docstring of its own (it shows
  the one of `ValueError`). The page uses its own sentence.
- `Model.add_reaction` refers to "the module documentation" for the counting
  of identical reactants. The module docstring of `crn_craft` does not have it.
- `crn_craft.ensemble`, `prepare` and `loop_source` refer to
  `crn_craft.rules.ensemble`, `rules.fsp` and `rules.loop_source` for the
  arguments, which are not public names of `crn_craft`. The argument
  descriptions on the page are the script's own table.
- `crn_craft.rules.ensemble` says `mtdna_direct` is about 3.3 and 3.9 times
  faster than Static Direct; the manual (12 November measurement) says 2.8 and
  3.1. The docstring is stale.
- `import_sbml` and `import_bng_net` do not say what they return (the
  `crn-craft.import/1.0` report; the module docstring does).
- `Result` documents `TERMINATIONS`, which is not a public name of `crn_craft`
  (it is `crn_craft.rules.TERMINATIONS`).

Not changed, outside this task: crate and script sources still name
`docs/rules-language.md`, `docs/mtdna-rules.md` and the other old paths
(`crates/crn-craft/src/lib.rs`, `rules/model.rs`, `rules/atp.rs`, `unified.rs`,
`validation/package_*.py`, `render_acceleration_report.py`). The stubs keep
the URLs working. The `check_rust_docs.py` script is not in CI.

Deviations: stub of `mtdna-rules.md` points to `reference/methods.md` only
(the stub format allows one link; the table of contents also names guide page
12). Links in `README.md`, `CHANGELOG.md` (link targets only), `AGENTS.md`,
`docs/tutorial.md`, `docs/mtdna-user-manual.md` and others were changed from
the three stubs to their targets, as `check_links.py` requires.
