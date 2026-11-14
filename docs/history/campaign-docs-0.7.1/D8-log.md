# D8 log: one checker for the guide; the two old guides retired

Date: 13 November 2026. Branch `docs-0.7.1`, main checkout, on top of D5
(`3d80c6d`). Nothing pushed. Specification: sections 1.7, 2.4 and 9 of the
table of contents (D1, second version), the D3 and D4 logs, the D3 chapter
driver.

## Commits

| Commit | Content |
| --- | --- |
| `18ae89a` | `check_guide.py`, `guide-blocks.json`, CI, `check_rust_docs.py` on the crate READMEs, stubs, link and script updates; `check_tutorial.py` removed |
| this commit | this log |

## What was produced

- **`validation/check_guide.py`**, from the D3 driver (with D4's smoke mode
  and `file` binding). Default set: every `docs/guide/*.md`, `README.md`,
  `python/README.md` and `docs/index.md` (26 files). The eleven points of
  section 1.7:

  | Point | How |
  | --- | --- |
  | 1. One file, one session | each file in a child process and a fresh, empty scratch directory |
  | 2. Setup block | the first `python` block of a file |
  | 3. Section-local state | at every level-2 heading the namespace is replaced and the setup block runs again; the scratch directory is kept for the whole file (chapters 10, 11 and 14 pass files between sections), as the D3 driver does |
  | 4. Files made in sight | `<!-- fixtures: ... -->` copies repository files under their paths; nothing else is in the directory |
  | 5. Fence languages | `python`, `bash` run; `text` must follow a runnable block; `json` strict (duplicate keys, NaN) and `crn-craft check` for models and requests; `sh` `bash -n`; `r` only in a `check: none` page; any other language, or none, is an error. A `check: none` page still syntax-checks its `sh` and parses its `json` blocks |
  | 6. Output comparison | exact after stripping; timings (`<number> s`) normalized on lines that hold one |
  | 7. Standalone pages | `<!-- check: standalone -->`: each block in a fresh interpreter from the repository root; `<!-- file: PATH -->` binds a block to a file byte for byte; the README's example and output equal `docs/index.md`'s, byte for byte |
  | 8. Long pages | `<!-- check: long -->`: smoke mode (`CRN_CRAFT_SMOKE=1`, output not compared) unless `--long` |
  | 9. Cargo | `compile_loop` blocks: skipped to the end of the file without cargo; `--require-cargo` fails instead |
  | 10. Count | `validation/guide-blocks.json`, keyed by repository path; a page with fewer runnable blocks fails, a checked page without an entry fails, an entry without a page fails (default set). Skipped-for-cargo blocks count as present |
  | 11. Options | `--binary`, `--max-replicates` (caps `crn.ensemble` and `crn.rules.ensemble`, outputs not compared; standalone blocks through a launcher), `--require-cargo`, `--output`; the report has no timings and is identical under `python -O` |

  Also kept for maintainers: `--write` (regenerate outputs), `--long`,
  `--keep`, `--counts`, `--root`. Each report row names its level-2 section.
- **`validation/guide-blocks.json`**: 108 runnable blocks (73 in D3's
  pages, 31 in D4's, 4 in the READMEs and `docs/index.md`); the numbers
  are the blocks each page has today.
- **`validation/check_tutorial.py` removed.** What it enforced and where it
  is now: exact output comparison (point 6; only timing lines are
  normalized, and no tutorial block printed one), strict JSON and
  `crn-craft check` (point 5), `bash -n` on `sh` (point 5), unknown fence
  languages (point 5), at least 15 runnable blocks (replaced by the
  per-page counts), the native CLI first on PATH (`--binary`).
- **`validation/check_manual.py`** stays for the frozen manual. Its default
  `--manual` is now `validation/compat/manual-0.6.0.md` (was the manual
  that is now a stub, which would have passed with no blocks) and its
  docstring says so; the code is otherwise unchanged. The CI step for the
  frozen manual is unchanged.
- **`validation/check_rust_docs.py`** checks `developer/rust-api.md` and
  the three crate READMEs: complete `rust` blocks compile and run;
  `<!-- check: fragment -->` before a fence marks a fragment, which is
  listed and not compiled (the `crn-craft-rules` README block); `toml`
  blocks parse and their `crn-craft` version matches the workspace version;
  `sh` blocks pass `bash -n`; other languages fail; at least one program
  must compile.
- **CI** (`.github/workflows/ci.yml`): in the `python` job the steps of
  `check_manual.py` on the 0.7.0 manual and on `python/README.md` and of
  `check_tutorial.py` are replaced by one step that runs `check_guide.py`
  normally and under `-O` and diffs the reports (`CRN_CRAFT_SOURCE` and
  `CARGO_TARGET_DIR` as for the compiled-loop step); pip installs
  `matplotlib==3.11.2 fcmaes==2.1.0 scipy==1.18.1 loguru==0.8.0` (the
  versions of the campaign venv); timeout 40 → 60 minutes (one run takes
  about 4 minutes here with warm compiled loops; CI builds two loops cold).
  The `contracts` job runs `check_rust_docs.py`. Every job still runs real
  checks; the frozen-manual steps are untouched.
- **Stubs.** `docs/tutorial.md` and `docs/mtdna-user-manual.md` are short
  pages with the table of section 3 (old section → guide section, linked
  with anchors), a link to the version at tag `v0.7.0`, and for the manual
  the case study, page 14a, Migrating and the frozen 0.6.0 copy. They are
  not one-sentence stubs in the sense of `check_links.py`, so their links
  are checked like any page's.
- **Links to the old manual retargeted**: `developer/architecture.md`
  (→ Migrating), `developer/specialized-kernels.md` (→ chapter 5, method
  reference), `evidence/mtdna-performance-analysis.md` three links (→ case
  study, 14a, method reference). `examples/mtdna/README.md`: says that the
  numbered scripts are the 0.7.0 manual's complete programs on
  `crn_craft.rules` and the model files, and maps each to its guide page;
  the sentence on `check_manual.py` is gone. `developer/extending.md`: the
  check list names `check_guide.py` and `check_rust_docs.py`, with a
  paragraph on `guide-blocks.json`, `--write` and fragments.
  `pyproject.toml`: `Documentation` URL → the site.
- **`examples/tutorial/workdir.py`**: `prepare_workdir` (used only by the
  removed checker) is gone. `make_assets.py` was run into a scratch folder:
  the six figures and `numbers.json` are byte-identical to the committed
  ones, and the twelve terminal transcripts were written.
- `validation/mutation_check.py` names neither checker; no other script
  under `validation/`, `scripts/` or the workflows names `check_tutorial.py`.

## Gate G1

Every block that `check_tutorial.py` and `check_manual.py` ran on `main` at
0.7.0 (33 tutorial, 26 manual, 2 `python/README.md`) has a destination that
`check_guide.py` ran in the final run. The table is produced by
`/home/claude/docs071/d8/g1.py` (outside the repository) from the final
report: for each ledger destination it requires the page to pass and every
block of the named language in the named level-2 section to have status
`ran`. M9, M11 and M24 have two destinations each (section 2.4). 61 of 61
ids, 0 missing, 0 not run. No blocker.

| Ledger id | File | Section | Runs |
| --- | --- | --- | --- |
| T1 | `02-first-model.md` | Species, parameters and reactions | runs: yes (python) |
| T2 | `14-model-files.md` | What the builder writes | runs: yes (python) |
| T3 | `14-model-files.md` | Packaged examples | runs: yes (python) |
| T4 | `02-first-model.md` | Checking a model | runs: yes (python) |
| T5 | `10-command-line.md` | A request file | runs: yes (python) |
| T6 | `02-first-model.md` | Checking a model | runs: yes (python) |
| T7 | `03-running.md` | Many trajectories in one call | runs: yes (python) |
| T8 | `10-command-line.md` | Python and the command line give the same result | runs: yes (python) |
| T9 | `04-results.md` | What a result holds | runs: yes (python) |
| T10 | `04-results.md` | Quantiles and the median | runs: yes (python) |
| T11 | `04-results.md` | Save a run in one file | runs: yes (python) |
| T12 | `11-analysis-import-export.md` | Oscillations | runs: yes (python) |
| T13 | `07-methods.md` | Exact distributions with finite state projection | runs: yes (python) |
| T14 | `07-methods.md` | Approximate leaping | runs: yes (python) |
| T15 | `03-running.md` | Budgets | runs: yes (python) |
| T16 | `03-running.md` | Stopping a run | runs: yes (python) |
| T17 | `08-reproducibility.md` | What makes a run reproducible | runs: yes (python) |
| T18 | `11-analysis-import-export.md` | Conservation laws and structure | runs: yes (python) |
| T19 | `11-analysis-import-export.md` | First passage from saved paths | runs: yes (python) |
| T20 | `11-analysis-import-export.md` | Convert a network model of release 0.6.0 | runs: yes (python) |
| T21 | `11-analysis-import-export.md` | Import SBML and BioNetGen | runs: yes (python) |
| T22 | `07-methods.md` | Static Direct and its selectors | runs: yes (python) |
| B1 | `10-command-line.md` | Check a model and a request | runs: yes (bash) |
| B2 | `10-command-line.md` | Check a model and a request | runs: yes (bash) |
| B3 | `10-command-line.md` | Run, archive and verify | runs: yes (bash) |
| B4 | `10-command-line.md` | Python and the command line give the same result | runs: yes (bash) |
| B5 | `04-results.md` | Save a run in one file | runs: yes (python) |
| B6 | `10-command-line.md` | Requests for fsp and tau | runs: yes (bash) |
| B7 | `10-command-line.md` | Requests for fsp and tau | runs: yes (bash) |
| B8 | `10-command-line.md` | Stop a run with a cancel file | runs: yes (bash) |
| B9 | `10-command-line.md` | Replay | runs: yes (bash) |
| B10 | `11-analysis-import-export.md` | Convert a network model of release 0.6.0 | runs: yes (bash) |
| B11 | `11-analysis-import-export.md` | Import SBML and BioNetGen | runs: yes (bash) |
| M1 | `12a-mtdna-takeover.md` | A first ensemble | runs: yes (python) |
| M2 | `12a-mtdna-takeover.md` | The mean hides the typical cell | runs: yes (python) |
| M3 | `04-results.md` | Save a run in one file | runs: yes (python) |
| M4 | `03-running.md` | Seeds and threads | runs: yes (python) |
| M5 | `12a-mtdna-takeover.md` | Takeover fraction and time to threshold | runs: yes (python) |
| M6 | `12b-mtdna-events.md` | Class appearances and extinctions | runs: yes (python) |
| M7 | `12b-mtdna-events.md` | Mutation events in a log | runs: yes (python) |
| M8 | `12b-mtdna-events.md` | Stop at the first advantaged mutant | runs: yes (python) |
| M9 | `12d-mtdna-methods.md` | The same cells under every method | runs: yes (python) |
| M9 | `14a-mtdna-model-files.md` | The fast path | runs: yes (python) |
| M10 | `07-methods.md` | Static Direct and its selectors | runs: yes (python) |
| M11 | `14-model-files.md` | A model document for three strains | runs: yes (python) |
| M11 | `14-model-files.md` | Run it | runs: yes (python) |
| M12 | `12d-mtdna-methods.md` | A compiled model | runs: yes (python) |
| M13 | `05-speed.md` | Compiled models | runs: yes (python) |
| M14 | `14a-mtdna-model-files.md` | Every lineage with the clone model | runs: yes (python) |
| M15 | `12c-mtdna-atp.md` | ATP as a continuous variable | runs: yes (python) |
| M16 | `09-continuous-variables.md` | A variable between events | runs: yes (python) |
| M17 | `09-continuous-variables.md` | Independent variables | runs: yes (python) |
| M18 | `14a-mtdna-model-files.md` | The stored ATP block and thinning | runs: yes (python) |
| M19 | `12e-mtdna-large-studies.md` | Parameter sweeps | runs: yes (python) |
| M20 | `02-first-model.md` | A second model with two reactants | runs: yes (python) |
| M21 | `06-stops-trackers-logs.md` | Stop conditions | runs: yes (python) |
| M22 | `12-mtdna.md` | The model | runs: yes (python) |
| M23 | `07-methods.md` | Static Direct and its selectors | runs: yes (python) |
| M24 | `12a-mtdna-takeover.md` | Variants | runs: yes (python) |
| M24 | `14a-mtdna-model-files.md` | The fast path | runs: yes (python) |
| M25 | `03-running.md` | Budgets | runs: yes (python) |
| M26 | `12e-mtdna-large-studies.md` | From the command line | runs: yes (python) |
| R1 | `python/README.md` | Example | runs: yes (python) |
| R2 | `python/README.md` | Example | runs: yes (python) |

File names are under `docs/guide/` unless a path is given. The section is
the level-2 section of the block (M23 is under the level-3 heading "Larger
models"; M16's model is the setup block of chapter 9).

## Checks run

All on the 0.7.0 wheel of this checkout in `/home/claude/m07/venv-m10a`,
with `target/release/crn-craft`, each heavy run under the CPU lock.

| Check | Result |
| --- | --- |
| `check_guide.py --binary target/release/crn-craft --require-cargo` | passed: 26 files, 108 runnable blocks ran, 72 outputs compared, 6 `sh` blocks syntax-checked, 3 fitting blocks in smoke mode (not compared), 4 `r` blocks not run; README and landing-page example equal; 3 min 49 s |
| the same under `python -O` | passed; report byte-identical (SHA-256 `93686df7…`) |
| `check_guide.py --max-replicates 8` on four pages | passed; 11 outputs "not compared: capped"; 12a in 1.2 s instead of 17.5 s |
| `check_guide.py` on chapter 5 without cargo on PATH | the compiled-models block and its output listed as skipped, page passes; with `--require-cargo` it fails |
| `check_manual.py --manual validation/compat/manual-0.6.0.md --binary ... --require-cargo`, normal and `-O` | 22 of 22 passed; reports identical; the frozen file is unchanged |
| `check_rust_docs.py` | 2 rust programs compiled and ran (`rust-api.md`); 1 fragment; 2 `toml`, 1 `sh` checked |
| `check_links.py` | 934 links, 36 stubs, 53 SUMMARY entries, 0 broken |
| `mdbook build` and `check_links.py --site` | no warnings; 3,524 references, 0 broken |
| `check_installed.py` from `/tmp` with `python -I`, and with `-I -O` | 18 checks passed; reports identical |
| `ruff check validation python examples`; `ruff format --check validation` | passed |
| `cargo fmt --all -- --check` | passed |
| `make_assets.py` into a scratch folder | figures and `numbers.json` byte-identical |

Negative tests (`/home/claude/docs071/d8/negative.sh`, outside the
repository), each on a scratch page, each fails with the expected message:
a wrong output, a name used from another section, a `toml` fence, a fence
without a language, an `r` block outside a `check: none` page, a `text`
block before any runnable block, a `json` block with a duplicate key, with
`NaN`, a model document that fails `crn-craft check`, an `sh` syntax error,
a failing `bash` block, fewer blocks than counted, an unknown `check`
option, a `file` binding that differs, a page without an entry in the
counts file, and a README example that differs from the landing page's. A
control page passes. For `check_rust_docs.py`: an unmarked fragment fails
to compile, and `crn-craft = "0.6"` fails the version check.

## Deviations

1. `check_manual.py` is not byte-unchanged: its default `--manual` now
   points at the frozen copy, and its docstring says that this is its one
   use. The prompt's "the file is unchanged" was read as the frozen manual
   (G2), which is unchanged; its CI step is unchanged.
2. The report schema is `crn-craft.guide-check/1.0` and the counts file is
   keyed by repository path, not file name (the driver used file names;
   `README.md` occurs twice).
3. Block counts count `skipped: needs cargo` blocks as present, so a local
   run without cargo still passes the count; CI uses `--require-cargo`.
4. The python job's timeout is raised to 60 minutes as a margin for the
   two cold compiled-loop builds and the fitting smoke run.
5. `pyproject.toml`'s `Documentation` URL now names the site, which exists
   only after GitHub Pages is enabled.

## Open questions

1. `ruff format --check examples` reports `examples/mtdna/fit_takeover.py`
   and `fit_heteroplasmy.py` (from D4; already so before this task). CI
   formats only `validation`. Formatting them means changing the two blocks
   of chapter 13 that are bound to the files; their outputs would not
   change. Fix in D9 or the release preparation?
2. The crate doc comments `crates/crn-craft/src/lib.rs` (line 122) and
   `src/rules/mod.rs` (line 198) still link `docs/tutorial.md` and
   `docs/mtdna-user-manual.md` on GitHub. Both now lead to the stubs'
   tables; changing them touches Rust sources in a docs-only release.
3. The numbered scripts in `examples/mtdna/` (the 0.7.0 manual's complete
   programs) are not run by any checker, as before this task.
4. `docs/campaign/README.md` still links `../tutorial.md` (D5's note); it
   now reaches the stub table, so it works.
