# D4 log: the case study, fitting, R and migrating

Date: 13 November 2026. Branch `docs-0.7.1`, on top of D2, D3 and DR
(`108ee0d`). Nothing pushed. Specification: the table of contents of task
D1, second version (`docs-toc.md` of the campaign, sections 1.5, 1.7, 2.3,
2.4, 4.3, 4.4, 4.7, 4.11 and 9), and the rules of the D3 prompt.

## What was written

| Page | Lines | Budget | Difference | Runnable blocks | Count in the table of contents |
| --- | ---: | ---: | ---: | ---: | ---: |
| `12-mtdna.md` | 150 | 120 | +25% | 2 | 1 |
| `12a-mtdna-takeover.md` | 212 | 230 | −8% | 5 | 5 |
| `12b-mtdna-events.md` | 186 | 190 | −2% | 4 | 4 |
| `12c-mtdna-atp.md` | 123 | 120 | +3% | 2 | 2 |
| `12d-mtdna-methods.md` | 155 | 150 | +3% | 3 | 3 |
| `12e-mtdna-large-studies.md` | 186 | 160 | +16% | 5 | 4 |
| `13-fitting.md` | 588 | 400 | +47% | 3 (standalone, long) | 3 |
| `14a-mtdna-model-files.md` | 250 | 240 | +4% | 5 | 5 |
| `15-from-r.md` | 150 | 150 | 0% | 0 (`check: none`, 4 `r` blocks) | 0 |
| `migrating.md` | 201 | 150 | +34% | 2 | 2 |

Over the margin of 15 percent:

- `12-mtdna.md`: the setup block, the function `three_class`, is 25 lines
  where the table of contents planned 16 (the parameter dictionary, three
  degradation reactions in a loop, the stop), and every section has its
  "Applies" line. Prose is short; the table of pages and questions stays,
  because it is the navigation of the case study.
- `12e`: one more block than planned. With `--archive`, `crn-craft run`
  prints a check of the archive, not the result, so a short Python block
  reads the three JSON outputs and the archived result.
- `13-fitting.md`: the three scripts are 299 lines (each repeats
  `three_class`, 22 lines, and has a three-line smoke switch), and their
  outputs 25; the block budget was 290. Prose follows the old tutorial's
  sections, with sections 7 and 8 cut to six points each as specified.
  Blocks are complete programs bound to their files and cannot be cut.
- `migrating.md`: the five migration tables of the sources (the changelog's
  Python and command-line tables, the manual's tables for 0.5.0, 0.4.0 and
  0.2.0) are 95 lines of the 201. The table of contents assigns all of them
  to this page; none was cut.

Also: `examples/mtdna/fit_takeover.py`, `fit_heteroplasmy.py` and
`noise_and_parallel.py` rewritten to the builder model and the top-level
API; `examples/r/sir.R` on the top-level module; `examples/mtdna/README.md`
(the fitting rows, one sentence); permanent stubs at
`docs/mtdna-optimization-tutorial.md` and `docs/using-from-r.md`, with the
links in `README.md`, the manual, `mtdna-rules.md` and the performance
analysis pointing to the guide; a declared timing of the builder form
(below).

## Commits

| Commit | Content |
| --- | --- |
| `7310c60` | declaration: timing of the builder form of the mtDNA model |
| `5eb42d5` | results of that timing (`validation/merge/manual-ratios/builder-form.json`) |
| `b23b403` | the three fitting scripts on the builder model; `ruff format` of the timing script |
| `6cf43e0` | chapter 12, the index page |
| `5c4c42c` | 12a |
| `ff80833` | 12b |
| `1cc82f4` | 12c |
| `892fce7` | 12d |
| `0ce391c` | Migrating |
| `51627ed` | chapter 15 and `examples/r/sir.R` |
| `b30d144` | the two stubs and the links to chapters 13 and 15; one sentence of chapter 15 |
| `c0af7a7` | 12e |
| `a1bf8d6` | 14a |
| `44579e6` | chapter 13 |
| this commit | this log; final verification |

## The builder form is slower than the table of contents assumed

The table of contents and the coordinator's answer to question 1 assumed
that the builder model under Static Direct costs what the model file under
Static Direct costs, 2.8 times the fast path. It costs more. A first run of
the index page's model took 5.3 s for the 200 cells that the fast path runs
in 0.5 s. Before quoting a number, I declared a measurement in
`validation/merge/RESULTS.md` ("Docs 0.7.1: the builder form of the mtDNA
model", `7310c60`) and ran it under the CPU lock (`5eb42d5`). Medians of
three on the workload of the manual's first example, 13 November 2026:

| Arm | Median s | Ratio |
| --- | ---: | --- |
| A model file, `mtdna_direct` | 0.543 | 1 |
| B model file, `static_direct` | 1.497 | B / A 2.76 |
| E builder model, `static_direct` | 5.373 | E / A 9.89; E / B 3.59 |
| F E with the compiled loop | 0.759 | E / F 7.08; F / A 1.40 |

The builder writes whole-rule propensities, which the interpreter evaluates
as expressions at every event (eight per event here, since every event
changes every rate); the file form uses per-instance rates. Writing the
degradation as numeric constants, or simplifying the expressions, changed
little (5.0 against 5.5 s in a first look). The compiled loop removes most
of the cost, so page 12d shows a gain of 7.1, not 2.1.

Consequences, stated in the pages: chapter 12 (Conventions) quotes 9.9 and
3.6; 12d quotes 7.1 and 1.4; 14a quotes 2.8 (12 November), 9.9 and 1.4;
chapter 13 says what one evaluation costs and that the file under the fast
path would take about a tenth. Chapter 5 (D3) still quotes the compiled
gain of 2.1 for the file form, which is correct as stated.

**The fitting run.** One evaluation of 800 cells took about 21 s instead of
2.2 s, so the full run of the three scripts took 2 h 12 min (2,375 s,
4,473 s and 1,071 s) instead of about 12 minutes. It ran once, under the
CPU lock, one lock hold per script. `fit_heteroplasmy.py` was stopped twice
(once by me, to restart the runner within the two-hour limit of a
background task, once by a restart of the environment) and run again from
its start; `fit_takeover.py`'s output was kept, after checking that the
SHA-256 of all three scripts still matched `scripts.sha256`, recorded
before the first run. The scripts matched again after the run.

What the full run shows (chapter 13 says it): takeover data alone give
`mu = 2.0·10⁻⁵`, `sm = 1.58` (true 3·10⁻⁵ and 1.5); the two-stage fit on
heteroplasmy distributions recovers `mu = 2.95·10⁻⁵`, `sm = 1.485`, with a
held-out loss of 30.1 against 41.2 for the true parameters and 143.2 for
stage 1. The prose of the old tutorial was rewritten to these numbers.

## How the outputs were produced

- Wheel: the 0.7.0 wheel that D3 built from this checkout
  (`/home/claude/docs071/d3/dist/`), installed in the campaign venv; release
  CLI `target/release/crn-craft`. fcmaes 2.1.0, scipy 1.18.1 and loguru
  0.8.0 were installed into the venv for chapter 13.
- Session pages: the D3 driver with `--write --binary
  target/release/crn-craft --require-cargo`, under the CPU lock. No output
  was typed. Prose was edited after the write runs only outside the fenced
  blocks; the final verification below compares every block again.
- Chapter 13: each script ran in full, under the CPU lock, from the
  repository root (`python examples/mtdna/SCRIPT.py`), its stdout captured
  to a file, resumable per script; the `text` blocks are those files,
  inserted by a script, unedited. The scripts' SHA-256 were recorded before
  the run and checked after it.
- Chapter 15: R 4.3.3 is installed on this host, with Ubuntu's reticulate
  1.35.0. CRAN, Posit's package manager and r-universe are refused by the
  egress proxy (403), so a newer reticulate could not be installed.
  reticulate 1.35 failed on the campaign venv (NumPy 2.3, pandas 3.1) with
  "C stack usage is too close to the limit". In a separate venv with
  NumPy 1.26.4 and pandas 2.3.3 and the same 0.7.0 wheel, the example ran
  with `import("crn_craft")`, and its output is the page's output (the same
  numbers as the 0.6.0 page). The event-log lines and three pitfall
  messages were also run. The page keeps `<!-- check: none -->`.

## Driver changes

`/home/claude/docs071/run_chapter.py` (outside the repository; the copy
before D4 is `/home/claude/docs071/d4/run_chapter.py.before-D4`):

1. **Smoke mode for long pages.** Without `--long`, a block of a
   `<!-- check: long -->` page runs with `CRN_CRAFT_SMOKE=1` and its output
   is not compared ("not compared: smoke"); with `--long` it runs in full
   and its output is compared or written. Before, such blocks did not run
   at all without `--long`.
2. **`<!-- file: PATH -->`** before a block binds it to a repository file:
   the block must equal the file byte for byte, and a standalone page runs
   the file itself (`python PATH` from the repository root), so that
   `noise_and_parallel.py` can import `fit_heteroplasmy`.

The fitting scripts read `CRN_CRAFT_SMOKE`: with `1`, every ensemble has
1% of its cells (at least 16) and every optimizer three evaluations. A
smoke run of the page takes about 10 seconds.

## Block ledger

Every block that the ledger (section 2.4) assigns to these pages. "ran: yes"
means the block ran in the final verification and its shown output matched.

| Ledger id | File | Section | Model | Ran |
| --- | --- | --- | --- | --- |
| M22 | `12-mtdna.md` | The model (setup block `three_class`, and a run of 100 cells) | builder | yes |
| M1 | `12a-mtdna-takeover.md` | A first ensemble | builder | yes |
| M2 | `12a-mtdna-takeover.md` | The mean hides the typical cell | builder | yes |
| M5 | `12a-mtdna-takeover.md` | Takeover fraction and time to threshold | builder | yes |
| M24, parts 1 and 3 | `12a-mtdna-takeover.md` | Variants | builder | yes |
| M6 | `12b-mtdna-events.md` | Class appearances and extinctions | builder | yes |
| M7 | `12b-mtdna-events.md` | Mutation events in a log | builder | yes |
| M8 | `12b-mtdna-events.md` | Stop at the first advantaged mutant | builder | yes |
| M15 | `12c-mtdna-atp.md` | ATP as a continuous variable | builder (`add_ode`) | yes |
| M9, part 1 | `12d-mtdna-methods.md` | The same cells under every method (seven methods) | builder | yes |
| M12 | `12d-mtdna-methods.md` | A compiled model (last block) | builder | yes, with cargo |
| M19 | `12e-mtdna-large-studies.md` | Parameter sweeps | builder | yes |
| M26 | `12e-mtdna-large-studies.md` | From the command line (the Python block) | builder, `to_json()` | yes |
| M26, CLI | `12e-mtdna-large-studies.md` | From the command line (`bash`: check, run, verify; not run before) | file from the builder | yes |
| — | `12e-mtdna-large-studies.md` | From the command line (reads the three outputs) | | yes |
| O1 | `13-fitting.md` | A first attempt with takeover fractions (`fit_takeover.py`) | builder | yes: full run once; smoke in the checker |
| O2 | `13-fitting.md` | A better objective with heteroplasmy distributions (`fit_heteroplasmy.py`) | builder | yes: full run once; smoke in the checker |
| O3 | `13-fitting.md` | Noise and the shape of the objective (`noise_and_parallel.py`) | builder | yes: full run once; smoke in the checker |
| M9, part 2 | `14a-mtdna-model-files.md` | The fast path (file under `mtdna_direct` and `static_direct`, builder under `static_direct`) | document | yes |
| M24, part 2 | `14a-mtdna-model-files.md` | The fast path (two fallbacks: an unguarded stop, a builder model) | document and builder | yes |
| M14 | `14a-mtdna-model-files.md` | Every lineage with the clone model | document | yes |
| M18 | `14a-mtdna-model-files.md` | The stored ATP block and thinning | document | yes |
| N11 | `migrating.md` | From the 0.6.0 module to the top level | builder (`crn.Model`, `crn.rules.Model`) | yes |
| R example | `15-from-r.md` | Example: SIR model | builder | not run by the checker; run in R by hand |

With D3's ledger, every manual block M1 to M26 now has a destination: M3,
M4, M10, M11, M13, M16, M17, M20, M21, M23, M25 in D3's pages, the rest
here. Setup blocks of 12a to 12e, 14a and Migrating are not in the ledger;
each is the model of the index page (or the strain model of chapter 7)
and runs before every section.

## Blocks that changed beyond the specification

1. **`three_class`** has `delta` as a parameter and degradation written as
   `delta * W` (the manual's M22 used the number), so that
   `bind_parameters` can change the half-life. It takes `m0`, the initial
   neutral mutants, for the variants.
2. **M1** prints the first five takeover days rounded to 0.1.
3. **M2** prints the median of `A` from the statistics and the median and
   97.5% quantile of `M + A` from the kept paths; the check against NumPy
   is on `A`.
4. **M6**: the first-passage tracker follows `W` with the predicate over
   the three species; the time-averaged count is rounded.
5. **M24, part 3**: the error is a misspelled species in an added reaction
   (`delta * WW`), caught when the model is used.
6. **M15**: the ATP model is a function `atp_cell` in the section block,
   with `add_ode` and the default `on_bound="stop"`; the lines that printed
   the file's fields are gone, as specified.
7. **M9, part 1** prints "mutant fraction" for the final mean of `M + A`
   over all molecules; 14a's part 2 runs three arms (the file under both
   methods and the builder model).
8. **M19** runs 200 cells per candidate instead of 400 (the specification
   allowed fewer cells); the block takes about a minute.
9. **M26**: 200 cells instead of 1,000, `static_direct`, a tracker on `A`;
   the archive is `takeover-archive` in the working folder (the parent of
   an archive must exist).
10. **M14** prints the median lifetime rounded and no longer prints the
    final instances of cell 0.
11. **N11** compares the two modules on twenty strains (same sampler, the
    top level adds `median` and `quantiles`) and shows the refusal of a
    hundred strains by `crn.rules.Model`, and the top level running them.
12. **O1 to O3**: the model is the builder model, the columns are `W`, `M`,
    `A`; the mean mutant fraction is `(M + A)` over the sum of the columns;
    the heteroplasmy summary reads `M + A` and the sum over `variable` from
    the kept paths; the smoke switch is new. The optimizer settings, seeds,
    cells and evaluations are those of the old scripts.

## Checks run

- Final run of the ten pages, under the CPU lock: with the release CLI and
  `--require-cargo`; the same under `python -O`; and with the `crn-craft`
  command of the wheel. All passed, and the three reports are
  byte-identical: 31 runnable blocks ran (28 in session pages, the three
  fitting scripts in smoke mode in 20 s), 19 outputs compared, four `r`
  blocks listed as not run.
- Negative test: one changed character in the `fit_takeover.py` block of
  chapter 13 fails the file binding.
- The fitting scripts in smoke mode by hand, and in full once (above).
- `validation/check_links.py`: no broken link or anchor from these pages
  except the forward links below; 739 links and 33 stubs checked, 5 broken,
  all to pages of other tasks (`comparison.md`, `howto.md`,
  `reference/python-api.md`, `developer/rust-api.md` twice).
- `ruff check` on `examples/mtdna`, the timing script and the driver;
  `ruff format --check validation`.
- The R example and the event-log lines in R 4.3.3 (above).

## Forward links

Links in these pages whose targets other tasks write:
`../developer/rust-api.md` (Migrating, chapter 13; D6), and the anchors of
D3's pages, which resolve. `python/README.md` still links to the old R page
by an absolute URL (D5's file).

## Open questions

1. **The builder form costs 9.9 times the fast path, not 2.8.** The pages
   say so with the measured ratios. For sweeps and fits, where candidates
   run interpreted, it is the full factor: chapter 13's full run took
   2 h 12 min. Keep the builder form for chapter 13, as decided, or move the
   fitting scripts to the model file and `mtdna_direct` (the alternative of
   question 1)? A faster evaluation of whole-rule propensities, or a
   certificate of the fast path for builder models, would remove the
   question; both are engine work for a later release.
2. The task prompt names `14-fitting.md` and `15-from-r.md`; the table of
   contents, `SUMMARY.md` and the coordinator's note name `13-fitting.md`.
   I used `13-fitting.md`.
3. Chapter 15 was run with reticulate 1.35 and NumPy 1.26, because the
   proxy refuses CRAN. A run with a current reticulate and NumPy 2 would
   confirm the page's sentence on version 1.47.
