# D3 log: the merged user guide, chapters 1 to 11 and 14, and Troubleshooting

Date: 13 November 2026. Branch `docs-0.7.1`, on top of the D2 commits
(`5bf52dc`). Nothing pushed. Specification: the table of contents of task D1,
second version (`docs-toc.md` of the campaign, sections 1.7, 2.1, 2.2, 2.4
and 9).

## What was written

Thirteen pages of `docs/guide/`: chapters 1 to 11, chapter 14 and
`troubleshooting.md`. Every model is built with the builder, except where the
table of contents puts a model document (chapter 14) or a file of another
format (chapter 11). Every page uses `import crn_craft as crn`.

Chapter 14 is not named in the scope sentence of the D3 task ("the basics,
the topic chapters, troubleshooting"), and the D4 task excludes it (D4 takes
the page that follows it). It would have had no writer, so D3 wrote it.

## Commits

| Commit | Content |
| --- | --- |
| `876f1bf` | chapter 1, Install |
| `82199f6` | chapter 2, A first model |
| `a98611e` | chapter 3, Running ensembles |
| `cde7eb4` | chapter 4, Results |
| `942996b` | chapter 5, Speed: threads, parameter sweeps and compiled models |
| `5d88650` | chapter 6, Stop conditions, trackers and event logs |
| `99c9bf5` | chapter 7, Choosing a method |
| `0a90b4c` | chapter 8, Reproducing runs and published results |
| `32d28b2` | chapter 9, Continuous variables |
| `528b270` | chapter 10, The command line, request files and archives; `make_assets.py`, `workdir.py`, `examples/guide/dimer.xml`, six screenshots, `numbers.json` |
| `fac5e79` | chapter 11, Analysis, import and export; two screenshots |
| `a531ace` | chapter 14, Model files and species with attributes |
| `97e6d34` | Troubleshooting |
| this commit | two spellings, this log |

## How the outputs were produced

- Wheel: built from this checkout with maturin 1.9.6 under the CPU lock
  (`crn_craft-0.7.0-cp310-abi3-manylinux_2_34_x86_64.whl`) and installed into
  the campaign's virtual environment. Release CLI: `target/release/crn-craft`
  of the same checkout (cargo reported it up to date).
- Driver: `run_chapter.py` of the campaign folder, outside the repository. It
  implements section 1.7 for one file per process: a fresh scratch directory,
  the fixtures comment, the fence languages, exact comparison with timings
  normalized, `standalone`, `none` and `long` pages, the cargo rule and a
  count file. With `--write` it replaces every `text` block by the real
  output of the block before it. No output in these pages was typed or
  edited by hand.
- Section-local state, as the driver reads point 3 of section 1.7: at every
  level-2 heading the Python namespace is replaced by a fresh one in which
  the setup block (the first `python` block of the file) has run again. The
  scratch directory is kept for the whole file. Chapters 10, 11 and 14 pass
  files, not names, from one section to a later one, and say so.
- Final run of all thirteen pages, under the CPU lock: with the release CLI
  and `--require-cargo`; the same under `python -O`; and with the `crn-craft`
  command of the wheel instead of the release CLI. All three passed, and the
  three reports are byte-identical: 73 runnable blocks ran, 49 outputs
  compared, 2 `sh` blocks syntax-checked.
- Negative tests of the driver: a wrong output, a name used from another
  section, an unknown fence language, a `text` block without a runnable
  block, a `json` block with a duplicate key and a count below the expected
  number all fail.
- `validation/check_tutorial.py` still passes on `docs/tutorial.md` (59
  blocks) after the change to `workdir.py`.
- `validation/check_links.py`: no broken link between the pages of this
  task, anchors included. 22 links are broken until their targets exist; see
  "Forward links" below.

## Lines per page against the budget

| Page | Lines | Budget | Difference | Runnable blocks | Count in the table of contents |
| --- | ---: | ---: | ---: | ---: | ---: |
| `01-install.md` | 97 | 85 | +14% | 1 | 1 |
| `02-first-model.md` | 192 | 175 | +10% | 5 | 5 |
| `03-running.md` | 208 | 200 | +4% | 6 | 6 |
| `04-results.md` | 226 | 195 | +16% | 7 | 7 |
| `05-speed.md` | 194 | 200 | −3% | 4 | 4 |
| `06-stops-trackers-logs.md` | 219 | 205 | +7% | 5 | 5 |
| `07-methods.md` | 291 | 250 | +16% | 6 | 6 |
| `08-reproducibility.md` | 181 | 200 | −10% | 4 | 3 |
| `09-continuous-variables.md` | 174 | 160 | +9% | 3 | 2 |
| `10-command-line.md` | 293 | 250 | +17% | 15 | 14 |
| `11-analysis-import-export.md` | 229 | 200 | +15% | 10 | 9 |
| `14-model-files.md` | 289 | 250 | +16% | 7 | 7 |
| `troubleshooting.md` | 69 | 60 | +15% | 0 | 0 |
| Total | 2,662 | 2,430 | +10% | 73 | 69 |

Four pages are one or two points over the 15 percent margin: chapters 4, 7,
10 and 14. In all four the fenced blocks alone exceed the block budget of
the table of contents (chapter 7 has the 34-line strain block of M10 with
its five-line output; chapter 10 has fifteen blocks, nine screenshots and
two tables; chapter 14 has the two model documents; chapter 4 has the new
export-and-plot block). Prose was cut in each of them; no checked block was
dropped.

The four extra runnable blocks: one in chapter 8 (it prints the fields that
identify a run), the setup block of chapter 9 (the table of contents did not
count one), the block that writes `long.run.json` in chapter 10 (in the
section that uses it), and a `bash` block for `crn-craft export` in chapter
11, which no checker ran before.

## Block ledger

Every block of section 2.4 that belongs to these pages. "ran: yes" means the
block ran in the final run and, where an output is shown, the output
matched.

| Ledger id | File | Section | Ran |
| --- | --- | --- | --- |
| S1 | `01-install.md` | Install the Python package | `sh`, syntax-checked |
| S2 | `01-install.md` | Build from source | `sh`, syntax-checked |
| N1 | `01-install.md` | Check the installation | yes |
| T1 | `02-first-model.md` | Species, parameters and reactions (setup block) | yes |
| N13 | `02-first-model.md` | A first run | yes |
| T4 | `02-first-model.md` | Checking a model | yes |
| T6 | `02-first-model.md` | Checking a model | yes |
| M20 | `02-first-model.md` | A second model with two reactants | yes |
| T7 | `03-running.md` | Many trajectories in one call | yes |
| M4 | `03-running.md` | Seeds and threads | yes |
| T15 | `03-running.md` | Budgets | yes |
| M25 | `03-running.md` | Budgets (the wall-time limit) | yes |
| T16 | `03-running.md` | Stopping a run | yes |
| T9 | `04-results.md` | What a result holds | yes |
| T10 | `04-results.md` | Quantiles and the median | yes |
| N2 | `04-results.md` | Single trajectories | yes |
| N3 | `04-results.md` | Export and plot | yes |
| T11 with M3 | `04-results.md` | Save a run in one file | yes |
| B5 | `04-results.md` | Save a run in one file (a `python` block) | yes |
| N14 | `05-speed.md` | Threads | yes |
| N8 | `05-speed.md` | Parameter sweeps | yes |
| M13 | `05-speed.md` | Compiled models (last block of the page) | yes, with cargo |
| M21 | `06-stops-trackers-logs.md` | Stop conditions | yes |
| N4 | `06-stops-trackers-logs.md` | Trackers | yes |
| N5 | `06-stops-trackers-logs.md` | Event logs | yes |
| N6 | `06-stops-trackers-logs.md` | Callbacks | yes |
| T22 | `07-methods.md` | Static Direct and its selectors | yes |
| M10 | `07-methods.md` | Static Direct and its selectors | yes |
| M23 | `07-methods.md` | Larger models | yes |
| T13 | `07-methods.md` | Exact distributions with finite state projection | yes |
| T14 | `07-methods.md` | Approximate leaping | yes |
| T17 | `08-reproducibility.md` | What makes a run reproducible | yes |
| N10 | `08-reproducibility.md` | Reproducing a published simulation | yes |
| M16 | `09-continuous-variables.md` | setup block (the model) and A variable between events (the run) | yes |
| M17 | `09-continuous-variables.md` | Independent variables | yes |
| T5 | `10-command-line.md` | A request file | yes |
| B1 | `10-command-line.md` | Check a model and a request | yes |
| B2 | `10-command-line.md` | Check a model and a request | yes |
| B3 | `10-command-line.md` | Run, archive and verify | yes |
| N9 | `10-command-line.md` | Run, archive and verify | yes |
| B9 | `10-command-line.md` | Replay | yes |
| B4 | `10-command-line.md` | Python and the command line give the same result | yes |
| T8 | `10-command-line.md` | Python and the command line give the same result | yes |
| J1 | `10-command-line.md` | Requests for fsp and tau (printed output, compared; the file is run by B6) | yes |
| B6 | `10-command-line.md` | Requests for fsp and tau | yes |
| B7 | `10-command-line.md` | Requests for fsp and tau | yes |
| B8 | `10-command-line.md` | Stop a run with a cancel file | yes |
| T18 | `11-analysis-import-export.md` | Conservation laws and structure | yes |
| T12 | `11-analysis-import-export.md` | Oscillations | yes |
| T19 | `11-analysis-import-export.md` | First passage from saved paths | yes |
| B10 | `11-analysis-import-export.md` | Convert a network model of release 0.6.0 | yes |
| T20 | `11-analysis-import-export.md` | Convert a network model of release 0.6.0 | yes |
| B11 | `11-analysis-import-export.md` | Import SBML and BioNetGen | yes |
| T21 | `11-analysis-import-export.md` | Import SBML and BioNetGen | yes |
| T2 | `14-model-files.md` | What the builder writes | yes |
| N15 | `14-model-files.md` | Save and load a model file | yes |
| T3 | `14-model-files.md` | Packaged examples | yes |
| M11, part 1 | `14-model-files.md` | A model document for three strains | yes |
| M11, part 2 | `14-model-files.md` | Run it | yes |
| N7 | `14-model-files.md` | New identities and random attributes | yes |

All 33 runnable tutorial blocks (T1 to T22, B1 to B11) are in these pages.
Of the manual's blocks, M3, M4, M10, M11, M13, M16, M17, M20, M21, M23 and
M25 are here. The others (M1, M2, M5 to M9, M12, M14, M15, M18, M19, M22,
M24, M26), N11, the fitting scripts and N12 belong to D4 and D5.

Outputs that the tutorial already showed came out unchanged on the builder
models: 88,037 events and the mean 20.09 (T7, T3), the quantiles 20, 12 and
30 (T10), the projection (T13), the leaping counts (T14), the analysis of
the enzyme mechanism (T18), `oscillatory 40.1` (T12), the first-passage
bounds (T19) and the imports (T21).

## Blocks that changed beyond the specification

1. **N1** prints the version and the schema of `crn.describe()`. Its output
   shows `0.7.0`. When the version is raised, chapter 1 must be rewritten
   with the driver.
2. **T6** prints one more line, the hint of the diagnostic.
3. **M20** prints "events per replicate" where the manual printed "events
   per cell".
4. **M4** keeps its five comparisons. Its last three lines, which printed
   the shape of the paths and one replicate's series, are not carried:
   chapter 4 shows path selections with N2.
5. **M25**, the wall-time limit: the run is so long that no replicate can
   complete (end 10⁷, limit 0.2 s), and the block prints the set of
   terminations. The manual's output listed `completed` beside
   `wall_time_limit`, which depends on the speed of the host.
6. **T9** prints one more line, the standard error, which the section on
   the standard error then reads. That section has no block of its own.
7. **B5** prints the first 23 lines of the report, to the end of the
   parameter values, not 30.
8. **M21** prints two more lines: the share of epidemics over by day 100 and
   three quantiles of the stop time, the two outputs the section names.
9. **M23** is under a level-3 heading, "Larger models", inside "Static Direct
   and its selectors". It reuses the function `strains` of M10, and names
   are local to a level-2 section. The anchor `#larger-models` is the same.
10. **M16**: the model is the setup block of chapter 9 and is named
    `resource`; the section's block runs it.
11. **M17** does not print the schema. Its refusal is shown on a second
    builder model, because dictionaries appear only in chapter 14.
12. **M11**: the document has the schema `crn-craft.model/1.0`. Part 1 writes
    it to `three-strains.model.json` and part 2 loads that file, because the
    two parts are in two level-2 sections. The rules are written more
    compactly; the document is the same.
13. **Chapter 10**: the setup block also makes the folder `runs`, because
    `--archive runs/first` needs its parent. The imports and `model_file` of
    T5 are in the setup block. J1 is printed one field per line. N9 has five
    lines and also shows that the archive written from Python holds the same
    result bytes as the one the command wrote.
14. **Chapter 11**: T18 writes `enzyme-mechanism.model.json`, which the new
    `bash` block of the export section reads. B10 and T20 read
    `models/legacy/immigration-death.model.json` and T21 reads
    `examples/guide/dimer.xml`, both under their repository paths, as the
    fixtures comment names them.
15. **N10** uses five report times, 10 to 50.

## Screenshots and figures

`examples/tutorial/workdir.py` has a new function `prepare_guide_workdir`. It
writes the files of chapters 10 and 11 with the chapters' builder lines and
`to_json()`. `make_assets.py` uses it and builds the models of the figures
with the same builder lines. `prepare_workdir`, the directory of the old
tutorial, is kept for `check_tutorial.py` until D8 retires that checker; it
now copies `examples/guide/dimer.xml`, the one source of the SBML text.

All assets were regenerated with the release CLI and rendered with
`render_terminals.py --browser /usr/bin/google-chrome`.

| Files | Result |
| --- | --- |
| `02-check`, `03-diagnostics`, `04-simulate`, `06-replay`, `08-tau`, `09-cancel` | regenerated and committed: the model hash is now `4dc25f88…`; the replay session runs the chapter's request |
| `11-export`, `12-convert` | regenerated and committed: the export reads the builder enzyme model; the convert reads `models/legacy/immigration-death.model.json` |
| `05-tamper`, `07-fsp` | regenerated, byte-identical to the committed files |
| `01-help`, `10-import` | byte-identical |
| the six figures | byte-identical: the builder models give the figures' numbers |
| `numbers.json` | only the three thread digests changed: they now come from the immigration–death model, as in chapter 8 |

Each of the eight committed screenshots was looked at before the commit.

## Chapter 8: what was verified

- The citation: Evans, Gillespie and Wilkinson, "The SBML discrete
  stochastic models test suite", Bioinformatics 24(2), 285–286 (2008),
  doi:10.1093/bioinformatics/btm566, confirmed with the Crossref record.
- The mean `100 exp(-0.01 t)` and the variance
  `2100 exp(-0.01 t) (1 - exp(-0.01 t))` of the birth–death process with
  rates 0.1 and 0.11 and 100 molecules: derived, and confirmed with
  `Model.fsp` on 601 states (the mean agrees to 5·10⁻⁹, the variance to
  3·10⁻⁷, with a mass defect of 9·10⁻¹¹). No data file of the suite is
  copied.
- Not verified at the source: the suite's own pages could not be opened from
  this host. The model's values (rates 0.1 and 0.11, 100 molecules) are
  taken from `validation/check_fsp.py`. The sentence "The test suite defines
  these two statistics" (Z and Y) is written from memory of the paper. The
  chapter gives no acceptance range in the suite's name.

## Forward links

These links are in the pages and resolve once the other tasks have written
their pages. The headings must keep the names of the table of contents.

| Target | Used by |
| --- | --- |
| `../comparison.md` | chapter 1 |
| `../howto.md` | chapter 1 |
| `../reference/python-api.md` | chapter 3 |
| `12-mtdna.md` | chapter 5 |
| `12a-mtdna-takeover.md#the-mean-hides-the-typical-cell` | chapter 4 |
| `12b-mtdna-events.md` | chapter 6 |
| `12c-mtdna-atp.md` | chapter 9 |
| `12d-mtdna-methods.md#matching-an-hourly-scheme` | chapters 7 and 8 |
| `13-fitting.md` | chapters 1 and 5 |
| `14a-mtdna-model-files.md`, with `#the-fast-path` and `#every-lineage-with-the-clone-model` | chapters 7 and 14, Troubleshooting |
| `15-from-r.md` | chapter 1 |
| `migrating.md`, with `#what-replays-and-what-does-not` | chapters 1, 8 and 11 |

## Open questions

1. Chapter 14 was written here although the task's scope sentence does not
   name it. If another task was meant to write it, this page is the draft.
2. The version in chapter 1 (item 1 of the changed blocks): keep the printed
   version and regenerate at each release, or print something that does not
   change?
3. Chapter 8: should the sentence on the two statistics of the test suite be
   confirmed against the paper, or reworded without naming the suite?
4. Section-local state: the driver resets Python names at every level-2
   heading and keeps files. Chapters 10, 11 and 14 rely on files that an
   earlier section wrote. If the final checker also resets the scratch
   directory per section, those three chapters need helper blocks in more
   sections.
5. Chapter 7 says plainly that on the hundred-strain model the selector
   matters little, because every event changes every rate. A model on which
   a selector pays off would be a better example for "Larger models"; none
   of the checked blocks has one.
6. Chapter 5 quotes one measured ratio for compiled models, 2.1 for the file
   form of the mtDNA model (12 November 2026). The block's own output on the
   shared host shows a smaller gain for SIR. The timings in the shown
   outputs are real but not controlled measurements.
