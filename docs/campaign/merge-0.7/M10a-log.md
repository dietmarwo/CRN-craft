# M10a: tutorial

Branch `merge/M10a`, from `merge-0.7` at 7da5060 (M8b committed) merged with
`merge/M10c` (M10b manual and references, M10c architecture documents). The old
engine is deleted in this tree. Only documentation, the tutorial checker and
the tutorial's asset scripts changed; no Rust code.

## What changed

- `docs/tutorial.md` is rewritten on the unified API and command line
  (`import crn_craft as crn`, `crn-craft run|check|convert|describe|analyze|
  export|verify`, `crn-craft.model/1.0`, `crn-craft.run/1.0`). Sections: install,
  first model with the builder, checking and diagnostics, exact simulation,
  results as data frames, quantiles and the median, `Result.save`, oscillations,
  FSP, tau, budgets and cancellation, reproducibility, analysis (with sampled
  first passage), import/convert/export, choosing a method (selectors and
  `auto`), where to go next. The network and rules engines are mentioned once,
  in a short history note in section 1.
- Two sections are new: quantiles and the median (with a new figure
  `quantiles.png`) and saving a result.
- `validation/check_tutorial.py` now checks every fenced block, not only
  `bash` and `python`:
  - `text` blocks must equal the output of the runnable block before them
    (before, `text`, `json` and `sh` blocks were skipped);
  - `json` blocks must be strict JSON (duplicate keys and non-finite
    constants rejected), and a `model/` or `run/` document must pass
    `crn-craft check` (a request against the model it names, so a stale hash
    fails);
  - `sh` blocks (install commands) must pass `bash -n`; they are not run;
  - any other fence language is an error;
  - fences may have four or more backticks, so a block can show markdown that
    has fences (the saved report);
  - the minimum count of runnable blocks (15) is kept, and now counts only
    `bash` and `python` blocks.
  Nothing it checked before is dropped. The report now has one row per
  `text`, `json` and `sh` block too, and is identical under `-O`.
- `examples/tutorial/workdir.py` adds `dimer.xml` (SBML) and
  `immigration-death.network.json` (a copy of the 0.6.0 network model from
  `models/legacy/`), which the import and convert sections read.
- `examples/tutorial/make_assets.py`: new `quantile_band` figure and a
  `12-convert` terminal session.
- `examples/tutorial/render_terminals.py`: new option `--browser PATH` (and
  environment variable `CRN_TUTORIAL_BROWSER`), passed to Playwright as
  `executable_path`. The default behavior is unchanged.
- Images in `docs/images/tutorial/`: all figures and all twelve terminal
  screenshots regenerated (`06-cached.png` removed, `06-replay.png` and
  `12-convert.png`, `quantiles.png` added). Every image was viewed after
  rendering. Every image the tutorial references exists, and no image is
  unreferenced. `numbers.json` regenerated.

## Tests first

The new checker rules were written before the text: a `text` block with a
wrong output, a `json` block with a duplicate key, a model with an unknown
name and a request with a stale hash were each fed to the checker functions
and rejected (run by hand, not kept as a test file). The tutorial's first
draft had placeholders in every `text` block and failed until each was filled
from the real output. A request without `keep_instances` failed the checker,
and two statements of the draft were wrong (a Python result's metadata is not
equal to the command line's result, because the request hash and the
`replicates` records differ; `quantiles=()` and the file bytes fix the
comparison) and were corrected from the run.

## Checks

| Check | Result |
| --- | --- |
| `validation/check_tutorial.py` normally and under `-O` | pass, 59 blocks, reports identical |
| `validation/check_installed.py` from `/tmp` with `python -I`, normally and `-O` | pass, 18 checks, reports identical |
| `validation/check_links.py` (docs, READMEs, tutorials) | 438 links, 0 broken |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| `examples/tutorial/make_assets.py`, `render_terminals.py --browser /usr/bin/google-chrome` | pass; 12 PNGs rendered |

No Rust or other Python file changed, so the cargo set and I1 to I3 were not
rerun. The wheel and binary were built from this worktree (`maturin build
--release --locked` and `cargo build --release --workspace --locked`, target
directory `/home/claude/m07/wt-m10a/target`, removed at the end).

## Numbers that changed (new stream)

| Quantity | Old tutorial | New |
| --- | --- | --- |
| immigration-death events, mean at t=30 (400 paths, seed 42) | 87,508, 19.7 | 88,037, 20.09 |
| tau variance ratio at epsilon 0.03 / 0.1 | 1.18 / 1.68 | 1.09 / 1.83 |
| sampled first passage, grid 0.05: detected of 500 | 460 | 464 (upper bound 6.075) |

The tau text changed with the numbers. With 1,000 paths the variance ratio
has a standard error of about 6%, and the ratios at epsilon 0.003, 0.01 and 0.03
(1.11, 1.01, 1.09) lie within two standard errors of 1. The text says so, and
only the 83% at epsilon 0.1 is called an inflation. The old text claimed 18% at
0.03.

## Notes and limitations

- The `--version` screenshot (`01-help.png`) shows version 0.6.0, because the
  workspace version is bumped in M11. M11 must render `01-help` again (run
  `make_assets.py` and `render_terminals.py`).
- Printed outputs in `text` blocks include pandas and NumPy formatting, so
  they depend on the pinned package versions of CI.
- The tutorial has no `json` block of a model. The checker supports them and
  one `json` block (the FSP request) is checked.
- The installation `sh` blocks cannot be run by a checker; they are only
  syntax-checked.
- `crn-craft` in the checker's scratch directory is the native binary. The
  Python console script of the same name has the same commands.
- The terminal transcripts (`*.txt`, `*.terminal.html`) written by
  `make_assets.py` are not committed.

## Migration notes

None for users. Anyone running `check_tutorial.py` must have `pandas` and
`xarray` installed (as before, through the wheel).
