# D10 log: the documentation after E1

Date: 13 November 2026. Branch `docs-0.7.1`, main checkout, on top of the
merge of E1 (`8d4b7cb`). Nothing pushed. Task: bring the documentation in
line with E1 (faster builder models; `mtdna_direct` certifies the builder
form under `run/1.0`) and close the open items of D6 and D8.

## Setup

- Wheel of this checkout built with maturin under the CPU lock and
  installed into `/home/claude/m07/venv-m10a` (`uv pip install --reinstall
  --no-deps`); rebuilt twice after docstring edits. Release CLI
  `target/release/crn-craft` (`-p crn-craft-cli`).
- A first check on the new wheel: the builder model of chapter 12 under
  `crn.ensemble(..., method="mtdna_direct")` reports
  `chacha8-u53-ziggurat/rules-mtdna-direct/1` and no diagnostics, also with
  `keep_paths="all"` and for candidates of `bind_parameters`; so does
  `crn.rules.ensemble` with a `crn.Model`, since a `model/1.0` document
  always selects `run/1.0`.

## Commits

| Commit | Content |
| --- | --- |
| `a05345b` | fitting scripts on `mtdna_direct` (builder model kept); `ruff format` of `fit_takeover.py` and `fit_heteroplasmy.py`; chapter 13's bound blocks follow |
| `96baba0` | Rust doc comments link the guide and the new reference paths |
| `36619f6` | CHANGELOG |
| `37d146e` | chapter 13: outputs and prose of the full fitting run |
| `e18c470` | docstrings; `python-api.md` regenerated |
| `ce37d6c`, `70e8843` | guide pages and the methods reference; outputs of 12d and 14a regenerated; `guide-blocks.json` (12d: 4 blocks) |
| this commit | this log |

## 1. Numbers

New numbers, from `validation/merge/RESULTS.md`, "E1 ... results", 13
November 2026: E/A 3.6, E/B 1.40, E/F 2.6, F/A 1.4, G/A 0.97.

| Page | Before | Now |
| --- | --- | --- |
| `12-mtdna.md`, Conventions | fast path certifies only the file; 9.9 and 3.6 | the fast path certifies the builder model under `run/1.0` (`crn.ensemble`); builder on the fast path 0.97 of the file; under Static Direct 3.6 times the fast path and 1.4 times the file under Static Direct; "use the fast path for long runs". The case study keeps `static_direct` as its method (it runs every variant, the ATP model included); no output changed |
| `12-mtdna.md`, The model as a file | files add "the fast path" | removed from that list |
| `12d-mtdna-methods.md` | seven methods; "the fast path needs the model file"; compiled gain 7.1; "evaluates every rate at every event" | eight methods (`mtdna_direct` added to the first block); new section **The fast path** with a checked block (`crn.ensemble`, 200 cells, `static_direct` and `mtdna_direct`, prints `crn-craft.result/1.0` and no warnings); compiled gain 2.6, 1.4 of the fast path; "runs as a small program after every event that changes its species" |
| `13-fitting.md` | Static Direct, 21 s per evaluation, "load the file for the fast path" | scripts pass `method="mtdna_direct"` with the builder model; about 2 s per evaluation; Static Direct about 3.6 times as long; the last bullet of "Making each evaluation cheap" rewritten; every output and the prose that reads them from the new full run |
| `14a-mtdna-model-files.md` | "three things need the files"; certificate does not read the builder; 9.9; fallback block expects whole-rule fallback | two things need the files; the certificate's builder form described with a link to the methods reference; the 64-cell block gains the arm "builder, mtdna_direct" (same events and takeovers as the file); 2.8 (file, 12 November) kept, 0.97, 3.6 and 1.4 with 13 November; the fallback block rerun: the builder model now runs the fast path, prose says so |
| `12e-mtdna-large-studies.md` 78 | "this is where the builder form costs most" | sweeps of the plain model can pass `mtdna_direct`, about 3.6 times as fast |
| `07-methods.md` table and Other methods | "the three-class mtDNA model as a file" | built with `crn.Model` or a model file; links 12d and 14a |
| `14-model-files.md` 209–217 | fallback only | adds that the certificate reads whole-rule propensities under `run/1.0`, so the builder form runs the fast path |
| `reference/methods.md` (whole-rule propensities) | "`mtdna_direct` does not certify them" | certifies them only in the builder form of the three-class model under `run/1.0` |

Searched also: `README.md`, `python/README.md`, crate READMEs, `index.md`,
`howto.md`, `comparison.md`, `troubleshooting.md`, `05-speed.md`,
`migrating.md`, `reference/*`, `developer/*`. None quotes the builder
ratios or the builder fallback. Unchanged as specified: README 2.1×,
`05-speed.md` 2.1×, `comparison.md` 2.14 and 2.84. `howto.md` and
`troubleshooting.md` describe `CRN-R051` in general, which stays correct.
`evidence/` pages are dated measurements and were left as they are.

Byte identity seen in the rerun outputs: every unchanged arm of 12d and
14a printed the event counts and takeovers of the 0.7.0 outputs; only the
timings changed (Static Direct on the builder model 1.69 to 0.55 s in
12d's first block, 5.46 to 1.71 s interpreted in the compiled block).

## 2. Fitting chapter

- `examples/mtdna/fit_takeover.py` and `fit_heteroplasmy.py`:
  `method="mtdna_direct"` with a comment, the builder model kept;
  formatted with `ruff format` (only line wrapping changed otherwise).
  `noise_and_parallel.py` imports `fit_heteroplasmy` and needed no change.
  The three bound blocks of chapter 13 were replaced by the files' text.
- Smoke runs of the three scripts passed before the full run.
- Full run with D4's runner `/home/claude/docs071/d4/fit/run_full.sh`,
  unchanged, under the CPU lock, one lock hold per script. D4's outputs were
  moved to `/home/claude/docs071/d4/fit/before-E1/` first so that the
  runner's resume markers did not skip the scripts. SHA-256 of the scripts
  recorded before (`d4/fit/scripts.sha256`) and checked after the run: OK.
  Times: 233 s, 441 s and 105 s (13 minutes; was 2 h 12 min).
- The outputs were inserted unedited. What changed in the story: the
  takeover-only fit is now worse (`mu = 1.87·10⁻⁵`, `sm = 2.27`; held out
  13.3% taken over against 10.8%), which shows the confounding more
  clearly; the two-stage fit gives `mu = 3.06·10⁻⁵`, `sm = 1.475` with
  held-out losses 68.6 (stage 1), 9.8 (stage 2) and 6.8 (true), so the
  paragraph that explained stage 2 beating the truth was replaced; the
  noise scan's minimum is at `sm = 1.60` (was 1.40), the same point where
  stage 1 landed; both ways of using two cores give about 2 s per
  evaluation.

## 3. Rust doc comments

Every `//!` link to a moved page now names its new path: `lib.rs`
(user guide, architecture, model language, contracts, archives, technical
report, diagnostics), `rules/mod.rs` (model reference, methods and the
mtDNA case study instead of `mtdna-rules.md` and the manual, kernel
checklist, architecture, diagnostics), `rules/model.rs`, `rules/atp.rs`
(method reference, ATP section), `generated.rs`, `generated_loop.rs`,
`archive.rs`, `diagnostic.rs`, `export.rs`, and `crn-craft-rules/src/lib.rs`.
Not changed: three runtime strings name `docs/rules-language.md`
(`rules/model.rs` 581, `rules/mod.rs` 2871, `unified.rs` 138). They are
written into diagnostics and `describe` output, so changing them changes
outputs; the stub at that path keeps them working. Candidate for 0.7.2.

## 4. Docstrings (D6's list)

| D6 item | Fix |
| --- | --- |
| `PreparedEnsemble` has no docstring | docstring added (the page's former own sentence) |
| `importers.ImportFailure` has no docstring | docstring added |
| `Model.add_reaction` refers to "the module documentation" | the counting rule is stated in place |
| `rules.ensemble` says 3.3× and 3.9× | 2.8 and 3.1 times less time per event (12 November, file form), and the builder form under `run/1.0` |
| `import_sbml`, `import_bng_net` do not say what they return | the `crn-craft.import/1.0` report and `ImportFailure` |
| `Result` names `TERMINATIONS` | `crn_craft.rules.TERMINATIONS` and its values; also `FspResult` named `fsp` (lint: "reference not on this page"), now `crn_craft.rules.fsp` |
| `ensemble`, `prepare`, `loop_source` refer to `crn_craft.rules.*` | not changed: `crn_craft.rules` is a documented, importable module, and the page states that `crn_craft.rules.ensemble` takes the same arguments |

`MISSING_DOC` in `make_python_reference.py` is now empty (the mechanism
stays). `--lint` reports no problem. `python-api.md` regenerated with
`--write`.

## 5. CHANGELOG

The E1 entries match the results: 186.7 to 68.4 ns per event is 2.7 times;
E/B 1.40; G/A 0.97; law check 120 comparisons, 0 failures. The first
sentence names both things. Added: the ns figures, "under Static Direct",
the 0.97 and the law check; one Documentation bullet for the fitting
scripts, docstrings and Rust doc links.

## 6. Checks

See the report below.

## Report

Checks on the final tree (`70e8843` plus this log), wheel of this checkout,
release CLI, each heavy step under the CPU lock (runner
`/home/claude/docs071/d10/checks/run.sh`, logs beside it):

| Check | Result |
| --- | --- |
| `check_guide.py --binary target/release/crn-craft --require-cargo` | passed: 26 files, 109 runnable blocks ran (was 108; 12d has 4), 73 outputs compared, 6 `sh` syntax-checked, 3 fitting blocks in smoke mode, 4 `r` not run |
| the same under `python -O` | passed; report byte-identical |
| frozen manual, `check_manual.py --manual validation/compat/manual-0.6.0.md --require-cargo`, normal and `-O` | 22 of 22; reports identical; the file is unchanged |
| `check_links.py` | 948 links, 36 stubs, 53 SUMMARY entries, 0 broken |
| `make_python_reference.py --check` | up to date; `--lint` empty |
| `check_rust_docs.py` | 2 rust programs compiled and ran, 0 failed |
| `mdbook build` (mdBook 0.5.0) and `check_links.py --site site` | no warnings; 3,554 references, 0 broken |
| `ruff check` and `ruff format --check` on `validation python examples` | clean (124 files) |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | clean |
| `cargo test --release --workspace --locked` | 406 passed, 0 failed |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | clean |

Deviations:

1. The case-study pages 12a to 12e keep `method="static_direct"`; only 12d
   (new section and the methods table), 14a and chapter 13 run
   `mtdna_direct` on the builder model. Switching every page would change
   every output of the case study and drop the general method from pages
   whose variants (ATP, other stops) need it.
2. D4's runner was used unchanged; its earlier outputs were moved to
   `d4/fit/before-E1/` so that its resume markers did not skip the run.
3. Three runtime strings that name `docs/rules-language.md` stay (they are
   program output); the stub resolves them.
4. `evidence/` pages keep their dated numbers.

Open questions:

1. Should the case study switch to `mtdna_direct` as its default method in
   a later release (faster pages, all outputs regenerated)?
2. The numbered scripts in `examples/mtdna/` are still not run by a checker
   (D8's note).
