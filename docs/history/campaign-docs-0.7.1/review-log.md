# D9 review log: the documentation restructure and the engine change E1

Date: 13 November 2026. Branch `docs-0.7.1`. Reviewed state: `c436e1a`
(every documentation task, E1 and D10 merged), then the merge of the
case-study switch to the fast path (`9623068`, task FP), which arrived
during the review and was reviewed after it. Nothing pushed. The reviewer
did not write any of the reviewed work.

Method: step 1 as a first-time user in a clean virtual environment, before
reading the plan and the table of contents; step 2 as a reviewer of the
documentation and of the engine change E1. The walkthrough notes are
outside the repository (`/home/claude/docs071/D9-walkthrough.md`), as are
the scripts and reports named below (`/home/claude/docs071/d9/`).

## Result

The five requests of the researcher are met and the gates G1 to G7 hold
after the fixes below. The engine change is correct as far as it was
tested here: no near-miss of the builder certificate was accepted, and
every existing result is unchanged, with one exception that was not
documented and is now (finding 2). Fifteen findings were fixed in nine
commits; seven design matters are listed for the owner.

## Step 1: walkthrough

Wheel of this checkout (maturin, under the CPU lock) in a fresh venv
(uv, Python 3.12.3; numpy 2.5.4, pandas 3.1.1, xarray 2026.9.0); site built
with mdBook 0.5.0, no warnings.

- README, first screen: what the tool does and how to install it are clear
  within a minute, and the example ran as printed. The comparison is one
  click away, but on GitHub the link is on the second screen: the banner
  and the note on the related repositories come first, as the table of
  contents specifies.
- Chapters 1 to 4: every block ran and every output matched, except the
  plot block of chapter 4, which stopped without matplotlib. Four points
  where the reader had to guess or read another page (findings 5 and 6).
- Ten tasks from the how-to page: nine were answered one click away. The
  tenth, "run the same model from the command line and archive it", had
  two dead ends (finding 1).

## Findings

Severity: **high** means a reader is misled or a documented step fails;
**medium** means a gap that weakens a gate or a page; **low** is wording.

| # | Severity | Finding | Fix |
| --- | --- | --- | --- |
| 1 | high | Chapter 10, Replay: "run its `model.json` with its `request.json`" is refused with `CRN-R020` for every archive, also the one the chapter's own commands write. The archived request keeps its original model path. The Archives reference had it right. Nothing told the reader how to repeat from the shell a run that `crn.ensemble` prepared: the request it writes names the model `model.rules.json` and hashes other bytes than `model.to_json()` gives | `8e1738a`: the working steps as a checked block (`cmp` of the two results), and the route through `run.archive`. Checked by hand: an archive written from Python with `crn.ensemble(...).archive()` replays with the CLI to the same `result.json` bytes |
| 2 | high | One class of 0.7.0 results does not replay byte for byte, and no page said so. A `run/1.0` request that names `mtdna_direct` for a builder-built three-class model ran the Static Direct fallback in 0.7.0 (`chacha8-u53/rules-static-direct/1`, warning `CRN-R051`) and runs the fast path now. The request hash is the same (seen on both wheels: 651 against 675 events for one request) | `79343cc`: changelog, Migrating, chapter 8 and the methods reference say it, with the method that gives the 0.7.0 numbers again (`static_direct` with inverse waits and per-event trackers; equal records, paths, statistics and tracker columns on both wheels). Whether a patch release may change such a result is for the owner (D6) |
| 3 | medium | E1, tests of the certificate: nine of its checks had no test that fails without them. Found with 19 additional mutants: 10 caught, 9 survived. The checks themselves are in the code and refuse the models | `a321f71`: one near-miss per check in a new test; the nine mutants join `mutation_check.py` (67 mutants). All nine are caught now |
| 4 | medium | `docs/comparison.md`, Speed: the 34% outlier was attributed to GillesPy2 (it is a CRN-craft run; no GillesPy2 run is more than 7% from its median); "rebop repeats up to 1.6 times" (the report has 1.96 for one row); the 5.8 to 6.1 thread figures stood under "`mtdna_direct`, 2,000 cells" (they are for Static Direct, next reaction and the tree method) | `eab51bf`, with a correction note in `validation/merge/RESULTS.md` |
| 5 | medium | Chapter 4, Export and plot: the block needs matplotlib and said so only in chapter 1; on a fresh install pandas 3.1 prints a deprecation warning for `join(..., rsuffix=)` | `853efa3`: the sentence at the block; `add_suffix` |
| 6 | medium | `method="static_direct"` is in every call from the README on and was first explained in chapter 3. After the FP merge, page 12 used `mtdna_direct` two sections before explaining it | `853efa3`, `8fd6f6f` |
| 7 | medium | Page 12, Conventions (after FP): "a variant outside that structure runs Static Direct with `CRN-R051`" is wrong for the ATP model, which `mtdna_direct` refuses with `CRN-R020` | `8d56ddd` |
| 8 | medium | Page 12d: the 200-cell comparison prints 26 and 13 takeovers for the two methods and says only that they agree in distribution | `8d56ddd`: the expected count and its standard deviation. A separate run of 4,000 cells per method gave 10.75% and 10.80% |
| 9 | low | Three details of the 0.7.0 guides were on no page: `omitted_duration` and `omitted_peak` of a tracker, and the snapshot fields of an FSP result (`projected_first_moment`). Found by listing the 362 code identifiers of the old user documents; 350 are in the new pages, and the other nine are covered under another spelling or are placeholders | `1e92974` |
| 10 | low | Two-engine wording outside history and Migrating: "every simulation method of the one engine" (methods reference); "the network engine" in chapter 8 | `79343cc` |
| 11 | low | Chapter 7: the selector timings of the twenty-strain table differ by noise, and the text said "none is faster than the scan" under a row that is | `853efa3` |
| 12 | low | Chapter 2 gave no time unit for the SIR model; chapter 3 named "instances" without saying what they are | `8d56ddd` |
| 13 | low | Spelling: "towards", "cancelled" in prose | `1e92974` |
| 14 | low | The how-to entry on the fast path linked only page 14a | `1e92974` |
| 15 | low | The changelog did not say that the case study runs the fast path (after FP) | `d6b4bd3` |

Seen and left:

- `docs/mtdna-rules.md` (a stub) points to the methods reference only; the
  model description it once held is in chapter 12 and page 14a. The stub
  format allows one target.
- The E1 log and `RESULTS.md` give the stored-request tests as
  "`stored_requests.rs` 9, `static_direct_defaults.rs` 2"; it is 2 and 9.
- Older changelog entries and logs under `docs/campaign/` name a
  contributor. They predate this campaign and the plan leaves them as they
  are. No page written in this campaign has a personal name.
- `tutorials/optimization/` and `RELEASING.md` have `bash` blocks that no
  checker runs, and the numbered scripts of `examples/mtdna/` are not run
  either. All three predate this campaign and are outside the guide.
- README, "Why CRN-craft": the bullets carry no links by design, so the
  mtDNA model is named before the reader can know it.

## The researcher's five requests

| Request | Evidence |
| --- | --- |
| One guide | `docs/guide/`, 23 pages in the sidebar; `docs/tutorial.md` and `docs/mtdna-user-manual.md` are tables from each old section to its new place |
| Simple model first, complicated cases last | chapters 1 to 4 use two builder models and no shell, file or hash; attributes, model files and the stored forms are chapter 14 and page 14a |
| Models built with `add_species` and `add_reaction` everywhere; files as the late special case | a search of the README, the landing page and chapters 1 to 13: no `crn.load`, `load_example` or model dictionary. Chapter 10 and page 12e write their file from the built model; chapter 11 reads other formats, which is its subject |
| No engine history up front; a working example first | README and landing page: one sentence, install, the example with its output; the history is one closing sentence. The guide names the old engines in Migrating only (after finding 10) |
| An early, fair comparison with GillesPy2 and rebop that brings out ensembles and compiled models | `docs/comparison.md`, third sidebar entry, linked from the README, the landing page, chapter 1 and the how-to page. It says where the others are better (rebop's macro was faster in the one measured optimization; selectors that lose on W4), and "Two things to try first" are ensembles and compiled models |

## Gates

| Gate | Result | Evidence |
| --- | --- | --- |
| G1 every block that ran at 0.7.0 still runs | holds | the ledger of the table of contents (61 blocks) against the final `check_guide.py` report: 61 of 61 have a destination whose blocks ran. 110 runnable blocks run now |
| G2 frozen manual | holds | `git diff a43f7ca -- validation/compat/manual-0.6.0.md` is empty; 22 of 22 blocks, normal and `-O`, reports equal |
| G3 links and old paths | holds | `check_links.py`: 957 links, 36 stubs, 53 SUMMARY entries, 0 broken. All 144 paths under `docs/` at 0.7.0 exist |
| G4 site | holds | mdBook 0.5.0 build without warnings; `check_links.py --site`: 0 broken |
| G5 README links | holds | seven documentation links, the banner and the link to the repository overview |
| G6 names, spelling, unchecked blocks | holds after findings 10 and 13 | searches described above; every `python` and `bash` block of the guide, the READMEs and the landing page is in the checker's set; the reference pages have `json` and `text` only (the five `json` blocks parse, and the one model passes `crn.check`) |
| G7 no dead end in the walkthrough | holds after finding 1 | step 1 |

## The comparison page

Five claims about other tools, drawn at random from the 34 on the page
(seed 20261113), against the source distributions of GillesPy2 1.8.3 and
rebop 0.10.0 from PyPI: GillesPy2's model classes and its two kinds of
rate (`core/reaction.py`); SBML import and export through libSBML, and
JSON (`core/model.py`, `sbml/`, `core/jsonify.py`); rate rules and
assignment rules in the hybrid solver (`get_supported_features` of
`tau_hybrid_solver.py`); Python as its only language; rebop's Python API
needs no compilation and its macro supports mass action only (README).
All five hold. Also confirmed from PyPI: the versions, dates and licenses
of GillesPy2, rebop and StochPy, and CRN-craft's first release on
4 November 2026.

Every speed number was recomputed from its file: the GillesPy2 ratios and
the 28 agreement tests (`gillespy2-comparison-0.7.json`), the rebop ratios
and 8 tests (`rebop-comparison-0.7.json`), 3.44, 1.72, 0.94 to 1.00, 2.5%
and 15 times (parameter optimization page), 2.14 and 2.84
(`manual-ratios/timings.json`), the twelve compiled networks
(`RESULTS.md`), 2.1 and 5.8 to 6.1 (performance analysis). The numbers
hold; three statements around them did not (finding 4). The builder
ratios of the guide (3.6, 1.4, 2.6, 0.97) equal `builder-speed.json`.

## E1

**Certificate.** Read: `Reader::whole`, `Reader::stop`, `builder_form` and
`certify` in `rules/mtdna.rs`. Each propensity is split into factors; the
count must be the reactant's own count, the vacancy `max(0, K - N)` once,
the denominator one linear form without a constant, every coefficient
exact; laws must agree across the eight channels; the shared checks of the
file form follow. `max` has two arguments in the parser, so indexing them
is safe. Tested from Python on the wheel (`d9/e1/nearmiss.py`): 62
near-misses (7 parameter sets out of range, 15 altered vacancy terms and
denominators, 8 altered single rules, 6 altered deaths, 9 altered
topologies, a continuous variable, 16 altered stops) all fall back with
`CRN-R051` or are refused; none ran the fast path. One model first listed
as a near-miss, a replication written `sw * (2 * M)`, is accepted, and
rightly: it is the same structure with the weight 2 `sw`.

**Law.** `d9/e1/agree.py`: the builder form under `mtdna_direct` against
the same model under `static_direct`, 100,000 cells each, in seven
regimes that the declared check does not cover (extinction at `K = 5` and
`K = 1`, a stop, `mu = 1`, `mu = 0`, initial mutants, numeric death
rates): 95 comparisons of means, stopped fractions, stop times and event
counts, largest |z| 2.05. At `K = 1000`, 4,000 cells each: takeover
fractions 0.1075 and 0.1080. Trackers, an event log, a callback, kept
paths and `bind_parameters` run on the fast path; a bound candidate
outside the structure falls back. The declared check was run again: 120
comparisons, 0 failures, every model entry equal to
`laws/e1-builder-mtdna_direct.json`.

**Byte identity.** `d9/e1/identity.py`: 390 cases on this wheel and on
crn-craft 0.7.0 from PyPI: nine builder models and 40 random expression
models with shared subexpressions and nine kinds of stop, under Static
Direct with every selector and inverse waits, next reaction, tau, the tree
method, binomial, hybrid, trackers with predicates, event logs, bound
parameters, FSP, generated loop sources, and every model file of
`models/rules/` under six methods. A SHA-256 over metadata, records, rule
events, paths and observations, or over the error: all 390 lines equal
(361 results, 95 of them with failed replicates, and 29 errors). The
stored requests replay (2 and 9 tests). The full law gate with the release
CLI: every model entry equal to `validation/merge/laws/m8b-cli.json`,
4,814 comparisons, 0 failures, largest |z| 3.88.

**Tests and the mutation gate.** `cargo test --release --workspace
--locked`: 407 passed. The 19 additional mutants and the nine that joined
the gate are finding 3. The whole gate was then run once: 67 of 67
mutants caught.

## Checks on the final tree

Run on `8fd6f6f` (every fix but this log and the history index, which
change no checked block), with the wheel of this checkout in the campaign
venv and the release CLI, each heavy step under the CPU lock (runner and
logs: `/home/claude/docs071/d9/checks/`):

| Check | Result |
| --- | --- |
| `check_guide.py --binary target/release/crn-craft --require-cargo` | passed: 26 files, 110 runnable blocks ran (97 `python`, 13 `bash`), 74 outputs compared, 6 `sh` syntax-checked, 3 fitting blocks in smoke mode, 3 `r` blocks and their output not run |
| the same under `python -O` | passed; report byte-identical |
| frozen manual, `check_manual.py --manual validation/compat/manual-0.6.0.md --require-cargo`, normal and `-O` | 22 of 22; reports identical |
| `check_links.py` | 966 links, 36 stubs, 53 SUMMARY entries, 0 broken (with this log and the index entry) |
| `make_python_reference.py --check` | up to date |
| `mdbook build` (0.5.0) and `check_links.py --site` | no warnings; 3,574 references, 0 broken |
| `check_rust_docs.py` | 2 rust programs compiled and ran, 0 failed |
| `ruff check` and `ruff format --check` on `validation python examples` | clean (124 files) |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | clean |
| `cargo test --release --workspace --locked` | 407 passed, 0 failed |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | clean |
| `check_installed.py` from `/tmp` with `python -I`, normal and `-O`, in the fresh venv of step 1 | 18 checks passed; reports identical |
| `cargo package --list` of `crn-craft`, `crn-craft-cli`, `crn-craft-rules` against a worktree at `a43f7ca` | 56, 9 and 8 files; the three lists are identical to 0.7.0 |
| law gate, `check_laws.py --binary`, `static_direct`, twelve models | every model entry equal to `laws/m8b-cli.json`; 4,814 comparisons, 0 failures |
| E1 law check, builder form under `mtdna_direct` | 120 comparisons, 0 failures; equal to the committed report |
| `mutation_check.py --lock`, all 67 mutants (on `a321f71` and later; the engine sources are unchanged since) | 67 of 67 caught, complete |

## For the owner: design matters, not changed

- **D1. The default method.** `crn.ensemble` without `method` runs the
  tree method on a builder model, slower and without a hint, so the guide
  has to say "name `static_direct` in every call". A default that picks
  Static Direct for a closed model would remove the sentence; it changes
  results of existing calls, so it needs a new request version or a
  warning first.
- **D2. Two hashes for one model.** `crn.ensemble(model)` hashes compact
  bytes, `model.to_json()` writes indented ones, so a Python run and a
  command-line run of "the same model" name it differently (chapter 14
  explains it). A `to_json` that can write the hashed bytes would let the
  two agree.
- **D3. Replaying an archive takes a manual copy,** because the archive
  stores `model.json` and the request keeps its path; requests written by
  `crn.ensemble` name the model `model.rules.json`. A `crn-craft replay
  DIR`, or `run` accepting the archive's own file name, would make
  finding 1 impossible.
- **D4. Program output still speaks of rules and engines:** "rules request
  does not identify the supplied model", "generated rules loops require
  Static Direct", "random product attributes require the dynamic engine",
  the progress line `crn_craft.rules.compile_loop: ...`, and
  `docs/rules-language.md` in diagnostics. Changing them changes outputs.
- **D5. Two vocabularies in one API:** the builder says species and
  reactions; results say `rule_events`, `rule_names`, `event["rule"]`,
  `observables` and `variable`.
- **D6. Finding 2 as a policy question.** The fast path for the builder
  form changes the result of an unchanged `run/1.0` request between 0.7.0
  and 0.7.1. It is documented now. The alternatives are an option that
  asks for the builder certificate, or leaving it as it is and naming it in
  the release notes.
- **D7. The builder cannot declare a sum** such as "all mutants", so
  statistics of a total need kept paths (chapter 12 says so), and it
  cannot write attributes (chapter 14).

For the release preparation: the wheel and the pages still say 0.7.0
while three pages say "since release 0.7.1" (the version line of the
README, the output of chapter 1, the header of the Python reference and
the comparison page follow the version bump); the PyPI README and
`pyproject.toml` link the site, which exists once GitHub Pages is enabled;
a compiled loop builds against the crate of the wheel's version, so the
crate has to be on crates.io before the wheel is used.

## Deviations

1. The table of contents and the plan were read only after step 1, as the
   task asked; `SUMMARY.md` was read in step 1 as the site's sidebar.
2. The five claims were checked against the source distributions on PyPI,
   not against the tools' web pages.
3. Pages 12 to 12e were not edited before the FP merge, as the coordinator
   asked; findings 6, 7 and 8 on them were fixed after it.
4. One new runnable block (chapter 10, Replay); `guide-blocks.json` counts
   16 for that page.
