# D5 log: how-to page, landing page and READMEs

Date: 13 November 2026. Branch `docs-0.7.1`, main checkout. Inputs: the brief, sections 5, 7 and 9 of the table of contents (D1), every page of `docs/guide/`, `docs/comparison.md`.

## Produced

- `docs/howto.md`: the fifty tasks of section 5 in its eleven groups, one line each, 75 links. Every task has an answering section; no task was left without a link.
- `docs/index.md`: the landing page of section 7.2: sentence, install, the example with its output, "Why CRN-craft" with one link per line, the comparison line, "Where to start", an eight-row "What is included" table (section 4.9), "The rest of the site", the history sentence and the related-repositories line. 90 lines against a budget of about 55: the eight-row table of section 4.9 and the wrapped lines account for the difference. Checked standalone.
- `README.md`: rewritten to section 7.1 (80 lines; budget about 75). First screen: title, banner, the related-repositories line, one sentence, install, the example and its output. Seven documentation links (comparison, site, guide, how-to, case study, reference, changelog); `LICENSE` and `THIRD_PARTY_NOTICES.md` named, not linked. Checked standalone (`<!-- check: standalone -->`).
- `python/README.md`: one sentence on what it is, `sh` install block (it was `bash`, which the chapter driver would run), the two examples unchanged, the first item of "What is included" without the history, the merge paragraph cut to one sentence at the end, seven documentation links to the site (the old absolute link to `docs/using-from-r.md` is gone). Checked standalone.
- Crate READMEs (`crates/crn-craft`, `crn-craft-cli`, `crn-craft-rules`): links to the site, the user guide, and where it fits the Rust API page or chapter 10. In `crn-craft` the dependency block moved up to the first screen.

## Content of the old README and where it is now

| Old section | Now |
| --- | --- |
| Version paragraph, "start with" | one sentence in the README; the history is one sentence at the end |
| One engine (format, JSON rule, method table, reproducibility) | `reference/model-language.md`, chapter 7 "Which method for which model", chapter 3 and 8, Migrating (unchanged, already there) |
| Optimization tutorials | chapter 13 "Beyond this example" links the tutorials; `tutorials/optimization/README.md` links the shared catalog |
| Documentation lists: start here, using, model cards | `docs/SUMMARY.md` (every target is a sidebar entry) |
| Documentation lists: campaign, plans and decisions, implementation logs and validation results | **moved** to `docs/history/index.md`, section "Lists by campaign" (links rebased) |
| Quick start (packaged example, dimerization, command list, diagnostic JSON) | chapters 2, 10 and 14; `reference/diagnostics.md` |
| What is included (table) | `comparison.md` feature table; eight rows in `docs/index.md`; evidence numbers in `evidence/technical-report.md`; strict input in chapter 2 and `reference/contracts.md` (LF and `core.autocrlf`) |
| Release history | `CHANGELOG.md` (0.2.0 to 0.7.0) |
| Developer checks | **moved** to `docs/developer/extending.md`, new section "Checks" (fence changed from `bash` to `sh`) |
| Campaign and remaining validation | **moved**: open validation work and the ecosystem boundary to the end of "Checks" in `extending.md`; the campaign pointer to `docs/history/index.md` |

## Checks

- How-to anchors: `/home/claude/docs071/d5/verify_howto.py` (outside the repository) parses the 50 task lines and checks each of the 75 links: the file exists, the fragment is the slug of a real heading, and the link text equals that heading. 50 tasks, 75 links, 0 errors; two link texts are deliberate short forms of a page title ("12b, Mutant classes", "15, Using CRN-craft from R").
- `validation/check_links.py`: 825 relative links, 36 stubs, 53 SUMMARY entries, 0 broken (before D5: 5 broken, all `index.md` and `howto.md`).
- Chapter driver on the 0.7.0 wheel (venv-m10a), normal and `python -O`: `README.md` 1 block, `docs/index.md` 1 block, `python/README.md` 2 blocks; all outputs equal. `validation/check_manual.py --manual python/README.md` (the CI step today): 2 of 2 passed.
- The README's example and output are byte-equal to those of `docs/index.md` (and its means equal those of `python/README.md`, 12.623 and 17.305).
- `sh` blocks of the crate READMEs pass `bash -n`.
- mdBook 0.5.0 build: no warnings; `check_links.py --site`: 3,512 references, 0 broken. Screenshots of the landing page, the how-to page and a pandoc rendering of the README in `/home/claude/docs071/d5/`.

## Deviations and open points

- The crate READMEs hold `toml` and `rust` blocks that no checker runs, before and after this task; the `rust` block of `crn-craft-rules` is a fragment (it uses `model_bytes` without defining it). The chapter driver accepts neither fence. Proposal for D8: compile the crate READMEs' `rust` blocks like `check_rust_docs.py`, or mark them as fragments.
- The site links in `python/README.md` and the crate READMEs point at `https://dietmarwo.github.io/CRN-craft/`, which exists only after the owner enables GitHub Pages.
- D8 must add the rule that the example of `docs/index.md` equals the README's (section 7.2); today it is checked by hand in this log.
- "Why CRN-craft" quotes the compiled mtDNA model at 2.1 times (12 November 2026), which is the file form of the model; chapter 5 says so. The builder form compiles to a different ratio (D4 finding).
- `docs/campaign/README.md` still links `../tutorial.md` (a stub); `docs/campaign/` is not edited by this campaign and not link-checked.
