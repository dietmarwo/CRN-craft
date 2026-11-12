# Final review of the 0.7.0 engine merge

The review before the pull request `merge-0.7` → `main` (the plan's "one
Opus review before the PR"), on the integration checkout, branch
`merge-0.7`, from `1d86b12` (all tasks M0–M10 merged). The reviewer wrote
none of the reviewed work. Evidence of every run: `/home/claude/m07/review/`
(`gates/`, `ci/` for the state at `1d86b12`; `final/` and `ci2/` for the
state after the fixes; `manual-ratios/`).

## Commits

| Commit | Title |
| --- | --- |
| `a203b68` | Review: Declare the re-measurement of the manual's speed ratios |
| `f286215` | Review: Describe the one engine in the crate READMEs and facade docs |
| `a947cb0` | Review: Name the one engine and the canonical commands in messages and API docs |
| `3bcb038` | Review: Match the diagnostics reference to the codes the code emits |
| `1093420` | Review: Correct documentation that disagrees with the code |
| `a4f824c` | Review: Correct claims of the CHANGELOG's Unreleased section |
| `91f040c` | Review: Mark the M1b page and the model cards as history |
| `205d776` | Review: Re-measure the manual's speed ratios on its own workload |
| `585561d` | Review: Drop the deleted tau test from the mutation gate's test list |
| `1a99975` | Review: Update the test, line and mutation counts of the technical report |
| (this log) | Review: Add the review log |

## Findings

No blocker was found: every gate passes on the merged tree, before and after
the fixes (section "Gates"). Severity: **blocker** (must not merge),
**should-fix** (wrong or misleading, fixed here unless stated), **note**.

### 1. Whole-docs consistency

Searched every markdown file outside `docs/campaign/` and the
`validation/*/RESULTS.md` files for two-engine wording, removed commands,
removed Python names, removed contracts, `fast`, `cached_direct`,
`direct_cr`, `sorting_direct` and `xoshiro`. Hits in the CHANGELOG's
migration tables, the migration tables of `alternative-ssa.md`,
`reusable-parameters.md`, `network-api.md`, `contracts.md` and the manual,
the dated decision records and iteration logs, the comparison pages with
their dated 0.6.0 notes, the history sections of the technical report and
`comparison-gillespy2.md`, and the evidence files under `validation/` are
history or migration and stay. Claims about what is read, written or
deleted were run against the release CLI and the installed wheel.

| # | Severity | Finding | Fix |
| --- | --- | --- | --- |
| D1 | should-fix | `crates/crn-craft/README.md` described the network engine (`cached_direct`, `ssa-run/0.3`, the `fast` policy, parameter templates) and sent attributed models to `crn-craft-rules`; the CLI README called the engine "the experimental rules engine"; the facade's README and crate docs presented it as an engine next to "the network engine in `crn_craft`" and said `compile_loop` builds against it. | `f286215`: the three READMEs describe the one engine; the facade's crate docs and package description call it the deprecated 0.6.0 name of `crn_craft::rules` (the 0.6.0 quickstart stays as its doctest). |
| D2 | should-fix | `CRN-E032` (`crn-craft convert` given a `model/1.0` or rules document) said "`crn-craft.model/1.0` is a rules-engine contract" and named `crn-craft rules-check` and `crn-craft rules`. The FSP refusal hint and the API docs of `run_files`, `solve_files`, `check_model_file` and `archive::v1` named the alias `crn-craft rules`; the crate docs counted thinning as approximate (`describe()` lists it as exact); the `rules` module docs said "the rules engine" and "zero, one or two reactant patterns". | `a947cb0`; codes and paths unchanged (the tests check both). |
| D3 | should-fix | `docs/diagnostics.md` was written from the plan's list (M10b log), not from the code: `CRN-E022`, `E042`, `E061`, `E062` and `E064` were listed as current but no code emits them; `CRN-E014`, `E032` and `W001`–`W003` were listed as retired but the network model reader, `convert` and `analyze` emit them; `CRN-E044` is still emitted by `crn-diffsol`. | `3bcb038`: each code placed by a grep of `crates`, `bindings`, `python` and the integrations. |
| D4 | should-fix | `docs/direct-ssa.md` gave the inverse wait as `-ln(u)/total` (known item). | `1093420`: `-ln(1 - u)/total` with `u` redrawn when 0, as `rules/flat.rs` and `flat_path.rs` do. `rules-language.md`'s `E = −ln u` for `hybrid` is right (`hybrid.rs`). |
| D5 | should-fix | Five pages (CHANGELOG, manual, model reference, tutorial, `python/README.md`) said `crn_craft.rules` is "the same module" as the top level. It is not: `crn_craft is crn_craft.rules` is false and `Model`, `ensemble`, `prepare`, `compile_loop`, `describe` differ (0.6.0 defaults against 1.0 defaults). | `1093420`, `a4f824c`. |
| D6 | should-fix | `docs/contracts.md` listed `state/0.1` among the replaced contracts (known item). M8b kept `state` (log, "Deviations"), `analyze MODEL STATE` reads it, and `describe` lists it. | `1093420`: the row no longer lists it. `architecture.md` already matched M8b. |
| D7 | should-fix | `docs/result-contract.md` did not describe `started` or the pre-cancel case (M8b note). Run: a pre-cancelled ensemble of 8 replicates on 2 threads recorded replicate 0 as started and `cancelled` with 0 events and 0 words, the others as not started; `boundaries_v1::unstarted_replicates_are_accounted_for` pins at most one started replicate per worker. | `1093420`: the `replicates` row says so. The archive pages and the Rust docs of `Control` and `Record` were already consistent; the tutorial's pre-cancelled archive verifies. |
| D8 | should-fix | `docs/reference-export.md` said Hill and Michaelis–Menten laws and orders above two are refused with `CRN-E050`; since M5f they export as custom propensities (`crn-craft export models/examples/michaelis-menten.model.json gillespy2` succeeds). | `1093420`, from the M5f log and `check_export.py`. |
| D9 | should-fix | README: thinning listed as approximate, "five selectors" (the crate docs say four, `auto` chooses among them). `python/README.md` listed ODE right-hand sides among the wheel's features (the wheel has no ODE API). C4's title and the worker note of `validation/parameter-optimization/README.md` said "the rules engine"; the optimization tutorial named `crn_craft_rules::execute`. | `1093420`. |
| D10 | should-fix | `docs/m1b-analysis.md`, linked from the README as current, described analysis 0.2, `analyze-v1`, `CRN-E042` and `analyze_linkage(network, …)`; the Direct and rate-law model cards described the network engine's codes. | `91f040c`: dated history notes, the linkage signature corrected. |
| D11 | should-fix | `docs/technical-report.md`: 538 Rust tests (before M8b) and the line count only by reference (known item). | `1a99975`: 397 tests (29 doctests) by a fresh `cargo test`; 29,372 source lines after the review (29,482 at `1d86b12`, as in the M8b log, before the facade docs were shortened); the mutation gate from the run below. |
| D12 | note | The crate README dependency snippets still say `"0.6"`, and the `--version` screenshot of the tutorial shows 0.6.0. | Left to M11 (version bump, snippets `"0.7"`, `01-help` rendered again), as the plan says. |
| D13 | note | No personal names in the documentation written during the merge (the diff since `43105b2`, searched for the names of people; host paths such as `/home/claude/m07/…` name evidence locations). | none |

### 2. Known review items

| Item | Result |
| --- | --- |
| `architecture.md`: module names, converter path, archive versions of `verify`, `state/0.1` | Correct against the tree: every named module exists, `crn_craft::rules::convert::from_network` is public, `verify` reads `archive/1.0` only (refuses 0.1 with `CRN-E065`, checked by `check_rules_python.py`), `analyze MODEL STATE` reads `state/0.1`. `contracts.md` fixed (D6). |
| `technical-report.md` counts | D11. |
| `direct-ssa.md` inverse wait | D4. |
| Manual ratios (sections 3, 8, 16) | Re-measured below; prose updated (`205d776`). |
| CHANGELOG `[Unreleased]` and README | Below. |
| M8b notes | RELEASING.md describes no removed engine (it lists the three crates, which stay); crate READMEs, facade docs and `python/README.md`: D1, D9. `validation/merge/kinetics-reference` says in its header that it builds only against the 0.6.0 crate: fine. Pre-cancel: D7. |
| Personal names | D13. |

**CHANGELOG `[Unreleased]` and README** (`a4f824c`, `1093420`). Every
Python name of the tables exists or is gone as stated (`hasattr` on the
installed wheel; `crn.rules.model_bytes` exists); every removed command is
a usage error (`CRN-E031`, exit 1); `check` of a network `run/0.1` answers
`check/1.0` with `CRN-E033`; `describe` equals `crn.describe()` and
`python -m crn_craft describe`; the quick-start blocks and commands of the
README run; the README's diagnostic example prints as shown; the speed
range 0.66–1.20 and the 1.55 noise figure are those of the M4d rerun;
`Kinetics::species`, `initial`, `parameter_values`, `sampling::poisson_with`
and the new struct fields exist. Corrected:

- the list of new contracts lacked `archive-check/1.0` and `report/1.0`;
- "limits 16× larger": types and rules grew 64 times, observables,
  parameters, tuples and channels 16 times;
- "the mtDNA manual's workload" was the three-class CME workload of M10c;
- the clock rule applies to Static Direct under `run/1.0`
  (`lenient_clock` in `flat_path.rs` and `loop_host.rs`), not to every method;
- the law-gate summary did not say that the first RSSA-CR run failed;
- the contract migration line named a `result/1` that never existed
  (`result-v1` was the schema of `result/0.1` and `0.2`), only the first
  `ssa-run` and `ssa-result` versions, and missed `analysis` and `describe`.

Note: "5–18 % slower" is the pooled figure of both M4d runs; run 2 alone
gives 6–20 %.

### Manual speed ratios

Declared in `validation/merge/RESULTS.md` (`a203b68`) before the first
timing, measured on the manual's section 3 example (200 cells, two
threads), medians of three alternating runs under the CPU lock; every
declared check passed (`validation/merge/manual-ratios/timings.json`).

| Ratio | Measured | Manual before |
| --- | ---: | ---: |
| Static Direct / `mtdna_direct` (wall time; per event 2.81) | 2.84 | 3.3 |
| Static Direct with inverse waits / `mtdna_direct` | 3.13 | 3.9 |
| interpreted Static Direct / compiled loop (from Python) | 2.14 | 2.4 |

Sections 3, 8 and 16 now say 2.8, 3.1 and 2.1 with the date; no code block
changed. The optimization tutorial and the performance analysis quoted the
same old ratios and now cite the new ones. Note: the compiled-loop figures
of section 8 (2.21× and 1.86× on the 1000-cell workload) and the 0.2.0
sections of the performance analysis were not re-measured; they are
labeled with their source.

### 3. Code review of the merge

| # | Severity | Finding | Action |
| --- | --- | --- | --- |
| C1 | should-fix | Messages and API docs naming the alias commands and two engines | D2 (`a947cb0`) |
| C2 | should-fix | `validation/mutation_check.py` listed the deleted test file `tau` (skipped silently) | `585561d` |
| C3 | should-fix, left | `crn-craft run` or `check` with a `model/1.0` model and a network `crn-craft.run/0.1` request reports `CRN-E003` "missing field `max_sampling_steps_per_replicate`" at `/budget`, not that `run/0.1` is no longer read. Correct code, unhelpful message. | For the owner: a dedicated message at `/schema` is a behavior change of an error path and wants its own test. |
| C4 | note | The mutation gate's test list does not include `network_models` (the converter's strict parsing), as M8b noted. Adding it changes the gate's test set. | Owner's decision. |
| C5 | note | `describe()["legacy"]["rules"]` keeps the 0.6.0 entry, with `"status": "experimental"`, `"commands": ["rules-check", "rules"]` and `contracts_1_0.archive.command` naming the alias `rules`; the `rules` module docs keep "Status: experimental". Aliases are not removed commands. | Whether the one engine is still experimental at 0.7.0 is the owner's call. |
| C6 | note | Public API: rustdoc with `-D warnings` and `#![warn(missing_docs)]` pass; every entry point named in the crate docs exists; every Python name the docs use as current exists. Dead code: `cargo clippy -- -W dead_code -W unused` and the bindings' clippy report nothing. The 15 files in `schemas/` equal `crn-craft describe KIND`, and `describe` lists each. Every file under `models/` is read by a test or script; every script in `validation/` and `validation/merge/` imports. | none |
| C7 | note | CI: every job runs checks with nonzero counts (below); no step passes by running nothing. CI has no rustdoc step; I5 rests on local runs. | Owner's decision. |
| C8 | note | The plan's rust matrix includes Rust 1.85.0 and 1.88.0; the proxy here refuses `static.rust-lang.org`, so neither toolchain could be installed. Stable clippy's `incompatible_msrv` lint (rust-version 1.85) passes and no let chains are used. | CI decides. |

Process note: a `--help` sweep over the validation scripts (to find scripts
that no longer import) ran some report scripts that ignore `--help` and
rewrote tracked figures and evidence (`docs/images/julia-*`,
`docs/comparison-*.md`, `validation/parameter-optimization/…`). They were
restored with `git checkout` before any commit and before the first suite
read them; no commit contains them. Every script in `validation/` and
`validation/merge/` imports.

## Gates

Rust (stable 1.97.0) and Python 3.12.3 in `/home/claude/fb2/engine-venv`
(NumPy 2.3.5, pandas 3.1.1, xarray 2026.9.0, jsonschema 4.26.0, SciPy
1.17.0, GillesPy2 1.8.3, rebop 0.10.0, ruff 0.11.13), DSMTS at `e79fd58`,
everything under the CPU lock. The CI-equivalent suite ran twice: at
`1d86b12` with the CLI and wheel of `1d86b12` (all checks passed and all
normal and `-O` report pairs were identical; the diff steps of that run had
a script error and were redone by `cmp`, and its last document checks read
files that were being edited), and at `585561d` with a fresh build, which is
the record (table). The commits after `585561d` change only the technical
report and this log; `check_links.py` passed again on them (434 links,
normally and under `-O`).

| Check | `585561d` |
| --- | --- |
| `cargo fmt --all -- --check`; clippy `-D warnings` (workspace and `bindings/python`) | pass |
| `cargo test --release --workspace --locked` | 397 passed, 0 failed, 30 result groups; I1 `stored_requests` 2, `static_direct_defaults` 9; I3 `generated_loop` 28, `generated_kernel` 3; facade 3 + doctest |
| rustdoc `-D warnings` | pass (only cargo's known output-filename collision note) |
| clippy `-W dead_code -W unused` | no warning |
| schemas against `describe KIND` | 15 of 15 equal |
| ruff format (`validation`, and `.` with `--exclude '*.md'`), ruff check | pass |
| `check_dependencies.py` | pass, reports identical |
| `check_contracts.py` | 139 of 139 |
| `check_network.py`, `check_structure.py` | 219 of 219; 309 of 309 |
| `reference_rng.py` | 9 streams, 360 words |
| `check_ssa.py`, `check_cme.py` | 17 of 17; 14 of 14 |
| `check_export.py` | 9 of 9, 15,756 rate comparisons |
| `check_fsp.py` (DSMTS) | 16 rows, largest error 4.1e-14 |
| `check_tau.py` | pass, 30,000 paths |
| `check_rules.py`, `check_mtdna.py` | 47 of 47; 151 checks, 32,768 paths |
| `check_python.py`, `check_rules_python.py` | 17 checks; 27 checks |
| `check_compiled_loop.py --require-cargo` | 59 rows |
| I2 `check_manual.py --require-cargo` | 26 of 26 blocks |
| I2b `check_manual.py --manual validation/compat/manual-0.6.0.md` | 22 of 22 blocks |
| `check_manual.py --manual python/README.md` | 2 of 2 blocks |
| `check_links.py` | 434 links, 0 broken |
| `check_builder.py`, `check_hybrid.py` | 18 tests; 11 checks + 10 independent |
| `check_alternative_contracts.py`, `check_alternative_ssa.py` | 9 checks; 168 of 168 |
| `check_importers.py`, `check_tutorial.py` | 49 rows; 59 blocks |
| `check_installed.py` (`python -I`, outside the checkout) | 18 checks |
| every Python check under `python -O` | reports identical to the normal run |
| line-endings job (fresh clone with `core.autocrlf=true`) | no file without LF; 192 tests of the job's targets pass |
| `mutation_check.py` (CI's `contracts` job) | 47 of 47 mutants caught, complete, unmutated baseline passed |

## For the owner and M11

- C3, C4, C5 and C7 above.
- M11: version 0.7.0 everywhere, crate README snippets `"0.7"`, the
  tutorial's `01-help` screenshot, the facade's dependencies and lockfile
  entries (M8c note), and the Rust 1.85.0 and 1.88.0 jobs of CI.
