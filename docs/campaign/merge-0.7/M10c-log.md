# M10c log: architecture and reports

Task M10c of [the plan](../../history/plan-0.7.0-engine-merge.md): the architecture,
agent and report documents describe the one engine, and the comparison
documents say what they measured. Branch `merge/M10c`, from `merge/M10b` at
`d419347`. The documents are written for the state after M8b (network
engine, `xoshiro`, the `fast` variants, the old contracts, the hidden
commands and `crn_craft._legacy` deleted), which ran in parallel and is not
in this branch.

## What changed

| Document | Change |
| --- | --- |
| `docs/architecture.md` | Rewritten: one engine in `crn-craft` (`crn_craft::rules`, the flat view and the features on it, selectors, the converter, archives), the facade crate, the CLI and the Python package with `crn_craft.rules`; the contracts table (1.0 contracts, the `rules-*` contracts still read, what is no longer read); preparation and methods, Static Direct and its selectors, random streams, parallelism, budgets and cancellation, generated code. Diagrams in mermaid, as before. |
| `docs/ai-context.md`, `AGENTS.md` | What an agent needs to use and extend the one engine: crate and module, contracts, methods and selectors, Python names, where the tests, stored requests, fixtures and gates are, and the rules for changes. |
| `docs/technical-report.md` | Rewritten for 0.7.0: methods of the one engine, the law gate (N = 4,814, every selector and `auto`, the failed first RSSA-CR run kept), the reference gates of FSP, analysis, kinetics, export, higher-order rules, tau and quantiles, the speed gate with its noise note, the mtDNA workload table, parameter workers (failed `binding12` row kept), tau on large networks, GillesPy2. The 0.6.0 numbers are in a section marked as history. |
| `docs/extending.md` | One crate, 1.0-only fields, the flat view, the stored requests and fixtures under `crates/crn-craft/tests/`, the law gate. |
| `docs/specialized-kernels.md` | Test paths `crates/crn-craft/tests/generated_{kernel,loop}.rs`, 19 loop fixtures (higher-order and clock-rule fixtures added to the list), `compile_loop` builds against `crn-craft`, the generator examples live in the facade crate, the deferred-projection bound of `model/1.0`, selectors in the checklist, and a note that the timings predate 0.7.0's faster Static Direct loop. |
| `docs/comparison-gillespy2.md` | Re-measured on the unified engine (below); the 2 November measurement of the 0.2.0 network engine kept in a separate history section with its own table and figure. |
| `docs/comparison-julia-8-threads.md`, `comparison-julia-gap10.md`, `comparison-catalyst-jumpprocesses.md`, `comparison-catalyst-jumpprocesses-before-acceleration.md` | A dated note at the top: measured with the 0.6.0 network engine, which 0.7.0 replaced; Julia is not available here; pointer to the speed gate. Tables unchanged. |
| `docs/comparison-parameter-optimization.md` | The same note, plus the M4d worker rerun (`validation/merge/parameter-workers-m4d.json`) with its failed `binding12` `reuse` row (1.55 against the 1.25 target). |
| `docs/mtdna-performance-analysis.md` | New section "Release 0.7.0" with the 0.6.0 against 0.7.0 table; the older sections are marked as release 0.2.0 measurements; CLI, `compile_loop` and parameter-binding wording of the one engine. |
| `docs/mtdna-rules.md` | One engine; the Python blocks use `import crn_craft as crn`; `model/1.0` and `run/1.0` in the version history; higher-order rules per method; lowering limits of both model families; `selector` in the Static Direct options and in the Rust struct literal; dated timings. |
| `validation/merge/RESULTS.md` | The GillesPy2 declaration (committed before the run) and results; the mtDNA workload results, marked as measured without a declared target. |
| `validation/merge/mtdna-workload/` | `bench.py`, `run.sh`, `timings.jsonl`, `summary.json`, copied unchanged from `/home/claude/m07/mtdna-speed/`. |
| `validation/merge/gillespy2-comparison-0.7.json`, `docs/images/gillespy2-time-per-event-0.7.png` | The merged GillesPy2 report and its figure. |

## Measurements

### GillesPy2 (declared)

Declared in `RESULTS.md` (`cec91ed`) before the first run: the script of
M8a unchanged, the CLI and wheel of this branch, four parts (`--only W1` to
`W4`) under the CPU lock, the script's repeats, and an agreement rule
(Bonferroni over 28 tests at family-wise 0.01: a p-value below 3.6e-4
fails; W4 species stay equal). Environment: `/home/claude/m07/venv-m10c`
(the package set of `venv-m6` without `crn-craft`: GillesPy2 1.8.3, NumPy
2.3.5, SciPy 1.18.1, matplotlib 3.11.2, SCons 4.11.1, Python 3.12.3), g++
13.3, `CARGO_TARGET_DIR` inside the worktree.

Each part ran once and completed (finished 13:11, 13:13, 13:16, 13:20 on
2026-11-12); the merge checked that the non-workload fields of the four
parts are equal. **Agreement: 28 of 28 passed** (smallest KS p-value 0.209,
smallest chi-square p-value 0.044; W4 species equal in both simulators).
Times (no target): GillesPy2's C++ solver took 2.9, 10.6, 3.6 and 9.3 times
as long as CRN-craft's default request on W1 to W4, and 2.8, 5.6, 3.7 and
1.3 times as long as the linear scan with inverse waits; the NumPy solver
134 to 358 times the default; C++ compile 2.1 to 2.9 s. Composition–rejection
and sorting are slower than GillesPy2 on W4 (0.88 and 0.93).

### mtDNA workload (no declared target)

Measured by the coordinator; copied unchanged. Recomputing the medians from
`timings.jsonl` gives the six rows of the table in the task prompt, and every
row of both versions has one digest and one sampler per method and thread
count. Both wheels report the version string 0.6.0; the `python` field of a
row names its environment (0.7.0: the integration wheel with the Rust
sources of `21f71b8`).

## Checks

| Check | Result |
| --- | --- |
| `validation/check_links.py` (relative links and heading fragments in `docs/` and the READMEs) | 429 links, 0 broken |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass (the copied `bench.py` included) |
| The two Python blocks of `docs/mtdna-rules.md` after the switch to `crn` | run with 20 replicates: Static Direct sampler `chacha8-u53-ziggurat/rules-static-direct/1`, `mtdna_direct` sampler as stated, `bind_parameters` works |
| `auto` on the three-class models | picks the linear scan (sampler without selector suffix), as the documents say |
| Rust standard set | not run: no Rust, binding or package code changed |

## Deviations

- The GillesPy2 figure is a new file (`gillespy2-time-per-event-0.7.png`);
  the old figure stays with the history section.
- The generator examples stay in the facade crate (the M8c log offered to
  move them); the documents say where they are.
- `docs/mtdna-rules.md` switched its Python blocks to `import crn_craft as
  crn`, as the manual did in M10b.

## For the final review: details that depend on M8b

- `architecture.md`: module names (`fsp.rs`, `analysis.rs`, `cone.rs`,
  `conservation.rs`, `linkage.rs`, `kinetics.rs`, `export.rs`, `sampling.rs`,
  `strict_json`, `diagnostic`, `unified`, `archive`), the converter path
  `crn_craft::rules::convert::from_network` with the network validation
  inside it, `verify` reading `archive/1.0` only, `state/0.1` still read by
  `analyze`, and the list of contracts no longer read ("their schema files
  are gone"). No file name is given for the next-reaction clock heap.
- `technical-report.md`: the Rust test count (538 at `d419347`, before M8b
  removes the network tests), the source line count (pointed to the M8b
  log), the mutation gate (47 of 47 per M8a).
- `specialized-kernels.md`: 19 loop fixtures and 1 kernel.
- `comparison-gillespy2.md`: measured on a tree that still had the network
  code. `validation/benchmark_campaign.py` reads
  `models/examples/repressilator.model.json` as a network model; after M8b
  it must read `models/legacy/`, or the reproduce command breaks.

## Disagreements with documents owned by others

- `docs/direct-ssa.md` gives the inverse wait as `-ln(u) / total`; the code
  (`rules/flat.rs`) and the M4f and M8a logs use `-ln(1 - u) / total` with
  `u` redrawn when 0. The technical report follows the code.
- The draft CHANGELOG (`/home/claude/m07/m10/changelog-unreleased.md`,
  "Added", archives) says "`crn-craft verify` reads both archive versions";
  `docs/contracts.md`, `docs/archive-contract.md` and the plan say
  `archive/0.1` is no longer read after M8b. The architecture follows the
  latter.
- `docs/mtdna-user-manual.md` sections 3, 8 and 16 give `mtdna_direct` as
  about 3.3 times faster than Static Direct and the generated loop as about
  2.4 times faster than interpreted Static Direct. Both ratios were measured
  before 0.7.0's faster Static Direct loop; on the mtDNA workload table the
  per-cell ratio is about 2.3 (another workload). The performance analysis
  and the technical report say so; the manual does not.
- `docs/contracts.md` keeps `crn-craft.state/0.1`, which the plan's M8b list
  deletes; the M5b and M7 logs keep it for `analyze MODEL STATE`.

## Migration notes

None: documentation and evidence only.
