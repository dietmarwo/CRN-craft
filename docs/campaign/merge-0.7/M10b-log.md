# M10b: manual and references

Branch `merge/M10b`, from `merge-0.7` at 342c322 (M8a, M8c and M9 merged).
The references are written for the state after M8b (network engine deleted,
as the plan's section 3 and M8b list it). M8b was not merged here, so the
statements about removed contracts, retired codes and `crn_craft._legacy`
could not be run against the final tree. See "Open points".

## What changed

1. **I2b first** (commit `eddadc8`).
   - `validation/compat/manual-0.6.0.md` is `docs/mtdna-user-manual.md` at the
     tag `v0.6.0` (the tag exists; `cmp` against the manual of the branch
     point showed no difference). A short `README.md` states its origin.
   - `validation/check_manual.py` accepts fences of more than three backticks
     (a block ends at a fence of at least as many) and caps both
     `crn.ensemble` and `crn.rules.ensemble` under `--max-replicates`.
     `--manual PATH` already existed.
   - CI runs the frozen copy normally and under `-O` and diffs the two
     reports.
   - Result: all 22 blocks of the frozen copy print the output shown, normally
     and under `-O`, with identical reports. No block, including the two with
     a cache path or a timing, needed the "variable" allowance beyond the
     checker's normalization of timings.
2. **Manual** (`docs/mtdna-user-manual.md`).
   - All blocks use `import crn_craft as crn`. The introduction says once that
     `crn_craft.rules` is the same module under its 0.6.0 name; section 18 has
     a table of the differences between the two spellings.
   - New short sections: quantiles and the median (section 3), saving results
     (section 3, with a four-backtick output block that shows part of a
     markdown report), selectors for Static Direct with the rule that
     `crn_craft.rules` keeps the linear scan (section 8), and larger models via
     `model/1.0` (section 12). The CLI section uses `crn-craft.run/1.0`,
     `check`, `run`, `--archive` and `verify`. The generated-loop section names
     `crn-craft` and the one-time rebuild of 0.6.0 caches.
   - Expected outputs were produced by running the blocks
     (`/home/claude/m07/m10b/tools/update_outputs.py`): only blocks whose
     normalized output differed were rewritten. Every pre-existing block kept
     its numbers under the top-level names except two: the Gompertz block
     prints `crn-craft.model/1.0` instead of `rules-model/0.4`, and the CLI
     block's result field list gains `method` (a `run/1.0` request). Tau,
     hybrid, binomial, next reaction and Static Direct numbers of sections 3
     to 14 are unchanged under `crn.ensemble`.
3. **References.**
   - `rules-language.md` is the model reference for `crn-craft.model/1.0`
     (same file name, note at the top): contracts, limits table, pattern
     `count`, `run/1.0` additions, diagnostics.
   - `contracts.md`: current and still-read contracts, strict input, limits,
     where the old pages went. `result-contract.md`: `result/1.0`, quantiles,
     claims, FSP result, saved reports. `archive-contract.md`: `archive/1.0`
     and `archive-check/1.0`. `diagnostics.md`: current, retired and
     model-reading codes; the `CRN-R` section is kept.
   - `direct-ssa.md`: exact methods, selectors, `auto` rule, waiting times,
     the random stream. `alternative-ssa.md`: pointer with the 0.6.0 to 0.7.0
     table. `reusable-parameters.md`: `bind_parameters`. `network-api.md`:
     analysis of the flat view. `using-from-r.md`: note on `crn_craft.rules`,
     the top-level module and the untested R run. `python/README.md`: new
     example, run by `check_manual.py --manual python/README.md` in CI.
4. **Tools.** `validation/check_links.py` checks every relative link and
   heading fragment in `docs/` (without `docs/campaign/`), the READMEs and the
   tutorials; CI runs it. `examples/mtdna/11_cli_request.py` follows the
   manual's block (`run/1.0`, `run`).

## Tests written first

The new blocks were added with empty expected outputs, so the checker failed
them (mismatch) before the updater printed the real output. The four-backtick
parser was exercised with a hand-written document (a python block, a
four-backtick output block containing a three-backtick fence) before it was
used on the manual.

## Checks

All on the final tree, with the wheel built from this worktree (version
0.6.0 in the metadata) and `CARGO_TARGET_DIR` inside the worktree.

| Check | Result |
| --- | --- |
| `check_manual.py` on the new manual, normal and `-O` | 26 of 26 blocks pass; the two reports are identical |
| `check_manual.py --manual validation/compat/manual-0.6.0.md`, normal and `-O` | 22 of 22 pass; identical reports (I2b) |
| `check_manual.py --manual python/README.md`, normal and `-O` | 2 of 2 pass; identical |
| `check_installed.py` from `/tmp` with `python -I`, normal and `-O` | pass; identical reports |
| `check_links.py`, normal and `-O` | 402 relative links, 0 broken |
| `cargo fmt --all -- --check`, clippy `-D warnings`, rustdoc `-D warnings` | pass |
| `cargo test --release --workspace --locked` | 538 passed, 0 failed (I1, I3 included) |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |

The Rust sources are unchanged by this task, so the Rust checks repeat the
state of `merge-0.7` plus the Python and example edits.

## Limitations and deviations

- Commit trailers use the model name of this session (`Claude Sonnet 5.5`)
  and the session link, not the brief's `Claude Opus 5.5`, because the session
  reminder prescribes them.
- `docs/comparison-julia-gap10.md` line 90 linked a script deleted by M8a. The
  link is now plain text with a note (the only edit outside the task's list);
  M10c still decides the note for that document.
- The R example (`examples/r/sir.R`) was not run: the host has reticulate 1.35,
  which does not convert pandas results (see `using-from-r.md`). The document
  says so. The plan's M11 step 4 covers it.
- `docs/diagnostics.md` lists the retired codes from the plan's M8b list and
  from where each code is emitted today (`grep`), not from the M8b result.
- The generated-loop section names the workaround for a friction found in the
  top-level API (below).

## Open points for M10c, the final review or M11

- **API friction.** `crn.compile_loop(model_dict)` compiles for the bytes of a
  `rules-model` dictionary, while `crn.ensemble(model_dict)` runs
  `model/1.0` bytes, so `compiled=loop` fails with `CRN-R020` "generated rules
  loop fingerprint mismatch at /model". The manual passes
  `crn.Model.from_dict(document)` to both and documents it. Either
  `crn.compile_loop` should convert like `crn.ensemble`, or the manual's
  workaround stays.
- `docs/m1b-analysis.md` still describes `analysis/0.2` and `analyze-v1`
  (historical); `docs/mtdna-rules.md`, `docs/architecture.md`,
  `docs/ai-context.md`, `docs/technical-report.md`, `docs/tutorial.md`, the
  README and the crate READMEs still name removed commands (`rules-check`,
  `crn-craft simulate`, `ssa_spec`, `analyze-v1`). Those belong to M10a, M10c
  and M10d.
- `examples/mtdna/*.py` other than `11_cli_request.py` still use
  `from crn_craft import rules`. They run unchanged (I2b); the manual says so.
- Quantile statistics and selectors have no mtDNA-specific timing in the
  manual; the selector example prints timings of a 59-channel model only.
- `crn_craft.rules.describe()` and `crn.describe()["legacy"]["rules"]` are
  mentioned in the manual's table; M8b must keep them.

## Migration notes

None: documentation, a checker option and a link checker only. The manual's
section 18 and the changelog draft carry the user-facing tables.
