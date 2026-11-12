# M1 log: one crate, no behavior change

Task M1 of [the plan](../../plan-0.7.0-engine-merge.md): move the rules
engine from the crate `crn-craft-rules` into the module `crn_craft::rules`
and turn `crn-craft-rules` into a facade, without changing behavior.

## What changed

- `git mv crates/crn-craft-rules/src/*.rs crates/crn-craft/src/rules/`,
  with `lib.rs` renamed to `mod.rs` (21 files, 14,205 lines). `pub mod
  rules;` is added to `crates/crn-craft/src/lib.rs`.
- Path edits inside the moved files, applied mechanically:
  - `crate::` became `crate::rules::` (code, doc links and the string
    paths in `deserialize_with = "..."` and `schemars(with = "...")`);
    `crate::rules::` is used everywhere rather than `super::`, because the
    `#[path]` submodules of `engine.rs` sit one level deeper;
  - `crn_craft::` became `crate::` in code and unit tests; doctests keep
    `crn_craft::sha256`, which is the right path from a doctest;
  - doctest paths `crn_craft_rules::` became `crn_craft::rules::`. The
    emitted source of the kernel and loop generators (`generated.rs`,
    `generated_loop.rs`, `expr.rs`) and the doctest assertions on that
    source keep `crn_craft_rules::` (invariant I3);
  - `include_str!`/`include_bytes!` paths of the unit tests gained one
    `../` (`../../../../models/rules/...`);
  - `cargo fmt` then rewrapped the lines that the longer paths pushed over
    the width (adding trailing commas where it split argument lists) and
    reordered `use` lines alphabetically.
- Crate-level items of the old `lib.rs`:
  - `#![warn(missing_docs)]` stays as an inner attribute of `rules/mod.rs`,
    where a lint attribute is valid; the `crn-craft` root also has it, and
    `#![forbid(unsafe_code)]` plus the workspace lint `unsafe_code =
    "forbid"` cover the module. The old crate had no other crate-level
    attributes.
  - The crate docs became the module docs of `crn_craft::rules`. Two
    intra-doc links named the crate `crn_craft`, which rustdoc cannot
    resolve inside that crate; they now link to `crate` and
    `crate::diagnostic::Diagnostic` with unchanged text.
- `crates/crn-craft-rules/src/lib.rs` is the facade: the 0.6.0 crate docs
  unchanged (so its README, docs.rs page and quickstart doctest stay), then
  `pub use crn_craft::rules::*;`. Its `Cargo.toml` (name, version,
  dependencies) is unchanged.
- New test `crates/crn-craft-rules/tests/facade.rs`: every public path of
  0.6.0 resolves through the facade and is the engine's own item (same
  `TypeId`, same function item, equal constants, facade-trait
  implementations accepted by the engine's generic entry points). The list
  is the rustdoc item index of 0.6.0 (`all.html`, 78 items) plus the root
  re-export `Model` and the hidden `libm`.
- `crates/crn-craft-rules/tests/` (stored requests, generated fixtures,
  defaults, independent variables) stays and now runs through the facade.

Outside `crates/` only this log is new: no schema, Python, binding, CLI
or lockfile edit.

## Public paths that needed attention in the facade

The single glob re-export suffices for every path; no path needed an
explicit re-export to compile. Checked individually:

| Path | Why it needed checking | Outcome |
| --- | --- | --- |
| `crn_craft_rules::{atp, continuous, model, observe}` | public modules | glob re-exports modules |
| `crn_craft_rules::Model` | root `pub use model::Model` | glob re-exports re-exports |
| `crn_craft_rules::{LoopHalt, LoopHost}` | re-exported from the private `engine::loop_host` | resolves; used by every generated loop |
| `crn_craft_rules::{RulesKernel, RulesLoop, generate_rules_*, execute_with_rules_*}` | re-exported from private modules | resolves; the two `execute_with_rules_*` take `&impl Trait`, so the test checks them by coercing to `fn` pointers with stub implementations of the facade traits |
| `crn_craft_rules::libm` | `#[doc(hidden)] pub use libm;`, called by generated rates and stop predicates as `crn_craft_rules::libm::{pow, exp, log, sqrt, floor, ceil}` | resolves through the glob (`#[doc(hidden)]` does not affect name resolution); see the rustdoc limitation below |
| `crn_craft_rules::observe::EventObserver` | has an associated constant, so it is not dyn-compatible | test implements it through the facade and passes the type to the engine's generic bound |
| `crn_craft_rules::observe::Event<'a>`, `LoopHost<'h, 'm>` | lifetime parameters | compared at `'static` |

The plan's extra line `pub use crn_craft::rules::{atp, continuous, model,
observe, libm};` was not added: those names come with the glob, and the test
proves it. As a check that the test discriminates, a facade with only that
line (no glob) fails the test with 54 unresolved paths.

## Tests written first

The move itself is guarded by tests that already existed: the 69 rules
unit tests (now in `crn-craft`'s lib tests), the 14 rules doctests (now
`crn-craft` doctests) and the facade crate's integration tests (I1, I3).
The new `facade.rs` was written once the facade compiled; its failing case
is the incomplete facade above.

## Checks

All on this host (2 CPUs, Rust 1.97.0, Python 3.12 venv), heavy steps under
the CPU lock, on the final tree.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 108 (39 network + 69 moved rules unit tests), its untouched integration tests, facade crate: `facade` 3, `generated_kernel` 3, `generated_loop` 25, `independent` 7, `static_direct_defaults` 9, `stored_requests` 2, CLI 12; doctests `crn_craft` 33 (19 + 14 moved), `crn_craft_rules` 1 (the facade's crate quickstart) |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (the existing cargo warning about the `crn-craft` binary and library sharing a doc file name is unchanged) |
| I1: stored requests replay byte-identically | pass (`stored_requests`, `static_direct_defaults`) |
| I3: generated fixtures match | pass (`generated_kernel`, `generated_loop`) |
| `crn-craft describe`, `describe KIND` for all 20 kinds, `--version` | byte-identical to the output of the 0.6.0 binary built before the move |
| Rustdoc item index of `crn_craft::rules` | the same pages as the 0.6.0 crate `crn_craft_rules` |
| Moved files against 0.6.0 | equal after the mechanical substitutions and `rustfmt` (checked by a script), except the two doc links above and the lint attribute on `Model::continuous` (see Limitations): 11 changed lines |
| `maturin build --release --locked --out dist`, wheel installed in the venv | pass |
| I2: `check_manual.py --binary target/release/crn-craft --require-cargo`, `CRN_CRAFT_SOURCE` set | pass, 22 blocks; `-O` report identical |
| `check_compiled_loop.py --require-cargo`, `CRN_CRAFT_SOURCE` and `CARGO_TARGET_DIR` as in CI | pass, 50 rows (31 identical to the interpreter, 18 refusals, 0 differences); `-O` report identical |
| `check_contracts.py` | pass, 86 checks; `-O` identical |
| `check_rules.py` | pass, 47 checks; `-O` identical |
| `check_installed.py` (`python -I`, outside the tree) | pass, 16 checks; `-O` identical |
| `check_rules_python.py`, `check_builder.py`, `check_hybrid.py`, `check_mtdna.py` | pass; `-O` reports identical |
| `check_dependencies.py` | pass (36 packages, policy unchanged); `-O` identical |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| `scripts/check_release.py` | pass |
| `cargo publish -p crn-craft -p crn-craft-rules -p crn-craft-cli --dry-run --locked --allow-dirty` | pass; each crate verified from its archive (`crn-craft` 56 files including the 21 in `src/rules/`, `crn-craft-rules` 8 files); the index access worked, with the expected warnings that 0.6.0 is already published |
| Optional workspaces `integrations/jobs`, `integrations/diffsol`, `validation/parameter-optimization/worker`, and `bindings/python` | `cargo metadata --locked` passes in all four; `cargo build --locked` passes in the first three (temporary target directory, removed afterward); no lockfile changed |

No timings were measured; M1 has no declared targets.

## Limitations

- rustdoc renders the facade a little differently from 0.6.0. When it
  inlines a glob re-export from another crate it ignores the
  `#[doc(hidden)]` on the re-exported `libm`, so the facade's docs now list a
  `libm` module; an explicit `#[doc(hidden)] pub use crn_craft::rules::libm;`
  in the facade, before or after the glob, did not hide it, so it was not
  kept. `Model` also gets its own page at the facade root, next to
  `model::Model`. Paths and behavior are unaffected.
- The doc of `Model::continuous` keeps the text `See [`crate::continuous`].`
  because schemars copies it into the `rules-model` schema description;
  rewriting it to `crate::rules::continuous` changed `crn-craft describe
  rules-model`. rustdoc cannot resolve the old path inside `crn_craft`, so
  the field carries `#[allow(rustdoc::broken_intra_doc_links)]` and the link
  renders as plain text. A later schema version (M3) can fix the text.
- The facade keeps the 0.6.0 dependency list, although its library uses
  only `crn-craft` (and its doctest and tests `serde_json`). Trimming it
  would change the `crn-craft-rules` entries of `Cargo.lock` and
  `bindings/python/Cargo.lock`; it is left for M8c or M11, when the facade
  shrinks to a smoke test.
- `pub(crate)` items of the rules engine are now visible to all of
  `crn-craft`. This widens nothing outside the crate.
- `python/crn_craft/_rules_loop.py` hashes `crates/crn-craft/src` into the
  compiled-loop cache key, so the move invalidates compiled-loop caches built
  from a 0.6.0 checkout once (a rebuild, no behavior change).
- The repository docs that describe the crate layout (`docs/architecture.md`,
  `docs/specialized-kernels.md`, the crate READMEs and the `crn-craft` crate
  docs, which still call `crn-craft-rules` a separate crate) are not changed
  here; M10 rewrites them.
- MSRV toolchains (1.85, 1.88) of the CI matrix are not installed on this
  host; fmt, clippy and tests ran with 1.97.

## Deviations from the plan

- The facade is the glob alone, without the explicit
  `{atp, continuous, model, observe, libm}` line (see above).
- New test file `crates/crn-craft-rules/tests/facade.rs` (not in the plan).
- Gate wording "`git diff --stat` shows renames and path edits only": plain
  `git diff -M --stat` shows `crates/crn-craft/src/rules/mod.rs` as a new
  file and the facade `lib.rs` as rewritten, because the facade reuses the
  path `crates/crn-craft-rules/src/lib.rs`. `git diff -M -B --stat` pairs
  `lib.rs → rules/mod.rs` as a rename with 64 changed lines.

## Migration notes

- Rust: `crn_craft_rules::X` is now `crn_craft::rules::X`. The old paths
  keep compiling through the facade crate, which re-exports the module.
- Generated kernels and loops still name `crn_craft_rules::` (I3), and the
  Python loop driver still depends on `crn-craft-rules`; M8c switches both.
- `crn_craft_rules::describe()` still reports `"crate": "crn-craft-rules"`
  (M7 changes it).
