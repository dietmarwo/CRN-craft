# Extending the project

1. State the supported mathematical and contract scope and the milestone. Read the [architecture](architecture.md) and the plan of the task, and write a model card before implementing a numerical algorithm.
2. Add typed input fields or variants with strict unknown-field handling, validation, codes and field paths. New fields go into `crn-craft.model/1.0` or `crn-craft.run/1.0` (or a later version); the `rules-*` versions refuse them with `CRN-R020` at `/schema`. Update a delivered schema version when its shape changes. Generate snapshots from `crn-craft describe KIND`; do not hand-edit JSON Schema.
3. Add an independently justified acceptance fixture, including failure cases. A parser change must keep the duplicate-key, numeric and depth checks. A rate change must cover unavailable reactants, dependency reads and the flat view (`FlatView::propensity` returns the engine's value), so that FSP, analysis, kinetics and export see the same rates.
4. Keep numerical compilation and algorithms in the `crn-craft` crate (`crates/crn-craft/src/rules/` for the engine); the CLI handles files and JSON output, and Python and the job adapters call the same core.
5. Keep the stored requests replaying: `crates/crn-craft/tests/stored_requests.rs` and `static_direct_defaults.rs` must pass unchanged, and `crates/crn-craft/tests/generated/` changes only with the generators (regenerate with the facade crate's examples `generate_rules_kernel` and `generate_rules_loop`).
6. Run the documented checks and retain evidence with toolchain identity. Never introduce a reference simulator as a core dependency.

For a new rate form, add exact small-count checks against hand-computed propensities and overflow and domain tests (`CRN-R030`). For a new exact method or Static Direct selector, give it a sampler identity, then add stream vectors, analytic marginal tests, report-time boundaries, absorbing states and budget and cancellation tests, and run the law gate (`validation/merge/check_laws.py`) and the exact-law checks (`validation/check_alternative_ssa.py`, `validation/check_cme.py`) before interpreting performance.

The dependency policy is intentionally closed: a new transitive package fails `check_dependencies.py` until reviewed and listed. Check package sources and `links` metadata before updating the allowlist. This gate is a review aid, not a comprehensive supply-chain audit.

## Checks

CI (`.github/workflows/ci.yml`) is the complete list. The core of it:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
cargo build --release --locked -p crn-craft-cli
python3 -m pip install -r validation/requirements.txt
python3 -m ruff format --check validation python
python3 -m ruff check validation python
python3 validation/check_contracts.py
python3 validation/check_ssa.py --output ssa.json
python3 -O validation/check_ssa.py --output ssa-optimized.json
diff ssa.json ssa-optimized.json
python3 validation/check_cme.py --output cme.json
python3 validation/check_rules.py --output rules.json
python3 validation/check_mtdna.py --binary target/release/crn-craft --output mtdna.json
python3 validation/check_links.py
python3 validation/check_rust_docs.py
# After installing the wheel (maturin build --release --locked --out dist; pip install dist/*.whl):
python3 validation/check_python.py --binary target/release/crn-craft --output python.json
python3 validation/check_rules_python.py --output rules-python.json
python3 validation/check_guide.py --binary target/release/crn-craft --require-cargo --output guide.json
python3 validation/check_manual.py --binary target/release/crn-craft --manual validation/compat/manual-0.6.0.md --output manual-0.6.0.json
(cd /tmp && python3 -I "$OLDPWD/validation/check_installed.py" --output installed.json)
python3 validation/mutation_check.py
```

Each Python check also runs under `python -O`, and CI compares the two
reports.

`check_guide.py` runs every page of `docs/guide/`, the two READMEs and
`docs/index.md`; its docstring gives the rules for blocks. A new page needs
an entry in `validation/guide-blocks.json`, the least number of runnable
blocks it has. After a change, `--write` replaces each shown output by the
real one; never type an output by hand. The fitting page runs in smoke mode
in CI; `--long` runs it in full (over two hours) before a release. The Rust
blocks of `developer/rust-api.md` and the crate READMEs are checked by
`check_rust_docs.py`; a block that is not a complete program is marked
`<!-- check: fragment -->`. The engine-merge gates (`validation/merge/check_laws.py`, the speed
gate and the worker comparison) are run by hand and recorded in
[validation/merge/RESULTS.md](../../validation/merge/RESULTS.md). GitHub
Actions runs fmt, clippy and release tests on Rust 1.85.0, 1.88.0 and
stable, and the installed-wheel acceptance on Linux, Windows and macOS. The
default build has no C/C++ simulator, Julia or Python runtime dependency; it
pins pure-Rust libm for rate evaluation.

Optional integrations have their own lockfiles and minimum Rust versions.
Release wheels pass the installed-package acceptance on Linux, Windows and
macOS ([RELEASING.md](../../RELEASING.md)). A Julia reference run and broad
biological validation remain open work for the owner's platforms. The
circuit demonstrations are synthetic, and detector diagnostics do not
establish universal sensitivity. The ecosystem boundary: CRN-craft supplies
numerical models, fcmaes searches parameters, optimization-jobs executes
jobs, and experiment-craft coordinates campaigns.
