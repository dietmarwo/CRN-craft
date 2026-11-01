# CRN-craft iteration 1: implementation log

**Date:** 1 November 2026  
**Version:** 0.1.0-dev.1  
**Milestone:** M0 — workspace, contracts, discovery and preflight validation  
**Delivery:** full new project source, revised plan and this log

## Outcome

Implemented and tested M0. The project builds as a Rust library and a `crn-craft` CLI. It can describe its actual capabilities and validate model/run JSON, including byte identity and filesystem references. It has no numerical simulator yet. The next implementation is M1a; Direct SSA starts in M2.

The attached 31 October plan is preserved unchanged. A new reviewed plan is in `docs/implementation-plan-2026-11-01.md`, also delivered separately. The project ZIP includes source, Cargo.lock, generated schemas, examples, invalid fixtures, documentation, CI and local evidence; build products and cached dependencies are excluded.

## Review changes

- Defined the stochastic count convention and deterministic Ω/combinatorial mapping; nonlinear ODE trajectories are not claimed as exact stochastic means.
- Separated population safety stops from finite-state truncation. An unbounded CRN can still be approximated with a finite FSP domain and leakage accounting.
- Corrected reproducibility scope to completed fixed-work numerical payloads, excluding time/provenance and deadline-dependent completion.
- Replaced impossible cyclic-tree testing with malformed-reference, duplicate-key and depth-limit tests.
- Corrected conservation-basis comparisons, correlated statistical testing, uncertainty/censoring claims and FSP error bookkeeping.
- Made later aggregators and tau-leaping depend on measured workloads; the first useful circuit migration can use verified Direct without waiting for every aggregator.
- Replaced unsupported exclusivity/performance statements and calendar promises with explicit milestone evidence gates.

## Implemented files and behavior

| Area | Main paths | Behavior |
|---|---|---|
| Workspace | Cargo.toml, Cargo.lock, two crate manifests | Edition 2024, MSRV 1.85, project unsafe code forbidden, publishing disabled |
| Contracts | crates/crn-craft/src/contracts.rs | Model/run/reserved-result types; tagged rates/methods; strict fields |
| JSON parser | crates/crn-craft/src/strict_json.rs | Duplicate-key rejection before map conversion; 16 MiB/depth limits; exact integer parsing; coded errors |
| Semantic checks | crates/crn-craft/src/validate.rs | Names, references, domains, counts, stoichiometry, time/ensemble/budget/output rules, FSP truncations |
| File/API boundary | crates/crn-craft/src/lib.rs | Raw SHA-256; bounded regular-file reads; relative model reference identity; structured reports |
| CLI | crates/crn-craft-cli/src/main.rs | describe/check/help/version; JSON output and defined exit codes; unimplemented commands fail |
| Examples | models/examples/ | Four complete models; four Direct-run contracts and one FSP-run contract |
| Negative fixtures | models/invalid/expected.json and JSON files | 13 named invalid cases with expected diagnostic codes and pointers |
| Verification | crate integration tests; validation/*.py | Rust regression tests, independent schema/CLI acceptance, dependency audit |
| Documentation and CI | docs/, README, AGENTS, RELEASING, .github/workflows/ci.yml | Agent context, formulas, scope, statistical policy, extension cards and local/hosted validation distinction |

The repressilator now includes all three transcription/translation/degradation modules: six species and twelve reactions. The enzyme-mechanism model and reduced Michaelis–Menten model remain separate examples. They are syntax/domain fixtures, not numerically validated biological models.

Every check result states `execution_available: false`. `describe` labels Direct/FSP as contract-only. Generated result types are reserved; no fake run output is emitted. Model/run input format changes after this delivery require versioning.

## Tests actually executed

| Validation | Observed result |
|---|---|
| Release tests on Rust 1.85.0 | 36 passed, 0 failed |
| Release tests on Rust 1.88.0 | 36 passed, 0 failed |
| Release tests on Rust 1.99.1 (installed stable) | 36 passed, 0 failed |
| fmt and clippy `-D warnings` | Passed on all three toolchains |
| Independent schema/CLI acceptance | 43 passed, 0 failed |
| Same acceptance under Python `-O` | 43 passed; byte-identical report |
| Python formatting/lint | Passed |
| Dependency audit, normally and under `-O` | Passed; identical reports; five injected bad dependency cases rejected |

The 36 Rust tests consist of 30 library/contract and 6 CLI tests. Repeating them on three toolchains is compatibility evidence, not three different test suites. The 43 Python checks are similarly repeated under `-O`. The full reports are in `validation/evidence`; `validation/RESULTS.md` explains their coverage.

Testing used Linux x86_64, Python 3.12.14, jsonschema 4.26.0 and Ruff 0.11.13. Rust versions were read from the actual tools. The core was built from the shipped lockfile with cached crates; external dependencies are not vendored in the archive. A new machine needs ordinary Cargo dependency access once.

## Archive verification

The staged project ZIP passed its CRC check and every extracted file matched `SOURCE_MANIFEST.json`. A fresh extraction built successfully with Rust 1.99.1 using `--release --workspace --locked --offline`, and its independent schema/CLI suite passed all 43 checks. This evidence is retained in `validation/evidence/archive-check.json` and `archive-build.log`. The final archive adds those records and this note; runtime source and schema snapshots are unchanged. Its manifest and extraction checks are repeated before delivery.

## Issues found and resolved

The first build exposed a Schemars API mismatch for an exclusive minimum attribute. A schema helper now generates the positive-domain constraint, and independent schema checks verify it. One negative fixture expected a less precise path than the parser produced; it now pins `/reactions/0/rate/kind`. Clippy caught one test-formatting warning. Dependency-gate review added explicit rejection of unreviewed local packages and non-registry sources. All final gates passed after these corrections; no statistical threshold was relaxed.

The remote repository URL was not successfully inspected: the public lookup was unavailable and unauthenticated `git ls-remote` could not resolve it without prompting for credentials. That does not establish whether a private repository exists. The attached plan identified the repository as not yet created; this task was fulfilled as a new local project and ZIP, without claiming a GitHub merge, push or hosted-CI run.

## How to reproduce

After extracting the ZIP, enter `CRN-craft`:

```bash
cargo build --release --workspace --locked
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
python3 -m pip install -r validation/requirements.txt
python3 -m ruff format --check validation
python3 -m ruff check validation
python3 validation/check_dependencies.py
python3 validation/check_contracts.py
python3 -O validation/check_contracts.py
./target/release/crn-craft check models/examples/repressilator.model.json models/examples/repressilator.run.json
```

For explicit compatibility checks, use `cargo +1.85.0`, `cargo +1.88.0` or the installed `+stable` after installing their fmt/clippy components. The checked-in CI config does this on GitHub when the project is uploaded; it has not run there as part of this delivery.

## Remaining work and claims

**Next: M1a.** Compile indexed stoichiometry and dependencies; implement checked mass-action/Hill/MM propensities and reactant gates; expose `analyze`; verify exact small-count values and 200 seeded small networks independently. Complete conservation/linkage work in M1a/M1b, with checked overflow and rational-span tests.

**Then: M2.** Implement Direct SSA, pinned samplers, ensembles, statistics, budgets, cancellation and verifiable run archives. M3 provides FSP/DSMTS/independent numerical checks. M4a establishes a reference speed table; M5 builds the minimal Python interface and circuit study.

No solver speedup, performance ratio, exact-SSA correctness, FSP guarantee, cross-platform stochastic identity, optimization integration or biological conclusion is claimed for this iteration. M0 is the validated foundation needed to make those later claims testable.
