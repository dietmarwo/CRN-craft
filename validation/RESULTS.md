# Iteration 1 validation results

**Date:** 1 November 2026. **Scope:** contracts, validation and CLI only. No simulator benchmark or statistical solver campaign was performed.

| Check | Result | Evidence |
|---|---|---|
| Rust 1.85.0 MSRV release suite | 36 tests passed, 0 failed | evidence/rust-1.85-tests.log |
| Rust 1.88.0 release suite | 36 tests passed, 0 failed | evidence/rust-1.88-tests.log |
| Installed stable, Rust 1.99.1 release suite | 36 tests passed, 0 failed | evidence/rust-stable-tests.log |
| fmt and clippy with warnings denied | Passed on all three toolchains | evidence/rust-*-fmt.log and rust-*-clippy.log |
| Independent Python JSON Schema + CLI suite | 43 checks passed, 0 failed | evidence/contracts.json |
| Same suite under Python -O | 43 checks passed; report byte-identical | evidence/contracts-optimized.json |
| Python formatting/lint | Passed | evidence/python-format.log and python-lint.log |
| Dependency gate | Passed; 5 deliberately bad dependency probes rejected | evidence/dependencies.json |
| Dependency gate under Python -O | Passed; report byte-identical | evidence/dependencies-optimized.json |

Rust totals comprise 30 contract/library integration tests and 6 CLI integration tests, on Linux x86_64. They cover the positive examples and all 13 persisted negative fixtures, plus malformed/deep/oversized JSON, duplicate escaped keys, values above 2^53, count overflow, unknown nested fields, numeric domains, invalid budgets, time grids, population limits, hash/path identity and escaping symlinks. The three runs are toolchain compatibility checks of the same 36 tests, not 108 distinct tests.

Python independently parses generated schemas with jsonschema, validates examples and check reports, compares schema snapshots, verifies hashes with hashlib, checks malformed-input behavior and confirms that unimplemented commands fail. Its 43 checks are repeated under `-O`; assertions are not used for validation.

The dependency gate reviews 27 external package versions (26 distinct names, including two syn versions) plus two workspace crates. It rejects unreviewed names, local packages, non-registry sources, known native build dependencies and Cargo `links` metadata. No Julia, Python runtime, C/C++ solver, RNG or ODE package is in the core tree. This is a dependency-policy check, not a proof that dependencies use no unsafe code or platform instructions.

The staged source archive also passed CRC and manifest verification, built from a fresh extraction with Rust 1.99.1, and passed all 43 independent contract checks from that extracted project. See `evidence/archive-check.json` and `evidence/archive-build.log`. Final packaging adds these evidence records without changing runtime source or schemas.

## Development findings retained

1. The first compile used `schemars` range syntax unsupported by pinned 0.8.22. Replaced it with a type-derived positive-number schema helper; generated-schema tests verify zero Ω is rejected.
2. The unsupported-expression fixture initially expected `/reactions/0/rate`; parsing correctly reported `/reactions/0/rate/kind`. The fixture was corrected to require the more precise field pointer. No validation rule was weakened.
3. Clippy found an uninlined format argument in a CLI test; corrected it, then reran the gate.
4. Explicit toolchain inspection distinguished the default Rust 1.88.0 from the installed stable alias, Rust 1.99.1. Both and MSRV 1.85.0 were tested; documentation records actual versions.
5. Reviewing the dependency gate found that a same-name git or local dependency needed separate rejection. Added source/local allowlists and corresponding mutation probes.

## Repository review, 1 November 2026

The delivery was imported unchanged as the first commit of `dietmarwo/CRN-craft`, including `SOURCE_MANIFEST.json`, whose 86 file hashes all matched. The changes after that commit are listed in CHANGELOG.md. Results on the reviewed tree, on Linux x86_64 with Rust 1.97.0 (the stable toolchain available there), Python 3.13.18 and the pinned jsonschema and Ruff:

| Check | Result |
|---|---|
| fmt and clippy with warnings denied | Passed |
| Release tests | 41 passed, 0 failed: 35 contract/library and 6 CLI tests |
| Independent Python JSON Schema + CLI suite | 43 checks passed; normal and `-O` reports byte-identical, and identical to `evidence/contracts.json` |
| Dependency gate | Passed; normal and `-O` reports identical, and identical to `evidence/dependencies.json` |
| Checkout with `core.autocrlf=true` | Without `.gitattributes`, 82 files were converted to CRLF and the example-run test failed with hash mismatches. With it, no file was converted and the contract suite passed |
| Decimal float parsing | Of 10000 shortest round-trip literals of random doubles, 2983 parsed to a different double before `float_roundtrip` was enabled, 0 after |
| Mutation probe | 12 mutants of parser and validation rules: the delivered tests caught 9; the three survivors are caught by the added tests |

Rust 1.85.0 and 1.88.0 could not be installed in the review environment; the GitHub Actions matrix runs both. The `evidence/` directory holds the delivery's records and was not regenerated.

## Limits

GitHub Actions runs on Linux only; Windows checkouts are covered by the `core.autocrlf=true` job, not by a Windows runner. Other operating systems/architectures, SSA accuracy, FSP error estimates, stochastic reproducibility, throughput, Python bindings and ecosystem workers are untested because their implementations are outside M0. No rate formula has been evaluated by CRN-craft in this iteration.

The environment is recorded in `evidence/environment.json`; input-document hashes in `evidence/review-inputs.json`. These are local test results and do not characterize the user's home hardware.
