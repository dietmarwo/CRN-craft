# Changelog

## [Unreleased]

### Added

- M0 Rust workspace with strict model/run contracts and reserved result envelope.
- JSON duplicate-key, unknown-field, domain, identity, budget and truncation checks.
- Structured `describe` and `check` CLI commands with explicit unavailable simulation status.
- Generated schemas, complete examples, negative fixtures and regression/acceptance suites.
- Local validation evidence, CI configuration, revised plan and agent/release documentation.
- Tests for unsupported run versions, methods and aggregators, and for both sides of the JSON nesting limit.
- Tests for both sides of the population-limit, FSP-bound and minimum-output-budget rules, which three probe mutants had passed.
- A CI job that checks out the repository with `core.autocrlf=true` and runs the contract suite on it.

### Changed

- CI cancels superseded runs and runs the dependency gate once, normally and under Python `-O`, comparing the two reports.
- `SOURCE_MANIFEST.json` from the delivery archive is not tracked; git history identifies the files.

### Fixed

- Decimal literals parse to the correctly rounded double (serde_json `float_roundtrip`). The default parser read about a quarter of random literals one ulp away, including shortest round-trip forms such as `3.983236103458967e-14`.
- Text files keep LF line endings in every checkout (`.gitattributes`), so exact-byte model hashes hold under Git for Windows' `core.autocrlf=true`.
- The nesting-limit hint states the rule the parser applies: depth at most 48, with the root at depth 0.

### Limitations

- No propensity evaluator, network analyzer, SSA, FSP, ODE, bindings, importers or ecosystem worker yet.
- No simulator correctness, throughput or biological validation is claimed.
