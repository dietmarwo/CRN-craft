# Changelog

## [Unreleased] — iteration 3 / M1b (`0.1.0-dev.3`)

### Added

- Exact nonnegative conservation-cone rays with deterministic work/ray caps, checked i128 arithmetic and explicit unavailable statuses.
- Effective initial-class bounds that resolve the reviewed mixed-sign-basis gap when cone enumeration completes.
- Reaction complexes, linkage/strong/terminal classes, weak reversibility and structural deficiency, with iterative graph walks.
- GillesPy2 1.8.3 elementary mass-action Python export, explicit homodimer conversion, collision-free aliases and unsupported-feature diagnostics.
- Independent support-enumeration/reachability oracle for 300 seeded networks, real target propensity/SSA-smoke checks and seven additional source mutants.

### Changed

- Analysis advances to 0.2; `analyze-v1` preserves 0.1. Discovery advances to 0.3; export starts at 0.1. Model/run/check/result/state schemas remain unchanged.
- W001 suppression uses exact cone bounds when available and retains sound basis bounds on unavailable outcomes.
- CI covers structure/oracle checks under normal Python and `-O`, 14 mutants after a passing baseline, and optional pinned reference validation.

### Limits

- Cone enumeration can grow exponentially; size alone does not ensure completion. No partial cone is labelled complete.
- Export supports mass-action reactant order at most two, no population safety limits and exactly representable initial counts up to 2^53. It is for reference SSA, not target ODE semantics or CRN-craft budget/range equivalence.
- Direct SSA and all CRN-craft simulation remain M2 or later.

## Iteration 2 / M1a (`0.1.0-dev.2`)

### Added

- Immutable compiled networks, sparse stoichiometry, reactant gates and propensity dependency graphs.
- Checked mass-action, Hill activation/repression and Michaelis–Menten evaluation using pinned pure-Rust libm 0.2.15; atomic i64 reaction updates.
- `analyze MODEL [STATE]`, new state/analysis 0.1 schemas, exact checked small-network conservation with explicit cap/overflow status and initial-class bounds.
- Independent 200-network Fraction/Decimal oracle, source mutation tests, numerical edge/atomicity regressions, documentation and retained evidence.
- Tests for the `log1p` Hill path near x = k, for weight division in conservation bounds and for binomial symmetry at high reactant orders, which three probe mutants had passed.

### Changed

- Discovery is `crn-craft.describe/0.2` and advertises implemented rate evaluation; model/run/check/result 0.1 contracts are unchanged.
- Analyzer growth warnings are suppressed only by proven initial-class bounds. Structural `check` semantics are unchanged.
- CI includes the independent network oracle normally and under `-O`, source mutation checks and network tests in the LF checkout job.
- The mutation check runs the unmutated suite first and requires it to pass, then runs seven mutants.

### Limitations

- No scheduler, RNG, SSA/FSP/ODE, linkage classes, expression AST, exporters or ecosystem integrations yet.
- Conservation uses checked i128 within 64 species / 256 reactions; its nonnegative basis heuristic is not a complete bound solver.
- f64 propensity evaluation is tested against declared tolerances, not certified correctly rounded or cross-platform identical.

## Iteration 1 / M0 (`0.1.0-dev.1`)

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
