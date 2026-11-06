# Input and analysis contracts in M1b

Canonical definitions are `crates/crn-craft/src/contracts.rs`, `analysis.rs` and their referenced types, plus the generated files in `schemas/` (20 at present). `describe` includes all schemas and capability status. Syntax/shape validation precedes semantic checks. Generated schemas express selected local bounds; Rust adds constraints that schemas do not fully express.

## Model

`schema` is exactly `crn-craft.model/0.1`; `convention` is `stochastic_counts`. Model/species/parameter/reaction names are 1–64 ASCII identifier characters. Species and parameters share a disjoint namespace; reaction names are unique within reactions. Counts are non-negative i64, including integers larger than 2^53 without round-tripping through f64. Stoichiometry entries are integers from 1 through 2147483647. Omit zeros; empty sides represent sources/sinks. Zero reactions are allowed for an absorbing model.

Parameters must be finite. Their domain depends on where they are used: coefficients, vmax and basal are non-negative; k, km and n are positive. Ω is finite and positive. The rate-card formulas define meanings. `check` validates arguments; `analyze` also evaluates rates.

`population_limit` is optional, non-negative and at least initial_count. It is enforced by state validation and checked reaction updates; it is not a reflecting wall or automatic FSP truncation. A structural growth warning is not a boundedness theorem.

## Run

`schema` is exactly `crn-craft.run/0.1`. `model.path` is relative to the run's parent directory and must stay inside it after resolving symlinks. No empty components, absolute paths, `.`, `..`, drive prefixes, backslashes or NUL bytes. The supplied model must resolve to that same path. Input files are read once into bounded byte buffers and hashes refer to these buffers. The directory check is not a defense against a concurrently hostile filesystem that changes paths during checking.

`time.end` is finite and positive; `save_times` contains 2–1000000 finite, strictly increasing values, beginning at 0 and ending exactly at end. Time units are seconds. The future save rule is state after events at or before a save time.

SSA/direct requires an ensemble with positive replicate count and a positive per-replicate event cap. `root_seed` is u64. At least two replicates are needed if sample variance is requested. FSP requires inclusive truncation `bounds` for every species, containing the initial state, positive `max_states` and a finite tolerance between 0 and 1. FSP does not accept an ensemble, event cap or retained trajectories. Its state cap constrains reachable enumeration; it need not cover the Cartesian product of the box.

Both methods remain **contract-only** in M1b. FSP acceptance of a bounded truncation does not implement uniformization or establish any error bound.

Output statistics are `mean` and/or `variance`, without duplicates. Kept trajectories cannot exceed SSA replicates. Some output must be requested. Threads range from 1 to 1024; wall-time limits, if present, are finite and positive. The output byte budget is positive and must cover at least `8 * save_times * (1 + species * (statistics + kept_trajectories))`. This is a minimum numeric payload, not a serialized archive-size guarantee; live enforcement belongs to M2. FSP distribution output gets a new contract when implemented.

## Parsing, limits and API

- Input file size: at most 16 MiB, including whitespace.
- JSON nesting depth: at most 48, root at depth 0; no reference-based expressions.
- Species: 1–4096; parameters: at most 16384; reactions: at most 65536.
- Duplicate keys are rejected even when JSON escaping gives two spellings of the same key.
- Unknown fields are rejected in nested objects and tagged variants.
- `NaN`/`Infinity` are invalid JSON; overflowing numeric literals such as `1e400` fail with a numeric diagnostic.
- Decimal literals parse to the correctly rounded double, as Rust's `str::parse` and Python's `json` do, so the shortest round-trip text of a double reads back to the same bits (serde_json `float_roundtrip`).

`parse_model`/`parse_run` reject syntax and shape errors but are **not** semantic validation. Use `check_bytes` for byte/hash and semantic checks; use `check_files` when filesystem reference identity matters. `CompiledNetwork::compile` always calls `validate_model`, including for directly constructed Rust models and non-finite numbers that cannot appear in JSON. Run validation remains separate. `analyze_bytes`/`analyze_files` compile, evaluate and analyze without requiring a run document.

Rust float serialization may normalize literals, so hashes are never taken from re-serialized models to claim original-byte identity. Line endings are part of the bytes: the repository's `.gitattributes` keeps LF endings in every checkout, and CI runs the contract suite on a checkout made with `core.autocrlf=true`. BTreeMap ordering and ordered vectors make diagnostics deterministic on repeat checks.

## New state and analysis contracts

`crn-craft.state/0.1` contains exactly `schema` and `counts`, a map from every declared species name to an i64 count. Missing or extra names, negative counts and exceeded population limits are errors. Duplicate keys and unknown fields use the existing strict parser. No partial state or default count is inferred.

`crn-craft.analysis/0.1` contains `success`, `simulation_available: false`, raw `model_sha256`, optional raw `state_sha256`, `analysis` and diagnostics. Input/numerical failure has `success: false` and `analysis: null`. Successful analysis includes species in declaration order, initial and evaluation counts, reaction structures, dependencies, propensities and conservation. A conservation cap/overflow leaves analysis success true with an explicit unavailable conservation status and warning. No partial basis or floating-point rank guess is emitted. Read [network semantics](network-api.md) for matrix orientation, integer-string coefficients, limits and bound interpretation.

## Migration from iteration 1

> Since release 0.2.0 every generated schema, the legacy snapshots included,
> carries the Rust API documentation as `description` text. Shapes are
> unchanged, so the statements below about unchanged snapshots hold for shape;
> the bytes now differ by those descriptions.

The four model/run/result/check 0.1 schemas are byte-for-byte unchanged. No model or run migration is needed. Discovery advances from `crn-craft.describe/0.1` to `crn-craft.describe/0.2`: clients must accept the new version, `analyze`, `rate_evaluation_available`, analysis limits and the two new schemas. Simulation methods remain unavailable. Programs reading only `check` retain its existing structural/domain semantics and conservative warnings; `analyze` may suppress a growth warning when a conservation law proves an initial-class bound.

## Migration from iteration 2

Current `analyze` emits `crn-craft.analysis/0.2`. It keeps the original fields and adds `nonnegative_conservation`, `linkage` and `initial_class_upper_bounds` inside `analysis`. The last field selects complete cone bounds when the cone is exact, otherwise the previous sound basis-derived bounds. `conservation.initial_class_upper_bounds` retains its basis-only meaning. The new cone has independent status, limits, work count, extreme rays, totals and bounds. Successful structure/rate evaluation can coexist with an explicitly unavailable cone; inspect its status.

`analyze-v1 MODEL [STATE]` and `analysis_v1::analyze_bytes` preserve the original analysis 0.1 shape and warning policy. Its generated schema is `describe analysis-v1` / `schemas/analysis-v1.schema.json`, byte-identical to the old analysis snapshot. Model/run/check/result/state 0.1 snapshots are unchanged. No model, run or state migration is needed.

Discovery advances to `crn-craft.describe/0.3`, adding the legacy-analysis/export commands, exact-cone limits and the pinned export target. `crn-craft.export/0.1` reports raw model SHA-256, success, an optional Python artifact, target identity, aliases, limitations and diagnostics. An export failure has `artifact: null`; unsupported semantics are never silently removed. CLI output remains JSON and exits 0/1. See [reference export](reference-export.md).

## M2 execution and result migration

`simulate MODEL RUN NEW_DIRECTORY [--cancel-file PATH]` executes Direct SSA;
`verify DIRECTORY` checks its archive. The original model/run shapes remain
0.1. Result 0.2 replaces the never-emitted reserved result for actual execution;
`describe result-v1` retains its old schema. Archive and archive-check begin at
0.1; discovery is 0.10, or 0.12 from the native CLI (`alternative_exact_ssa`
lists the `ssa-run`/`ssa-result` methods, the four accepted request versions
and the `fast` default of each: none in 0.1, off in 0.2, on in 0.3 and 0.4,
where `next_reaction` runs the checked policy only; the top-level
`aggregators` list still describes legacy Direct only). FSP executes via
`fsp`. Existing check/analysis/export
reports remain unchanged; query discovery for current executable methods.

New runtime limits supplement input-domain validation. Invalid/unsupported
execution creates no archive. Incomplete trajectories cannot acquire an exact
ensemble assessment or unconditional statistics. See [Direct SSA](direct-ssa.md)
and the [result/archive contract](archive-contract.md) for numerical mappings,
state/sample semantics, output accounting, completion markers and limitations.
