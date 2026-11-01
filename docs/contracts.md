# Input contracts in M0

Canonical definitions are `crates/crn-craft/src/contracts.rs` and the four generated files in `schemas/`. `describe` includes all schemas plus capability status. Syntax/shape validation precedes semantic checks. Generated schemas express selected local bounds; `check` adds constraints that schemas do not fully express.

## Model

`schema` is exactly `crn-craft.model/0.1`; `convention` is `stochastic_counts`. Model/species/parameter/reaction names are 1–64 ASCII identifier characters. Species and parameters share a disjoint namespace; reaction names are unique within reactions. Counts are non-negative i64, including integers larger than 2^53 without round-tripping through f64. Stoichiometry entries are integers from 1 through 2147483647. Omit zeros; empty sides represent sources/sinks. Zero reactions are allowed for an absorbing model.

Parameters must be finite. Their domain depends on where they are used: coefficients, vmax and basal are non-negative; k, km and n are positive. Ω is finite and positive. The rate-card formulas define meanings; M0 checks their arguments without evaluating the formulas.

`population_limit` is optional, non-negative and at least initial_count. It denotes a future execution safety stop, not a reflecting wall or automatic FSP truncation. A structural growth warning is not a boundedness theorem.

## Run

`schema` is exactly `crn-craft.run/0.1`. `model.path` is relative to the run's parent directory and must stay inside it after resolving symlinks. No empty components, absolute paths, `.`, `..`, drive prefixes, backslashes or NUL bytes. The supplied model must resolve to that same path. Input files are read once into bounded byte buffers and hashes refer to these buffers. The directory check is not a defense against a concurrently hostile filesystem that changes paths during checking.

`time.end` is finite and positive; `save_times` contains 2–1000000 finite, strictly increasing values, beginning at 0 and ending exactly at end. Time units are seconds. The future save rule is state after events at or before a save time.

SSA/direct requires an ensemble with positive replicate count and a positive per-replicate event cap. `root_seed` is u64. At least two replicates are needed if sample variance is requested. FSP requires inclusive truncation `bounds` for every species, containing the initial state, positive `max_states` and a finite tolerance between 0 and 1. FSP does not accept an ensemble, event cap or retained trajectories. Its state cap constrains reachable enumeration; it need not cover the Cartesian product of the box.

Both methods are **contract-only** in M0. FSP acceptance of a bounded truncation does not implement uniformization or establish any error bound.

Output statistics are `mean` and/or `variance`, without duplicates. Kept trajectories cannot exceed SSA replicates. Some output must be requested. Threads range from 1 to 1024; wall-time limits, if present, are finite and positive. The output byte budget is positive and must cover at least `8 * save_times * (1 + species * (statistics + kept_trajectories))`. This is a minimum numeric payload, not a serialized archive-size guarantee; live enforcement belongs to M2. FSP distribution output gets a new contract when implemented.

## Parsing, limits and API

- Input file size: at most 16 MiB, including whitespace.
- JSON nesting depth: at most 48, root at depth 0; no reference-based expressions.
- Species: 1–4096; parameters: at most 16384; reactions: at most 65536.
- Duplicate keys are rejected even when JSON escaping gives two spellings of the same key.
- Unknown fields are rejected in nested objects and tagged variants.
- `NaN`/`Infinity` are invalid JSON; overflowing numeric literals such as `1e400` fail with a numeric diagnostic.
- Decimal literals parse to the correctly rounded double, as Rust's `str::parse` and Python's `json` do, so the shortest round-trip text of a double reads back to the same bits (serde_json `float_roundtrip`).

`parse_model`/`parse_run` reject syntax and shape errors but are **not** semantic validation. Use `check_bytes` for byte/hash and semantic checks; use `check_files` when filesystem reference identity matters. Directly constructed Rust values still require `validate_model` and `validate_run`; those functions check non-finite values that cannot appear in valid JSON. No validated/compiled network object exists until M1.

Rust float serialization may normalize literals, so hashes are never taken from re-serialized models to claim original-byte identity. BTreeMap ordering and ordered vectors make diagnostics deterministic on repeat checks.
