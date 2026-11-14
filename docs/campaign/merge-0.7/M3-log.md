# M3 log: unified contract names

Task M3 of [the plan](../../history/plan-0.7.0-engine-merge.md): add
`crn-craft.model/1.0`, `crn-craft.run/1.0`, `crn-craft.result/1.0` and
`crn-craft.check/1.0` as new versions of the rules contracts, without new
features, and dispatch by schema string.

## What changed

- Versions (`rules/model.rs`, `rules/mod.rs`): `ModelVersion::V10`
  (`crn-craft.model/1.0`), `RunVersion::V10` (`crn-craft.run/1.0`),
  `ResultVersion::V10` (`crn-craft.result/1.0`), `CheckVersion::V10`
  (`crn-craft.check/1.0`). Each is the last variant of its ordered enum,
  so every `<`/`>=` version test gives a 1.0 document the meaning of the
  newest 0.x version: a `model/1.0` model means what `rules-model/0.4`
  means, a `run/1.0` request what `rules-run/0.6` means (Ziggurat waits and
  automatic sparse trackers when absent, `method` and `continuous` echoed),
  and `run/1.0` returns `result/1.0` with the fields of `rules-result/0.3`.
  The variant docs use `concat!`, as the existing ones do.
- Dispatch:
  - `prepare` refuses a `model/1.0` model with a `rules-run` request:
    `CRN-R020` at `/schema` ("runs only with a crn-craft.run/1.0
    request"), so old requests cannot reach 1.0 features. `run/1.0` accepts
    every rules-model version and `model/1.0`.
  - `check_model` (`crn-craft rules-check`, `crn_craft.rules.load`)
    reports `crn-craft.check/1.0` for a `model/1.0` document, also when it
    does not parse (the top-level `schema` string is read leniently), and
    `rules-check/0.1` otherwise.
  - Network readers refuse rules-engine documents with the new
    `CRN-E032` at `/schema`: `parse_model`, `parse_run` and the `ssa-run`
    and `tau-request` parsers go through `parse_network`, which turns a
    strict-parse `CRN-E003` into `CRN-E032` when the document's `schema` is
    a `crn-craft.rules-*` or a `…/1.0` contract. Successful parses and other
    errors are unchanged. So `check`, `analyze`, `export`, `simulate`,
    `simulate-cached`, `fsp`, `tau`, `ssa` and the Python network functions
    refuse `model/1.0` models and `run/1.0` requests clearly instead of
    reporting the first unknown field. `docs/diagnostics.md` lists
    `CRN-E032` and the new `CRN-R020` case.
  - **`crn-craft check MODEL` with a `model/1.0` model refuses it**
    (`CRN-E032`, hint naming `rules-check`) instead of dispatching to the
    rules check. Reason: `check` keeps one report contract (`check/0.1`)
    until M7 defines the unified `check MODEL [REQUEST]` with
    `rules-check` as its alias; dispatching now would give one command two
    report shapes for one release step. Reversible in M7.
- The converter writes `crn-craft.model/1.0` (`CONVERTED_VERSION`).
- Schemas (`rules::schema`): new kinds `model-1.0`, `run-1.0`,
  `result-1.0`, `check-1.0`, generated from the same Rust types as the
  `rules-*` kinds. Each schema lists only its own versions in the version
  definition's `enum` (the 0.x versions for `rules-*`, the 1.0 version for
  the 1.0 kinds), so the four `rules-*` schemas stay byte-identical. The
  1.0 schemas get their own top-level and version descriptions, and the
  `continuous` description of `model-1.0` reads
  "See `crn_craft::rules::continuous`." instead of the 0.6.0 rustdoc path
  noted in the M1 log; the `rules-model` schema keeps the old text. Files
  `schemas/{model,run,result,check}-1.0.schema.json` were written by
  `crn-craft describe KIND`. Struct-level rustdoc about 1.0 is written with
  `concat!`, which the 0.x schemas do not see.
- Discovery: `crn_craft::rules::describe()` gains `contracts_1_0`
  (`model`, `request`, `result`, `check`, the 0.x versions they mean, the
  models a `run/1.0` request accepts, and the schema kinds); `convert.to`
  is now `crn-craft.model/1.0`. `accepted_models` and `accepted_requests`
  keep their 0.x lists: `static_direct_defaults.rs`
  (`discovery_lists_the_new_version`, part of the I1 check) asserts six
  request versions. The native discovery stays `describe/0.14` (set in M2,
  unreleased) and its help lists the new `describe` kinds.
- Python stays on 0.x: `rules.Model.to_dict()` and every request built for
  a 0.x model are unchanged (I2). Two places had to know the 1.0 gate:
  `rules.ensemble` writes `crn-craft.run/1.0` when the model document is
  `crn-craft.model/1.0` (only that version fits), and the generator request
  of `compile_loop`/`loop_source` does the same (the loop does not depend on
  the request). Without them a converted model could not run through
  `crn_craft.rules`.
- `validation/check_contracts.py` checks the four new schema snapshots.

## Tests written first

`crates/crn-craft/tests/contracts_v1.rs` (7 tests) was written before the
implementation and run against the M2 code: 7 of 7 failed. One failure was a
test bug (it expected at least 18 stored requests; there are 17), fixed by
comparing with the keys of `stored/expected.json` as `stored_requests.rs`
does. The tests:

- `every_stored_request_as_run_1_0_gives_the_rules_run_0_6_result`: the
  equivalence test. For each of the 17 stored requests (never edited), the
  request rewritten as `run/1.0` gives the result of the same request
  rewritten as `rules-run/0.6` except `schema` and `run_sha256`; the model
  rewritten as `model/1.0` (only its schema string changed) gives it except
  `schema`, `run_sha256` and `model_sha256`. The baseline is the 0.6
  rewrite, not the stored bytes, because `run/1.0` means 0.6: stored 0.2–0.5
  Static Direct requests without options run inverse waits under their own
  version and Ziggurat waits under 0.6 and 1.0, and older result versions
  lack `method` and `continuous`. The stored bytes themselves are I1.
- `run_1_0_accepts_every_rules_model_and_rules_run_refuses_model_1_0`
  (also: `bind_parameters` candidates keep the 1.0 versions).
- `checks_of_model_1_0_report_check_1_0`, `the_converter_writes_model_1_0`.
- `network_entry_points_refuse_1_0_and_rules_documents`: `check_bytes`,
  `parse_model`, `parse_run`, `analyze_bytes`, `export_gillespy2`,
  `simulation::simulate_bytes`, `ssa::prepare`, `tau::solve`, `fsp::solve`;
  other shape and syntax errors keep `CRN-E003` and `CRN-E001`.
- `schemas_match_their_snapshots_and_list_only_their_versions` (Rust
  snapshot test of the eight rules-family schemas) and
  `discovery_names_the_1_0_contracts_and_keeps_the_rules_lists`.

`crates/crn-craft-cli/tests/cli.rs` gained
`network_commands_refuse_1_0_documents_and_rules_commands_run_them`
(`check` with and without a run, `analyze`, `export`, `simulate` (no
directory created), `fsp`, `ssa`, `tau` → `CRN-E032`; `rules-check` →
`check/1.0`; `rules` with `run/1.0` → `result/1.0`; `describe` of the four
kinds equals the files). `tests/convert.rs` now expects `model/1.0`, runs
the converted model with `run/1.0`, and expects `CRN-E032` when given a
rules model or a converted model.

## Checks

On the final tree, under the CPU lock.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 111, `contracts_v1` 7 (new), `convert` 10, CLI 14 (1 new), facade 3, generated kernels 3 and loops 25 (I3), `independent` 7, `static_direct_defaults` 9 and `stored_requests` 2 (I1), all other network tests, doctests 34 + 1 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (only the known output-filename collision warning) |
| Every committed schema against `crn-craft describe KIND` | 24 of 24 byte-identical (20 old, 4 new) |
| `cargo build --release --workspace --locked`; wheel built from the root and installed | pass |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` (`CRN_CRAFT_SOURCE` set) | pass, 22 blocks; `-O` identical |
| `check_contracts.py` | 90 of 90 (86 + 4 new 1.0 schema snapshots); `-O` identical |
| `check_rules.py` | 47 of 47; `-O` identical |
| `check_compiled_loop.py --require-cargo` (`CRN_CRAFT_SOURCE`, `CARGO_TARGET_DIR` as in CI) | pass, 50 rows: 31 identical to the interpreter, 18 refusals, 1 cancellation, as in M1; `-O` identical |
| `check_rules_python.py` (22), `check_builder.py`, `check_python.py`, `check_alternative_contracts.py` (20), `check_importers.py` (48 rows), `check_network.py`, `check_structure.py`, `check_tutorial.py` (17 blocks) | pass; `-O` reports identical |
| `check_installed.py` (`python -I`, outside the tree) | 16 checks pass; `-O` identical |
| `check_dependencies.py` | pass; `-O` identical |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| Law gate carried over: `check_laws.py` with the M3 build on W1, enzyme-mechanism, immigration-death and michaelis-menten (200,000 replicates each; converted models now `model/1.0`, requests `run/1.0`) | power sums, rows and event counts equal to the M2 gate report for all four |

## Limitations

- `crn-craft check` refuses `model/1.0` instead of checking it (decision
  above); use `crn-craft rules-check` until M7.
- `rules-check` and `rules` still report `CRN-E003` for a network
  `model/0.1` document; M7 decides whether the unified commands convert it.
- The 1.0 schemas are the 0.x schemas with a filtered version list and
  replaced descriptions, done in `rules::schema`; a later change to a
  rules type changes both, as intended until the 1.0 contracts are
  delivered.
- Python's builder and request helpers keep writing 0.x documents for 0.x
  models (I2); top-level 1.0 output is M6.

## Deviations from the plan

- Besides the three schemas the plan names, `schemas/check-1.0.schema.json`
  is added, since `check/1.0` is an emitted contract.
- `accepted_models`/`accepted_requests` in the rules discovery keep their
  0.x lists (an I1 test pins them); the 1.0 versions are listed under the
  new key `contracts_1_0`.
- `CRN-E032` and the network-side refusal come from the coordinator's
  instructions, not the plan; `rules.ensemble` and the loop generator
  request select `run/1.0` for `model/1.0` models (needed for converted
  models to run from Python at all).

## Migration notes

- New contracts: `crn-craft.model/1.0`, `crn-craft.run/1.0`,
  `crn-craft.result/1.0`, `crn-craft.check/1.0`; schemas via
  `crn-craft describe model-1.0` (`run-1.0`, `result-1.0`, `check-1.0`).
  A 0.x rules document keeps its meaning and output byte for byte.
- `crn-craft convert` and `crn_craft.convert` now write
  `crn-craft.model/1.0`; run it with a `crn-craft.run/1.0` request
  (`crn_craft.rules.ensemble` selects it automatically).
- A `crn-craft.model/1.0` model with a `rules-run` request is refused
  (`CRN-R020` at `/schema`). Network commands refuse rules-engine documents
  with `CRN-E032` at `/schema`; before, they reported `CRN-E003` at the
  first unknown field.
