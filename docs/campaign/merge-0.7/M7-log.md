# M7 log: one command line

Task M7 of [the plan](../../history/plan-0.7.0-engine-merge.md) (section 3,
"Unified command line"), done with M6a on the branch `merge/M6`; step 0
(the check of the coordinator's merge) is in the [M6a log](M6a-log.md).

## What changed

### Commands (`crates/crn-craft-cli/src/main.rs`)

```text
crn-craft describe [KIND]
crn-craft check MODEL [REQUEST]
crn-craft run MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]
crn-craft verify DIRECTORY
crn-craft analyze MODEL [STATE]
crn-craft convert MODEL
crn-craft export MODEL gillespy2
crn-craft --version
```

- `rules` and `rules-check` are aliases of `run` and `check` (same code
  path, same bytes; `rules-check` now also takes a request).
- `run` is the former `rules` command: a simulation request prints the
  result as compact JSON with sorted keys, a `run/1.0` request with method
  `fsp` prints `fsp-result/1.0`, `--archive` writes `archive/1.0` and prints
  the `archive-check/1.0` report. A network model gets `CRN-E033`.
- `analyze MODEL [STATE]` keeps the optional state document of M5b (the
  plan lists `analyze MODEL`; `check_export.py` and `compare_julia.py` need
  rates at a given state, as the M5b log says).
- `simulate`, `simulate-cached`, `ssa`, `tau` and `fsp` (the network
  engine's commands) keep working until M8b but are hidden: absent from
  `--help`, from the usage message (`CRN-E031` hint) and from `describe`.
  In the code they are one function `legacy` and a constant `LEGACY`,
  both marked "LEGACY: delete in M8b".
- `--help` lists the visible commands, the aliases, which models and
  requests each reads, the `describe` kinds (1.0 first, then "still read")
  and the exit codes. `--version` prints `{"version": ..., "milestone":
  "unified-engine"}` (was `"SSA-acceleration"`).

### Discovery (`crates/crn-craft/src/unified.rs`)

`crn_craft::describe()` is the new `crn-craft.describe/1.0`, printed by
`crn-craft describe` and returned by Python's `crn_craft.describe()` (the
native CLI and the core no longer differ, as `describe/0.15` and `0.11`
did). Top-level fields describe the unified engine: `crate`
(`"crn-craft"`), `version`, `milestone`, `commands` (the seven above),
`aliases`, `python`, `contracts` (model, request, result, check, FSP
result, analysis, archive, archive check, export, state, import and report
versions), `schemas` (the `describe` kinds of the 1.0 contracts and current
reports), `methods` (with `fsp`), `approximate_methods`, `default_method`,
`selectors`, `quantiles` (with the Python default), `fsp`, `archives`,
`reference_exporters`, `limits` (the `model/1.0` limits and the run
limits), `sampler`, `reference`. Under `legacy`: the still-read model,
request, result, check and archive versions (rules-model 0.1–0.4, the
network `model/0.1`, rules-run 0.1–0.6, rules-result, rules-check,
`archive/0.1`), their schema kinds (`rules-model`, `rules-run`,
`rules-result`, `rules-check`, `model`, `archive`, `archive-check`) and
`legacy.rules`, the rules engine's 0.6.0 discovery entry
(`crn_craft::rules::describe()`), whose `crate` field is now `crn-craft`.

Every JSON object is printed with sorted keys, so "the 1.0 contracts
first" means: they are the top-level fields and the start of `--help`'s
kind list, and the older contracts are inside `legacy`.

The network document is `crn_craft::describe_network()`
(`crn-craft.describe/0.11`, marked legacy), read by
`crn_craft._legacy.describe()`; three Rust tests of it were pointed at the
new name. `describe KIND` still serves every kind, the hidden network ones
(`run`, `result`, `result-v1`, `check`, `fsp-result`, `tau-request`,
`tau-result`, `ssa-run`, `ssa-result`) included, until M8b.

### Check (`crn_craft::check_any`, `check_any_files`)

`check MODEL [REQUEST]` reads every model version and dispatches:

| Input | Report |
| --- | --- |
| rules-model alone | `crn-craft.rules-check/0.1`, byte-identical to 0.6.0's `rules-check` |
| `crn-craft.model/1.0` alone | `crn-craft.check/1.0` (as `rules-check` since M3) |
| network `crn-craft.model/0.1` alone | `crn-craft.check/1.0` of its conversion: valid when `convert` accepts it; otherwise the converter's diagnostics (the network validation errors with their 0.6.0 codes and paths, or `CRN-R016`) |
| any rules or 1.0 model with a request | `crn-craft.check/1.0` with `run_sha256`: the model's check, then exactly the errors `run` reports before running (`rules::prepare`, or the new `fsp::v1::check` for an FSP request), and with files the `model.path` reference (`CRN-R020` at `/model/path`, as `run`) |
| network model with a network `run/0.1` request | the network engine's `crn-craft.check/0.1` (legacy, for the hidden commands; M8b deletes it) |
| network model with any other request | `check/1.0` with `CRN-E033` |

**Decision for `model/0.1`: the converted check.** Reasons: M8b deletes
the network validator's report (`check/0.1`) and keeps the validation only
inside the converter, so the converted check is the behavior that lasts;
for an invalid model it reports the same codes and paths as before (the
converter returns the network validation errors); and it answers what a
user of 0.7.0 needs to know, whether the unified engine can take the model.
What changes: the report is `check/1.0` (no `scope`,
`execution_available` or `run_sha256` fields), the network warnings
`CRN-W001`–`W003` are no longer in a check (`analyze` reports them), and a
valid network model the converter refuses (population limits, no
reactions) is now invalid with `CRN-R016`.

`RulesCheck` gained `run_sha256` (skipped when absent); the `rules-check`
schema removes it (`ONLY_1_0_ROOT`), so `schemas/rules-check.schema.json`
is unchanged; `schemas/check-1.0.schema.json` was regenerated with
`crn-craft describe check-1.0` (new property, new top-level and
`diagnostics` descriptions).

### `CRN-E033`

A network model given where a rules or 1.0 model is read (`run`, `check
MODEL REQUEST`, `crn_craft.rules.prepare`/`ensemble`/`load`, FSP) is now
refused with `CRN-E033` at `/schema` ("convert it first"), instead of
`CRN-E003` at its first network field (for example `/convention`). It is
the shape error of the strict rules-model parse turned into the new code
when the document's `schema` is `crn-craft.model/0.1`, the mirror of M3's
`CRN-E032`. `docs/diagnostics.md` lists it.

### Python command line (`python/crn_craft/__main__.py`)

`python -m crn_craft` (the `crn-craft` script of the wheel) has the same
command set, aliases, hidden legacy commands and usage report; argument
parsing follows the native CLI (no argparse, so the usage message cannot
list hidden commands). `run` prints the exact bytes of the native `run`
(native `run_files_json`); `check` uses the native `check_any_files_json`;
the other reports are the native reports, pretty-printed with sorted keys
as before. `--version` prints `{"version": ..., "interface": "python",
"milestone": "unified-engine"}`.

## Tests

The implementation was drafted first. `crates/crn-craft-cli/tests/cli.rs`
gained four tests (`help_usage_and_describe_list_the_unified_commands_only`,
`run_and_check_are_rules_and_rules_check`,
`check_reads_every_model_version`,
`hidden_network_commands_keep_working_until_m8b`) and four were updated
(the describe assertions moved to `describe/1.0`; `check` left the list of
commands that refuse `model/1.0`). Their first run failed in two tests,
both on the `check-1.0` schema snapshot, which then was regenerated; the
first regeneration also changed the `rules-check` schema (a field doc), so
the request text moved to a 1.0-only field description and the
`rules-check` schema is byte-identical again.

`validation/check_contracts.py` was rewritten for the new commands and run
first against the step-0 binary (`f3fcb6b`): it stopped at the first
example model (`jsonschema.ValidationError`: a `check/0.1` report where
`check/1.0` is required); against the new binary: 126 of 126. It keeps
every strictness check and adds:

- schema snapshots of every kind (1.0, still read, hidden legacy);
- `describe/1.0` fields, no hidden command or schema listed, `--help`,
  `--version`;
- network examples: `check` (converted, `check/1.0`), `convert` (valid
  against `model-1.0`, `check` and `rules-check` agree), `analyze`,
  `analyze-v1` gone, `export`; network run files (LEGACY, `check/0.1`);
- rules models: `check` equals `rules-check` (`rules-check/0.1`);
- the four rules-model/0.3 runs through `run` and `rules` with identical
  bytes, and `check MODEL REQUEST` agreeing with `run` (the refused
  continuous model gives the same diagnostics);
- FSP through `run` and `check`, the missing bound refused identically;
- strict requests through `run` and `check`: unknown field (`CRN-E003`),
  duplicate key (`CRN-E002`), trailing document (`CRN-E001`), above 16 MiB
  and too deep (`CRN-E004`), every refusal exit 1;
- `CRN-E033` for a network model with a `run/1.0` request;
- the invalid network fixtures, model alone against `check-1.0`, with a
  run against the legacy `check` schema;
- schema rejections of the network `model` and of `model-1.0` (unknown
  top-level and pattern fields, wrong and rules-model versions);
- strict model files through `check`: trailing document, non-UTF-8, depth,
  size, duplicate key, unknown field;
- usage errors (`run`, `rules`, `report`, `export`, `check`, no argument),
  whose hint names no hidden command; deterministic reports.

`check_python.py` (`command_line`) checks the Python CLI against the native
one: `--help` commands, `describe` equal to the native document and to
`crn.describe()`, `run` and `rules` printing the native bytes, `check` and
`rules-check` reports equal, `convert` bytes equal, the usage error of a
hidden command without the hidden names, `--version`.

## Checks

See the [M6a log](M6a-log.md#checks): one gate covers both tasks.

## Deviations from the plan

- `analyze` keeps its optional STATE argument.
- `check` of a network model with a network `run/0.1` request keeps the
  network `check/0.1` report until M8b (the hidden commands' contract).
- The plan says the old commands are "removed"; as instructed they are
  hidden until M8b.
- New diagnostic `CRN-E033` (not in the plan).
- `describe` is one document for the CLI and both Python names
  (`describe/1.0`); the core `describe/0.11` stays as
  `crn_craft::describe_network()` for `crn_craft._legacy`.

## Migration notes (command line)

| 0.6.0 | 0.7.0 |
| --- | --- |
| `crn-craft rules MODEL REQUEST [--archive DIR] [--cancel-file PATH]` | `crn-craft run ...` (`rules` stays an alias) |
| `crn-craft rules-check MODEL` | `crn-craft check MODEL [REQUEST]` (`rules-check` stays an alias); same report for a rules model |
| `crn-craft check MODEL` (network, `check/0.1`) | `crn-craft check MODEL`: `check/1.0` of the conversion |
| `crn-craft check MODEL RUN` (network) | unchanged until M8b (hidden contract); then convert the model and write a `run/1.0` request |
| `crn-craft simulate MODEL RUN DIR`, `simulate-cached` | `crn-craft convert MODEL > model.json`, a `run/1.0` request, `crn-craft run model.json request.json --archive DIR`; hidden until M8b |
| `crn-craft ssa MODEL REQUEST` | `crn-craft run` with `"method": {"kind": "static_direct", "selector": ...}` or `next_reaction`; hidden until M8b |
| `crn-craft tau MODEL REQUEST` | `crn-craft run` with `"method": {"kind": "tau", "epsilon": ...}`; hidden until M8b |
| `crn-craft fsp MODEL RUN` | `crn-craft run MODEL REQUEST` with a `run/1.0` request whose method is `fsp`; hidden until M8b |
| `crn-craft describe` (`describe/0.15` native, `0.11` Python) | `crn-craft.describe/1.0`, the same document from both; old fields of the rules engine under `legacy.rules` |
| `crn-craft --version` milestone `SSA-acceleration` | `unified-engine` |
| `CRN-E003` for a network model given to `rules`/`rules-check` | `CRN-E033` at `/schema` |
| Rust `crn_craft::describe()` (`describe/0.11`) | `crn_craft::describe()` is `describe/1.0`; the old one is `crn_craft::describe_network()` until M8b |
| — | Rust `crn_craft::check_any`, `check_any_files`, `AnyCheck`, `crn_craft::unified::{MILESTONE, COMMANDS, SCHEMAS, LEGACY_SCHEMAS}`, `crn_craft::fsp::v1::check`, `RulesCheck::run_sha256` |
| Rust `RulesCheck { schema, valid, model_sha256, diagnostics }` | struct literals add `run_sha256: None` |
