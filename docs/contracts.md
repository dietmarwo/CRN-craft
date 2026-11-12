# Contracts

Every document that CRN-craft reads or writes is versioned JSON with a
`schema` field. The Rust types define the shapes. The files in `schemas/` and
`crn-craft describe KIND` are generated from them, so they cannot disagree.
`crn-craft describe` without a kind prints the discovery document
(`crn-craft.describe/1.0`): commands, contracts, methods, selectors,
quantiles and limits. The Python call `crn_craft.describe()` returns the same
document.

## Current contracts

| Contract | Role | `describe` kind |
| --- | --- | --- |
| `crn-craft.model/1.0` | the model: [model reference](rules-language.md) | `model-1.0` |
| `crn-craft.run/1.0` | the request: simulation methods, or `fsp` | `run-1.0` |
| `crn-craft.result/1.0` | the result of a simulation request ([result contract](result-contract.md)) | `result-1.0` |
| `crn-craft.fsp-result/1.0` | the result of an `fsp` request | `fsp-result-1.0` |
| `crn-craft.check/1.0` | the report of `check`: whether `run` would accept the documents | `check-1.0` |
| `crn-craft.analysis/1.0` | structure, conservation laws and bounds of a model | `analysis-1.0` |
| `crn-craft.archive/1.0`, `crn-craft.archive-check/1.0` | a sealed run and its verification ([archive contract](archive-contract.md)) | `archive-1.0`, `archive-check-1.0` |
| `crn-craft.import/1.0` | the report of the SBML and BioNetGen importers; its `model` is a `model/1.0` document | none |
| `crn-craft.export/0.1` | the report of `export MODEL gillespy2` | `export` |
| `crn-craft.state/0.1` | an optional state for `analyze MODEL STATE`, by species label | `state` |
| `crn-craft.report/1.0` | the JSON file written by `Result.save`: model, request and result in one document. Written and read by Python only, so it has no schema file | none |
| `crn-craft.describe/1.0` | the discovery document | none |

## Contracts that are still read

| Contract | How it is read |
| --- | --- |
| `crn-craft.rules-model/0.1` to `0.4` | as before. A `rules-model/0.4` document means what `model/1.0` means, with the older limits. |
| `crn-craft.rules-run/0.1` to `0.6` | as before; the results are `rules-result/0.1` to `0.3` and replay byte for byte. |
| `crn-craft.rules-check/0.1` | the report of `check MODEL` for a rules model alone |
| network `crn-craft.model/0.1` | through the converter: `check`, `analyze`, `convert` and `export` read it, and `convert` writes it as `model/1.0`. `run` refuses it with `CRN-E033`. |

A 1.0 model needs a 1.0 request: a `model/1.0` model with a `rules-run`
request is refused (`CRN-R020` at `/schema`). A `rules-model` model runs
under either request family. The `rules-*` request versions refuse the
fields that arrived with 1.0 (`method.selector`, `outputs.quantiles`, the
`fsp` method), also with `CRN-R020` at `/schema`.

Documents of the removed network engine are not read any more: network
`run/0.1`, `ssa-run/0.1` to `0.4`, `tau-request/0.1`, `result/0.1`,
`ssa-result`, `tau-result`, `fsp-result/0.1`, `analysis/0.1` and `0.2`,
`archive/0.1` and `archive-check/0.1`. Their schema files are gone. The
[changelog](../CHANGELOG.md) has the migration tables.

## Strict input

- A file is at most 16 MiB. JSON nesting is at most 48 deep, with the root
  at depth 0.
- Duplicate keys are rejected, also when JSON escapes give two spellings of
  one key. Unknown fields are rejected at every depth. `NaN` and `Infinity`
  are not JSON, and a literal such as `1e400` fails with a numeric
  diagnostic.
- Decimal literals parse to the correctly rounded double, as Rust's
  `str::parse` and Python's `json` do, so the shortest round-trip text of a
  double reads back to the same bits.
- Syntax and shape are checked first. Semantic checks follow and report
  every error they find, each with a stable code and a JSON path
  ([diagnostics](diagnostics.md)).
- Hashes are SHA-256 of the exact bytes consumed, never of re-serialized
  JSON. Line endings are part of the bytes: `.gitattributes` keeps LF
  endings in every checkout, and CI runs the contract suite on a checkout
  made with `core.autocrlf=true`.
- A request names its model by a path relative to the request's directory
  and by the model's SHA-256. The path must stay inside that directory after
  resolving symbolic links: no empty components, absolute paths, `.`, `..`,
  drive prefixes, backslashes or NUL bytes.

## Limits

`describe()["limits"]` lists the limits of the current contracts. Models of
`model/1.0` accept 4,096 types, 16,384 parameters, 4,096 observables and
65,536 rules; the static methods accept 4,096 instance tuples and 65,536
channels. The `rules-model` versions keep 64, 1,024, 256, 1,024, 256 and
4,096. See the [model reference](rules-language.md#limits).

## Changing a contract

A change to the shape of a delivered contract needs a new version, and the
old version keeps its meaning. A new field in a new version is refused by
the versions before it. Schemas are generated, never edited by hand.

## Where the old pages went

The earlier pages on the network model, its run request, the reserved result
envelope and the separate analysis reports described contracts that the
unified engine replaced. Their successors are:

| Earlier topic | Now |
| --- | --- |
| network `model/0.1`, `run/0.1`, `state/0.1` | [model reference](rules-language.md); `crn-craft convert` writes an old model as `model/1.0` |
| `analysis/0.1`, `0.2` and `analyze-v1` | `crn-craft.analysis/1.0` ([network analysis](network-api.md)) |
| `result/0.1` and `0.2`, `archive/0.1` | [result contract](result-contract.md), [archive contract](archive-contract.md) |
| `ssa-run`, `tau-request` | `crn-craft.run/1.0` methods `static_direct`, `next_reaction` and `tau` ([methods](direct-ssa.md)) |
