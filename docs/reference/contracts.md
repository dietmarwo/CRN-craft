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
| `crn-craft.model/1.0` | the model: [model reference](model-language.md) | `model-1.0` |
| `crn-craft.run/1.0` | the request: simulation methods, or `fsp` | `run-1.0` |
| `crn-craft.result/1.0` | the result of a simulation request ([results](results.md)) | `result-1.0` |
| `crn-craft.fsp-result/1.0` | the result of an `fsp` request | `fsp-result-1.0` |
| `crn-craft.check/1.0` | the report of `check`: whether `run` would accept the documents | `check-1.0` |
| `crn-craft.analysis/1.0` | structure, conservation laws and bounds of a model | `analysis-1.0` |
| `crn-craft.archive/1.0`, `crn-craft.archive-check/1.0` | a sealed run and its verification ([archive contract](archives.md)) | `archive-1.0`, `archive-check-1.0` |
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
[changelog](../../CHANGELOG.md) has the migration tables.

## Versions of the request and the model

The request version decides which fields are allowed:

| Version | Adds |
| --- | --- |
| `crn-craft.run/1.0` | `method.selector` for `static_direct`, `outputs.quantiles`, the method `fsp`; Static Direct defaults to the selector `auto`; the result is `crn-craft.result/1.0` |
| `crn-craft.rules-run/0.1` | the tree method only |
| `crn-craft.rules-run/0.2` | `method`, `trackers`, `event_log`; the result is rules-result/0.2, with a sampler identity and optional ordered `observations` |
| `crn-craft.rules-run/0.3` | `mtdna_direct` |
| `crn-craft.rules-run/0.4` | Static Direct `waiting_time` and `sparse_trackers` |
| `crn-craft.rules-run/0.5` | `hybrid` with its `rtol` and `atol`; the result is rules-result/0.3, which adds `method` and `continuous` |
| `crn-craft.rules-run/0.6` | no field; Static Direct defaults to Ziggurat waits and sparse trackers; the result is rules-result/0.3 |

`crn-craft.run/1.0` has the shape of `rules-run/0.6` and means what it means,
plus the features of 1.0. Python writes the versions as follows.
`crn.ensemble` always writes `run/1.0`. `crn.rules.ensemble` picks the lowest
`rules-run` version that fits its arguments, and rules-run/0.6 for
`static_direct`; it writes `run/1.0` only when an argument needs it
(`quantiles=`, `selector=`) or the model is `crn-craft.model/1.0`.
`crn.prepare(model, request)` takes a request you wrote yourself.

Model files carry their own version:

| Version | Adds |
| --- | --- |
| `crn-craft.rules-model/0.1` | the base language |
| `crn-craft.rules-model/0.2` | the legacy `atp` block |
| `crn-craft.rules-model/0.3` | whole-rule `propensity` and `continuous` variables; from 0.3 a rule gives exactly one of `rate` and `propensity` |
| `crn-craft.rules-model/0.4` | `"independent": true` for continuous variables whose derivative reads only parameters and other independent variables; `hybrid` integrates them once per prepared run and shares the solution among all replicates |
| `crn-craft.model/1.0` | what 0.4 means, with larger limits and the pattern `count` |

A 0.4 document without the field means what it means as 0.3. `crn.Model`
writes `model/1.0`, and `crn.rules.Model` writes 0.3 (0.4 with an independent
ODE). A `rules-model` file runs under either request family. Stored
rules-model 0.1 to 0.3 files and rules-run 0.1 to 0.6 requests keep their
meaning and replay byte for byte. Each backend has its own replay identity.
The worker count preserves a completed numerical payload apart from its
request hash.

## The run request

`crn-craft.run/1.0` has the shape of `rules-run/0.6` and adds three things:
`method.selector` for `static_direct` (`linear`, `composition_rejection`,
`sorting`, `rssa_cr` or `auto`, the default; see [exact methods](methods.md)),
`outputs.quantiles` (up to 16 probabilities in (0, 1); the result's
`statistics.quantiles` holds the quantile of every report time and observable,
and the median is the quantile at 0.5; see the [result contract](results.md))
and the method `fsp` (finite state projection: `bounds`, `max_states` and
`tolerance`; the result is `crn-craft.fsp-result/1.0`). Its result is
`crn-craft.result/1.0`. A `rules-run` request refuses these fields (`CRN-R020`
at `/schema`) and keeps its meaning, its numbers and its hashes. The rest of
this section describes the common core, starting with
`crn-craft.rules-run/0.1`, which binds the model by relative path and SHA-256 and
gives `time.end`, strictly increasing `time.report_times` in [0, end],
`ensemble.replicates` and `root_seed`, `outputs` (`statistics` from mean and
variance, `keep_paths`, `keep_instances`) and an explicit `budget`
(`max_events_per_replicate`, `max_sampling_steps_per_replicate` counting every
random word, `max_instances` alive at once, `threads` 1–64,
`max_output_bytes`, optional `wall_time_s`). Before running, the request is
admitted only if a conservative reservation for records, in-flight and kept
series and kept final instances (`keep_instances × max_instances`) stays below
2 GiB, and if `max_output_bytes` can hold at least 160 bytes per replicate
record (`CRN-R020`).

Trackers and event logs (rules-run/0.2 and later) reserve output per
replicate: 256 bytes per tracker plus 64 per kept episode, and `max_entries`
times a bound on one log entry. The entry bound counts every number at its
longest compact JSON form (25 bytes with its comma) and uses the model's own
sizes: 160 bytes for the keys, time and rule, 25 per observable column, and
per participant of an exact step 62 bytes plus 25 per attribute, where the
participants are the largest number of reactant patterns plus products of
any rule. Leap entries (`tau`, `binomial`) list firings instead, 45 bytes per
lowered channel. Up to 0.2.0 the bound assumed two reactants and 64 products
for every rule; the smaller of the two bounds applies, so no request admitted
before is refused now. For the three-class mtDNA model an entry reserves 721
bytes, so 50 entries for 1,000 cells take 36 MB of the 64 MiB Python default.
Replicates times this reservation must fit `max_output_bytes` (and 2 GiB). The
`CRN-R020` message at `/budget` itemizes the reservation, and its hint names
the smallest `max_output_bytes` that admits the request.

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
4,096. See the [model reference](model-language.md#limits).

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
| network `model/0.1`, `run/0.1` | [model reference](model-language.md); `crn-craft convert` writes an old model as `model/1.0`. `state/0.1` stays current (above) |
| `analysis/0.1`, `0.2` and `analyze-v1` | `crn-craft.analysis/1.0` ([network analysis](analysis.md)) |
| `result/0.1` and `0.2`, `archive/0.1` | [result contract](results.md), [archive contract](archives.md) |
| `ssa-run`, `tau-request` | `crn-craft.run/1.0` methods `static_direct`, `next_reaction` and `tau` ([methods](methods.md)) |
