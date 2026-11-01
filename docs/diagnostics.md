# Diagnostics

Every diagnostic contains `severity`, `code`, `path`, `message`, `hint`. Paths use RFC 6901 (`~1` for slash, `~0` for tilde); the empty path denotes the root. Type errors may identify the containing object when the parser cannot name a missing/unknown child. Inside the tagged `rate` and `method` variants, a missing field or an invalid field value (for example `"aggregator": "nrm"`) is reported at the variant object, such as `/method`; an unknown `kind` is reported at `.../kind`. Messages provide detail; code and relevant pointer are the machine contract.

| Code | Meaning | Typical repair |
|---|---|---|
| CRN-E001 | Invalid JSON, encoding or trailing content | Supply one UTF-8 JSON document |
| CRN-E002 | Duplicate JSON key at any depth | Retain one value per decoded key |
| CRN-E003 | Wrong type, unknown/missing field, schema version or enum variant | Use the generated schema |
| CRN-E004 | Input depth/size/model collection limit or no species | Reduce size/depth; provide at least one species |
| CRN-E010 | Invalid identifier | Use the documented ASCII identifier format |
| CRN-E011 | Duplicate name or species/parameter collision | Rename the conflicting entry |
| CRN-E012 | Unknown species | Declare or correct it |
| CRN-E013 | Non-positive or oversized integer stoichiometry | Use 1–2147483647; omit zeros |
| CRN-E014 | Unknown parameter | Declare it or provide a numeric literal |
| CRN-E015 | Non-finite or out-of-domain numeric argument | Use the required finite positive/non-negative domain |
| CRN-E016 | Negative population or insufficient safety limit | Correct counts/limits |
| CRN-E020 | Invalid time grid, budget, ensemble, output or method combination | Follow the field hint |
| CRN-E021 | Malformed or mismatched SHA-256 | Review model bytes and compute their lowercase digest |
| CRN-E022 | Unsafe/non-relative path, missing target or mismatch with supplied model | Put the model below the run directory and pass that file |
| CRN-E023 | Missing/invalid FSP bounds, capacity or tolerance | Supply a finite truncation containing the initial state |
| CRN-E030 | Input/output filesystem or stream failure | Supply a readable regular file or writable output stream |
| CRN-E031 | Unsupported command or wrong arguments | Use `--help`; describe/check/analyze/analyze-v1/export execute in M1b |
| CRN-E040 | Dependency graph exceeds 1000000 edges | Reduce network coupling or size; edges are never silently dropped |
| CRN-E041 | Invalid state length, missing/extra name, count or safety limit | Provide all species with valid non-negative i64 counts |
| CRN-E042 | Evaluated rate overflows, is non-finite/negative, or a positive rate becomes zero | Revise/rescale the model; do not treat underflow as an absorbing channel |
| CRN-E043 | Invalid reaction index or disabled reaction update | Select an in-range channel with positive propensity and sufficient reactants |
| CRN-E044 | Reaction update overflows i64 or exceeds a population limit | Stop or explicitly revise the safety limit; state remains unchanged |
| CRN-E050 | Export target cannot preserve a requested feature | Use the documented elementary mass-action subset; keep safety semantics in the source |
| CRN-E051 | Homodimer coefficient cannot be halved exactly in f64 | Rescale the reference fixture; no source is emitted |
| CRN-W001 | Net-producing reaction without a net-consuming one | Review growth and budgets; no proof of unboundedness is implied |
| CRN-W002 | Zero mass-action coefficient | Confirm the channel is intentionally disabled |
| CRN-W003 | Reaction has no net population effect | Remove it unless event counts have a stated purpose |
| CRN-W010 | Exact conservation analysis exceeds 64 species / 256 reactions | Use a smaller network or independent analyzer; no absence of laws is implied |
| CRN-W011 | Checked i128 conservation arithmetic exceeds its range | Use arbitrary precision; no partial basis or rank is reported |
| CRN-W012 | Cone exceeds 64 species / 256 reactions | Use an independent analyzer; no absence of laws is implied |
| CRN-W013 | Cone exhausts intermediate-ray or deterministic-work cap | Inspect cone status/limits; only sound basis bounds remain |
| CRN-W014 | Checked i128 cone arithmetic overflows | Use arbitrary precision; no partial rays or cone bounds are emitted |
| CRN-W015 | Internal exact cone residual verification fails | Report a defect; no cone result is trusted |

The 13 persisted negative fixtures in `models/invalid/expected.json` pin representative codes and paths. Rust tests cover additional boundary and parser cases. Parsing stops on the first syntax/type error per document; semantic validation accumulates errors after a successful parse. A failed stdout write can only report to stderr and exit 1; it cannot guarantee a JSON response on a broken stream.

State paths such as `/counts/A` refer to the supplied state, while rate and reaction paths refer to the model. `check` retains structural CRN-W001 warnings. `analyze` suppresses one only when its species has a proven initial-class conservation bound. Conservation warnings appear both in the conservation report and the top-level diagnostic list for direct consumers.

## M2 runtime and archive diagnostics

| Code | Meaning | Response |
|---|---|---|
| CRN-E060 | Unsupported execution or resource/accounting reservation exceeded | Use Direct SSA; reduce workers, replicates or retained output; budget room for accounting |
| CRN-E061 | Invalid/unresolved total or categorical probability | Rescale or revise dynamic range; no positive channel is silently dropped |
| CRN-E062 | Invalid exponential wait or stagnant/nonfinite event time | Rescale rates/time; do not clip or force time forward |
| CRN-E064 | Actual serialized output exceeds the cap | Increase max_output_bytes; optional numerical output is removed and accounting retained when it fits |
| CRN-E065 | Archive I/O, integrity, size or relational failure | Preserve evidence; use a new output directory and verify the complete file set |

CRN-E042 remains a numerical failure for positive-to-zero propensity underflow.
CRN-E044 retains its checked update meaning: population-limit attempts stop
without applying the event; integer overflow is a failed trajectory. Ordinary
event, deadline and cancellation stops are typed termination states, not
fabricated numerical errors. A valid incomplete archive remains incomplete.

Campaign codes: CRN-E070 covers FSP capacity, work, normalization, output and
interruption failures; CRN-E080 covers optional deterministic domains/unsupported
smooth extensions/solver failures; CRN-E085 covers approximate Poisson-leap
requests, budgets and numerical failures; CRN-E090 covers importer syntax,
unsupported semantics, explicit-convention and resource failures. None denotes
biological validation or a certified floating-point error bound.
