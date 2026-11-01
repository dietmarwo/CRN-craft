# Diagnostics

Every diagnostic contains `severity`, `code`, `path`, `message`, `hint`. Paths use RFC 6901 (`~1` for slash, `~0` for tilde); the empty path denotes the root. Type errors may identify the containing object when the parser cannot name a missing/unknown child. Messages provide detail; code and relevant pointer are the machine contract.

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
| CRN-E031 | Unsupported command or wrong arguments | Use `--help`; only describe/check execute in M0 |
| CRN-W001 | Net-producing reaction without a net-consuming one | Review growth and budgets; no proof of unboundedness is implied |
| CRN-W002 | Zero mass-action coefficient | Confirm the channel is intentionally disabled |
| CRN-W003 | Reaction has no net population effect | Remove it unless event counts have a stated purpose |

The 13 persisted negative fixtures in `models/invalid/expected.json` pin representative codes and paths. Rust tests cover additional boundary and parser cases. Parsing stops on the first syntax/type error per document; semantic validation accumulates errors after a successful parse. A failed stdout write can only report to stderr and exit 1; it cannot guarantee a JSON response on a broken stream.
