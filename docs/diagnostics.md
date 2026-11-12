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
| CRN-E014 | Unknown parameter; for parameter templates, a wrong value count or name set | Declare it or provide a numeric literal; bind exactly the declared names |
| CRN-E015 | Non-finite or out-of-domain numeric argument | Use the required finite positive/non-negative domain |
| CRN-E016 | Negative population or insufficient safety limit | Correct counts/limits |
| CRN-E020 | Invalid time grid, budget, ensemble, output or method combination; at the API path `/kernel`, a generated kernel for another topology or kernel generation for non-mass-action laws | Follow the field hint; regenerate the kernel for this model |
| CRN-E021 | Malformed or mismatched SHA-256 | Review model bytes and compute their lowercase digest |
| CRN-E022 | Unsafe/non-relative path, missing target or mismatch with supplied model | Put the model below the run directory and pass that file |
| CRN-E023 | Missing/invalid FSP bounds, capacity or tolerance | Supply a finite truncation containing the initial state |
| CRN-E030 | Input/output filesystem or stream failure | Supply a readable regular file or writable output stream |
| CRN-E031 | Unsupported command or wrong arguments | Use `--help` for the command list |
| CRN-E032 | A rules-engine document (`crn-craft.model/1.0`, `crn-craft.run/1.0` or a `crn-craft.rules-*` contract) given to a network command or reader, reported at `/schema` | Use `rules-check` and `rules`; `convert` writes a network model as `crn-craft.model/1.0` |
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
| CRN-E067 | Sampling-word budget exhausted in an alternative SSA method (`ssa`; recorded as `event_limit` for that replicate) | Increase max_sampling_steps_per_replicate; incomplete ensembles withhold statistics |

CRN-E042 remains a numerical failure for positive-to-zero propensity underflow.
CRN-E044 retains its checked update meaning: population-limit attempts stop
without applying the event; integer overflow is a failed trajectory. Ordinary
event, deadline and cancellation stops are typed termination states, not
fabricated numerical errors. A valid incomplete archive remains incomplete.
CRN-E043 also reports a reaction that would make a count negative: in the
Direct loop this is a scheduler defect, never a population-limit stop. The
optional `crn-diffsol` integration reports a species exceeding its
`population_limit` between or at requested times with CRN-E044.

Campaign codes: CRN-E070 covers FSP capacity, work, normalization, output and
interruption failures, for the network `fsp` command and for `crn-craft.run/1.0`
requests with method `fsp`; CRN-E080 covers optional deterministic domains/unsupported
smooth extensions/solver failures; CRN-E085 covers approximate Poisson-leap
requests, budgets and numerical failures; CRN-E090 covers importer syntax,
unsupported semantics, explicit-convention and resource failures. None denotes
biological validation or a certified floating-point error bound.

## Rules engine (`CRN-R…`)

The experimental rules engine for dynamic, attributed species
([rules language](rules-language.md)) shares `CRN-E001`–`CRN-E004` for JSON
input and the read error `CRN-E030`, and has its own codes: `CRN-R010` names,
`CRN-R011` unknown types, `CRN-R012` expressions and scopes, `CRN-R013` rule
shapes and product attributes (in `crn-craft.model/1.0` also a pattern
`count` outside 1 to 2,147,483,647 and a higher-order pattern over an
attributed type or with `where` or `factor`), `CRN-R014` initial state, `CRN-R015` model
limits (64 types, 1,024 parameters, 256 observables and 1,024 rules up to
`rules-model/0.4`; 4,096 types, 16,384 parameters, 4,096 observables and
65,536 rules in `crn-craft.model/1.0`), `CRN-R016` a network model that the
converter (`crn-craft convert`, `crn_craft.convert`) cannot write as a rules
model yet (no reaction, more than 64 product molecules in a reaction, an
initial count above 2^53, or a population limit; the path points into the
network model), `CRN-R017` a model that has no flat view, which finite
state projection (and from 0.7.0 analysis, ODE and export) needs: a model
whose products can create attribute tuples during a run, with random
product attributes, a whole-rule propensity over a type with several tuples
or above the static lowering limits (at `/rules`), with an `atp` block (at
`/atp`) or with continuous variables (at `/continuous`),
`CRN-R020` run requests (including method, tracker, log and version
gates, a request with method `fsp` where a simulation is expected (at
`/method`), the `fsp` request's version (at `/schema`), `max_states`,
`tolerance` and bounds (a missing bound, a bound below the initial count,
or an unknown species label, at `/method/bounds/LABEL`), a `crn-craft.model/1.0` model with a `rules-run` request, a
higher-order rule with `tree`, `mtdna_direct`, `binomial` or a generated
rates kernel, `outputs.quantiles` in a request before `crn-craft.run/1.0`
(at `/schema`), more than 16 or not strictly increasing quantile
probabilities in (0, 1), a quantile buffer above the 2 GiB workspace limit
(at `/outputs/quantiles`), `method.selector` in a request before
`crn-craft.run/1.0` (at `/schema`), a generated event loop or rates kernel
with a Static Direct selector other than `linear`, also when `auto` picks one
(at `/method/selector`), and a generated-kernel fingerprint mismatch),
`CRN-R030` numerical
failures during a replicate (with the selectors of `crn-craft.run/1.0` also
a non-finite total, a composition or acceptance draw rounded to its bound,
and a propensity outside its RSSA-CR bounds), `CRN-R031` weights below
53-bit resolution (in the tree a rule or instance weight; in Static Direct
and its fallbacks, with every selector, a channel below 2⁻⁵³ of the total or
lost in the total's pairwise sums; with RSSA-CR the smallest positive
propensity of a channel's population box, so use `tree` or rescale),
`CRN-R040`–`CRN-R045` event, sampling, instance, wall-time, cancellation and
output limits, and the warnings `CRN-R050` for an incomplete ensemble and
`CRN-R051` for a requested specialization that did not apply (an
`mtdna_direct` fallback, or an explicit `sparse_trackers: true` without an
eligible tracker; the rules-run/0.6 default does not warn).
