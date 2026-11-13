# Diagnostics

Every diagnostic contains `severity`, `code`, `path`, `message` and `hint`.
Paths use RFC 6901 (`~1` for a slash, `~0` for a tilde); the empty path
denotes the root. A type error may name the containing object when the parser
cannot name a missing or unknown child. Messages give detail; the code and the
path are the machine contract. Codes are stable: a number is never reused for
another meaning.

Codes with the prefix `CRN-E` belong to documents, files, the sampler and the
tools around the engine. Codes with the prefix `CRN-R` belong to models and
run requests; they are listed below the first tables. A failed `check`
reports every error it finds. Parsing stops at the first syntax or type error
of a document, and semantic validation accumulates.

## Documents, files and tools (`CRN-E…`)

| Code | Meaning | Typical repair |
| --- | --- | --- |
| CRN-E001 | Invalid JSON, encoding or trailing content | Supply one UTF-8 JSON document |
| CRN-E002 | Duplicate JSON key at any depth | Keep one value per decoded key |
| CRN-E003 | Wrong type, unknown or missing field, schema version or enum variant | Use the generated schema (`crn-craft describe KIND`) |
| CRN-E004 | Input depth or size limit, or a model collection limit | Reduce size or depth |
| CRN-E015 | Non-finite or out-of-domain number | Use the required finite domain |
| CRN-E030 | Input or output filesystem or stream failure | Supply a readable regular file or a writable stream |
| CRN-E031 | Unsupported command or wrong arguments | Use `--help` |
| CRN-E033 | A network model (`crn-craft.model/0.1`) given where a rules or 1.0 model is read (`run`, `check MODEL REQUEST`, `crn_craft.rules.prepare`, finite state projection), reported at `/schema` | Convert it with `crn-craft convert MODEL` or `crn_craft.convert`; `crn_craft.ensemble`, `crn_craft.prepare` and `Model.from_dict` convert it themselves |
| CRN-E041 | Invalid state for `analyze` or the projection: a missing or extra name, or an invalid count | Name every species with a valid count |
| CRN-E044 | `crn-diffsol` only: the ODE trajectory of a network model passes a species' `population_limit` | Stop earlier or revise the limit; the ODE never clips a state |
| CRN-E050 | The export target cannot preserve a requested feature | Keep to the documented subset |
| CRN-E051 | An exported coefficient cannot be halved exactly in f64, reported at `/rules/i/rate` | Rescale the model; no source is emitted |
| CRN-E065 | Archive I/O, integrity, size or relational failure | Keep the evidence; use a new directory and verify the complete file set |
| CRN-E070 | Finite state projection: capacity, work, normalization, output or interruption | Lower the bounds or raise `max_states` within the cap |
| CRN-E080 | ODE right-hand side: domain, unsupported smooth extension or solver failure; a negative propensity | Revise the model; see `crn-diffsol` |
| CRN-E090 | Importer: syntax, unsupported semantics, missing explicit convention or resource failure | Follow the message; pass `convention=` |

None of these codes denotes biological validation or a certified
floating-point error bound.

### Reading a network model

`crn-craft check`, `analyze`, `convert` and `export` read a network
`crn-craft.model/0.1` through the converter. Its checks report these codes:

| Code | Meaning |
| --- | --- |
| CRN-E010 | Invalid identifier |
| CRN-E011 | Duplicate name, or a species and parameter with one name |
| CRN-E012 | Unknown species name |
| CRN-E013 | Stoichiometry outside 1 to 2147483647 |
| CRN-E014 | Unknown parameter name |
| CRN-E016 | Negative population or insufficient population limit |
| CRN-E032 | The document given to the converter (`crn-craft convert`, `crn_craft.convert`) is a `crn-craft.model/1.0` or `crn-craft.rules-*` document, which needs no conversion; reported at `/schema` |

The converter reports a model that the unified model cannot hold with
`CRN-R016`.

### Warnings of analysis

| Code | Meaning |
| --- | --- |
| CRN-W001 | A species has a net-producing reaction and no net-consuming one (at `/species/i`); `analyze` drops it for a species with a proven conservation bound. It does not prove unboundedness |
| CRN-W002 | A reaction's rate is the constant 0 (at `/reactions/i`) |
| CRN-W003 | A reaction has no net population effect (at `/reactions/i`) |
| CRN-W010 | Exact conservation analysis exceeds 64 species or 256 reactions: no laws are claimed or excluded |
| CRN-W011 | Checked i128 conservation arithmetic exceeds its range: no partial basis is reported |
| CRN-W012 | The nonnegative cone exceeds 64 species or 256 reactions |
| CRN-W013 | The cone exhausts its intermediate-ray or work cap: only sound basis bounds remain |
| CRN-W014 | Checked i128 cone arithmetic overflows: no partial rays or bounds are emitted |
| CRN-W015 | Internal verification of the cone failed: report a defect, no cone result is trusted |

### Retired codes

The network engine is gone, and so are the codes that only it reported. They
are not reused: `CRN-E020`, `CRN-E021`, `CRN-E022`, `CRN-E023`, `CRN-E040`,
`CRN-E042`, `CRN-E043`, `CRN-E060`, `CRN-E061`, `CRN-E062`, `CRN-E064`,
`CRN-E067` and `CRN-E085`. `CRN-E044` is left only in `crn-diffsol`. Their
replacements are the `CRN-R` codes: a run request, budget or method problem is `CRN-R020`, a
numerical failure `CRN-R030` and `CRN-R031`, a budget stop `CRN-R040` to
`CRN-R045`, and a document of the wrong family `CRN-E033` or `CRN-R020` at
`/schema`.

## Models and run requests (`CRN-R…`)

The engine for dynamic, attributed species ([model reference](model-language.md))
shares `CRN-E001`–`CRN-E004` for JSON input and the read error `CRN-E030`, and
has its own codes: `CRN-R010` names,
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
and a propensity outside its RSSA-CR bounds; at `/time` a clock that fails to
advance: every `rules-run` version fails at the first positive wait below
half an ulp of the time, Static Direct under `crn-craft.run/1.0` accepts such
a wait and fails only after 1,000,000 consecutive steps without an advance), `CRN-R031` weights below
53-bit resolution (in the tree a rule or instance weight; in Static Direct
and its fallbacks, with every selector, a channel below 2⁻⁵³ of the total or
lost in the total's pairwise sums; with RSSA-CR the smallest positive
propensity of a channel's population box, so use `tree` or rescale),
`CRN-R040`–`CRN-R045` event, sampling, instance, wall-time, cancellation and
output limits, and the warnings `CRN-R050` for an incomplete ensemble and
`CRN-R051` for a requested specialization that did not apply (an
`mtdna_direct` fallback, or an explicit `sparse_trackers: true` without an
eligible tracker; the rules-run/0.6 default does not warn).

### The `CRN-R` codes at a glance

| Code | Meaning |
|---|---|
| CRN-E001–E004 | malformed JSON, duplicate keys, unknown fields or wrong types (`E003`), size limits |
| CRN-R010 | invalid, reserved or duplicate name; observable shadowing a parameter |
| CRN-R011 | unknown type |
| CRN-R012 | expression syntax, unknown name, out-of-scope name or random draw outside products; an independent variable's derivative that reads an observable or a dependent variable |
| CRN-R013 | unsupported rule shape, missing or unknown attribute, bad `copy`; both or neither of `rate` and `propensity`, `where` or `factor` on a propensity rule's reactant; `lower_bound` without `on_bound` or the reverse |
| CRN-R014 | invalid initial count or attribute; a continuous variable's initial value or bound that is not finite, or an initial value below its bound |
| CRN-R013 and R020 | for `model/1.0`: a pattern `count` outside 1 to 2,147,483,647, or on an attributed type or with `where` or `factor` (`R013`); a higher-order rule under a method that refuses it (`R020` at `/method`) |
| CRN-R015 | model size limit (the [limits](model-language.md#limits); also 64 continuous variables); a field used before its model version (`atp` in 0.1, `propensity` and `continuous` in 0.1 or 0.2, `independent` before 0.4); `atp` together with `continuous` |
| CRN-R016, R017 | a network model the converter cannot write (`R016`); a model without a flat view, which finite state projection, analysis, ODE export and GillesPy2 export need (`R017`) |
| CRN-R020 | invalid run request (hash, path, grid, ensemble, budget, method, selector, quantiles, tracker, log or version gate; generated-kernel fingerprint mismatch) |
| CRN-R030 | negative or non-finite rate, factor, weight or product attribute; overflow; under `run/1.0`, 1,000,000 consecutive waits that do not advance the clock |
| CRN-R031 | weight below 53-bit resolution (Static Direct: a channel below 2⁻⁵³ of the total) |
| CRN-R040–R045 | event, sampling, instance (including a count above 2⁵³), wall-time, cancellation and output limits |
| CRN-R050 | warning: incomplete ensemble, statistics withheld |
| CRN-R051 | warning: a requested specialization did not apply (`mtdna_direct` could not certify the model and ran Static Direct, or an explicit `sparse_trackers: true` found no eligible tracker; the rules-run/0.6 default does not warn); the message gives the reason |
