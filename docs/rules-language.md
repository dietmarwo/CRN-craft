# Dynamic, attributed species: the rules language

**Status:** experimental (K1). Contracts `crn-craft.rules-model/0.1`,
`crn-craft.rules-check/0.1`, `crn-craft.rules-run/0.1`,
`crn-craft.rules-result/0.1`; crate `crn-craft-rules`; CLI `rules-check` and
`rules`. Protocol and gates: `docs/decisions/rules-k0-k1-2026-11-02.md`.

The network engine needs every species listed in advance. Some models cannot
do that: each new mitochondrial mutant is its own clone, polymers have a
length, clusters have a size, cells carry a lineage. The rules engine
simulates such systems exactly. Species *types* carry numeric *attributes*;
an *instance* is a type plus one value per attribute and a copy number.
Instances appear when a rule produces a new attribute vector and disappear
when their count reaches zero. Rules are schemata: one rule stands for the
reactions of every matching instance.

## A first model

Kowald's fixed-population mtDNA model: molecules degrade, replicate into the
vacancies below a copy-number target `K`, and a wild-type replication can
produce a new mutant clone, advantaged (`sm`) with probability `fracAdv`.

```json
{
  "schema": "crn-craft.rules-model/0.1",
  "name": "mtdna_fixed_population",
  "parameters": {"K": 20, "delta": 1.0, "sw": 1.0, "sm": 1.5, "mu": 0.1, "fracAdv": 0.5, "w0": 20},
  "types": [{"name": "mtDNA", "attributes": ["clone", "s", "wt"]}],
  "initial": [{"type": "mtDNA", "attributes": {"clone": 0, "s": "sw", "wt": 1}, "count": "w0"}],
  "observables": {
    "N": {"type": "mtDNA"},
    "D": {"type": "mtDNA", "weight": "x.wt + sm * (1 - x.wt)"},
    "mutants": {"type": "mtDNA", "where": "x.wt == 0"},
    "mutant_clones": {"type": "mtDNA", "where": "x.wt == 0", "reduce": "distinct"}
  },
  "rules": [
    {"name": "degradation", "reactants": [{"label": "x", "type": "mtDNA"}], "rate": "delta", "products": []},
    {"name": "replication",
     "reactants": [{"label": "x", "type": "mtDNA", "factor": "x.s * (1 - mu * x.wt)"}],
     "rate": "max(0, K - N) / D",
     "products": [{"type": "mtDNA", "copy": "x"}, {"type": "mtDNA", "copy": "x"}]},
    {"name": "mutation",
     "reactants": [{"label": "x", "type": "mtDNA", "where": "x.wt == 1", "factor": "x.s * mu"}],
     "rate": "max(0, K - N) / D",
     "products": [{"type": "mtDNA", "copy": "x"},
                  {"type": "mtDNA", "attributes": {"clone": "new_id()", "s": "choice(fracAdv, sm, sw)", "wt": 0}}]}
  ]
}
```

The full file is `models/rules/mtdna-fixed.rules.json`; the aggregation,
immigration–death and immigration–dimerization models sit beside it.

```sh
crn-craft rules-check models/rules/mtdna-fixed.rules.json
crn-craft rules MODEL REQUEST [--cancel-file PATH]
crn-craft describe rules-model     # JSON Schema; also rules-check, rules-run, rules-result
```

## Model

| Field | Meaning |
|---|---|
| `parameters` | Finite numbers by name. Parameters and observables share one namespace. |
| `types` | `{name, attributes}`. Attribute values are finite f64; `-0` is stored as `0`. |
| `initial` | `{type, attributes, count}`: every attribute, and a nonnegative integer count, as numbers or parameter expressions. Equal attribute vectors merge. |
| `observables` | Named reductions over the instances of one type (below). |
| `rules` | `{name, reactants, rate, products}` with zero, one or two reactant patterns. |
| `stop` | Optional condition over parameters and observables. Once true, the state is held to the end of the report grid and the replicate is `stopped`. |

Names are ASCII identifiers (`[A-Za-z_][A-Za-z0-9_]*`, at most 64 bytes) and
cannot be function names. Unknown fields and duplicate keys are rejected.

### Observables

`{type, label?, where?, weight?, reduce?}`. The label defaults to `x`, the
weight to 1. Over the present instances of `type` that satisfy `where`:

| `reduce` | Value |
|---|---|
| `sum` (default) | Σ weight × count |
| `distinct` | Σ weight, each instance once (weight 1 counts distinct instances, e.g. clones) |
| `max`, `min` | largest or smallest weight; 0 when no instance matches |

Weights must be finite; they may be negative.

### Rules and propensities

A reactant pattern is `{label, type, where?, factor?}`. `factor` is a local
factor h(x) ≥ 0 (default 1); an instance that fails `where` has h = 0. The
rule's `rate` is a global factor g ≥ 0 over parameters and observables. With
n_x the count of instance x:

| Reactants | Propensity |
|---|---|
| none | g |
| one pattern | g × Σ_x h(x) n_x |
| two identical patterns (same type, `where` and `factor`) | g × Σ over unordered pairs of distinct molecules h(x) h(y) = g × ½[(Σ h n)² − Σ h² n] |
| two patterns over different types | g × (Σ h_x n_x)(Σ h_y n_y) |

For one instance this is the stochastic mass-action convention n(n−1)/2 of the
network engine (a `2A` channel with coefficient c). Two *different* patterns
over the same type are rejected in 0.1. Patterns are compared after their
labels are resolved, so `a.size` and `b.size` are identical factors.

`g` is evaluated only when the reactant weight is positive, so a rate such as
`(K - N) / D` is never evaluated after extinction. A negative or non-finite
rate or factor fails the replicate (`CRN-R030`).

Reactants are consumed. A product is `{type, copy?, attributes}`: `copy` names
a reactant of the same type whose attributes are copied, and `attributes`
overrides some (with `copy`) or gives all (without). A product equal to a
reactant cancels its consumption, so `x → x + x` only adds one copy.

### Expressions

Operators, loosest first: `||`, `&&`, `== !=`, `< <= > >=`, `+ -`, `* /`,
unary `- !`, `^` (right associative; `-x^2` is `-(x^2)`). Comparisons and
logic return 1 or 0; nonzero is true; `&&`, `||` and `if` short-circuit.
Functions: `min max exp log sqrt abs floor ceil if(c, a, b)`. Numbers are
decimal with an optional exponent.

| Field | Names it may use |
|---|---|
| `initial` attributes and counts | parameters |
| pattern `where`, `factor`; observable `where`, `weight` | parameters, `label.attribute` of that instance |
| rule `rate`, `stop` | parameters, observables |
| product `attributes` | parameters, `label.attribute` of the rule's reactants, random draws |

Random draws exist only in product attributes: `bernoulli(p)`,
`choice(p, a, b)` (a with probability p), `uniform(a, b)`, `normal(m, s)`
(Box–Muller) and `new_id()`, which returns 1 + the largest initial attribute
value (floored) on its first call and counts up by one, so new clones never
reuse an identity of the initial state or of each other.

## Run request and result

`crn-craft.rules-run/0.1` binds the model by relative path and SHA-256 and
gives `time.end`, strictly increasing `time.report_times` in [0, end],
`ensemble.replicates` and `root_seed`, `outputs` (`statistics` from mean and
variance, `keep_paths`, `keep_instances`) and an explicit `budget`
(`max_events_per_replicate`, `max_sampling_steps_per_replicate` counting every
random word, `max_instances` alive at once, `threads` 1–64,
`max_output_bytes`, optional `wall_time_s`).

Report values are right-continuous: an event at a report time is included.
`crn-craft.rules-result/0.1` holds the type, rule and observable names, one
record per replicate (termination, events, stop time, events per rule, random
words, rejected pair draws, peak and created instances, diagnostic), the
observable series of the kept paths, final instances of the kept replicates
(sorted by type and attributes) and, when every replicate completed or
stopped, the mean, sample variance and standard error per report time and
observable. Terminations: `completed`, `stopped`, `event_limit`,
`sampling_limit`, `instance_limit`, `wall_time_limit`, `cancelled`, `failed`,
`output_limit`; the ensemble reports `completed` when every replicate completed
or stopped.

## Numerical policy

- Gillespie's direct method. Rule selection scans the rule propensities in
  order; instance selection descends a pairwise sum tree per distinct pattern.
  Tree totals are recomputed from their children on each update (no
  accumulated drift); their rounding depends on the slot layout, which is a
  deterministic function of the replicate's history.
- Identical-pattern pairs draw two molecules from the tree and reject the same
  molecule (probability 1/n when both draws hit one instance). This is exact;
  it is slow when one instance dominates the weight (multiplicative kernels
  after gelation) and is counted in `rejected_pairs`.
- A positive rule propensity or instance weight below 2⁻⁵³ of its total cannot
  be selected exactly by a 53-bit uniform and fails with `CRN-R031`.
- Each replicate draws from its own ChaCha8 stream keyed by SHA-256 of the
  domain `crn-craft.rules/0.1` and the root seed, with the replicate index as
  stream id. Results do not depend on the thread count.
- Counts and identities stay below 2⁵³.

## Diagnostics

| Code | Meaning |
|---|---|
| CRN-E001–E004 | malformed JSON, duplicate keys, unknown fields or wrong types (`E003`), size limits |
| CRN-R010 | invalid, reserved or duplicate name; observable shadowing a parameter |
| CRN-R011 | unknown type |
| CRN-R012 | expression syntax, unknown name, out-of-scope name or random draw outside products |
| CRN-R013 | unsupported rule shape, missing or unknown attribute, bad `copy` |
| CRN-R014 | invalid initial count or attribute |
| CRN-R015 | model size limit |
| CRN-R020 | invalid run request (hash, path, grid, ensemble, budget) |
| CRN-R030 | negative or non-finite rate, factor, weight or product attribute; overflow |
| CRN-R031 | weight below 53-bit resolution |
| CRN-R040–R045 | event, sampling, instance, wall-time, cancellation and output limits |
| CRN-R050 | warning: incomplete ensemble, statistics withheld |

## Validation (K1)

`validation/check_rules.py` compares 100,000 replicates per fixture with exact
laws (Poisson; master equations by uniformization over truncated or finite
state spaces): immigration–death, immigration–dimerization (rules and network
engines), Kowald's fixed-population mtDNA model, and aggregation with constant
and multiplicative kernels over the 77 partitions of 12. Thirty marginal CDFs
share one DKW family at alpha 0.001. Results: `validation/rules/RESULTS.md`.

## Limits of 0.1

- No time-dependent rates or ODE coupling (Kowald's variable-population model
  with ATP needs these; planned as K2), no compartments.
- Two different patterns over one type, three or more reactants, and counts
  inside expressions (largest clone *size*) are not supported yet.
- No Python builder; models are JSON.
- Performance is not yet benchmarked; propensities of all rules are
  recomputed after every event, which suits models with few rules.
