# Model reference: `crn-craft.model/1.0`

> This file keeps its former name, "rules language", so that old links
> resolve. It is the reference for the model contract `crn-craft.model/1.0`,
> and for the older `rules-model` versions that mean the same thing.

**Contracts.** The model is `crn-craft.model/1.0` (`crn-craft describe
model-1.0`). The versions `crn-craft.rules-model/0.1` to `0.4` are still read
and run; a `rules-model/0.4` document means what a `model/1.0` document means,
with smaller limits ([limits](#limits)). The request is `crn-craft.run/1.0`
or `rules-run/0.1` to `0.6`, and the result `crn-craft.result/1.0` or
`rules-result/0.1` to `0.3` ([contracts](contracts.md)). The crate is
`crn-craft` (module `crn_craft::rules`), the command line `crn-craft check`
and `crn-craft run`, and the Python package `crn_craft`. This page describes
the language, starting from the `tree` method; [methods reference](methods.md)
covers the other methods. Protocol and
gates: `docs/decisions/rules-k0-k1-2026-11-02.md`.

Species are not a fixed list. Some models cannot list their species in
advance: each new mitochondrial mutant is its own clone, polymers have a
length, clusters have a size, cells carry a lineage. The engine simulates
such systems exactly. Species *types* carry numeric *attributes*; an
*instance* is a type plus one value per attribute and a copy number.
Instances appear when a rule produces a new attribute vector and disappear
when their count reaches zero. Rules are schemata: one rule stands for the
reactions of every matching instance. A model whose types have no attributes
is an ordinary reaction network, written with `crn_craft.Model`, and a network
`model/0.1` file converts to one with `crn-craft convert`.

## A first model

the fixed-population mtDNA model: molecules degrade, replicate into the
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

```text
crn-craft check models/rules/mtdna-fixed.rules.json
crn-craft run MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]
crn-craft describe model-1.0       # JSON Schema; also run-1.0, result-1.0, check-1.0
```

The file above is a `rules-model/0.1` document. Writing its `schema` as
`crn-craft.model/1.0` gives the same model. `rules-check` and `rules` are
aliases of `check` and `run`.

## Model

| Field | Meaning |
|---|---|
| `parameters` | Finite numbers by name. Parameters and observables share one namespace. |
| `types` | `{name, attributes}`. Attribute values are finite f64; `-0` is stored as `0`. |
| `initial` | `{type, attributes, count}`: every attribute, and a nonnegative integer count, as numbers or parameter expressions. Equal attribute vectors merge. |
| `observables` | Named reductions over the instances of one type (below). |
| `rules` | `{name, reactants, rate, products}` with reactant patterns; in 0.3 `propensity` may replace `rate`. Up to `rules-model/0.4`, at most two patterns. In `model/1.0`, any number, and a pattern may carry a `count` ([higher-order rules](#higher-order-rules-model10)). |
| `stop` | Optional condition over parameters and observables. Once true, the state is held to the end of the report grid and the replicate is `stopped`. |
| `continuous` | rules-model/0.3: continuous variables, ODEs between events that rates may read (below); rules-model/0.4 adds `independent` variables. |

Names, including the model's `name`, are ASCII identifiers
(`[A-Za-z_][A-Za-z0-9_]*`, at most 64 bytes) and cannot be function names. Unknown fields and duplicate keys are rejected.

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

For one instance this is the stochastic mass-action convention n(n−1)/2 (a
`2A` channel with coefficient c). Two *different* patterns over the same type
are rejected in 0.1. Patterns are compared after their
labels are resolved, so `a.size` and `b.size` are identical factors.

#### Counting identical reactants, and rebop

A `rate` counts *distinct combinations of reactant molecules*, Gillespie's
h_μ (Gillespie 1977, DOI
[10.1021/j100540a008](https://doi.org/10.1021/j100540a008)): with rate k,
`A` fires at k·n, `A + B` at k·n_A·n_B and `A + A` at k·n(n−1)/2. rebop 0.10
uses the same expression for distinct species but leaves out the 1/m! for m
identical reactants: its numeric rate k on `["A", "A"]` fires at k·n(n−1)
(its documentation equates `4.2` with `"4.2 * A * (A-1)"`). A rebop model
with a numeric rate k on `A + A` corresponds to rate 2k here; models without
identical reactants agree. String rates are full propensities in both. The
Python builder `crn_craft.Model` follows the convention of this
engine, so its JSON equals a hand-written model; see
`validation/check_builder.py` for the statistical comparison with rebop.

#### Whole-rule propensities (rules-model/0.3)

A rule may give `propensity`, its total rate P as an expression over
parameters and observables, instead of `rate`. Exactly one of the two is
present (`CRN-R013`); `propensity` in 0.1 or 0.2 is `CRN-R015`, and a 0.1 or
0.2 rule without `rate` stays `CRN-E003`.

| Reactants | Propensity |
|---|---|
| none | P |
| one pattern of type T | P if T has a molecule, else 0 |
| two patterns of type T | P if T has two molecules, else 0 |
| patterns of types T and U | P if T and U each have a molecule, else 0 |

The reactant patterns take no `where` or `factor` (`CRN-R013`), and P cannot
read their labels (`CRN-R012`): the expression already says how the rate
depends on the state, for example `"vmax * S / (Km + S)"` over an observable
`S`. P is evaluated only while it applies; a negative or non-finite value
fails the replicate (`CRN-R030` at `/rules/i/propensity`). A firing picks
reactant molecules uniformly: one of the type's molecules, or an unordered
pair of distinct molecules, whatever their attributes. The `tree` method
samples this for any type. The static methods lower each pattern to its
tuples, so they accept a propensity rule only when each reactant type has a
single tuple, as for attribute-free species (`CRN-R020` at `/method`
otherwise); `mtdna_direct` falls back to Static Direct with `CRN-R051`, and
generated rates kernels refuse such models (`CRN-R020`), while generated
event loops emit them with the interpreter's arithmetic. For one
reactant, `propensity: "k * A"` and `rate: "k"` give the same floating-point
value, so they reproduce each other's trajectories exactly.

`g` is evaluated only when the reactant weight is positive, so a rate such as
`(K - N) / D` is never evaluated after extinction. A negative or non-finite
rate or factor fails the replicate (`CRN-R030`).

Reactants are consumed. A product is `{type, copy?, attributes}`: `copy` names
a reactant of the same type whose attributes are copied, and `attributes`
overrides some (with `copy`) or gives all (without). A product equal to a
reactant cancels its consumption, so `x → x + x` only adds one copy.

### Higher-order rules (`model/1.0`)

A reactant pattern of an attribute-free type may carry `count`, an integer
from 1 to 2,147,483,647 (default 1). It stands for `count` molecules of that
type. A rule has any number of patterns, so `A + A + A → ...` is a pattern with
`count: 3`, or three identical patterns, which means the same. Patterns with
`count` take no `where` or `factor`, and their type has no attributes
(`CRN-R013`). The rule's `rate` gives

rate × Π C(n, k)

over the patterns, with n the count of the type's molecules and k the pattern's
`count`; this is the combinatorial convention of every other rate. For a
whole-rule `propensity`, the value applies while every pattern has enough
molecules. A product list may be any length up to 64 molecules.

Which methods run higher-order rules:

| Method | Higher-order rules |
| --- | --- |
| `static_direct`, `next_reaction`, `tau`, `hybrid`, generated event loops | yes |
| `tree` (the default when `method` is absent), `mtdna_direct`, `binomial`, generated rates kernels | refused: `CRN-R020` at `/method` |

`crn_craft.Model.add_reaction` takes three or more reactant molecules and
writes the `count` for you. `crn-craft convert` writes a network reaction of
three or more molecules the same way, with labels `x`, `y` or `x1`, `x2`, ….

### Expressions

Operators, loosest first: `||`, `&&`, `== !=`, `< <= > >=`, `+ -`, `* /`,
unary `- !`, `^` (right associative; `-x^2` is `-(x^2)`). Comparisons and
logic return 1 or 0; nonzero is true; `&&`, `||`, `if` and `choice`
short-circuit. An undefined (NaN) operand of `!`, `&&`, `||`, `if`, `min` or
`max` gives NaN, and a NaN `where` or `stop` fails the replicate (`CRN-R030`)
instead of counting as true; comparisons with NaN are false (IEEE). Functions:
`min max exp log sqrt abs floor ceil if(c, a, b)`. Numbers are decimal with an
optional exponent; only ASCII whitespace separates tokens.

| Field | Names it may use |
|---|---|
| `initial` attributes and counts | parameters |
| pattern `where`, `factor`; observable `where`, `weight` | parameters, `label.attribute` of that instance |
| rule `rate`, `stop` | parameters, observables |
| product `attributes` | parameters, `label.attribute` of the rule's reactants, random draws |

Random draws exist only in product attributes: `bernoulli(p)`,
`choice(p, a, b)` (a with probability p; only the chosen branch is evaluated),
`uniform(a, b)`, `normal(m, s)` (Box–Muller) and `new_id()`. `new_id()`
returns 1 + max(0, largest initial attribute value), floored, on its first call
and counts up by one, so new clones never reuse an identity of the initial
state or of each other. Values written by constants or parameters in products
are not tracked: keep such identities below the initial maximum. Identities
must stay below 2⁵³: a model that calls `new_id()` with an initial attribute of
2⁵³ − 1 or more is rejected (`CRN-R014`), and reaching 2⁵³ during a run fails
the replicate (`CRN-R030`).

## Continuous variables (rules-model/0.3)

A model may declare continuous variables in `continuous`, each
`{name, initial, derivative, lower_bound?, on_bound?}`:

```json
{
  "continuous": [
    {"name": "ATP", "initial": "ATP0",
     "derivative": "f * wild_type - v1 * ATP - v2 * N",
     "lower_bound": 0, "on_bound": "stop"}
  ]
}
```

| Field | Meaning |
|---|---|
| `name` | shares one namespace with parameters and observables (`CRN-R010`); also the result column |
| `initial` | number or parameter expression, finite and not below the bound (`CRN-R014`) |
| `derivative` | dx/dt between events, over parameters, observables and continuous variables; a non-finite value fails the replicate (`CRN-R030` at `/continuous/i/derivative`) |
| `lower_bound`, `on_bound` | together or not at all (`CRN-R013`): a number or parameter expression, and `stop` or `clamp` |

Rule rates and propensities may read the variables, so rates change between
events. The model's `stop`, observable filters and weights, reactant
patterns and products cannot read them (`CRN-R012`; a stop is checked at
events only), and neither can run-request trackers (`CRN-R020`). Results
list the variables as observable columns after the model's observables, in
declaration order, and rules-result/0.3 names them in `continuous`.

At a lower bound, `stop` ends the replicate as `stopped` when the variable
reaches the bound while its derivative there is not positive (also at time 0
and after an event); the variable reads the bound on every later report and
the other observables hold their values. This is the rule of the legacy ATP
block at zero. `clamp` holds the variable at the bound while its derivative
there is negative, and it follows the ODE again once the derivative is
positive; rates and reports never read it below the bound.

`models/rules/mtdna-three-class-atp-continuous.rules.json` is the ATP model
of `mtdna-three-class-atp.rules.json` written this way: the parameters `f`,
`v1`, `v2`, `c` and `ATP0` replace the `atp` block, and each coupled
replication rate is multiplied by `c / (c + ATP)` in its expression. A model
cannot have both `atp` and `continuous` (`CRN-R015`); the legacy block and
`thinning` keep working unchanged.

### Independent variables (rules-model/0.4)

`crn-craft.rules-model/0.4` is 0.3 plus one optional field of a continuous
variable, `"independent": true` (default `false`; `null` is rejected with
`CRN-E003`, and the field in 0.1–0.3 with `CRN-R015` at
`/continuous/i/independent`). An independent variable's ODE does not depend
on the stochastic state: its `derivative` may read only parameters and
other independent variables. A read of an observable or of a variable
without `independent: true` is refused with `CRN-R012` at
`/continuous/i/derivative`. `initial` and `lower_bound` read parameters, as
for every variable. Rule rates and propensities, and the derivatives of the
other (dependent) variables, read independent variables like any continuous
variable, and results report them like the others.

```json
{
  "continuous": [
    {"name": "m", "initial": "m0", "derivative": "g * m", "independent": true},
    {"name": "x", "initial": 0, "derivative": "m * A - x"}
  ]
}
```

Without `independent`, a rules-model/0.4 document means exactly what the
same document means as 0.3. Stored 0.1–0.3 documents are unchanged.

Bounds keep their meaning. `clamp` holds an independent variable at its
bound as above. `stop` ends the replicate when the variable reaches its
bound with a derivative that is not positive; because the variable does not
depend on the replicate, that happens at the same time in every replicate
that has not stopped, failed or reached the horizon before.

### The `hybrid` method

Only `hybrid` runs a model with continuous variables; every other method
refuses one with `CRN-R020` at `/method`. Its tolerances `rtol` and `atol`,
its integration scheme and its error control are in the
[methods reference](methods.md#the-hybrid-method).

## Run request and result

The request (`crn-craft.run/1.0` and the `rules-run` versions) is in the
[contracts](contracts.md#the-run-request), with its budgets and the admission
rule. The result is in the [results reference](results.md). The numerical
policy of the `tree` method is in the
[methods reference](methods.md#the-tree-method-and-its-numerical-policy).

## Diagnostics

Every code has an entry in the [diagnostics reference](diagnostics.md), which
also has a table of the `CRN-R` codes of models and run requests.

## Validation

`validation/check_rules.py` compares the engine with exact laws. The results
are in `validation/rules/RESULTS.md`, and the
[technical report](../evidence/technical-report.md) has the account.

## Limits

| | `rules-model/0.1` to `0.4` | `crn-craft.model/1.0` |
| --- | ---: | ---: |
| Types | 64 | 4,096 |
| Parameters | 1,024 | 16,384 |
| Observables | 256 | 4,096 |
| Rules | 1,024 | 65,536 |
| Instance tuples (static methods) | 256 | 4,096 |
| Channels (static methods) | 4,096 | 65,536 |

Both versions share the other limits: 64 attributes per type, 64 products
per rule, expressions of at most 4,096 bytes, 100,000 initial entries,
documents of at most 16 MiB and JSON nesting of at most 48
(`describe()["limits"]`). Exceeding a model limit is `CRN-R015`. To use the
larger limits, change only the `schema` string to `crn-craft.model/1.0`.

## What is not supported

- ODE coupling: rules-model/0.2 adds one fixed case, the ATP flow of
  the variable-population mtDNA model, solved with exact thinning
  ([ATP and thinning](methods.md#the-atp-block-and-thinning)); rules-model/0.3 adds general
  continuous variables, simulated approximately by [`hybrid`](methods.md#the-hybrid-method), and
  rules-model/0.4 independent variables, which `hybrid` integrates once per
  run. An independent variable can carry explicit time dependence (a clock
  `t' = 1`, or a forcing such as a seasonal oscillator). No compartments.
- Two different patterns over one type, higher-order rules over attributed
  types, and counts inside expressions (largest clone *size*) are not
  supported yet.
- Models are JSON. `crn_craft` prepares and runs ensembles from Python
  ([Python API](python-api.md)), and `crn_craft.Model` builds
  `model/1.0` documents over attribute-free species from
  `add_reaction(rate, reactants, products)` calls: a numeric rate is a
  mass-action constant with the counting above, a string the whole
  propensity.
- The tree method recomputes the propensities of all rules after every
  event, which suits models with few rules. Closed models can use the
  static methods; benchmarks are in `validation/mtdna/`.
