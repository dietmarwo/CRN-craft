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
the language, starting from the `tree` method; [mtDNA ensembles](mtdna-rules.md)
covers the other methods, trackers and the older request versions. Protocol and
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

```sh
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
"continuous": [
  {"name": "ATP", "initial": "ATP0",
   "derivative": "f * wild_type - v1 * ATP - v2 * N",
   "lower_bound": 0, "on_bound": "stop"}
]
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
"continuous": [
  {"name": "m", "initial": "m0", "derivative": "g * m", "independent": true},
  {"name": "x", "initial": 0, "derivative": "m * A - x"}
]
```

Without `independent`, a rules-model/0.4 document means exactly what the
same document means as 0.3. Stored 0.1–0.3 documents are unchanged.

Bounds keep their meaning. `clamp` holds an independent variable at its
bound as above. `stop` ends the replicate when the variable reaches its
bound with a derivative that is not positive; because the variable does not
depend on the replicate, that happens at the same time in every replicate
that has not stopped, failed or reached the horizon before.

### The `hybrid` method (rules-run/0.5)

Only `hybrid` runs a model with continuous variables; every other method
refuses one with `CRN-R020` at `/method`. The request gives its tolerances:

```json
"method": {"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}
```

with `rtol` in [1e-12, 1e-2] and `atol` in (0, 1e-2] (`CRN-R020` at
`/method/rtol` or `/method/atol`); `hybrid` needs rules-run/0.5 (`CRN-R020`
at `/schema`). The Python helper writes these values for `method="hybrid"`.

The model is a piecewise-deterministic Markov process (Davis 1984, DOI
[10.1111/j.2517-6161.1984.tb01308.x](https://doi.org/10.1111/j.2517-6161.1984.tb01308.x)),
simulated by integrated propensities (Alfonsi et al. 2005, DOI
[10.1051/proc:2005001](https://doi.org/10.1051/proc:2005001)):

1. After each event draw a unit exponential threshold E from one open
   uniform, E = −ln u.
2. Integrate the variables and the integrated total propensity
   Λ(t) = ∫ a₀ ds, with Λ = 0 at the event, by the Dormand–Prince 5(4)
   pair (DOI [10.1016/0771-050X(80)90013-3](https://doi.org/10.1016/0771-050X(80)90013-3)).
   Steps follow the standard error control of Hairer, Nørsett and Wanner
   (*Solving Ordinary Differential Equations I*, 2nd ed., 1993, II.4): the
   root-mean-square of the error estimate over all components, each scaled
   by `atol + rtol·max(|y₀|, |y₁|)`, must not exceed 1; the next step is
   h·min(10, max(0.2, 0.9·err^(−1/5))). The first step uses their starting
   rule, later ones the previous prediction, and each step is at most four
   times (E − Λ)/a₀ at its start (but not below 16 rounding units of the
   time).
3. The event happens where Λ reaches E, located by the Illinois method on
   the fourth-order dense output (II.6, after Shampine 1986, DOI
   [10.1090/S0025-5718-1986-0815836-3](https://doi.org/10.1090/S0025-5718-1986-0815836-3)).
   It selects channel j with probability a_j/a₀ evaluated at that point,
   from one uniform, with the resolution rule of Static Direct
   (`CRN-R031`).
4. Report values and bound crossings come from the same dense output; a
   bound crossing that comes first ends the step there.

Rates that do not read a variable are computed once per interval.

Independent variables (rules-model/0.4) are integrated once per prepared
run, by the same Dormand–Prince pair, step control and bound handling,
over [0, end] at the request's `rtol` and `atol`, with their own error
norm. The accepted steps and their fourth-order dense output are kept, and
every replicate and worker reads the variables from them, so results stay
independent of `threads`. Per replicate, the method then integrates only
the dependent variables and Λ, reading the independent variables from the
dense output at each stage time; the error norm covers those components
only. A stop bound of an independent variable ends the shared solution, and
every replicate that reaches its time stops there (a tie with the
threshold goes to the stop, as above). If the shared solve fails (a
non-finite derivative, a step that underflows, or more attempted steps than
`max_sampling_steps_per_replicate`), it ends at the start of the failing
step, and every replicate that reaches that time fails with the solve's
diagnostic (`CRN-R030`, or `CRN-R041` for the step budget). The solve runs
in `prepare` (and again for every `bind_parameters` candidate); the step
budget bounds it.

The method lowers the model like the static methods, so it refuses models with
dynamic instances; it refuses ATP models (use `thinning`); and it also runs
models without continuous variables, as an exact jump process up to
rounding. It is approximate at its tolerances. Its sampler is
`chacha8-u53/rules-approximate-hybrid-dopri54/1`; rules-result/0.3 records
the request's `method`, tolerances included. Integration steps, accepted or
rejected, count against `max_sampling_steps_per_replicate` like leap
iterations (`CRN-R041`), and a step below four rounding units of the time
fails the replicate (`CRN-R030`). Validation:
[`validation/hybrid/RESULTS.md`](../validation/hybrid/RESULTS.md).

## Run request and result

`crn-craft.run/1.0` has the shape of `rules-run/0.6` and adds three things:
`method.selector` for `static_direct` (`linear`, `composition_rejection`,
`sorting`, `rssa_cr` or `auto`, the default; see [exact methods](direct-ssa.md)),
`outputs.quantiles` (up to 16 probabilities in (0, 1); the result's
`statistics.quantiles` holds the quantile of every report time and observable,
and the median is the quantile at 0.5; see the [result contract](result-contract.md))
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

Report values are right-continuous: an event at a report time is included.
`crn-craft.rules-result/0.1` holds the type, rule and observable names, one
record per replicate (termination, events, stop time, events per rule, random
words, peak instances (including the one that exceeded `max_instances`) and
created instances, diagnostic), the
observable series of the kept paths, final instances of the kept replicates
(sorted by type and attributes) and, when every replicate completed or
stopped, the mean, sample variance and standard error per report time and
observable. Statistics are accumulated in replicate order (Welford) while the
ensemble runs; only kept series stay in memory. A path's random words include
product draws, so the event that hits the sampling budget can take
`rng_words` slightly past it. Terminations: `completed`, `stopped`, `event_limit`,
`sampling_limit`, `instance_limit`, `wall_time_limit`, `cancelled`, `failed`,
`output_limit`; the ensemble reports `completed` when every replicate completed
or stopped.

## Numerical policy

- Gillespie's direct method. Rule selection scans the rule propensities in
  order; instance selection descends a pairwise sum tree per distinct pattern.
  Tree totals are recomputed from their children on each update (no
  accumulated drift); their rounding depends on the slot layout, which is a
  deterministic function of the replicate's history.
- Identical-pattern pairs use a pair tree: each node stores its pattern
  weight S and the weight P of pairs of distinct molecules inside it,
  P = P_left + P_right + S_left·S_right (a leaf holds h²·n(n−1)/2). The
  propensity is g × P at the root, a sum of nonnegative terms without the
  cancellation of ½[(Σ h n)² − Σ h² n]. Selection descends by P into a child or
  stops at the cross term, then picks one molecule on each side by S; no
  rejection (sampler `chacha8-u53/rules-direct/2`; version 1, used only for the
  first K1 gate run, drew two molecules and rejected the same one).
- A positive rule propensity or instance weight below 2⁻⁵³ of its total cannot
  be selected exactly by a 53-bit uniform and fails with `CRN-R031`. For pair
  rules the check is conservative: the smallest factor squared must reach
  2⁻⁵³ of the pair total.
- Overflowing weights, pair totals and observable totals fail the replicate
  (`CRN-R030`) rather than becoming 0 or `null`.
- Each replicate draws from its own ChaCha8 stream keyed by SHA-256 of the
  domain `crn-craft.rules/0.1` and the root seed, with the replicate index as
  stream id. Results do not depend on the thread count.
- Counts stay at most 2⁵³ (`CRN-R042` otherwise) and identities below 2⁵³.

## Diagnostics

| Code | Meaning |
|---|---|
| CRN-E001–E004 | malformed JSON, duplicate keys, unknown fields or wrong types (`E003`), size limits |
| CRN-R010 | invalid, reserved or duplicate name; observable shadowing a parameter |
| CRN-R011 | unknown type |
| CRN-R012 | expression syntax, unknown name, out-of-scope name or random draw outside products; an independent variable's derivative that reads an observable or a dependent variable |
| CRN-R013 | unsupported rule shape, missing or unknown attribute, bad `copy`; both or neither of `rate` and `propensity`, `where` or `factor` on a propensity rule's reactant; `lower_bound` without `on_bound` or the reverse |
| CRN-R014 | invalid initial count or attribute; a continuous variable's initial value or bound that is not finite, or an initial value below its bound |
| CRN-R013 and R020 | for `model/1.0`: a pattern `count` outside 1 to 2,147,483,647, or on an attributed type or with `where` or `factor` (`R013`); a higher-order rule under a method that refuses it (`R020` at `/method`) |
| CRN-R015 | model size limit (the table above; also 64 continuous variables); a field used before its model version (`atp` in 0.1, `propensity` and `continuous` in 0.1 or 0.2, `independent` before 0.4); `atp` together with `continuous` |
| CRN-R016, R017 | a network model the converter cannot write (`R016`); a model without a flat view, which finite state projection, analysis, ODE export and GillesPy2 export need (`R017`) |
| CRN-R020 | invalid run request (hash, path, grid, ensemble, budget, method, selector, quantiles, tracker, log or version gate; generated-kernel fingerprint mismatch) |
| CRN-R030 | negative or non-finite rate, factor, weight or product attribute; overflow; under `run/1.0`, 1,000,000 consecutive waits that do not advance the clock |
| CRN-R031 | weight below 53-bit resolution (Static Direct: a channel below 2⁻⁵³ of the total) |
| CRN-R040–R045 | event, sampling, instance (including a count above 2⁵³), wall-time, cancellation and output limits |
| CRN-R050 | warning: incomplete ensemble, statistics withheld |
| CRN-R051 | warning: a requested specialization did not apply (`mtdna_direct` could not certify the model and ran Static Direct, or an explicit `sparse_trackers: true` found no eligible tracker; the rules-run/0.6 default does not warn); the message gives the reason |

## Validation (K1)

`validation/check_rules.py` compares 100,000 replicates per fixture with exact
laws (Poisson; master equations by uniformization over truncated or finite
state spaces): immigration–death, immigration–dimerization, the fixed-population mtDNA model, and aggregation with constant
and multiplicative kernels over the 77 partitions of 12. Thirty marginal CDFs
share one DKW family at alpha 0.001. Results: `validation/rules/RESULTS.md`.

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
  ([mtDNA ensembles](mtdna-rules.md)); rules-model/0.3 adds general
  continuous variables, simulated approximately by `hybrid` (below), and
  rules-model/0.4 independent variables, which `hybrid` integrates once per
  run. An independent variable can carry explicit time dependence (a clock
  `t' = 1`, or a forcing such as a seasonal oscillator). No compartments.
- Two different patterns over one type, higher-order rules over attributed
  types, and counts inside expressions (largest clone *size*) are not
  supported yet.
- Models are JSON. `crn_craft` prepares and runs ensembles from Python
  ([mtDNA ensembles](mtdna-rules.md)), and `crn_craft.Model` builds
  `model/1.0` documents over attribute-free species from
  `add_reaction(rate, reactants, products)` calls: a numeric rate is a
  mass-action constant with the counting above, a string the whole
  propensity.
- The tree method recomputes the propensities of all rules after every
  event, which suits models with few rules. Closed models can use the
  static methods; benchmarks are in `validation/mtdna/`.

## Prepared Python ensembles and static methods

Python prepares these runs through the top-level `crn_craft`, which writes
`crn-craft.model/1.0` and `crn-craft.run/1.0`, or through `crn_craft.rules`,
the API of release 0.6.0 under its old name, which writes the `rules-*`
versions unless a call needs a 1.0 feature. Both have event observation, Static Direct and next reaction,
approximate leaps, and ATP thinning. See the
[mtDNA rules guide](mtdna-rules.md) for rules-run/0.2–0.4, rules-result/0.2 and rules-model/0.2, examples,
supported model shapes and validation commands. rules-run/0.5 adds `hybrid`
and returns rules-result/0.3 (0.2 plus `method` and `continuous`); see
[continuous variables](#continuous-variables-rules-model03). rules-run/0.6
adds no field: `static_direct` defaults to Ziggurat waits and sparse
integration of eligible count trackers, `waiting_time: "inverse"` with
`sparse_trackers: false` restores the 0.5 default, and the result is
rules-result/0.3
([Static Direct options](mtdna-rules.md#generic-static-direct-options)).
