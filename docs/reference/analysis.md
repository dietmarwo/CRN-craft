# Network analysis

Every closed model has a flat view: its species, its channels, the change
each channel makes and the channels whose rates each change affects.
`crn-craft analyze MODEL [STATE]` and `Model.analyze()` (Python) report the
structure of that view as `crn-craft.analysis/1.0` and run no simulation.
Finite state projection, the ODE right-hand sides of `crn-diffsol` and the
GillesPy2 export read the same view. A model has a flat view when its set of
species is closed: no product creates a new attribute tuple during a run, no
product attribute is random, and there is no `atp` block or continuous
variable. Any other model is refused with `CRN-R017`.

This page replaces the former description of compiled network models, which
were removed with the network engine in release 0.7.0. Network models of
`crn-craft.model/0.1` are still read: they are converted first
(`crn-craft convert`), and then every statement below applies to the converted
model.

## Indices and labels

Species and channels are indexed from zero. The species are the labels of the
flat view: the declared name of an attribute-free type, and `T{a=1,b=0.5}`
for an attributed type. For a converted network model they are the network's
species in declaration order, and the channels are its reactions in
declaration order.

Every reaction lists sorted sparse `reactants`, `products` and `net_changes`,
each term `{species, coefficient}`. Net coefficients are products minus
reactants with zeros omitted. They are the columns of the stoichiometric
matrix N (species by reactions), not rows and not propensity multipliers.
`reads` lists the species that the channel's rate reads. `dependencies[j]`
is the sorted, unique list of the channels whose rates can change when channel
j fires: the channels that read a species j changes. It describes rate
invalidation only; a scheduler must also redraw the channel that fired.

## The report

`crn-craft.analysis/1.0` has `success`, `simulation_available` (always
`false`), `model_sha256`, an optional `state_sha256`, `analysis` and
`diagnostics`. An input or numerical failure has `success: false` and
`analysis: null`. A successful analysis holds:

| Field | Content |
| --- | --- |
| `species`, `reactions`, `dependencies` | the structure above |
| `initial_counts`, `evaluation_counts`, `propensities` | the propensity of every channel, evaluated at the evaluation state: the initial counts, or the counts of a supplied `crn-craft.state/0.1` document, which names species by label |
| `conservation` | the linear conservation laws (below) |
| `nonnegative_conservation` | the nonnegative conservation cone: extreme rays, work and limits |
| `linkage` | complexes, linkage classes, deficiency and weak reversibility |
| `initial_class_upper_bounds` | the best available bound of each species in the initial conservation class |

`model_sha256` and `state_sha256` hash the raw supplied bytes. The evaluation
state need not lie in the initial conservation class, and the report makes no
reachability claim.

## Linear conservation laws

The analyzer computes the rational left nullspace of N by checked integer
row elimination with gcd reduction, then checked rational back substitution.
Each primitive integer basis row l satisfies lᵀN = 0 and is verified against
the sparse columns. A limited sequence of elementary combinations exposes
nonnegative rows where available. The rows span the rational conservation
space, not necessarily the integer lattice of laws. Coefficients and the
initial total lᵀx₀ are decimal strings, so 128-bit integers survive JSON
readers that use doubles.

| Status | Meaning |
| --- | --- |
| `exact` | a complete rational basis, with exact residuals and totals in checked i128 arithmetic |
| `size_limit` | more than 64 species or 256 channels: no basis or rank is computed (`CRN-W010`) |
| `arithmetic_overflow` | a checked i128 operation overflowed: nothing partial is kept (`CRN-W011`) |

A size or overflow status leaves the structure and the propensities
available, and `success` stays true. It does not mean that no law exists.
A model without reactions has the identity basis, and a full-rank network
an empty one.

For a nonnegative row, each positive coefficient gives
`floor(initial_total / coefficient)` as an upper bound for that species in the
initial class. A signed row gives no count bound. A `null` bound means
"not proved", not "unbounded". The enzyme example has the laws E + ES and
S + ES + P, which bound all four species.

## The nonnegative cone

The cone of nonnegative conservation laws is computed separately, by exact
extreme-ray enumeration with limits of 64 species, 256 channels, 4,096
intermediate rays and 2,000,000 work units. Its status, limits and work count
are in the report. When the cone is exact, `initial_class_upper_bounds`
uses its bounds, and otherwise the sound bounds of the basis.
`conservation.initial_class_upper_bounds` keeps the basis-only meaning.
Exhausting a limit is `CRN-W012` to `CRN-W014`, and a failed internal check
`CRN-W015`. The linkage analysis adds the deficiency of the network; see
[the analysis notes](../history/m1b-analysis.md).

## Moved

| Earlier | Now |
| --- | --- |
| `CompiledNetwork::compile`, `analyze_bytes` on a network model | `crn_craft::analyze_bytes` and `analyze_files` read every model version; the analysis functions take a `&crn_craft::rules::FlatView` |
| `crn-craft analyze-v1`, `analysis/0.1` and `0.2` | `crn-craft analyze`, `analysis/1.0` |
| `topology_sha256`, `mass_action_coefficient`, `propensity_in_scheduler_state`, generated kernels | [parameter sweeps](../reusable-parameters.md); `crn_craft.compile_loop` |
| the validation of the conservation analysis | `validation/check_structure.py` and `validation/merge/check_analysis_references.py` |
