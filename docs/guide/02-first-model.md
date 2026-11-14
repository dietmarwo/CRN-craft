# 2. A first model

A model has species, parameters and reactions. This chapter builds one in
Python, runs it once, and shows what happens when a model has a mistake.

## Species, parameters and reactions

The builder `crn.Model` writes a model call by call. In this one, molecules
of `A` arrive at rate 4, and each molecule decays at rate 0.2. It is called
immigration–death, and it is the running example of this guide:

```python
import numpy as np
import crn_craft as crn

model = crn.Model("immigration_death")
model.add_parameter({"birth": 4.0, "death": 0.2})
model.add_species("A", 0)
model.add_reaction("birth", [], ["A"], name="birth")
model.add_reaction("death * A", ["A"], [], name="death")
print(model)
print(model.reaction_names)
```

```text
Model('immigration_death': 1 species, 2 reactions)
('birth', 'death')
```

`add_species(name, initial)` declares a species and its count at time 0.
`add_reaction(rate, reactants, products)` takes the species of each side as
lists. An empty list means that nothing is consumed or produced. The blocks
of this chapter continue from this one.

## A first run

`crn.ensemble` prepares a run and `simulate()` executes it. Here 400
independent trajectories run to time 30, and the result holds their mean at
three report times. `method="static_direct"` names the exact method for a
reaction network, Gillespie's direct method. Pass it in every call, as this
guide does; [chapter 3](03-running.md#many-trajectories-in-one-call) says
why:

```python
result = crn.ensemble(
    model, end=30.0, report_times=[0.0, 10.0, 30.0], replicates=400, seed=42,
    method="static_direct",
).simulate()
print(result.statistics["mean"])
print("exact:", (20 * (1 - np.exp(-0.2 * np.array([0.0, 10.0, 30.0])))).round(2))
```

```text
            A
time         
0.0    0.0000
10.0  17.2450
30.0  20.0925
exact: [ 0.   17.29 19.95]
```

The count of `A` has a known law: it is Poisson with mean
`20 (1 - exp(-0.2 t))`. The simulated means agree with it within their
sampling error, which chapter 4 shows how to read. An ensemble is one call
for any number of trajectories ([chapter 3](03-running.md)), and threads and
compiled models make it fast ([chapter 5](05-speed.md)).

## Rates as numbers and as strings

The first argument of `add_reaction` is the rate, and it has two forms.

- **A number is a mass-action constant.** The reaction fires at the constant
  times the number of distinct combinations of reactant molecules: `A` at
  `k * A`, `A + B` at `k * A * B`, and `A + A` at `k * A * (A - 1) / 2`.
  Write a repeated species once per molecule: `["A", "A"]` is `2 A`.
- **A string is the whole propensity.** It is an expression over parameters
  and species counts, with the operators `+ - * / ^`, comparisons, and the
  functions `min max exp log sqrt abs floor ceil if`. Michaelis–Menten
  kinetics is `"vmax * S / (Km + S)"` on `["S"]`. A string rate applies while
  each reactant has enough molecules and is 0 otherwise.

So `model.add_reaction(0.2, ["A"], [])` is the same reaction as the death
reaction above. Counts are integers. A reaction may take any number of
reactant molecules: `["A", "A", "A"]` is a third-order reaction. Tools differ
in how they count `A + A`, which matters when you take a rate constant from
a paper or from another program
([chapter 8](08-reproducibility.md#reproducing-a-published-simulation)).

## Checking a model

`crn.check` validates a model without running it. It returns a report with
the field `valid` and a list of diagnostics:

```python
report = crn.check(model)
print(report["valid"], report["diagnostics"])
```

```text
True []
```

Every problem is a diagnostic with a stable code, the path of the field in
the model or the request, a message and a hint. A function that cannot go on
raises `crn.CrnError`, a `ValueError` that carries the diagnostics. Here a
report time lies after the end of the run:

```python
try:
    crn.ensemble(model, end=1.0, report_times=[0.0, 2.0], method="static_direct")
except crn.CrnError as error:
    print(error)
    print(error.diagnostics[0]["path"])
    print(error.diagnostics[0]["hint"])
```

```text
CRN-R020: report times must be strictly increasing within [0, end] at /time/report_times/1
/time/report_times/1
Sort the grid and keep it inside the horizon.
```

Validation is strict. A misspelled parameter in a string rate, an unknown
field and a duplicate key are errors, never warnings, and nothing is guessed.
[Diagnostics](../reference/diagnostics.md) lists every code.

## A second model with two reactants

The SIR model of an epidemic has what immigration–death lacks: a reaction of
two molecules. A susceptible `S` meets an infected `I` and becomes infected,
and an infected recovers to `R`. Time is in days. The model is the second
example of this guide:

```python
sir = crn.Model("sir")
sir.add_species({"S": 490, "I": 10, "R": 0})  # or one at a time: add_species("S", 490)
sir.add_parameter("gamma", 0.1)
# A number is a mass-action constant: S + I fires at 5e-4 * S * I.
sir.add_reaction(5e-4, ["S", "I"], ["I", "I"], name="infection")
# A string is the whole propensity, over parameters and species counts.
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
print(sir)
run = crn.ensemble(
    sir,
    end=60.0,
    report_times=[0.0, 20.0, 40.0, 60.0],
    replicates=1000,
    seed=1,
    threads=2,
    method="static_direct",
)
result = run.simulate()
print(result.statistics["mean"].round(1))  # species columns, sorted by name
print("events per replicate:", result.rule_events.mean().to_dict())
```

```text
Model('sir': 3 species, 2 reactions)
         I      R      S
time                    
0.0   10.0    0.0  490.0
20.0  90.9   86.1  323.0
40.0  87.4  295.0  117.6
60.0  28.7  404.4   66.9
events per replicate: {'infection': 423.065, 'recovery': 404.412}
```

The epidemic peaks near day 20 and has infected most of the population by
day 60. Each replicate fires about 423 infections and 404 recoveries.

## The builder calls

| Call | Effect |
| --- | --- |
| `crn.Model(name="model")` | a new, empty model |
| `add_parameter(name, value)` or `add_parameter(mapping)` | parameters that string rates may read |
| `add_species(name, initial=0)` or `add_species(mapping)` | species and their initial counts |
| `set_initial(mapping)` or `set_initial(S=490, I=10)` | new initial counts of declared species |
| `add_reaction(rate, reactants, products, name=None)` | a reaction; names default to `r0`, `r1`, ... |
| `set_stop(condition)` | a condition that ends a replicate ([chapter 6](06-stops-trackers-logs.md)) |
| `add_ode(name, initial, derivative)` | a continuous variable ([chapter 9](09-continuous-variables.md)) |

The `add_*` calls return the model, so they can be chained. A mapping, as in
`add_parameter({"gamma": 0.1, "beta": 5e-4})`, declares its entries in order.
If one entry is invalid, the call raises and declares none of them. Species
that a reaction names without a declaration start at 0.

A `Model` is not a dictionary: `model["stop"] = ...` raises `TypeError`, and
`set_stop` is the call to use. `crn.ensemble` and `crn.check` take the
`Model` itself. What the builder writes, and how to save and load a model as
a file, is the subject of [chapter 14](14-model-files.md).

See also: [Running ensembles](03-running.md),
[Model language](../reference/model-language.md),
[Diagnostics](../reference/diagnostics.md).
