# 9. Continuous variables

Some quantities are not counts of molecules: an energy level, a drug
concentration, a temperature, a hazard that grows with age. A model can
carry such a quantity as a *continuous variable*. It follows an ordinary
differential equation between the events, and reaction rates may read it.
The setup block builds a population `W` whose replication is slowed by a
resource `ATP` that the population itself produces and uses. The blocks of
this chapter continue from it:

```python
import math

import crn_craft as crn

cell = crn.Model("resource")
cell.add_parameter({"f": 0.1, "v1": 0.2, "v2": 0.01, "c": 1.3034153747670567})
cell.add_species("W", 1000)
cell.add_ode("ATP", 450.0, "f * W - v1 * ATP - v2 * W", lower_bound=0)
cell.add_reaction("c / (c + ATP) * W", ["W"], ["W", "W"], name="replication")
cell.add_reaction(0.0028881132523331052, ["W"], [], name="degradation")
```

## A variable between events

`add_ode(name, initial, derivative)` declares the variable. Above, `ATP`
starts at 450 and changes at the rate `f * W - v1 * ATP - v2 * W`, which
reads the count `W`. The replication rate reads `ATP` in turn. A model with
a continuous variable runs under the method `hybrid`:

```python
result = crn.ensemble(
    cell,
    end=2_400.0,
    report_times=[0.0, 1_200.0, 2_400.0],
    replicates=16,
    seed=3,
    threads=2,
    method={"kind": "hybrid", "rtol": 1e-7, "atol": 1e-9},
).simulate()
print("method:", result.metadata["method"])
print(result.statistics["mean"].round(2))
```

```text
method: {'atol': 1e-09, 'kind': 'hybrid', 'rtol': 1e-07}
              W     ATP
time                   
0.0     1000.00  450.00
1200.0  1000.06  449.55
2400.0  1005.75  452.43
```

Results report the variable as a column after the species. The arguments of
`add_ode`:

| Argument | Meaning |
| --- | --- |
| `name` | the name; rates and derivatives read it like a species count |
| `initial` | the value at time 0, a number or an expression over parameters |
| `derivative` | dx/dt between events, an expression over parameters, species counts and continuous variables |
| `lower_bound`, `on_bound` | optional: a bound, and `"stop"` (the default) or `"clamp"` |
| `independent` | `True` for a variable that does not depend on the molecules (below) |

## The hybrid method

Only `method="hybrid"` runs a model with continuous variables. The other
methods refuse it with `CRN-R020`. Between events, `hybrid` integrates the
variables together with the total rate, with an adaptive Runge–Kutta method
(Dormand–Prince 5(4)), and places the next event where the integrated rate
reaches an exponentially distributed threshold.

The method is approximate at its tolerances. `method="hybrid"` uses
`rtol=1e-6` and `atol=1e-9`, and the dictionary form above sets others. The
result records them in `metadata["method"]`, and its sampler has
`approximate` in its name. To see whether the tolerances matter, run again
with tighter ones and compare. Each event costs an integration step, so
`hybrid` is slower than Static Direct: use it when a rate reads a continuous
variable, not otherwise.
[Methods and their options](../reference/methods.md) has the algorithm.

## Bounds

A variable may have a lower bound, and `on_bound` says what happens there:

- `"stop"` ends the replicate as `stopped` when the variable reaches the
  bound while its derivative there is not positive. The variable reads the
  bound at every later report time, and the species keep their counts. In
  the model above, a population that runs out of `ATP` ends.
- `"clamp"` holds the variable at the bound and lets it follow the equation
  again once its derivative turns positive.

A bound is the way to stop on a continuous variable. The stop condition of
`set_stop` reads species counts only.

## Independent variables

Some variables do not depend on the molecules at all: a daily rhythm, a dose
that decays on its own, a hazard that grows with age. Declare such a
variable with `independent=True`. Its `initial` and `derivative` may then
read only parameters and other independent variables. `hybrid` integrates
it once for the whole run, and every replicate reads that one solution. In
this model the death rate `m` grows as dm/dt = g m, whatever the molecules
do (Gompertz mortality):

```python
gompertz = crn.Model("gompertz")
gompertz.add_parameter({"g": 0.5, "m0": 0.1})
gompertz.add_species({"A": 20})
gompertz.add_ode("m", "m0", "g * m", independent=True)
gompertz.add_reaction("m * A", ["A"], [], name="death")
result = crn.ensemble(
    gompertz,
    end=4.0,
    report_times=[0.0, 2.0, 4.0],
    replicates=4000,
    seed=1,
    threads=2,
    method="hybrid",
).simulate()
print(result.metadata["continuous"], result.metadata["sampler"])
print(result.statistics["mean"].round(3))
survival = math.exp(-0.1 * (math.exp(2.0) - 1.0) / 0.5)
print(f"exact mean at t = 4: {20 * survival:.3f}")

# An independent derivative may read only parameters and independent variables.
wrong = crn.Model("wrong")
wrong.add_parameter("g", 0.5)
wrong.add_species("A", 20)
wrong.add_ode("m", 0.1, "g * m * A", independent=True)
wrong.add_reaction("m * A", ["A"], [], name="death")
try:
    crn.ensemble(wrong, end=4.0, report_times=[4.0], method="hybrid")
except crn.CrnError as error:
    print(error.diagnostics[0]["code"], error.diagnostics[0]["path"])
    print(error.diagnostics[0]["message"])
```

```text
['m'] chacha8-u53/rules-approximate-hybrid-dopri54/1
           A      m
time               
0.0   20.000  0.100
2.0   14.154  0.272
4.0    5.585  0.739
exact mean at t = 4: 5.573
CRN-R012 /continuous/0/derivative
the derivative of independent variable `m` reads observable `A`
```

Each molecule survives to time t with probability exp(−∫m), so the count at
t = 4 is binomial with mean 20 × 0.2786 = 5.573. The 4,000 replicates give
5.585.

An independent variable costs one integration per run, not one per
replicate, which matters when the forcing is large or stiff. A bound keeps
its meaning: with `"stop"`, every replicate that is still running stops at
the same time, because the variable is the same in all of them.

## Rules for continuous variables

- Rates may read them. The stop condition of `set_stop` may not
  (`CRN-R012`): stop on a variable with its `lower_bound`.
- Trackers follow species only (`CRN-R020` for a continuous variable). Event
  logs and callbacks see the variables' values at each event.
- A bound is checked along the integration steps. A variable that dips below
  its bound and returns within one step can be missed, and tighter
  tolerances shorten the steps.
- Compiled models and finite state projection do not take continuous
  variables.

See also: [Choosing a method](07-methods.md),
[Variable population: ATP](12c-mtdna-atp.md),
[Model language](../reference/model-language.md#continuous-variables-rules-model03).
