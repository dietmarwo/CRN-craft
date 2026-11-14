# Variable population: ATP

In the model of the index page a fixed target `K` limits replication. In
this variant an energy balance limits it: a pool of ATP, which wild-type
molecules produce and every molecule consumes. The population is then free
to grow or shrink. The setup block is the model of the
[index page](12-mtdna.md#the-model); this page builds its variant from the
same parts. The blocks of this page continue from it:

```python
import math

import numpy as np
import crn_craft as crn

END = 25_920.0  # three years (1,080 days), in hours


def three_class(K=1000, w0=1000, m0=0, mu=3e-5, sm=1.5, half_life_days=10, takeover=None):
    """The three-class mtDNA model: wild type W, neutral M, advantaged A."""
    model = crn.Model("mtdna_three_class")
    model.add_parameter({"K": K, "delta": math.log(2) / (24 * half_life_days),
                         "sw": 1.0, "sm": sm, "mu": mu, "fracAdv": 0.2})
    model.add_species({"W": w0, "M": m0, "A": 0})
    free = "max(0, K - W - M - A) / (W + sm * (M + A))"
    model.add_reaction(f"{free} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt")
    model.add_reaction(f"{free} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt")
    model.add_reaction(f"{free} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv")
    model.add_reaction(f"{free} * sw * M", ["M"], ["M", "M"], "replication_mt")
    model.add_reaction(f"{free} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
    for label, species in (("wt", "W"), ("mt", "M"), ("mtAdv", "A")):
        model.add_reaction(f"delta * {species}", [species], [], f"degradation_{label}")
    if takeover is not None:
        model.set_stop(f"W + M + A > 0 && M + A >= {takeover} * (W + M + A)")
    return model
```

## ATP as a continuous variable

Applies: [chapter 9, A variable between events](09-continuous-variables.md#a-variable-between-events).

ATP is a continuous quantity. Between events it follows

dATP/dt = f · W − v1 · ATP − v2 · (W + M + A),

and every replication rate is multiplied by c / (c + ATP) instead of the
free places below `K`. `add_ode` declares the variable with a lower bound of
0, and the string rates read it. Only the method `hybrid` runs such a model.
The fast path of the other pages does not apply, because this model is not
the three-class structure:

```python
def atp_cell(f=0.1, v1=0.2, v2=0.01, c=1.3034153747670567, atp0=450.0):
    """The three-class model with replication limited by an ATP balance."""
    model = crn.Model("mtdna_three_class_atp")
    model.add_parameter({"delta": math.log(2) / 240, "sw": 1.0, "sm": 1.5, "mu": 3e-5,
                         "fracAdv": 0.2, "f": f, "v1": v1, "v2": v2, "c": c})
    model.add_species({"W": 1000, "M": 0, "A": 0})
    model.add_ode("ATP", atp0, "f * W - v1 * ATP - v2 * (W + M + A)", lower_bound=0)
    energy = "c / (c + ATP)"
    model.add_reaction(f"{energy} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt")
    model.add_reaction(f"{energy} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt")
    model.add_reaction(f"{energy} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv")
    model.add_reaction(f"{energy} * sw * M", ["M"], ["M", "M"], "replication_mt")
    model.add_reaction(f"{energy} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
    for label, species in (("wt", "W"), ("mt", "M"), ("mtAdv", "A")):
        model.add_reaction(f"delta * {species}", [species], [], f"degradation_{label}")
    return model


result = crn.ensemble(atp_cell(), end=END, report_times=24.0 * np.arange(0, 1081, 180),
                      replicates=64, seed=7, threads=2, method="hybrid").simulate()
meta = result.metadata
print("columns:", list(result.statistics["mean"].columns))
print("sampler:", meta["sampler"], "method:", meta["method"])
print("stopped (ATP exhausted):", (result.records["termination"] == "stopped").sum())
for hours, row in result.statistics["mean"].iterrows():
    print(f"day {hours / 24:4.0f}: molecules {row['W'] + row['M'] + row['A']:7.1f}  "
          f"mutants {row['M'] + row['A']:6.1f}  ATP {row['ATP']:7.2f}")
```

```text
columns: ['A', 'M', 'W', 'ATP']
sampler: chacha8-u53/rules-approximate-hybrid-dopri54/1 method: {'atol': 1e-09, 'kind': 'hybrid', 'rtol': 1e-06}
stopped (ATP exhausted): 10
day    0: molecules  1000.0  mutants    0.0  ATP  450.00
day  180: molecules  1006.8  mutants    4.5  ATP  450.32
day  360: molecules  1826.1  mutants  771.8  ATP  438.88
day  540: molecules  2963.6  mutants 1835.6  ATP  424.91
day  720: molecules  4570.4  mutants 3328.0  ATP  409.96
day  900: molecules  5191.7  mutants 3908.4  ATP  402.75
day 1080: molecules  7011.7  mutants 5611.1  ATP  379.08
```

`ATP` is a column after the species. The population is no longer held near
1,000: where mutants expand, it grows, because mutants consume ATP but do
not produce it, and the replication factor c / (c + ATP) rises as ATP falls.
`metadata["method"]` records the tolerances of `hybrid`, which is
approximate at them.

## What ends a cell

Applies: [chapter 9, Bounds](09-continuous-variables.md#bounds).

`lower_bound=0` with the default `on_bound="stop"` ends a cell when ATP
reaches 0 while its derivative there is not positive: the energy of the cell
is exhausted. The cell is then `stopped`, ATP reads 0 at every later report
time, and the species keep their counts. With `on_bound="clamp"` ATP would
instead stay at 0 until its derivative turns positive, and the cell would
live on. On the stored form of this model, `hybrid` took about eleven times
as long as exact thinning, because every event needs an integration step;
its error has not been measured as a function of the tolerances
([validation](https://github.com/dietmarwo/CRN-craft/blob/main/validation/hybrid/RESULTS.md)).

## The stored form

Applies: [chapter 7, Other methods](07-methods.md#other-methods).

The repository also has this model as a stored file with an older `atp`
block. Such files run exactly under the method `thinning`;
[page 14a](14a-mtdna-model-files.md#the-stored-atp-block-and-thinning)
shows it.

See also: [Continuous variables](09-continuous-variables.md),
[Methods and speed](12d-mtdna-methods.md).
