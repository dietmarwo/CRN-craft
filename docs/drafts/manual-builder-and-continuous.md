# Draft chapters for the mtDNA user manual (0.3.0)

These chapters are written for `docs/mtdna-user-manual.md` and are merged
into it in I7. Each `python` block runs on its own from the repository root.
Outputs were produced on the 2-vCPU development container with the wheel of
the engine stream; rerun them after the merge, since the result objects
change in I3. The blocks read results through `result.metadata`, the full
result document, which I3 keeps.

Placement: "Building models in Python" goes after chapter 2 (the models) or
before the chapter on model variants; "Continuous variables" replaces
chapter 10 ("Variable population: the ATP model"). The troubleshooting entry
replaces the existing `observation reservation exceeds budget` item.

---

## Building models in Python

`rules.Model` writes a rules model from Python calls, in the style of
rebop's `Gillespie` class. Every species is an attribute-free type, and
every species count is an observable with the species name.

```python
"""An SIR epidemic written with rules.Model."""

import numpy as np
from crn_craft import rules

sir = rules.Model("sir")
sir.add_species("S", 490).add_species("I", 10).add_species("R")
sir.add_parameter("gamma", 0.1)
# A number is a mass-action constant: S + I fires at 5e-4 * S * I.
sir.add_reaction(5e-4, ["S", "I"], ["I", "I"], name="infection")
# A string is the whole propensity, over parameters and species counts.
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
print(sir)
run = rules.ensemble(
    sir,
    end=60.0,
    report_times=[0.0, 20.0, 40.0, 60.0],
    replicates=1000,
    seed=1,
    threads=2,
    method="static_direct",
)
result = run.simulate()
meta = result.metadata
mean = np.asarray(meta["statistics"]["mean"])
names = meta["observables"]  # species, sorted by name
for t, row in zip(meta["report_times"], mean, strict=True):
    cells = "  ".join(f"{n} {v:6.1f}" for n, v in zip(names, row, strict=True))
    print(f"t {t:4.0f}: {cells}")
```

```text
Model('sir': 3 species, 2 reactions)
t    0: I   10.0  R    0.0  S  490.0
t   20: I   91.2  R   86.4  S  322.4
t   40: I   87.7  R  296.7  S  115.7
t   60: I   28.4  R  405.7  S   66.0
```

The calls:

| Call | Effect |
| --- | --- |
| `rules.Model(name="model")` | a new rules-model/0.3 document |
| `add_parameter(name, value)` | a parameter that string rates may read |
| `add_species(name, initial=0)` | a species type, its count observable and its initial count |
| `set_initial(mapping)` or `set_initial(S=490, I=10)` | initial counts of declared species |
| `add_reaction(rate, reactants, products, name=None)` | a rule; species are lists of names, repeated for multiplicity; undeclared species start at 0; names default to `r0`, `r1`, ... |
| `to_dict()`, `to_json()` | the JSON document |

`rules.ensemble` and `rules.prepare` accept the `Model` itself. Every method
runs builder models; `tree` is never needed, since the species have no
attributes.

**Rates.** A numeric rate k is a mass-action constant, counted like the rest
of CRN-craft: k times the number of distinct combinations of reactant
molecules. `A` fires at k·A, `A + B` at k·A·B and `A + A` at k·A(A − 1)/2.
rebop drops the 1/2 for identical reactants, so its rate k on `["A", "A"]`
equals `2 * k` here; reactions without identical reactants agree. A string
rate is the reaction's whole propensity (rules-model/0.3 `propensity`), an
expression over parameters and species counts with the operators
`+ - * / ^`, comparisons and the functions `min max exp log sqrt abs floor
ceil if`. It applies while each reactant has enough molecules (two for
`["A", "A"]`) and is 0 otherwise, so `"vmax * S / (Km + S)"` on `["S"]` stops
when `S` is used up. A reaction takes at most two reactant molecules.
`validation/check_builder.py` compares builder models with exact means and
with rebop 0.10.0 ([results](../validation/builder/RESULTS.md)).

The three-class model can be written the same way, with one species per
class and the vacancy law in the propensities. `examples/mtdna/12_builder.py`
builds it and compares it with the JSON file:

```python
"""The three-class model built with rules.Model."""

from crn_craft import rules

model = rules.Model("mtdna_three_species")
for name, value in {"K": 1000, "delta": 0.0028881132523331052, "sw": 1.0,
                    "sm": 1.5, "mu": 3e-05, "fracAdv": 0.2}.items():
    model.add_parameter(name, value)
model.add_species("W", 1000).add_species("M").add_species("A")
vacancy = "max(0, K - W - M - A) / (W + sm * (M + A))"
model.add_reaction(f"{vacancy} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt")
model.add_reaction(f"{vacancy} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt")
model.add_reaction(f"{vacancy} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv")
model.add_reaction(f"{vacancy} * sw * M", ["M"], ["M", "M"], "replication_mt")
model.add_reaction(f"{vacancy} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
for species in ("W", "M", "A"):
    model.add_reaction(0.0028881132523331052, [species], [], f"degradation_{species}")
result = rules.ensemble(
    model, end=25_920.0, report_times=[0.0, 25_920.0], replicates=100,
    seed=1, threads=2, method="mtdna_direct",
).simulate()
print(result.metadata["sampler"])
print(result.metadata["diagnostics"][0]["message"])
```

```text
chacha8-u53/rules-static-direct/1
mtdna_direct ran Static Direct: the model has whole-rule propensities
```

`mtdna_direct` certifies only the JSON form of the model; for builder
models it runs Static Direct and says so. Species with attributes (clone
identities, drawn fitness) need the JSON language.

---

## 10. Variable population: continuous variables

*Replaces "Variable population: the ATP model".*

In the variable-population model replication is limited by an ATP pool
instead of the vacancy `K − N`. ATP is a continuous quantity: between events
it follows

dA/dt = f · wt − v1 · A − v2 · N,

and the replication rules are multiplied by c/(c + A). Since 0.3.0 this is
written in the model file itself, as a *continuous variable*:
`models/rules/mtdna-three-class-atp-continuous.rules.json` (rules-model/0.3)
declares `ATP` with its initial value, its derivative and a lower bound,
and each replication rate reads it.

```python
"""Variable population with ATP as a continuous variable (hybrid method)."""

import json

import numpy as np
from crn_craft import rules

data = rules.load("models/rules/mtdna-three-class-atp-continuous.rules.json")
model = json.loads(data)
print("continuous:", model["continuous"])
print("replication_wt rate:", model["rules"][1]["rate"])

days = np.arange(0, 1081, 180)
run = rules.ensemble(
    data,
    end=25_920.0,
    report_times=24.0 * days,
    replicates=64,
    seed=7,
    threads=2,
    method="hybrid",  # rtol 1e-6, atol 1e-9 unless a method dict says otherwise
)
result = run.simulate()
meta = result.metadata
names = meta["observables"]
print("observables:", names, "continuous columns:", meta["continuous"])
print("sampler:", meta["sampler"], "method:", meta["method"])
print("stopped (ATP exhausted):", meta["stopped_replicates"])
mean = np.asarray(meta["statistics"]["mean"])
for day, row in zip(days, mean, strict=True):
    print(
        f"day {day:4d}: N {row[names.index('N')]:7.1f}  "
        f"mutants {row[names.index('mutants')]:6.1f}  ATP {row[names.index('ATP')]:7.2f}"
    )
```

```text
continuous: [{'name': 'ATP', 'initial': 'ATP0', 'derivative': 'f * wild_type - v1 * ATP - v2 * N', 'lower_bound': 0, 'on_bound': 'stop'}]
replication_wt rate: (sw * (1 - mu)) * c / (c + ATP)
observables: ['N', 'advantaged_mutants', 'mutants', 'neutral_mutants', 'wild_type', 'ATP'] continuous columns: ['ATP']
sampler: chacha8-u53/rules-approximate-hybrid-dopri54/1 method: {'atol': 1e-09, 'kind': 'hybrid', 'rtol': 1e-06}
stopped (ATP exhausted): 4
day    0: N  1000.0  mutants    0.0  ATP  450.00
day  180: N   998.6  mutants    0.8  ATP  449.13
day  360: N  1005.9  mutants    4.3  ATP  450.41
day  540: N  1607.2  mutants  561.5  ATP  445.71
day  720: N  2778.9  mutants 1660.0  ATP  428.70
day  900: N  3449.8  mutants 2282.5  ATP  422.58
day 1080: N  3493.0  mutants 2325.0  ATP  420.51
```

A continuous variable has

| Field | Meaning |
| --- | --- |
| `name` | its name; rates and derivatives read it like an observable, and results report it as an observable column after the model's observables |
| `initial` | its value at time 0, a number or a parameter expression |
| `derivative` | dx/dt between events: an expression over parameters, observables (here `wild_type` and `N`) and continuous variables |
| `lower_bound`, `on_bound` | optional, together: a bound and `"stop"` or `"clamp"` |

With `"stop"` a cell ends as `stopped` when the variable reaches its bound
while its derivative there is not positive, as in the ATP model when ATP is
exhausted; the variable then reads the bound on every later report, and the
other observables keep their values. With `"clamp"` the variable is held at
the bound and follows the ODE again once its derivative turns positive.

**The `hybrid` method.** Only `method="hybrid"` runs models with continuous
variables; the other methods refuse them with `CRN-R020`. Between events it
integrates the variables together with the integrated total propensity, with
the adaptive Dormand–Prince 5(4) Runge–Kutta pair, and places the next event
where the integrated propensity reaches an exponential threshold. It is
approximate at its tolerances `rtol` and `atol`, which the request states and
`metadata["method"]` records. The defaults `rtol=1e-6` and `atol=1e-9` apply
unless you pass `method={"kind": "hybrid", "rtol": ..., "atol": ...}`. The
request uses rules-run/0.5, and the result is rules-result/0.3 with the
`method` and the `continuous` column names. `validation/check_hybrid.py`
checks `hybrid` against exact laws and against exact thinning on this model
([results](../validation/hybrid/RESULTS.md)). On the test host it took about
eleven times as long as thinning on the ATP model, since every event needs an
integration step.

**The legacy file.** `mtdna-three-class-atp.rules.json` (rules-model/0.2)
writes the same model with a fixed `atp` block, which `method="thinning"`
samples exactly. It keeps working unchanged; a model cannot have both an
`atp` block and `continuous` variables.

**Your own variables.** `rules.Model.add_continuous(name, initial,
derivative, lower_bound=None, on_bound="stop")` adds a variable to a builder
model; string rates may read it:

```python
"""A continuous variable in a builder model."""

import numpy as np
from crn_craft import rules

cell = rules.Model("atp_wild_type")
for name, value in {"f": 0.1, "v1": 0.2, "v2": 0.01, "c": 1.3034153747670567}.items():
    cell.add_parameter(name, value)
cell.add_species("W", 1000)
cell.add_continuous("ATP", 450.0, "f * W - v1 * ATP - v2 * W", lower_bound=0)
cell.add_reaction("c / (c + ATP) * W", ["W"], ["W", "W"], name="replication")
cell.add_reaction(0.0028881132523331052, ["W"], [], name="degradation")
result = rules.ensemble(
    cell,
    end=2_400.0,
    report_times=[0.0, 1_200.0, 2_400.0],
    replicates=16,
    seed=3,
    threads=2,
    method={"kind": "hybrid", "rtol": 1e-7, "atol": 1e-9},
).simulate()
print("method:", result.metadata["method"])
print("observables:", result.metadata["observables"])
print("mean:", np.round(np.asarray(result.metadata["statistics"]["mean"]), 2).tolist())
```

```text
method: {'atol': 1e-09, 'kind': 'hybrid', 'rtol': 1e-07}
observables: ['W', 'ATP']
mean: [[1000.0, 450.0], [1000.06, 449.55], [1005.75, 452.43]]
```

Rules for continuous variables:

- Rule rates and propensities may read them; the model's `stop`, observable
  filters and weights, reactant patterns and products may not (`CRN-R012`).
  Stop on a variable with its `lower_bound`.
- Trackers follow discrete observables only (`CRN-R020` for a continuous
  one). Event logs and callbacks see the variables' values at each event.
- `hybrid` lowers the model to fixed tuples like `static_direct`, so models
  whose instances change identity (the clone model) need `tree`, which has
  no continuous variables.

For the method table in chapter 8, add:

| `method` | Kind | Use it for |
| --- | --- | --- |
| `"hybrid"` or `{"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}` | approximate at its tolerances | models with continuous variables (section 10) |

---

## Troubleshooting entry (replaces the reservation item)

- **`observation reservation exceeds budget`.** Event logs reserve space for
  every cell before the run. The message itemizes the reservation (cells ×
  bytes per cell, log entries × bytes per entry), and its hint names the
  smallest `max_output_bytes` that admits the request. Raise
  `max_output_bytes` to that value, or use fewer cells or entries. For the
  three-class model an entry takes 721 bytes, so 50 entries for 1,000 cells
  fit the 64 MiB default.
