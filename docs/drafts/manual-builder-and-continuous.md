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
