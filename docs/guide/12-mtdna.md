# 12. Case study: mtDNA mutants

## What this case study covers

Applies: [chapter 1, How to read this guide](01-install.md#how-to-read-this-guide).

A cell holds about a thousand copies of its mitochondrial DNA. Copies
degrade and are replaced by replication, and a replication can produce a
mutant. Some mutants replicate faster than the wild type. Over years, a
mutant class can take over a cell. The case study asks the questions a
researcher asks of such a model and answers each with a chapter of this
guide. It explains no API a second time: each section names the chapter it
applies.

| Page | Question | Applies |
| --- | --- | --- |
| [Takeover in a fixed population](12a-mtdna-takeover.md) | How many cells are taken over, and when? | chapters 3, 4, 6 |
| [Mutant classes](12b-mtdna-events.md) | How often do mutant classes appear and vanish? Which mutations happened when? | chapter 6 |
| [Variable population: ATP](12c-mtdna-atp.md) | What if an energy balance, not a fixed target, limits replication? | chapter 9 |
| [Methods and speed](12d-mtdna-methods.md) | Which method, and how fast? | chapters 5, 7 |
| [Large studies](12e-mtdna-large-studies.md) | How do parameters change takeover? How to run a study as files? | chapters 3, 5, 8, 10 |
| [Fitting parameters to data](13-fitting.md) | Which mutation probability and advantage explain measured cells? | chapter 5 |
| [The mtDNA model files](14a-mtdna-model-files.md) | How to follow every lineage, and run the fast path? | chapter 14 |

## The model

Applies: [chapter 2, Species, parameters and reactions](02-first-model.md#species-parameters-and-reactions).

One cell has a copy-number target `K`. Every molecule degrades at a
constant rate. Every molecule replicates into the free places below `K`:
the fewer free places, the slower. When a wild-type molecule replicates,
the copy is a mutant with a small probability `mu`. A share `fracAdv` of new
mutants is *advantaged*: it replicates with the weight `sm` instead of 1.

The model has three species: `W` (wild type), `M` (neutral mutants) and `A`
(advantaged mutants). The function `three_class` builds it with the builder
([chapter 2](02-first-model.md)). Every page of the case study starts with
this block:

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

A hundred cells for three years:

```python
model = three_class()
print(model)
result = crn.ensemble(model, end=END, report_times=[0.0, END], replicates=100,
                      seed=1, threads=2, method="static_direct").simulate()
final = result.statistics["mean"].loc[END]
print(result.metadata["sampler"], result.metadata["termination"])
print(f"mean mutant fraction after three years: {(final['M'] + final['A']) / final.sum():.4f}")
```

```text
Model('mtdna_three_class': 3 species, 8 reactions)
chacha8-u53-ziggurat/rules-static-direct/1 completed
mean mutant fraction after three years: 0.1108
```

## Parameters and rates

Applies: [chapter 2, Rates as numbers and as strings](02-first-model.md#rates-as-numbers-and-as-strings).

| Argument of `three_class` | Parameter | Default | Meaning |
| --- | --- | ---: | --- |
| `K` | `K` | 1000 | copy-number target |
| `half_life_days` | `delta` | ln 2 / 240 h⁻¹ | degradation rate per molecule (half-life 10 days) |
| | `sw` | 1.0 | replication weight of wild type and neutral mutants |
| `sm` | `sm` | 1.5 | replication weight of advantaged mutants |
| `mu` | `mu` | 3·10⁻⁵ | mutation probability per wild-type replication |
| | `fracAdv` | 0.2 | share of new mutants that are advantaged |
| `w0`, `m0` | | 1000, 0 | initial wild-type and neutral-mutant molecules |
| `takeover` | | none | the takeover threshold, a fraction of all molecules |

Every class degrades at rate `delta`. A molecule of weight `s` replicates
at rate `max(0, K − N) / D · s`, where `N = W + M + A` counts all molecules
and `D = W + sm (M + A)` is the competition denominator. A wild-type
replication yields an unchanged copy with probability 1 − `mu`, a neutral
mutant with probability `mu` (1 − `fracAdv`) and an advantaged mutant with
probability `mu` · `fracAdv`. The eight reactions write this out channel by
channel.

**Totals are sums of columns.** The builder gives one column per species
and no other observable. All molecules are `W + M + A`, all mutants
`M + A`. Write them out in stop conditions and first-passage predicates,
and add the columns in pandas for means. The variance or a quantile of a
total needs the kept paths ([chapter 4](04-results.md#single-trajectories)).

**Takeover** is the stop condition that `takeover=0.9` sets: a cell is
taken over when 90% of its molecules are mutants,
`W + M + A > 0 && M + A >= 0.9 * (W + M + A)`. Other thresholds are other
values of `takeover`.

## Conventions

Applies: [chapter 7, Which method for which model](07-methods.md#which-method-for-which-model).

- Time is in hours. `END`, 25,920 hours, is 1,080 days, about three years.
  Report times are written as `24.0 * days`.
- One replicate is one cell. A cell at `K = 1000` has about 150,000 events
  in three years.
- The method is `static_direct`, which runs every variant on these pages.
  The three-class model also has a fast path, `mtdna_direct`, and since
  release 0.7.1 it certifies the builder model as well as the model file
  when the request is `crn-craft.run/1.0`, which `crn.ensemble` writes
  ([page 12d](12d-mtdna-methods.md#the-fast-path)). On 13 November 2026, on
  a two-CPU host, the builder model under the fast path took the time of
  the model file under it (0.97 times), and under Static Direct 3.6 times
  as long, and 1.4 times as long as the model file under Static Direct
  ([measurement](https://github.com/dietmarwo/CRN-craft/blob/main/validation/merge/RESULTS.md)).
  For long runs of the plain three-class model, use the fast path.

## The model as a file

Applies: [chapter 14, Save and load a model file](14-model-files.md#save-and-load-a-model-file).

The repository also has this model as JSON files, in which every molecule
carries attributes. Those files add a model that follows every mutant
lineage and a stored form of the ATP model;
[page 14a](14a-mtdna-model-files.md), after the chapter on model files,
shows them.

See also: [Takeover in a fixed population](12a-mtdna-takeover.md),
[mtDNA performance analysis](../evidence/mtdna-performance-analysis.md).
