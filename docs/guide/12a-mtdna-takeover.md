# Takeover in a fixed population

How many cells are taken over by mutants, and when? This page runs the
three-class model of the [index page](12-mtdna.md#the-model) with a fixed
copy-number target. The setup block is that model. The blocks of this page
continue from it:

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

## A first ensemble

Applies: [chapter 3, Many trajectories in one call](03-running.md#many-trajectories-in-one-call), and [chapter 4, What a result holds](04-results.md#what-a-result-holds).

Two hundred cells for three years, with the takeover stop at 90% mutants.
Every cell is one replicate, so one call runs them all:

```python
import time

model = three_class(takeover=0.9)
days = np.arange(0, 1081, 90)
run = crn.ensemble(model, end=END, report_times=24.0 * days, replicates=200,
                   seed=1, threads=2, method="mtdna_direct")
start = time.perf_counter()
result = run.simulate()
seconds = time.perf_counter() - start

meta = result.metadata
records = result.records  # one row per cell
stopped = records["termination"] == "stopped"
print(meta["sampler"], meta["termination"])
print(f"{meta['events']:,} events in {seconds:.2f} s")
print(f"taken over: {stopped.sum()} of {len(records)} cells")
print(records[["termination", "stop_time", "events"]].head(3))
first = records.loc[stopped, "stop_time"].nsmallest(5)
print("takeover days (first five):", (first.to_numpy() / 24).round(1))
mean = result.statistics["mean"]  # one row per report time
fraction = (mean["M"] + mean["A"]) / mean.sum(axis=1)
for hours, value in fraction.items():
    print(f"day {hours / 24:5.0f}: mean mutant fraction {value:.4f}")
```

```text
chacha8-u53-ziggurat/rules-mtdna-direct/1 completed
28,536,814 events in 0.50 s
taken over: 23 of 200 cells
          termination  stop_time  events
replicate                               
0           completed    25920.0  148726
1           completed    25920.0  149476
2           completed    25920.0  150111
takeover days (first five): [272.7 300.8 304.3 358.8 426. ]
day     0: mean mutant fraction 0.0000
day    90: mean mutant fraction 0.0004
day   180: mean mutant fraction 0.0055
day   270: mean mutant fraction 0.0156
day   360: mean mutant fraction 0.0242
day   450: mean mutant fraction 0.0384
day   540: mean mutant fraction 0.0460
day   630: mean mutant fraction 0.0516
day   720: mean mutant fraction 0.0620
day   810: mean mutant fraction 0.0838
day   900: mean mutant fraction 0.1011
day   990: mean mutant fraction 0.1065
day  1080: mean mutant fraction 0.1097
```

A cell has about 150,000 events in three years. The cells that were taken
over stopped early and hold their state at every later report time. The
mean mutant fraction is the mean of `M + A` over the mean of all molecules
at each report time.

## The mean hides the typical cell

Applies: [chapter 4, Quantiles and the median](04-results.md#quantiles-and-the-median), and [Single trajectories](04-results.md#single-trajectories).

Most cells stay near zero mutants, while a few are taken over. The mean then
sits far from the typical cell. The quantiles of one species come from
`statistics`. The quantiles of all mutants, `M + A`, are quantiles of a sum,
so they come from the kept paths:

```python
days = np.arange(0, 1081, 360)
result = crn.ensemble(
    three_class(), end=END, report_times=24.0 * days, replicates=200, seed=1,
    threads=2, method="mtdna_direct",
    quantiles=[0.025, 0.25, 0.5, 0.75, 0.975], keep_paths="all",
).simulate()
statistics = result.statistics
mean = statistics["mean"]
print("quantile levels:", list(statistics["quantiles"]))
print("mean mutants:     ", (mean["M"] + mean["A"]).round(1).tolist())
print("median advantaged:", statistics["median"]["A"].tolist())
paths = result.paths  # cells x report times x species
mutants = (paths.sel(variable="M") + paths.sel(variable="A")).to_numpy()
print("median mutants:   ", np.median(mutants, axis=0).tolist())
print("97.5% of mutants: ", np.quantile(mutants, 0.975, axis=0).round(1).tolist())
last = paths.sel(time=END, variable="A").to_numpy()
print("equal to numpy:", np.quantile(last, 0.975) == statistics["quantiles"][0.975]["A"].iloc[-1])
```

```text
quantile levels: [0.025, 0.25, 0.5, 0.75, 0.975]
mean mutants:      [0.0, 25.2, 67.0, 120.7]
median advantaged: [0.0, 0.0, 0.0, 0.0]
median mutants:    [0.0, 0.0, 0.0, 0.0]
97.5% of mutants:  [0.0, 336.4, 995.0, 998.0]
equal to numpy: True
```

The median is zero at every report time, while the mean rises: the typical
cell has no mutant, and the mean is carried by the cells that were taken
over. The 97.5% quantile shows how far a cell can go. The quantiles that
the run computes equal NumPy's over the kept paths, so both routes give the
same numbers.

## Takeover fraction and time to threshold

Applies: [chapter 6, Stop conditions](06-stops-trackers-logs.md#stop-conditions).

Two outputs of a stop condition answer the question of this page: the share
of cells whose stop time lies before an age, and the distribution of the
stop times:

```python
result = crn.ensemble(three_class(takeover=0.9), end=END, report_times=[0.0, END],
                      replicates=200, seed=1, threads=2, method="mtdna_direct").simulate()
stopped = result.records["termination"] == "stopped"
for age in (360, 720, 1080):
    taken = (stopped & (result.records["stop_time"] <= 24.0 * age)).mean()
    print(f"taken over by day {age:4d}: {taken:.3f}")
days = result.records.loc[stopped, "stop_time"] / 24
print(f"days to 90% mutants: first {days.min():.0f}, median {days.median():.0f}")
```

```text
taken over by day  360: 0.020
taken over by day  720: 0.055
taken over by day 1080: 0.115
days to 90% mutants: first 273, median 726
```

These are the 200 cells of the first ensemble: the same seed gives the same
cells, whatever the report times. Another threshold is another value of
`takeover`, for example 0.6 for a 60% deficiency threshold. To record when a
threshold is passed without ending the cell, use a first-passage tracker
([page 12b](12b-mtdna-events.md#class-appearances-and-extinctions)).

## Variants

Applies: [chapter 2, The builder calls](02-first-model.md#the-builder-calls), and [Checking a model](02-first-model.md#checking-a-model).

A variant of the model is a call of `three_class` with other arguments:
copy number `K`, half-life, initial heteroplasmy (`w0` wild-type and `m0`
neutral-mutant molecules) and the threshold. Here a smaller cell with a
longer half-life starts with 20% neutral mutants, and a cell counts as taken
over at 60%:

```python
model = three_class(K=500, w0=400, m0=100, half_life_days=20, takeover=0.6)
result = crn.ensemble(model, end=END, report_times=[0.0, 8640.0, 17_280.0, END],
                      replicates=300, seed=9, threads=2, method="mtdna_direct").simulate()
stopped = result.records["termination"] == "stopped"
print(f"cells past 60% within three years: {stopped.mean():.2%}")

# A misspelled name is reported with a code and a path before anything runs.
broken = three_class()
broken.add_reaction("delta * WW", ["W"], [], "extra_loss")
try:
    crn.ensemble(broken, end=1.0, report_times=[0.0, 1.0])
except crn.CrnError as error:
    print("rejected:", error)
```

```text
cells past 60% within three years: 4.33%
rejected: CRN-R012: unknown name `WW`; expected a parameter or observable in `delta * WW` at /rules/8/propensity
```

A builder model is checked when it is used, and a mistake raises
`crn.CrnError` with a code, a path and the offending expression. To change a
parameter of many runs without building the model again, bind it
([page 12e](12e-mtdna-large-studies.md#parameter-sweeps)).

See also: [Mutant classes](12b-mtdna-events.md),
[Large studies](12e-mtdna-large-studies.md).
