# Mutant classes: trackers, event logs and callbacks

How often does a mutant class appear and vanish in a cell? When did the
mutations happen? These questions are about events between report times.
This page answers them with the tools of
[chapter 6](06-stops-trackers-logs.md). The setup block is the model of the
[index page](12-mtdna.md#the-model). The blocks of this page continue from
it:

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

## Class appearances and extinctions

Applies: [chapter 6, Trackers](06-stops-trackers-logs.md#trackers).

One tracker follows the neutral class `M`, one the advantaged class `A`.
A third records the first time mutants reach 60% of the molecules, written
over the three species, without stopping the cell:

```python
sixty = "W + M + A > 0 && M + A >= 0.6 * (W + M + A)"
result = crn.ensemble(
    three_class(takeover=0.9), end=END, report_times=[0.0, END], replicates=500,
    seed=3, threads=2, method="mtdna_direct",
    trackers=[
        {"observable": "M", "keep_episodes": 3},
        {"observable": "A", "keep_episodes": 3},
        {"observable": "W", "first_passage": sixty},
    ],
).simulate()
neutral, advantaged, passage = (result.tracker_columns(i) for i in range(3))
print(advantaged[["appearances", "extinctions", "maximum", "integral"]].head(3))
print("neutral classes arising per cell:", neutral["appearances"].mean())
print("neutral classes lost per cell:   ", neutral["extinctions"].mean())
print("advantaged arising/lost:         ",
      advantaged["appearances"].mean(), advantaged["extinctions"].mean())
print("largest advantaged class seen:   ", advantaged["maximum"].max())
# The integral runs until the cell's stop or the horizon.
load = advantaged["integral"] / result.records["stop_time"]
print("time-averaged advantaged count, mean over cells:", round(load.mean(), 3))
reached = passage["first_passage"].dropna()  # NaN: never reached 60%
print(f"cells reaching 60% mutants: {len(reached)}; median day {reached.median() / 24:.0f}")
first = next(o for o in result.observations if o["trackers"][1]["episodes"])
for episode in first["trackers"][1]["episodes"]:
    print("cell", first["id"], "advantaged episode:", episode)
```

```text
           appearances  extinctions  maximum   integral
replicate                                              
0                    0            0      0.0   0.000000
1                    1            1      1.0  64.109016
2                    0            0      0.0   0.000000
neutral classes arising per cell: 1.566
neutral classes lost per cell:    1.472
advantaged arising/lost:          0.404 0.266
largest advantaged class seen:    900.0
time-averaged advantaged count, mean over cells: 17.69
cells reaching 60% mutants: 52; median day 494
cell 1 advantaged episode: {'duration': 64.10901603704042, 'end': 5246.85151399302, 'peak': 1.0, 'start': 5182.742497955979}
```

An appearance of a class is the birth of a mutant lineage *while its class
is absent*. A new mutant that arises while its class is present merges into
it, so a class count is not a lineage count. The integral divided by the
observed time is the time-averaged count: here the average number of
advantaged molecules per hour lived. The first passage at 60% is recorded
for cells that went on to 90% as well as for cells that stayed below it.

## Mutation events in a log

Applies: [chapter 6, Event logs](06-stops-trackers-logs.md#event-logs).

An event log keeps up to `max_entries` events of the listed reactions for
every cell, with the species counts after each event. The two mutation
reactions are the ones to keep:

```python
result = crn.ensemble(
    three_class(), end=END, report_times=[0.0, END], replicates=100, seed=4,
    threads=2, method="mtdna_direct",
    event_log={"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]},
).simulate()
log = result.event_columns()  # cells x entries x species, NaN-padded
names = log["rule_name"].to_series().dropna()  # padding has no name
print("logged mutation events:", names.value_counts().sort_index().to_dict())
print("omitted (log full):", sum(o["omitted_log_entries"] for o in result.observations))
print(result.event_frame(0)[["time", "rule_name", "W", "M", "A"]])
days = result.event_columns(0)["time"] / 24  # cell 0, without padding
print("cell 0 mutation days:", days.round(1).values.tolist())
```

```text
logged mutation events: {'mutation_mt': 158, 'mutation_mtAdv': 39}
omitted (log full): 0
               time       rule_name      W    M      A
entry                                                 
0       1442.924093     mutation_mt  994.0  1.0    0.0
1       7746.367689     mutation_mt  996.0  1.0    0.0
2      17740.279613  mutation_mtAdv  997.0  0.0    1.0
3      21692.680185     mutation_mt  744.0  1.0  252.0
cell 0 mutation days: [60.1, 322.8, 739.2, 903.9]
```

No log was full: 50 entries per cell are more than three years of mutations
at `mu = 3·10⁻⁵`. Each entry holds the time, the reaction and the counts of
`W`, `M` and `A` after the event.

## Stop at the first advantaged mutant

Applies: [chapter 6, Callbacks](06-stops-trackers-logs.md#callbacks).

A callback sees every event of a cell in order and may end the cell. This
one counts events by reaction and stops a cell at its first advantaged
mutation. A callback holds Python's interpreter lock, so it suits a few
cells:

```python
from collections import Counter

run = crn.ensemble(three_class(), end=END, report_times=[0.0, END], replicates=20,
                   seed=4, threads=1, method="mtdna_direct")
seen = Counter()


def watch(replicate, event):
    """Count events by reaction; stop a cell at its first advantaged mutant."""
    name = run.rule_names[event["rule"]]
    seen[name] += 1
    return name != "mutation_mtAdv"


result = run.simulate(callback=watch)
stopped = result.records["termination"] == "stopped"
print("events seen by the callback:", sum(seen.values()))
print("cells stopped at their first advantaged mutation:", stopped.sum())
days = result.records.loc[stopped, "stop_time"] / 24
print("their stop days:", days.round(1).to_numpy())
```

```text
events seen by the callback: 2675817
cells stopped at their first advantaged mutation: 5
their stop days: [739.2 672.8 574.3 938.7 230.9]
```

Both this run and the log above use the fast path with seed 4, so their
first 20 cells are the same cells. The advantaged mutation of cell 0 on day
739.2 in the log reappears here as the first stop day: an observer reads the
events, it does not change them.

## Single lineages

Applies: [chapter 14, New identities and random attributes](14-model-files.md#new-identities-and-random-attributes).

A class is not a lineage. Following every mutant lineage, with its birth,
its death and its largest size, needs every new mutant to carry its own
identity. That is a species with attributes, which the builder cannot
write. [Page 14a](14a-mtdna-model-files.md#every-lineage-with-the-clone-model)
runs the clone model, in which every lineage has its own identity.

See also: [Stop conditions, trackers and event logs](06-stops-trackers-logs.md),
[Methods and speed](12d-mtdna-methods.md).
