# Simulating mtDNA mutant accumulation: user manual

This manual is for researchers who simulate how mitochondrial DNA mutants
arise and spread in single cells. It shows how to do it with the rules
engine of CRN-craft, from Python, from the command line and from Rust. Every
Python block is a complete program: copy it and run it from the repository
root. Most are also in `examples/mtdna/`. The examples were run as printed.
Each one's output is shown below it, timings from a two-CPU Linux host.
`validation/check_manual.py` runs every block and compares its output with
the manual.

Companion documents:

- [Optimization tutorial](mtdna-optimization-tutorial.md): fitting model
  parameters to data with fcmaes.
- [Performance analysis](mtdna-performance-analysis.md): what each method
  and option costs, and how to size a study.
- [Rules language](rules-language.md) and [mtDNA ensembles](mtdna-rules.md):
  the reference pages behind this manual.

Contents:

1. [Installation](#1-installation)
2. [The models](#2-the-models)
3. [A first ensemble](#3-a-first-ensemble)
4. [Seeds, threads and reproducibility](#4-seeds-threads-and-reproducibility)
5. [Takeover, thresholds and stop conditions](#5-takeover-thresholds-and-stop-conditions)
6. [Event-level bookkeeping with trackers](#6-event-level-bookkeeping-with-trackers)
7. [Exact event records: logs and callbacks](#7-exact-event-records-logs-and-callbacks)
8. [Choosing a method](#8-choosing-a-method)
9. [Individual mutant lineages: the clone model](#9-individual-mutant-lineages-the-clone-model)
10. [Variable population: continuous variables](#10-variable-population-continuous-variables)
11. [Parameter sweeps](#11-parameter-sweeps)
12. [Building models in Python](#12-building-models-in-python)
13. [Model variants](#13-model-variants)
14. [Budgets, wall time and cancellation](#14-budgets-wall-time-and-cancellation)
15. [The command line and JSON requests](#15-the-command-line-and-json-requests)
16. [Rust](#16-rust)
17. [Checking results](#17-checking-results)
18. [Reference](#18-reference)
19. [Troubleshooting](#19-troubleshooting)

## 1. Installation

```sh
pip install crn-craft                    # Python package (abi3 wheels for CPython 3.10+)
pip install "crn-craft[xarray]"          # optional: results as xarray Datasets
cargo install crn-craft-cli --locked     # optional: the native command-line tool
```

The Python package needs NumPy and pandas, which pip installs with it;
pandas is a required dependency since 0.3.0, because results are pandas
DataFrames.

The Python package contains the native engine, so the CLI is optional. Every
feature in this manual is in release 0.4.0; features that exist only in
Rust are marked where they appear.

Both packages install a command named `crn-craft`. The rules commands
(`rules-check`, `rules`, `describe rules-run`) exist only in the native CLI
from `cargo install`; the Python package's `crn-craft` command covers the
fixed-network tools. If both are installed, the one that comes first on your
`PATH` wins. `crn-craft --version` tells them apart (the Python command
prints `"interface": "python"`); call the native binary by its full path,
by default `~/.cargo/bin/crn-craft`, if needed. Section 15
shows the rules commands.

To use the current source instead of the release, build the wheel from a
clone of the repository:

```sh
pip install maturin
maturin build --release --locked --out dist
pip install --force-reinstall dist/*.whl
```

Run the examples from the repository root. They read the model files under
`models/rules/`:

```sh
python examples/mtdna/01_first_run.py
```

The optimization tutorial also needs `pip install fcmaes scipy`.

## 2. The models

The engine simulates *rules models*: species types carry numeric
attributes, and rules act on every instance that matches a pattern. The
mtDNA models describe one cell. Its molecules degrade at a constant rate
and replicate into the free places below a copy-number target `K`.
Replication can produce a mutant. The repository has five model files:

| File | What it models | Methods |
| --- | --- | --- |
| `mtdna-three-class.rules.json` | fixed population target; wild type, neutral mutants and advantaged mutants as three classes | every generic method, and the opt-in fast path `mtdna_direct` |
| `mtdna-three-class-cme.rules.json` | the same with K = 20 and fast rates, for exact checks against the master equation | as above |
| `mtdna-three-class-atp-continuous.rules.json` | variable population: replication limited by a cellular ATP balance, written as a continuous variable (rules-model/0.3) | `hybrid` |
| `mtdna-three-class-atp.rules.json` | the same model in the stored rules-model/0.2 form, with a fixed `atp` block | `thinning` |
| `mtdna-fixed.rules.json` | every new mutant is its own clone with an identity | `tree` |

Section 12 shows how to write models like these from Python with
`rules.Model` instead of JSON.

### The three-class model

Each molecule is an `mtDNA` instance with three attributes: `clone` (0 for
wild type, 1 neutral, 2 advantaged), `s` (its replication weight) and `wt`
(1 for wild type, 0 for mutants). The parameters are:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `K` | 1000 | copy-number target |
| `delta` | ln 2 / 240 h⁻¹ | degradation rate per molecule (half-life 10 days) |
| `sw` | 1.0 | replication weight of wild type and neutral mutants |
| `sm` | 1.5 | replication weight of advantaged mutants |
| `mu` | 3·10⁻⁵ | mutation probability per wild-type replication |
| `fracAdv` | 0.2 | share of new mutants that are advantaged |
| `w0` | 1000 | initial wild-type molecules |

The observables are sums over the molecules:

| Observable | Value |
| --- | --- |
| `N` | all molecules |
| `wild_type` | molecules with `wt == 1` |
| `mutants` | molecules with `wt == 0` |
| `neutral_mutants` | `clone == 1` |
| `advantaged_mutants` | `clone == 2` |
| `D` | the competition denominator, Σ (`wt` + `sm` · (1 − `wt`)) over all molecules |

Every class degrades at rate `delta`. A molecule of weight `s` replicates
at rate `max(0, K − N) / D · s`. A wild-type replication yields:

- an unchanged copy with probability 1 − `mu`;
- a neutral mutant with probability `mu` · (1 − `fracAdv`);
- an advantaged mutant with probability `mu` · `fracAdv`.

The eight rules (`degradation_wt`, `replication_wt`, `mutation_mt`,
`mutation_mtAdv`, `degradation_mt`, `replication_mt`, `degradation_mtAdv`,
`replication_mtAdv`) write this out channel by channel. The model has no
stop condition. The examples add the takeover rule used throughout this
manual, `N > 0 && mutants >= 0.9 * N`: a cell is taken over when 90% of its
molecules are mutants.

Results are pandas DataFrames whose columns carry the observable names, so
select a column by name, as in `result.statistics["mean"]["N"]`. Observable
columns are sorted by name, not in the order of the file, and continuous
variables follow them (section 10); the raw matrices in `result.metadata`
use the same order.

## 3. A first ensemble

An *ensemble* is a set of independent replicates of the model. Here each
replicate is one cell. `rules.ensemble` builds and validates a run request
and prepares it; `simulate()` runs the cells on native threads and returns
pandas DataFrames.

Most examples pass `method="mtdna_direct"`, an opt-in fast path for the
three-class structure of this model (section 8). `method="static_direct"`
runs the same model, and any other closed model, about 3.3× slower; the
numbers then differ cell by cell but not in distribution.

```python
"""Three-class model: 200 cells for three years with the takeover stop."""

import json
import time

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
days = np.arange(0, 1081, 90)
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=24.0 * days,
    replicates=200,
    seed=1,
    threads=2,
    method="mtdna_direct",
)
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
print("takeover days (first five):", first.to_numpy() / 24)
mean = result.statistics["mean"]  # one row per report time
fraction = mean["mutants"] / mean["N"]
for hours, value in fraction.items():
    print(f"day {hours / 24:5.0f}: mean mutant fraction {value:.4f}")
```

```text
chacha8-u53-ziggurat/rules-mtdna-direct/1 completed
28,536,814 events in 0.53 s
taken over: 23 of 200 cells
          termination  stop_time  events
replicate                               
0           completed    25920.0  148726
1           completed    25920.0  149476
2           completed    25920.0  150111
takeover days (first five): [272.73509567 300.76054915 304.31896137 358.80591137 426.00462126]
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

What the result holds:

- `result.records`: a DataFrame with one row per cell, indexed by
  `replicate`, the cell's index. The columns are `started`, `termination`,
  `stop_time`, `events`, `rng_words`, `peak_instances`, `created_instances`
  and `saved_points`. `termination` is categorical: compare it with a name
  such as `"stopped"` (section 18 lists them).
- `result.rule_events`: events per cell and rule, one column per rule name.
- `result.statistics`: a dictionary of DataFrames `mean`, `sample_variance`
  and `mean_standard_error`, indexed by report `time` with one column per
  observable, and `samples`, the number of cells. They are withheld (`None`)
  if any cell ended early for a budget reason.
- `result.paths`: the kept series of single cells (next section).
- `result.metadata`: the full result document (sampler, termination,
  totals, statistics, kept paths, observations, diagnostics).

`.to_numpy()` turns any DataFrame or column into a NumPy array.

Report times are right-continuous: an event exactly at a report time is
included. A cell that stops holds its state at every later report time.

## 4. Seeds, threads and reproducibility

Every cell draws from its own ChaCha8 random stream, keyed by the root seed
and the cell's index. Results therefore depend on the seed and the method,
but not on the number of threads. A cell's trajectory is also the same in a
larger ensemble.

```python
"""Seeds, threads and reproducibility."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"


def run(seed, threads):
    ensemble = rules.ensemble(
        model,
        end=25_920.0,
        report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
        replicates=64,
        seed=seed,
        threads=threads,
        method="mtdna_direct",
        keep_paths="all",  # every cell's series, as result.paths
    )
    return ensemble.simulate()


one, two, other = run(7, 1), run(7, 2), run(8, 2)
same = one.records.equals(two.records) and one.rule_events.equals(two.rule_events)
print("1 vs 2 threads, same seed, identical records:", same)
print("same statistics:", one.statistics["mean"].equals(two.statistics["mean"]))
print("same paths:", one.paths.equals(two.paths))
print(
    "seed 7 vs 8, identical stop times:",
    one.records["stop_time"].equals(other.records["stop_time"]),
)
# Replicate i always uses the same stream, whatever the ensemble size.
small = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
    replicates=8,
    seed=7,
    threads=1,
    method="mtdna_direct",
).simulate()
print(
    "first 8 of 64 equal an 8-cell run:",
    small.records["events"].equals(one.records["events"].head(8)),
)
print(one.paths[["N", "mutants"]].head(5))
```

```text
1 vs 2 threads, same seed, identical records: True
same statistics: True
same paths: True
seed 7 vs 8, identical stop times: False
first 8 of 64 equal an 8-cell run: True
                        N  mutants
replicate time                    
0         0.0      1000.0      0.0
          8640.0    998.0      1.0
          17280.0   997.0      0.0
          25920.0   999.0      0.0
1         0.0      1000.0      0.0
```

`keep_paths` decides which cells keep their observable series: a count
keeps cells 0, 1, … up to it, and `"all"` (or `True`) keeps every cell.
`result.paths` holds them as one DataFrame indexed by (`replicate`,
`time`), with one column per observable. Some common forms:

- `result.paths.loc[3]`: cell 3's series;
- `result.paths.xs(8640.0, level="time")`: every kept cell at day 360;
- `result.paths["mutants"].unstack("time")`: cells × report times;
- `result.to_xarray()`: an xarray Dataset with dimensions `replicate` and
  `time`, one variable per observable (`pip install "crn-craft[xarray]"`).

A cell that ended on a budget has rows only up to its last saved report
time; `to_xarray()` fills the rest with NaN. Kept paths count against
`max_output_bytes`; the
[performance analysis](mtdna-performance-analysis.md#8-sizing-a-study)
gives their size.

Each method has its own sampler identity in `metadata["sampler"]`. The same
seed gives different trajectories under different methods. Results are
comparable across methods in distribution, not path by path.

## 5. Takeover, thresholds and stop conditions

A `stop` condition ends a cell the first time it is true after an event.
The cell's termination is `stopped` and `stop_time` is the event time.
The two common outputs are:

- **Takeover fraction by age:** the share of cells whose stop time lies
  before that age.
- **Time-to-threshold distribution:** the stop times of the stopped cells.

```python
"""Takeover fraction by age and time to the 90% threshold."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
result = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=200,
    seed=1,
    threads=2,
    method="mtdna_direct",
).simulate()
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

These are the 200 cells of section 3: the same seed gives the same cells
whatever the report times.

Any expression over parameters and observables can be a stop, for example
`mutants >= 0.6 * N` for a 60% COX-deficiency threshold, or
`advantaged_mutants > 0` to stop at the first advantaged mutant. Every
method accepts any stop. The fast path `mtdna_direct` keeps its speed only
with no stop or the guarded fraction `N > 0 && mutants >= f * N`; any other
stop makes it fall back to Static Direct with a warning (sections 8 and
13). To record a threshold without ending the cell, use a first-passage
tracker (next section).

## 6. Event-level bookkeeping with trackers

Trackers follow one observable at every event, at full native speed, and
report per cell:

- `appearances` and `extinctions`: how often the observable rose from 0, or
  fell to 0;
- `last_appearance` and `last_extinction`: when that last happened;
- `minimum` and `maximum`;
- `integral`: the time integral up to the cell's stop or horizon;
- `first_passage`: the first time an optional predicate was true;
- `episodes`: the first `keep_episodes` excursions above zero, each with
  start, end, peak and duration. Later ones are summarized as
  `omitted_episodes`, `omitted_duration` and `omitted_peak`.

```python
"""Event-level bookkeeping with trackers."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=500,
    seed=3,
    threads=2,
    method="mtdna_direct",
    trackers=[
        {"observable": "neutral_mutants", "keep_episodes": 3},
        {"observable": "advantaged_mutants", "keep_episodes": 3},
        # Time at which mutants first reach 60% of the molecules.
        {"observable": "mutants", "first_passage": "N > 0 && mutants >= 0.6 * N"},
    ],
)
result = run.simulate()
print("sampler:", result.metadata["sampler"])
# One DataFrame per tracker, one row per cell.
neutral = result.tracker_columns(0)
advantaged = result.tracker_columns(1)
passage = result.tracker_columns(2)
print(advantaged[["appearances", "extinctions", "maximum", "integral"]].head(3))
print("neutral lineages arising per cell:", neutral["appearances"].mean())
print("neutral lineages lost per cell:   ", neutral["extinctions"].mean())
print(
    "advantaged lineages arising/lost: ",
    advantaged["appearances"].mean(),
    advantaged["extinctions"].mean(),
)
print("largest advantaged class seen:    ", advantaged["maximum"].max())
# The integral runs until the replicate's stop or the horizon.
load = advantaged["integral"] / result.records["stop_time"]
print("time-averaged advantaged count, mean over cells:", load.mean())
reached = passage["first_passage"].dropna()  # NaN: never reached 60%
print(
    f"cells reaching 60% mutants: {len(reached)}; median day "
    f"{reached.median() / 24:.0f}"
)
# Episodes: each excursion of a class above zero, with start, end, peak.
first = next(o for o in result.observations if o["trackers"][1]["episodes"])
for episode in first["trackers"][1]["episodes"]:
    print("cell", first["id"], "advantaged episode:", episode)
```

```text
sampler: chacha8-u53-ziggurat/rules-mtdna-direct/1
           appearances  extinctions  maximum   integral
replicate                                              
0                    0            0      0.0   0.000000
1                    1            1      1.0  64.109016
2                    0            0      0.0   0.000000
neutral lineages arising per cell: 1.566
neutral lineages lost per cell:    1.472
advantaged lineages arising/lost:  0.404 0.266
largest advantaged class seen:     900.0
time-averaged advantaged count, mean over cells: 17.689950061057512
cells reaching 60% mutants: 52; median day 494
cell 1 advantaged episode: {'duration': 64.10901603704042, 'end': 5246.85151399302, 'peak': 1.0, 'start': 5182.742497955979}
```

In the three-class model a class appearance is the birth of a mutant
lineage *when its class is absent*. A new mutant that arises while its class
is present merges into it. For the fate of every single lineage, use the
clone model (section 9).

The integral divided by the observed time is the time-averaged count. For
`advantaged_mutants` above, it is the average number of advantaged
molecules per hour lived.

For plain count trackers (an unweighted class count without a
first-passage predicate), `"static_direct"` and `mtdna_direct` integrate
sparsely: an integral is advanced only when its class changes. The last
bits of an integral can then differ from per-event integration.
`{"kind": "static_direct", "sparse_trackers": False}` integrates at every
event, as Static Direct did in release 0.3.0.

## 7. Exact event records: logs and callbacks

An *event log* keeps, for every cell, up to `max_entries` events of the
listed rules (all rules if the list is empty). Each entry holds the time,
the rule, the participating instances and the observables after the event.

```python
"""Exact event records: a bounded event log for every cell."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=100,
    seed=4,
    threads=2,
    method="static_direct",
    event_log={"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]},
)
result = run.simulate()
# Every logged event, indexed by (replicate, entry); rule_name names the rule.
log = result.event_columns()
counts = log.groupby("rule_name", observed=True).size()
print("logged mutation events:", counts.to_dict())
print("omitted (log full):", sum(o["omitted_log_entries"] for o in result.observations))
cell = log.loc[0]  # the log of cell 0
print(cell[["time", "rule", "rule_name", "N", "mutants"]])
print("cell 0 mutation days:", (cell["time"] / 24).round(1).tolist())
```

```text
logged mutation events: {'mutation_mt': 166, 'mutation_mtAdv': 39}
omitted (log full): 0
               time  rule       rule_name      N  mutants
entry                                                    
0       1442.924093     2     mutation_mt  995.0      1.0
1       7745.901298     2     mutation_mt  997.0      1.0
2      17740.229245     3  mutation_mtAdv  998.0      1.0
cell 0 mutation days: [60.1, 322.7, 739.2]
```

A *live callback* is a Python function called for every committed event of
a cell, in order. Returning `False` stops that cell. Each event carries the
rule index; `run.rule_names[event["rule"]]` gives the rule's name:

```python
"""Exact event records: a live callback that stops cells."""

import json
from collections import Counter

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=20,
    seed=4,
    threads=1,
    method="static_direct",
)
seen = Counter()


def watch(replicate, event):
    """Count events by rule; stop a cell at its first advantaged mutant."""
    name = run.rule_names[event["rule"]]
    seen[name] += 1
    return name != "mutation_mtAdv"


# Callbacks hold the GIL, so use them on few cells.
result = run.simulate(callback=watch)
stopped = result.records["termination"] == "stopped"
print("events seen by the callback:", sum(seen.values()))
print("cells stopped at their first advantaged mutation:", stopped.sum())
days = result.records.loc[stopped, "stop_time"] / 24
print("their stop days:", days.round(1).to_numpy())
```

```text
events seen by the callback: 2770371
cells stopped at their first advantaged mutation: 5
their stop days: [739.2 672.8 574.3 938.7 914.1]
```

Both runs use Static Direct with seed 4. Cell 0's advantaged mutation in
the log, on day 739.2, reappears as its stop day in the callback run, which
stopped it at that first advantaged mutation. An observer reads the same
events; it does not change them.

Rules for observation:

- Observers read state and may stop a replicate. They never change the
  state or draw random numbers, so attaching one leaves results unchanged.
- A callback receives each event as a parsed JSON document and holds
  Python's interpreter lock. It cost about 14 µs per event on the test host,
  a hundred times the simulation itself, and it serializes the threads. Use
  it on a few cells, and use trackers or logs for whole ensembles.
- Log retention counts against the admission budget. Large logs need small
  cohorts and a larger `max_output_bytes`.
- `batch_callback=fn` replays the retained log after the run in chunks.
- `result.event_columns()` returns every logged event as one DataFrame
  indexed by (`replicate`, `entry`): `time`, `rule` (the index),
  `rule_name`, then the observables after the event.
  `result.event_columns(id)` returns one cell's log. Participating
  instances are in `result.observations`.

## 8. Choosing a method

The generic methods run every model within their scope. None of them
depends on mtDNA or on the names in a model:

| `method` | Kind | Scope | Use it for |
| --- | --- | --- | --- |
| `"static_direct"` | exact | closed models (a fixed set of attribute vectors) | the default for most work; Ziggurat waiting times and sparse count trackers (rules-run/0.6) |
| `{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False}` | exact | closed models | the Static Direct of release 0.3.0: inverse waiting times and per-event trackers; reproduces its numbers |
| `"next_reaction"` | exact | closed models | Gibson–Bruck; pays off in large sparse networks, not here, where every event changes every rate |
| `None` (tree) | exact | every model without an ATP block or continuous variables | models whose instances change identity, such as the clone model |
| `{"kind": "tau", "epsilon": 0.03}` | approximate | closed models | adaptive Poisson leaps; at K = 1000 the vacancy limits leaps, so it runs exact events |
| `{"kind": "binomial", "step": 1.0}` | approximate | closed models with single-parent unit births and deaths | fixed hourly splitting, similar to a binomial hourly scheme in Java |
| `"hybrid"` or `{"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}` | approximate at its tolerances | closed models, with or without continuous variables | models with continuous variables (section 10), which no other method runs; without them it runs exact events up to rounding, slower than Static Direct |

Two methods serve one model structure each:

| `method` | Kind | Use it for |
| --- | --- | --- |
| `"mtdna_direct"` | exact | an opt-in fast path for the three-class fixed-population structure, about 3.3× faster than Static Direct; other models fall back to Static Direct with inverse waits (below) |
| `"thinning"` | exact | stored models with the legacy ATP block (section 10), and only those |

The same 64 cells under each method that runs the three-class model:

```python
"""The same 64 cells under every method that runs the three-class model."""

import json
import time

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
methods = {
    "tree": None,
    "static_direct": "static_direct",
    "static_direct, inverse": {
        "kind": "static_direct",
        "waiting_time": "inverse",
        "sparse_trackers": False,
    },
    "next_reaction": "next_reaction",
    "tau (epsilon 0.03)": {"kind": "tau", "epsilon": 0.03},
    "binomial (1 h steps)": {"kind": "binomial", "step": 1.0},
    "hybrid": "hybrid",
    "mtdna_direct (opt-in)": "mtdna_direct",
}
for label, method in methods.items():
    run = rules.ensemble(
        model,
        end=25_920.0,
        report_times=[0.0, 12_960.0, 25_920.0],
        replicates=64,
        seed=5,
        threads=2,
        method=method,
    )
    start = time.perf_counter()
    result = run.simulate()
    seconds = time.perf_counter() - start
    meta = result.metadata
    stopped = result.records["termination"] == "stopped"
    final = result.statistics["mean"].iloc[-1]  # the means at the last report time
    fraction = final["mutants"] / final["N"]
    print(
        f"{label:24s} {seconds:6.2f} s  {meta['events']:>11,} events  "
        f"taken over {stopped.sum():2d}  final mutant fraction {fraction:.3f}  "
        f"{meta['sampler']}"
    )
```


```text
tree                       1.56 s    9,167,872 events  taken over  6  final mutant fraction 0.097  chacha8-u53/rules-direct/2
static_direct              0.37 s    8,861,769 events  taken over 10  final mutant fraction 0.164  chacha8-u53-ziggurat/rules-static-direct/1
static_direct, inverse     0.40 s    9,180,176 events  taken over  8  final mutant fraction 0.144  chacha8-u53/rules-static-direct/1
next_reaction              0.62 s    9,226,649 events  taken over  7  final mutant fraction 0.105  chacha8-u53/rules-next-reaction/1
tau (epsilon 0.03)         2.34 s    9,225,966 events  taken over  5  final mutant fraction 0.108  chacha8-u53/rules-approximate-cgp-poisson/2
binomial (1 h steps)       0.39 s    9,158,734 events  taken over  6  final mutant fraction 0.088  chacha8-u53/rules-approximate-binomial-split/2
hybrid                     1.56 s    9,181,653 events  taken over  8  final mutant fraction 0.144  chacha8-u53/rules-approximate-hybrid-dopri54/1
mtdna_direct (opt-in)      0.12 s    8,969,870 events  taken over  8  final mutant fraction 0.122  chacha8-u53-ziggurat/rules-mtdna-direct/1
```

The takeover counts differ because each method uses the random stream
differently. With 64 cells, 5 to 10 takeovers are within sampling error of
each other. `"static_direct"` draws Ziggurat waiting times;
`static_direct, inverse` is the Static Direct of release 0.3.0 and gives
its numbers. The approximate methods are statistically close here but are
labeled as approximations in their sampler identity. Without continuous
variables, `hybrid` uses the random stream as Static Direct with inverse
waits does: the same cells take over at nearly the same times, and its
event times differ only by rounding, so paths drift apart over many
events. It costs more per
event, since it integrates between events. The
[performance analysis](mtdna-performance-analysis.md) compares the methods
at full scale.

### The opt-in fast path `mtdna_direct`

`mtdna_direct` is an exact direct method fused for one model structure:
three classes in a fixed population, as in section 2. It reads that
structure from the lowered rules; the model file and its names do not
matter. At preparation it lowers the model to a fixed network, as Static
Direct does, and certifies that:

1. **Three classes, eight channels.** The network has exactly three
   attribute vectors (classes) and eight channels. Each channel consumes
   one molecule and changes one class by one molecule: every class
   degrades and replicates, and one class, the wild type, also yields one
   molecule of each other class (the two mutation channels). The roles are
   read from the channels.
2. **Degradation.** All three classes degrade at one constant rate.
3. **Replication.** Every replication and mutation rate is
   `max(0, K - N) / D` times constant factors. All five share one `K`, an
   expression of parameters whose value is an integer from 1 to 1,000,000
   and at least the initial population. `N` is an observable that counts
   every molecule, and `D` an observable that sums the classes with finite
   positive weights.
4. **Observables.** Every observable is a finite weighted sum over the
   three classes.
5. **Stop.** The stop is absent, or `N > 0 && mutants >= f * N` with a
   number `f` in [0, 1], where `mutants` counts the two classes other than
   the wild type.
6. **Resolution.** Bounds from `K`, the weights and the rate constants
   show that every hazard stays finite, and every positive hazard at least
   16 · 2⁻⁵³ of the total, in every reachable state.

The certificate holds for one model and one parameter set:
`bind_parameters` certifies each candidate again. If a check fails, the
run uses generic Static Direct with inverse waits and per-event trackers.
`metadata["sampler"]` is then `chacha8-u53/rules-static-direct/1`, the
numbers equal those of `{"kind": "static_direct", "waiting_time":
"inverse", "sparse_trackers": False}`, and a `CRN-R051` warning in
`metadata["diagnostics"]` names the failed check. A model that cannot be
lowered at all, such as the clone model, is rejected with `CRN-R020`, as
by every method except the tree.

`mtdna_direct` samples the same process as Static Direct with a different
random stream: each cell's ChaCha8 words go to Ziggurat waiting times and
to a selection over grouped hazards. Its sampler identity is
`chacha8-u53-ziggurat/rules-mtdna-direct/1`. Results therefore agree in
distribution, not cell by cell; section 17 lists the checks. Trackers,
event logs and callbacks work as under Static Direct, and count trackers
integrate sparsely (section 6).

What it gains: it computes the eight hazards from three counts with one
division and reads its stop from the counts. On the three-class model
that is about 3.3× less time per event than Static Direct with its
default Ziggurat waits and sparse trackers, and 3.9× less than Static
Direct with inverse waits
([performance analysis](mtdna-performance-analysis.md#2-methods)).

The structure need not be called mtDNA. Below, three strains share `room`
places, and the resident strain gives rise to two variants. Every name
differs from the model file, and `mtdna_direct` still certifies the model.
A strain that dies faster breaks check 2:

```python
"""mtdna_direct recognizes a model structure, not a file or names."""

import copy
import time

from crn_craft import rules

# Three strains share `room` places. The resident strain (0) gives rise to
# the two variants; every strain dies at one rate and grows into free room.
grow = "max(0, room - total) / crowding"
model = {
    "schema": "crn-craft.rules-model/0.1",
    "name": "three_strains",
    "parameters": {"room": 200, "loss": 0.01, "g": 0.5, "eps": 0.0002},
    "types": [{"name": "Cell", "attributes": ["strain"]}],
    "initial": [{"type": "Cell", "attributes": {"strain": 0}, "count": 200}],
    "observables": {
        "total": {"type": "Cell"},
        "crowding": {"type": "Cell", "weight": "1 + (x.strain == 2)"},
        "variants": {"type": "Cell", "where": "x.strain > 0"},
    },
    "stop": "total > 0 && variants >= 0.5 * total",
    "rules": [],
}
for strain, rate in ((0, "g"), (1, "g"), (2, "2 * g")):
    cell = [{"label": "c", "type": "Cell", "where": f"c.strain == {strain}"}]
    model["rules"] += [
        {"name": f"die_{strain}", "reactants": cell, "rate": "loss", "products": []},
        {
            "name": f"grow_{strain}",
            "reactants": cell,
            "rate": f"{grow} * {rate}",
            "products": [{"type": "Cell", "copy": "c"}, {"type": "Cell", "copy": "c"}],
        },
    ]
for strain in (1, 2):
    model["rules"].append(
        {
            "name": f"switch_to_{strain}",
            "reactants": [{"label": "c", "type": "Cell", "where": "c.strain == 0"}],
            "rate": f"{grow} * eps",
            "products": [
                {"type": "Cell", "copy": "c"},
                {"type": "Cell", "attributes": {"strain": strain}},
            ],
        }
    )


def run(model, method):
    start = time.perf_counter()
    result = rules.ensemble(
        model,
        end=2000.0,
        report_times=[0.0, 2000.0],
        replicates=4000,
        seed=1,
        threads=2,
        method=method,
    ).simulate()
    seconds = time.perf_counter() - start
    taken = (result.records["termination"] == "stopped").mean()
    print(
        f"{method:14s} {seconds:5.2f} s  variants took over: {taken:.3f}  "
        f"{result.metadata['sampler']}"
    )
    return result


run(model, "static_direct")
run(model, "mtdna_direct")
# A strain that dies faster breaks the shared degradation rate, so the run
# falls back to Static Direct and says why.
faster = copy.deepcopy(model)
next(r for r in faster["rules"] if r["name"] == "die_1")["rate"] = "2 * loss"
result = run(faster, "mtdna_direct")
print([d["message"] for d in result.metadata["diagnostics"]])
```

```text
static_direct   1.09 s  variants took over: 0.421  chacha8-u53-ziggurat/rules-static-direct/1
mtdna_direct    0.34 s  variants took over: 0.427  chacha8-u53-ziggurat/rules-mtdna-direct/1
mtdna_direct    1.11 s  variants took over: 0.434  chacha8-u53/rules-static-direct/1
['mtdna_direct ran Static Direct: degradation is not one constant rate shared by all classes']
```

The first two runs agree within sampling error: with 4000 cells each, the
standard error of their difference is about 0.011. The third ran Static
Direct, at its speed.

### Optional: compiled event loop

*Opt-in, off by default; needs a Rust toolchain. New in 0.4.0.*

`rules.compile_loop(model)` compiles the whole Static Direct event of one
model into a small native program; `compiled=` then runs Static Direct
requests in it. Results are byte-identical to the interpreter's: the same
records, statistics, paths, trackers, event logs and sampler identity, so
replays and stored results keep their meaning. Unlike `mtdna_direct`, it
keeps Static Direct's use of the random stream and is not limited to the
three-class structure.

```python
"""Static Direct with a compiled event loop (needs cargo)."""

import json
import time

from crn_craft import CrnError, rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
loop = rules.compile_loop(model)  # about a minute the first time, then cached
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 12_960.0, 25_920.0],
    replicates=200,
    seed=5,
    threads=2,
    method="static_direct",
    trackers=[{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}],
)
for label, compiled in (("interpreted", None), ("compiled", loop)):
    start = time.perf_counter()
    result = run.simulate(compiled=compiled)
    seconds = time.perf_counter() - start
    stopped = (result.records["termination"] == "stopped").sum()
    print(
        f"{label:12s} {seconds:5.2f} s  {result.metadata['events']:>11,} events  "
        f"taken over {stopped:2d}  {result.metadata['sampler']}"
    )
# Runs the request both ways; raises unless the results are byte-identical.
loop.verify(run)
print("byte-identical")

# Without cargo: does the generator support a model?
other = json.loads(json.dumps(model))
other["rules"][0]["rate"] = "delta * (1 + mutants / K)"
try:
    rules.loop_source(other)
except CrnError as error:
    print(error.diagnostics[0]["code"], error.diagnostics[0]["message"])
```

```text
interpreted   1.19 s   28,523,982 events  taken over 20  chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1
compiled      0.51 s   28,523,982 events  taken over 20  chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1
byte-identical
CRN-R020 generated rules loops support constant or vacancy rates and sum observables
```

On the 1000-cell workload of the
[performance analysis](mtdna-performance-analysis.md#2-methods), timed
from Python including the subprocess, the compiled loop took 2.6 s
against 5.8 s for interpreted Static Direct with its defaults (2.21×),
and 3.8 s against 7.0 s with inverse waits and per-event trackers (1.86×),
medians of three on the two-CPU host
(`validation/mtdna/compiled-loop-python.json`).
`mtdna_direct` is faster still where it applies.

**Requirements.** Rust 1.85 or newer, the MSRV of the crates (install it
with [rustup](https://rustup.rs)), and network access to crates.io: the
program depends on `crn-craft-rules` at exactly the installed version.
Offline, or with a version that is not on crates.io (a wheel built from a
clone), pass `source=` a CRN-craft checkout of that version, or set
`CRN_CRAFT_SOURCE`. `cargo=` names another cargo executable.

**Build time and cache.** The first build compiles `crn-craft-rules` and
its dependencies: about 40 s on the two-CPU host after the downloads, several
minutes on slower machines. `compile_loop` prints cargo's progress and the
path of the build log to stderr (`progress=False` silences it). Builds are
cached in the user cache directory (`~/.cache/crn-craft/rules-loops` on
Linux; `cache_dir=` or `CRN_CRAFT_CACHE_DIR` choose another), keyed by the
loop's fingerprint, the installed version, `rustc -vV`, the target and
Rust flags. A cached call took 0.06 s. All builds share one Cargo target
directory (`CARGO_TARGET_DIR` if set), so a second model built in 3 s.
The cache takes about 160 MB; delete the directory to reclaim it.

**Limitations.**

- Only `method="static_direct"`, with any options. Other methods are
  refused with `CRN-R020` at `/method`.
- The generator needs constant or vacancy rates (`max(0, K - N) / D` times
  constants) and weighted-sum observables, and at most 1 MiB of source.
  It refuses expression rates, `propensity` rules, max or distinct
  observables, continuous variables, ATP blocks and models whose instances
  change identity with `CRN-R020` at `/method`, before cargo runs.
  `rules.loop_source(model)` returns the source, or that error, without
  cargo.
- Parameters are compiled in, so another parameter set needs another loop.
  A loop passed for other model bytes, or for a `bind_parameters`
  candidate, is refused with `CRN-R020` at `/model`; candidates run
  interpreted.
- `callback` and `batch_callback` are refused with `ValueError`: the loop
  runs in its own process. Trackers and event logs work.
- Threads, report times, stops, budgets, `wall_time_s` and `cancel_file`
  work as interpreted. Ctrl-C stops the program and raises
  `KeyboardInterrupt`.
- On networks with many constant-rate channels the loop can be slower than
  the interpreter
  ([details](specialized-kernels.md#generated-rules-event-loops)).
- Missing cargo, a failed build (`error.log` is the build log) and a
  result that differs raise `rules.CompiledLoopError`, whose `reason` is
  `"cargo_missing"`, `"build_failed"` or `"mismatch"`.

**Confirming identity.** Every `compile_loop` call, cached or not, runs
short requests through the program and the interpreter, with both waiting
times, both tracker policies and an event log, and returns the loop only if
the metadata and every record column are byte-identical. `loop.verify(run)`
does the same for your own request, as above, and returns the result. On a
mismatch both results are saved next to the build for inspection.

## 9. Individual mutant lineages: the clone model

`mtdna-fixed.rules.json` gives every new mutant its own `clone` identity
(`new_id()`) and draws its weight (`choice(fracAdv, sm, sw)`). The number of
instances therefore changes during a run, and only the `tree` method can
simulate it. The observable `mutant_clones` counts distinct mutant clones.
The full event log then gives every lineage's birth, death and peak size:

```python
"""Per-lineage fates with the clone model (tree method)."""

import json

import numpy as np
from crn_craft import rules

# Every new mutant is its own clone (attribute `clone` from new_id()).
model = json.loads(rules.load("models/rules/mtdna-fixed.rules.json"))
model["parameters"].update({"K": 100, "w0": 100, "mu": 0.01, "delta": 0.05})
run = rules.ensemble(
    model,
    end=400.0,
    report_times=[0.0, 100.0, 200.0, 400.0],
    replicates=20,
    seed=6,
    threads=2,
    keep_instances=3,
    trackers=[{"observable": "mutant_clones"}],
    event_log={"max_entries": 5000, "rules": []},
    max_output_bytes=1 << 30,
)
result = run.simulate()
print("sampler:", result.metadata["sampler"])
mean = result.statistics["mean"]
print("mean number of distinct mutant clones:", mean["mutant_clones"].to_numpy())
print(
    "peak instances (clones + wild type), max over cells:",
    result.records["peak_instances"].max(),
)
print("instances created per cell:", result.records["created_instances"].mean())
omitted = sum(o["omitted_log_entries"] for o in result.observations)
print("log entries omitted:", omitted)  # 0: every event is in the log

# Lineage lifetimes from the full event log: replay each cell's net changes
# per clone id (reactants are consumed, products created).
born, died, peak = {}, {}, {}
for observation in result.observations:
    counts = {}
    for entry in observation["log"]:
        change = {}
        for sign, records in ((-1, entry["reactants"]), (1, entry["products"])):
            for record in records:
                clone, _, wt = record["attributes"]
                if wt == 0:  # mutant clones only
                    key = (observation["id"], clone)
                    change[key] = change.get(key, 0) + sign * record["count"]
        for key, delta in change.items():
            counts[key] = counts.get(key, 0) + delta
            born.setdefault(key, entry["time"])
            peak[key] = max(peak.get(key, 0), counts[key])
            if counts[key] == 0:
                died[key] = entry["time"]
life = np.array([died[k] - born[k] for k in born if k in died])
print(f"{len(born)} clones born, {len(life)} lost within the horizon")
print("median lifetime of lost clones (h):", np.median(life))
print("largest clone ever (molecules):", max(peak.values()))
print("final instances of cell 0:", result.metadata["final_instances"][0])
```

```text
sampler: chacha8-u53/rules-direct/2
mean number of distinct mutant clones: [0.   1.9  2.6  1.85]
peak instances (clones + wild type), max over cells: 7
instances created per cell: 10.1
log entries omitted: 0
202 clones born, 165 lost within the horizon
median lifetime of lost clones (h): 18.433420576806014
largest clone ever (molecules): 100
final instances of cell 0: {'id': 0, 'instances': [{'attributes': [0.0, 1.0, 1.0], 'count': 7, 'type_index': 0}, {'attributes': [4.0, 1.5, 0.0], 'count': 17, 'type_index': 0}, {'attributes': [5.0, 1.5, 0.0], 'count': 68, 'type_index': 0}, {'attributes': [8.0, 1.0, 0.0], 'count': 1, 'type_index': 0}]}
```

The example uses a small cell (K = 100) and a high mutation probability so
that every event fits in the log. For K = 1000 over three years, a cell has
about 150,000 events. Track lineages with `mutant_clones`, `peak_instances`
and `created_instances`, and log only the mutation and degradation rules of
a few cells. Per cell, `created_instances` counts new instances, here one
per lineage.

## 10. Variable population: continuous variables

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
run = rules.ensemble(
    data,
    end=25_920.0,
    report_times=24.0 * np.arange(0, 1081, 180),
    replicates=64,
    seed=7,
    threads=2,
    method="hybrid",  # rtol 1e-6, atol 1e-9 unless a method dict says otherwise
)
result = run.simulate()
meta = result.metadata
print("columns:", list(result.statistics["mean"].columns))
print("continuous columns:", meta["continuous"])
print("sampler:", meta["sampler"], "method:", meta["method"])
stopped = result.records["termination"] == "stopped"
print("stopped (ATP exhausted):", stopped.sum())
mean = result.statistics["mean"]
for hours, row in mean.iterrows():
    print(
        f"day {hours / 24:4.0f}: N {row['N']:7.1f}  "
        f"mutants {row['mutants']:6.1f}  ATP {row['ATP']:7.2f}"
    )
```

```text
continuous: [{'name': 'ATP', 'initial': 'ATP0', 'derivative': 'f * wild_type - v1 * ATP - v2 * N', 'lower_bound': 0, 'on_bound': 'stop'}]
replication_wt rate: (sw * (1 - mu)) * c / (c + ATP)
columns: ['N', 'advantaged_mutants', 'mutants', 'neutral_mutants', 'wild_type', 'ATP']
continuous columns: ['ATP']
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
integration step. The error has not yet been measured as a function of the
tolerances.

**Your own variables.** `rules.Model.add_continuous(name, initial,
derivative, lower_bound=None, on_bound="stop")` adds a variable to a builder
model (section 12); string rates may read it:

```python
"""A continuous variable in a builder model."""

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

Rules for continuous variables:

- Rule rates and propensities may read them; the model's `stop`, observable
  filters and weights, reactant patterns and products may not (`CRN-R012`).
  Stop on a variable with its `lower_bound`.
- Trackers follow discrete observables only (`CRN-R020` for a continuous
  one). Event logs and callbacks see the variables' values at each event.
- `hybrid` lowers the model to fixed tuples like `static_direct`, so models
  whose instances change identity (the clone model) need `tree`, which has
  no continuous variables.
- A bound is checked on the integrator's dense output at each step. A
  variable that dips below its bound and returns within one step can be
  missed; tighter tolerances shorten the steps.

### Stored models: the ATP block and `thinning`

`mtdna-three-class-atp.rules.json` (rules-model/0.2) writes the same model
with a fixed `atp` block, which `method="thinning"` samples exactly against
a bound on each interval. Stored files and requests of this form keep
working unchanged and replay byte-identically. A model cannot have both an
`atp` block and `continuous` variables; for new models, use continuous
variables.

```python
"""A stored rules-model/0.2 file with the legacy ATP block (exact thinning)."""

import json

from crn_craft import rules

data = rules.load("models/rules/mtdna-three-class-atp.rules.json")
print("ATP block:", json.loads(data)["atp"])
result = rules.ensemble(
    data,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=64,
    seed=7,
    threads=2,
    method="thinning",
).simulate()
print("sampler:", result.metadata["sampler"])
print("stopped (ATP exhausted):", (result.records["termination"] == "stopped").sum())
print(result.statistics["mean"][["N", "mutants", "ATP"]].round(2))
```

```text
ATP block: {'wild_type': 'wild_type', 'total': 'N', 'initial': 450.0, 'production': 0.1, 'decay': 0.2, 'maintenance': 0.01, 'scale': 1.3034153747670567, 'replication_rules': ['replication_wt', 'mutation_mt', 'mutation_mtAdv', 'replication_mt', 'replication_mtAdv']}
sampler: chacha8-u53/rules-atp-thinning/1
stopped (ATP exhausted): 8
               N  mutants     ATP
time                             
0.0      1000.00     0.00  450.00
25920.0  5919.92  4587.14  394.85
```

The final means of the two 64-cell runs differ (`N` 3493 against 5920; 4
and 8 cells stopped). Mutant expansion is rare and lets the population
grow, so the final `N` is heavy-tailed, and 64 cells are too few to compare
methods. `validation/check_hybrid.py` compares the two methods on
400 cells each, and optionally 4,000, with predeclared two-sample tests.

## 11. Parameter sweeps

`template.bind_parameters(values)` prepares a candidate with new parameter
values and exactly the same request settings. Each candidate is validated
and lowered once, in about 0.3 ms.

```python
"""A parameter sweep that reuses one prepared template."""

import json
import time

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
template = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=400,
    seed=8,
    threads=2,
    method="mtdna_direct",
)
start = time.perf_counter()
print("takeover fraction after three years")
print("   mu \\ sm " + "".join(f"{sm:8.2f}" for sm in (1.0, 1.25, 1.5, 2.0)))
for mu in (1e-5, 3e-5, 1e-4):
    row = []
    for sm in (1.0, 1.25, 1.5, 2.0):
        result = template.bind_parameters({"mu": mu, "sm": sm}).simulate()
        row.append((result.records["termination"] == "stopped").mean())
    print(f"{mu:9.0e}  " + "".join(f"{v:8.3f}" for v in row))
print(f"12 candidates in {time.perf_counter() - start:.1f} s")
# The same seed gives every candidate the same random streams (common
# random numbers), so differences between cells of the table are not noise
# from different seeds alone.
```

```text
takeover fraction after three years
   mu \ sm     1.00    1.25    1.50    2.00
    1e-05     0.000   0.018   0.035   0.072
    3e-05     0.000   0.045   0.098   0.190
    1e-04     0.000   0.142   0.292   0.482
12 candidates in 14.5 s
```

With equal weights (`sm = 1`), no cell reached 90% mutants within three
years at K = 1000: drift alone is too slow. The takeover fraction grows
roughly in proportion to `mu` and steeply with `sm`. The
[optimization tutorial](mtdna-optimization-tutorial.md) turns this into a
fit. A candidate's `model_sha256` and `run_sha256` hash the bound model and
request, so every result names its parameters.

## 12. Building models in Python

`rules.Model` writes a rules model from Python calls, in the style of
rebop's `Gillespie` class. Every species is an attribute-free type, and
every species count is an observable with the species name.

```python
"""An SIR epidemic written with rules.Model."""

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
print(result.statistics["mean"].round(1))  # species columns, sorted by name
print("events per cell:", result.rule_events.mean().to_dict())
```

```text
Model('sir': 3 species, 2 reactions)
         I      R      S
time                    
0.0   10.0    0.0  490.0
20.0  90.9   86.1  323.0
40.0  87.4  295.0  117.6
60.0  28.7  404.4   66.9
events per cell: {'infection': 423.065, 'recovery': 404.412}
```

The calls:

| Call | Effect |
| --- | --- |
| `rules.Model(name="model")` | a new rules-model/0.3 document |
| `add_parameter(name, value)` | a parameter that string rates may read |
| `add_species(name, initial=0)` | a species type, its count observable and its initial count |
| `set_initial(mapping)` or `set_initial(S=490, I=10)` | initial counts of declared species |
| `add_reaction(rate, reactants, products, name=None)` | a rule; species are lists of names, repeated for multiplicity; undeclared species start at 0; names default to `r0`, `r1`, ... |
| `add_continuous(name, initial, derivative, lower_bound=None, on_bound="stop")` | a continuous variable (section 10) |
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
parameters = {"K": 1000, "delta": 0.0028881132523331052, "sw": 1.0, "sm": 1.5}
parameters |= {"mu": 3e-05, "fracAdv": 0.2}
for name, value in parameters.items():
    model.add_parameter(name, value)
model.add_species("W", 1000).add_species("M").add_species("A")
vacancy = "max(0, K - W - M - A) / (W + sm * (M + A))"
model.add_reaction(f"{vacancy} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt")
model.add_reaction(f"{vacancy} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt")
model.add_reaction(f"{vacancy} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv")
model.add_reaction(f"{vacancy} * sw * M", ["M"], ["M", "M"], "replication_mt")
model.add_reaction(f"{vacancy} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
for species in ("W", "M", "A"):
    model.add_reaction(parameters["delta"], [species], [], f"degradation_{species}")
result = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=100,
    seed=1,
    threads=2,
    method="mtdna_direct",
).simulate()
print(result.metadata["sampler"])
print(result.metadata["diagnostics"][0]["message"])
final = result.statistics["mean"].loc[25_920.0]
print(f"mean mutant fraction after three years: {(final['M'] + final['A']) / final.sum():.4f}")
```

```text
chacha8-u53/rules-static-direct/1
mtdna_direct ran Static Direct: the model has whole-rule propensities
mean mutant fraction after three years: 0.0828
```

`mtdna_direct` certifies only the JSON form of the model; for builder
models it runs Static Direct and says so in a `CRN-R051` warning. Species
with attributes (clone identities, drawn fitness) need the JSON language
(section 2 and the [rules language](rules-language.md)).

## 13. Model variants

The model file is plain JSON. Load it with `json.loads(rules.load(path))`,
change it, and pass the dictionary to `rules.ensemble`.

- **Copy number and half-life:** set `K`, `w0` and `delta = ln 2 / half-life`.
- **Initial heteroplasmy:** add an `initial` entry for the mutant class.
- **Threshold:** change the fraction in `stop`.
- **Recognition:** see what happens when the fast path `mtdna_direct`
  cannot certify the model.
- **Errors:** see how mistakes are reported.

```python
"""Model variants: copy number, half-life, initial heteroplasmy, threshold."""

import json
import math

from crn_craft import CrnError, rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["parameters"].update({"K": 500, "w0": 400, "delta": math.log(2) / (20 * 24)})
# Start heteroplasmic: 400 wild type and 100 neutral mutants (class 1).
model["initial"].append(
    {"type": "mtDNA", "attributes": {"clone": 1, "s": "sw", "wt": 0}, "count": 100}
)
model["stop"] = "N > 0 && mutants >= 0.6 * N"  # 60% threshold
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
    replicates=300,
    seed=9,
    threads=2,
    method="mtdna_direct",
)
result = run.simulate()
stopped = result.records["termination"] == "stopped"
print(result.metadata["sampler"])
print(f"cells past 60% within three years: {stopped.mean():.2%}")

# mtdna_direct certifies the model at preparation. A model it does not
# recognize runs Static Direct and says why in a CRN-R051 warning.
odd = json.loads(json.dumps(model))
odd["stop"] = "mutants > 250"
result = rules.ensemble(
    odd, end=100.0, report_times=[0.0, 100.0], replicates=2, method="mtdna_direct"
).simulate()
print(result.metadata["sampler"])
print([d["message"] for d in result.metadata["diagnostics"]])

# Mistakes are reported with a code and a JSON path before anything runs.
broken = json.loads(json.dumps(model))
broken["rules"][1]["rate"] = "max(0, K - N) / DD * sw"
try:
    rules.ensemble(broken, end=1.0, report_times=[0.0, 1.0])
except CrnError as error:
    print("rejected:", error)
```

```text
chacha8-u53-ziggurat/rules-mtdna-direct/1
cells past 60% within three years: 4.33%
chacha8-u53/rules-static-direct/1
['mtdna_direct ran Static Direct: stop is neither absent nor the guarded mutant fraction']
rejected: CRN-R012: unknown name `DD`; expected a parameter or observable in `max(0, K - N) / DD * sw` at /rules/1/rate
```

Validation is strict:

- unknown fields, duplicate keys and unknown names are errors;
- every diagnostic has a stable code (`CRN-R0xx`) and a JSON path;
- `rules.load(path)` checks a model file before any run;
- `crn-craft rules-check MODEL` does the same from the command line.

## 14. Budgets, wall time and cancellation

Every run has explicit limits:

| Limit | Python argument | Default | Termination when reached |
| --- | --- | --- | --- |
| events | `max_events` | 1,000,000 | `event_limit` (`CRN-R040`) |
| random words | `max_sampling_steps` | 8,000,000 | `sampling_limit` (`CRN-R041`) |
| live instances | `max_instances` | 100,000 | `instance_limit` (`CRN-R042`) |
| wall time | `wall_time_s` | none | `wall_time_limit` (`CRN-R043`) |
| result size | `max_output_bytes` | 64 MiB | `output_limit` (`CRN-R045`) |

Defaults are per cell, except wall time, which covers the whole run.
`max_events` must exceed the events of one cell: about 150,000 for three
years at K = 1000, and 2–3 times more with a doubled `K` or horizon. An
ensemble that has a cell ending on a budget withholds statistics
(`CRN-R050`), so increase the limit rather than average over truncated
cells.

```python
"""Budgets, wall time and cancellation."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))


def summary(result):
    terms = sorted(set(result.records["termination"]))
    codes = [d["code"] for d in result.metadata["diagnostics"]]
    return result.metadata["termination"], terms, codes


# Too few events for three years: every cell ends with event_limit and the
# ensemble withholds its statistics.
small = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=4,
    seed=1,
    max_events=10_000,
    method="mtdna_direct",
)
result = small.simulate()
print(summary(result), "statistics:", result.statistics)

# A wall-time budget stops the ensemble cooperatively.
slow = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=20_000,
    seed=1,
    threads=2,
    wall_time_s=1.0,
    method="mtdna_direct",
)
print(summary(slow.simulate()))
```

```text
('event_limit', ['event_limit'], ['CRN-R050']) statistics: None
('wall_time_limit', ['completed', 'wall_time_limit'], ['CRN-R050'])
```

Ctrl-C cancels a Python run cooperatively. `simulate(cancel_file=path)`
cancels when the file appears, which suits batch systems. The engine polls
every 1024 events.

## 15. The command line and JSON requests

The CLI runs the same requests as Python. A request binds the model file by
path and SHA-256, so a result always names the exact model it came from.

```python
"""Write a model and a run request for the command line, and run them.

Usage: python examples/mtdna/11_cli_request.py [DIRECTORY [CRN_CRAFT_BINARY]]

DIRECTORY defaults to `takeover` in the temporary directory, /tmp/takeover
on Linux. With the native CLI binary, the script also runs the request.
"""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

if len(sys.argv) > 1:
    folder = Path(sys.argv[1])
else:
    folder = Path(tempfile.gettempdir()) / "takeover"
folder.mkdir(parents=True, exist_ok=True)
model = json.loads(Path("models/rules/mtdna-three-class.rules.json").read_text())
model["stop"] = "N > 0 && mutants >= 0.9 * N"
data = (json.dumps(model, indent=2) + "\n").encode()
(folder / "takeover.rules.json").write_bytes(data)
request = {
    "schema": "crn-craft.rules-run/0.3",
    # The request binds the model by relative path and exact SHA-256.
    "model": {
        "path": "takeover.rules.json",
        "sha256": hashlib.sha256(data).hexdigest(),
    },
    "time": {"end": 25920.0, "report_times": [0.0, 8640.0, 17280.0, 25920.0]},
    "ensemble": {"replicates": 1000, "root_seed": 11},
    "outputs": {
        "statistics": ["mean", "variance"],
        "keep_paths": 0,
        "keep_instances": 0,
    },
    "budget": {
        "threads": 2,
        "max_events_per_replicate": 10_000_000,
        "max_sampling_steps_per_replicate": 100_000_000,
        "max_instances": 10,
        "max_output_bytes": 64 << 20,
    },
    "method": {"kind": "mtdna_direct"},
    "trackers": [{"observable": "advantaged_mutants"}],
}
(folder / "takeover.run.json").write_text(json.dumps(request, indent=2) + "\n")
print("model sha256:", request["model"]["sha256"])
if len(sys.argv) > 2:
    done = subprocess.run(
        [sys.argv[2], "rules", "takeover.rules.json", "takeover.run.json"],
        cwd=folder,
        capture_output=True,
        check=True,
    )
    result = json.loads(done.stdout)
    print(result["termination"], result["sampler"])
    print(
        "cells:",
        result["completed_replicates"],
        "stopped:",
        result["stopped_replicates"],
    )
    print("events:", result["events"])
    print("result fields:", sorted(result))
```

These commands need the native CLI (`cargo install crn-craft-cli --locked`),
not the command of the same name from the Python package; see section 1.

```sh
python examples/mtdna/11_cli_request.py /tmp/takeover
crn-craft rules-check /tmp/takeover/takeover.rules.json
crn-craft rules /tmp/takeover/takeover.rules.json /tmp/takeover/takeover.run.json > result.json
crn-craft describe rules-run      # the JSON Schema of the request
crn-craft describe rules-result   # and of the result
```

```text
model sha256: f9eaaa59e4e449a55bd991350c748fdd8118b7c79dd81d6bcb0fc57bf025c0c2
completed chacha8-u53-ziggurat/rules-mtdna-direct/1
cells: 1000 stopped: 107
events: 142356193
result fields: ['completed_replicates', 'diagnostics', 'events', 'final_instances', 'model_sha256', 'observables', 'observations', 'paths', 'replicates', 'report_times', 'requested_replicates', 'root_seed', 'rules', 'run_sha256', 'sampler', 'schema', 'statistics', 'stopped_replicates', 'termination', 'types']
```

The request versions decide which fields are allowed:

| Version | Adds |
| --- | --- |
| `crn-craft.rules-run/0.1` | the tree method only |
| `crn-craft.rules-run/0.2` | `method`, `trackers`, `event_log` |
| `crn-craft.rules-run/0.3` | `mtdna_direct` |
| `crn-craft.rules-run/0.4` | Static Direct `waiting_time` and `sparse_trackers` |
| `crn-craft.rules-run/0.5` | `hybrid` with its `rtol` and `atol`; the result is rules-result/0.3, which adds `method` and `continuous` |
| `crn-craft.rules-run/0.6` | no field; Static Direct defaults to Ziggurat waits and sparse trackers; the result is rules-result/0.3 |

`rules.ensemble` picks the lowest version that fits its arguments, and
rules-run/0.6 for `static_direct`.
`rules.prepare(model, request)` takes a request you wrote yourself. Stored
requests of every earlier version keep their meaning and replay
byte-identically.

Model files carry their own version: rules-model/0.1 is the base language,
0.2 adds the legacy `atp` block, and 0.3 adds whole-rule `propensity` and
`continuous` variables. `rules.Model` writes 0.3.

## 16. Rust

The crate `crn-craft-rules` is the engine behind Python and the CLI.
`prepare` validates exact model and request bytes, and `execute` runs them:

```rust
use crn_craft_rules::{Control, execute, prepare};

let prepared = prepare(&model_bytes, &request_bytes)?;
let result = execute(&prepared, &Control::default())?;
println!("{} cells stopped", result.stopped_replicates);
```

- `execute_with_observer(&prepared, &control, &observer)` calls an
  `EventObserver` for every committed event on the native workers. It
  returns `ObserverAction::Stop` to end that cell. An observer that always
  continues compiles away and leaves results byte-identical.
- *Opt-in:* `generate_rules_loop(&prepared)` emits Rust source
  for the whole Static Direct event of one prepared model. It requires a
  `static_direct` request, constant or vacancy rates and sum observables,
  and at most 1 MiB of source; it refuses anything else with `CRN-R020`.
  Parameters are compiled in, so another parameter set needs another loop.
  Compile it into your program and run it with `execute_with_rules_loop`,
  which runs only a loop generated for the same model bytes. Results are
  byte-identical to `execute`, and on the three-class model it is about
  2.4× faster than interpreted Static Direct. On networks with many
  constant-rate channels it can be slower, so measure first
  ([details](specialized-kernels.md#generated-rules-event-loops)). From
  Python, `rules.compile_loop` builds such a loop into a separate program
  and runs it (section 8).

## 17. Checking results

Exactness is checked against independent references, and the checks can be
rerun:

- **Master equation.** `validation/check_mtdna.py` compares the
  `mtdna-three-class-cme` model (K = 20) under every exact method, the fast
  path `mtdna_direct` included, with an exact solution of the chemical
  master equation, using DKW bounds at family α = 0.001. Approximate
  methods are compared with the exact ones.
- **Lumping.** The clone model is checked against the same master-equation
  law for its class totals.
- **Observation.** Trackers equal their recomputation from complete event
  logs, and results with and without observers are byte-identical.
- **Workers.** Results are byte-identical for one and many threads.
- **Java reference.** On 10,000 cells, `mtdna_direct` agreed with a Java
  implementation of the fixed-population model. The takeover-time
  distributions were within the declared bound, and the mutation counts
  within 1% (`validation/mtdna/java-gap.json`).

```sh
python validation/check_mtdna.py --binary target/release/crn-craft --output mtdna-gates.json
python validation/check_rules_python.py --output mtdna-python.json
```

A numerical check confirms that the simulator solves the model; it does not
validate the model against biology.

## 18. Reference

`rules.ensemble` arguments:

| Argument | Meaning |
| --- | --- |
| `model` | model bytes from `rules.load`, a dictionary, or a `rules.Model` (section 12) |
| `end`, `report_times` | horizon in hours; increasing report times in [0, end] |
| `replicates`, `seed`, `threads` | number of cells, root seed, 1–64 native threads |
| `keep_paths` | how many cells keep their full observable series (`result.paths`); `"all"` or `True` for every cell |
| `keep_instances` | how many cells keep their final instances (`metadata["final_instances"]`) |
| `method` | a method name or dictionary (section 8) |
| `compiled` | a loop from `rules.compile_loop` for this model, run by `simulate()` (section 8) |
| `trackers` | list of `{"observable", "first_passage"?, "keep_episodes"?}` |
| `event_log` | `{"max_entries": n, "rules": [...]}` |
| budgets | `max_events`, `max_sampling_steps`, `max_instances`, `max_output_bytes`, `wall_time_s` |

The result, `run.simulate()`:

| Attribute | Content |
| --- | --- |
| `records` | DataFrame, one row per cell, index `replicate`: `started` (bool), `termination` (categorical), `stop_time` (float), `events`, `rng_words`, `peak_instances`, `created_instances`, `saved_points` (unsigned integers) |
| `rule_events` | DataFrame, one row per cell, one column per rule name: events of each rule |
| `statistics` | `None` if withheld; else `samples` and the DataFrames `mean`, `sample_variance`, `mean_standard_error`, index `time`, one column per observable |
| `paths` | DataFrame of the kept series, index (`replicate`, `time`), one column per observable; `to_xarray()` gives a Dataset with dimensions `replicate` and `time` |
| `tracker_columns(i)` | DataFrame of tracker `i`, one row per observed cell; NaN where a time or peak is absent |
| `event_columns()`, `event_columns(id)` | DataFrame of every logged event, index (`replicate`, `entry`), or of one cell's log, index `entry` |
| `rule_names` | rule names in declaration order, also on the prepared run |
| `observations`, `record_diagnostics`, `metadata` | trackers and logs as documents, per-cell errors by id, the full result document |

Every DataFrame and column has `.to_numpy()`. Code written for the NumPy
results of release 0.2.0 changes as follows:

| 0.2.0 | Now |
| --- | --- |
| `records["termination"] == rules.TERMINATIONS.index("stopped")` | `records["termination"] == "stopped"`; the old form is silently all `False` |
| `records["id"]` | `records.index` |
| `records["rule_events"]` | `result.rule_events` |
| `records[name]` (an array) | `records[name].to_numpy()`; integer termination codes are `records["termination"].cat.codes` |
| `statistics["mean"][:, names.index("N")]` | `statistics["mean"]["N"]` |
| `statistics["mean"]` (an array) | `statistics["mean"].to_numpy()` |
| `metadata["paths"]` rows | `result.paths`, or `result.paths.to_numpy()` |
| `tracker_columns(i)["id"]` and fields | `tracker_columns(i).index` and its columns |
| `event_columns(id)["observables"]` | `event_columns(id)[list(metadata["observables"])].to_numpy()` |
| `keep_paths=True` (one cell) | `keep_paths=1`; `True` now keeps every cell |

Result metadata:

| Key | Content |
| --- | --- |
| `sampler`, `termination`, `diagnostics` | what ran, how the ensemble ended, warnings and errors |
| `events`, `completed_replicates`, `stopped_replicates` | totals |
| `rules`, `observables`, `types`, `report_times` | names and grid, in result order; rule names also as `run.rule_names` and `result.rule_names` |
| `statistics` | `mean`, `sample_variance`, `mean_standard_error`, `samples` (also `result.statistics`) |
| `paths`, `final_instances`, `observations` | kept series, kept instances, and (with trackers or a log) trackers and logs |
| `model_sha256`, `run_sha256`, `root_seed` | provenance |
| `method`, `continuous` | rules-result/0.3 (`hybrid`, `static_direct`): the request's method, with the tolerances of `hybrid`, and the names of the continuous columns |

Terminations, the categories of `records["termination"]` in the order of
`rules.TERMINATIONS`: `completed`, `stopped`, `event_limit`,
`sampling_limit`, `instance_limit`, `wall_time_limit`, `cancelled`,
`failed`, `output_limit`.

Frequent diagnostics:

| Code | Meaning |
| --- | --- |
| `CRN-R012` | expression error or unknown name |
| `CRN-R020` | invalid request: hash, budget, method, tracker, log or version |
| `CRN-R030` | invalid rate, overflow or undefined stop in a cell |
| `CRN-R031` | a positive channel below 53-bit resolution of the total; use `tree` or rescale |
| `CRN-R040`–`R045` | a budget was reached |
| `CRN-R050` | incomplete ensemble, statistics withheld |
| `CRN-R051` | a requested specialization did not apply; the message says why |

## 19. Troubleshooting

- **`event_limit` in every cell.** Raise `max_events` above the events of
  one cell (section 14).
- **`mtdna_direct` is not faster.** Check `metadata["sampler"]` and the
  `CRN-R051` warning, which names the failed check of section 8. A stop
  other than the guarded fraction, unequal degradation rates or a rate
  outside `max(0, K - N) / D` times constants make it fall back to Static
  Direct.
- **`observation reservation exceeds budget`.** Event logs reserve space for
  every cell before the run. The message itemizes the reservation (cells ×
  bytes per cell, log entries × bytes per entry), and its hint
  (`error.diagnostics[0]["hint"]`) names the smallest `max_output_bytes`
  that admits the request. Raise `max_output_bytes` to that value, or use
  fewer cells or entries. For the three-class model an entry takes 721
  bytes, so 50 entries for 1,000 cells fit the 64 MiB default.
- **`CRN-R020` at `/method` for a model with continuous variables.** Only
  `method="hybrid"` runs them (section 10); stored models with an `atp`
  block need `thinning`.
- **A Python callback is slow.** It holds the interpreter lock. Use trackers
  or event logs for whole ensembles.
- **`rules.compile_loop` fails.** `CompiledLoopError.reason` says why:
  `"cargo_missing"` (install Rust with rustup), `"build_failed"` (read
  `error.log`; for a version that is not on crates.io pass `source=` a
  checkout), or `"mismatch"`. A `CrnError` with `CRN-R020` at `/method`
  means the generator does not support the model (section 8).
- **Different results with another method.** Each method uses the random
  stream differently. Compare distributions, not paths.
- **`CRN-R031`.** A rate is below 2⁻⁵³ of the total, for example a tiny
  mutation probability in a huge population. Static methods select in one
  level; the tree method resolves rules and instances separately.
- **No cell stopped, though some should have.** Compare `termination` with
  a name, `records["termination"] == "stopped"`. A comparison with the
  integer `rules.TERMINATIONS.index("stopped")`, as in release 0.2.0, is
  `False` for every cell.
- **Column positions.** Select DataFrame columns by observable name. In the
  raw matrices of `result.metadata`, observables are sorted by name; use
  `metadata["observables"].index(name)`.
