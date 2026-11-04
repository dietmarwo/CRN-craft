# Simulating mtDNA mutant accumulation: user manual

This manual is for researchers who simulate how mitochondrial DNA mutants
arise and spread in single cells. It shows how to do it with the rules
engine of CRN-craft, from Python, from the command line and from Rust. Every
scenario comes with a complete example in `examples/mtdna/`. The examples
were run as printed. Each one's output is shown below it, timings from a
two-CPU Linux host.

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
10. [Variable population: the ATP model](#10-variable-population-the-atp-model)
11. [Parameter sweeps](#11-parameter-sweeps)
12. [Model variants](#12-model-variants)
13. [Budgets, wall time and cancellation](#13-budgets-wall-time-and-cancellation)
14. [The command line and JSON requests](#14-the-command-line-and-json-requests)
15. [Rust](#15-rust)
16. [Checking results](#16-checking-results)
17. [Reference](#17-reference)
18. [Troubleshooting](#18-troubleshooting)

## 1. Installation

```sh
pip install crn-craft                    # Python package (abi3 wheels for CPython 3.10+)
cargo install crn-craft-cli --locked     # optional: the native command-line tool
```

The Python package contains the native engine, so the CLI is optional. Every
Python feature in this manual is in release 0.2.0. Features that are newer
than the release are Rust only and are marked where they appear. To use the
current source instead of the release, build the wheel from a clone of the
repository:

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
Replication can produce a mutant. The repository has four model files:

| File | What it models | Methods |
| --- | --- | --- |
| `mtdna-three-class.rules.json` | fixed population target; wild type, neutral mutants and advantaged mutants as three classes | every static method, `mtdna_direct`, `tree` |
| `mtdna-three-class-cme.rules.json` | the same with K = 20 and fast rates, for exact checks against the master equation | as above |
| `mtdna-three-class-atp.rules.json` | variable population: replication limited by a cellular ATP balance | `thinning` |
| `mtdna-fixed.rules.json` | every new mutant is its own clone with an identity | `tree` |

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

The order of observables in results is sorted by name, not the order in the
file. Look up a column with `result.metadata["observables"].index("N")`
instead of assuming a position.

## 3. A first ensemble

An *ensemble* is a set of independent replicates of the model. Here each
replicate is one cell. `rules.ensemble` builds and validates a run request
and prepares it; `simulate()` runs the cells on native threads and returns
NumPy columns.

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
records = result.records
stopped = records["termination"] == rules.TERMINATIONS.index("stopped")
print(meta["sampler"], meta["termination"])
print(f"{meta['events']:,} events in {seconds:.2f} s")
print(f"taken over: {stopped.sum()} of {len(stopped)} cells")
print("takeover days (first five):", np.sort(records["stop_time"][stopped])[:5] / 24)
names = meta["observables"]
mean = result.statistics["mean"]
fraction = mean[:, names.index("mutants")] / mean[:, names.index("N")]
for day, value in zip(days, fraction, strict=True):
    print(f"day {day:5d}: mean mutant fraction {value:.4f}")
```

```text
chacha8-u53-ziggurat/rules-mtdna-direct/1 completed
28,536,814 events in 0.53 s
taken over: 23 of 200 cells
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

- `result.metadata`: the full result document (sampler, termination,
  totals, statistics, kept paths, observations, diagnostics).
- `result.records`: one entry per cell. The columns are `id`, `started`,
  `termination`, `stop_time`, `events`, `rng_words`, `peak_instances`,
  `created_instances`, `saved_points` and `rule_events` (cells × rules).
- `result.statistics`: `mean`, `sample_variance` and `mean_standard_error`
  (report times × observables), and `samples`. They are withheld if any
  cell ended early for a budget reason.

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

import numpy as np
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
        keep_paths=64,
    )
    return ensemble.simulate()


one, two, other = run(7, 1), run(7, 2), run(8, 2)
same = all(np.array_equal(one.records[k], two.records[k]) for k in one.records)
print("1 vs 2 threads, same seed, identical records:", same)
print(
    "same statistics:", np.array_equal(one.statistics["mean"], two.statistics["mean"])
)
print(
    "seed 7 vs 8, identical stop times:",
    np.array_equal(one.records["stop_time"], other.records["stop_time"]),
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
    np.array_equal(small.records["events"], one.records["events"][:8]),
)
```

```text
1 vs 2 threads, same seed, identical records: True
same statistics: True
seed 7 vs 8, identical stop times: False
first 8 of 64 equal an 8-cell run: True
```

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
stopped = result.records["termination"] == rules.TERMINATIONS.index("stopped")
ages = 24.0 * np.array([360, 720, 1080])
taken = [(stopped & (result.records["stop_time"] <= a)).mean() for a in ages]
```

Any expression over parameters and observables can be a stop, for example
`mutants >= 0.6 * N` for a 60% COX-deficiency threshold, or
`advantaged_mutants > 0` to stop at the first advantaged mutant.
`mtdna_direct` runs at full speed with no stop or with the guarded fraction
form `N > 0 && mutants >= f * N`. Any other stop makes it fall back to Static
Direct with a warning (section 12). To record a threshold without ending
the cell, use a first-passage tracker (next section).

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

import numpy as np
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
neutral = result.tracker_columns(0)
advantaged = result.tracker_columns(1)
passage = result.tracker_columns(2)
stop = result.records["stop_time"]
print("neutral lineages arising per cell:", neutral["appearances"].mean())
print("neutral lineages lost per cell:   ", neutral["extinctions"].mean())
print(
    "advantaged lineages arising/lost: ",
    advantaged["appearances"].mean(),
    advantaged["extinctions"].mean(),
)
print("largest advantaged class seen:    ", advantaged["maximum"].max())
# The integral runs until the replicate's stop or the horizon.
load = advantaged["integral"] / stop
print("time-averaged advantaged count, mean over cells:", load.mean())
hit = ~np.isnan(passage["first_passage"])
print(
    f"cells reaching 60% mutants: {hit.sum()}; median day "
    f"{np.median(passage['first_passage'][hit]) / 24:.0f}"
)
# Episodes: each excursion of a class above zero, with start, end, peak.
first = next(o for o in result.observations if o["trackers"][1]["episodes"])
for episode in first["trackers"][1]["episodes"]:
    print("cell", first["id"], "advantaged episode:", episode)
```

```text
sampler: chacha8-u53-ziggurat/rules-mtdna-direct/1
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

For plain count trackers, Static Direct offers sparse integration
(`{"kind": "static_direct", "sparse_trackers": True}`), and `mtdna_direct`
does it automatically: an integral is advanced only when its class changes.
The last bits of an integral can then differ from per-event integration.

## 7. Exact event records: logs and callbacks

An *event log* keeps, for every cell, up to `max_entries` events of the
listed rules (all rules if the list is empty). Each entry holds the time,
the rule, the participating instances and the observables after the event.
A *live callback* is a Python function called for every committed event. It
may stop a cell by returning `False`.

```python
"""Exact event records: the event log and live callbacks."""

import json
from collections import Counter

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
names = [rule["name"] for rule in model["rules"]]

# 1. A bounded log of selected rules for every cell.
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
total = Counter()
for observation in result.observations:
    for entry in observation["log"]:
        total[names[entry["rule"]]] += 1
print("logged mutation events:", dict(total))
print("omitted (log full):", sum(o["omitted_log_entries"] for o in result.observations))
columns = result.event_columns(result.observations[0]["id"])
print(
    "cell",
    result.observations[0]["id"],
    "mutation days:",
    np.round(columns["time"] / 24, 1),
)
observables = result.metadata["observables"]
print("N after each logged event:", columns["observables"][:, observables.index("N")])

# 2. A live callback sees every committed event in order; returning False
#    stops that replicate. Callbacks hold the GIL, so use them on few cells.
seen = Counter()


def watch(replicate, event):
    seen[event["rule"]] += 1
    # Stop a cell as soon as an advantaged mutant arises.
    return event["rule"] != names.index("mutation_mtAdv")


run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=20,
    seed=4,
    threads=1,
    method="static_direct",
)
result = run.simulate(callback=watch)
stopped = result.records["termination"] == rules.TERMINATIONS.index("stopped")
print("events seen by the callback:", sum(seen.values()))
print("cells stopped at their first advantaged mutation:", stopped.sum())
print("their stop days:", np.round(result.records["stop_time"][stopped] / 24, 1))
```

```text
logged mutation events: {'mutation_mt': 157, 'mutation_mtAdv': 40}
omitted (log full): 0
cell 0 mutation days: [323.4 744.1]
N after each logged event: [997. 999.]
events seen by the callback: 2512850
cells stopped at their first advantaged mutation: 7
their stop days: [ 744.1  364.2  914.2  232.4   47.8 1061.3  747.2]
```

Both runs use Static Direct with seed 4. Cell 0's mutation days in the log
(323.4 and 744.1) reappear as the stop day of a cell stopped at its first
advantaged mutation. An observer reads the same events; it does not change
them.

Rules for observation:

- Observers read state and may stop a replicate. They never change the
  state or draw random numbers, so attaching one leaves results unchanged.
- A callback receives each event as a parsed JSON document and holds
  Python's interpreter lock. It cost about 14 µs per event on the test host,
  a hundred times the simulation itself, and it serializes the threads. Use
  it on a few cells, and use trackers or logs for whole ensembles.
- Log retention counts against the admission budget. Large logs need small
  cohorts and a larger `max_output_bytes`.
- `batch_callback=fn` replays the retained log after the run in chunks, and
  `result.event_columns(id)` returns time, rule and observables as arrays.

## 8. Choosing a method

| `method` | Kind | Use it for |
| --- | --- | --- |
| `"mtdna_direct"` | exact | the three-class model: certified once at preparation, about 3× faster than Static Direct |
| `"static_direct"` | exact | any closed model (fixed set of attribute vectors); the default for most work |
| `{"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": True}` | exact | Static Direct with faster waiting times and sparse count trackers; its own replay identity |
| `"next_reaction"` | exact | Gibson–Bruck; pays off in large sparse networks, not here, where every event changes every rate |
| `None` (tree) | exact | models whose instances change identity, such as the clone model |
| `{"kind": "tau", "epsilon": 0.03}` | approximate | adaptive Poisson leaps; at K = 1000 the vacancy limits leaps, so it runs exact events |
| `{"kind": "binomial", "step": 1.0}` | approximate | fixed hourly splitting, similar to a binomial hourly scheme in Java |
| `"thinning"` | exact | the ATP model (section 10) |

The same 64 cells under each method:

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
    "static_direct, ziggurat": {"kind": "static_direct", "waiting_time": "ziggurat"},
    "mtdna_direct": "mtdna_direct",
    "next_reaction": "next_reaction",
    "tau (epsilon 0.03)": {"kind": "tau", "epsilon": 0.03},
    "binomial (1 h steps)": {"kind": "binomial", "step": 1.0},
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
    stopped = result.records["termination"] == rules.TERMINATIONS.index("stopped")
    mean = result.statistics["mean"]
    names = meta["observables"]
    fraction = mean[-1, names.index("mutants")] / mean[-1, names.index("N")]
    print(
        f"{label:24s} {seconds:6.2f} s  {meta['events']:>11,} events  "
        f"taken over {stopped.sum():2d}  final mutant fraction {fraction:.3f}  "
        f"{meta['sampler']}"
    )
```

```text
tree                       2.58 s    9,167,872 events  taken over  6  final mutant fraction 0.097  chacha8-u53/rules-direct/2
static_direct              0.72 s    9,180,176 events  taken over  8  final mutant fraction 0.144  chacha8-u53/rules-static-direct/1
static_direct, ziggurat    0.56 s    8,861,769 events  taken over 10  final mutant fraction 0.164  chacha8-u53-ziggurat/rules-static-direct/1
mtdna_direct               0.18 s    8,969,870 events  taken over  8  final mutant fraction 0.122  chacha8-u53-ziggurat/rules-mtdna-direct/1
next_reaction              1.03 s    9,226,649 events  taken over  7  final mutant fraction 0.105  chacha8-u53/rules-next-reaction/1
tau (epsilon 0.03)         4.08 s    9,225,966 events  taken over  5  final mutant fraction 0.108  chacha8-u53/rules-approximate-cgp-poisson/2
binomial (1 h steps)       0.68 s    9,158,734 events  taken over  6  final mutant fraction 0.088  chacha8-u53/rules-approximate-binomial-split/2
```

The takeover counts differ because each method uses the random stream
differently. With 64 cells, 5 to 10 takeovers are within sampling error of
each other. The approximate methods are statistically close here but are
labelled as approximations in their sampler identity. The
[performance analysis](mtdna-performance-analysis.md) compares the methods
at full scale.

`mtdna_direct` checks at preparation that the model has the three-class
structure, that every positive hazard stays resolvable in every reachable
state, and that the stop is absent or the guarded fraction. If any check
fails, the run uses Static Direct, `metadata["sampler"]` says so, and a
`CRN-R051` warning in `metadata["diagnostics"]` names the reason.

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
names = result.metadata["observables"]
print("mean number of distinct mutant clones:", mean[:, names.index("mutant_clones")])
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

## 10. Variable population: the ATP model

`mtdna-three-class-atp.rules.json` (rules-model/0.2) lets the population
float. Replication is limited by an ATP pool instead of the vacancy
`K − N`. ATP follows dA/dt = f · wt − v1 · A − v2 · N between events, solved
exactly, and the replication rules are multiplied by c/(c + A). The
`thinning` method samples events exactly against a bound on each interval.
ATP appears as an extra continuous observable, `ATP`. A path also stops when
ATP is exhausted.

```python
"""Variable population with the ATP balance (exact thinning)."""

import json

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class-atp.rules.json"))
print("ATP block:", model["atp"])
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=24.0 * np.arange(0, 1081, 180),
    replicates=64,
    seed=7,
    threads=2,
    method="thinning",
)
result = run.simulate()
meta = result.metadata
names = meta["observables"]
print("observables:", names)
print("sampler:", meta["sampler"], "stopped:", meta["stopped_replicates"])
mean = result.statistics["mean"]
for day, row in zip(np.arange(0, 1081, 180), mean, strict=True):
    print(
        f"day {day:4d}: N {row[names.index('N')]:7.1f}  "
        f"mutants {row[names.index('mutants')]:6.1f}  ATP {row[names.index('ATP')]:7.2f}"
    )
```

```text
ATP block: {'wild_type': 'wild_type', 'total': 'N', 'initial': 450.0, 'production': 0.1, 'decay': 0.2, 'maintenance': 0.01, 'scale': 1.3034153747670567, 'replication_rules': ['replication_wt', 'mutation_mt', 'mutation_mtAdv', 'replication_mt', 'replication_mtAdv']}
observables: ['N', 'advantaged_mutants', 'mutants', 'neutral_mutants', 'wild_type', 'ATP']
sampler: chacha8-u53/rules-atp-thinning/1 stopped: 8
day    0: N  1000.0  mutants    0.0  ATP  450.00
day  180: N   998.9  mutants    0.2  ATP  449.16
day  360: N  1576.6  mutants  540.5  ATP  441.89
day  540: N  3599.7  mutants 2418.4  ATP  422.84
day  720: N  4765.5  mutants 3510.3  ATP  407.80
day  900: N  5352.7  mutants 4055.7  ATP  401.87
day 1080: N  5919.9  mutants 4587.1  ATP  394.85
```

Event trackers can follow the discrete observables of the ATP model, but
not `ATP` itself.

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
        stopped = result.records["termination"] == rules.TERMINATIONS.index("stopped")
        row.append(stopped.mean())
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

## 12. Model variants

The model file is plain JSON. Load it with `json.loads(rules.load(path))`,
change it, and pass the dictionary to `rules.ensemble`.

- **Copy number and half-life:** set `K`, `w0` and `delta = ln 2 / half-life`.
- **Initial heteroplasmy:** add an `initial` entry for the mutant class.
- **Threshold:** change the fraction in `stop`.
- **Recognition:** see what happens when `mtdna_direct` does not recognize
  the model.
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
stopped = result.records["termination"] == rules.TERMINATIONS.index("stopped")
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

## 13. Budgets, wall time and cancellation

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
    terms = [rules.TERMINATIONS[t] for t in result.records["termination"]]
    codes = [d["code"] for d in result.metadata["diagnostics"]]
    return result.metadata["termination"], sorted(set(terms)), codes


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

## 14. The command line and JSON requests

The CLI runs the same requests as Python. A request binds the model file by
path and SHA-256, so a result always names the exact model it came from.

```python
"""Write a model and a run request for the command line, and run them.

Usage: python examples/mtdna/11_cli_request.py DIRECTORY [CRN_CRAFT_BINARY]
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

folder = Path(sys.argv[1])
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

`rules.ensemble` picks the lowest version that fits its arguments.
`rules.prepare(model, request)` takes a request you wrote yourself.

## 15. Rust

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
- *Newer than 0.2.0:* `generate_rules_loop(&prepared)` emits Rust source
  for the whole Static Direct event of a closed model with constant or
  vacancy rates, including its stop. Compile it into your program and run
  it with `execute_with_rules_loop`. Results are byte-identical to
  `execute`, and on the three-class model it is about 2.4× faster than
  interpreted Static Direct
  ([details](specialized-kernels.md#generated-rules-event-loops)). From
  Python, use `mtdna_direct` for the same speed.

## 16. Checking results

Exactness is checked against independent references, and the checks can be
rerun:

- **Master equation.** `validation/check_mtdna.py` compares the
  `mtdna-three-class-cme` model (K = 20) under every exact method with an
  exact solution of the chemical master equation, using DKW bounds at family
  α = 0.001. Approximate methods are compared with the exact ones.
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

## 17. Reference

`rules.ensemble` arguments:

| Argument | Meaning |
| --- | --- |
| `model` | model bytes from `rules.load`, or a dictionary |
| `end`, `report_times` | horizon in hours; increasing report times in [0, end] |
| `replicates`, `seed`, `threads` | number of cells, root seed, 1–64 native threads |
| `keep_paths` | how many cells keep their full observable series (`metadata["paths"]`) |
| `keep_instances` | how many cells keep their final instances (`metadata["final_instances"]`) |
| `method` | a method name or dictionary (section 8) |
| `trackers` | list of `{"observable", "first_passage"?, "keep_episodes"?}` |
| `event_log` | `{"max_entries": n, "rules": [...]}` |
| budgets | `max_events`, `max_sampling_steps`, `max_instances`, `max_output_bytes`, `wall_time_s` |

Result metadata:

| Key | Content |
| --- | --- |
| `sampler`, `termination`, `diagnostics` | what ran, how the ensemble ended, warnings and errors |
| `events`, `completed_replicates`, `stopped_replicates` | totals |
| `rules`, `observables`, `types`, `report_times` | names and grid, in result order |
| `statistics` | `mean`, `sample_variance`, `mean_standard_error`, `samples` (also `result.statistics`) |
| `paths`, `final_instances`, `observations` | kept series, kept instances, and (with trackers or a log) trackers and logs |
| `model_sha256`, `run_sha256`, `root_seed` | provenance |

Terminations, indexed by `rules.TERMINATIONS`: `completed`, `stopped`,
`event_limit`, `sampling_limit`, `instance_limit`, `wall_time_limit`,
`cancelled`, `failed`, `output_limit`.

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

## 18. Troubleshooting

- **`event_limit` in every cell.** Raise `max_events` above the events of
  one cell (section 13).
- **`mtdna_direct` is not faster.** Check `metadata["sampler"]` and the
  `CRN-R051` warning. A stop other than the guarded fraction, or another
  model change, makes it fall back.
- **`observation reservation exceeds budget`.** Event logs reserve space for
  every cell. Use fewer cells or entries, or raise `max_output_bytes`.
- **A Python callback is slow.** It holds the interpreter lock. Use trackers
  or event logs for whole ensembles.
- **Different results with another method.** Each method uses the random
  stream differently. Compare distributions, not paths.
- **`CRN-R031`.** A rate is below 2⁻⁵³ of the total, for example a tiny
  mutation probability in a huge population. Static methods select in one
  level; the tree method resolves rules and instances separately.
- **Column positions.** Observables are sorted by name. Use
  `metadata["observables"].index(name)`.
