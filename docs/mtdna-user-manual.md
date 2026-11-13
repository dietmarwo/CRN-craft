# Simulating mtDNA mutant accumulation: user manual

This manual is for researchers who simulate how mitochondrial DNA mutants
arise and spread in single cells. It shows how to do it with CRN-craft, from
Python, from the command line and from Rust. Every Python block is a complete
program: copy it and run it from the repository root. Most have a counterpart in
`examples/mtdna/`, written with `crn_craft.rules`. The examples were run as printed. Each one's output is
shown below it, timings from a two-CPU Linux host.
`validation/check_manual.py` runs every block and compares its output with
the manual.

The blocks use `import crn_craft as crn`. The module `crn_craft.rules` keeps
the API of release 0.6.0 under its old name: code written for 0.6.0 with
`from crn_craft import rules` keeps running with the same documents and the
same results. The top level has the same names with the defaults of the 1.0
contracts; section 18 lists the differences.

Companion documents:

- [Optimization tutorial](mtdna-optimization-tutorial.md): fitting model
  parameters to data with fcmaes.
- [Performance analysis](evidence/mtdna-performance-analysis.md): what each method
  and option costs, and how to size a study.
- [Model reference](reference/model-language.md) and [mtDNA ensembles](reference/methods.md):
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
cargo install crn-craft-cli --locked     # optional: the native command-line tool
```

The Python package needs NumPy, pandas and xarray, which pip installs with
it. Results are pandas DataFrames (pandas is required since 0.3.0), and
kept paths and event logs are xarray DataArrays (xarray is required since
0.5.0). The extra `crn-craft[xarray]` of earlier releases still installs.

From R, the Python package works through reticulate; see
[Using CRN-craft from R](using-from-r.md).

The Python package contains the native engine, so the CLI is optional. Every
feature in this manual is in release 0.7.0; features that exist only in
Rust are marked where they appear.

Both packages install a command named `crn-craft`, with the same commands and
the same reports (`describe`, `check`, `run`, `verify`, `analyze`, `convert`,
`export`). If both are installed, the one that comes first on your `PATH`
wins. `crn-craft --version` tells them apart (the Python command prints
`"interface": "python"`); call the native binary by its full path, by default
`~/.cargo/bin/crn-craft`, if needed. Section 15 shows the commands.

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

The engine simulates models whose species types carry numeric
attributes, and whose rules act on every instance that matches a pattern. The
mtDNA models describe one cell. Its molecules degrade at a constant rate
and replicate into the free places below a copy-number target `K`.
Replication can produce a mutant. The repository has five model files:

| File | What it models | Methods |
| --- | --- | --- |
| `mtdna-three-class.rules.json` | fixed population target; wild type, neutral mutants and advantaged mutants as three classes | every generic method, and the opt-in fast path `mtdna_direct` |
| `mtdna-three-class-cme.rules.json` | the same with K = 20 and fast rates, for exact checks against the master equation | as above |
| `mtdna-three-class-atp-continuous.rules.json` | variable population: replication limited by a cellular ATP balance, written as a continuous variable (`rules-model/0.3`) | `hybrid` |
| `mtdna-three-class-atp.rules.json` | the same model in the stored rules-model/0.2 form, with a fixed `atp` block | `thinning` |
| `mtdna-fixed.rules.json` | every new mutant is its own clone with an identity | `tree` |

Section 12 shows how to write models like these from Python with
`crn.Model` instead of JSON.

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
select a column by name, as in `result.statistics["mean"]["N"]`. Kept paths
and event logs are xarray DataArrays whose `variable` coordinate carries the
names, as in `result.paths.sel(variable="N")`. Observables are sorted by
name, not in the order of the file, and continuous variables follow them
(section 10); the raw matrices in `result.metadata` use the same order.

## 3. A first ensemble

An *ensemble* is a set of independent replicates of the model. Here each
replicate is one cell. `crn.ensemble` builds and validates a run request
(`crn-craft.run/1.0`) and prepares it; `simulate()` runs the cells on native
threads and returns pandas DataFrames.

Most examples pass `method="mtdna_direct"`, an opt-in fast path for the
three-class structure of this model (section 8). `method="static_direct"`
runs the same model, and any other closed model, about 2.8× slower on this
example (measured on 12 November 2026); the numbers then differ cell by
cell but not in distribution.

```python
"""Three-class model: 200 cells for three years with the takeover stop."""

import json
import time

import numpy as np
import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
days = np.arange(0, 1081, 90)
run = crn.ensemble(
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
  observable, `median` and `quantiles` (next section), and `samples`, the
  number of cells. They are withheld (`None`) if any cell ended early for a
  budget reason.
- `result.paths`: the kept series of single cells, an xarray DataArray
  (next section).
- `result.metadata`: the full result document (sampler, termination,
  totals, statistics, kept paths, observations, diagnostics).

`.to_numpy()` turns any DataFrame, column or DataArray into a NumPy array.

Report times are right-continuous: an event exactly at a report time is
included. A cell that stops holds its state at every later report time.

### Quantiles and the median

Means hide the shape of the distribution: in the three-class model most cells
stay near zero mutants while a few are taken over, so the mean sits far from
the typical cell. `crn.ensemble` therefore also computes quantiles of the
cells' values at every report time, for every observable. By default it asks
for the probabilities 0.025, 0.5 and 0.975, so `result.statistics["median"]`
is the median and `result.statistics["quantiles"]` is a dictionary from each
probability to a DataFrame shaped like `mean`. The values are NumPy's default
quantiles (linear interpolation, type 7), equal bit for bit to
`numpy.quantile` over the cells. Pass `quantiles=[0.05, 0.25, 0.5, 0.75,
0.95]` for others, up to 16 strictly increasing probabilities in (0, 1), or
`quantiles=()` for none. Quantiles keep one value per cell, report time and
observable in memory while the run lasts.

```python
"""Quantiles and the median of the mutant count."""

import numpy as np
import crn_craft as crn

model = crn.load("models/rules/mtdna-three-class.rules.json")
days = np.arange(0, 1081, 360)
result = crn.ensemble(
    model,
    end=25_920.0,
    report_times=24.0 * days,
    replicates=200,
    seed=1,
    threads=2,
    method="mtdna_direct",
    quantiles=[0.025, 0.25, 0.5, 0.75, 0.975],
).simulate()
statistics = result.statistics
print("quantile levels:", list(statistics["quantiles"]))
print("mean mutants:  ", statistics["mean"]["mutants"].round(1).tolist())
print("median mutants:", statistics["median"]["mutants"].tolist())
print("97.5% quantile:", statistics["quantiles"][0.975]["mutants"].round(1).tolist())
paths = crn.ensemble(
    model,
    end=25_920.0,
    report_times=24.0 * days,
    replicates=200,
    seed=1,
    threads=2,
    method="mtdna_direct",
    keep_paths="all",
).simulate()
last = paths.paths.sel(time=25_920.0, variable="mutants").to_numpy()
print("equal to numpy:", np.quantile(last, 0.975) == statistics["quantiles"][0.975]["mutants"].iloc[-1])
```

```text
quantile levels: [0.025, 0.25, 0.5, 0.75, 0.975]
mean mutants:   [0.0, 25.2, 67.0, 120.7]
median mutants: [0.0, 0.0, 0.0, 0.0]
97.5% quantile: [0.0, 336.4, 995.0, 998.0]
equal to numpy: True
```

The median is zero at every report time, while the mean rises to 121: the
typical cell has no mutant, and the mean is carried by the cells that were
taken over. The 97.5% quantile shows how far a cell can go. Quantiles need all
cells: like the other statistics they are withheld when a cell ended on a
budget.

### Saving results

`result.save(file_name)` writes the model, the simulation settings and the
results to one file. The extension chooses the format: `.md` writes a markdown
report that a person can read, `.json` writes `crn-craft.report/1.0`, a JSON
document with the model, the request and the results, and `crn.load_report`
reads it back. `paths=False` leaves the kept trajectories out. The same result
always gives the same bytes: a report holds no timestamp, host name or timing.
It names the package version and both SHA-256 hashes, so a report identifies
the exact model and request. A report is a readable record of a run. For a
sealed record that `crn-craft verify` can check, use an archive (section 15).

```python
"""Save a result as a markdown report and as a JSON report."""

import tempfile
from pathlib import Path

import crn_craft as crn

model = crn.load("models/rules/mtdna-three-class.rules.json")
result = crn.ensemble(
    model,
    end=8_640.0,
    report_times=[0.0, 4_320.0, 8_640.0],
    replicates=50,
    seed=2,
    threads=2,
    method="mtdna_direct",
).simulate()
with tempfile.TemporaryDirectory() as folder:
    markdown = result.save(Path(folder) / "takeover.md")
    saved = result.save(Path(folder) / "takeover.json", paths=False)
    print(markdown.name, saved.name)
    report = crn.load_report(saved)
    print(report["schema"], sorted(report))
    print("same sampler:", report["result"]["sampler"] == result.metadata["sampler"])
    # The markdown report: one section per statistic.
    text = markdown.read_text()
    print(text[text.index("### median") : text.index("### 97.5% percentile")].rstrip())
```

````text
takeover.md takeover.json
crn-craft.report/1.0 ['crn_craft_version', 'model', 'request', 'result', 'schema']
same sampler: True
### median

```text
D: 1000, 997, 997.5
N: 1000, 997, 997
advantaged_mutants: 0, 0, 0
mutants: 0, 0, 0
neutral_mutants: 0, 0, 0
wild_type: 1000, 997, 997
```
````

`result.model` and `result.request` give the documents a result came from, as
new dictionaries; `result.save` raises `ValueError` for a result that does not
know them (one built by hand from a metadata document).

## 4. Seeds, threads and reproducibility

Every cell draws from its own ChaCha8 random stream, keyed by the root seed
and the cell's index. Results therefore depend on the seed and the method,
but not on the number of threads. A cell's trajectory is also the same in a
larger ensemble.

```python
"""Seeds, threads and reproducibility."""

import json

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"


def run(seed, threads):
    ensemble = crn.ensemble(
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
small = crn.ensemble(
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
paths = one.paths  # cells x report times x observables, an xarray DataArray
print(dict(paths.sizes))
print(paths.sel(replicate=0, variable=["N", "mutants"]).to_pandas())
```

```text
1 vs 2 threads, same seed, identical records: True
same statistics: True
same paths: True
seed 7 vs 8, identical stop times: False
first 8 of 64 equal an 8-cell run: True
{'replicate': 64, 'time': 4, 'variable': 6}
variable       N  mutants
time                     
0.0       1000.0      0.0
8640.0     998.0      1.0
17280.0    997.0      0.0
25920.0    999.0      0.0
```

`keep_paths` decides which cells keep their observable series: a count
keeps cells 0, 1, … up to it, and `"all"` (or `True`) keeps every cell.
`result.paths` holds them as one xarray DataArray with dimensions
(`replicate`, `time`, `variable`): cells × report times × observables, the
order of `event_columns()` (section 7). The coordinates are the cell
indices, the report times and the observable names (continuous variables
last); the values are floats. `paths.attrs`
names the `sampler`, the `method` (as JSON text, `"null"` for results
that report none; `json.loads` turns it back into a dictionary), the
`root_seed`, `model_sha256` and `run_sha256`, all of which netCDF can store. Some common
forms:

- `result.paths.sel(replicate=3)`: cell 3's series, report times ×
  observables;
- `result.paths.sel(time=8640.0)`: every kept cell at day 360;
- `result.paths.sel(variable="mutants")`: cells × report times;
- `result.paths.sel(variable="mutants").mean("replicate")`: the mean over
  the kept cells at each report time;
- `result.paths.transpose("time", "variable", "replicate")`: the
  dimension order of release 0.5.0;
- `result.paths_frame()`: the DataFrame of release 0.4.0, indexed by
  (`replicate`, `time`), one column per observable;
- `result.to_xarray()`: an xarray Dataset with dimensions `replicate` and
  `time`, one variable per observable: the same values as
  `paths.to_dataset("variable")`.

A cell that ended on a budget has values only up to its last saved report
time and NaN after it; `paths_frame()` has rows only up to that time. Kept
paths count against `max_output_bytes`; the
[performance analysis](evidence/mtdna-performance-analysis.md#8-sizing-a-study)
gives their size.

Each method has its own sampler identity in `metadata["sampler"]`. The same
seed gives different trajectories under different methods. Results are
comparable across methods in distribution, not path by path.

## 5. Takeover, thresholds and stop conditions

A `stop` condition ends a cell the first time it is true after an event.
The cell's termination is `stopped` and `stop_time` is the event time.
Here it is a field of the model dictionary; a model built with `crn.Model`
sets it with `model.set_stop(...)` (section 12).
The two common outputs are:

- **Takeover fraction by age:** the share of cells whose stop time lies
  before that age.
- **Time-to-threshold distribution:** the stop times of the stopped cells.

```python
"""Takeover fraction by age and time to the 90% threshold."""

import json

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
result = crn.ensemble(
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

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
run = crn.ensemble(
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

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
run = crn.ensemble(
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
# Every logged event as a DataArray: cells x entries x observables after the
# event, NaN-padded; time, rule and rule_name are coordinates of each entry.
log = result.event_columns()
names = log["rule_name"].to_series().dropna()  # padding has no name
print("logged mutation events:", names.value_counts().sort_index().to_dict())
print("omitted (log full):", sum(o["omitted_log_entries"] for o in result.observations))
table = result.event_frame(0)  # the log of cell 0 as a DataFrame
print(table[["time", "rule", "rule_name", "N", "mutants"]])
days = result.event_columns(0)["time"] / 24  # cell 0, without padding
print("cell 0 mutation days:", days.round(1).values.tolist())
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

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
run = crn.ensemble(
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
- `result.event_columns()` returns every log as one xarray DataArray with
  dimensions (`replicate`, `entry`, `variable`): the observables after
  each event, NaN past a cell's last entry. Over (`replicate`, `entry`) it
  has the coordinates `time` (NaN as padding), `rule` (the rule index, −1
  for a leap, −2 as padding) and `rule_name` (NaN for a leap and as
  padding); over `replicate`, `length`, the number of entries each cell
  logged (0 for a cell without logged events). `result.event_columns(id)`
  returns one cell's log with dimensions (`entry`, `variable`) and no
  padding.
- `result.event_frame()` and `result.event_frame(id)` return the same
  entries as the DataFrames of release 0.4.0, indexed by (`replicate`,
  `entry`) or by `entry`: `time`, `rule`, `rule_name`, then the
  observables. Participating instances are in `result.observations`.

## 8. Choosing a method

The generic methods run every model within their scope. None of them
depends on mtDNA or on the names in a model:

| `method` | Kind | Scope | Use it for |
| --- | --- | --- | --- |
| `"static_direct"` | exact | closed models (a fixed set of attribute vectors) | the default for most work; Ziggurat waiting times and sparse count trackers; a selector chooses the next event (below) |
| `{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False}` | exact | closed models | the Static Direct of release 0.3.0: inverse waiting times and per-event trackers; reproduces its numbers |
| `"next_reaction"` | exact | closed models | Gibson–Bruck; pays off in large sparse networks, not here, where every event changes every rate |
| `None` (tree) | exact | every model without an ATP block or continuous variables | models whose instances change identity, such as the clone model |
| `{"kind": "tau", "epsilon": 0.03}` | approximate | closed models | adaptive Poisson leaps; at K = 1000 the vacancy limits leaps, so it runs exact events |
| `{"kind": "binomial", "step": 1.0}` | approximate | closed models with single-parent unit births and deaths | fixed hourly splitting, similar to a binomial hourly scheme in Java |
| `"hybrid"` or `{"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}` | approximate at its tolerances | closed models, with or without continuous variables | models with continuous variables (section 10), which no other method runs; without them it runs exact events up to rounding, slower than Static Direct |

Two methods serve one model structure each:

| `method` | Kind | Use it for |
| --- | --- | --- |
| `"mtdna_direct"` | exact | an opt-in fast path for the three-class fixed-population structure, about 2.8× faster than Static Direct; other models fall back to Static Direct with inverse waits (below) |
| `"thinning"` | exact | stored models with the legacy ATP block (section 10), and only those |

The same 64 cells under each method that runs the three-class model:

```python
"""The same 64 cells under every method that runs the three-class model."""

import json
import time

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
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
    run = crn.ensemble(
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
[performance analysis](evidence/mtdna-performance-analysis.md) compares the methods
at full scale.

### Selectors for Static Direct

Static Direct chooses the next event among the channels of the lowered model.
A plain scan is the fastest choice for the eight channels of the three-class
model, and a slow one for a model with hundreds of channels. A *selector*
changes how the channel is chosen, not the process: every selector samples the
same law, uses the random stream in its own way, and puts its name into
`metadata["sampler"]`, so the same seed gives different paths.

| `selector` | How it chooses |
| --- | --- |
| `"linear"` | scans the channels in order |
| `"composition_rejection"` | groups the channels by the binary exponent of their rates, then draws a group and a member |
| `"sorting"` | scans in an order that moves each fired channel one place forward |
| `"rssa_cr"` | draws candidates from upper bounds of the rates over population intervals and accepts them against lower bounds or the exact rate |
| `"auto"` | picks one of the first, three, or `linear`, at preparation (below) |

`crn.ensemble(..., method="static_direct", selector="rssa_cr")` and
`method={"kind": "static_direct", "selector": "rssa_cr"}` are the same request.
`crn.ensemble` leaves the choice to `auto`. It picks RSSA-CR from 8 channels
when no rate reads two or more species, and from 128 channels when the species
that start with molecules average fewer than 8; otherwise composition–rejection
from 128 channels, and otherwise the linear scan. It decides from the initial
counts, so a model whose populations later grow or shrink far from them can end
up with a slower selector than another choice would give.
`crn.describe()["selectors"]["auto"]` lists the thresholds. The three-class
model has eight channels whose rates read the molecule counts, so `auto` runs
the linear scan there, and every example above gives the same numbers as
before.

**`crn_craft.rules` keeps the linear scan.** A request that `crn.rules.ensemble`
writes always scans, as in every earlier release, so adding `quantiles=` to a
0.6.0 script never changes its paths. Only `crn.ensemble` leaves the choice to
`auto`, and `selector=` overrides it in either module. `crn.prepare` of a
`rules-run` request keeps that request's linear scan.

The example builds twenty strains that compete for room. Each strain grows into
the free room, dies, and gives rise to the next strain. Every rate reads all
the strains, so `auto` runs the linear scan on these 59 channels:

```python
"""Selectors on a twenty-strain competition model."""

import time

import crn_craft as crn


def strains(count, builder=crn.Model):
    """Strains S0, S1, ... share room K; strain i gives rise to strain i + 1."""
    names = [f"S{i}" for i in range(count)]
    model = builder(f"strains_{count}")
    model.add_parameter({"K": 500.0, "g": 0.02, "mu": 0.001})
    model.add_species({name: 200 if i == 0 else 0 for i, name in enumerate(names)})
    free = f"max(0, K - ({' + '.join(names)})) / K"
    for i, name in enumerate(names):
        model.add_reaction(f"g * {name} * {free}", [name], [name, name], f"grow_{i}")
        model.add_reaction(0.01, [name], [], f"die_{i}")
        if i + 1 < count:
            model.add_reaction(f"mu * {name}", [name], [name, names[i + 1]], f"mutate_{i}")
    return model


model = strains(20)
for selector in ("linear", "composition_rejection", "sorting", "rssa_cr", "auto"):
    run = crn.ensemble(
        model,
        end=300.0,
        report_times=[0.0, 300.0],
        replicates=200,
        seed=1,
        threads=2,
        method="static_direct",
        selector=selector,
    )
    start = time.perf_counter()
    result = run.simulate()
    seconds = time.perf_counter() - start
    total = result.statistics["mean"].iloc[-1].sum()
    print(
        f"{selector:22s} {seconds:5.2f} s  {result.metadata['events']:>9,} events  "
        f"molecules {total:6.1f}  {result.metadata['sampler']}"
    )
# The 0.6.0 module scans, whatever the model's size.
old = crn.rules.ensemble(
    model, end=300.0, report_times=[0.0, 300.0], replicates=200, seed=1,
    threads=2, method="static_direct",
).simulate()
print("crn.rules:", old.metadata["sampler"])
```

```text
linear                  0.10 s    312,169 events  molecules  270.5  chacha8-u53-ziggurat/rules-static-direct/1
composition_rejection   0.14 s    313,771 events  molecules  270.8  chacha8-u53-ziggurat/rules-static-direct-cr/1
sorting                 0.13 s    312,179 events  molecules  270.7  chacha8-u53-ziggurat/rules-static-direct-sorting/1
rssa_cr                 0.20 s    313,771 events  molecules  270.8  chacha8-u53-ziggurat/rules-static-direct-rssa-cr/1
auto                    0.10 s    312,169 events  molecules  270.5  chacha8-u53-ziggurat/rules-static-direct/1
crn.rules: chacha8-u53-ziggurat/rules-static-direct/1
```

The molecule counts agree within sampling error, and each selector names itself
in the sampler. Every rate here reads all strains, so RSSA-CR pins its bounds to
the current counts and draws the same events as composition–rejection. With 59
channels none is faster than the scan. Selectors pay off on larger models: see
[Larger models](#larger-models). Generated event loops (below) run only the
linear scan.

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
division and reads its stop from the counts. On the example of section 3
that is about 2.8× less time per event than Static Direct with its
default Ziggurat waits and sparse trackers, and 3.1× less than Static
Direct with inverse waits and per-event trackers (medians of three runs on
the two-CPU host, 12 November 2026, `validation/merge/manual-ratios/`; the
[performance analysis](evidence/mtdna-performance-analysis.md#2-methods) has older
measurements).

The structure need not be called mtDNA. Below, three strains share `room`
places, and the resident strain gives rise to two variants. Every name
differs from the model file, and `mtdna_direct` still certifies the model.
A strain that dies faster breaks check 2:

```python
"""mtdna_direct recognizes a model structure, not a file or names."""

import copy
import time

import crn_craft as crn

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
    result = crn.ensemble(
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

*Opt-in, off by default; needs a Rust toolchain. New in 0.4.0; string
rates since 0.5.0.*

`crn.compile_loop(model)` compiles the whole Static Direct event of one
model into a small native program; `compiled=loop`, on `crn.ensemble` or
on `simulate()`, then runs Static Direct requests in it. Results are byte-identical to the interpreter's: the same
records, statistics, paths, trackers, event logs and sampler identity, so
replays and stored results keep their meaning. Unlike `mtdna_direct`, it
keeps Static Direct's use of the random stream and is not limited to the
three-class structure.

```python
"""Static Direct with a compiled event loop (needs cargo)."""

import json
import time

import crn_craft as crn

document = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
document["stop"] = "N > 0 && mutants >= 0.9 * N"
model = crn.Model.from_dict(document)  # model/1.0, the form crn.ensemble runs
loop = crn.compile_loop(model)  # about a minute the first time, then cached
run = crn.ensemble(
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
```

```text
interpreted   1.19 s   28,523,982 events  taken over 20  chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1
compiled      0.51 s   28,523,982 events  taken over 20  chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1
byte-identical
```

On the 1000-cell workload of the
[performance analysis](evidence/mtdna-performance-analysis.md#2-methods), timed
from Python including the subprocess, the compiled loop took 2.6 s
against 5.8 s for interpreted Static Direct with its defaults (2.21×),
and 3.8 s against 7.0 s with inverse waits and per-event trackers (1.86×),
medians of three on the two-CPU host
(`validation/mtdna/compiled-loop-python.json`).
`mtdna_direct` is faster still where it applies.

`compile_loop` also compiles string rates. A string rate in `crn.Model`
is a whole-rule propensity, and the generator emits it, like any rate
expression, with the interpreter's operations in the interpreter's order:
the same `libm` functions, NaN propagation through comparisons, logic and
`if`, short-circuit `&&` and `||`, and `min` and `max` that keep NaN. The
results stay byte-identical to interpreted Static Direct, which `verify`
checks:

```python
"""String rates in a compiled loop: the builder's SIR model (needs cargo)."""

import time

import crn_craft as crn

sir = crn.Model("sir")
sir.add_species({"S": 9_990, "I": 10, "R": 0})
sir.add_parameter({"beta": 0.3, "gamma": 0.1, "N": 10_000.0})
sir.add_reaction("beta * S * I / N", ["S", "I"], ["I", "I"], name="infection")
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
loop = crn.compile_loop(sir)  # string rates compile like constant ones
run = crn.ensemble(
    sir,
    end=400.0,
    report_times=[0.0, 100.0, 200.0, 400.0],
    replicates=500,
    seed=3,
    threads=2,
    method="static_direct",
)
for label, compiled in (("interpreted", None), ("compiled", loop)):
    start = time.perf_counter()
    result = run.simulate(compiled=compiled)
    seconds = time.perf_counter() - start
    print(f"{label:12s} {seconds:5.2f} s  {result.metadata['events']:>11,} events")
loop.verify(run)
print("byte-identical")

# Without cargo: the generator still refuses models only another method runs.
atp = crn.load("models/rules/mtdna-three-class-atp-continuous.rules.json")
try:
    crn.loop_source(atp)
except crn.CrnError as error:
    print(error.diagnostics[0]["code"], error.diagnostics[0]["message"])
```

```text
interpreted   0.61 s    9,400,096 events
compiled      0.22 s    9,400,096 events
byte-identical
CRN-R020 continuous variables require the hybrid method
```

On the same SIR model with 2,000 replicates (37.6 million events), timed
from Python including the subprocess, the compiled loop took 1.64 s
against 3.35 s with inverse waits (2.00×) and 1.36 s against 2.97 s with
Ziggurat waits (2.18×), medians of three on the two-CPU host
(`validation/rules/compiled-loop-string-rates.json`).

**`compiled=` on `crn.ensemble` or on `simulate()`.** A loop can be
given in two places, and both run the same program with the same results:

- `crn.ensemble(..., compiled=loop)` (or `crn.prepare(model, request,
  compiled=loop)`) attaches the loop to the prepared run. Preparation
  checks once that the loop was compiled for these model bytes and that
  the request selects `static_direct`, and raises `crn.CrnError` `CRN-R020`
  (at `/model` or `/method`) otherwise, before anything runs. Every
  `run.simulate()` then uses the loop without further arguments. Use it
  when all runs of a request should be compiled: repeated calls, or code
  such as a fitting loop that only calls `simulate()`.
- `run.simulate(compiled=loop)` uses the loop for that call only, on a run
  prepared without one; the first such call makes the same checks. The
  same prepared run can then also run interpreted, as the example above
  does to compare the two.

A run with an attached loop always runs compiled: `simulate(callback=...)`
raises `ValueError`, so prepare the request without `compiled` to watch
events. `run.bind_parameters(...)` returns a candidate without the loop,
which runs interpreted, since a loop serves only the parameter set it was
compiled for.

**Requirements.** Rust 1.85 or newer, the MSRV of the crates (install it
with [rustup](https://rustup.rs)), and network access to crates.io: the
program depends on `crn-craft` at exactly the installed version.
Offline, or with a version that is not on crates.io (a wheel built from a
clone), pass `source=` a CRN-craft checkout of that version, or set
`CRN_CRAFT_SOURCE`. `cargo=` names another cargo executable.

**Build time and cache.** The first build compiles `crn-craft` and
its dependencies: about 40 s on the two-CPU host after the downloads, several
minutes on slower machines. `compile_loop` prints cargo's progress and the
path of the build log to stderr (`progress=False` silences it). Builds are
cached in the user cache directory (`~/.cache/crn-craft/rules-loops` on
Linux; `cache_dir=` or `CRN_CRAFT_CACHE_DIR` choose another), keyed by the
loop's fingerprint, the installed version, `rustc -vV`, the target and
Rust flags. A cached call took 0.06 s. All builds share one Cargo target
directory (`CARGO_TARGET_DIR` if set), so a second model built in 3 s.
The cache takes about 160 MB; delete the directory to reclaim it. Loops cached
by release 0.6.0 are built again once on first use, with the same results, and
their entries stay in the cache until you delete it.

**Limitations.**

- Only `method="static_direct"`, with any options, and only the linear scan:
  `crn.ensemble` leaves the selector to `auto`, so for a model whose `auto`
  choice is another selector, name `selector="linear"` together with
  `compiled=` (`CRN-R020` at `/method/selector` otherwise). Other methods are
  refused with `CRN-R020` at `/method`.
- The generator emits every model that `method="static_direct"` runs:
  constant, vacancy and expression rates, string rates (`propensity`
  rules) and `sum`, `distinct`, `max` and `min` observables. It refuses,
  with `CRN-R020` at `/method` and before cargo runs, models that only
  another method runs (continuous variables, ATP blocks, instances that
  change identity), and loops that would exceed 1 MiB of source, which
  densely coupled networks with many channels reach.
  `crn.loop_source(model)` returns the source, or that error, without
  cargo.
- The loop is compiled for the exact model bytes. `crn.compile_loop`
  builds for the document `crn.ensemble` runs (`model/1.0`), so the same
  dictionary or `crn.Model` works for both calls. `crn.rules.compile_loop`
  builds for the bytes as given, which `crn.rules.ensemble` runs; a loop
  compiled for other bytes is refused with `CRN-R020` at `/model`.
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
  ([details](developer/specialized-kernels.md#generated-rules-event-loops)).
- Missing cargo, a failed build (`error.log` is the build log) and a
  result that differs raise `crn.CompiledLoopError`, whose `reason` is
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
import crn_craft as crn

# Every new mutant is its own clone (attribute `clone` from new_id()).
model = json.loads(crn.load("models/rules/mtdna-fixed.rules.json"))
model["parameters"].update({"K": 100, "w0": 100, "mu": 0.01, "delta": 0.05})
run = crn.ensemble(
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
import crn_craft as crn

data = crn.load("models/rules/mtdna-three-class-atp-continuous.rules.json")
model = json.loads(data)
print("continuous:", model["continuous"])
print("replication_wt rate:", model["rules"][1]["rate"])
run = crn.ensemble(
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
request is `run/1.0` (`rules-run/0.5` under `crn.rules.ensemble` for a
`rules-model` document), and the result records the `method` and the
`continuous` column names. `validation/check_hybrid.py`
checks `hybrid` against exact laws and against exact thinning on this model
([results](../validation/hybrid/RESULTS.md)). On the test host it took about
eleven times as long as thinning on the ATP model, since every event needs an
integration step. The error has not yet been measured as a function of the
tolerances.

**Your own variables.** `crn.Model.add_ode(name, initial, derivative, *,
independent=False, lower_bound=None, on_bound="stop")` adds a variable to a
builder model (section 12); string rates may read it. `independent=True`
marks an ODE that reads only parameters and other independent variables
(from `rules-model/0.4` on). `add_ODE` is the same method, and `add_continuous`, its
name before 0.5.0, still works with the same arguments but emits a
`DeprecationWarning`:

```python
"""A continuous variable in a builder model."""

import crn_craft as crn

cell = crn.Model("atp_wild_type")
cell.add_parameter({"f": 0.1, "v1": 0.2, "v2": 0.01, "c": 1.3034153747670567})
cell.add_species("W", 1000)
cell.add_ode("ATP", 450.0, "f * W - v1 * ATP - v2 * W", lower_bound=0)
cell.add_reaction("c / (c + ATP) * W", ["W"], ["W", "W"], name="replication")
cell.add_reaction(0.0028881132523331052, ["W"], [], name="degradation")
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

### Independent ODEs

*New in 0.5.0 (`rules-model/0.4`).*

Some continuous variables do not depend on the molecules at all: a forcing
such as a daily rhythm, a drug concentration that decays on its own, or a
hazard that grows with age. Declare such a variable with
`add_ode(..., independent=True)`, which writes `"independent": true`. Its `initial` and `derivative` may then read only
parameters and other independent variables; anything else is refused with
`CRN-R012` at the derivative's path. Rates, propensities and the
derivatives of the other variables read it as before, and results report it
like any continuous variable.

`hybrid` integrates the independent variables once per prepared run, over
the whole horizon at the request's `rtol` and `atol`, and every replicate
and thread reads them from that one solution, so results do not depend on
`threads`. Each replicate integrates only the other variables and the
integrated propensity. A large or stiff independent system is then no
longer integrated again in every replicate: on an eight-compartment pulse
driving arrivals, `independent=True` was 2.8 times as fast
(`validation/hybrid/RESULTS.md`, section E). When events dominate, as with
a two-variable daily rhythm on the three-class model
(`examples/mtdna/14_independent_ode.py`), the run time hardly changes, but
the rhythm stays more accurate because it is not restarted at every event.

Bounds keep their meaning. With `on_bound="clamp"` the variable is held at
its bound. With `"stop"`, reaching the bound ends the replicate, and since
the variable is the same in every replicate, every replicate that is still
running stops at the same time. The shared solution is computed while the
request is prepared, so `cancel_file` and `wall_time_s` do not interrupt it;
the step budget `max_sampling_steps` bounds it. Without `independent`, a
rules-model/0.4 document means what it means as 0.3.

```python
"""An independent ODE: Gompertz mortality (hybrid)."""

import math

import crn_craft as crn

# The hazard m grows as dm/dt = g m whatever the molecules do.
gompertz = crn.Model("gompertz")
gompertz.add_parameter({"g": 0.5, "m0": 0.1})
gompertz.add_species({"A": 20})
gompertz.add_ode("m", "m0", "g * m", independent=True)
gompertz.add_reaction("m * A", ["A"], [], name="death")
print("schema:", gompertz.to_dict()["schema"])
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
model = gompertz.to_dict()
model["continuous"][0]["derivative"] = "g * m * A"
try:
    crn.ensemble(model, end=4.0, report_times=[4.0], method="hybrid")
except crn.CrnError as error:
    print(error.diagnostics[0]["code"], error.diagnostics[0]["path"])
    print(error.diagnostics[0]["message"])
```

```text
schema: crn-craft.model/1.0
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

Each molecule survives to time t with probability exp(−∫m), so the count
at t = 4 is binomial with mean 20 × 0.2786 = 5.573; 4,000 replicates give
5.585.

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

import crn_craft as crn

data = crn.load("models/rules/mtdna-three-class-atp.rules.json")
print("ATP block:", json.loads(data)["atp"])
result = crn.ensemble(
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

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
template = crn.ensemble(
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

`crn.Model` writes a model from Python calls, in the style of
rebop's `Gillespie` class. Every species is an attribute-free type, and
every species count is an observable with the species name.

```python
"""An SIR epidemic written with crn.Model."""

import crn_craft as crn

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
| `crn.Model(name="model")` | a new `crn-craft.model/1.0` document (`crn.rules.Model` writes `rules-model/0.3`, or 0.4 with an independent ODE) |
| `add_parameter(name, value)` or `add_parameter(mapping)` | parameters that string rates may read; `add_parameters` is the same method |
| `add_species(name, initial=0)` or `add_species(mapping)` | species types, their count observables and their initial counts |
| `set_initial(mapping)` or `set_initial(S=490, I=10)` | initial counts of declared species |
| `add_reaction(rate, reactants, products, name=None)` | a rule; species are lists of names, repeated for multiplicity; undeclared species start at 0; names default to `r0`, `r1`, ... |
| `add_ode(name, initial, derivative, *, independent=False, lower_bound=None, on_bound="stop")` | a continuous variable (section 10); the other arguments are keyword-only. `independent=True` marks an ODE that reads only parameters and other independent variables. `add_ODE` is the same method; `add_continuous`, the name before 0.5.0, warns |
| `set_stop(condition)` | the model's stop condition, an expression over parameters and species counts such as `"I == 0"`; `None` removes it |
| `to_dict()`, `to_json()` | the JSON document |

A mapping, as in `add_parameter({"gamma": 0.1, "beta": 5e-4})`, declares
its entries in order with the checks of one call each. If an entry is
invalid, the call raises and declares none of them. A mapping takes no
second argument: `add_species({"S": 490}, 10)` is a `TypeError`.

`crn.ensemble` and `crn.prepare` accept the `Model` itself, and `Model.from_dict`
builds one from a model dictionary or model bytes. Every method
runs builder models; `tree` is never needed, since the species have no
attributes.

**Stop conditions.** A `Model` is not a dictionary, so `model["stop"] = ...`
raises `TypeError: 'Model' object does not support item assignment`. Use
`set_stop`, which writes the same `stop` field (section 5): once the
condition is true, checked at time 0 and after every event, the cell keeps
its state to the end of the report grid, its termination is `stopped` and
`stop_time` is the event time. The condition reads parameters and species
counts, not continuous variables (a `lower_bound` of `add_ode` stops on
those), and `set_stop(None)` removes it. `model.to_dict()` returns a
dictionary that can be edited like a model file, for fields the builder
does not cover:

```python
"""A stop condition in a builder model: the end of an SIR epidemic."""

import crn_craft as crn

sir = crn.Model("sir")
sir.add_species({"S": 490, "I": 10, "R": 0})
sir.add_parameter({"beta": 5e-4, "gamma": 0.1})
sir.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
sir.set_stop("I == 0")  # the epidemic is over
result = crn.ensemble(
    sir,
    end=1000.0,
    report_times=[0.0, 1000.0],
    replicates=1000,
    seed=1,
    threads=2,
    method="static_direct",
).simulate()
print(result.records["termination"].value_counts().to_dict())
print(f"mean end of the epidemic: day {result.records['stop_time'].mean():.1f}")
print(result.statistics["mean"].round(1))
```

```text
{'stopped': 1000, 'completed': 0, 'event_limit': 0, 'sampling_limit': 0, 'instance_limit': 0, 'wall_time_limit': 0, 'cancelled': 0, 'failed': 0, 'output_limit': 0}
mean end of the epidemic: day 110.5
           I      R      S
time                      
0.0     10.0    0.0  490.0
1000.0   0.0  446.1   53.9
```

A stopped cell keeps its state, so the means at day 1000 are the final
sizes of the epidemics.

**Rates.** A numeric rate k is a mass-action constant, counted like the rest
of CRN-craft: k times the number of distinct combinations of reactant
molecules. `A` fires at k·A, `A + B` at k·A·B and `A + A` at k·A(A − 1)/2.
rebop drops the 1/2 for identical reactants, so its rate k on `["A", "A"]`
equals `2 * k` here; reactions without identical reactants agree. A string
rate is the reaction's whole propensity (the `propensity` field of the
model), an expression over parameters and species counts with the operators
`+ - * / ^`, comparisons and the functions `min max exp log sqrt abs floor
ceil if`. It applies while each reactant has enough molecules (two for
`["A", "A"]`) and is 0 otherwise, so `"vmax * S / (Km + S)"` on `["S"]` stops
when `S` is used up. A reaction may take any number of reactant molecules:
`["A", "A", "A"]` is a third-order reaction, written as a pattern with
`count: 3` (the model reference describes
[higher-order rules](reference/model-language.md#higher-order-rules-model10)).
`crn.rules.Model` keeps the limit of two.
`validation/check_builder.py` compares builder models with exact means and
with rebop 0.10.0 ([results](../validation/builder/RESULTS.md)).

The three-class model can be written the same way, with one species per
class and the vacancy law in the propensities. `examples/mtdna/12_builder.py`
builds it and compares it with the JSON file:

```python
"""The three-class model built with crn.Model."""

import crn_craft as crn

model = crn.Model("mtdna_three_species")
parameters = {"K": 1000, "delta": 0.0028881132523331052, "sw": 1.0, "sm": 1.5}
parameters |= {"mu": 3e-05, "fracAdv": 0.2}
model.add_parameter(parameters)
model.add_species({"W": 1000, "M": 0, "A": 0})
vacancy = "max(0, K - W - M - A) / (W + sm * (M + A))"
model.add_reaction(f"{vacancy} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt")
model.add_reaction(f"{vacancy} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt")
model.add_reaction(f"{vacancy} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv")
model.add_reaction(f"{vacancy} * sw * M", ["M"], ["M", "M"], "replication_mt")
model.add_reaction(f"{vacancy} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
for species in ("W", "M", "A"):
    model.add_reaction(parameters["delta"], [species], [], f"degradation_{species}")
result = crn.ensemble(
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
(section 2 and the [rules language](reference/model-language.md)).

### Larger models

A model written with `crn.Model` is `crn-craft.model/1.0`, which accepts 4,096
species types, 16,384 parameters, 4,096 observables and 65,536 rules, sixteen
times the limits of the `rules-model` versions (64, 1,024, 256 and 1,024). The
static methods accept 4,096 instance tuples and 65,536 channels. To lift the
limits of a model file, change its `schema` string to `crn-craft.model/1.0`:
the two mean the same thing. The next block builds the competition model of
section 8 with a hundred strains, which `crn.rules.Model` refuses. With 299
channels, `auto` runs composition–rejection:

```python
"""A hundred strains: beyond the limits of the 0.6.0 model, run with auto."""

import crn_craft as crn


def strains(count, builder=crn.Model):
    """Strains S0, S1, ... share room K; strain i gives rise to strain i + 1."""
    names = [f"S{i}" for i in range(count)]
    model = builder(f"strains_{count}")
    model.add_parameter({"K": 500.0, "g": 0.02, "mu": 0.001})
    model.add_species({name: 200 if i == 0 else 0 for i, name in enumerate(names)})
    free = f"max(0, K - ({' + '.join(names)})) / K"
    for i, name in enumerate(names):
        model.add_reaction(f"g * {name} * {free}", [name], [name, name], f"grow_{i}")
        model.add_reaction(0.01, [name], [], f"die_{i}")
        if i + 1 < count:
            model.add_reaction(f"mu * {name}", [name], [name, names[i + 1]], f"mutate_{i}")
    return model


settings = dict(end=100.0, report_times=[0.0, 100.0], replicates=20, seed=1, threads=2)
old = strains(100, crn.rules.Model)
print("rules.Model schema:", old.to_dict()["schema"])
try:
    crn.rules.ensemble(old, method="static_direct", **settings)
except crn.CrnError as error:
    print("refused:", error.diagnostics[0]["code"], error.diagnostics[0]["path"])
model = strains(100)
print("Model schema:", model.to_dict()["schema"])
result = crn.ensemble(model, method="static_direct", **settings).simulate()
print(result.metadata["sampler"], result.metadata["termination"])
print(f"{result.metadata['events']:,} events;", len(result.rule_names), "rules")
```

```text
rules.Model schema: crn-craft.rules-model/0.3
refused: CRN-R015 /types
Model schema: crn-craft.model/1.0
chacha8-u53-ziggurat/rules-static-direct-cr/1 completed
10,163 events; 299 rules
```

## 13. Model variants

The model file is plain JSON. Load it with `json.loads(crn.load(path))`,
change it, and pass the dictionary to `crn.ensemble`.

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

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
model["parameters"].update({"K": 500, "w0": 400, "delta": math.log(2) / (20 * 24)})
# Start heteroplasmic: 400 wild type and 100 neutral mutants (class 1).
model["initial"].append(
    {"type": "mtDNA", "attributes": {"clone": 1, "s": "sw", "wt": 0}, "count": 100}
)
model["stop"] = "N > 0 && mutants >= 0.6 * N"  # 60% threshold
run = crn.ensemble(
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
result = crn.ensemble(
    odd, end=100.0, report_times=[0.0, 100.0], replicates=2, method="mtdna_direct"
).simulate()
print(result.metadata["sampler"])
print([d["message"] for d in result.metadata["diagnostics"]])

# Mistakes are reported with a code and a JSON path before anything runs.
broken = json.loads(json.dumps(model))
broken["rules"][1]["rate"] = "max(0, K - N) / DD * sw"
try:
    crn.ensemble(broken, end=1.0, report_times=[0.0, 1.0])
except crn.CrnError as error:
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
- `crn.load(path)` checks a model file before any run;
- `crn-craft check MODEL [REQUEST]` does the same from the command line, and
  with a request it checks that `run` would accept both.

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

import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))


def summary(result):
    terms = sorted(set(result.records["termination"]))
    codes = [d["code"] for d in result.metadata["diagnostics"]]
    return result.metadata["termination"], terms, codes


# Too few events for three years: every cell ends with event_limit and the
# ensemble withholds its statistics.
small = crn.ensemble(
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
slow = crn.ensemble(
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

Usage: python SCRIPT [DIRECTORY [CRN_CRAFT_BINARY]]

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
    "schema": "crn-craft.run/1.0",
    # The request binds the model by relative path and exact SHA-256.
    "model": {
        "path": "takeover.rules.json",
        "sha256": hashlib.sha256(data).hexdigest(),
    },
    "time": {"end": 25920.0, "report_times": [0.0, 8640.0, 17280.0, 25920.0]},
    "ensemble": {"replicates": 1000, "root_seed": 11},
    "outputs": {
        "statistics": ["mean", "variance"],
        "quantiles": [0.5],
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
        [sys.argv[2], "run", "takeover.rules.json", "takeover.run.json"],
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

These commands run with the native CLI (`cargo install crn-craft-cli --locked`)
or with the `crn-craft` command of the Python package, which has the same
commands and reports (section 1). `rules-check` and `rules` are aliases of
`check` and `run`.

```sh
python examples/mtdna/11_cli_request.py /tmp/takeover
crn-craft check /tmp/takeover/takeover.rules.json /tmp/takeover/takeover.run.json
crn-craft run /tmp/takeover/takeover.rules.json /tmp/takeover/takeover.run.json > result.json
crn-craft run /tmp/takeover/takeover.rules.json /tmp/takeover/takeover.run.json --archive /tmp/takeover/archive
crn-craft verify /tmp/takeover/archive
crn-craft describe run-1.0       # the JSON Schema of the request
crn-craft describe result-1.0    # and of the result
```

`check` reports whether `run` would accept the two documents, without running
them. `run --archive NEW_DIRECTORY` also writes the model, the request and the
result with a manifest (`crn-craft.archive/1.0`), and `verify` checks such a
directory later; from Python, `run.archive(directory)` and `crn.verify_archive`
do the same. See the [archive contract](reference/archives.md).

```text
model sha256: f9eaaa59e4e449a55bd991350c748fdd8118b7c79dd81d6bcb0fc57bf025c0c2
completed chacha8-u53-ziggurat/rules-mtdna-direct/1
cells: 1000 stopped: 107
events: 142356193
result fields: ['completed_replicates', 'diagnostics', 'events', 'final_instances', 'method', 'model_sha256', 'observables', 'observations', 'paths', 'replicates', 'report_times', 'requested_replicates', 'root_seed', 'rules', 'run_sha256', 'sampler', 'schema', 'statistics', 'stopped_replicates', 'termination', 'types']
```

The request versions decide which fields are allowed:

| Version | Adds |
| --- | --- |
| `crn-craft.run/1.0` | `method.selector` for `static_direct`, `outputs.quantiles`, the method `fsp`; Static Direct defaults to the selector `auto`; the result is `crn-craft.result/1.0` |
| `crn-craft.rules-run/0.1` | the tree method only |
| `crn-craft.rules-run/0.2` | `method`, `trackers`, `event_log` |
| `crn-craft.rules-run/0.3` | `mtdna_direct` |
| `crn-craft.rules-run/0.4` | Static Direct `waiting_time` and `sparse_trackers` |
| `crn-craft.rules-run/0.5` | `hybrid` with its `rtol` and `atol`; the result is rules-result/0.3, which adds `method` and `continuous` |
| `crn-craft.rules-run/0.6` | no field; Static Direct defaults to Ziggurat waits and sparse trackers; the result is rules-result/0.3 |

`crn.ensemble` always writes `run/1.0`. `crn.rules.ensemble` picks the lowest
`rules-run` version that fits its arguments, and rules-run/0.6 for
`static_direct`; it writes `run/1.0` only when an argument needs it
(`quantiles=`, `selector=`) or the model is `crn-craft.model/1.0`.
`crn.prepare(model, request)` takes a request you wrote yourself. Stored
requests of every earlier version keep their meaning and replay
byte-identically.

Model files carry their own version: `rules-model/0.1` is the base language,
0.2 adds the legacy `atp` block, 0.3 adds whole-rule `propensity` and
`continuous` variables, and 0.4 adds independent variables. `crn-craft.model/1.0`
means what 0.4 means, with larger limits and the pattern `count`.
`crn.Model` writes `model/1.0`, and `crn.rules.Model` writes 0.3 (0.4 with an
independent ODE). A `rules-model` file runs under either request family.

## 16. Rust

The crate `crn-craft` is the engine behind Python and the CLI; its module
`crn_craft::rules` holds the functions below. `prepare` validates exact model
and request bytes, and `execute` runs them. The crate `crn-craft-rules` is a
deprecated facade that re-exports the module, so `crn_craft_rules::X` still
compiles:

```rust
use crn_craft::rules::{Control, execute, prepare};

let prepared = prepare(&model_bytes, &request_bytes)?;
let result = execute(&prepared, &Control::default())?;
println!("{} cells stopped", result.stopped_replicates);
```

- `execute_with_observer(&prepared, &control, &observer)` calls an
  `EventObserver` for every committed event on the native workers. It
  returns `ObserverAction::Stop` to end that cell. An observer that always
  continues compiles away and leaves results byte-identical.
- *Opt-in:* `generate_rules_loop(&prepared)` emits Rust source
  for the whole Static Direct event of one prepared model, for every model
  that Static Direct runs (expression and string rates included). It
  requires a `static_direct` request and at most 1 MiB of source, and
  refuses anything else with `CRN-R020`.
  Parameters are compiled in, so another parameter set needs another loop.
  Compile it into your program and run it with `execute_with_rules_loop`,
  which runs only a loop generated for the same model bytes. Results are
  byte-identical to `execute`, and on the example of section 3 it is about
  2.1× faster than interpreted Static Direct (timed from Python with
  `compile_loop`, 12 November 2026). On networks with many
  constant-rate channels it can be slower, so measure first
  ([details](developer/specialized-kernels.md#generated-rules-event-loops)). From
  Python, `crn.compile_loop` builds such a loop into a separate program
  and runs it (section 8).

## 17. Checking results

Exactness is checked against independent references, and the checks can be
rerun:

- **Quantiles and selectors.** `validation/merge/check_quantiles.py` checks
  the quantile statistics against NumPy, bit for bit, and
  `validation/merge/check_selectors.py` checks every selector against the
  exact laws, with declared tests (`validation/merge/RESULTS.md`).
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

`crn.ensemble` arguments (`crn.rules.ensemble` has the same, with the defaults
of release 0.6.0):

| Argument | Meaning |
| --- | --- |
| `model` | model bytes from `crn.load`, a dictionary, or a `crn.Model` (section 12) |
| `end`, `report_times` | horizon in hours; increasing report times in [0, end] |
| `replicates`, `seed`, `threads` | number of cells, root seed, 1–64 native threads |
| `keep_paths` | how many cells keep their full observable series (`result.paths`); `"all"` or `True` for every cell |
| `keep_instances` | how many cells keep their final instances (`metadata["final_instances"]`) |
| `method` | a method name or dictionary (section 8) |
| `selector` | the Static Direct selector (section 8); absent means `auto` for `crn.ensemble` and the linear scan for `crn.rules.ensemble` |
| `quantiles` | probabilities for `statistics["quantiles"]` and `["median"]`: 0.025, 0.5 and 0.975 by default in `crn.ensemble`, none in `crn.rules.ensemble`; `()` asks for none (section 3) |
| `compiled` | a loop from `crn.compile_loop` for this model, checked here and run by every `simulate()`; `simulate(compiled=loop)` uses one for a single call instead (section 8) |
| `trackers` | list of `{"observable", "first_passage"?, "keep_episodes"?}` |
| `event_log` | `{"max_entries": n, "rules": [...]}` |
| budgets | `max_events`, `max_sampling_steps`, `max_instances`, `max_output_bytes`, `wall_time_s` |

The result, `run.simulate()`:

| Attribute | Content |
| --- | --- |
| `records` | DataFrame, one row per cell, index `replicate`: `started` (bool), `termination` (categorical), `stop_time` (float), `events`, `rng_words`, `peak_instances`, `created_instances`, `saved_points` (unsigned integers) |
| `rule_events` | DataFrame, one row per cell, one column per rule name: events of each rule |
| `statistics` | `None` if withheld; else `samples` and the DataFrames `mean`, `sample_variance`, `mean_standard_error`, index `time`, one column per observable; with quantiles also `quantiles` (a dictionary from probability to DataFrame) and, for 0.5, `median` |
| `model`, `request` | the model and request documents that ran, as new dictionaries |
| `save(path, *, paths=True)` | writes a markdown (`.md`) or JSON (`.json`) report; see `crn.load_report` (section 3) |
| `paths` | DataArray of the kept series, dimensions (`replicate`, `time`, `variable`), NaN after a budget stop; `attrs` with `sampler`, `method`, `root_seed`, `model_sha256`, `run_sha256` |
| `paths_frame()` | DataFrame of the kept series, index (`replicate`, `time`), one column per observable (release 0.4.0's `paths`) |
| `to_xarray()` | Dataset of the kept series, dimensions `replicate` and `time`, one variable per observable |
| `tracker_columns(i)` | DataFrame of tracker `i`, one row per observed cell; NaN where a time or peak is absent |
| `event_columns()`, `event_columns(id)` | DataArray of every log, dimensions (`replicate`, `entry`, `variable`), NaN-padded, with the coordinates `time`, `rule`, `rule_name` and `length`; or of one cell's log, dimensions (`entry`, `variable`) |
| `event_frame()`, `event_frame(id)` | DataFrame of every logged event, index (`replicate`, `entry`), or of one cell's log, index `entry` (release 0.4.0's `event_columns`) |
| `rule_names` | rule names in declaration order, also on the prepared run |
| `observations`, `record_diagnostics`, `metadata` | trackers and logs as documents, per-cell errors by id, the full result document |

Every DataFrame, column and DataArray has `.to_numpy()`. Code written for
release 0.6.0 with `from crn_craft import rules` runs unchanged: the same
names, defaults, documents, hashes and numbers (CI runs the 0.6.0 manual to
check it). Code that moves to the top-level names changes as follows:

| `crn.rules` (0.6.0) | `crn` (0.7.0) |
| --- | --- |
| `rules.Model` writes `rules-model/0.3` | `crn.Model` writes `crn-craft.model/1.0`: larger limits, three or more reactant molecules |
| `rules.ensemble` writes the lowest `rules-run` version, and Static Direct scans | `crn.ensemble` writes `crn-craft.run/1.0`; Static Direct uses the selector `auto`, so its paths can differ for models with many channels (`selector="linear"` keeps them) |
| `statistics` has `mean`, `sample_variance` and `mean_standard_error` | also `median` and `quantiles` by default |
| `result.metadata["schema"]` is `crn-craft.rules-result/0.x` | `crn-craft.result/1.0`, with new hashes |
| `crn_craft.rules.describe()` | `crn.describe()`: `crn-craft.describe/1.0`; the 0.6.0 fields are under `legacy.rules` |
| the network functions of the top level (`run_spec`, `simulate`, `ssa_spec`, `simulate_ssa`, `prepare_ssa`, `prepare_parameters`, `tau`, `fsp`, `analyze` on network documents) | removed: `crn.ensemble`, `Model.fsp(...)`, `Model.analyze()`; the changelog lists each |

Code written for release 0.5.0 changes as follows:

| 0.5.0 | Now |
| --- | --- |
| `result.paths` with dimensions (`time`, `variable`, `replicate`) | dimensions (`replicate`, `time`, `variable`); `result.paths.transpose("time", "variable", "replicate")` gives the old order |
| `result.paths.values[t, v, r]` | `result.paths.values[r, t, v]`; code that selects by name (`sel`, `mean("replicate")`) is unchanged |

Code written for release 0.4.0 changes as follows:

| 0.4.0 | Now |
| --- | --- |
| `result.paths`, a DataFrame | `result.paths_frame()`; `result.paths` is a DataArray |
| `result.paths.loc[3]` | `result.paths.sel(replicate=3)`, or `result.paths_frame().loc[3]` |
| `result.paths["mutants"].unstack("time")` | `result.paths.sel(variable="mutants").to_pandas()` |
| `result.event_columns()`, `result.event_columns(id)`, DataFrames | `result.event_frame()`, `result.event_frame(id)`; `event_columns` returns DataArrays |
| `pip install "crn-craft[xarray]"` | `pip install crn-craft`: xarray is required; the extra still installs |
| `Model.add_continuous(...)` | `Model.add_ode(...)` with the same arguments; `add_continuous` still works and emits a `DeprecationWarning` |

Code written for the NumPy results of release 0.2.0 changes as follows:

| 0.2.0 | Now |
| --- | --- |
| `records["termination"] == crn.rules.TERMINATIONS.index("stopped")` | `records["termination"] == "stopped"`; the old form is silently all `False` |
| `records["id"]` | `records.index` |
| `records["rule_events"]` | `result.rule_events` |
| `records[name]` (an array) | `records[name].to_numpy()`; integer termination codes are `records["termination"].cat.codes` |
| `statistics["mean"][:, names.index("N")]` | `statistics["mean"]["N"]` |
| `statistics["mean"]` (an array) | `statistics["mean"].to_numpy()` |
| `metadata["paths"]` rows | `result.paths_frame()`, or `result.paths_frame().to_numpy()` |
| `tracker_columns(i)["id"]` and fields | `tracker_columns(i).index` and its columns |
| `event_columns(id)["observables"]` | `event_columns(id).to_numpy()` |
| `keep_paths=True` (one cell) | `keep_paths=1`; `True` now keeps every cell |

Result metadata:

| Key | Content |
| --- | --- |
| `sampler`, `termination`, `diagnostics` | what ran, how the ensemble ended, warnings and errors |
| `events`, `completed_replicates`, `stopped_replicates` | totals |
| `rules`, `observables`, `types`, `report_times` | names and grid, in result order; rule names also as `run.rule_names` and `result.rule_names` |
| `statistics` | `mean`, `sample_variance`, `mean_standard_error`, `samples` and, when asked, `quantiles` as a list of `{p, values}` (also `result.statistics`, where `quantiles` is a dictionary) |
| `paths`, `final_instances`, `observations` | kept series, kept instances, and (with trackers or a log) trackers and logs |
| `model_sha256`, `run_sha256`, `root_seed` | provenance |
| `method`, `continuous` | rules-result/0.3 (`hybrid`, `static_direct`): the request's method, with the tolerances of `hybrid`, and the names of the continuous columns |

Terminations, the categories of `records["termination"]` in the order of
`crn.rules.TERMINATIONS`: `completed`, `stopped`, `event_limit`,
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
- **`'Model' object does not support item assignment`.** A `crn.Model` is
  not a dictionary: set a stop with `model.set_stop(...)`, or edit
  `model.to_dict()` and pass that dictionary (section 12).
- **`crn.compile_loop` fails.** `CompiledLoopError.reason` says why:
  `"cargo_missing"` (install Rust with rustup), `"build_failed"` (read
  `error.log`; for a version that is not on crates.io pass `source=` a
  checkout), or `"mismatch"`. A `crn.CrnError` with `CRN-R020` at `/method`
  means the generator does not support the model (section 8).
- **Different results with another method.** Each method uses the random
  stream differently. Compare distributions, not paths.
- **`CRN-R031`.** A rate is below 2⁻⁵³ of the total, for example a tiny
  mutation probability in a huge population. Static methods select in one
  level; the tree method resolves rules and instances separately.
- **No cell stopped, though some should have.** Compare `termination` with
  a name, `records["termination"] == "stopped"`. A comparison with the
  integer `crn.rules.TERMINATIONS.index("stopped")`, as in release 0.2.0, is
  `False` for every cell.
- **Column positions.** Select DataFrame columns by observable name. In the
  raw matrices of `result.metadata`, observables are sorted by name; use
  `metadata["observables"].index(name)`.
