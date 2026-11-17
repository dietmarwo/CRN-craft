# 3. Running ensembles

An *ensemble* is a set of independent trajectories of one model, called
*replicates*. This chapter runs ensembles, makes them reproducible, and
shows how a run can end. The setup block builds the immigration–death model
of [chapter 2](02-first-model.md). The blocks of this chapter continue from
it:

```python
import numpy as np
import crn_craft as crn

model = crn.Model("immigration_death")
model.add_parameter({"birth": 4.0, "death": 0.2})
model.add_species("A", 0)
model.add_reaction("birth", [], ["A"], name="birth")
model.add_reaction("death * A", ["A"], [], name="death")
```

## Many trajectories in one call

One call runs any number of independent trajectories on native threads and
returns their statistics. There is no loop to write and no result to collect
by hand. `crn.ensemble` builds a run request, checks it and prepares it.
`simulate()` then runs the replicates:

```python
run = crn.ensemble(
    model, end=30.0, report_times=list(np.linspace(0, 30, 31)),
    replicates=400, seed=42, keep_paths=8, method="static_direct",
)
result = run.simulate()
assert result.metadata["termination"] == "completed"

t = np.array(result.metadata["report_times"])
mean = result.statistics["mean"]["A"].to_numpy()
exact = 20 * (1 - np.exp(-0.2 * t))
print(result.metadata["events"], "events")
print(f"mean at t=30: {mean[-1]:.2f}, exact {exact[-1]:.2f}")
assert abs(mean[-1] - exact[-1]) < 4 * np.sqrt(exact[-1] / 400)
```

```text
88037 events
mean at t=30: 20.09, exact 19.95
```

The immigration–death law is Poisson with mean `20 (1 - exp(-0.2 t))`, so
the mean and the variance agree:

![Eight paths and the ensemble mean; variance equals the mean](../images/tutorial/immigration-death.png)

`method="static_direct"` is the exact method for a reaction network
(Gillespie's direct method). Name it in every call: without `method`, a
slower exact method for models with changing species runs
([chapter 7](07-methods.md)).

## The everyday arguments

| Argument | Meaning |
| --- | --- |
| `end` | the time at which every replicate ends |
| `report_times` | increasing times in [0, `end`] at which the state is recorded |
| `replicates` | the number of trajectories; 1 by default |
| `seed` | the root seed; 0 by default |
| `threads` | 1 to 64 native threads; 1 by default |
| `keep_paths` | how many replicates keep their whole series: a count, or `"all"` ([chapter 4](04-results.md#single-trajectories)) |
| `method` | a method name, or a dictionary with options ([chapter 7](07-methods.md)) |
| `quantiles` | probabilities for the quantile statistics ([chapter 4](04-results.md#quantiles-and-the-median)) |

Report times are right-continuous: an event exactly at a report time is
included. Statistics exist only at the report times, so choose them for the
question. [Python API](../reference/python-api.md) lists every argument.

## Seeds and threads

Every replicate draws from its own random stream, keyed by the root seed and
the replicate's index. A result therefore depends on the seed and the
method, and not on the number of threads:

```python
def run(seed, threads, replicates=64):
    return crn.ensemble(
        model, end=30.0, report_times=[0.0, 10.0, 20.0, 30.0],
        replicates=replicates, seed=seed, threads=threads, keep_paths="all",
        method="static_direct",
    ).simulate()


one, two, other = run(7, 1), run(7, 2), run(8, 2)
same = one.records.equals(two.records) and one.rule_events.equals(two.rule_events)
print("1 and 2 threads, same seed, identical records:", same)
print("identical statistics:", one.statistics["mean"].equals(two.statistics["mean"]))
print("identical paths:", one.paths.equals(two.paths))
print("seeds 7 and 8, identical events:", one.records["events"].equals(other.records["events"]))
# Replicate i always uses the same stream, whatever the size of the ensemble.
small = run(7, 1, replicates=8)
print("first 8 of 64 equal an 8-replicate run:",
      small.records["events"].equals(one.records["events"].head(8)))
```

```text
1 and 2 threads, same seed, identical records: True
identical statistics: True
identical paths: True
seeds 7 and 8, identical events: False
first 8 of 64 equal an 8-replicate run: True
```

The same seed gives the same result again, on any machine and any number of
threads. Replicate *i* is the same trajectory in a larger ensemble, so a run
can be extended without changing what it already gave. Another seed gives
other trajectories of the same law. [Chapter 5](05-speed.md#threads) times
the threads.

## Check the termination

A run does not always complete. `result.metadata["termination"]` says how
the ensemble ended, and `result.records["termination"]` says it for every
replicate, with one of nine names:

| Termination | Meaning |
| --- | --- |
| `completed` | the replicate reached `end` |
| `stopped` | a stop condition ended it ([chapter 6](06-stops-trackers-logs.md)); this counts as complete |
| `event_limit`, `sampling_limit`, `instance_limit`, `wall_time_limit`, `output_limit` | a budget was reached (next section) |
| `cancelled` | the run was stopped from outside (last section) |
| `failed` | a rate or a stop could not be evaluated; `result.record_diagnostics` says why |

Always check the termination before using statistics. If any replicate ended
on a budget, was canceled or failed, `result.statistics` is `None`: a
partial ensemble is never averaged as if it were complete.

## Budgets

Every run has explicit limits. A replicate that reaches one stops there.

| Limit | Argument | Default | Termination | Code |
| --- | --- | --- | --- | --- |
| events per replicate | `max_events` | 100,000,000 | `event_limit` | `CRN-R040` |
| random numbers per replicate | `max_sampling_steps` | 16 × `max_events` | `sampling_limit` | `CRN-R041` |
| live instances per replicate, for species with attributes ([chapter 14](14-model-files.md#when-species-are-not-a-fixed-list)) | `max_instances` | 100,000 | `instance_limit` | `CRN-R042` |
| wall time of the whole run | `wall_time_s` | none | `wall_time_limit` | `CRN-R043` |
| size of the result | `max_output_bytes` | 64 MiB | `output_limit` | `CRN-R045` |

The event and sampling limits are a safety net against a model that runs
away, for example a rate with a sign error that makes the population
explode. They count work, not seconds, so a run stops at the same event on
every machine. The default lets a replicate run 10⁸ events, about ten
seconds of work at ten million events per second. `max_sampling_steps` follows `max_events` unless you give
it: an event draws one to about three random numbers, so 16 per event leaves
room, and a long run needs only `max_events` raised.

Here 50 events are too few for a run to time 30:

```python
short = crn.ensemble(model, end=30.0, report_times=[0.0, 30.0], replicates=8,
                     seed=1, max_events=50, method="static_direct").simulate()
print(short.metadata["termination"], short.statistics)
print(short.records["termination"].value_counts()["event_limit"], "replicates hit the limit")
```

```text
event_limit None
8 replicates hit the limit
```

A replicate that ends early is not an error, but `simulate()` warns with
`crn.IncompleteEnsembleWarning`, which Python prints on stderr, every time,
also for each batch of a loop. Here the warning is caught to show its
message:

```python
import warnings

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    crn.ensemble(model, end=30.0, report_times=[0.0, 30.0], replicates=8, seed=1,
                 max_events=50, method="static_direct").simulate()
print(caught[0].category.__name__)
print(caught[0].message)
```

```text
IncompleteEnsembleWarning
8 of 8 replicates ended early: 8 event_limit (raise max_events); result.statistics is None. result.records["termination"] gives the termination of every replicate. Pass strict=True to raise IncompleteEnsembleError instead.
```

The message counts the replicates per termination and names the argument
to raise for each. `strict=True` turns the warning into an error,
`crn.IncompleteEnsembleError`, whose `result` still holds the records:

```python
try:
    crn.ensemble(model, end=30.0, report_times=[0.0, 30.0], replicates=8, seed=1,
                 max_events=50, method="static_direct", strict=True).simulate()
except crn.IncompleteEnsembleError as error:
    print(error)
    print(error.counts, len(error.result.records))
```

```text
8 of 8 replicates ended early: 8 event_limit (raise max_events); result.statistics is None. result.records["termination"] gives the termination of every replicate.
{'event_limit': 8} 8
```

The remedy is a larger limit, not an average over cut-off replicates.
`max_events=None` removes the event and sampling limits: the request then
records the largest possible value, 2⁶⁴ − 1, so it still replays exactly.
Without a limit, a runaway model runs until it exhausts memory, so give such
a run a wall-time limit. It ends a run that takes longer than planned, at a
point that depends on the machine:

```python
slow = crn.ensemble(model, end=1e7, report_times=[0.0, 1e7], replicates=4,
                    seed=1, threads=2, max_events=None, wall_time_s=0.2,
                    method="static_direct").simulate()
print(slow.metadata["termination"], slow.statistics)
print([d["code"] for d in slow.metadata["diagnostics"]])
print(sorted(set(slow.records["termination"])))
```

```text
wall_time_limit None
['CRN-R050']
['wall_time_limit']
```

The warning `CRN-R050` in the result says that the ensemble is incomplete
and its statistics are withheld. The records still account for every
replicate. A study that expects early ends, such as a parameter search that
gives up on extreme candidates, can silence the Python warning with
`warnings.simplefilter("ignore", crn.IncompleteEnsembleWarning)`. A
canceled run ([next section](#stopping-a-run)) and a replicate that a stop
condition ended are not reported.

## Stopping a run

A long run can be stopped from outside by creating a file. The run ends with
the termination `cancelled`:

```python
import pathlib

pathlib.Path("stop").touch()
long_run = crn.ensemble(model, end=1e6, report_times=[0.0, 1e6], replicates=8,
                        max_events=None, method="static_direct")
stopped = long_run.simulate(cancel_file="stop")
print(stopped.metadata["termination"], stopped.statistics)
```

```text
cancelled None
```

Here the file exists before the run starts, so the run ends at once. On a
batch system, another job or a time limit creates the file. The engine looks
for it every 1,024 events. **Ctrl-C** cancels a running call in the same
way and raises `KeyboardInterrupt` within about 50 ms.

See also: [Results](04-results.md),
[Speed](05-speed.md),
[Troubleshooting](troubleshooting.md#runs-that-end-early).
