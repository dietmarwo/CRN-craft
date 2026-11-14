# 6. Stop conditions, trackers and event logs

Statistics at report times do not answer every question. When did the
epidemic end? How high was its peak? When was a threshold first reached?
These need what happens between report times. This chapter shows four tools
that watch every event. The setup block builds the SIR model of
[chapter 2](02-first-model.md) with both rates as strings. The blocks of
this chapter continue from it:

```python
from collections import Counter

import crn_craft as crn

sir = crn.Model("sir")
sir.add_species({"S": 490, "I": 10, "R": 0})
sir.add_parameter({"beta": 5e-4, "gamma": 0.1})
sir.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
```

## Stop conditions

`set_stop` gives the model a condition over parameters and species counts.
It is checked at time 0 and after every event. The first time it is true,
the replicate ends: its termination is `stopped`, and its `stop_time` is the
time of that event. Here an epidemic is over when nobody is infected:

```python
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
over = result.records["termination"] == "stopped"
print("over by day 100:", (over & (result.records["stop_time"] <= 100.0)).mean())
print(result.records.loc[over, "stop_time"].quantile([0.05, 0.5, 0.95]).round(1).to_dict())
```

```text
{'stopped': 1000, 'completed': 0, 'event_limit': 0, 'sampling_limit': 0, 'instance_limit': 0, 'wall_time_limit': 0, 'cancelled': 0, 'failed': 0, 'output_limit': 0}
mean end of the epidemic: day 110.5
           I      R      S
time                      
0.0     10.0    0.0  490.0
1000.0   0.0  446.1   53.9
over by day 100: 0.3
{0.05: 85.8, 0.5: 107.2, 0.95: 143.2}
```

A stopped replicate keeps its state to the end of the report times, so the
means at day 1000 are the final sizes of the epidemics. The last two lines
are the two usual outputs of a stop: the share of replicates stopped by a
given time, and the distribution of the stop times.

Any expression can be a stop, for example `"I >= 100"` or
`"R > 0.5 * (S + I + R)"`. A stopped replicate counts as complete, so the
statistics are there. `set_stop(None)` removes the condition. To record a
threshold without ending the replicate, use a tracker.

## Trackers

A tracker follows one species at every event, at native speed, and reports
per replicate:

- `appearances` and `extinctions`: how often the count rose from 0 or fell
  to 0, and `last_appearance` and `last_extinction`: when it last did;
- `minimum` and `maximum`;
- `integral`: the time integral of the count up to the stop or the end;
- `first_passage`: the first time an optional condition was true;
- with `keep_episodes`, the first excursions above zero, each with start,
  end, peak and duration. Later excursions are summarized: their number in
  `omitted_episodes`, their total duration in `omitted_duration` and their
  highest peak in `omitted_peak`.

```python
run = crn.ensemble(
    sir, end=1000.0, report_times=[0.0, 1000.0], replicates=1000, seed=1,
    threads=2, method="static_direct",
    trackers=[{"observable": "I", "first_passage": "I >= 100"}],
)
result = run.simulate()
tracker = result.tracker_columns(0)  # tracker 0: one row per replicate
print(tracker[["maximum", "integral", "first_passage", "last_extinction"]].head(3).round(1))
print("mean peak of I:", tracker["maximum"].mean())
print("mean infected days per epidemic:", tracker["integral"].mean().round(1))
reached = tracker["first_passage"].dropna()  # NaN: never reached 100
print(f"reached 100 infected: {len(reached) / 1000:.3f}, median day {reached.median():.1f}")
```

```text
           maximum  integral  first_passage  last_extinction
replicate                                                   
0            115.0    4336.9           29.3            119.5
1            133.0    5046.4           14.7            107.4
2            122.0    4529.4           14.9             89.4
mean peak of I: 128.23
mean infected days per epidemic: 4457.5
reached 100 infected: 0.948, median day 20.7
```

The peak is exact, where the mean at report times would miss it. The
integral of `I` is the total time spent infected, so the integral over the
observed time is the time-averaged count. `first_passage` is the exact time
of the event that made the condition true. `last_extinction` equals the stop
time of the section before: the same seed gives the same replicates, and a
tracker does not change them.

## Event logs

An event log keeps, for every replicate, up to `max_entries` events of the
listed reactions, or of all reactions if the list is empty. Each entry holds
the time, the reaction and the species counts after the event. This log
keeps the first five infections of every replicate:

```python
run = crn.ensemble(
    sir, end=1000.0, report_times=[0.0, 1000.0], replicates=200, seed=1,
    threads=2, method="static_direct",
    event_log={"max_entries": 5, "rules": ["infection"]},
)
result = run.simulate()
print(result.event_frame(0))  # the log of replicate 0 as a data frame
log = result.event_columns()  # every log as one array
print(dict(log.sizes))
fifth = log["time"].sel(entry=4).to_series()
print(f"fifth infection: median day {fifth.median():.2f}")
print("later infections, not logged:", sum(o["omitted_log_entries"] for o in result.observations))
```

```text
           time  rule  rule_name     I    R      S
entry                                             
0      1.859423     0  infection   9.0  2.0  489.0
1      2.095936     0  infection   9.0  3.0  488.0
2      2.139428     0  infection  10.0  3.0  487.0
3      2.653214     0  infection  11.0  3.0  486.0
4      3.866311     0  infection   9.0  6.0  485.0
{'replicate': 200, 'entry': 5, 'variable': 3}
fifth infection: median day 1.88
later infections, not logged: 86268
```

`result.event_frame()` without an argument gives all logs as one frame,
indexed by replicate and entry. `result.event_columns()` gives them as an
xarray array with the dimensions `replicate`, `entry` and `variable` and the
coordinates `time`, `rule` and `rule_name`. A replicate with fewer events is
padded with `NaN`.

Logs reserve their space before the run: entries × replicates count against
`max_output_bytes`. A request that does not fit is refused with `CRN-R020`.
The message itemizes the reservation, and its hint names the smallest
`max_output_bytes` that admits the request. Use few entries, few replicates
or a larger limit.

## Callbacks

A callback is a Python function that the engine calls for every event of a
replicate, in order. Returning `False` stops that replicate. This one stops
each epidemic at its 50th infection:

```python
run = crn.ensemble(sir, end=1000.0, report_times=[0.0, 1000.0], replicates=20,
                   seed=1, method="static_direct")
infections = Counter()


def watch(replicate, event):
    """Count the infections of each replicate; stop it at the 50th."""
    if run.rule_names[event["rule"]] == "infection":
        infections[replicate] += 1
    return infections[replicate] < 50


result = run.simulate(callback=watch)
stopped = result.records["termination"] == "stopped"
print("stopped at the 50th infection:", stopped.sum(), "of 20")
print("mean day of the 50th infection:", result.records.loc[stopped, "stop_time"].mean().round(2))
```

```text
stopped at the 50th infection: 20 of 20
mean day of the 50th infection: 9.81
```

The event is a dictionary with the `time`, the index of the reaction in
`rule` and the species counts after the event in `observables`, sorted by
name. Three rules apply:

- A callback reads and may stop. It never changes the state or draws random
  numbers, so a run with a callback has the same events as one without.
- It is slow. Every event crosses into Python and holds the interpreter
  lock, which costs about a hundred times the simulation itself and makes
  the threads take turns. Use it on a few replicates.
- An exception in the callback cancels the run and is raised by
  `simulate()`. The order of calls across replicates is not fixed.

## Which one to use

| Question | Tool | Cost | Result |
| --- | --- | --- | --- |
| When does a condition first hold, if the replicate may end there? | stop condition | none | exact event time |
| Peak, time integral, first passage, appearances and extinctions | tracker | small | exact, one row per replicate |
| Which events happened, when, and what was the state after each? | event log | memory per entry, reserved before the run | exact, up to `max_entries` |
| Anything that needs Python code at every event | callback | slow, one thread at a time | exact |
| The state at chosen times | report times and kept paths ([chapter 4](04-results.md#single-trajectories)) | memory per report time | exact at those times; what happens between them is not seen |

A first passage read from saved paths is a bound, not the time itself
([chapter 11](11-analysis-import-export.md#first-passage-from-saved-paths)).

See also: [Results](04-results.md),
[Methods and their options](../reference/methods.md),
[trackers and logs on the mtDNA model](12b-mtdna-events.md).
