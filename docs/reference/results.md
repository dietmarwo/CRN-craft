# Results and claims

A simulation request (`crn-craft.run/1.0`, or an older `rules-run` version)
produces one result document. `crn-craft.result/1.0` has the shape of
`rules-result/0.3` plus `statistics.quantiles`. A request with method `fsp`
produces `crn-craft.fsp-result/1.0` instead (below). `crn-craft check`
produces `crn-craft.check/1.0`, a preflight report: `valid: true` means no
input diagnostic of severity error, and nothing was simulated.

The schemas are `crn-craft describe result-1.0`, `fsp-result-1.0` and
`check-1.0`. The Python object `Result` (`crn_craft.rules.Result`) wraps the document: see the
[Python API](python-api.md#result).

## What a result holds

| Field | Content |
| --- | --- |
| `schema`, `sampler` | the contract, and the sampler identity of the method that ran. A replay compares equal only under the same identity. |
| `model_sha256`, `run_sha256`, `root_seed` | provenance: the hashes of the exact bytes consumed |
| `types`, `rules`, `observables`, `report_times` | names and grid; observables are sorted by name, continuous variables follow |
| `requested_replicates`, `completed_replicates`, `stopped_replicates`, `events` | totals |
| `termination` | `completed` when every replicate completed or stopped; otherwise the most severe replicate termination |
| `replicates` | one record per replicate: whether it started, termination, events, stop time, events per rule, random words, peak and created instances, optional diagnostic. A replicate that a cancellation or the wall-time limit kept from starting has `started: false` and zero counters. A cancellation or wall-time limit that comes before the run can still start one replicate per worker thread: such a replicate is recorded as started, ends at its first poll (`cancelled` or `wall_time_limit`) and has 0 events and 0 random words |
| `paths`, `final_instances`, `observations` | kept series, kept final instances, and tracker and log results |
| `statistics` | `mean`, `sample_variance`, `mean_standard_error`, `samples` and `quantiles` (below). Withheld unless every replicate completed or stopped. |
| `method`, `continuous` | the request's method with its tolerances, and the names of the continuous columns |
| `diagnostics` | ensemble-level diagnostics (`CRN-R050`, `CRN-R051`, `CRN-R045`) |

The terminations, in the order of `crn_craft.rules.TERMINATIONS`, are `completed`,
`stopped`, `event_limit`, `sampling_limit`, `instance_limit`,
`wall_time_limit`, `cancelled`, `failed` and `output_limit`.

## Quantiles

A request may list up to 16 probabilities in `outputs.quantiles`, strictly
increasing, each in (0, 1). The result's `statistics.quantiles` then holds,
for each probability in request order, `{"p": ..., "values": ...}`: the
quantile over all replicates for every report time and observable, indexed
like `mean`. The definition is NumPy's default (linear interpolation,
Hyndman and Fan type 7), and the values equal `numpy.quantile` bit for bit.
The median is the quantile at 0.5. Quantiles are withheld with the other
statistics. They keep one value per replicate, report time and observable in
memory during the run, which counts against the 2 GiB workspace limit.

The Python helpers write the request field. `crn_craft.ensemble` asks for
0.025, 0.5 and 0.975 unless `quantiles=` says otherwise, and the result's
`statistics["median"]` and `statistics["quantiles"]` hold DataFrames. Earlier
request versions refuse `outputs.quantiles` (`CRN-R020` at `/schema`).

## Claims

- A result states what the engine computed, with Monte Carlo error. An exact
  method samples the model's process; it does not remove sampling error.
  `tau`, `binomial` and `hybrid` are approximate at their declared
  parameters, and their sampler identities say so.
- Statistics describe all requested replicates or none. A truncated or
  failed ensemble withholds them (`CRN-R050`), because state-dependent
  truncation biases the completed subset. Report censored paths as such.
- Timings, thread counts and file paths are not part of numerical
  comparisons. For a fixed request and sampler identity the numbers are
  byte-identical for every thread count.
- A successful numerical run is not biological validation.

## Finite state projection

`fsp-result/1.0` holds the probability distribution of the model's species
over the states the projection enumerated: `species` (the flat species
labels), `states`, the `snapshots` per report time, the work and sparse edge
counts, the uniformization rate and `numerical_error_certified`, which is
`false`: the tolerance bounds the Poisson-series truncation, not
floating-point error. A snapshot has its `time`, the `probabilities` of the
states, the error terms `mass_defect`, `physical_leakage_estimate`,
`series_tail_upper_bound`, `discarded_weight_estimate` and
`roundoff_mass_residual`
([chapter 7](../guide/07-methods.md#exact-distributions-with-finite-state-projection)),
and `projected_first_moment` and `projected_second_raw_moment`, the moments
of each species over the projected states only. See `Model.fsp` in the
[Python API](python-api.md#modelfsp).

## Saved reports

`Result.save("name.md")` writes a markdown report and `Result.save("name.json")`
writes `crn-craft.report/1.0`, a JSON document with the model, the request and
the results. `crn_craft.load_report(path)` reads the JSON form back. The
guide teaches both ([chapter 4](../guide/04-results.md#save-a-run-in-one-file)),
and the [Python API](python-api.md#resultsave) lists the arguments.

## The metadata of the Python result

`result.metadata` is the result document. Its keys:

| Key | Content |
| --- | --- |
| `sampler`, `termination`, `diagnostics` | what ran, how the ensemble ended, warnings and errors |
| `events`, `completed_replicates`, `stopped_replicates` | totals |
| `rules`, `observables`, `types`, `report_times` | names and grid, in result order; rule names also as `run.rule_names` and `result.rule_names` |
| `statistics` | `mean`, `sample_variance`, `mean_standard_error`, `samples` and, when asked, `quantiles` as a list of `{p, values}` (in `result.statistics`, `quantiles` is a dictionary from probability to DataFrame) |
| `paths`, `final_instances`, `observations` | kept series, kept instances, and (with trackers or a log) trackers and logs |
| `model_sha256`, `run_sha256`, `root_seed` | provenance |
| `method`, `continuous` | the request's method, with the tolerances of `hybrid`, and the names of the continuous columns (`rules-result/0.3` and later) |

In the raw matrices of `result.metadata`, observables are sorted by name; use
`metadata["observables"].index(name)`. Select the columns of a DataFrame by
observable name instead. The attributes and methods of the Python object are
in the [Python API](python-api.md#results).

## Accumulation and report values

Report values are right-continuous: an event at a report time is included.
`crn-craft.rules-result/0.1` holds the type, rule and observable names, one
record per replicate (termination, events, stop time, events per rule, random
words, peak instances (including the one that exceeded `max_instances`) and
created instances, diagnostic), the
observable series of the kept paths, final instances of the kept replicates
(sorted by type and attributes) and, when every replicate completed or
stopped, the mean, sample variance and standard error per report time and
observable. Statistics are accumulated in replicate order (Welford) while the
ensemble runs; only kept series stay in memory. A path's random words include
product draws, so the event that hits the sampling budget can take
`rng_words` slightly past it. Terminations: `completed`, `stopped`, `event_limit`,
`sampling_limit`, `instance_limit`, `wall_time_limit`, `cancelled`, `failed`,
`output_limit`; the ensemble reports `completed` when every replicate completed
or stopped.

## Event logs and tracker columns

`result.event_columns(replicate)` returns the bounded log of one replicate as
an xarray DataArray with dimensions (`entry`, `variable`): the observables
after each event, and the coordinates `time`, `rule` (the index, -1 for a
leap) and `rule_name`. Without an argument it returns every log with
dimensions (`replicate`, `entry`, `variable`), padded with NaN (and `rule` -2)
past each replicate's `length`. `result.event_frame()` returns the DataFrames
of earlier releases: one row per logged event, indexed by `entry` or by
(`replicate`, `entry`). `result.tracker_columns(index)` returns the scalar
tracker summaries as a DataFrame indexed by replicate id, with NaN for
missing times or peaks. Ragged participant identities, leap multiplicities
and episodes remain in `result.observations`.

Semantics of observation. Trackers capture zero crossings, last crossing
times, extrema, time integrals and first passage. A tracked interval runs
until the replicate's stop or horizon. A predicate that evaluates to NaN does
not count as a hit. Episodes keep the configured prefix and summarize omitted
duration and peak. An event log retains up to `max_entries` selected events
per path and counts omissions; its entries include the rule index,
participant type and attributes, and the post-event observables. Retention
counts enter the admission checks, so use small cohorts for long logs.

A live callback receives `fn(replicate_id, event)` on native workers.
Returning `False` stops that replicate. Each path preserves event order.
Callbacks acquire the interpreter lock, and the interleaving across
replicates depends on the scheduler. An exception cancels the run and
propagates. A `batch_callback` receives chunks of the retained log after the
simulation. The Rust `execute_with_observer` entry point exposes the same
committed events and supports a per-path stop; the default observer is
compiled away. Ctrl-C and cancel files use cooperative cancellation, and an
interrupted ensemble withholds its statistics.
