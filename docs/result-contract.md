# Results and claims

A simulation request (`crn-craft.run/1.0`, or an older `rules-run` version)
produces one result document. `crn-craft.result/1.0` has the shape of
`rules-result/0.3` plus `statistics.quantiles`. A request with method `fsp`
produces `crn-craft.fsp-result/1.0` instead (below). `crn-craft check`
produces `crn-craft.check/1.0`, a preflight report: `valid: true` means no
input diagnostic of severity error, and nothing was simulated.

The schemas are `crn-craft describe result-1.0`, `fsp-result-1.0` and
`check-1.0`. The Python object `Result` (`crn_craft.rules.Result`) wraps the document: see the
[user manual](mtdna-user-manual.md#18-reference).

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
floating-point error. See `Model.fsp` in the
[user manual](mtdna-user-manual.md#18-reference).

## Saved reports

`Result.save("name.md")` writes a markdown report and `Result.save("name.json")`
writes `crn-craft.report/1.0`, a JSON document with the model, the request and
the results. `crn_craft.load_report(path)` reads the JSON form back. The
user manual describes both ([saving results](mtdna-user-manual.md#saving-results)).
