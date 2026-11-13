# Troubleshooting

Symptoms, causes and where the guide explains them.

## Runs that end early

- **`statistics` is `None`.** A replicate ended on a budget, was canceled
  or failed, so the statistics are withheld (`CRN-R050`). Look at
  `result.records["termination"]` ([chapter 3](03-running.md#check-the-termination)).
- **`event_limit` or `sampling_limit` in every replicate.** Raise
  `max_events` and `max_sampling_steps` ([chapter 3](03-running.md#budgets)).
- **`observation reservation exceeds budget`.** An event log reserves space
  for every replicate before the run. The hint of the error names the
  smallest `max_output_bytes` that admits the request; or use fewer
  replicates or entries ([chapter 6](06-stops-trackers-logs.md#event-logs)).

## Methods and models

- **A run is much slower than expected.** Name `method="static_direct"`:
  without `method` the tree method runs ([chapter 7](07-methods.md)).
- **Different results with another method or selector.** Each uses the
  random numbers differently. Compare distributions, not paths
  ([chapter 7](07-methods.md#same-law-different-random-numbers)).
- **`CRN-R020` at `/method` for a model with continuous variables.** Only
  `method="hybrid"` runs them ([chapter 9](09-continuous-variables.md#the-hybrid-method)).
- **`CRN-R031`.** A rate is below 2⁻⁵³ of the total rate, for example a tiny
  mutation probability in a huge population. Rescale the model, or run
  without `method`: the tree method resolves rules and instances separately.
- **`mtdna_direct` is not faster.** Read `metadata["sampler"]` and the
  `CRN-R051` warning, which names the check that failed
  ([page 14a](14a-mtdna-model-files.md#the-fast-path)).

## Python results

- **`'Model' object does not support item assignment`.** A `crn.Model` is
  not a dictionary. Use `model.set_stop(...)`, or edit `model.to_dict()` and
  pass the dictionary ([chapter 14](14-model-files.md#save-and-load-a-model-file)).
- **No replicate stopped, though some should have.** Compare the termination
  with a name: `records["termination"] == "stopped"`.
- **A column holds another species than expected.** Select columns by name.
  Species are sorted by name, also in the raw matrices of `result.metadata`
  ([chapter 4](04-results.md#select-by-name)).
- **A Python callback is slow.** It holds the interpreter lock. Use trackers
  or event logs for whole ensembles
  ([chapter 6](06-stops-trackers-logs.md#which-one-to-use)).

## Compiled models

- **`crn.compile_loop` raises `CompiledLoopError`.** Its `reason` says why:
  `"cargo_missing"` (install Rust with rustup), `"build_failed"` (read
  `error.log`; for a version that is not on crates.io, pass `source=` a
  checkout), or `"mismatch"`.
- **`CRN-R020` at `/method` or `/method/selector`.** The model needs another
  method, or `auto` chose another selector than the linear scan: name
  `selector="linear"` ([chapter 5](05-speed.md#limits)).

## Frequent diagnostic codes

| Code | Meaning |
| --- | --- |
| `CRN-R012` | an expression error or an unknown name |
| `CRN-R020` | an invalid request: hash, budget, method, tracker, log or version |
| `CRN-R030` | an invalid rate, an overflow or an undefined stop in a replicate |
| `CRN-R031` | a positive rate too small to select; see above |
| `CRN-R040` to `CRN-R045` | a budget was reached |
| `CRN-R050` | an incomplete ensemble; statistics withheld |
| `CRN-R051` | a requested fast path did not apply; the message says why |

[Diagnostics](../reference/diagnostics.md) lists every code.
