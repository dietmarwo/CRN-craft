# Limits in CRN-craft, and the lessons of the study

Author: Axel Kowald

Part of [tutorial S2](S2-published-study.md).

```python
import math
import crn_craft as crn

degRate = math.log(2) / 240.0
```

## Limits and when to raise them

CRN-craft puts explicit budgets on every run. A buggy or explosive model
then stops with a clear reason instead of exhausting memory or running
forever. The defaults suit typical models, and since 0.8 they also suit
long horizons like the ones of this study. A replicate that hits a budget
stops, and its reason is recorded in `result.records["termination"]`
(`completed`, `stopped`, `event_limit`, `sampling_limit`, `instance_limit`,
`wall_time_limit`, `cancelled`, `failed`).

Always make sure that every replicate completed before computing
statistics. An early stop leaves truncated paths, and dropping those cells
silently biases the result towards fast-ending, often atypical,
trajectories. `simulate()` warns when a replicate did not complete
(`crn.IncompleteEnsembleWarning`), and with `strict=True` it raises
instead:

```python
model = crn.Model("decay")
model.add_reaction(degRate, ["mtDNA"], [], name="degradation")
model.add_species({"mtDNA": 1000})
run = crn.ensemble(model, end=1e6, report_times=[0.0, 1e6], replicates=4, seed=1,
                   max_events=100, strict=True)
try:
    run.simulate()
except crn.IncompleteEnsembleError as error:
    print(error.counts)
```

```text
{'event_limit': 4}
```

| Limit (keyword) | What it bounds | In this study |
| --- | --- | --- |
| `max_events` | reaction events per replicate; default 10⁸ | a 120-year clone cell has about 6·10⁶ events; the default is enough (0.7.1 needed 5·10⁷) |
| `max_sampling_steps` | random-number draws per replicate; by default 16 × `max_events` | follows `max_events`, nothing to set (0.7.1 needed 5·10⁸) |
| `max_instances` | instances (here clones) alive at the same time in one replicate; default 100,000 | the default is enough: about 40 clones at most (0.7.1 needed 5,000 because of the memory reservation of `keep_instances`) |
| `max_output_bytes` | size of the returned result; default 64 MiB | `keep_paths="all"` with 1,441 monthly report times for 1,000 cells is about 23 MB |
| batch size (script) | memory held by kept trajectories | 1,000 cells per `crn.ensemble` call |
| `wall_time_s` | wall-clock time per run | not needed; use it as the guard when `max_events=None` |
| FSP work limit (`CRN-E070`) | state-space work of `Model.fsp` | the state space was too large, so sampling was used instead |

How to choose the limits:

1. Estimate the event count from the model: total propensity × simulated
   time. The default `max_events` of 10⁸ covers about 15 seconds of
   simulation per replicate at typical speeds.
2. Run a prototype with a few hundred replicates at the longest horizon
   and the highest rate, then inspect `records["termination"]`,
   `records["peak_instances"]` and `records["events"]`.
3. Raise only the limit that was hit, and keep it finite. A finite budget
   still guards against a model error that makes the event rate or the
   number of instances explode. `max_events=None` removes the event limit;
   combine it with `wall_time_s`.
4. If a reservation error mentions bytes (`CRN-R020`), its message lists
   the parts of the reservation and its hint names the largest one; reduce
   that one, for example the replicates per call. Do not raise memory
   limits blindly.

### What changed since 0.7.1

The tutorial was first written for CRN-craft 0.7.1. Version 0.8 changed:

- `max_events` defaults to 10⁸, and `max_sampling_steps` follows it, so
  the long runs need no raised limits.
- `keep_instances` no longer reserves `keep_instances × max_instances`
  instances before the run, so `max_instances` need not be lowered.
- `simulate()` warns about incomplete replicates; `strict=True` makes it an
  error.
- The builder writes species with attributes, `crn.copy`, `crn.new`,
  `crn.mass_action` and `add_observable`, so the clone model needs no JSON.
- `reduce="max_count"` gives the largest clone as an observable.
- Pair rules give the labels `x` and `y` at random, so the exact Moran step
  can be one rule.

## Summary of lessons

- **Lump when you can.** If only counts matter, an exact two-reaction model
  with expression propensities plus `static_direct` and `compile_loop` is
  fastest.
- **Use one seed per grid point** in parameter sweeps that feed a
  statistical model, so the samples stay independent.
- **Use species with attributes and `new_id()`** when individual lineages
  matter, together with `reduce="distinct"` and `reduce="max_count"`
  observables and the tree method. Use `keep_instances` and
  `metadata["final_instances"]` only for quantities that observables cannot
  express.
- **Check every replicate's termination status** (or use `strict=True`),
  and raise budgets deliberately, after a prototype run at the hardest
  settings.
- **Record alternative definitions as separate observables** (here types
  with and without wild type). They are cheap, and they can resolve
  ambiguities in a paper's text.

Back to [the start of the chapter](S2-published-study.md).
