# 5. Speed: threads, parameter sweeps and compiled models

Three tools make a study fast: threads for one ensemble, parameter binding
for many ensembles of one model, and a compiled event loop for a model that
runs often. None of them changes a result. The setup block builds the SIR
model of [chapter 2](02-first-model.md) with 10,000 individuals. The blocks
of this chapter continue from it:

```python
import time

import crn_craft as crn

sir = crn.Model("sir")
sir.add_species({"S": 9_990, "I": 10, "R": 0})
sir.add_parameter({"beta": 0.3, "gamma": 0.1, "N": 10_000.0})
sir.add_reaction("beta * S * I / N", ["S", "I"], ["I", "I"], name="infection")
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
```

## Threads

`threads=n` divides the replicates among `n` native threads. Every replicate
has its own random stream, so the result is the same for any number of
threads:

```python
def timed(threads):
    run = crn.ensemble(
        sir, end=400.0, report_times=[0.0, 100.0, 200.0, 400.0], replicates=500,
        seed=3, threads=threads, method="static_direct",
    )
    start = time.perf_counter()
    result = run.simulate()
    return result, time.perf_counter() - start


one, seconds_one = timed(1)
two, seconds_two = timed(2)
assert one.statistics["mean"].equals(two.statistics["mean"])
print(f"1 thread  {seconds_one:5.2f} s  {one.metadata['events']:,} events")
print(f"2 threads {seconds_two:5.2f} s  identical statistics")
```

```text
1 thread   0.77 s  9,400,096 events
2 threads  0.39 s  identical statistics
```

The timings are from a host with two CPUs. Use as many threads as the
machine has free CPUs, up to 64. More threads than CPUs give nothing. The
[performance analysis](../evidence/mtdna-performance-analysis.md#3-threads)
has measurements on up to 20 threads.

## What output costs

Simulation time grows with the number of events. Output adds to it:

- Kept paths take memory for species × report times × kept replicates. Keep
  the paths of a few replicates to look at, and all of them only for a
  statistic that needs them.
- Quantiles hold one value per replicate, report time and species while the
  run lasts. `quantiles=()` turns them off.
- Every report time is a row in every statistic. A fine grid for a plot is
  cheap. A grid of 100,000 times for a million replicates is not.

A result larger than `max_output_bytes`, 64 MiB by default, ends as
`output_limit` without statistics and paths: raise the limit or keep less.
The [performance analysis](../evidence/mtdna-performance-analysis.md#8-sizing-a-study)
sizes a large study.

## Parameter sweeps

A sweep runs one model for many parameter values. Prepare the run once, then
bind new values with `bind_parameters`. Here the infection rate `beta`
varies:

```python
template = crn.ensemble(
    sir, end=400.0, report_times=[0.0, 400.0], replicates=200, seed=3,
    threads=2, method="static_direct",
)
for beta in (0.08, 0.15, 0.2, 0.3, 0.5):
    result = template.bind_parameters({"beta": beta}).simulate()
    recovered = result.statistics["mean"]["R"].iloc[-1]
    print(f"beta {beta:4.2f}  mean final R {recovered:7.1f}  "
          f"model {result.metadata['model_sha256'][:8]}")
```

```text
beta 0.08  mean final R    45.5  model 50926733
beta 0.15  mean final R  5752.8  model 35da249c
beta 0.20  mean final R  7964.3  model a1fcd763
beta 0.30  mean final R  9404.3  model 52d43047
beta 0.50  mean final R  9930.9  model 212bd588
```

Below `beta = gamma = 0.1` the epidemic dies out. Above it, most of the
population is infected.

`bind_parameters` takes a mapping from parameter names to numbers and
returns a new prepared run. Only the named parameters change. The request,
the seed, the method, the initial counts and the reactions stay as they are;
to change one of those, prepare a new run. The template and earlier
candidates stay intact. An unknown name or a value that is not a finite
number is refused with `CRN-R020`. Each candidate has its own model hash, so
every result names its parameters.

## Common random numbers

The same seed gives every candidate the same random streams. Differences
between two candidates then come from the parameters and not from different
random numbers, which makes a sweep smoother and helps an optimizer that
compares candidates. For an independent evaluation, such as a final check of
a fitted value, prepare a run with another `seed`.

## Compiled models

`crn.compile_loop(model)` compiles the event loop of one model into a small
native program. The results are byte-identical to the interpreted run, the
run is faster, and you stay in Python:

```python
loop = crn.compile_loop(sir)  # about a minute the first time, then cached
run = crn.ensemble(
    sir, end=400.0, report_times=[0.0, 100.0, 200.0, 400.0], replicates=500,
    seed=3, threads=2, method="static_direct",
)
for label, compiled in (("interpreted", None), ("compiled", loop)):
    start = time.perf_counter()
    result = run.simulate(compiled=compiled)
    seconds = time.perf_counter() - start
    print(f"{label:12s} {seconds:5.2f} s  {result.metadata['events']:>11,} events")
loop.verify(run)  # runs both ways; raises unless the results are byte-identical
print("byte-identical")

# Without a build: loop_source says whether a model can be compiled.
dosed = crn.Model("dosed")
dosed.add_species("A", 10)
dosed.add_ode("drug", 1.0, "-0.1 * drug")
dosed.add_reaction("drug * A", ["A"], [], name="death")
try:
    crn.loop_source(dosed)
except crn.CrnError as error:
    print(error.diagnostics[0]["code"], error.diagnostics[0]["message"])
```

```text
interpreted   0.45 s    9,400,096 events
compiled      0.29 s    9,400,096 events
byte-identical
CRN-R020 continuous variables require the hybrid method
```

A loop can be given in two places. `run.simulate(compiled=loop)` uses it for
one call, as above. `crn.ensemble(..., compiled=loop)` attaches it to the
prepared run, and every `simulate()` then runs compiled; use this when code
that knows nothing of loops calls `simulate()`.

What a compiled model needs:

- A Rust toolchain, 1.85 or newer, from [rustup](https://rustup.rs), and
  access to crates.io for the first build.
- Time for the first build, about a minute. Builds are cached, and a cached
  call returns in a fraction of a second.
- `method="static_direct"`. A model that another method must run, such as
  one with continuous variables, is refused, as the last lines show.

The gain depends on the model. On 12 November 2026, on a two-CPU host, the
compiled mtDNA model of [chapter 12](12-mtdna.md), in its file form, ran 2.1
times as fast as the interpreted one
([measurement](https://github.com/dietmarwo/CRN-craft/blob/main/validation/merge/RESULTS.md)).
A network with many constant-rate reactions can run slower compiled, so time
both. [Methods and their options](../reference/methods.md) has the cache,
the options and every limit.

## Limits

- **Binding prepares each candidate again.** This takes about 0.2 ms for a
  small model. On very short runs it can be a visible share of a sweep.
- **A candidate runs interpreted.** A loop is compiled for one parameter
  set, so `bind_parameters` returns a run without it. Compiling pays for
  many runs with fixed parameters, not for a sweep.
- **Compiled models use the linear scan only.** For a model whose `auto`
  choice is another selector ([chapter 7](07-methods.md)), name
  `selector="linear"` together with `compiled=`.
- **Very large models are refused.** A loop whose source would exceed 1 MiB
  is not generated.
- **No callbacks.** A compiled run takes no Python callback. Trackers and
  event logs work.

See also: [Choosing a method](07-methods.md),
[Fitting parameters to data](13-fitting.md),
[mtDNA performance analysis](../evidence/mtdna-performance-analysis.md).
