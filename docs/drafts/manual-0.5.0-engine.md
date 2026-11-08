# Draft manual text for 0.5.0: engine stream

User-facing text for `docs/mtdna-user-manual.md`, written by the engine
stream (I4 and I6) for the coordinator to merge in I7. Each Python block
runs from the repository root and prints the `text` block after it
(`check_manual.py` normalizes timings). The blocks use only the builder API
of 0.4.0 and JSON documents, so they do not depend on the Python stream.

## 1. Compiled loops with string rates (I4)

Replaces, in "Optional compiled event loop", the example's last part
("Without cargo: does the generator support a model?"), whose expected
output `CRN-R020 generated rules loops support constant or vacancy rates
and sum observables` no longer occurs, and the second bullet of
"Limitations". Insert the new example after the three-class example.

### Text: after the three-class example

`compile_loop` also compiles string rates. A string rate in `rules.Model`
is a whole-rule propensity, and the generator emits it, like any rate
expression, with the interpreter's operations in the interpreter's order:
the same `libm` functions, NaN propagation through comparisons, logic and
`if`, short-circuit `&&` and `||`, and `min` and `max` that keep NaN. The
results stay byte-identical to interpreted Static Direct, which `verify`
checks:

```python
import time

from crn_craft import CrnError, rules

sir = rules.Model("sir")
sir.add_species("S", 9_990).add_species("I", 10).add_species("R")
sir.add_parameter("beta", 0.3).add_parameter("gamma", 0.1).add_parameter("N", 10_000.0)
sir.add_reaction("beta * S * I / N", ["S", "I"], ["I", "I"], name="infection")
sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
loop = rules.compile_loop(sir)  # string rates compile like constant ones
run = rules.ensemble(
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

# What the generator still refuses: models that only another method runs.
atp = rules.load("models/rules/mtdna-three-class-atp-continuous.rules.json")
try:
    rules.loop_source(atp)
except CrnError as error:
    print(error.diagnostics[0]["code"], error.diagnostics[0]["message"])
```

```text
interpreted   0.85 s    9,400,096 events
compiled      0.37 s    9,400,096 events
byte-identical
CRN-R020 continuous variables require the hybrid method
```

On the same SIR model with 2,000 replicates (37.6 million events), timed
from Python including the subprocess, the compiled loop took 1.64 s
against 3.35 s with inverse waits (2.00×) and 1.36 s against 2.97 s with
Ziggurat waits (2.18×), medians of three on the two-CPU host
(`validation/rules/compiled-loop-string-rates.json`).

### Text: replaces the second bullet of "Limitations"

- The generator emits every model that `method="static_direct"` runs:
  constant, vacancy and expression rates, string rates (`propensity`
  rules) and `sum`, `distinct`, `max` and `min` observables. It refuses,
  with `CRN-R020` at `/method` and before cargo runs, models that only
  another method runs (continuous variables, ATP blocks, instances that
  change identity), and loops that would exceed 1 MiB of source, which
  densely coupled networks with many channels reach.
  `rules.loop_source(model)` returns the source, or that error, without
  cargo.

### Also to update

- The docstring of `rules.loop_source` in `python/crn_craft/_rules_loop.py`
  still lists expression rates, `propensity` rules and max or distinct
  observables among the refusals. The engine stream left that file to the
  Python stream; the text above is the replacement.
