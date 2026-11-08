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

## 2. Independent ODEs (I6)

A new subsection of "Continuous variables and `hybrid`", after the ATP
example. The block writes the model as JSON; once the Python stream's
`Model.add_ode(name, initial, derivative, independent=True)` is merged, the
coordinator may show the builder instead: it writes the same
`"independent": true` and selects rules-model/0.4.

### Text

Some continuous variables do not depend on the molecules at all: a forcing
such as a daily rhythm, a drug concentration that decays on its own, or a
hazard that grows with age. Mark such a variable `"independent": true`
(rules-model/0.4). Its `derivative` may then read only parameters and other
independent variables; anything else is refused with `CRN-R012` at the
derivative's path. Rates, propensities and the derivatives of the other
variables read it as before, and results report it like any continuous
variable.

`hybrid` integrates the independent variables once per prepared run, over
the whole horizon at the request's `rtol` and `atol`, and every replicate
and thread reads them from that one solution, so results do not depend on
`threads`. Each replicate integrates only the other variables and the
integrated propensity. A large or stiff independent system is no longer
integrated again in every replicate: on an eight-compartment pulse driving
arrivals, `independent: true` was 2.8 times as fast
(`validation/hybrid/RESULTS.md`, section E). When events dominate, as with
a two-variable rhythm on the three-class model
(`examples/mtdna/14_independent_ode.py`), the run time hardly changes, but
the rhythm stays more accurate because it is not restarted at every event.

Bounds keep their meaning. With `on_bound: "clamp"` the variable is held at
its bound. With `"stop"`, reaching the bound ends the replicate, and since
the variable is the same in every replicate, every replicate that is still
running stops at the same time. Without `independent`, a rules-model/0.4
document means what it means as 0.3.

```python
import math

from crn_craft import CrnError, rules

# Gompertz mortality: the hazard m grows as dm/dt = g m whatever the
# molecules do, so m is an independent ODE (rules-model/0.4).
model = {
    "schema": "crn-craft.rules-model/0.4",
    "name": "gompertz",
    "parameters": {"g": 0.5, "m0": 0.1},
    "types": [{"name": "A"}],
    "initial": [{"type": "A", "count": 20}],
    "observables": {"A": {"type": "A"}},
    "rules": [
        {
            "name": "death",
            "reactants": [{"label": "a", "type": "A"}],
            "rate": "m",
            "products": [],
        }
    ],
    "continuous": [
        {"name": "m", "initial": "m0", "derivative": "g * m", "independent": True}
    ],
}
result = rules.ensemble(
    model,
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
model["continuous"][0]["derivative"] = "g * m * A"
try:
    rules.ensemble(model, end=4.0, report_times=[4.0], method="hybrid")
except CrnError as error:
    print(error.diagnostics[0]["code"], error.diagnostics[0]["path"])
    print(error.diagnostics[0]["message"])
```

```text
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

### Also to update

- The examples table of `examples/mtdna/README.md`: add
  `14_independent_ode.py` under section 10, continuous variables and
  `hybrid`.
- The contract table of the manual (if it lists model versions): add
  `crn-craft.rules-model/0.4`, 0.3 plus `independent`.
