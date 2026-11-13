# Methods and speed

Which method should run the mtDNA model, and how fast can it run? This page
runs the same cells under every method that takes the builder model, and
then compiles the model. The setup block is the model of the
[index page](12-mtdna.md#the-model). The blocks of this page continue from
it:

```python
import math

import numpy as np
import crn_craft as crn

END = 25_920.0  # three years (1,080 days), in hours


def three_class(K=1000, w0=1000, m0=0, mu=3e-5, sm=1.5, half_life_days=10, takeover=None):
    """The three-class mtDNA model: wild type W, neutral M, advantaged A."""
    model = crn.Model("mtdna_three_class")
    model.add_parameter({"K": K, "delta": math.log(2) / (24 * half_life_days),
                         "sw": 1.0, "sm": sm, "mu": mu, "fracAdv": 0.2})
    model.add_species({"W": w0, "M": m0, "A": 0})
    free = "max(0, K - W - M - A) / (W + sm * (M + A))"
    model.add_reaction(f"{free} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt")
    model.add_reaction(f"{free} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt")
    model.add_reaction(f"{free} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv")
    model.add_reaction(f"{free} * sw * M", ["M"], ["M", "M"], "replication_mt")
    model.add_reaction(f"{free} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
    for label, species in (("wt", "W"), ("mt", "M"), ("mtAdv", "A")):
        model.add_reaction(f"delta * {species}", [species], [], f"degradation_{label}")
    if takeover is not None:
        model.set_stop(f"W + M + A > 0 && M + A >= {takeover} * (W + M + A)")
    return model
```

## The same cells under every method

Applies: [chapter 7, Which method for which model](07-methods.md#which-method-for-which-model).

Sixty-four cells with the takeover stop, under seven methods. The tree
method runs when `method` is `None`:

```python
import time

model = three_class(takeover=0.9)
methods = {
    "tree": None,
    "static_direct": "static_direct",
    "static_direct, inverse": {"kind": "static_direct", "waiting_time": "inverse",
                               "sparse_trackers": False},
    "next_reaction": "next_reaction",
    "tau (epsilon 0.03)": {"kind": "tau", "epsilon": 0.03},
    "binomial (1 h steps)": {"kind": "binomial", "step": 1.0},
    "hybrid": "hybrid",
}
for label, method in methods.items():
    run = crn.ensemble(model, end=END, report_times=[0.0, 12_960.0, END], replicates=64,
                       seed=5, threads=2, method=method)
    start = time.perf_counter()
    result = run.simulate()
    seconds = time.perf_counter() - start
    stopped = (result.records["termination"] == "stopped").sum()
    final = result.statistics["mean"].iloc[-1]  # the means at the last report time
    fraction = (final["M"] + final["A"]) / final.sum()
    print(f"{label:22s} {seconds:6.2f} s  {result.metadata['events']:>10,} events  "
          f"taken over {stopped:2d}  mutant fraction {fraction:.3f}  "
          f"{result.metadata['sampler']}")
```

```text
tree                     3.59 s   9,234,548 events  taken over  7  mutant fraction 0.130  chacha8-u53/rules-direct/2
static_direct            1.69 s   9,231,158 events  taken over  7  mutant fraction 0.113  chacha8-u53-ziggurat/rules-static-direct/1
static_direct, inverse   1.92 s   9,120,187 events  taken over  6  mutant fraction 0.092  chacha8-u53/rules-static-direct/1
next_reaction            1.98 s   9,329,541 events  taken over  4  mutant fraction 0.059  chacha8-u53/rules-next-reaction/1
tau (epsilon 0.03)       9.72 s   9,080,008 events  taken over  8  mutant fraction 0.155  chacha8-u53/rules-approximate-cgp-poisson/2
binomial (1 h steps)     0.93 s   9,158,734 events  taken over  6  mutant fraction 0.088  chacha8-u53/rules-approximate-binomial-split/2
hybrid                   3.74 s   9,123,754 events  taken over  6  mutant fraction 0.092  chacha8-u53/rules-approximate-hybrid-dopri54/1
```

The takeover counts differ because each method uses the random numbers in
its own way ([chapter 7](07-methods.md#same-law-different-random-numbers)).
With 64 cells, 4 to 8 takeovers are within sampling error of each other. The approximate methods, `tau` and `binomial`, are labeled as
approximate in their sampler names. `hybrid` integrates between events even
without a continuous variable, so it costs more per event. At `K = 1000` the
free places below `K` are few, so `tau` finds no safe leap and runs exact
events, at a higher cost per event; leaping pays off when every rate changes
slowly. Static Direct is the method to use. The fast path `mtdna_direct`
needs the model file ([page 14a](14a-mtdna-model-files.md#the-fast-path)).

## Matching an hourly scheme

Applies: [chapter 7, Approximate leaping](07-methods.md#approximate-leaping), and [chapter 8, Reproducing a published simulation](08-reproducibility.md#reproducing-a-published-simulation).

Many published mtDNA simulations advance the cell in fixed steps of one
hour. To compare with such a result, run the model exactly first, under
Static Direct, and compare the distributions, not single cells. The method
`{"kind": "binomial", "step": 1.0}` then imitates the hourly scheme: in
each step it removes degraded molecules, freezes the replication rates, and
splits every parent's outcomes by their competing-risk probabilities. Its
steps lie on a fixed lattice, so the report times do not change a path.
Frozen rates can put the population a few molecules above `K`. A difference
between the binomial and the exact results above is the error of the step,
and a difference between the binomial result and the published one is a
difference of the schemes, such as how simultaneous new mutants are
grouped.

## A compiled model

Applies: [chapter 5, Compiled models](05-speed.md#compiled-models).

`crn.compile_loop(model)` compiles the event loop of the builder model. The
results stay byte-identical, and the run takes a fraction of the time. The
first build takes about a minute and needs a Rust toolchain; later calls
read the cache:

```python
import time

model = three_class(takeover=0.9)
loop = crn.compile_loop(model)  # about a minute the first time, then cached
run = crn.ensemble(model, end=END, report_times=[0.0, 12_960.0, END], replicates=200,
                   seed=5, threads=2, method="static_direct",
                   trackers=[{"observable": "M"}, {"observable": "A"}])
for label, compiled in (("interpreted", None), ("compiled", loop)):
    start = time.perf_counter()
    result = run.simulate(compiled=compiled)
    seconds = time.perf_counter() - start
    stopped = (result.records["termination"] == "stopped").sum()
    print(f"{label:12s} {seconds:5.2f} s  {result.metadata['events']:>11,} events  "
          f"taken over {stopped:2d}  {result.metadata['sampler']}")
loop.verify(run)  # runs the request both ways; raises unless byte-identical
print("byte-identical")
```

```text
interpreted   5.46 s   28,452,894 events  taken over 26  chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1
compiled      0.76 s   28,452,894 events  taken over 26  chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1
byte-identical
```

For this model the compiled loop pays off more than for most: the builder
writes every rate as an expression, which the interpreter evaluates at
every event, and the compiled loop computes it natively. In the measurement
of 13 November 2026 on a two-CPU host, the compiled builder model ran 7.1
times as fast as the interpreted one, and took 1.4 times as long as the
model file under the fast path
([measurement](https://github.com/dietmarwo/CRN-craft/blob/main/validation/merge/RESULTS.md)).
A loop is compiled for one parameter set, so it does not speed up a sweep
or a fit ([chapter 5, Limits](05-speed.md#limits)).

See also: [Choosing a method](07-methods.md),
[The fast path](14a-mtdna-model-files.md#the-fast-path),
[mtDNA performance analysis](../evidence/mtdna-performance-analysis.md).
