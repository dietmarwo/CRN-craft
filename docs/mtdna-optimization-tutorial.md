# Fitting mtDNA models to data: an optimization tutorial

This tutorial estimates two parameters of the three-class mtDNA model:

- the mutation probability `mu`;
- the replication advantage `sm` of advantaged mutants.

They are fitted to single-cell heteroplasmy data with the derivative-free
optimizers of [fcmaes](https://github.com/dietmarwo/fast-cma-es). The
objective is a stochastic simulation, so the tutorial also covers what
that changes: Monte Carlo noise, common random numbers, staged refinement,
held-out checks and how to spend the cores of a machine. Every script is in
`examples/mtdna/` and was run as shown. Timings are from a two-CPU Linux
host.

Read the [user manual](mtdna-user-manual.md) first for the model, the
ensemble API and the result columns. The
[performance analysis](mtdna-performance-analysis.md) helps size the
evaluations.

```sh
pip install crn-craft fcmaes scipy
```

## 1. The problem

Single-cell measurements, such as COX staining or single-cell sequencing,
give the mutant fraction of individual cells at a few ages. The model
predicts, for any `mu` and `sm`, the distribution of mutant fractions across
cells at each age. Fitting means finding the parameters whose predicted
distributions match the observed ones.

Each prediction is an ensemble of simulated cells. A three-year cell at
K = 1000 has about 150,000 events. An objective evaluation with 800 cells
takes about 2.2 s on two cores with `mtdna_direct`, the opt-in fast path
for this model's three-class structure, and a reused prepared template.
Two hundred evaluations are therefore minutes, not hours. Static Direct,
which runs any closed model, takes about 2.8× longer with its default
Ziggurat waits and sparse trackers (user manual, section 8, measured on
12 November 2026): about 20 minutes for the same fit.

To be able to check the answer, the tutorial makes its own data from known
parameters (`mu = 3·10⁻⁵`, `sm = 1.5`), with a seed that the fit never sees.
With real data, replace `OBSERVED` by the measured fractions.

## 2. A first attempt: takeover fractions only

Takeover data are often all that is available: the fraction of cells past a
threshold, for example 90% mutants, at a few ages. The first script fits
`mu` and `sm` to the takeover fractions at 360, 720 and 1080 days, plus the
mean mutant fraction. It runs 100 BiteOpt evaluations.

```python
"""Fit the mutation probability and the replication advantage to takeover data."""

import json
import sys
import time

import numpy as np
from crn_craft import rules
from fcmaes import bitecpp
from scipy.optimize import Bounds

DAYS = np.array([360, 720, 1080])
MODEL = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
MODEL["stop"] = "N > 0 && mutants >= 0.9 * N"


def template(cells, seed, threads=2):
    return rules.ensemble(
        MODEL,
        end=25_920.0,
        report_times=24.0 * np.concatenate(([0], DAYS)),
        replicates=cells,
        seed=seed,
        threads=threads,
        method="mtdna_direct",
    )


def summaries(result):
    """Fraction of cells taken over by each age, and mean mutant fraction."""
    stop = result.records["stop_time"]
    stopped = result.records["termination"] == "stopped"
    taken = np.array([(stopped & (stop <= 24.0 * d)).mean() for d in DAYS])
    mean = result.statistics["mean"].iloc[1:]  # the report times after day 0
    load = (mean["mutants"] / mean["N"]).to_numpy()
    return taken, load


# 1. Synthetic observations from known parameters and an independent seed.
TRUE = {"mu": 3e-5, "sm": 1.5}
if __name__ == "__main__":
    start = time.perf_counter()
    observed = summaries(template(4000, seed=2024).bind_parameters(TRUE).simulate())
    print("observed takeover:", np.round(observed[0], 4))
    print("observed mutant fraction:", np.round(observed[1], 4))
    print(f"data: {time.perf_counter() - start:.1f} s")

SEARCH = template(800, seed=7)  # fixed seed: common random numbers


def objective(x, prepared=SEARCH, observed_data=None):
    mu, sm = 10.0 ** x[0], x[1]
    taken, load = summaries(prepared.bind_parameters({"mu": mu, "sm": sm}).simulate())
    taken_obs, load_obs = observed_data or OBSERVED
    # Binomial standard errors of the observed takeover fractions; a floor
    # keeps empty early bins from dominating.
    se = np.sqrt(np.maximum(taken_obs * (1 - taken_obs), 0.01) / 800)
    return float(
        np.sum(((taken - taken_obs) / se) ** 2)
        + np.sum(((load - load_obs) / (0.1 * load_obs + 1e-3)) ** 2)
    )


if __name__ == "__main__":
    OBSERVED = observed
    evaluations = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    bounds = Bounds([-5.5, 1.0], [-3.5, 2.5])
    start = time.perf_counter()
    found = bitecpp.minimize(
        objective,
        bounds=bounds,
        max_evaluations=evaluations,
        rg=np.random.default_rng(42),
    )
    seconds = time.perf_counter() - start
    mu, sm = 10.0 ** found.x[0], found.x[1]
    print(
        f"BiteOpt: {found.nfev} evaluations in {seconds:.0f} s, "
        f"loss {found.fun:.2f}, mu {mu:.2e}, sm {sm:.3f}"
    )
    # 2. Held out: the fitted and the true parameters with fresh cells.
    check = template(4000, seed=99)
    for label, values in (("fitted", {"mu": mu, "sm": sm}), ("true", TRUE)):
        taken, load = summaries(check.bind_parameters(values).simulate())
        print(
            f"{label:6s} held-out takeover {np.round(taken, 4)} load {np.round(load, 4)}"
        )
```

```text
observed takeover: [0.0125 0.062  0.1078]
observed mutant fraction: [0.0254 0.0675 0.1119]
data: 12.4 s
BiteOpt: 100 evaluations in 219 s, loss 0.90, mu 1.87e-05, sm 2.269
fitted held-out takeover [0.0335 0.083  0.133 ] load [0.0374 0.0823 0.1264]
true   held-out takeover [0.0112 0.0605 0.1078] load [0.0226 0.0682 0.1097]
```

The fit reaches a low loss on its own 800 cells, but with the wrong
parameters: `mu` too small and `sm` too large. On fresh cells it predicts the
data worse than the true parameters. Two things went wrong:

- **Confounding.** Takeover fractions mainly measure how many lineages
  succeed per unit time. Fewer mutations with a larger advantage give
  nearly the same curve as more mutations with a smaller one.
- **Overfitting to noise.** With one fixed set of 800 cells, the optimizer
  can find parameters that fit those cells' accidents. The held-out check
  catches it.

The remedy for the first is more informative data. The remedy for the
second is more cells in a second stage and an independent check.

## 3. A better objective: heteroplasmy distributions

The distribution of mutant fractions separates the two parameters:

- `mu` sets how many cells carry any mutants at all;
- `sm` sets how fast they move towards high mutant fractions.

The objective compares, at each age, the fraction of cells whose mutant
fraction reaches 1%, 10%, 50% and 90%. Each difference is scaled by its
binomial standard error.

```python
"""Fit the mutation probability and the replication advantage to
single-cell heteroplasmy distributions at three ages."""

import json
import sys
import time

import numpy as np
from crn_craft import rules
from fcmaes import bitecpp
from scipy.optimize import Bounds

DAYS = np.array([360, 720, 1080])
EDGES = np.array([0.01, 0.1, 0.5, 0.9])  # mutant fraction thresholds
MODEL = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
MODEL["stop"] = "N > 0 && mutants >= 0.9 * N"
TRUE = {"mu": 3e-5, "sm": 1.5}


def template(cells, seed, threads=2):
    return rules.ensemble(
        MODEL,
        end=25_920.0,
        report_times=24.0 * DAYS,
        replicates=cells,
        seed=seed,
        threads=threads,
        method="mtdna_direct",
        keep_paths="all",
        max_output_bytes=256 << 20,
    )


def heteroplasmy(result):
    """Fraction of cells whose mutant fraction reaches each edge, per age."""
    paths = result.paths  # cells x report times x observables
    fraction = (paths.sel(variable="mutants") / paths.sel(variable="N")).to_numpy()
    return (fraction[:, :, None] >= EDGES).mean(axis=0)  # ages x edges


def loss(simulated, observed, cells):
    # Binomial standard errors of the observed fractions, with a floor.
    se = np.sqrt(np.maximum(observed * (1 - observed), 0.005) / cells)
    return float(np.sum(((simulated - observed) / se) ** 2))


SEARCH_CELLS = 800
SEARCH = template(SEARCH_CELLS, seed=7)  # one seed for every candidate


def objective(x):
    candidate = SEARCH.bind_parameters({"mu": 10.0 ** x[0], "sm": x[1]})
    return loss(heteroplasmy(candidate.simulate()), OBSERVED, SEARCH_CELLS)


if __name__ == "__main__":
    start = time.perf_counter()
    OBSERVED = heteroplasmy(template(4000, seed=2024).bind_parameters(TRUE).simulate())
    print(f"observed ({time.perf_counter() - start:.0f} s), rows are days {DAYS}:")
    print(np.round(OBSERVED, 4))
    evaluations = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    start = time.perf_counter()
    found = bitecpp.minimize(
        objective,
        bounds=Bounds([-5.5, 1.0], [-3.5, 2.5]),
        max_evaluations=evaluations,
        rg=np.random.default_rng(42),
    )
    mu, sm = 10.0 ** found.x[0], found.x[1]
    print(
        f"stage 1: {found.nfev} evaluations in {time.perf_counter() - start:.0f} s, "
        f"loss {found.fun:.1f}, mu {mu:.2e}, sm {sm:.3f}"
    )
    # Stage 2: more cells and a new seed in a box around the first estimate.
    SEARCH_CELLS = 2400
    SEARCH = template(SEARCH_CELLS, seed=8)
    start = time.perf_counter()
    lower = [found.x[0] - 0.2, max(1.0, found.x[1] - 0.3)]
    upper = [found.x[0] + 0.2, found.x[1] + 0.3]
    refined = bitecpp.minimize(
        objective,
        bounds=Bounds(lower, upper),
        x0=found.x,
        max_evaluations=evaluations // 3,
        rg=np.random.default_rng(43),
    )
    first = {"mu": mu, "sm": sm}
    mu, sm = 10.0 ** refined.x[0], refined.x[1]
    print(
        f"stage 2: {refined.nfev} evaluations in {time.perf_counter() - start:.0f} s, "
        f"loss {refined.fun:.1f}, mu {mu:.2e}, sm {sm:.3f}"
    )
    check = template(4000, seed=99)
    for label, values in (("stage 1", first),):
        held = heteroplasmy(check.bind_parameters(values).simulate())
        print(f"{label:7s} held-out loss {loss(held, OBSERVED, 4000):7.1f}")
    for label, values in (("stage 2", {"mu": mu, "sm": sm}), ("true", TRUE)):
        held = heteroplasmy(check.bind_parameters(values).simulate())
        print(f"{label:7s} held-out loss {loss(held, OBSERVED, 4000):7.1f}")
```

```text
observed (11 s), rows are days [ 360  720 1080]:
[[0.0648 0.036  0.025  0.0125]
 [0.122  0.0828 0.0712 0.062 ]
 [0.1798 0.1352 0.1212 0.1078]]
stage 1: 100 evaluations in 209 s, loss 1.5, mu 2.99e-05, sm 1.601
stage 2: 33 evaluations in 208 s, loss 5.7, mu 3.06e-05, sm 1.475
stage 1 held-out loss    68.6
stage 2 held-out loss     9.8
true    held-out loss     6.8
```

How the script is built:

- **One prepared template.** `template(cells, seed)` prepares the ensemble
  once. Every candidate is `SEARCH.bind_parameters({...}).simulate()`, which
  binds new values without rebuilding the request.
- **Common random numbers.** All candidates of a stage share one seed, so
  every candidate sees the same random streams. Differences between
  candidates are then due to the parameters, not to different random cells.
- **A log scale for `mu`.** The optimizer works on `log10(mu)` in
  [−5.5, −3.5], which spans two orders of magnitude evenly.
- **Two stages.** Stage 1 (800 cells, 100 evaluations) finds the region.
  Stage 2 (2400 cells, a new seed, 33 evaluations) refines inside a box around
  it. The box is ±0.2 in `log10(mu)` and ±0.3 in `sm`.
- **A held-out check.** 4000 fresh cells with a third seed score stage 1,
  stage 2 and the true parameters. Only this number says how good a fit is.

Stage 2 recovers `mu = 3.06·10⁻⁵` and `sm = 1.475`, against the true 3·10⁻⁵
and 1.5. Its held-out loss (9.8) is close to that of the true parameters
(6.8), while stage 1's estimate scores 68.6.

## 4. Noise and the shape of the objective

```python
"""Monte Carlo noise of the objective, and parallel evaluation with fcmaes."""

import time

import numpy as np
from fcmaes import de
from loguru import logger
from scipy.optimize import Bounds

import fit_heteroplasmy as fit

fit.OBSERVED = fit.heteroplasmy(
    fit.template(4000, seed=2024).bind_parameters(fit.TRUE).simulate()
)
x_true = np.array([np.log10(3e-5), 1.5])

# 1. The same parameters under five search seeds: the loss is noisy.
values = []
for seed in range(5):
    fit.SEARCH = fit.template(800, seed=seed)
    values.append(fit.objective(x_true))
print("loss at the true parameters, five seeds:", np.round(values, 1))

# 2. With one fixed seed the loss is a deterministic function of the
#    parameters: a scan over sm is repeatable, but not smooth.
fit.SEARCH = fit.template(800, seed=7)
for sm in np.linspace(1.3, 1.7, 9):
    print(f"sm {sm:.2f}: loss {fit.objective([x_true[0], sm]):6.1f}")


# 3. Parallel candidates: each process prepares its own template once.
_LOCAL = {}


def parallel_objective(x):
    if "search" not in _LOCAL:
        _LOCAL["search"] = fit.template(800, seed=7, threads=1)
    fit.SEARCH = _LOCAL["search"]
    return fit.objective(x)


if __name__ == "__main__":
    logger.remove()  # fcmaes logs every improvement through loguru
    start = time.perf_counter()
    found = de.minimize(
        parallel_objective,
        bounds=Bounds([-5.5, 1.0], [-3.5, 2.5]),
        popsize=8,
        max_evaluations=32,
        workers=2,
        rg=np.random.default_rng(5),
    )
    print(
        f"DE, 2 workers x 1 thread: {found.nfev} evaluations in "
        f"{time.perf_counter() - start:.0f} s, best loss {found.fun:.1f}"
    )
```

```text
loss at the true parameters, five seeds: [24.1  3.5  6.2 10.6 15.5]
sm 1.30: loss  166.5
sm 1.35: loss  101.1
sm 1.40: loss   68.1
sm 1.45: loss   32.3
sm 1.50: loss   20.9
sm 1.55: loss    4.9
sm 1.60: loss    1.7
sm 1.65: loss   10.5
sm 1.70: loss   16.3
DE, 2 workers x 1 thread: 32 evaluations in 71 s, best loss 263.0
```

The first line shows the Monte Carlo noise. At the true parameters, five
sets of 800 cells give losses between 3.5 and 24. A difference of that size
between two candidates means nothing unless both used the same cells.

With one fixed seed the loss is a deterministic function of the parameters.
The scan over `sm` is repeatable and has a clear valley. It is not smooth,
though: a small parameter change can move an event and change a cell's fate.
This is why the tutorial uses derivative-free optimizers. BiteOpt and
differential evolution handle such objectives; gradient methods and tight
convergence tolerances do not.

Two consequences:

- Do not trust the minimum of one seed. Refine with more cells, and check on
  a seed the optimizer never saw.
- With few cells, the minimum sits near the truth but where that seed's
  cells happen to fit best. Stage 1 above landed at `sm = 1.60`, where this
  seed's scan is lowest.

## 5. Using all cores

An evaluation can use cores in two ways:

- **Threads inside one simulation.** `threads=k` in `rules.ensemble` runs the
  cells of one evaluation in parallel. It needs nothing from the optimizer,
  and results do not depend on `k`.
- **Parallel candidates.** fcmaes evaluates several candidates at once in
  worker processes: `de.minimize(..., workers=k)`, or parallel restarts with
  `retry.minimize`. Each worker runs its simulations with `threads=1`.

On this two-CPU host both gave the same throughput: about 2.2 s per
evaluation with 2 threads in sequence, and 32 evaluations in 71 s with 2
workers × 1 thread. The DE run only measures throughput: 32 evaluations
are far too few for DE to converge, hence its best loss of 263. A real DE
fit needs a few thousand evaluations. On a larger machine the choice depends on the cells per
evaluation:

- With thousands of cells, threads inside the simulation scale well.
- With a few hundred cells, the fixed costs of preparing and collecting an
  ensemble weigh more, and parallel candidates scale better.
- The two can be mixed, for example 4 workers × 4 threads on 16 cores.

Each worker process prepares its own template on first use, through the
cached `_LOCAL` dictionary in the script, because prepared ensembles cannot
be sent between processes. The script relies on Linux's `fork` to give the
workers the observed data. On macOS and Windows, which spawn workers, load
the data from a file inside the worker instead.

fcmaes logs improvements through `loguru`; `logger.remove()` silences it.

## 6. Choosing an optimizer and a budget

fcmaes offers, among others:

- **BiteOpt** (`bitecpp.minimize`): robust with few evaluations and noisy
  objectives. A good default for 2–6 parameters and 100–500 evaluations.
- **Differential evolution** (`de.minimize`): parallel population
  evaluation with `workers`. Suited to larger budgets.
- **CMA-ES** (`cmaescpp.minimize`) and parallel restarts
  (`retry.minimize`): for more parameters or multimodal problems.

Practical sizes for this model:

- Start with 100 evaluations per parameter pair at a few hundred cells.
  Refine with 3–4 times the cells and a third of the evaluations.
- Use enough cells that a bin's expected count is at least about 10. For a
  bin holding 1% of cells, that means 1000 cells.
- Put an event budget on candidates (`max_events` in `rules.ensemble`) and a
  `wall_time_s` on the whole run. Extreme parameters, such as a very large
  `K` or `mu`, can otherwise make one evaluation very long. Treat a candidate
  that hits a budget as a high loss.

## 7. Making each evaluation cheap

- **Opt in to `mtdna_direct`** where the model has the three-class
  structure ([user manual, section 8](mtdna-user-manual.md#the-opt-in-fast-path-mtdna_direct)).
  Check `metadata["sampler"]`: a stop other than `N > 0 && mutants >= f * N`
  makes it fall back to Static Direct with inverse waits, which is about
  3.1× slower (user manual, section 8). Record
  other thresholds with first-passage trackers, or compute them from kept
  paths as here. Other models run on the generic methods, Static Direct
  first.
- **Keep only what the loss needs.** `keep_paths="all"` keeps 3 report
  times × 6 observables per cell, as one `result.paths` DataArray (cells ×
  report times × observables). Avoid event logs and Python callbacks
  inside an objective.
- **Reuse the template.** `bind_parameters` costs about 0.3 ms, as does a
  fresh request. Both are negligible, but a template guarantees that every
  candidate uses exactly the same settings.
- **Stop cells early when you can.** A takeover stop ends a cell's
  simulation. For takeover-only objectives this saves the remaining events
  of every stopped cell.

## 8. Beyond this example

- **More parameters.** Add `fracAdv`, `K` or the half-life (`delta`) to the
  vector. `bind_parameters` accepts any model parameter. `K` and `w0` must
  stay integers for `mtdna_direct`. Expect more confounding, and plan for
  more data summaries.
- **Likelihoods.** With counts of cells per heteroplasmy bin, a multinomial
  log-likelihood replaces the scaled squares. The same noise rules apply.
- **Other models.** The same loop works for the ATP model
  (`method="thinning"`) and for the clone model (tree method), at the costs
  given in the [performance analysis](mtdna-performance-analysis.md).
- **Rust.** For very large studies, an optimizer in Rust can call
  `crn_craft::rules::execute` directly, or a generated event loop
  ([user manual, section 16](mtdna-user-manual.md#16-rust)).

## 9. Checklist

1. **Simulate the data first.** Fit to data simulated with parameters you
   know, before fitting real data.
2. **Use summaries that separate the parameters.** Distributions separate
   them better than takeover curves alone.
3. **Fix the seed within a stage.** Change it between stages.
4. **Refine, then hold out.** Refine with more cells, then score on a held-out
   seed.
5. **Report the held-out loss.** Report it with the estimate, along with the
   spread of the loss at the estimate over a few seeds.
