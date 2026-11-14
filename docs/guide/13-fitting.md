# 13. Fitting parameters to data

<!-- check: standalone, long -->

This chapter estimates two parameters of the three-class mtDNA model of
[chapter 12](12-mtdna.md): the mutation probability `mu` and the
replication advantage `sm` of advantaged mutants. They are fitted to
single-cell data with the derivative-free optimizers of
[fcmaes](https://github.com/dietmarwo/fast-cma-es). The objective is a
stochastic simulation, so the chapter also covers what that changes: Monte
Carlo noise, common random numbers, staged refinement, held-out checks and
how to spend the cores of a machine. The fitting scripts need two more
packages:

```sh
pip install crn-craft fcmaes scipy
```

Each block of this chapter is a complete script, the file named above it in
[`examples/mtdna/`](https://github.com/dietmarwo/CRN-craft/tree/main/examples/mtdna).
Run it from that folder or from the repository root. The outputs are from
one full run on 13 November 2026 on a two-CPU host. The guide's checker runs the
scripts in a short form, with few cells and three evaluations, when
`CRN_CRAFT_SMOKE=1` is set; the full run takes about thirteen minutes (4, 7
and 2 minutes for the three scripts).

## The problem

Applies: [chapter 12, The model](12-mtdna.md#the-model).

Single-cell measurements, such as COX staining or single-cell sequencing,
give the mutant fraction of individual cells at a few ages. The model
predicts, for any `mu` and `sm`, the distribution of mutant fractions
across cells at each age. Fitting means finding the parameters whose
predicted distributions match the observed ones.

Each prediction is an ensemble of simulated cells. A three-year cell at
`K = 1000` has about 150,000 events. Each script builds the model with the
`three_class` function of chapter 12 and runs it on the fast path,
`method="mtdna_direct"`, which certifies the builder model
([page 12d](12d-mtdna-methods.md#the-fast-path)). One objective evaluation
with 800 cells then takes about 2 s on two cores, so a hundred evaluations
take between three and four minutes. Under Static Direct, the same
evaluation takes about 3.6 times as long
([chapter 12, Conventions](12-mtdna.md#conventions)).

To be able to check the answer, the chapter makes its own data from known
parameters (`mu = 3·10⁻⁵`, `sm = 1.5`), with a seed that the fit never sees.
With real data, replace the observations by the measured fractions.

## A first attempt with takeover fractions

Applies: [page 12a, Takeover fraction and time to threshold](12a-mtdna-takeover.md#takeover-fraction-and-time-to-threshold).

Takeover data are often all that is available: the fraction of cells past a
threshold, for example 90% mutants, at a few ages. The first script fits
`mu` and `sm` to the takeover fractions at 360, 720 and 1080 days, plus the
mean mutant fraction. It runs 100 BiteOpt evaluations:

<!-- file: examples/mtdna/fit_takeover.py -->
```python
"""Fit the mutation probability and the replication advantage to takeover data."""

import math
import os
import sys
import time

import numpy as np
import crn_craft as crn
from fcmaes import bitecpp
from scipy.optimize import Bounds

SMOKE = os.environ.get("CRN_CRAFT_SMOKE") == "1"  # CI: few cells, 3 evaluations


def cells(count):
    return max(count // 100, 16) if SMOKE else count


def three_class(
    K=1000, w0=1000, m0=0, mu=3e-5, sm=1.5, half_life_days=10, takeover=None
):
    """The three-class mtDNA model: wild type W, neutral M, advantaged A."""
    model = crn.Model("mtdna_three_class")
    model.add_parameter(
        {
            "K": K,
            "delta": math.log(2) / (24 * half_life_days),
            "sw": 1.0,
            "sm": sm,
            "mu": mu,
            "fracAdv": 0.2,
        }
    )
    model.add_species({"W": w0, "M": m0, "A": 0})
    free = "max(0, K - W - M - A) / (W + sm * (M + A))"
    model.add_reaction(
        f"{free} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt"
    )
    model.add_reaction(
        f"{free} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt"
    )
    model.add_reaction(
        f"{free} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv"
    )
    model.add_reaction(f"{free} * sw * M", ["M"], ["M", "M"], "replication_mt")
    model.add_reaction(f"{free} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
    for label, species in (("wt", "W"), ("mt", "M"), ("mtAdv", "A")):
        model.add_reaction(f"delta * {species}", [species], [], f"degradation_{label}")
    if takeover is not None:
        model.set_stop(f"W + M + A > 0 && M + A >= {takeover} * (W + M + A)")
    return model


DAYS = np.array([360, 720, 1080])
MODEL = three_class(takeover=0.9)


def template(count, seed, threads=2):
    return crn.ensemble(
        MODEL,
        end=25_920.0,
        report_times=24.0 * np.concatenate(([0], DAYS)),
        replicates=cells(count),
        seed=seed,
        threads=threads,
        method="mtdna_direct",  # the fast path (guide page 12d)
    )


def summaries(result):
    """Fraction of cells taken over by each age, and mean mutant fraction."""
    stop = result.records["stop_time"]
    stopped = result.records["termination"] == "stopped"
    taken = np.array([(stopped & (stop <= 24.0 * d)).mean() for d in DAYS])
    mean = result.statistics["mean"].iloc[1:]  # the report times after day 0
    load = ((mean["M"] + mean["A"]) / mean.sum(axis=1)).to_numpy()
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
    evaluations = 3 if SMOKE else evaluations
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
data: 10.2 s
BiteOpt: 100 evaluations in 202 s, loss 0.90, mu 1.87e-05, sm 2.269
fitted held-out takeover [0.0335 0.083  0.133 ] load [0.0374 0.0823 0.1264]
true   held-out takeover [0.0112 0.0605 0.1078] load [0.0226 0.0682 0.1097]
```

The fit reaches a low loss on its own 800 cells, with `mu` more than a
third too small (1.9·10⁻⁵) and `sm` far too large (2.27). On 4000 fresh
cells its parameters predict 13.3% of cells taken over by day 1080, where
the data and the true parameters have 10.8%. Two things went wrong:

- **Confounding.** Takeover fractions mainly measure how many lineages
  succeed per unit time. Fewer mutations with a larger advantage give
  nearly the same curve as more mutations with a smaller one.
- **Overfitting to noise.** With one fixed set of 800 cells, the optimizer
  can find parameters that fit those cells' accidents. The held-out check
  shows it: the fit was best on its own cells, not on fresh ones.

The remedy for the first is more informative data. The remedy for the
second is more cells in a second stage and an independent check.

## A better objective with heteroplasmy distributions

Applies: [chapter 4, Single trajectories](04-results.md#single-trajectories).

The distribution of mutant fractions separates the two parameters:

- `mu` sets how many cells carry any mutants at all;
- `sm` sets how fast they move towards high mutant fractions.

The objective compares, at each age, the fraction of cells whose mutant
fraction reaches 1%, 10%, 50% and 90%. The mutant fraction of a cell is
`(M + A) / (W + M + A)`, read from the kept paths. Each difference is scaled
by its binomial standard error:

<!-- file: examples/mtdna/fit_heteroplasmy.py -->
```python
"""Fit the mutation probability and the replication advantage to
single-cell heteroplasmy distributions at three ages."""

import math
import os
import sys
import time

import numpy as np
import crn_craft as crn
from fcmaes import bitecpp
from scipy.optimize import Bounds

SMOKE = os.environ.get("CRN_CRAFT_SMOKE") == "1"  # CI: few cells, 3 evaluations


def cells(count):
    return max(count // 100, 16) if SMOKE else count


def three_class(
    K=1000, w0=1000, m0=0, mu=3e-5, sm=1.5, half_life_days=10, takeover=None
):
    """The three-class mtDNA model: wild type W, neutral M, advantaged A."""
    model = crn.Model("mtdna_three_class")
    model.add_parameter(
        {
            "K": K,
            "delta": math.log(2) / (24 * half_life_days),
            "sw": 1.0,
            "sm": sm,
            "mu": mu,
            "fracAdv": 0.2,
        }
    )
    model.add_species({"W": w0, "M": m0, "A": 0})
    free = "max(0, K - W - M - A) / (W + sm * (M + A))"
    model.add_reaction(
        f"{free} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt"
    )
    model.add_reaction(
        f"{free} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt"
    )
    model.add_reaction(
        f"{free} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv"
    )
    model.add_reaction(f"{free} * sw * M", ["M"], ["M", "M"], "replication_mt")
    model.add_reaction(f"{free} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
    for label, species in (("wt", "W"), ("mt", "M"), ("mtAdv", "A")):
        model.add_reaction(f"delta * {species}", [species], [], f"degradation_{label}")
    if takeover is not None:
        model.set_stop(f"W + M + A > 0 && M + A >= {takeover} * (W + M + A)")
    return model


DAYS = np.array([360, 720, 1080])
EDGES = np.array([0.01, 0.1, 0.5, 0.9])  # mutant fraction thresholds
MODEL = three_class(takeover=0.9)
TRUE = {"mu": 3e-5, "sm": 1.5}


def template(count, seed, threads=2):
    return crn.ensemble(
        MODEL,
        end=25_920.0,
        report_times=24.0 * DAYS,
        replicates=cells(count),
        seed=seed,
        threads=threads,
        method="mtdna_direct",  # the fast path (guide page 12d)
        keep_paths="all",
        max_output_bytes=256 << 20,
    )


def heteroplasmy(result):
    """Fraction of cells whose mutant fraction reaches each edge, per age."""
    paths = result.paths  # cells x report times x species
    mutants = paths.sel(variable="M") + paths.sel(variable="A")
    fraction = (mutants / paths.sum("variable")).to_numpy()
    return (fraction[:, :, None] >= EDGES).mean(axis=0)  # ages x edges


def loss(simulated, observed, count):
    # Binomial standard errors of the observed fractions, with a floor.
    se = np.sqrt(np.maximum(observed * (1 - observed), 0.005) / count)
    return float(np.sum(((simulated - observed) / se) ** 2))


SEARCH_CELLS = cells(800)
SEARCH = template(800, seed=7)  # one seed for every candidate


def objective(x):
    candidate = SEARCH.bind_parameters({"mu": 10.0 ** x[0], "sm": x[1]})
    return loss(heteroplasmy(candidate.simulate()), OBSERVED, SEARCH_CELLS)


if __name__ == "__main__":
    start = time.perf_counter()
    OBSERVED = heteroplasmy(template(4000, seed=2024).bind_parameters(TRUE).simulate())
    print(f"observed ({time.perf_counter() - start:.0f} s), rows are days {DAYS}:")
    print(np.round(OBSERVED, 4))
    evaluations = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    evaluations = 3 if SMOKE else evaluations
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
    SEARCH_CELLS = cells(2400)
    SEARCH = template(2400, seed=8)
    start = time.perf_counter()
    lower = [found.x[0] - 0.2, max(1.0, found.x[1] - 0.3)]
    upper = [found.x[0] + 0.2, found.x[1] + 0.3]
    refined = bitecpp.minimize(
        objective,
        bounds=Bounds(lower, upper),
        x0=found.x,
        max_evaluations=max(evaluations // 3, 1),
        rg=np.random.default_rng(43),
    )
    first = {"mu": mu, "sm": sm}
    mu, sm = 10.0 ** refined.x[0], refined.x[1]
    print(
        f"stage 2: {refined.nfev} evaluations in {time.perf_counter() - start:.0f} s, "
        f"loss {refined.fun:.1f}, mu {mu:.2e}, sm {sm:.3f}"
    )
    check = template(4000, seed=99)
    for label, values in (
        ("stage 1", first),
        ("stage 2", {"mu": mu, "sm": sm}),
        ("true", TRUE),
    ):
        held = heteroplasmy(check.bind_parameters(values).simulate())
        print(f"{label:7s} held-out loss {loss(held, OBSERVED, cells(4000)):7.1f}")
```

```text
observed (10 s), rows are days [ 360  720 1080]:
[[0.0648 0.036  0.025  0.0125]
 [0.122  0.0828 0.0712 0.062 ]
 [0.1798 0.1352 0.1212 0.1078]]
stage 1: 100 evaluations in 200 s, loss 1.5, mu 2.99e-05, sm 1.601
stage 2: 33 evaluations in 201 s, loss 5.7, mu 3.06e-05, sm 1.475
stage 1 held-out loss    68.6
stage 2 held-out loss     9.8
true    held-out loss     6.8
```

How the script is built:

- **One prepared template.** `template(cells, seed)` prepares the ensemble
  once. Every candidate is `SEARCH.bind_parameters({...}).simulate()`, which
  binds new values without rebuilding the request
  ([chapter 5](05-speed.md#parameter-sweeps)).
- **Common random numbers.** All candidates of a stage share one seed, so
  every candidate sees the same random streams. Differences between
  candidates are then due to the parameters, not to different random cells
  ([chapter 5](05-speed.md#common-random-numbers)).
- **A log scale for `mu`.** The optimizer works on `log10(mu)` in
  [−5.5, −3.5], which spans two orders of magnitude evenly.
- **Two stages.** Stage 1 (800 cells, 100 evaluations) finds the region.
  Stage 2 (2400 cells, a new seed, 33 evaluations) refines inside a box
  around it: ±0.2 in `log10(mu)` and ±0.3 in `sm`.
- **A held-out check.** 4000 fresh cells with a third seed score stage 1,
  stage 2 and the true parameters. Only this number says how good a fit is.

Stage 2 recovers `mu = 3.06·10⁻⁵` and `sm = 1.475`, against the true
3·10⁻⁵ and 1.5. Stage 1 found `mu` already, but `sm` too large (1.60). Its
held-out loss (68.6) shows it, while stage 2 scores 9.8, close to the true
parameters' 6.8. Even the true parameters do not reach a loss of zero: the
data are themselves a sample of 4000 cells, and the held-out cells are
another.

## Noise and the shape of the objective

Applies: [chapter 8, What makes a run reproducible](08-reproducibility.md#what-makes-a-run-reproducible).

The third script imports the second and looks at its objective:

<!-- file: examples/mtdna/noise_and_parallel.py -->
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
        max_evaluations=3 if fit.SMOKE else 32,
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
DE, 2 workers x 1 thread: 32 evaluations in 65 s, best loss 263.0
```

The first line shows the Monte Carlo noise. At the true parameters, five
sets of 800 cells give losses between 3.5 and 24.1. A difference of that
size between two candidates means nothing unless both used the same cells.

With one fixed seed the loss is a deterministic function of the parameters.
The scan over `sm` is repeatable and has a valley. It is not smooth,
though: a small parameter change can move an event and change a cell's
fate. This is why the chapter uses derivative-free optimizers. BiteOpt and
differential evolution handle such objectives; gradient methods and tight
convergence tolerances do not.

Two consequences:

- Do not trust the minimum of one seed. Refine with more cells, and check on
  a seed the optimizer never saw.
- With few cells, the minimum sits near the truth but where that seed's
  cells happen to fit best: with `mu` at its true value, this seed's scan is
  lowest at `sm = 1.60`, not 1.50. Stage 1 of the second script, with the
  same seed, landed there too.

## Using all cores

Applies: [chapter 5, Threads](05-speed.md#threads).

An evaluation can use cores in two ways:

- **Threads inside one simulation.** `threads=k` in `crn.ensemble` runs the
  cells of one evaluation in parallel. It needs nothing from the optimizer,
  and results do not depend on `k`.
- **Parallel candidates.** fcmaes evaluates several candidates at once in
  worker processes: `de.minimize(..., workers=k)`, or parallel restarts with
  `retry.minimize`. Each worker runs its simulations with `threads=1`.

On this two-CPU host both gave about the same throughput: about 2 s per
evaluation with 2 threads in sequence (stage 1 above), and 32 evaluations
in 65 s with 2 workers × 1 thread. The DE run only measures throughput:
32 evaluations are far too few for DE to converge, hence its best loss of
263.0. A real DE fit needs a few thousand evaluations. On a larger machine the choice depends on the cells per
evaluation:

- With thousands of cells, threads inside the simulation scale well.
- With a few hundred cells, the fixed costs of preparing and collecting an
  ensemble weigh more, and parallel candidates scale better.
- The two can be mixed, for example 4 workers × 4 threads on 16 cores.

Each worker process prepares its own template on first use, through the
cached `_LOCAL` dictionary in the script, because prepared ensembles cannot
be sent between processes. The script relies on Linux's `fork` to give the
workers the observed data. On macOS and Windows, which spawn workers, load
the data from a file inside the worker instead. fcmaes logs improvements
through `loguru`; `logger.remove()` silences it.

## Choosing an optimizer and a budget

Applies: [chapter 3, Budgets](03-running.md#budgets).

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
- Put an event budget on candidates (`max_events` in `crn.ensemble`) and a
  `wall_time_s` on the whole run. Extreme parameters, such as a very large
  `K` or `mu`, can otherwise make one evaluation very long. Treat a
  candidate that hits a budget as a high loss.

## Making each evaluation cheap

Applies: [chapter 5, What output costs](05-speed.md#what-output-costs).

- **Keep only what the loss needs.** Here that is three report times of
  three species per cell. Avoid event logs and callbacks inside an
  objective.
- **Stop cells early.** A takeover stop ends a cell's simulation and saves
  its remaining events.
- **Reuse the template**, so that every candidate uses exactly the same
  settings ([page 12e](12e-mtdna-large-studies.md#parameter-sweeps)).
- **A compiled model does not help here.** A loop is compiled for one
  parameter set, and every candidate has new parameters
  ([chapter 5, Limits](05-speed.md#limits)).
- **The fast path does.** The scripts pass `method="mtdna_direct"` with the
  builder model; it ran about 3.6 times as fast as Static Direct on the
  same model in the measurement of 13 November 2026
  ([page 12d](12d-mtdna-methods.md#the-fast-path)). It serves only the
  three-class structure: for another model, or a stop that the fast path
  does not certify, use Static Direct and check `metadata["sampler"]`.

## Beyond this example

Applies: [chapter 5, Parameter sweeps](05-speed.md#parameter-sweeps).

- **More parameters.** Add `fracAdv`, `K` or the half-life (`delta`) to the
  vector: `bind_parameters` accepts any parameter of the model. Expect more
  confounding, and plan for more data summaries.
- **Likelihoods.** With counts of cells per heteroplasmy bin, a multinomial
  log-likelihood replaces the scaled squares. The same noise rules apply.
- **Other models.** The same loop works for the ATP model of
  [page 12c](12c-mtdna-atp.md), at the cost of `hybrid`.
- **Generic lessons.** Four
  [optimization tutorials](https://github.com/dietmarwo/CRN-craft/tree/main/tutorials/optimization)
  calibrate other models with the same tools.
- **Rust.** For very large studies, an optimizer in Rust can call the
  engine directly ([Rust API](../developer/rust-api.md)).

## Checklist

Applies: [chapter 8, What to record for a paper](08-reproducibility.md#what-to-record-for-a-paper).

1. **Simulate the data first.** Fit to data simulated with parameters you
   know, before fitting real data.
2. **Use summaries that separate the parameters.** Distributions separate
   them better than takeover curves alone.
3. **Fix the seed within a stage.** Change it between stages.
4. **Refine, then hold out.** Refine with more cells, then score on a
   held-out seed.
5. **Report the held-out loss.** Report it with the estimate, along with the
   spread of the loss at the estimate over a few seeds.

See also: [Speed](05-speed.md),
[Large studies](12e-mtdna-large-studies.md),
[mtDNA performance analysis](../evidence/mtdna-performance-analysis.md).
