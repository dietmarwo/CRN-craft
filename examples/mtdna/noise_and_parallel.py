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
