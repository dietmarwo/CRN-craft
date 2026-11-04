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
        keep_paths=cells,
        max_output_bytes=256 << 20,
    )


def heteroplasmy(result):
    """Fraction of cells whose mutant fraction reaches each edge, per age."""
    names = result.metadata["observables"]
    values = np.array([path["values"] for path in result.metadata["paths"]])
    fraction = values[:, :, names.index("mutants")] / values[:, :, names.index("N")]
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
