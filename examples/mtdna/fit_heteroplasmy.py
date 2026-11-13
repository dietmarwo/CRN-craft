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
        method="static_direct",
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
    for label, values in (("stage 1", first), ("stage 2", {"mu": mu, "sm": sm}),
                          ("true", TRUE)):
        held = heteroplasmy(check.bind_parameters(values).simulate())
        print(f"{label:7s} held-out loss {loss(held, OBSERVED, cells(4000)):7.1f}")
