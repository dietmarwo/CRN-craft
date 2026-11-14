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
