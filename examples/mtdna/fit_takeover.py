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
