"""Run a time-step refinement of a source-scale fiber ensemble.

Births are Poisson draws using start-of-step feedback. Degradation and
hopping compete for the same copies, so a copy can depart only once per
step. The boundary acts as a fixed reservoir. This solver is approximate.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit, prange

from study import MANIFEST, ROOT, chain_initial, crossing, speeds


@njit
def binomial_clip(count, p):
    """Integer binomial draw, including pure boundaries."""
    return np.random.binomial(count, p) if count > 0 and p > 0 else 0


@njit(parallel=True)
def simulate(initial, n, seed, steps, dt, copies, delta, hop, death, c):
    sites = len(initial)
    result = np.empty((n, sites))
    events = death + 2 * hop
    loss = -math.expm1(-dt * events)
    terminal_loss = -math.expm1(-dt * (death + hop))
    for k in prange(n):
        np.random.seed(seed + k)
        state = initial.copy()
        future = initial.copy()
        for _ in range(steps):
            future[:, :] = state
            for genotype in range(2):
                future[1, genotype] += np.random.poisson(dt * hop * state[0, genotype])
            for i in range(1, sites):
                rate = max(
                    0.0, death + c * (copies - state[i, 0] - delta * state[i, 1])
                )
                for genotype in range(2):
                    count = state[i, genotype]
                    if count == 0:
                        continue
                    last = i + 1 == sites
                    departures = binomial_clip(count, terminal_loss if last else loss)
                    left = binomial_clip(
                        departures, hop / (death + hop if last else events)
                    )
                    right = (
                        0
                        if last
                        else binomial_clip(departures - left, hop / (death + hop))
                    )
                    if i > 1:
                        future[i - 1, genotype] += left
                    if not last:
                        future[i + 1, genotype] += right
                    future[i, genotype] += (
                        np.random.poisson(count * rate * dt) - departures
                    )
            state, future = future, state
        for i in range(sites):
            total = state[i, 0] + state[i, 1]
            result[k, i] = state[i, 1] / total if total else 0.0
    return result


def summarize(position, starting, terminal, dt, seed):
    """Smooth three adjacent mean bins before testing a descending crossing."""
    rng = np.random.default_rng(seed)
    n = len(terminal)
    mean = terminal.mean(axis=0)

    def location(profile, h):
        return crossing(
            position[1:-1], np.convolve(profile, [0.25, 0.5, 0.25], mode="valid"), h
        )

    dates = 1480
    archive = speeds(pd.read_csv(ROOT / "reference-2e.csv"))
    rows = []
    for h in MANIFEST["fig2_crossings"]:
        start = location(starting, h)
        end = location(mean, h)
        estimates, ambiguous = [], 0
        for counts in rng.multinomial(n, np.full(n, 1 / n), size=700):
            profile = counts @ terminal / n
            try:
                estimates.append((location(profile, h) - start) / dates)
            except ValueError:
                ambiguous += 1
        low, high = (
            np.quantile(estimates, [0.005, 0.995]) if estimates else (np.nan, np.nan)
        )
        ref = archive.loc[archive.threshold == h, "speed_um_day"].iloc[0]
        rows.append(
            {
                "step_day": dt,
                "fibers": n,
                "threshold": h,
                "speed_um_day": (end - start) / dates,
                "archived_speed": ref,
                "low_99": low,
                "high_99": high,
                "ambiguous": ambiguous,
            }
        )
    return mean, rows


def run(n, steps, output):
    p = MANIFEST["fig2_site_parameters"]
    x, initial = chain_initial(p)
    starting = initial[:, 1] / initial.sum(axis=1)
    output.mkdir(parents=True, exist_ok=True)
    rows, profiles = [], [(0.0, 0, pos, h) for pos, h in zip(x, starting, strict=True)]
    for dt in steps:
        file = output / f"approx-wave-{dt}-{n}.npy"
        seed = MANIFEST["seed"] + 6000 + 1000 * round(-math.log2(dt) - 1)
        if file.is_file():
            terminal = np.load(file, allow_pickle=False)
        else:
            terminal = simulate(
                initial,
                n,
                seed,
                round(1480 / dt),
                dt,
                p["copies"],
                p["delta"],
                p["hop"],
                p["death"],
                p["c"],
            )
            pending = file.with_suffix(".tmp.npy")
            np.save(pending, terminal)
            pending.replace(file)
        profile, values = summarize(x, starting, terminal, dt, seed + 5000)
        profiles.extend((dt, 1480, pos, h) for pos, h in zip(x, profile, strict=True))
        rows.extend(values)
        print(json.dumps(values, sort_keys=True), flush=True)
    pd.DataFrame(rows).to_csv(output / "approx-wave-speeds.csv", index=False)
    pd.DataFrame(profiles, columns=["step_day", "day", "position_um", "mean"]).to_csv(
        output / "approx-wave-profiles.csv", index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fibers", type=int, default=64)
    parser.add_argument("--output", type=Path, default=ROOT / "runs" / "wave")
    parser.add_argument(
        "--step", type=float, choices=(0.5, 0.25, 0.125), action="append"
    )
    args = parser.parse_args()
    run(args.fibers, args.step or (0.5, 0.25), args.output)
