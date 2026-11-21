"""Review checks of a finished full run: exact variance, timing, ensembles.

Run after `study.py --size full` and
`study.py --size full --free-final-degradation` (about 1 and 2 minutes
with Numba). Writes the `results/review-*.csv` files.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

from model import accept_distance, predictions
from study import FREE, ROOT, data, ess

THRESHOLD = 40.0
AGES = np.array([10.0, 20.0, 35.0, 50.0, 70.0, 100.0])


def states(run):
    """Stored walker states per ensemble after the first quarter."""
    out = []
    for path in sorted(run.glob("ensemble-360000-*.npz")):
        with np.load(path, allow_pickle=False) as saved:
            s = saved["states"]
        out.append(s[len(s) // 4 :])
    return out


def family_step(n, birth, death, hours, rng):
    """Exact linear birth-death law: zero or geometric per founder copy."""
    if hours <= 0:
        return n
    net = birth - death
    if net == 0:
        success = 1 / (1 + birth * hours)
        survival = success
    else:
        success = net / (birth * math.exp(net * hours) - death)
        survival = math.exp(net * hours) * success
    alive = rng.binomial(n, min(1.0, survival))
    extra = rng.negative_binomial(np.maximum(alive, 1), min(1.0, success))
    return alive + np.where(alive > 0, extra, 0)


def lineages(theta, ages, n, rng, h0=0.2):
    """Exact six-phase lineages; observation at a division sees the daughter."""
    ends = [203.0, 315.0, *(315.0 + 24 * np.cumsum(theta[13:16]))]
    divisions = [7.0 * (k + 1) for k in range(29)]
    divisions += [203.0 + 16.0 * (k + 1) for k in range(7)]
    stop = 24 * max(ages)
    cuts = sorted({0.0, *ends, *divisions, *(24 * a for a in ages)})
    cuts = [c for c in cuts if c <= stop]
    mutant = int(round(round(theta[0]) * h0))
    counts = np.tile([round(theta[0]) - mutant, mutant], (n, 1))
    seen, phase = {0.0: counts.copy()}, 0
    for a, b in zip(cuts[:-1], cuts[1:], strict=True):
        while phase < 5 and a >= ends[phase]:
            phase += 1
        death = theta[7 + phase]
        growth = -death if phase == 5 else theta[1 + phase]
        counts = family_step(counts, growth + death, death, b - a, rng)
        if b in divisions:
            counts = rng.binomial(counts, 0.5)
        seen[b] = counts.copy()
    rows = []
    for age in ages:
        c = seen[24 * age]
        total = c.sum(axis=1)
        alive = total > 0
        h = c[alive, 1] / total[alive]
        scale = h.mean() * (1 - h.mean())
        rows.append(
            {
                "mean_copies_simulated": total.mean(),
                "vprime_exact": h.var(ddof=1) / scale,
                "vprime_se": np.sqrt(np.var((h - h.mean()) ** 2) / alive.sum()) / scale,
                "extinct_fraction": 1 - alive.mean(),
                "fixed_or_lost_fraction": np.mean((h == 0) | (h == 1)),
            }
        )
    return rows


def distance(theta, ages, kinds, measured, vprime=None):
    theory = predictions(theta, ages)[0]
    if vprime is not None:
        theory[:, 1] = vprime
    picked = theory[np.arange(len(kinds)), kinds]
    weight = np.where(kinds == 1, 1000.0, 1.0)
    return float(np.sum(weight * (picked - measured) ** 2))


def exact_variance(flat, rng, n=20000):
    """E1: first-order normalized variance against exact lineages."""
    _, age, kinds, measured = data()
    pick = rng.choice(len(flat), 24, replace=False)
    late = np.array([predictions(x, np.array([100.0]))[0][0, 0] for x in flat[::50]])
    lowest = np.argsort(late)[:4] * 50
    grid = np.arange(0.25, 30.0, 0.25)
    rows = []
    for label, chosen in (("random", pick), ("lowest_day100_copies", lowest)):
        for j in chosen:
            theta = flat[j]
            minimum = float(grid[np.argmin(predictions(theta, grid)[0][:, 0])])
            ages = np.array(sorted({*AGES, minimum, *age[kinds == 1]}))
            sim = lineages(theta, ages, n, rng)
            theory = predictions(theta, ages)[0]
            lookup = {a: k for k, a in enumerate(ages)}
            vp = np.array(
                [
                    sim[lookup[a]]["vprime_exact"] if k else np.nan
                    for a, k in zip(age, kinds, strict=True)
                ]
            )
            se = math.sqrt(
                sum(
                    (
                        2000
                        * (sim[lookup[a]]["vprime_exact"] - m)
                        * sim[lookup[a]]["vprime_se"]
                    )
                    ** 2
                    for a, k, m in zip(age, kinds, measured, strict=True)
                    if k == 1
                )
            )
            first = distance(theta, age, kinds, measured)
            exact = distance(theta, age, kinds, measured, vp)
            for a in ages:
                k = lookup[a]
                rows.append(
                    {
                        "draw": int(j),
                        "set": label,
                        "day": a,
                        "copy_minimum": a == minimum,
                        "mean_copies": math.exp(theory[k, 0]),
                        **sim[k],
                        "vprime_first_order": theory[k, 1],
                        "relative_error": theory[k, 1] / sim[k]["vprime_exact"] - 1,
                        "distance_first_order": first,
                        "distance_exact": exact,
                        "distance_se": se,
                    }
                )
    return pd.DataFrame(rows)


def age_shift(flat, rng):
    """E2: the three ages restored to division times, against rounded labels."""
    _, age, kinds, measured = data()
    rounded = age.copy()
    for exact, label in ((7 / 24, 0.29), (14 / 24, 0.58), (21 / 24, 0.88)):
        rounded[np.isclose(age, exact)] = label
    order = np.argsort(rounded, kind="stable")
    sub = flat[rng.choice(len(flat), 20000, replace=False)]
    base = np.array([accept_distance(x, age, kinds, measured) for x in sub])
    moved = np.array(
        [accept_distance(x, rounded[order], kinds[order], measured[order]) for x in sub]
    )
    return pd.DataFrame(
        [
            {
                "states": len(sub),
                "max_distance_division_times": base.max(),
                "accepted_with_rounded_labels": np.mean(moved <= THRESHOLD),
                "median_change": np.median(moved - base),
                "change_q05": np.quantile(moved - base, 0.05),
                "change_q95": np.quantile(moved - base, 0.95),
            }
        ]
    )


def rhat(chains):
    """Rank-normalized split R-hat, maximum of bulk and folded."""
    n = min(len(c) for c in chains) // 2 * 2
    x = np.array([h for c in chains for h in (c[: n // 2], c[n // 2 : n])])

    def single(y):
        z = norm.ppf((rankdata(y).reshape(y.shape) - 0.375) / (y.size + 0.25))
        within = z.var(axis=1, ddof=1).mean()
        between = z.shape[1] * z.mean(axis=1).var(ddof=1)
        k = z.shape[1]
        return math.sqrt(((k - 1) / k * within + between / k) / within)

    return max(single(x), single(np.abs(x - np.median(x))))


def ensembles(per):
    """E3 and E5: agreement of the four ensembles and their effective sizes."""
    selected = np.array([10.0, 35.0, 70.0, 100.0])
    names = [f"theta{k}" for k in FREE]
    names += [f"{o}_{a:g}" for a in selected for o in ("log_copies", "vprime")]
    values, sizes = [], []
    for s in per:
        flat = s.reshape(-1, s.shape[-1])
        obs = np.asarray([predictions(x, selected)[0] for x in flat])
        obs = obs.reshape(len(s), s.shape[1], -1)
        values.append(np.concatenate([s[:, :, FREE], obs], axis=2))
        series = obs.mean(axis=1)
        marginal = obs.reshape(-1, obs.shape[2]).var(axis=0)
        sizes.append(
            [ess(v) * m / v.var() for v, m in zip(series.T, marginal, strict=True)]
        )
    rows = []
    for j, name in enumerate(names):
        chains = [v[::10, w, j] for v in values for w in range(v.shape[1])]
        row = {"quantity": name, "rhat_walkers": rhat(chains)}
        for e, v in enumerate(values):
            q = np.quantile(v[:, :, j], [0.025, 0.5, 0.975])
            row |= {f"q025_{e}": q[0], f"median_{e}": q[1], f"q975_{e}": q[2]}
            if j >= len(FREE):
                row[f"ess_{e}"] = sizes[e][j - len(FREE)]
        rows.append(row)
    return pd.DataFrame(rows)


def final_degradation(run, free_run):
    """E4: phase-six degradation fixed to zero against free on [0, 1]."""
    rows = []
    for label, folder in (("fixed_zero", run), ("free", free_run)):
        f = pd.read_csv(folder / "figure3-comparison.csv")
        hb = pd.read_csv(folder / "hb-comparison.csv")
        abc = pd.read_csv(folder / "abc.csv")
        flat = np.concatenate([s.reshape(-1, 16) for s in states(folder)])
        late = f[(f.observable == "copies") & (f.day == 100)].iloc[0]
        nu = np.quantile(flat[:, 12], [0.025, 0.5, 0.975])
        start = np.quantile(13.125 + flat[:, 13:16].sum(axis=1), [0.025, 0.5, 0.975])
        rows.append(
            {
                "phase6_degradation": label,
                "acceptance_min": abc.acceptance.min(),
                "acceptance_max": abc.acceptance.max(),
                "ess_sum_of_minima": abc.ess_min_observable.sum(),
                "mad_log_copies": f[f.observable == "copies"].difference.abs().median(),
                "mad_vprime": f[f.observable != "copies"].difference.abs().median(),
                "day100_copies_median": late["median"],
                "day100_copies_q025": late.q025,
                "day100_copies_q975": late.q975,
                "hb_in_sample_band": int(hb.in_sample_band.sum()),
                "nu6_q025": nu[0],
                "nu6_median": nu[1],
                "nu6_q975": nu[2],
                "phase6_start_day_median": start[1],
                "phase6_start_day_q025": start[0],
                "phase6_start_day_q975": start[2],
            }
        )
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, default=ROOT / "runs/full")
    p.add_argument("--free-run", type=Path, default=ROOT / "runs/full-free-final")
    p.add_argument("--output", type=Path, default=ROOT / "results")
    a = p.parse_args()
    per = states(a.run)
    flat = np.concatenate([s.reshape(-1, 16) for s in per])
    rng = np.random.default_rng(2101)
    a.output.mkdir(parents=True, exist_ok=True)
    age_shift(flat, rng).to_csv(a.output / "review-age-shift.csv", index=False)
    ensembles(per).to_csv(a.output / "review-ensembles.csv", index=False)
    final_degradation(a.run, a.free_run).to_csv(
        a.output / "review-final-degradation.csv", index=False
    )
    exact_variance(flat, rng).to_csv(a.output / "review-vprime.csv", index=False)


if __name__ == "__main__":
    main()
