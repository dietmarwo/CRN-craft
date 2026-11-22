"""Review experiments for the Figure 2E speed gap and the Figure 1C settings.

Each subcommand writes rows to results/review-experiments.csv (or --output).
Targets and decision rules are stated in the guide page (tutorial S4a).

  deterministic  mean-field chain without noise, and the mean of the step scheme
  cohort         approximate fibers (wave.simulate), checkpointed per 16 fibers
  noise          exact engine vs approximate scheme, local variance after 20 days
  equal          per-region spread at day 700 for the equal-density condition
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from study import MANIFEST, ROOT, chain_initial, fractions, paths

ARCHIVED = 0.18973014373038566
SPACING = 30.0


def front(x, h, t):
    """Descending crossing by linear interpolation, NaN if not unique."""
    i = np.flatnonzero((h[:-1] >= t) & (h[1:] < t))
    if len(i) != 1:
        return np.nan
    i = i[0]
    return x[i] + (x[i + 1] - x[i]) * (t - h[i]) / (h[i + 1] - h[i])


def smoothed(x, h, t):
    return front(x[1:-1], np.convolve(h, [0.25, 0.5, 0.25], mode="valid"), t)


def chain_rates(p, copies):
    mu, c, g, d = p["death"], p["c"], p["hop"], p["delta"]

    def flows(w, m, dt=None):
        birth = np.maximum(0.0, mu + c * (copies - w - d * m))
        out = []
        for x in (w, m):
            if dt is None:
                inflow = np.zeros_like(x)
                inflow[1:] += g * x[:-1]
                inflow[1:-1] += g * x[2:]
                outflow = 2 * g * x
                outflow[-1] = g * x[-1]
                change = (birth - mu) * x + inflow - outflow
            else:
                loss = -np.expm1(-dt * (mu + 2 * g))
                last = -np.expm1(-dt * (mu + g))
                gone = x * loss
                gone[-1] = x[-1] * last
                left = gone * g / (mu + 2 * g)
                left[-1] = gone[-1] * g / (mu + g)
                right = (gone - left) * g / (mu + g)
                change = x * birth * dt - gone
                change[2:] += right[1:-1]
                change[1:-1] += left[2:]
                change[1] += dt * g * x[0]
            change[0] = 0.0
            out.append(change)
        return out

    return flows


def deterministic():
    p = MANIFEST["fig2_site_parameters"]
    rows = []
    for copies in (3500, 3000):
        x, start = chain_initial({**p, "copies": copies})
        w0, m0 = start[:, 0].astype(float), start[:, 1].astype(float)
        h0 = m0 / (w0 + m0)
        flows = chain_rates(p, copies)

        def rhs(_, y, flows=flows, n=x.size):
            return np.concatenate(flows(y[:n], y[n:]))

        y = solve_ivp(
            rhs,
            (0, 1480),
            np.concatenate([w0, m0]),
            method="LSODA",
            rtol=1e-8,
            atol=1e-6,
        ).y[:, -1]
        cases = [("ode", None, y[: len(x)], y[len(x) :])]
        for dt in (0.5, 0.25, 0.125):
            w, m = w0.copy(), m0.copy()
            for _ in range(round(1480 / dt)):
                dw, dm = flows(w, m, dt)
                w, m = w + dw, m + dm
            cases.append(("mean-of-step-scheme", dt, w, m))
        for kind, dt, w, m in cases:
            h1 = m / (w + m)
            for t in (0.4, 0.5, 0.6):
                rows.append(
                    {
                        "experiment": "T1-deterministic",
                        "variant": f"{kind} Nss={copies}"
                        + ("" if dt is None else f" step={dt}"),
                        "threshold": t,
                        "value": (front(x, h1, t) - front(x, h0, t)) / 1480,
                        "unit": "um/day",
                    }
                )
    return rows


def cohort(name, copies, dt, fibers, seed, folder):
    from wave import simulate

    p = {**MANIFEST["fig2_site_parameters"], "copies": copies}
    x, start = chain_initial(p)
    h0 = start[:, 1] / start.sum(axis=1)
    folder = folder / name
    folder.mkdir(parents=True, exist_ok=True)
    for first in range(0, fibers, 16):
        file = folder / f"chunk{first:04d}.npy"
        if not file.exists():
            terminal = simulate(
                start,
                min(16, fibers - first),
                seed + first,
                round(1480 / dt),
                dt,
                p["copies"],
                p["delta"],
                p["hop"],
                p["death"],
                p["c"],
            )
            np.save(file.with_suffix(".tmp.npy"), terminal)
            file.with_suffix(".tmp.npy").replace(file)
    sample = np.concatenate([np.load(f) for f in sorted(folder.glob("chunk*.npy"))])
    n = len(sample)
    weights = np.random.default_rng(7).multinomial(n, np.full(n, 1 / n), 2000) / n
    rows = []
    for t in (0.4, 0.5, 0.6):
        start_front = smoothed(x, h0, t)
        boot = np.array(
            [(smoothed(x, w @ sample, t) - start_front) / 1480 for w in weights]
        )
        low, high = np.nanquantile(boot, [0.005, 0.995])
        rows.append(
            {
                "experiment": name,
                "variant": f"Nss={copies} step={dt}",
                "threshold": t,
                "fibers": n,
                "value": (smoothed(x, sample.mean(axis=0), t) - start_front) / 1480,
                "low_99": low,
                "high_99": high,
                "unit": "um/day",
                "label": label(low, high),
            }
        )
    area = (sample.sum(axis=1) - h0.sum()) * SPACING / 1480
    rows.append(
        {
            "experiment": name,
            "variant": f"Nss={copies} step={dt} area",
            "fibers": n,
            "value": area.mean(),
            "low_99": area.mean() - 2.576 * area.std(ddof=1) / np.sqrt(n),
            "high_99": area.mean() + 2.576 * area.std(ddof=1) / np.sqrt(n),
            "unit": "um/day",
        }
    )
    return rows


def label(low, high, margin=MANIFEST["fig2_half_height_margin_um_day"]):
    if ARCHIVED - margin <= low and high <= ARCHIVED + margin:
        return "numerical agreement"
    if high < ARCHIVED or low > ARCHIVED:
        return "unresolved discrepancy"
    return "inconclusive"


def noise(folder):
    from wave import simulate

    p = MANIFEST["fig2_site_parameters"]
    x, start = chain_initial(p)
    h0 = start[:, 1] / start.sum(axis=1)
    folder = folder / "noise"
    folder.mkdir(parents=True, exist_ok=True)
    exact = []
    for k in range(16):
        file = folder / f"exact{k:02d}.npy"
        if not file.exists():
            counts = paths(
                start,
                (p["copies"], p["delta"], p["hop"]),
                [0.0, 20.0],
                1,
                MANIFEST["seed"] + 11000 + k,
                1,
                death=p["death"],
                c=p["c"],
                fixed_left=True,
            )
            np.save(file, fractions(counts)[0, -1])
        exact.append(np.load(file))
    exact = np.array(exact)
    front_sites = (h0 > 0.05) & (h0 < 0.95)

    def spread(a):
        return a[:, front_sites].var(axis=0, ddof=1).mean()

    rng = np.random.default_rng(3)
    rows = []
    for dt in (0.5, 0.125):
        file = folder / f"approx{dt}.npy"
        if not file.exists():
            np.save(
                file,
                simulate(
                    start,
                    64,
                    MANIFEST["seed"] + 11500 + int(dt * 1000),
                    round(20 / dt),
                    dt,
                    p["copies"],
                    p["delta"],
                    p["hop"],
                    p["death"],
                    p["c"],
                ),
            )
        approx = np.load(file)
        boot = [
            spread(approx[rng.integers(0, 64, 64)])
            / spread(exact[rng.integers(0, 16, 16)])
            for _ in range(2000)
        ]
        low, high = np.quantile(boot, [0.005, 0.995])
        rows.append(
            {
                "experiment": "T4-noise",
                "variant": f"step={dt} vs exact, day 20",
                "fibers": 64,
                "value": spread(approx) / spread(exact),
                "low_99": low,
                "high_99": high,
                "unit": "variance ratio",
                "label": "no evidence of bias" if low <= 1 <= high else "biased",
            }
        )
    return rows


def equal():
    archived = {"neutral": (0.006117, 0.006090), "extra": (0.005973, 0.005990)}
    rows = []
    for copies, start in ((200, 100), (150, 75)):
        for treatment, extra in (("neutral", 0.0), ("extra", 0.00007)):
            counts = paths(
                [(start, start)] * 2,
                (copies, 1.0, 0.035),
                [0.0, 700.0],
                2000,
                MANIFEST["seed"] + 200 + int(extra > 0),
                1,
                extra=extra,
            )
            h = fractions(counts)[:, -1, :]
            ratio = (
                h.std(axis=0, ddof=1) / np.sqrt(2000) / np.array(archived[treatment])
            )
            for region, value in enumerate(ratio):
                rows.append(
                    {
                        "experiment": "E1-equal-density",
                        "variant": f"Nss={copies} start={start}/{start} {treatment} region {region}",
                        "fibers": 2000,
                        "value": value,
                        "unit": "SEM ratio fresh/archive",
                        "label": "match" if 0.95 <= value <= 1.05 else "mismatch",
                    }
                )
    return rows


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("deterministic", "cohort", "noise", "equal")
    )
    parser.add_argument("--work", type=Path, default=ROOT / "runs" / "review")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "results" / "review-experiments.csv"
    )
    args = parser.parse_args()
    os.environ.setdefault("NUMBA_NUM_THREADS", "1")
    if args.command == "deterministic":
        rows = deterministic()
    elif args.command == "noise":
        rows = noise(args.work)
    elif args.command == "equal":
        rows = equal()
    else:
        rows = []
        for name, copies, dt, fibers, offset in (
            ("T2-copies-3000", 3000, 0.125, 96, 9000),
            ("T3-step-0.0625", 3500, 0.0625, 32, 10000),
            ("T3b-step-0.5", 3500, 0.5, 32, 12000),
        ):
            rows += cohort(
                name, copies, dt, fibers, MANIFEST["seed"] + offset, args.work
            )
    table = pd.DataFrame(rows)
    if args.output.exists():
        old = pd.read_csv(args.output)
        old = old[~old.experiment.isin(table.experiment.unique())]
        table = pd.concat([old, table], ignore_index=True)
    table.to_csv(args.output, index=False, lineterminator="\n", float_format="%.6g")
    print(json.dumps({"rows": len(rows)}))


if __name__ == "__main__":
    main()
