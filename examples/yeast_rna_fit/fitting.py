"""Fit mean and histogram observations with the same physical model."""

import hashlib
from concurrent.futures import ThreadPoolExecutor
from time import perf_counter

import numpy as np
import pandas as pd
from fcmaes import decpp
from scipy.optimize import Bounds

from model import ROOT, compiled, exact_means, nonspatial, quick_means, spatial


def observations():
    return pd.read_csv(ROOT / "source-counts.csv", dtype={"replicate": str})


class Data:
    def __init__(self, frame, spatial_model=False):
        self.dose = sorted(frame.dose.unique())
        self.groups = []
        self.spatial = spatial_model
        for (dose, rep, minute), rows in frame.groupby(["dose", "replicate", "minute"]):
            if minute == 0:
                # Published fits start in S1 with zero RNA. Time zero is
                # displayed with the data and omitted from both objectives.
                continue
            weights = rows.cells.to_numpy()
            xs = (
                rows[["nuclear", "cytoplasm"]].to_numpy()
                if spatial_model
                else rows[["total"]].to_numpy()
            )
            n = int(weights.sum())
            mean = weights @ xs / n
            centered = xs - mean
            # Fixed empirical Gaussian covariance for the average of n cells.
            covariance = ((centered.T * weights) @ centered) / max(1, n * (n - 1))
            covariance += np.eye(xs.shape[1]) / (12 * n)
            inv = np.linalg.inv(covariance)
            self.groups.append((dose, rep, minute, n, xs, weights, mean, inv))

    def means_cost(self, predicts):
        score = 0.0
        for dose, _, minute, _, _, _, mean, inv in self.groups:
            error = mean - predicts[dose][minute]
            score += 0.5 * (error @ inv @ error)
        return float(score)

    def counts_cost(self, predicts, per_rep=False):
        grouped = {}
        for dose, rep, minute, n, xs, weights, _, _ in self.groups:
            law = predicts[dose][minute]
            if self.spatial:
                inside = (xs[:, 0] < law.shape[1]) & (xs[:, 1] < law.shape[0])
                pv = law[xs[inside, 1], xs[inside, 0]]
            else:
                inside = xs[:, 0] < len(law)
                pv = law[xs[inside, 0]]
            # Compare with explicit overflow as a coarse bin. References
            # verify cap expansion separately. A nonpositive occupied bin
            # invalidates the candidate at the current numerical resolution.
            overflow = weights[~inside].sum()
            if not np.all(pv > 0) or (overflow and law.sum() >= 1):
                return np.inf
            logloss = -np.dot(weights[inside], np.log(pv))
            if overflow:
                logloss -= overflow * np.log1p(-law.sum())
            initial = grouped.get((dose, rep), (0.0, 0))
            grouped[dose, rep] = (initial[0] + float(logloss), initial[1] + n)
        return grouped if per_rep else sum(x[0] for x in grouped.values())

    def n(self):
        return sum(x[3] for x in self.groups)


def nonspatial_laws(x, data, cap, step):
    return {dose: nonspatial(x, dose, cap, step).sum(-1) for dose in data.dose}


def spatial_laws(x, data, ncap, ccap, tolerance=2e-6):
    minutes = np.unique([0, *[t for _, _, t, *_ in data.groups]])
    out = {}
    for dose in data.dose:
        raw = spatial(x, dose, minutes, ncap, ccap, tolerance).sum(-1)
        out[dose] = dict(zip(minutes, raw, strict=True))
    return out


def exact_moments(x, data):
    return {
        dose: exact_means(x, dose, np.arange(56), data.spatial) for dose in data.dose
    }


def one_fit(cost, free, budget, seed):
    t0 = perf_counter()
    dimension = len(free)
    calls, invalid = 0, 0

    def evaluate(z):
        nonlocal calls, invalid
        calls += 1
        value = cost(z)
        if not np.isfinite(value):
            invalid += 1
            return 1e18
        return value

    initial = evaluate(np.zeros(dimension))
    bound = Bounds(-np.asarray(free), np.asarray(free))
    found = decpp.minimize(
        evaluate,
        bounds=bound,
        popsize=max(12, 2 * dimension),
        max_evaluations=budget,
        rg=np.random.default_rng(seed),
    )
    z = found.x if found.fun < initial else np.zeros(dimension)
    return z, {
        "initial": initial,
        "value": min(initial, found.fun),
        "evaluations": calls,
        "invalid": invalid,
        "seed": seed,
        "requested": budget,
        "seconds": perf_counter() - t0,
    }


def fit_nonspatial(data, anchor, preset, seed, *, workers, kind, extra=False):
    compiled()
    step, cap = preset["step"], preset["cap"]
    # The same multiplicative box is used for both objectives and preserves
    # the sign of the Hog1 coefficient from the author's parameter vector.
    bounds = np.log([3, 4, 4, 3, 3, 3, 3, 4, 3, 4, 4, 3, 1.6])

    def transform(z):
        x = anchor.copy()
        x[:13] *= np.exp(z)
        return x

    def cost(z):
        x = transform(z)
        if kind == "mean":
            return data.means_cost({s: quick_means(x, s, step) for s in data.dose})
        return data.counts_cost(nonspatial_laws(x, data, cap, step))

    budgets = [preset["evaluations"]] * preset["restarts"]
    if extra:
        budgets += [4 * preset["evaluations"]]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(
            pool.map(
                lambda job: one_fit(cost, bounds, job[1], seed + job[0]),
                enumerate(budgets),
            )
        )
    best_z, best = min(results, key=lambda r: r[1]["value"])
    return transform(best_z), best, [r[1] for r in results], transform


def refine_nonspatial(x, data, kind, cap, step):
    if kind == "mean":
        return data.means_cost(exact_moments(x, data))
    return data.counts_cost(nonspatial_laws(x, data, cap, step))


def spatial_grid(data, anchor, output, workers, budgets, seed):
    """Fit two spatial rates against a checked interpolation of sparse solves."""
    from scipy.interpolate import RectBivariateSpline

    width = np.log(1.35)
    coordinates = np.array([-width, 0.0, width])
    support = (21, 240)
    shape = (len(coordinates),) * 2
    key_data = (ROOT / "source-counts.csv").read_bytes()
    fingerprint = hashlib.sha256(
        anchor.tobytes() + np.asarray(support).tobytes() + key_data
    ).hexdigest()[:16]

    def transform(z):
        x = anchor.copy()
        x[[9, 11]] *= np.exp(z)
        return x

    def node(index):
        z = [coordinates[i] for i in index]
        x = transform(z)
        path = output / f"spatial-grid-{fingerprint}-{index[0]}-{index[1]}.npz"
        if path.exists():
            with np.load(path, allow_pickle=False) as bundle:
                values = bundle["value"]
        else:
            means = data.means_cost(exact_moments(x, data))
            histogram = data.counts_cost(spatial_laws(x, data, *support)) / data.n()
            values = np.array([means, histogram])
            pending = path.with_suffix(".tmp.npz")
            np.savez(pending, value=values)
            pending.replace(path)
        return index, values

    scores = np.zeros((*shape, 2))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for index, score in pool.map(node, list(np.ndindex(shape))):
            if np.any(~np.isfinite(score)):
                raise RuntimeError(
                    f"Invalid occupied bin on the declared spatial grid at {index}"
                )
            scores[index] = score
    grid = pd.DataFrame(
        [(coordinates[i], coordinates[j], *scores[i, j]) for i, j in np.ndindex(shape)],
        columns=["log_i3_factor", "log_decay_factor", "mean_cost", "nll_per_cell"],
    )
    grid.to_csv(output / "spatial-grid.csv", index=False)
    fits, rows, diagnostics = {}, [], []
    histogram_law = None
    for i, kind in enumerate(("mean", "histogram")):
        predicted_score = RectBivariateSpline(
            coordinates, coordinates, scores[..., i], kx=2, ky=2
        )

        def model_cost(z, interpolation=predicted_score):
            return float(interpolation(z[0], z[1], grid=False))

        trials = [
            one_fit(model_cost, [width, width], budget, seed + r)
            for r, budget in enumerate(budgets)
        ]
        z, _ = min(trials, key=lambda x: x[1]["value"])
        diagnostics.extend({"objective": kind, **r} for _, r in trials)
        x = transform(z)
        if kind == "mean":
            proposed = data.means_cost(exact_moments(x, data))
        else:
            proposed = data.counts_cost(spatial_laws(x, data, *support)) / data.n()
        node = np.unravel_index(np.argmin(scores[..., i]), shape)
        source = "optimized"
        if scores[node][i] < proposed:
            # Preserve the best evaluated candidate when the interpolant
            # exaggerates the gain at a point between grid nodes.
            z = np.array([coordinates[at] for at in node])
            x = transform(z)
            source = "grid"
        if source == "grid":
            coarse_hist = float(scores[node][1])
        elif kind == "histogram":
            coarse_hist = proposed
        else:
            coarse_hist = data.counts_cost(spatial_laws(x, data, *support)) / data.n()
        fine = spatial_laws(x, data, 28, 305, tolerance=4e-7)
        fine_hist = data.counts_cost(fine) / data.n()
        if kind == "histogram":
            histogram_law = fine
        fit_mean = data.means_cost(exact_moments(x, data))
        actual = fit_mean if kind == "mean" else fine_hist
        fits[kind] = x
        rows.append(
            {
                "objective": kind,
                "log_i3_factor": z[0],
                "log_decay_factor": z[1],
                "interpolated_objective": model_cost(z),
                "verified_selection": source,
                "optimizer_direct_coarse": proposed,
                "direct_objective": actual,
                "mean_cost": fit_mean,
                "direct_histogram_nll_per_cell": fine_hist,
                "coarse_refined_nats_per_cell": abs(fine_hist - coarse_hist),
                "max_physical_leakage_estimate": max(
                    1 - law.sum() for p in fine.values() for law in p.values()
                ),
            }
        )
    pd.DataFrame(rows).to_csv(output / "spatial-fits.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(output / "spatial-restarts.csv", index=False)
    return fits, histogram_law
