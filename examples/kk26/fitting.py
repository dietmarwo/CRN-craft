"""Preparation, calibration and uncertainty on synthetic cross-sectional data."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import Bounds, least_squares
from scipy.stats import qmc
from scipy.sparse.linalg import spsolve

from model import (
    AGES,
    BUILD_FLAGS,
    HERE,
    REGIMES,
    SAMPLER,
    SPLIT,
    cohort,
    forward,
    library,
    tails,
)

ROOT_SEED = 20261120
PARAMETERS = ("mutation", "advantage", "advantageous")
GRID = "local-cubic-spline/1"
CAPS = {"baseline": 192, "strong": 320}
STEP = 5.0
STORED = HERE / "results" / "reconstruction"
TOLERANCE = 0.001


def axes_for(regime):
    if regime == "baseline":
        return (
            np.log10([3.8e-6, 4.9e-6, 6.1e-6, 7.5e-6, 9.1e-6, 1.12e-5]),
            np.array([0.052, 0.072, 0.086, 0.1, 0.115, 0.136, 0.166, 0.21]),
            np.array([0.32, 0.54, 0.77, 1.0]),
        )
    if regime == "strong":
        return (
            np.log10([1.65e-5, 2.08e-5, 2.63e-5, 3.30e-5, 4.15e-5, 5.2e-5]),
            np.array([0.26, 0.35, 0.415, 0.465, 0.515, 0.57, 0.64, 0.73, 0.87]),
            np.array([0.065, 0.15, 0.295, 0.47]),
        )
    raise ValueError(f"unknown regime {regime}")


def seed_for(root, label, index):
    value = f"{root}:{label}:{index}".encode()
    return int.from_bytes(hashlib.blake2b(value, digest_size=8).digest(), "little")


def code_hash():
    return hashlib.sha256((HERE / "moran.cpp").read_bytes()).hexdigest()


def save_csv(path, rows):
    rows = list(rows)
    columns = list(dict.fromkeys(k for r in rows for k in r))
    with tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, newline="", delete=False
    ) as f:
        out = csv.DictWriter(f, fieldnames=columns, lineterminator="\n")
        out.writeheader()
        out.writerows(rows)
    path.with_name(f.name.rsplit("/", 1)[-1]).replace(path)


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _save_npz(path, **kwargs):
    with tempfile.NamedTemporaryFile(mode="wb", dir=path.parent, delete=False) as f:
        np.savez_compressed(f, **kwargs)
    path.with_name(f.name.rsplit("/", 1)[-1]).replace(path)


def from_x(model, x):
    return replace(model, mutation=10 ** x[0], advantage=x[1], advantageous=x[2])


def to_x(model):
    return np.array([math.log10(model.mutation), model.advantage, model.advantageous])


def _metadata(regime, model, ages, axes, cap, step):
    return json.dumps(
        {
            "schema": GRID,
            "forward": SPLIT,
            "code": code_hash(),
            "build": " ".join(BUILD_FLAGS),
            "case": regime,
            "copies": model.copies,
            "half_days": model.half_days,
            "ages": np.asarray(ages).tolist(),
            "axes": [x.tolist() for x in axes],
            "cap": cap,
            "step_days": step,
        },
        sort_keys=True,
    )


class Grid:
    def __init__(self, model, days, axes, values, cap, step, meta):
        self.model = model
        self.days = days
        self.axes = axes
        self.cap = cap
        self.step = step
        self.meta = meta
        self.values = values
        self.bounds = Bounds([a[0] for a in axes], [a[-1] for a in axes])
        self.spline = RegularGridInterpolator(
            axes, values[..., :4], method="cubic", solver=spsolve
        )

    def predict(self, x):
        return self.spline(x)[0]

    def direct(self, x):
        r = forward(from_x(self.model, x), self.days, cap=self.cap, step=self.step)
        return tails(r[:, :5]), r[:, 5:].max(axis=0)


def prepare_grid(
    regime,
    cache,
    days=AGES,
    threads=1,
    *,
    copies=1000,
    half_days=10.0,
    axes=None,
    cap=None,
    step=STEP,
):
    """Build or resume a grid. A complete matching grid is portable across paths."""
    model = replace(REGIMES[regime], copies=copies, half_days=half_days)
    cap = min(copies, CAPS[regime] if cap is None else cap)
    axes = axes or axes_for(regime)
    metadata = _metadata(regime, model, days, axes, cap, step)
    key = hashlib.sha256(metadata.encode()).hexdigest()[:16]
    filename = f"{regime}-{key}.npz"
    local = cache / filename
    expected = tuple(len(a) for a in axes) + (len(days), 6)
    cache.mkdir(parents=True, exist_ok=True)
    for candidate in (local, STORED / "grids" / filename):
        if candidate.is_file():
            with np.load(candidate, allow_pickle=False) as z:
                if str(z["metadata"]) != metadata or z["values"].shape != expected:
                    raise ValueError(
                        f"corrupt or mismatched forward grid {candidate.name}"
                    )
                return Grid(
                    model, days, axes, z["values"], cap, step, metadata
                ), candidate

    library()  # Build once before launching parallel solver calls.
    work = cache / (filename + ".points")
    work.mkdir(exist_ok=True)
    coords = list(np.ndindex(expected[:3]))
    values = np.empty(expected)

    def solve(index):
        path = work / ("-".join(map(str, index)) + ".npz")
        if path.is_file():
            with np.load(path, allow_pickle=False) as z:
                return index, z["value"]
        x = np.array([axis[j] for axis, j in zip(axes, index, strict=True)])
        data = forward(from_x(model, x), days, cap=cap, step=step)
        if data[:, 6].max() > 1e-8 or data[:, 5].min() < -1e-8:
            raise RuntimeError(f"{regime} node {index} has invalid forward masses")
        value = np.column_stack([tails(data[:, :5]), data[:, 5:]])
        _save_npz(path, value=value)
        return index, value

    with ThreadPoolExecutor(max_workers=threads) as pool:
        for i, (index, value) in enumerate(pool.map(solve, coords), 1):
            values[index] = value
            if i == len(coords) or i % 24 == 0:
                print(f"{regime}: {i}/{len(coords)} forward nodes", flush=True)
    _save_npz(local, metadata=np.array(metadata), values=values)
    shutil.rmtree(work)
    return Grid(model, days, axes, values, cap, step, metadata), local


def fit(grid, counts, seed, budget):
    from fcmaes.decpp import minimize

    if budget < 32:
        raise ValueError("the budget must allow a population of 32")
    counts = np.asarray(counts)
    if (
        counts.shape != (len(grid.days), 5)
        or (counts < 0).any()
        or (counts.sum(axis=1) <= 0).any()
    ):
        raise ValueError("one nonnegative, nonempty five-bin row is required per age")
    empirical = tails(counts) / counts.sum(axis=1, keepdims=True)

    def objective(x):
        residual = grid.predict(x) - empirical
        return float(np.sum(residual**2))

    result = minimize(
        objective,
        bounds=grid.bounds,
        popsize=32,
        max_evaluations=budget,
        rg=np.random.default_rng(seed),
        workers=1,
    )
    if result.x is None or result.status < 0 or not np.isfinite(result.fun):
        raise RuntimeError("fcmaes did not return a fit")
    model = from_x(grid.model, result.x)
    return {
        **{name: float(getattr(model, name)) for name in PARAMETERS},
        "ssq": float(objective(result.x)),
        "at_reference_ssq": objective(to_x(grid.model)),
        "evaluations": int(result.nfev),
        "optimizer_status": int(result.status),
        "optimizer_seed": seed,
        "fraction_upper_bound": bool(float(result.x[2]) >= grid.bounds.ub[2] - 1e-4),
    }


def fit_summary(regime, rows, bounds):
    """Mean, CV, linear-interpolation percentiles and fits near the box edge."""
    truth = REGIMES[regime]
    out = []
    for k, name in enumerate(PARAMETERS):
        data = np.array([float(row[name]) for row in rows])
        lo, hi = np.quantile(data, [0.025, 0.975])
        actual = getattr(truth, name)
        x = np.log10(data) if name == "mutation" else data
        span = bounds.ub[k] - bounds.lb[k]
        near = (x - bounds.lb[k] <= 0.01 * span) | (bounds.ub[k] - x <= 0.01 * span)
        box = (bounds.lb[k], bounds.ub[k])
        if name == "mutation":
            box = (10 ** box[0], 10 ** box[1])
        out.append(
            {
                "regime": regime,
                "parameter": name,
                "cohorts": len(data),
                "truth": actual,
                "mean": float(data.mean()),
                "median": float(np.median(data)),
                "cv_percent": float(100 * data.std(ddof=1) / data.mean())
                if len(data) > 1
                else None,
                "p2_5": float(lo),
                "p97_5": float(hi),
                "truth_in_interval": bool(lo <= actual <= hi),
                "box_low": float(box[0]),
                "box_high": float(box[1]),
                "fits_within_1pct_of_box_edge": int(near.sum()),
            }
        )
    return out


def diagnostics(grid, regime, fits, snapshots):
    # Off-grid tests are fixed independently of the selected observations.
    checks = {
        "baseline": (
            (5.55e-6, 0.093, 0.89),
            (7.05e-6, 0.107, 0.62),
            (8.8e-6, 0.145, 0.92),
        ),
        "strong": (
            (2.35e-5, 0.472, 0.265),
            (2.98e-5, 0.537, 0.214),
            (3.63e-5, 0.616, 0.32),
        ),
    }
    sequences = [
        ("off-grid", ix, to_x(with_parameters(grid.model, values)), None)
        for ix, values in enumerate(checks[regime], 1)
    ]
    # The whole search box, not only the neighborhood of the truth.
    box = qmc.scale(
        qmc.Sobol(3, scramble=True, seed=26).random_base2(4),
        grid.bounds.lb,
        grid.bounds.ub,
    )
    sequences += [("box", ix, point, None) for ix, point in enumerate(box, 1)]
    ids = sorted({0, len(fits) // 2, len(fits) - 1})
    for ix in ids:
        row = fits[ix]
        point = to_x(with_parameters(grid.model, [row[name] for name in PARAMETERS]))
        data = snapshots[ix]
        sequences.append(
            ("fitted", ix, point, tails(data) / data.sum(axis=1, keepdims=True))
        )
    out = []
    for label, ix, point, observation in sequences:
        expected, (leak, negative) = grid.direct(point)
        pred = grid.predict(point)
        gap = float(np.max(np.abs(pred - expected)))
        out.append(
            {
                "regime": regime,
                "kind": label,
                "index": ix,
                "max_curve_difference": gap,
                "signed_leakage": float(leak),
                "negative_mass": float(negative),
                "ssq_spline": None
                if observation is None
                else float(np.sum((pred - observation) ** 2)),
                "ssq_forward": None
                if observation is None
                else float(np.sum((expected - observation) ** 2)),
                "accepted": gap <= TOLERANCE,
            }
        )
    return out


def with_parameters(model, values):
    return replace(model, **dict(zip(PARAMETERS, values, strict=True)))


def fraction_profile(grid, counts, best, regime):
    """Refit p and s on the first cohort along a fixed advantageous fraction."""
    empirical = tails(counts) / counts.sum(axis=1, keepdims=True)
    lo, hi = grid.bounds.lb, grid.bounds.ub
    initial = to_x(with_parameters(grid.model, [best[x] for x in PARAMETERS]))[:2]
    fractions = np.unique(
        np.r_[
            np.linspace(lo[2], hi[2], 25), grid.model.advantageous, best["advantageous"]
        ]
    )
    for fraction in fractions:

        def residual(x, fraction=fraction):
            return (grid.predict([*x, fraction]) - empirical).ravel()

        result = least_squares(
            residual,
            initial,
            bounds=(lo[:2], hi[:2]),
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        yield {
            "regime": regime,
            "fixed_advantageous_fraction": float(fraction),
            "mutation": 10 ** float(result.x[0]),
            "advantage": float(result.x[1]),
            "profile_ssq": float(result.fun @ result.fun),
            "cohort_best_ssq": best["ssq"],
        }


def load_counts(path):
    rows = read_csv(path)
    try:
        ages = np.array([float(r["age_days"]) for r in rows])
        counts = np.array([[int(r[f"bin_{i}"]) for i in range(5)] for r in rows])
    except (KeyError, ValueError) as exc:
        raise ValueError("expected columns age_days and bin_0 through bin_4") from exc
    if (
        not len(ages)
        or not np.isfinite(ages).all()
        or ages[0] < 0
        or ages[-1] > 20000
        or np.any(np.diff(ages) <= 0)
        or np.any(counts < 0)
        or np.any(counts.sum(axis=1) == 0)
    ):
        raise ValueError("ages must increase and each age needs nonnegative bin counts")
    return ages, counts


def campaign(out, cache, regimes, reps, cells, budget, seed, threads):
    """Regenerate cohorts and fit each using its own deterministic seed."""
    out.mkdir(parents=True, exist_ok=True)
    fits, summaries, checks, profiles, timings = [], [], [], [], []
    manifest = {
        "root_seed": seed,
        "cases": regimes,
        "time_days": AGES.tolist(),
        "bootstrap_replicates": reps,
        "cells_per_replicate": cells,
        "max_evaluations": budget,
        "cohort_sampler": SAMPLER,
        "fit_sampler": "fcmaes.decpp",
        "model": SPLIT,
        "proxy": GRID,
        "moran_cpp_sha256": code_hash(),
        "inputs": {r: asdict(REGIMES[r]) for r in regimes},
        "grids": {},
    }
    work = out / "cohorts"
    work.mkdir(exist_ok=True)
    for regime in regimes:
        start = time.perf_counter()
        grid, path = prepare_grid(regime, cache, threads=threads)
        manifest["grids"][regime] = {
            **json.loads(grid.meta),
            "file": path.name,
            "max_node_leakage": float(grid.values[..., 4].max()),
            "max_node_negative_mass": float(grid.values[..., 5].max()),
        }
        grid_time = time.perf_counter() - start
        snapshots, current = [], []
        sim_time = fit_time = 0.0
        for i in range(reps):
            cell_seed = seed_for(seed, regime + "-data", i)
            search_seed = seed_for(seed, regime + "-de", i)
            file = work / f"{regime}-{cells}-{cell_seed:016x}-{code_hash()[:12]}.npz"
            tic = time.perf_counter()
            if file.is_file():
                with np.load(file, allow_pickle=False) as z:
                    count = z["counts"]
            else:
                count = cohort(grid.model, cells, cell_seed, threads=threads)[0]
                _save_npz(file, counts=count)
            snapshots.append(count)
            sim_time += time.perf_counter() - tic
            tic = time.perf_counter()
            row = {
                "regime": regime,
                "cohort": i,
                "cells": cells,
                "cohort_seed": cell_seed,
                **fit(grid, count, search_seed, budget),
            }
            current.append(row)
            fits.append(row)
            fit_time += time.perf_counter() - tic
            save_csv(out / "fits.csv", fits)
            if i + 1 == reps or (i + 1) % 10 == 0:
                print(f"{regime}: {i + 1}/{reps} fitted cohorts", flush=True)
        summaries.extend(fit_summary(regime, current, grid.bounds))
        checks.extend(diagnostics(grid, regime, current, snapshots))
        profiles.extend(fraction_profile(grid, snapshots[0], current[0], regime))
        save_csv(out / "fraction-profile.csv", profiles)
        first = (
            {"age_days": age, **{f"bin_{i}": row[i] for i in range(5)}}
            for age, row in zip(AGES, snapshots[0], strict=True)
        )
        save_csv(out / f"snapshots-{regime}.csv", first)
        save_csv(out / "bootstrap.csv", summaries)
        save_csv(out / "grid-checks.csv", checks)
        timings.append(
            {
                "regime": regime,
                "prepare_seconds": grid_time,
                "cohort_seconds": sim_time,
                "fit_seconds": fit_time,
            }
        )
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "timings.json").write_text(json.dumps(timings, indent=2) + "\n")
    return manifest
