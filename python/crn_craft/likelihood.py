"""Observation maps and grouped likelihoods for stochastic count models."""

from dataclasses import dataclass
from math import fsum, isfinite

import numpy as np


def measured_counts(states, species, observable):
    """Map projected states to one count, a species sum, or a joint readout.

    A string names one species. Each entry in a tuple is one measured axis.
    An axis can be a string or a list of species to sum.
    """
    labels = tuple(species)
    axes = (observable,) if isinstance(observable, str) else tuple(observable)
    if not axes:
        raise ValueError("at least one measured count is required")
    columns = []
    for axis in axes:
        names = [axis] if isinstance(axis, str) else tuple(axis)
        if not names:
            raise ValueError("empty measured axis")
        columns.append(
            np.asarray(states)[:, [labels.index(key) for key in names]].sum(axis=1)
        )
    return np.column_stack(columns)


def count_probability(states, probabilities, values, *, detection=1.0):
    """Compute the probability of each requested count or joint count."""
    states = np.asarray(states)
    bins = np.asarray(values)
    if states.ndim == 1:
        states = states[:, None]
    if bins.ndim == 1:
        bins = bins[:, None]
    if (
        states.ndim != 2
        or bins.ndim != 2
        or states.shape[1] != bins.shape[1]
        or np.any(states < 0)
        or np.any(bins < 0)
        or np.any(states != np.floor(states))
        or np.any(bins != np.floor(bins))
    ):
        raise ValueError("states and bins must be matching nonnegative count vectors")
    p = _probabilities(probabilities, len(states))
    efficiency = np.broadcast_to(detection, (states.shape[1],))
    if np.any((efficiency < 0) | (efficiency > 1)) or not np.isfinite(efficiency).all():
        raise ValueError("detection probabilities must be in [0, 1]")
    if np.all(efficiency == 1):
        result = {}
        for key, weight in zip(map(tuple, states), p, strict=True):
            result[key] = result.get(key, 0.0) + weight
        return np.array([result.get(tuple(key), 0.0) for key in bins])

    from scipy.stats import binom

    return np.array(
        [(binom.pmf(key, states, efficiency).prod(axis=1)) @ p for key in bins]
    )


def binned_probability(
    counts, probabilities, edges, *, gain=None, background=0.0, sigma=None
):
    """Map one count to intervals or to Gaussian fluorescence bins.

    Edges are increasing and bins are left-closed, right-open. Infinite
    edges represent underflow or overflow. Fluorescence has mean
    background + gain * count and standard deviation sigma.
    """
    counts, edges = np.asarray(counts), np.asarray(edges, dtype=float)
    if (
        counts.ndim != 1
        or not np.isfinite(counts).all()
        or edges.ndim != 1
        or len(edges) < 2
        or np.any(np.isnan(edges))
        or np.any(edges[1:] <= edges[:-1])
    ):
        raise ValueError("finite counts and strictly increasing bin edges are required")
    p = _probabilities(probabilities, len(counts))
    if gain is None and sigma is None:
        i = np.searchsorted(edges, counts, side="right") - 1
        selected = (i >= 0) & (i < len(edges) - 1)
        return np.bincount(i[selected], p[selected], minlength=len(edges) - 1)
    if (
        gain is None
        or sigma is None
        or not all(map(isfinite, (gain, sigma, background)))
        or gain < 0
        or sigma <= 0
    ):
        raise ValueError("fluorescence needs a nonnegative gain and positive noise")
    from scipy.special import ndtr

    mean = background + gain * counts
    result = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        left, right = (lo - mean) / sigma, (hi - mean) / sigma
        weights = np.where(
            left > 0, ndtr(-left) - ndtr(-right), ndtr(right) - ndtr(left)
        )
        result.append(float(p @ weights))
    return np.asarray(result)


def _probabilities(probabilities, size):
    p = np.asarray(probabilities, dtype=float)
    if p.shape != (size,) or not np.isfinite(p).all() or np.any(p < -1e-10):
        raise ValueError("invalid probability array")
    if p.sum() > 1 + 1e-5:
        raise ValueError("probability sum exceeds one")
    return np.maximum(p, 0.0)


@dataclass(frozen=True)
class HistogramGroup:
    """Counts and predicted observed bin probabilities for one sample group."""

    condition: str
    replicate: str
    time: float
    counts: object
    probabilities: object
    missing_mass: float = 0.0
    time_unit: str = ""
    readout: str = ""


@dataclass(frozen=True)
class HistogramScore:
    """Negative log likelihood and a bound allowing for uncaptured mass."""

    nats_per_cell: float
    optimistic_nats_per_cell: float
    cells: int
    groups: tuple
    impossible_bins: int


def histogram_score(groups):
    """Score grouped histograms and expose the cost uncertainty due to leakage."""
    rows, cells, impossible = [], 0, 0
    for group in groups:
        counts = np.asarray(group.counts)
        if (
            counts.ndim != 1
            or not np.isfinite(counts).all()
            or (np.any(counts < 0) or np.any(counts != np.floor(counts)))
        ):
            raise ValueError("histogram counts must be finite nonnegative integers")
        if not isfinite(group.time) or not 0 <= group.missing_mass <= 1:
            raise ValueError("invalid observation time or missing probability")
        q = _probabilities(group.probabilities, len(counts))
        use = counts > 0
        n = sum(map(int, counts))
        floor = np.minimum(1.0, q[use])
        ceiling = np.minimum(1.0, floor + group.missing_mass)
        with np.errstate(divide="ignore"):
            cost = float(-counts[use] @ np.log(floor))
            optimistic = float(-counts[use] @ np.log(ceiling))
        impossible += int(np.count_nonzero(ceiling == 0))
        rows.append((group.condition, group.replicate, group.time, n, cost, optimistic))
        cells += n
    if not cells:
        raise ValueError("the score needs at least one observed cell")
    return HistogramScore(
        fsum(row[-2] for row in rows) / cells,
        fsum(row[-1] for row in rows) / cells,
        cells,
        tuple(rows),
        impossible,
    )


def mean_score(observed, predicted, covariance):
    """Return the quadratic mean objective for a group covariance."""
    delta = np.atleast_1d(np.asarray(observed, dtype=float) - predicted)
    v = np.atleast_2d(covariance)
    if v.shape != (len(delta), len(delta)) or not np.allclose(v, v.T):
        raise ValueError("a symmetric covariance must match the observation axes")
    transformed = np.linalg.solve(np.linalg.cholesky(v), delta)
    return float(transformed @ transformed / 2)
