"""Forward sensitivities on a shared prepared finite state projection."""

from dataclasses import dataclass
from math import fsum
from operator import index

import numpy as np

from . import _native
from .likelihood import measured_counts
from .projection import _native_error, _vector


@dataclass(frozen=True)
class ProjectionSensitivities:
    """Snapshots include a final bin for mass and derivatives in the FSP sink."""

    times: tuple
    species: tuple
    states: np.ndarray
    probabilities: np.ndarray
    derivatives: np.ndarray
    work: int

    def observed(self, axes):
        """Aggregate hidden states before computing any score."""
        values = measured_counts(self.states, self.species, axes)
        bins, inverse = np.unique(values, return_inverse=True, axis=0)
        channels = np.concatenate(
            (self.probabilities[:, None], self.derivatives), axis=1
        )
        observed = np.zeros((*channels.shape[:2], len(bins)))
        for i, frame in enumerate(channels[:, :, :-1]):
            for j, distribution in enumerate(frame):
                observed[i, j] = np.bincount(inverse, distribution, minlength=len(bins))
        return bins, observed[:, 0], observed[:, 1:]


def forward_sensitivities(
    graph,
    segments,
    times,
    parameters,
    *,
    initial=None,
    initial_derivatives=None,
    logarithmic=False,
    difference=1e-4,
    tolerance=1e-6,
    max_work=300_000,
    max_factor_bytes=512 << 20,
):
    """Solve p'=Qp and s'=Qs+(dQ/dtheta)p over piecewise-constant bindings.

    Each segment is (end_time, bindings), measured from time zero. Requested
    times increase strictly. Differentiated names must be bound in every
    segment. Use centered differences of the generator on the same graph.
    With logarithmic=True the directions are natural log parameters.
    Initial derivatives have one row per parameter and one column per state,
    in graph order, without a sink column. The default initial law is fixed.
    Count interventions, stop states and continuously varying signals are
    outside this interface. Work counts augmented solves, not individual
    block substitutions. The caller controls the truncation.
    """
    names = tuple(parameters)
    times = tuple(map(float, times))
    segments = tuple((float(end), dict(binding)) for end, binding in segments)
    knots = tuple(end for end, _ in segments)
    if (
        not names
        or len(set(names)) != len(names)
        or len(names) > 16
        or not times
        or not knots
        or not np.isfinite(times + knots).all()
        or times[0] < 0
        or times[-1] > knots[-1]
        or knots[0] <= 0
        or any(a >= b for a, b in zip(times, times[1:], strict=False))
        or any(a >= b for a, b in zip(knots, knots[1:], strict=False))
        or not 0 < difference < 0.1
        or not 0 < tolerance < 1
        or index(max_work) < 1
        or index(max_factor_bytes) < 1
    ):
        raise ValueError(
            "invalid sensitivity times, segments, parameters or tolerances"
        )
    n = len(graph.states)
    if len(times) * (n + 1) * (len(names) + 1) > 8_000_000:
        raise ValueError("sensitivity snapshots exceed the output limit")
    state = np.zeros((len(names) + 1, n + 1), dtype="<f8")
    for label, probability in graph._law(initial).items():
        position = graph._index(label)
        if position is None:
            raise ValueError("initial law outside the sensitivity graph")
        state[0, position] = probability
    if initial_derivatives is not None:
        directions = np.asarray(initial_derivatives, dtype=float)
        if (
            directions.shape != (len(names), n)
            or not np.isfinite(directions).all()
            or np.max(np.abs(directions.sum(axis=1))) > 1e-9
        ):
            raise ValueError("invalid initial derivatives")
        state[1:, :-1] = directions
    if abs(fsum(state[0]) - 1) > 1e-12:
        raise ValueError("sensitivity initial probability must sum to one")
    snapshots = [state.copy()] if times[0] == 0 else []
    start = 0.0
    work = 0
    for end, binding in segments:
        generator, _, _ = graph._generator(binding, {})
        pairs = []
        for name in names:
            value = binding.get(name)
            if value is None or not np.isfinite(value):
                raise ValueError(f"missing finite bound parameter {name}")
            if logarithmic:
                if value <= 0:
                    raise ValueError(
                        "logarithmic sensitivities need positive parameters"
                    )
                lower, upper = value * np.exp(np.array([-difference, difference]))
                width = 2 * difference
            else:
                step = difference * max(abs(value), 1e-3)
                lower, upper = value - step, value + step
                width = upper - lower
            a, _, _ = graph._generator(binding | {name: float(upper)}, {})
            b, _, _ = graph._generator(binding | {name: float(lower)}, {})
            pairs.append((a, b, float(width)))
        if graph._actions:
            raise ValueError(
                "count interventions are unavailable in sensitivity segments"
            )
        interior = [t for t in times if start < t < end]
        try:
            reports, final, steps = generator.tangents(
                state.tobytes(),
                pairs,
                start,
                end,
                interior,
                tolerance,
                max_work - work,
                max_factor_bytes,
            )
        except _native.NativeError as error:
            raise _native_error(error) from None
        snapshots.extend(_vector(reports).reshape(len(interior), *state.shape))
        state = _vector(final).reshape(state.shape)
        if end in times:
            snapshots.append(state)
        work += steps
        start = end
    snapshots = np.array(snapshots)
    if (
        not np.isfinite(snapshots).all()
        or np.max(np.abs(snapshots[:, 0].sum(axis=1) - 1.0)) > 1e-7
        or np.max(np.abs(snapshots[:, 1:].sum(axis=2))) > 1e-7
        or np.min(snapshots[:, 0]) < -max(10 * tolerance * 1e-4, 1e-11)
    ):
        raise ValueError("sensitivity solver lost probability or derivative mass")
    return ProjectionSensitivities(
        times, graph.species, graph.states, snapshots[:, 0], snapshots[:, 1:], work
    )


def fisher_information(probabilities, derivatives, cells, *, probability_floor=0.0):
    """Independent snapshot Fisher information after observation aggregation.

    Use a positive floor only to match a documented benchmark convention.
    Return the retained information and excluded probability at each time.
    """
    p = np.asarray(probabilities, dtype=float)
    s = np.asarray(derivatives, dtype=float)
    if (
        p.ndim != 2
        or s.ndim != 3
        or (s.shape[0], s.shape[2]) != p.shape
        or not np.isfinite(p).all()
        or not np.isfinite(s).all()
        or not 0 <= probability_floor < 1
        or np.min(p) < -1e-8
    ):
        raise ValueError("invalid Fisher probabilities or derivatives")
    cells = np.broadcast_to(cells, (len(p),))
    if not np.isfinite(cells).all() or np.any(cells < 0):
        raise ValueError("invalid cell counts")
    keep = (p > 0) & (p >= probability_floor)
    root = np.sqrt(np.divide(cells[:, None], p, out=np.zeros_like(p), where=keep))
    weighted = s * root[:, None]
    info = np.einsum("tip,tjp->ij", weighted, weighted)
    if not np.isfinite(info).all():
        raise ValueError("nonfinite Fisher information")
    return info, np.sum(np.where(keep, 0.0, np.maximum(p, 0.0)), axis=1)
