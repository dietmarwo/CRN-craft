"""Forward sensitivities and Fisher information on a prepared projection.

:func:`forward_sensitivities` integrates the law ``p' = Q p`` of a prepared
finite state projection (:func:`crn_craft.prepare_fsp`) together with its
derivatives ``s_j' = Q s_j + (dQ/dtheta_j) p`` over piecewise-constant
parameter segments. :func:`fisher_information` turns observed probabilities
and derivatives into the Fisher information of independent single-cell
snapshots. Aggregate hidden states into measured bins first
(:meth:`ProjectionSensitivities.observed`).
"""

from dataclasses import dataclass
from math import fsum
from operator import index

import numpy as np

from . import _native
from .likelihood import measured_counts
from .projection import _native_error, _vector


@dataclass(frozen=True)
class ProjectionSensitivities:
    """Laws and their parameter derivatives at the requested times.

    ``probabilities`` has shape ``(times, states + 1)`` and ``derivatives``
    ``(times, parameters, states + 1)``, in the order of ``states``; the last
    column is the absorbing sink. Each row of ``probabilities`` sums to one
    and each derivative row to zero, sink included. ``work`` counts the
    augmented implicit solves.
    """

    times: tuple
    species: tuple
    states: np.ndarray
    probabilities: np.ndarray
    derivatives: np.ndarray
    work: int

    def observed(self, axes):
        """Sum states and derivatives into the bins of the measured counts.

        ``axes`` names the measured species. Returns the distinct measured
        count tuples, their probabilities ``(times, bins)`` and derivatives
        ``(times, parameters, bins)``; the sink is not a bin.
        """
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

    ``segments`` is a sequence of ``(end_time, bindings)``, measured from
    time zero, with strictly increasing ends; each binding maps parameter
    names to values for that segment and must bind every name in
    ``parameters`` (at most 16 names). ``times`` increase strictly, lie in
    ``[0, last end]``; time zero returns the initial state. The law and its
    derivatives pass continuously across segment ends, where the integrator
    restarts.

    ``dQ/dtheta`` is the centered difference of the generator at
    ``theta (1 +- difference)`` (or ``theta exp(+-difference)`` with
    ``logarithmic=True``, which differentiates in natural-log coordinates)
    on the same state graph. It is exact up to roundoff for rates linear in
    the parameter and second-order accurate otherwise, so parameters inside
    any rate expression, Hill functions included, can be differentiated;
    the time integration is the implicit variable-order BDF of the prepared
    projection with the derivatives inside its error norm (relative
    ``tolerance``, absolute ``1e-4 tolerance``), and one sparse
    factorization serves the law and every direction.

    ``initial`` is a law as in :meth:`PreparedFsp.solve`; by default the
    initial law given to :func:`crn_craft.prepare_fsp`. It does not depend
    on the parameters unless ``initial_derivatives`` (one row per parameter,
    one column per state in graph order, no sink column, rows summing to
    zero) says so. Scheduled count interventions, stop states and
    continuously varying signals are outside this interface.

    Raises ``ValueError`` for invalid arguments, a lost probability or
    derivative mass above ``1e-7`` or a negative probability beyond the
    integration tolerance, and the errors of the prepared projection for
    budgets (``max_work`` augmented solves, ``max_factor_bytes``) and
    interruption. The caller chooses the truncation: read the sink.
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
    law = graph._initial if initial is None else graph._law(initial)
    for label, probability in law.items():
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
    """Fisher information of independent snapshots of measured bins.

    ``probabilities`` ``(times, bins)`` and ``derivatives``
    ``(times, parameters, bins)`` describe MEASURED bins (aggregate hidden
    states first). ``cells`` is the number of cells measured at each time, a
    number or one per time. Returns ``(matrix, excluded)``: the symmetric
    positive semidefinite matrix ``sum_t n_t sum_b g g^T / p`` over bins with
    ``p > 0`` and ``p >= probability_floor``, and per time the probability of
    the excluded bins. A bin with ``p = 0`` is unreachable or has underflowed
    and carries no information; a positive floor reproduces a published
    convention that drops rare bins and lowers the information.
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
        or np.min(p, initial=0.0) < -1e-8
    ):
        raise ValueError("invalid Fisher probabilities or derivatives")
    try:
        cells = np.broadcast_to(np.asarray(cells, dtype=float), (len(p),))
    except ValueError:
        raise ValueError("invalid cell counts") from None
    if not np.isfinite(cells).all() or np.any(cells < 0):
        raise ValueError("invalid cell counts")
    keep = (p > 0) & (p >= probability_floor)
    root = np.sqrt(np.divide(cells[:, None], p, out=np.zeros_like(p), where=keep))
    weighted = s * root[:, None]
    info = np.einsum("tip,tjp->ij", weighted, weighted)
    info = 0.5 * (info + info.T)
    if not np.isfinite(info).all():
        raise ValueError("nonfinite Fisher information")
    return info, np.sum(np.where(keep, 0.0, np.maximum(p, 0.0)), axis=1)
