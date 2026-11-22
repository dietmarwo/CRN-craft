"""Sparse stationary laws on finite and progressively enlarged projections."""

from dataclasses import dataclass
from operator import index

import numpy as np

from . import _native
from .projection import _MAX_FACTOR_BYTES, _native_error, _vector, prepare_fsp


@dataclass(frozen=True)
class StationaryResult:
    """Stationary probabilities and empirical projection diagnostics.

    ``method`` is ``"finite-closed-class"`` for an exact finite generator
    and ``"stationary-fsp-return-state"`` for outward-rate redirection. It
    distinguishes this result from the transient absorbing-sink projection.
    """

    species: tuple
    states: np.ndarray
    probabilities: np.ndarray
    residual_l1: float
    outward_flux: float
    return_state: tuple | None
    marginal_cdf_change: float | None
    history: tuple
    converged: bool
    method: str

    def marginal(self, species):
        """Return the projected probability of each count of one species."""
        return np.bincount(
            self.states[:, self.species.index(species)], self.probabilities
        )


def _solve(graph, parameters, redirect, selected, max_factor_bytes):
    generator, _, stopped = graph._generator(parameters, {})
    if np.any(stopped) and redirect is not None:
        raise ValueError(
            "stationary FSP requires an ergodic model without a reachable stop"
        )
    escapes = _vector(generator.escapes())
    anchor = None
    if redirect is not None:
        anchor = graph._index(redirect)
        if anchor is None:
            raise ValueError("return state is outside the reachable projection")
    elif np.any(escapes):
        raise ValueError(
            "the finite generator has escaping transitions; supply a return state"
        )
    target = None
    if selected is not None:
        target = graph._index(selected)
        if target is None:
            raise ValueError("selected class state is outside the projection")
    try:
        data, residual, flux, _ = generator.stationary(anchor, target, max_factor_bytes)
    except _native.NativeError as error:
        raise _native_error(error) from None
    return _vector(data), residual, flux, escapes


def _cdf_difference(previous, current):
    gap = 0.0
    for p, q in zip(previous, current, strict=True):
        width = max(len(p), len(q))
        diff = np.pad(p, (0, width - len(p))) - np.pad(q, (0, width - len(q)))
        gap = max(gap, float(np.max(np.abs(np.cumsum(diff)))))
    return gap


def stationary_fsp(
    model,
    bounds,
    *,
    return_state=None,
    class_state=None,
    parameters=None,
    expansions=0,
    tolerance=1e-6,
    flux_tolerance=1e-6,
    max_states=100_000,
    max_factor_bytes=_MAX_FACTOR_BYTES,
):
    """Solve a finite stationary law or an outward-flux-redirected FSP.

    With a return state, expand faces according to boundary flux exposure
    and compare marginal CDFs across successive bounds. The resulting
    convergence indication is empirical for an infinite network. A finite
    generator with several closed classes needs an explicit class state.
    The native engine builds the generator, finds the closed communicating
    classes and solves the anchored balance equations of one class with a
    sparse LU, refused when its bound exceeds ``max_factor_bytes``. Without
    a return state, every transition must stay inside the bounds.
    """
    expansions = index(expansions)
    max_factor_bytes = index(max_factor_bytes)
    if (
        not 0 <= expansions <= 20
        or not 0 < tolerance < 1
        or not np.isfinite(flux_tolerance)
        or flux_tolerance < 0
        or (class_state is not None and return_state is not None)
        or (expansions and return_state is None)
        or max_factor_bytes < 0
    ):
        raise ValueError("invalid stationary approximation controls")
    if return_state is not None:
        return_state = tuple(index(v) for v in return_state)
    if class_state is not None:
        class_state = tuple(index(v) for v in class_state)
    bounds = dict(bounds)
    history, previous = [], None
    for step in range(expansions + 1):
        graph = prepare_fsp(model, bounds, max_states=max_states)
        p, residual, flux, escapes = _solve(
            graph, parameters or {}, return_state, class_state, max_factor_bytes
        )
        marginals = tuple(
            np.bincount(graph.states[:, i], p) for i in range(len(graph.species))
        )
        delta = None if previous is None else _cdf_difference(previous, marginals)
        history.append(
            {
                "bounds": dict(bounds),
                "states": len(p),
                "residual_l1": residual,
                "outward_flux": flux,
                "marginal_cdf_change": delta,
            }
        )
        closed = return_state is None or flux == 0.0
        converged = closed or (
            delta is not None and delta <= tolerance and flux <= flux_tolerance
        )
        if converged or step == expansions:
            p.setflags(write=False)
            return StationaryResult(
                graph.species,
                graph.states,
                p,
                residual,
                flux,
                return_state,
                delta,
                tuple(history),
                converged,
                "finite-closed-class"
                if return_state is None
                else "stationary-fsp-return-state",
            )
        previous = marginals
        out = escapes * p
        face_fluxes = []
        for i, label in enumerate(graph.species):
            # A jump may cross a face from within the largest positive jump.
            reach = max(
                (
                    change
                    for c in graph._channels
                    for j, change in c["changes"]
                    if i == j and change > 0
                ),
                default=1,
            )
            face_fluxes.append(
                float(out[graph.states[:, i] > bounds[label] - reach].sum())
            )
        largest = max(face_fluxes)
        for face, label in zip(face_fluxes, graph.species, strict=True):
            if face >= 0.1 * largest:
                bounds[label] += max(2, (bounds[label] + 3) // 4)
    raise RuntimeError("unreachable stationary iteration")
