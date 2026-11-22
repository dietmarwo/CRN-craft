"""Reusable finite state projections on the native CRN-craft engine."""

import json
from collections.abc import Mapping
from dataclasses import dataclass
from hashlib import sha256
from math import fsum, isfinite
from operator import index

import numpy as np

from . import _native
from ._common import model_bytes
from .protocol import Protocol
from .rules import _as_model_1_0, _crn_error, _unique_keys

_MAX_COUNT = 2**53
# Counts and increments beyond this behave like any value outside the bounds.
_CLAMP = 2**62
_OPERATIONS = {"set": 0, "change": 1, "partition": 2}


def _require_scipy(feature):
    """Fail with the name of the optional extra when SciPy is missing."""
    try:
        import scipy  # noqa: F401
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            f"{feature} needs SciPy; install it with "
            "pip install 'crn-craft[projection]'",
            name="scipy",
        ) from error


class Signal:
    """A real-time rate input. Tables use linear or right-continuous steps.

    Tables are evaluated natively by the midpoint rule. Callables are a
    Python-only option, called from the native step loop. Provide their
    change points and a reproducible description of their parameters in the
    calling experiment.
    """

    def __init__(self, function, *, change_times=()):
        if not callable(function):
            raise TypeError("a prescribed signal needs a time function")
        self.evaluate = function
        self._table_sha256 = None
        self._table = None
        self.change_times = tuple(float(x) for x in change_times)
        if any(not isfinite(x) for x in self.change_times):
            raise ValueError("signal breakpoints must be finite")

    @classmethod
    def table(cls, times, values, *, interpolation="linear"):
        """Create a table with held endpoint values outside its range."""
        t = np.asarray(tuple(times), dtype=float)
        v = np.asarray(tuple(values), dtype=float)
        if (
            t.ndim != 1
            or not len(t)
            or v.shape != t.shape
            or not np.isfinite(t).all()
            or not np.isfinite(v).all()
            or np.any(np.diff(t) <= 0)
        ):
            raise ValueError("signals need finite values on an increasing time grid")
        if interpolation == "linear":
            signal = cls(lambda at: np.interp(at, t, v))
        elif interpolation == "step":
            signal = cls(
                lambda at: v[max(0, np.searchsorted(t, at, side="right") - 1)],
                change_times=t,
            )
        else:
            raise ValueError("interpolation must be linear or step")
        signal._table = (
            signal.evaluate,
            interpolation == "step",
            t.tolist(),
            v.tolist(),
        )
        signal._table_sha256 = sha256(
            model_bytes(
                {
                    "times": t.tolist(),
                    "values": v.tolist(),
                    "interpolation": interpolation,
                }
            )
        ).hexdigest()
        return signal


@dataclass(frozen=True)
class ProjectionState:
    """A report-time law and its restart bookkeeping."""

    time: float
    probabilities: np.ndarray
    stopped_probabilities: np.ndarray
    physical_leakage: float
    series_tail_bound: float
    mass_defect: float
    minimum_probability: float
    graph_id: str
    actions_applied: int
    binding_id: str


@dataclass(frozen=True)
class ProjectionResult:
    """The law at requested times, with one absorbing truncation sink."""

    species: tuple
    states: np.ndarray
    snapshots: tuple
    work: int
    solver: str
    model_sha256: str
    protocol_sha256: str
    signal_table_sha256: dict

    def marginal(self, species, *, at=-1):
        """Return the count probabilities of one species, with sink excluded."""
        values = self.states[:, self.species.index(species)]
        return np.bincount(values, self.snapshots[at].probabilities)


class _Work:
    def __init__(self, limit):
        self.limit = index(limit)
        self.count = 0
        # The native counter is an int64; larger limits cannot be reached.
        self.native_limit = max(-(2**63), min(self.limit, 2**63 - 1))

    def charge(self, amount=1):
        self.count += amount
        if self.count > self.limit:
            raise RuntimeError(
                f"projection work limit of {self.limit} exceeded; raise max_work"
            )


def _native_error(error, work=None, signals=()):
    """The exception the projection raises for a native refusal.

    Budget, work and integrator failures (``CRN-E070``) are ``RuntimeError``,
    inputs the projection cannot represent (``CRN-E071``) ``ValueError``;
    other diagnostics (rate evaluation) are a ``CrnError``.
    """
    try:
        first = json.loads(str(error))[0]
        code, path, message = first["code"], first["path"], first["message"]
    except (ValueError, LookupError, TypeError):
        return _crn_error(error)
    if code == "CRN-E070":
        if path == "/projection/work" and work is not None:
            message = f"projection work limit of {work.limit} exceeded; raise max_work"
        return RuntimeError(message)
    if code == "CRN-E071":
        if path == "/projection/overflow":
            return OverflowError(message)
        if path.startswith("/signals/"):
            name = signals[int(path.rsplit("/", 1)[1])]
            message = f"{name!r} is not a rule with a count-independent rate"
        return ValueError(message)
    return _crn_error(error)


def _clamp(n):
    return max(-_CLAMP, min(n, _CLAMP))


def _vector(data):
    return np.frombuffer(data, dtype="<f8")


def _sparse(generator, signal=None):
    """One native generator matrix as a SciPy CSC matrix."""
    from scipy import sparse

    dim, indptr, indices, data = generator.csc(signal)
    return sparse.csc_matrix(
        (
            _vector(data),
            np.frombuffer(indices, dtype="<i8"),
            np.frombuffer(indptr, dtype="<i8"),
        ),
        shape=(dim, dim),
    )


def _source(signal):
    """A native table, or the Python callable of the signal."""
    table = signal._table
    if table is not None and table[0] is signal.evaluate:
        return table[1:]
    return signal.evaluate


def _view(document):
    try:
        native = _native.ProjectionView(model_bytes(document))
        return native, json.loads(native.layout())
    except _native.NativeError as error:
        raise _crn_error(error) from None


def _signature(layout):
    return (
        tuple(layout["species"]),
        tuple(
            (c["name"], c["needs"], c["changes"], c["whole"])
            for c in layout["channels"]
        ),
    )


def _action_indices(actions, species, parameters):
    compiled = []
    previous = 0.0
    for action in actions:
        t = action["time"]
        if not isfinite(t) or t < previous:
            raise ValueError("action times must be finite and ordered")
        previous = t
        changes = {}
        for key in ("change", "set", "partition"):
            for name, value in action.get(key, {}).items():
                if name not in species or species[name] in changes:
                    raise ValueError(f"unknown or conflicting intervention on {name!r}")
                if (key == "partition" and not 0 <= value <= 1) or (
                    key == "set" and not 0 <= value <= _MAX_COUNT
                ):
                    raise ValueError(f"invalid intervention on {name!r}")
                changes[species[name]] = (key, value)
        edits = action.get("parameters", {})
        if edits.keys() - parameters.keys() or not all(map(isfinite, edits.values())):
            raise ValueError("invalid kinetic parameter edit")
        compiled.append((t, changes, edits))
    return tuple(compiled)


class PreparedFsp:
    """Fixed projection topology shared by independent numerical solves."""

    def __init__(
        self,
        model,
        bounds,
        *,
        protocol=None,
        initial=None,
        max_states=100_000,
        max_edges=2_000_000,
    ):
        self._document = json.loads(
            _as_model_1_0(model, network=True), object_pairs_hook=_unique_keys
        )
        self._document.setdefault("parameters", {})
        self._native, self._layout = _view(self._document)
        self.species = tuple(self._layout["species"])
        if bounds.keys() != set(self.species):
            raise ValueError("bounds must cover every lowered species")
        self._bounds = tuple(index(bounds[s]) for s in self.species)
        self._max_states = index(max_states)
        self._max_edges = index(max_edges)
        if (
            not 1 <= self._max_states <= 250_000
            or not 1 <= self._max_edges <= 2_000_000
            or any(not 0 <= n <= _MAX_COUNT for n in self._bounds)
        ):
            raise ValueError("projection bounds or budgets out of range")
        if protocol is not None and not isinstance(protocol, Protocol):
            raise TypeError("protocol must be a Protocol")
        self._protocol = (
            protocol.to_dict() if protocol else {"actions": [], "delays": []}
        )
        if self._protocol["delays"]:
            raise ValueError("fixed-delay histories require trajectory simulation")
        self._actions = _action_indices(
            self._protocol["actions"],
            dict(zip(self.species, range(len(self.species)), strict=True)),
            self._document["parameters"],
        )
        self._initial = self._law(initial)
        self._channels = self._layout["channels"]
        interventions = [
            [
                (
                    i,
                    _OPERATIONS[kind],
                    _clamp(value) if kind != "partition" else 0,
                    value if kind == "partition" else 0.0,
                )
                for i, (kind, value) in changes.items()
            ]
            for _, changes, _ in self._actions
            if changes
        ]
        try:
            self._topology = _native.ProjectionTopology(
                self._native,
                list(self._bounds),
                [[_clamp(n) for n in state] for state in sorted(self._initial)],
                interventions,
                self._max_states,
                self._max_edges,
            )
        except _native.NativeError as error:
            raise _native_error(error) from None
        self.states = np.frombuffer(self._topology.states(), dtype="<i8").reshape(
            -1 if self.species else 1, len(self.species)
        )
        state_bytes = self.states.astype("<i8", copy=False).tobytes()
        self._id = sha256(
            model_bytes([self.species, self._protocol, self._document]) + state_bytes
        ).hexdigest()
        self.states.setflags(write=False)

    def _index(self, state):
        """The index of a count tuple, or None outside the projection."""
        if len(state) != len(self.species) or any(
            not 0 <= n <= _MAX_COUNT for n in state
        ):
            return None
        return self._topology.index(list(state))

    def _law(self, initial):
        if initial is None:
            return {tuple(self._layout["initial"]): 1.0}
        if not isinstance(initial, Mapping):
            raise TypeError("the initial law maps count tuples to probabilities")
        law = {tuple(index(n) for n in state): float(p) for state, p in initial.items()}
        if (
            not law
            or any(not isfinite(p) or p < 0 for p in law.values())
            or abs(fsum(law.values()) - 1.0) > 1e-12
        ):
            raise ValueError("initial probabilities must be nonnegative and sum to one")
        return {k: p for k, p in law.items() if p}

    def _binding(self, parameters):
        if not parameters:
            return self._native, self._layout
        eligible = self._document["parameters"]
        if parameters.keys() - eligible.keys() or not all(
            map(isfinite, parameters.values())
        ):
            raise ValueError("invalid projection parameter binding")
        source = {**self._document, "parameters": {**eligible, **parameters}}
        native, layout = _view(source)
        if (
            _signature(layout) != _signature(self._layout)
            or layout["initial"] != self._layout["initial"]
            or any(
                (c["rate"] is None) != (d["rate"] is None)
                for c, d in zip(layout["channels"], self._channels, strict=True)
            )
        ):
            raise ValueError(
                "binding changes topology or initial counts; prepare a new projection"
            )
        return native, layout

    def _generator(self, parameters, signals):
        """The native generator under a binding, the signals in its order,
        and the stop flags."""
        native, layout = self._binding(parameters)
        channels = layout["channels"]
        names = tuple(signals)
        ids = [[i for i, c in enumerate(channels) if c["name"] == n] for n in names]
        try:
            generator = self._topology.generator(native, ids)
        except _native.NativeError as error:
            raise _native_error(error, signals=names) from None
        stopped = np.frombuffer(generator.stopped(), dtype=bool)
        return generator, tuple(signals.values()), stopped

    def _scipy_generator(self, parameters):
        """The generator without signals as a SciPy CSC matrix (sink last),
        and the stop flags."""
        generator, _, stopped = self._generator(parameters, {})
        return _sparse(generator), stopped

    def _kick(self, state, changes, work):
        if not changes:
            return state
        operations = [
            (
                i,
                _OPERATIONS[kind],
                _clamp(value) if kind != "partition" else 0,
                value if kind == "partition" else 0.0,
            )
            for i, (kind, value) in changes.items()
        ]
        try:
            data, work.count = self._topology.kick(
                state.tobytes(), operations, work.count, work.native_limit
            )
        except _native.NativeError as error:
            raise _native_error(error, work) from None
        return _vector(data)

    def solve(
        self,
        end,
        report_times,
        *,
        initial=None,
        parameters=None,
        signals=None,
        tolerance=1e-8,
        max_step=np.inf,
        solver="bdf",
        max_work=200_000,
        max_output_bytes=256 << 20,
        on_report=None,
        keep_snapshots=True,
    ):
        """Propagate from a prepared law or a compatible report-time checkpoint.

        Constant segments use positive uniformization. Inhomogeneous
        segments use SciPy's BDF (``solver="bdf"``, the one solve that needs
        SciPy) or the native adaptive exponential midpoint rule. Series-tail,
        leakage, minimum probability and mass defect are retained separately.
        The series-tail bound covers uniformization only; BDF and midpoint
        segments add no integration error estimate, and BDF can return
        slightly negative probabilities. A checkpoint continues only under
        the same parameters and signals; change parameters at a time point
        with a protocol action. An optional on_report callback receives each
        snapshot in order. Pass keep_snapshots=False when a reduction suffices.
        """
        if on_report is not None and not callable(on_report):
            raise TypeError("on_report must be a callable")
        times = np.asarray(tuple(report_times), dtype=float)
        max_step = float(max_step)
        if (
            not isfinite(end)
            or end <= 0
            or not 1e-14 <= tolerance <= 1e-2
            or np.isnan(max_step)
            or max_step <= 0
            or solver not in ("bdf", "midpoint")
            or len(times) == 0
            or not np.isfinite(times).all()
            or times[0] < 0
            or times[-1] > end
            or np.any(np.diff(times) <= 0)
        ):
            raise ValueError(
                "invalid projection horizon, report grid or solver controls"
            )
        if 16 * len(self.states) * len(times) > index(max_output_bytes):
            raise RuntimeError("projection output budget exceeded")
        parameters, signals = dict(parameters or {}), dict(signals or {})
        if not all(isinstance(s, Signal) for s in signals.values()):
            raise TypeError("a time-dependent rate needs a Signal")
        model_sha = sha256(
            model_bytes(
                {
                    **self._document,
                    "parameters": {**self._document["parameters"], **parameters},
                }
            )
        ).hexdigest()
        # A Python callable can carry mutable state. A checkpoint is private to its
        # in-memory signal objects as well as the numeric model binding.
        table_hashes = {
            name: signal._table_sha256
            for name, signal in signals.items()
            if signal._table_sha256 is not None
        }
        binding = (
            model_sha
            + ":"
            + ",".join(
                f"{k}={table_hashes[k] if k in table_hashes else id(v)}"
                for k, v in sorted(signals.items())
            )
        )
        start, actions_applied = 0.0, 0
        state = np.zeros(len(self.states) + 1)
        frozen = np.zeros(len(self.states))
        tail = 0.0
        if isinstance(initial, ProjectionState):
            if initial.graph_id != self._id or initial.binding_id != binding:
                raise ValueError(
                    "checkpoint does not match this projection and binding"
                )
            start, actions_applied = initial.time, initial.actions_applied
            frozen[:] = initial.stopped_probabilities
            state[:-1] = initial.probabilities - frozen
            state[-1], tail = initial.physical_leakage, initial.series_tail_bound
        else:
            for key, value in (
                self._initial if initial is None else self._law(initial)
            ).items():
                at = self._index(key)
                if at is None:
                    raise ValueError(
                        "initial support requires a new prepared projection"
                    )
                state[at] = value
        if times[0] < start or start > end:
            raise ValueError("reports must follow the restart time")
        actions = self._actions
        if actions and actions[-1][0] > end:
            # A shorter horizon intentionally returns an intermediate checkpoint.
            actions = tuple(action for action in actions if action[0] <= end)
        for _, _, edits in actions[:actions_applied]:
            parameters.update(edits)
        if signals and solver == "bdf":
            _require_scipy("solver='bdf'")
        generator, forcing, stopped = self._generator(parameters, signals)
        matrices = None
        boundaries = {start, float(end)}
        boundaries.update(t for t, _, _ in actions[actions_applied:] if t >= start)
        boundaries.update(
            t
            for signal in signals.values()
            for t in signal.change_times
            if start < t < end
        )
        output = []
        work = _Work(max_work)
        report_index = 0

        def snapshot(t, vector):
            probabilities = vector[:-1] + frozen
            checkpoint_stopped = np.where(stopped, probabilities, frozen)
            probabilities.setflags(write=False)
            checkpoint_stopped.setflags(write=False)
            value = ProjectionState(
                t,
                probabilities,
                checkpoint_stopped,
                float(vector[-1]),
                tail,
                float(1.0 - (probabilities.sum() + vector[-1])),
                float(min(vector.min(), probabilities.min())),
                self._id,
                actions_applied,
                binding,
            )
            if on_report:
                on_report(value)
            if keep_snapshots:
                output.append(value)

        for t in sorted(boundaries):
            stop_reports = int(np.searchsorted(times, t, side="left"))
            interior = times[report_index:stop_reports]
            if t > start:
                sample_times = np.r_[interior, t]
                if forcing and solver == "bdf":
                    if matrices is None:
                        matrices = (
                            _sparse(generator),
                            tuple(_sparse(generator, k) for k in range(len(forcing))),
                        )
                    values = _bdf(
                        matrices,
                        forcing,
                        state,
                        start,
                        t,
                        tolerance,
                        max_step,
                        work,
                        sample_times,
                    )
                    for at, vector in zip(interior, values[:-1], strict=True):
                        snapshot(float(at), vector)
                    state = values[-1]
                else:
                    previous = start
                    for at in sample_times:
                        if forcing:
                            state = _midpoint(
                                generator,
                                forcing,
                                state,
                                previous,
                                at,
                                tolerance,
                                max_step,
                                work,
                            )
                        else:
                            state, increment = _uniform(
                                generator,
                                state,
                                at - previous,
                                tolerance / (len(times) + len(boundaries)),
                                work,
                            )
                            tail += increment
                        if at < t:
                            snapshot(float(at), state)
                        previous = at
            report_index = stop_reports
            while actions_applied < len(actions) and actions[actions_applied][0] == t:
                # Preserve previously stopped mass through interventions.
                if t > 0.0 or actions_applied:
                    frozen[stopped] += state[:-1][stopped]
                    state[:-1][stopped] = 0.0
                _, changes, edits = actions[actions_applied]
                state = self._kick(state, changes, work)
                if edits:
                    parameters.update(edits)
                    generator, forcing, stopped = self._generator(parameters, signals)
                    matrices = None
                actions_applied += 1
            if report_index < len(times) and times[report_index] == t:
                snapshot(t, state)
                report_index += 1
            start = t

        return ProjectionResult(
            self.species,
            self.states,
            tuple(output),
            work.count,
            solver if forcing else "uniformization",
            model_sha,
            sha256(model_bytes(self._protocol)).hexdigest(),
            table_hashes,
        )


def _uniform(generator, p, duration, tol, work):
    """Native uniformization of the constant part; see the Rust `uniform`."""
    try:
        data, tail, work.count = generator.uniform(
            p.tobytes(), float(duration), float(tol), work.count, work.native_limit
        )
    except _native.NativeError as error:
        raise _native_error(error, work) from None
    return _vector(data), tail


def _midpoint(generator, signals, p, start, end, tol, max_step, work):
    """Native adaptive exponential midpoint rule; see the Rust `midpoint`."""
    try:
        data, work.count = generator.midpoint(
            p.tobytes(),
            float(start),
            float(end),
            float(tol),
            float(max_step),
            [_source(s) for s in signals],
            work.count,
            work.native_limit,
        )
    except _native.NativeError as error:
        raise _native_error(error, work) from None
    return _vector(data)


def _bdf(matrices, signals, p, start, end, tol, max_step, work, reports):
    """SciPy's BDF on the native generator matrices."""
    from scipy.integrate import solve_ivp

    constant, forcing = matrices

    def weights(t):
        # Integrate to a discontinuity using the incoming hazard. The next
        # interval starts at the boundary with the right-continuous hazard.
        at = np.nextafter(end, start) if t == end else t
        values = [float(s.evaluate(at)) for s in signals]
        if any(v < 0 or not isfinite(v) for v in values):
            raise ValueError(
                "time-dependent reaction rate must be finite and nonnegative"
            )
        return values

    def matrix_at(t):
        q = constant
        for v, matrix in zip(weights(t), forcing, strict=True):
            q = q + v * matrix
        return q

    def rhs(t, vector):
        work.charge()
        result = constant @ vector
        for v, matrix in zip(weights(t), forcing, strict=True):
            result += v * (matrix @ vector)
        return result

    result = solve_ivp(
        rhs,
        (start, end),
        p,
        method="BDF",
        t_eval=reports,
        jac=lambda t, _: matrix_at(t),
        rtol=tol,
        atol=tol * 1e-4,
        max_step=max_step,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.y.T


def prepare_fsp(
    model,
    bounds,
    *,
    protocol=None,
    initial=None,
    max_states=100_000,
    max_edges=2_000_000,
):
    """Prepare a reusable projection for a closed model.

    The native engine lowers the model, enumerates the states, builds the
    generator and propagates the law; only ``solve(solver="bdf")`` with
    signals needs SciPy (the ``projection`` extra).
    Existing :meth:`Model.fsp` is the native uniformization for one fixed
    model without protocols, signals or initial mixtures. A prepared
    projection holds no per-solve state, so threads can share it. It cannot
    be pickled; prepare one per process.
    """
    return PreparedFsp(
        model,
        bounds,
        protocol=protocol,
        initial=initial,
        max_states=max_states,
        max_edges=max_edges,
    )
