"""Pure-Python/SciPy reference of the prepared projection (not shipped).

A pure-Python implementation of ``crn_craft.prepare_fsp`` with the same
API and error messages, kept as an independent oracle for
``validation/check_projection_native.py``. Only model lowering and rate
evaluation use the native extension; enumeration, generator assembly,
interventions and propagation are Python, NumPy and SciPy. ``Signal`` and
``ProjectionState`` are the package's own classes, so signals and
checkpoints are shared with the native implementation.
"""

import json
from collections import defaultdict
from collections.abc import Mapping
from hashlib import sha256
from itertools import product
from math import comb, fsum, isfinite, prod
from operator import index

import numpy as np

from crn_craft import _native
from crn_craft._common import model_bytes
from crn_craft.projection import ProjectionResult, ProjectionState, Signal
from crn_craft.protocol import Protocol
from crn_craft.rules import _as_model_1_0, _crn_error, _unique_keys

_MAX_COUNT = 2**53


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


class _Work:
    def __init__(self, limit):
        self.limit = index(limit)
        self.count = 0

    def charge(self, amount=1):
        self.count += amount
        if self.count > self.limit:
            raise RuntimeError(
                f"projection work limit of {self.limit} exceeded; raise max_work"
            )


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


def _target(state, changes, bounds):
    target = list(state)
    for i, delta in changes:
        target[i] += delta
        if target[i] > bounds[i]:
            return None
    return tuple(target)


def _outcomes(count, operation):
    kind, value = operation
    if kind == "partition":
        return range(count + 1) if 0 < value < 1 else (count if value else 0,)
    return (value if kind == "set" else count + value,)


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
        _require_scipy("prepare_fsp")
        from scipy import sparse

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
        self._initialize_graph()
        state_bytes = self.states.astype("<i8", copy=False).tobytes()
        self._id = sha256(
            model_bytes([self.species, self._protocol, self._document]) + state_bytes
        ).hexdigest()
        self.states.setflags(write=False)
        self._constant = tuple(
            i for i, c in enumerate(self._channels) if c["rate"] is not None
        )
        self._variable = tuple(
            i for i, c in enumerate(self._channels) if c["rate"] is None
        )
        self._sparse = sparse

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

    def _initialize_graph(self):
        states, identities = [], {}
        origins, destinations, edges, multiplicities = [], [], [], []
        sink = self._max_states
        width = len(self._bounds)
        touched = 0

        def admit(state):
            if len(state) != width or any(n < 0 for n in state):
                raise ValueError("invalid initial or reachable counts")
            if any(n > bound for n, bound in zip(state, self._bounds, strict=True)):
                return sink
            if state not in identities:
                if len(states) >= sink:
                    raise RuntimeError("reachable-state budget exceeded")
                # Reserve the state table and the indices used for lookups.
                if (len(states) + 1) * (128 + 40 * width) > 256 << 20:
                    raise RuntimeError(
                        "projection state table exceeds the workspace budget"
                    )
                identities[state] = len(states)
                states.append(state)
            return identities[state]

        for state in sorted(self._initial):
            if admit(state) == sink:
                raise ValueError("initial law lies outside the projection")

        cursor = 0
        for changes in [None, *(changes for _, changes, _ in self._actions if changes)]:
            if changes:
                # Every state reached before this action may be its input.
                for state in tuple(states):
                    indices = tuple(changes)
                    supports = [_outcomes(state[i], changes[i]) for i in indices]
                    if any(len(s) and s[0] < 0 for s in supports):
                        continue
                    if prod(map(len, supports)) > self._max_edges * 4 - touched:
                        raise RuntimeError("intervention expansion exceeds work budget")
                    for values in product(*supports):
                        touched += 1
                        if touched > self._max_edges * 4:
                            raise RuntimeError(
                                "intervention expansion exceeds work budget"
                            )
                        output = list(state)
                        for i, n in zip(indices, values, strict=True):
                            output[i] = n
                        admit(tuple(output))
            while cursor < len(states):
                state = states[cursor]
                for i, channel in enumerate(self._channels):
                    touched += 1
                    if touched > self._max_edges * 25:
                        raise RuntimeError("projection enumeration exceeds work budget")
                    needs = channel["needs"]
                    if any(state[j] < n for j, n in needs):
                        continue
                    dst = _target(state, channel["changes"], self._bounds)
                    if dst == state:
                        continue
                    if len(edges) >= self._max_edges:
                        raise RuntimeError("projection edge budget exceeded")
                    origins.append(cursor)
                    destinations.append(sink if dst is None else admit(dst))
                    edges.append(i)
                    multiplicities.append(
                        1.0 if channel["whole"] else _combinations(state, needs)
                    )
                cursor += 1

        self.states = np.asarray(states, dtype=np.int64)
        self._list = tuple(states)
        self._ids = identities
        self._src = np.asarray(origins, dtype=np.int64)
        self._dst = np.asarray(destinations, dtype=np.int64)
        self._dst[self._dst == sink] = len(states)
        self._edges = np.asarray(edges, dtype=np.int64)
        self._multiplicities = np.asarray(multiplicities, dtype=float)
        self._shape = (len(states) + 1,) * 2

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

    def _matrix(self, rates, selected):
        src, dst, data = self._src[selected], self._dst[selected], rates[selected]
        return self._sparse.coo_matrix(
            (np.r_[data, -data], (np.r_[dst, src], np.r_[src, src])),
            shape=self._shape,
        ).tocsc()

    def _generator(self, parameters, signals):
        native, layout = self._binding(parameters)
        n = len(self.states)
        channels = layout["channels"]
        constants = np.array(
            [0.0 if c["rate"] is None else c["rate"] for c in channels]
        )
        factors = np.array([c["factor"] for c in channels])
        basis = self._multiplicities * factors[self._edges]
        with np.errstate(over="ignore", invalid="ignore"):
            rates = constants[self._edges] * basis
        stopped = np.zeros(n, dtype=bool)
        if self._variable or layout["has_stop"]:
            if n * max(1, len(self._variable)) > 10_000_000:
                raise RuntimeError(
                    "state-dependent rate table exceeds the workspace budget"
                )
            try:
                flags, data = native.evaluate(self._list, self._variable)
            except _native.NativeError as error:
                raise _crn_error(error) from None
            stopped[:] = flags
            data = np.asarray(data).reshape(n, len(self._variable))
            for i, channel in enumerate(self._variable):
                select = self._edges == channel
                rates[select] = data[self._src[select], i]
        signal_matrices = []
        for name, signal in signals.items():
            ids = [i for i, c in enumerate(channels) if c["name"] == name]
            if not ids or any(i not in self._constant for i in ids):
                raise ValueError(
                    f"{name!r} is not a rule with a count-independent rate"
                )
            mask = np.isin(self._edges, ids) & ~stopped[self._src]
            signal_matrices.append((signal, self._matrix(basis, mask)))
            rates[np.isin(self._edges, ids)] = 0.0
        rates[stopped[self._src]] = 0.0
        if np.any(rates < 0) or not np.isfinite(rates).all():
            raise ValueError("invalid rate in the projection")
        keep = rates > 0
        return self._matrix(rates, keep), tuple(signal_matrices), stopped

    def _kick(self, state, changes, work):
        from scipy.stats import binom

        if not changes:
            return state
        mass = {self._list[j]: state[j] for j in np.flatnonzero(state[:-1])}
        escape = state[-1]
        for i, (kind, value) in changes.items():
            next_mass = defaultdict(float)
            for counts, amount in mass.items():
                possible = _outcomes(counts[i], (kind, value))
                work.charge(len(possible))
                if possible[0] < 0:
                    raise ValueError(
                        f"intervention removes more {self.species[i]} than available"
                    )
                weights = (
                    binom.pmf(np.arange(len(possible)), counts[i], value)
                    if len(possible) > 1
                    else (1.0,)
                )
                for count, weight in zip(possible, weights, strict=True):
                    if count > self._bounds[i]:
                        escape += amount * weight
                    elif weight:
                        target = counts[:i] + (count,) + counts[i + 1 :]
                        next_mass[target] += amount * weight
                        if len(next_mass) > self._max_states:
                            raise RuntimeError("intervention workspace budget exceeded")
            mass = next_mass
        new = np.zeros_like(state)
        new[-1] = escape
        for counts, value in mass.items():
            if value:
                try:
                    new[self._ids[counts]] = value
                except KeyError as error:
                    raise ValueError(
                        "intervention reaches an unprepared state"
                    ) from error
        return new

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
        segments use SciPy's BDF or adaptive exponential midpoint. Series-tail,
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
                try:
                    state[self._ids[key]] = value
                except KeyError as error:
                    raise ValueError(
                        "initial support requires a new prepared projection"
                    ) from error
        if times[0] < start or start > end:
            raise ValueError("reports must follow the restart time")
        actions = self._actions
        if actions and actions[-1][0] > end:
            # A shorter horizon intentionally returns an intermediate checkpoint.
            actions = tuple(action for action in actions if action[0] <= end)
        for _, _, edits in actions[:actions_applied]:
            parameters.update(edits)
        generator, forcing, stopped = self._generator(parameters, signals)
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
                    values = _forced(
                        generator,
                        forcing,
                        state,
                        start,
                        t,
                        tolerance,
                        max_step,
                        solver,
                        work,
                        reports=sample_times,
                    )
                    for at, vector in zip(interior, values[:-1], strict=True):
                        snapshot(float(at), vector)
                    state = values[-1]
                else:
                    previous = start
                    for at in sample_times:
                        if forcing:
                            state = _forced(
                                generator,
                                forcing,
                                state,
                                previous,
                                at,
                                tolerance,
                                max_step,
                                solver,
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


def _combinations(state, needs):
    # The native engine bounds the work for large combinatorial factors.
    # Two nontrivial arguments above this range already overflow f64.
    result = 1
    for i, n in needs:
        choose = min(n, state[i] - n)
        if choose > 1024:
            raise ValueError("reactant combinatorics exceed probability solver range")
        result *= comb(state[i], choose)
    return float(result)


def _uniform(matrix, p, duration, tol, work):
    from scipy.special import gammainc

    rate = float(np.max(-matrix.diagonal()))
    if rate == 0.0:
        return p, 0.0
    pieces = max(1, int(np.ceil(rate * duration / 16.0)))
    if pieces > work.limit - work.count:
        raise RuntimeError("uniformization exceeds the work budget")
    transition = matrix / rate
    transition.setdiag(1.0 + transition.diagonal())
    mu = rate * duration / pieces
    tail = 0.0
    for _ in range(pieces):
        term = p
        weight = float(np.exp(-mu))
        p = weight * term
        k = 0
        while (remaining := float(gammainc(k + 1, mu))) > tol / pieces:
            work.charge()
            k += 1
            term = transition @ term
            weight *= mu / k
            p += weight * term
        tail += remaining
    return p, tail


def _forced(
    constant, forcing, p, start, end, tol, max_step, solver, work, *, reports=None
):
    def weights(t):
        # Integrate to a discontinuity using the incoming hazard. The next
        # interval starts at the boundary with the right-continuous hazard.
        at = np.nextafter(end, start) if t == end else t
        values = [float(s.evaluate(at)) for s, _ in forcing]
        if any(v < 0 or not isfinite(v) for v in values):
            raise ValueError(
                "time-dependent reaction rate must be finite and nonnegative"
            )
        return values

    def matrix_at(t):
        q = constant
        for v, (_, matrix) in zip(weights(t), forcing, strict=True):
            q = q + v * matrix
        return q

    if solver == "bdf":
        from scipy.integrate import solve_ivp

        def rhs(t, vector):
            work.charge()
            result = constant @ vector
            for v, (_, matrix) in zip(weights(t), forcing, strict=True):
                result += v * (matrix @ vector)
            return result

        result = solve_ivp(
            rhs,
            (start, end),
            p,
            method="BDF",
            t_eval=[end] if reports is None else reports,
            jac=lambda t, _: matrix_at(t),
            rtol=tol,
            atol=tol * 1e-4,
            max_step=max_step,
        )
        if not result.success:
            raise RuntimeError(result.message)
        return result.y[:, -1] if reports is None else result.y.T

    from scipy.sparse.linalg import expm_multiply

    def advance(at, h, vector):
        work.charge()
        return expm_multiply(h * matrix_at(at + h / 2), vector)

    t, h = start, min(end - start, max_step)
    while t < end:
        # A step that reaches the end within rounding ends there, and a
        # remainder of less than two steps is halved.
        remaining = end - t
        last = remaining <= h * (1 + 1e-9)
        if last:
            h = remaining
        elif remaining < 2 * h:
            h = remaining / 2
        if t + h <= t:
            raise RuntimeError("midpoint step cannot advance the clock")
        full = advance(t, h, p)
        half = advance(t + h / 2, h / 2, advance(t, h / 2, p))
        gap = float(np.linalg.norm(full - half, 1))
        goal = tol * h / (end - start)
        if gap <= goal:
            p = half
            t = end if last else t + h
        elif goal <= 64 * np.finfo(float).eps:
            # The step error estimate cannot resolve this goal in floating point.
            raise RuntimeError(
                "midpoint tolerance is below the rounding level of its step "
                "error estimate; use a looser tolerance or solver='bdf'"
            )
        factor = 2.0 if gap == 0 else np.clip(0.9 * (goal / gap) ** (1 / 3), 0.2, 2.0)
        h = min(max_step, h * factor)
    return p


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

    The native core lowers the model and evaluates rates. Python enumerates
    the states, and SciPy (the ``projection`` extra) propagates the law.
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
