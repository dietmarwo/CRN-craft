"""Small byte-preserving JSON interface to the native CRN-craft core."""

from collections.abc import Mapping, Set
import numbers
import hashlib
import importlib.metadata
import json
from . import _native

try:
    __version__ = importlib.metadata.version("crn-craft")
except importlib.metadata.PackageNotFoundError:  # a source tree, not installed
    __version__ = _native.__version__


class CrnError(ValueError):
    """A request the native core rejected.

    ``diagnostics`` is the list of CRN-craft diagnostics (dictionaries with
    ``code``, ``path``, ``message`` and ``hint``). It is a ``ValueError``, so
    existing ``except ValueError`` handlers keep working.
    """

    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        first = diagnostics[0] if diagnostics else {}
        message = f"{first.get('code', 'CRN')}: {first.get('message', 'rejected')}"
        if first.get("path"):
            message += f" at {first['path']}"
        if len(diagnostics) > 1:
            message += f" (and {len(diagnostics) - 1} more)"
        super().__init__(message)


def _crn_error(error):
    text = str(error)
    try:
        value = json.loads(text)
    except ValueError:
        value = [{"code": "CRN-E000", "path": "", "message": text}]
    return CrnError(value if isinstance(value, list) else [value])


def _call(method, *args, **kwargs):
    try:
        return json.loads(method(*args, **kwargs))
    except _native.NativeError as error:
        raise _crn_error(error) from None


def model_bytes(model):
    if isinstance(model, bytes):
        return model
    if isinstance(model, str):
        return model.encode("utf-8")
    return (json.dumps(model, allow_nan=False, separators=(",", ":")) + "\n").encode()


def describe():
    return json.loads(_native.describe_json())


def convert(model):
    """Convert a network model (``crn-craft.model/0.1``) into a rules model.

    ``model`` is bytes, text or a dictionary. Returns the converted model's
    exact bytes (pretty-printed JSON with a trailing newline), the same
    bytes ``crn-craft convert MODEL`` prints; run them with
    :mod:`crn_craft.rules`. Every species becomes a type and an observable
    of the same name, every reaction a rule with the same propensity. A
    model the rules language cannot hold yet (more than 64 species, more
    than two reactant molecules in a reaction, a population limit, ...)
    raises ``CrnError`` with ``CRN-R016``; an invalid network model raises
    ``CrnError`` with its validation diagnostics.
    """
    try:
        return _native.convert_network(model_bytes(model))
    except _native.NativeError as error:
        raise _crn_error(error) from None


def check(model, run=None):
    return json.loads(
        _native.check_json(
            model_bytes(model), None if run is None else model_bytes(run)
        )
    )


def analyze(model):
    return json.loads(_native.analyze_json(model_bytes(model)))


def simulate(model, run, *, cached=False, cancel_file=None):
    """Return the numerical result; inspect termination before using statistics.

    Ctrl-C (in the main thread) cancels the run and raises KeyboardInterrupt.
    Creating ``cancel_file`` stops the run with termination ``cancelled``.
    Rejected requests raise ``CrnError``.
    """
    method = _native.simulate_cached_json if cached else _native.simulate_json
    return _call(method, model_bytes(model), model_bytes(run), _path(cancel_file))


class PreparedSimulation:
    """Validated immutable inputs reusable across independent executions.

    Repeated calls restart the same fixed seeds. Prepare a new request to change
    inputs. Each call owns its execution state, including when called concurrently.
    """

    def __init__(self, native):
        self._prepared = native

    def simulate(self, *, cancel_file=None):
        return _call(self._prepared.simulate_json, _path(cancel_file))


def prepare(model, run, *, cached=False):
    """Validate and prepare a legacy Direct request once for repeated runs."""
    try:
        return PreparedSimulation(
            _native.prepare_json(model_bytes(model), model_bytes(run), cached)
        )
    except _native.NativeError as error:
        raise _crn_error(error) from None


_DEFAULT_POLICY = object()


def ssa_spec(
    model,
    *,
    aggregator="direct_cr",
    waiting_time="inverse",
    max_sampling_steps=1000000000,
    fast=_DEFAULT_POLICY,
    **kwargs,
):
    """Construct an exact-SSA request bound to model bytes.

    The request is ``crn-craft.ssa-run/0.3``, or ``0.4`` for
    ``aggregator="next_reaction"``. ``fast`` (default ``True``) emulates
    JumpProcesses.jl: Xoshiro256++ words, running totals without lost-increment
    or 53-bit resolution checks, and stop polls per batch of events.
    ``fast=False`` selects the checked exact policy (ChaCha8 words, both checks,
    stop polls before every random word), whose paths equal those of
    ``crn-craft.ssa-run/0.1`` requests. Next reaction runs only the checked
    policy, so its default is ``False``. The request always states
    ``method.fast`` explicitly.
    """
    if fast is _DEFAULT_POLICY:
        fast = aggregator != "next_reaction"
    if not isinstance(fast, bool):
        raise TypeError("fast must be True or False")
    request = run_spec(model, **kwargs)
    request["schema"] = (
        "crn-craft.ssa-run/0.4"
        if aggregator == "next_reaction"
        else "crn-craft.ssa-run/0.3"
    )
    request["method"] = {
        "kind": "ssa",
        "aggregator": aggregator,
        "waiting_time": waiting_time,
        "fast": fast,
    }
    request["budget"]["max_sampling_steps_per_replicate"] = max_sampling_steps
    return request


def simulate_ssa(model, request, *, cancel_file=None):
    """Alternative exact SSA; inspect termination before using statistics."""
    return _call(
        _native.ssa_json, model_bytes(model), model_bytes(request), _path(cancel_file)
    )


def prepare_ssa(model, request):
    """Prepare immutable alternative-SSA inputs for repeated/concurrent execution."""
    try:
        return PreparedSimulation(
            _native.prepare_ssa_json(model_bytes(model), model_bytes(request))
        )
    except _native.NativeError as error:
        raise _crn_error(error) from None


class BoundParameterSimulation(PreparedSimulation):
    """Independent bound model/run snapshots and restartable prepared execution."""

    @property
    def model_bytes(self):
        return self._prepared.model_bytes()

    @property
    def request_bytes(self):
        return self._prepared.request_bytes()


class ParameterTemplate:
    """Fixed topology with independently bound parameter values.

    All declared parameters must be supplied, in ``parameter_names`` order or
    as an exact mapping. Initial populations/run settings are fixed. Each bind
    refreshes law memos/bounds and retains shared immutable network structure.
    """

    def __init__(self, native):
        self._template = native
        self._parameter_names = tuple(native.parameter_names())
        self._default_values = tuple(native.default_values())
        self._topology_sha256 = native.topology_sha256()

    @property
    def parameter_names(self):
        return self._parameter_names

    @property
    def default_values(self):
        return self._default_values

    @property
    def topology_sha256(self):
        return self._topology_sha256

    def bind(self, values, *, seed=None):
        """Bind one candidate to the fixed topology without recompiling it.

        ``values`` is a mapping with exactly the names in ``parameter_names``,
        or an ordered sequence (list, tuple, numpy array) in that order, which
        is sorted by name. Unordered collections such as sets are rejected.
        ``seed`` overrides the run's root seed. Returns a
        ``BoundParameterSimulation`` whose ``model_bytes`` and
        ``request_bytes`` replay the run. Raises ``CrnError``.
        """

        def invalid(message, code="CRN-E015", path="/parameters"):
            return CrnError(
                [
                    {
                        "code": code,
                        "path": path,
                        "message": message,
                        "hint": "Use every parameter_names entry exactly once with real numeric values.",
                    }
                ]
            )

        if seed is not None:
            if (
                isinstance(seed, bool)
                or not isinstance(seed, numbers.Integral)
                or not 0 <= seed < 2**64
            ):
                raise invalid(
                    "seed must be an unsigned 64-bit integer",
                    "CRN-E020",
                    "/ensemble/root_seed",
                )
            seed = int(seed)
        if isinstance(values, Mapping):
            if set(values) != set(self.parameter_names):
                missing = sorted(set(self.parameter_names) - set(values))
                unknown = sorted(map(str, set(values) - set(self.parameter_names)))
                raise invalid(
                    f"parameter names do not match template: missing {missing}, unknown {unknown}",
                    "CRN-E014",
                )
            values = [values[name] for name in self.parameter_names]
        elif isinstance(values, str | bytes | bytearray | Set):
            raise invalid(
                "parameter values must be an ordered numeric sequence or a mapping"
            )
        else:
            try:
                values = list(values)
            except TypeError as error:
                raise invalid(
                    "parameter values must be a sequence or mapping"
                ) from error
        if any(isinstance(x, bool) or not isinstance(x, numbers.Real) for x in values):
            raise invalid("parameter values must be real numbers, excluding booleans")
        try:
            floats = [float(x) for x in values]
        except (OverflowError, ValueError) as error:
            raise invalid(f"parameter value is not a finite float: {error}") from None
        try:
            return BoundParameterSimulation(self._template.bind(floats, seed))
        except _native.NativeError as error:
            raise _crn_error(error) from None

    def simulate(self, values, *, seed=None, cancel_file=None):
        """Bind ``values`` (see ``bind``) and run the ensemble; returns the result."""
        return self.bind(values, seed=seed).simulate(cancel_file=cancel_file)


def prepare_parameters(model, request):
    """Compile fixed topology once for a parameter search; no runtime compiler."""
    try:
        return ParameterTemplate(
            _native.prepare_parameters_json(model_bytes(model), model_bytes(request))
        )
    except _native.NativeError as error:
        raise _crn_error(error) from None


def mass_action_kernel_source(model):
    """Generate optional Rust propensity source; compile externally once per topology."""
    try:
        return _native.mass_action_kernel_source(model_bytes(model))
    except _native.NativeError as error:
        raise _crn_error(error) from None


def fsp(model, run, *, cancel_file=None):
    """Solve the finite state projection; rejected requests raise ``CrnError``."""
    return _call(
        _native.fsp_json, model_bytes(model), model_bytes(run), _path(cancel_file)
    )


def _path(value):
    return None if value is None else str(value)


def run_spec(
    model,
    *,
    end=100.0,
    save_times=None,
    replicates=64,
    seed=4242,
    threads=1,
    keep=0,
    max_events=1000000,
    max_output_bytes=10000000,
):
    """Construct a legacy-compatible SSA run bound to the exact model bytes."""
    return {
        "schema": "crn-craft.run/0.1",
        "model": {
            "path": "model.json",
            "sha256": hashlib.sha256(model_bytes(model)).hexdigest(),
        },
        "method": {"kind": "ssa", "aggregator": "direct"},
        "time": {
            "end": end,
            "save_times": save_times if save_times is not None else [0.0, end],
        },
        "ensemble": {"replicates": replicates, "root_seed": seed},
        "outputs": {
            "statistics": ["mean", "variance"] if replicates > 1 else ["mean"],
            "keep_trajectories": keep,
        },
        "budget": {
            "threads": threads,
            "max_events_per_replicate": max_events,
            "max_output_bytes": max_output_bytes,
        },
    }


def tau(model, request, *, cancel_file=None):
    """Approximate tau-leaping ensemble; rejected requests raise ``CrnError``."""
    return _call(
        _native.tau_json, model_bytes(model), model_bytes(request), _path(cancel_file)
    )


def load_example(name):
    from importlib.resources import files

    if name not in [
        "immigration-death",
        "enzyme-mechanism",
        "michaelis-menten",
        "repressilator",
        "vilar-noise",
    ]:
        raise ValueError("unknown example")
    return files("crn_craft").joinpath("models", name + ".model.json").read_bytes()
