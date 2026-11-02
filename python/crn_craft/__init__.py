"""Small byte-preserving JSON interface to the native CRN-craft core."""

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


def _call(method, *args, **kwargs):
    try:
        return json.loads(method(*args, **kwargs))
    except _native.NativeError as error:
        text = str(error)
        try:
            value = json.loads(text)
        except ValueError:
            value = [{"code": "CRN-E000", "path": "", "message": text}]
        raise CrnError(value if isinstance(value, list) else [value]) from None


def model_bytes(model):
    if isinstance(model, bytes):
        return model
    if isinstance(model, str):
        return model.encode("utf-8")
    return (json.dumps(model, allow_nan=False, separators=(",", ":")) + "\n").encode()


def describe():
    return json.loads(_native.describe_json())


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
    """Prepare a legacy Direct request once; retain its exact bytes and hashes."""
    try:
        return PreparedSimulation(
            _native.prepare_json(model_bytes(model), model_bytes(run), cached)
        )
    except _native.NativeError as error:
        raise CrnError(json.loads(str(error))) from None


def ssa_spec(
    model,
    *,
    aggregator="direct_cr",
    waiting_time="inverse",
    max_sampling_steps=1000000000,
    **kwargs,
):
    """Construct a versioned alternative exact-SSA request, bound to model bytes."""
    request = run_spec(model, **kwargs)
    request["schema"] = "crn-craft.ssa-run/0.1"
    request["method"] = {
        "kind": "ssa",
        "aggregator": aggregator,
        "waiting_time": waiting_time,
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
        raise CrnError(json.loads(str(error))) from None


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
