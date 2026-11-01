"""Small byte-preserving JSON interface to the native CRN-craft core."""

import hashlib
import json
from . import _native


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


def simulate(model, run, *, cached=False):
    """Return the numerical result; inspect termination before using statistics."""
    method = _native.simulate_cached_json if cached else _native.simulate_json
    return json.loads(method(model_bytes(model), model_bytes(run)))


def fsp(model, run):
    return json.loads(_native.fsp_json(model_bytes(model), model_bytes(run)))


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


def tau(model, request):
    return json.loads(_native.tau_json(model_bytes(model), model_bytes(request)))


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
