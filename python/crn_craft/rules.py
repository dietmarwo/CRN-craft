"""Prepared rules ensembles and read-only NumPy result columns."""

import hashlib
import json
from pathlib import Path

import numpy as np

from . import _crn_error, _native, _path, model_bytes

TERMINATIONS = (
    "completed", "stopped", "event_limit", "sampling_limit", "instance_limit",
    "wall_time_limit", "cancelled", "failed", "output_limit",
)


def load(path):
    """Load and validate exact model bytes for reuse in multiple requests."""
    data = Path(path).read_bytes()
    check = json.loads(_native.check_rules_json(data))
    if not check["valid"]:
        from . import CrnError

        raise CrnError(check["diagnostics"])
    return data


class Result:
    def __init__(self, metadata, columns):
        self.metadata = metadata["result"]
        self.record_diagnostics = dict(metadata["record_diagnostics"])
        self.records = {}
        for key, buffer in columns.items():
            dtype = {"started": "?", "termination": "u1", "stop_time": "<f8"}.get(
                key, "<u8"
            )
            self.records[key] = np.frombuffer(buffer, dtype=dtype)
        self.records["rule_events"] = self.records["rule_events"].reshape(
            len(self.records["id"]), len(self.metadata["rules"])
        )
        stats = self.metadata["statistics"]
        self.statistics = None if stats is None else {
            key: np.asarray(value) if isinstance(value, list) else value
            for key, value in stats.items()
        }


class PreparedEnsemble:
    def __init__(self, native):
        self._native = native

    def simulate(self, *, cancel_file=None):
        """Execute with native workers. Ctrl-C or cancel_file cancels the run."""
        try:
            metadata, columns = self._native.simulate(_path(cancel_file))
        except _native.NativeError as error:
            raise _crn_error(error) from None
        return Result(json.loads(metadata), columns)


def prepare(model, request):
    """Prepare a rules-run contract, preserving the supplied model bytes."""
    try:
        return PreparedEnsemble(
            _native.prepare_rules(model_bytes(model), model_bytes(request))
        )
    except _native.NativeError as error:
        raise _crn_error(error) from None


def ensemble(
    model, *, end, report_times, replicates=1, seed=0, threads=1,
    keep_paths=0, keep_instances=0, max_events=1_000_000,
    max_sampling_steps=8_000_000, max_instances=100_000,
    max_output_bytes=64 << 20, wall_time_s=None,
):
    """Build and prepare a rules request with explicit work/output budgets."""
    data = model_bytes(model)
    request = {
        "schema": "crn-craft.rules-run/0.1",
        "model": {"path": "model.rules.json", "sha256": hashlib.sha256(data).hexdigest()},
        "time": {"end": end, "report_times": list(report_times)},
        "ensemble": {"replicates": replicates, "root_seed": seed},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": keep_paths,
            "keep_instances": keep_instances,
        },
        "budget": {
            "threads": threads,
            "max_events_per_replicate": max_events,
            "max_sampling_steps_per_replicate": max_sampling_steps,
            "max_instances": max_instances,
            "max_output_bytes": max_output_bytes,
            "wall_time_s": wall_time_s,
        },
    }
    return prepare(data, request)
