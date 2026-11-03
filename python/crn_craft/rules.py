"""Prepared rules ensembles and read-only NumPy result columns."""

import hashlib
import json
import operator
from pathlib import Path

import numpy as np

from . import _crn_error, _native, _path, model_bytes

TERMINATIONS = (
    "completed",
    "stopped",
    "event_limit",
    "sampling_limit",
    "instance_limit",
    "wall_time_limit",
    "cancelled",
    "failed",
    "output_limit",
)


def describe():
    """Return the native rules capabilities and limits."""
    return json.loads(_native.describe_rules_json())


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
        self.observations = self.metadata.get("observations", [])
        stats = self.metadata["statistics"]
        self.statistics = (
            None
            if stats is None
            else {
                key: np.asarray(value) if isinstance(value, list) else value
                for key, value in stats.items()
            }
        )

    def event_columns(self, replicate):
        """Return time, rule and observable arrays for a retained path log.

        Rule index -1 denotes a leap. Participant identities and leap
        multiplicities are available in the corresponding observation entries.
        """
        for observation in self.observations:
            if observation["id"] == replicate:
                entries = observation["log"]
                return {
                    "time": np.array([e["time"] for e in entries], dtype=np.float64),
                    "rule": np.array(
                        [-1 if e["rule"] is None else e["rule"] for e in entries],
                        dtype=np.int64,
                    ),
                    "observables": np.array(
                        [e["observables"] for e in entries], dtype=np.float64
                    ).reshape(len(entries), len(self.metadata["observables"])),
                }
        raise KeyError(replicate)

    def tracker_columns(self, index=0):
        """Return scalar tracker fields as arrays, one entry per observed path.

        Paths that failed before tracking began carry no observation.
        """
        if index < 0:
            raise ValueError("tracker index must be nonnegative")
        count = max((len(ob["trackers"]) for ob in self.observations), default=0)
        if index >= count:
            raise IndexError(f"tracker index {index} with {count} trackers")
        rows = [(ob["id"], ob["trackers"][index]) for ob in self.observations]
        out = {"id": np.array([i for i, _ in rows], dtype=np.uint64)}
        for key in ("appearances", "extinctions", "omitted_episodes"):
            out[key] = np.array([t[key] for _, t in rows], dtype=np.uint64)
        for key in (
            "last_appearance",
            "last_extinction",
            "first_passage",
            "minimum",
            "maximum",
            "integral",
            "omitted_duration",
            "omitted_peak",
        ):
            out[key] = np.array([t[key] for _, t in rows], dtype=np.float64)
        return out


class PreparedEnsemble:
    def __init__(self, native):
        self._native = native

    def bind_parameters(self, overrides):
        """Prepare an independent candidate with the same ensemble settings.

        The candidate's model and run hashes cover the re-serialized bound
        model and request, not the original file.
        """
        try:
            return PreparedEnsemble(self._native.bind_parameters(dict(overrides)))
        except _native.NativeError as error:
            raise _crn_error(error) from None

    def simulate(
        self, *, cancel_file=None, callback=None, batch_callback=None, batch_size=256
    ):
        """Run native workers.

        A live callback(id, event) returning False stops that replicate; any
        other return value continues. Workers preserve per-replicate order,
        and calls from different replicates interleave. An exception in the
        callback cancels the run and is raised here; partial results are
        discarded. A callback that blocks also blocks Ctrl-C. A
        batch_callback(id, events) replays the retained log after execution
        and cannot affect the run.
        """
        if batch_size < 1:
            raise ValueError("batch_size must be positive")

        def dispatch(replicate, event):
            return callback(replicate, json.loads(event))

        try:
            metadata, columns = self._native.simulate(
                _path(cancel_file), None if callback is None else dispatch
            )
        except _native.NativeError as error:
            raise _crn_error(error) from None
        result = Result(json.loads(metadata), columns)
        if batch_callback is not None:
            for observation in result.observations:
                log = observation["log"]
                for start in range(0, len(log), batch_size):
                    batch_callback(observation["id"], log[start : start + batch_size])
        return result


def prepare(model, request):
    """Prepare a rules-run contract, preserving the supplied model bytes."""
    try:
        return PreparedEnsemble(
            _native.prepare_rules(model_bytes(model), model_bytes(request))
        )
    except _native.NativeError as error:
        raise _crn_error(error) from None


def ensemble(
    model,
    *,
    end,
    report_times,
    replicates=1,
    seed=0,
    threads=1,
    keep_paths=0,
    keep_instances=0,
    max_events=1_000_000,
    max_sampling_steps=8_000_000,
    max_instances=100_000,
    max_output_bytes=64 << 20,
    wall_time_s=None,
    trackers=(),
    event_log=None,
    method=None,
):
    """Build and prepare a rules request with explicit work/output budgets."""
    data = model_bytes(model)
    # NumPy scalars and arrays are converted; counts must be integral.
    count = operator.index
    request = {
        "schema": "crn-craft.rules-run/0.1",
        "model": {
            "path": "model.rules.json",
            "sha256": hashlib.sha256(data).hexdigest(),
        },
        "time": {"end": float(end), "report_times": [float(t) for t in report_times]},
        "ensemble": {"replicates": count(replicates), "root_seed": count(seed)},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": count(keep_paths),
            "keep_instances": count(keep_instances),
        },
        "budget": {
            "threads": count(threads),
            "max_events_per_replicate": count(max_events),
            "max_sampling_steps_per_replicate": count(max_sampling_steps),
            "max_instances": count(max_instances),
            "max_output_bytes": count(max_output_bytes),
            "wall_time_s": None if wall_time_s is None else float(wall_time_s),
        },
    }
    if trackers or event_log is not None or method is not None:
        # Present fields must hold values; omit the ones that were not given.
        request["schema"] = "crn-craft.rules-run/0.2"
        if trackers:
            request["trackers"] = list(trackers)
        if event_log is not None:
            request["event_log"] = event_log
        if method is not None:
            request["method"] = {"kind": method} if isinstance(method, str) else method
        if (
            isinstance(request.get("method"), dict)
            and request["method"].get("kind") == "mtdna_direct"
        ):
            request["schema"] = "crn-craft.rules-run/0.3"
        elif (
            isinstance(request.get("method"), dict)
            and request["method"].get("kind") == "static_direct"
            and {"waiting_time", "sparse_trackers"} & request["method"].keys()
        ):
            request["schema"] = "crn-craft.rules-run/0.4"
    return prepare(data, request)
