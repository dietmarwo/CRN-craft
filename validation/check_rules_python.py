"""Exercise the installed wheel's rules API, buffers and callback contracts."""

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np

from crn_craft import rules
from check_rules import require

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require("thinning" in rules.describe()["methods"], "rules discovery missing")
    model = rules.load(ROOT / "models/rules/mtdna-three-class-cme.rules.json")
    common = {
        "end": 3.0,
        "report_times": [0.0, 1.0, 3.0],
        "replicates": 20,
        "seed": 1234,
        "method": "static_direct",
        "trackers": [{"observable": "mutants", "keep_episodes": 128}],
        "event_log": {"max_entries": 512},
    }
    expected = rules.ensemble(model, threads=1, **common).simulate()
    prepared = rules.ensemble(model, threads=20, **common)
    actual = prepared.simulate()
    for key in expected.records:
        require(
            actual.records[key].dtype == expected.records[key].dtype
            and actual.records[key].tobytes() == expected.records[key].tobytes()
            and not actual.records[key].flags.writeable,
            f"record column mismatch: {key}",
        )
    require(expected.observations == actual.observations, "thread-dependent log")
    require(
        np.array_equal(expected.statistics["mean"], actual.statistics["mean"]),
        "thread-dependent statistics",
    )
    require(actual.records["rule_events"].shape == (20, 8), "rule buffer shape")
    first = actual.observations[0]["log"]
    events = actual.event_columns(0)
    require(
        np.array_equal(events["time"], [e["time"] for e in first]), "event time column"
    )
    trackers = actual.tracker_columns()
    require(
        np.array_equal(
            trackers["extinctions"],
            [o["trackers"][0]["extinctions"] for o in actual.observations],
        ),
        "tracker count column",
    )
    received = [[] for _ in range(20)]
    batches = [[] for _ in range(20)]
    with_callback = prepared.simulate(
        callback=lambda i, e: received[i].append(e),
        batch_callback=lambda i, entries: batches[i].extend(entries),
        batch_size=13,
    )
    require(
        received == batches == [ob["log"] for ob in actual.observations],
        "callback sequence mismatch",
    )
    require(
        with_callback.metadata == actual.metadata, "observer changed numerical result"
    )
    stopped = prepared.simulate(callback=lambda _i, _event: False)
    require(
        np.all(stopped.records["events"] == 1)
        and stopped.metadata["stopped_replicates"] == 20,
        "callback did not stop the paths",
    )

    class CallbackError(Exception):
        pass

    def fail(_i, _event):
        raise CallbackError("callback error")

    try:
        prepared.simulate(callback=fail)
    except CallbackError:
        pass
    else:
        raise RuntimeError("callback exception was swallowed")

    with tempfile.TemporaryDirectory() as directory:
        flag = Path(directory) / "cancel"
        flag.touch()
        stopped = prepared.simulate(cancel_file=flag)
    require(
        stopped.metadata["termination"] == "cancelled" and stopped.statistics is None,
        "cancel file failed",
    )
    fixed = common | {"method": "mtdna_direct", "event_log": None}
    template = rules.ensemble(model, threads=8, **fixed)
    candidate = template.bind_parameters({"mu": 0.3, "sm": 1.8})
    a = candidate.simulate()
    b = candidate.simulate()
    fresh = json.loads(model)
    fresh["parameters"].update(mu=0.3, sm=1.8)
    c = rules.ensemble(fresh, threads=8, **fixed).simulate()
    require("ziggurat" in a.metadata["sampler"], "bound candidate fell back")
    for key in a.records:
        require(
            a.records[key].tobytes()
            == b.records[key].tobytes()
            == c.records[key].tobytes(),
            f"bound parameter replay mismatch: {key}",
        )
    require(a.observations == c.observations, "bound parameter observation mismatch")
    require(
        template.simulate().metadata["model_sha256"] != a.metadata["model_sha256"],
        "binding overwrote the template",
    )
    from crn_craft import CrnError

    for overrides in ({"unknown": 1.0}, {"mu": float("nan")}):
        try:
            template.bind_parameters(overrides)
        except CrnError:
            pass
        else:
            raise RuntimeError("invalid rules parameter accepted")
    fast = rules.ensemble(
        model,
        threads=8,
        **(
            common
            | {
                "method": {
                    "kind": "static_direct",
                    "waiting_time": "ziggurat",
                    "sparse_trackers": True,
                },
                "trackers": [
                    {"observable": "neutral_mutants"},
                    {"observable": "advantaged_mutants"},
                ],
            }
        ),
    )
    plain = fast.simulate()
    observed = fast.simulate(callback=lambda _id, _event: True)
    require(
        plain.metadata == observed.metadata,
        "generic Continue observer changed metadata",
    )
    require(
        "static-direct-sparse-trackers" in plain.metadata["sampler"],
        "generic opts missing",
    )
    for key in plain.records:
        require(
            plain.records[key].tobytes() == observed.records[key].tobytes(),
            f"generic Continue observer changed {key}",
        )
    # NumPy scalars and arrays are accepted for the grid, counts and seeds;
    # fractional counts are rejected rather than truncated.
    grid = rules.ensemble(
        model,
        end=np.float64(3.0),
        report_times=np.arange(4),
        replicates=np.int64(4),
        seed=np.uint64(7),
        threads=np.int32(2),
    )
    require(grid.simulate().metadata["termination"] == "completed", "NumPy arguments")
    try:
        rules.ensemble(model, end=3.0, report_times=[0.0, 3.0], replicates=2.5)
    except TypeError:
        pass
    else:
        raise RuntimeError("fractional replicate count accepted")
    try:
        actual.tracker_columns(1)
    except IndexError:
        pass
    else:
        raise RuntimeError("tracker index beyond the requested trackers")
    args.output.write_text(
        json.dumps(
            {
                "passed": True,
                "checks": [
                    "1/20 typed column bytes",
                    "ordered statistics",
                    "read-only buffers",
                    "event and batch delivery",
                    "event and tracker array views",
                    "continue identity",
                    "observer stop",
                    "callback exception propagation",
                    "cancel file",
                    "specialized parameter bind, replay and rejection",
                    "generic Ziggurat sparse tracker callback replay",
                    "NumPy arguments and strict counts",
                    "tracker index bounds",
                ],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
