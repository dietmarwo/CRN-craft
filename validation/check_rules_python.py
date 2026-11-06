"""Exercise the installed wheel's rules API, buffers and callback contracts."""

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from crn_craft import rules
from check_rules import require

ROOT = Path(__file__).resolve().parents[1]
RECORD_DTYPES = {
    "started": "?",
    "termination": "u1",
    "stop_time": "<f8",
    "events": "<u8",
    "rng_words": "<u8",
    "peak_instances": "<u8",
    "created_instances": "<u8",
    "saved_points": "<u8",
}


def frame_bytes(frame):
    """Index, columns, dtypes and exact value bytes of a DataFrame."""
    parts = [tuple(frame.index.names), frame.index.to_numpy().tobytes()]
    for name in frame.columns:
        column = frame[name]
        if isinstance(column.dtype, pd.CategoricalDtype):
            values = column.cat.codes.to_numpy()
        else:
            values = column.to_numpy()
        parts.append((name, str(column.dtype), values.tobytes()))
    return parts


def result_bytes(result):
    """Every DataFrame of a result, for byte-for-byte replay comparisons."""
    out = [frame_bytes(result.records), frame_bytes(result.rule_events)]
    if result.statistics is not None:
        for key, value in sorted(result.statistics.items()):
            out.append(
                (
                    key,
                    value
                    if value is None or isinstance(value, int)
                    else frame_bytes(value),
                )
            )
    return out


def check_frames(result, metadata, columns):
    """Check a Result's DataFrames against the raw native metadata/buffers."""
    meta = metadata["result"]
    observables = meta["observables"]
    times = meta["report_times"]
    ids = np.frombuffer(columns["id"], dtype="<u8")
    records = result.records
    require(
        list(records.columns) == list(RECORD_DTYPES)
        and records.index.name == "replicate"
        and records.index.dtype == np.int64
        and np.array_equal(records.index.to_numpy(), ids),
        "records index and columns",
    )
    for key, dtype in RECORD_DTYPES.items():
        raw = np.frombuffer(columns[key], dtype=dtype)
        column = records[key]
        if key == "termination":
            require(
                isinstance(column.dtype, pd.CategoricalDtype)
                and tuple(column.cat.categories) == rules.TERMINATIONS
                and np.array_equal(column.cat.codes.to_numpy(), raw)
                and list(column.astype(str)) == [rules.TERMINATIONS[i] for i in raw],
                "categorical termination",
            )
        else:
            require(
                column.dtype == np.dtype(dtype)
                and column.to_numpy().tobytes() == raw.tobytes(),
                f"record column {key}",
            )
    counts = np.frombuffer(columns["rule_events"], dtype="<u8")
    require(
        result.rule_events.shape == (len(ids), len(meta["rules"]))
        and list(result.rule_events.columns) == meta["rules"]
        and result.rule_events.index.equals(records.index)
        and set(result.rule_events.dtypes) == {np.dtype(np.uint64)}
        and result.rule_events.to_numpy().ravel().tobytes() == counts.tobytes(),
        "rule event frame",
    )
    stats = meta["statistics"]
    if stats is None:
        require(result.statistics is None, "withheld statistics")
    else:
        require(result.statistics["samples"] == stats["samples"], "statistics samples")
        for key in ("mean", "sample_variance", "mean_standard_error"):
            frame = result.statistics[key]
            if stats[key] is None:
                require(frame is None, f"absent statistic {key}")
                continue
            require(
                frame.index.name == "time"
                and list(frame.index) == times
                and list(frame.columns) == observables
                and frame.shape == (len(times), len(observables))
                and np.array_equal(frame.to_numpy(), np.array(stats[key])),
                f"statistic frame {key}",
            )
    paths = result.paths
    require(
        list(paths.index.names) == ["replicate", "time"]
        and list(paths.columns) == observables
        and len(paths) == sum(len(p["values"]) for p in meta["paths"]),
        "path frame shape",
    )
    for kept in meta["paths"]:
        rows = paths.loc[kept["id"]]
        require(
            np.array_equal(
                rows.to_numpy(), np.array(kept["values"]).reshape(-1, len(observables))
            )
            and list(rows.index) == times[: len(kept["values"])],
            f"path rows of replicate {kept['id']}",
        )
    for observation in meta.get("observations", []):
        if observation["log"]:
            frame = result.event_columns(observation["id"])
            log = observation["log"]
            require(
                list(frame.columns) == ["time", "rule", "rule_name", *observables]
                and list(frame["time"]) == [e["time"] for e in log]
                and list(frame["rule"])
                == [-1 if e["rule"] is None else e["rule"] for e in log]
                and [None if pd.isna(n) else n for n in frame["rule_name"]]
                == [
                    None if e["rule"] is None else meta["rules"][e["rule"]] for e in log
                ]
                and np.array_equal(
                    frame[observables].to_numpy(), [e["observables"] for e in log]
                ),
                "event log frame",
            )
    observed = meta.get("observations", [])
    everything = result.event_columns()
    require(
        list(everything.index.names) == ["replicate", "entry"]
        and len(everything) == sum(len(o["log"]) for o in observed)
        and all(
            everything.loc[o["id"]].equals(result.event_columns(o["id"]))
            for o in observed
            if o["log"]
        ),
        "event frame of every replicate",
    )
    for index, _ in enumerate(meta.get("observations", [{}])[0].get("trackers", [])):
        frame = result.tracker_columns(index)
        observed = meta["observations"]
        require(
            list(frame.index) == [o["id"] for o in observed]
            and all(
                frame[key].dtype == np.uint64
                for key in ("appearances", "extinctions", "omitted_episodes")
            )
            and all(
                np.array_equal(
                    frame[key].to_numpy(),
                    np.array(
                        [o["trackers"][index][key] for o in observed],
                        dtype=frame[key].dtype,
                    ),
                    equal_nan=frame[key].dtype.kind == "f",
                )
                for key in frame.columns
            ),
            "tracker frame",
        )


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
    require(result_bytes(actual) == result_bytes(expected), "thread-dependent frames")
    require(expected.observations == actual.observations, "thread-dependent log")
    require(
        np.array_equal(expected.statistics["mean"], actual.statistics["mean"]),
        "thread-dependent statistics",
    )
    require(actual.rule_events.shape == (20, 8), "rule frame shape")
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
    # Rule names: declaration order on the prepared run and the result, so a
    # callback names a rule as run.rule_names[event["rule"]].
    declared = tuple(rule["name"] for rule in json.loads(model)["rules"])
    require(
        prepared.rule_names == actual.rule_names == declared
        and tuple(actual.metadata["rules"]) == declared
        and prepared.bind_parameters({"mu": 0.3}).rule_names == declared,
        "rule names",
    )
    fired = {prepared.rule_names[e["rule"]] for log in received for e in log}
    require(fired and fired <= set(declared), "rule names of callback events")
    require(
        with_callback.metadata == actual.metadata, "observer changed numerical result"
    )
    stopped = prepared.simulate(callback=lambda _i, _event: False)
    require(
        bool((stopped.records["events"] == 1).all())
        and bool((stopped.records["termination"] == "stopped").all())
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
    require(
        result_bytes(a) == result_bytes(b) == result_bytes(c),
        "bound parameter replay mismatch",
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
    require(
        result_bytes(plain) == result_bytes(observed),
        "generic Continue observer changed the frames",
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
    # Trackers without a method run the tree under rules-run/0.2.
    tree = rules.ensemble(
        model,
        end=3.0,
        report_times=[0.0, 1.0, 3.0],
        replicates=4,
        trackers=[{"observable": "mutants"}],
    ).simulate()
    require(tree.metadata["termination"] == "completed", "trackers without method")
    try:
        actual.tracker_columns(1)
    except IndexError:
        pass
    else:
        raise RuntimeError("tracker index beyond the requested trackers")
    # DataFrames against the raw native metadata and buffers: a full run with
    # every path, a log and a first-passage tracker, and a run whose
    # replicates end on the event budget with partial paths.
    full = rules.ensemble(
        model,
        threads=2,
        keep_paths="all",
        **(
            common
            | {
                "trackers": [
                    {"observable": "mutants", "keep_episodes": 4},
                    {"observable": "N", "first_passage": "N < 7"},
                ]
            }
        ),
    )
    metadata, columns = full._native.simulate(None, None)
    metadata = json.loads(metadata)
    framed = rules.Result(metadata, columns)
    check_frames(framed, metadata, columns)
    require(
        len(framed.paths) == 20 * 3
        and framed.paths.index.get_level_values(0).nunique() == 20,
        "keep_paths='all' keeps every replicate",
    )
    require(
        framed.tracker_columns(1)["first_passage"].isna().any()
        and framed.tracker_columns(1)["first_passage"].notna().any(),
        "first passage NaN where absent",
    )
    short = rules.ensemble(
        model,
        end=3.0,
        report_times=[0.0, 1.0, 2.0, 3.0],
        replicates=6,
        seed=5,
        max_events=20,
        keep_paths=True,
        method="static_direct",
    )
    metadata, columns = short._native.simulate(None, None)
    metadata = json.loads(metadata)
    partial = rules.Result(metadata, columns)
    check_frames(partial, metadata, columns)
    saved = partial.paths.groupby(level="replicate").size()
    require(
        partial.statistics is None
        and set(partial.records["termination"]) == {"event_limit"}
        and saved.max() < 4
        and np.array_equal(
            saved.to_numpy(), partial.records["saved_points"].to_numpy()
        ),
        "partial paths of budget-ended replicates",
    )
    # Leaps of an approximate method are logged with rule -1 and no name.
    leaping = rules.ensemble(
        model,
        end=3.0,
        report_times=[0.0, 3.0],
        replicates=3,
        seed=2,
        method={"kind": "binomial", "step": 0.5},
        event_log={"max_entries": 64},
    )
    metadata, columns = leaping._native.simulate(None, None)
    metadata = json.loads(metadata)
    check_frames(rules.Result(metadata, columns), metadata, columns)
    require(
        any(
            e["rule"] is None
            for o in metadata["result"]["observations"]
            for e in o["log"]
        ),
        "leap entries in the log",
    )
    try:
        import xarray  # noqa: F401
    except ImportError:
        try:
            partial.to_xarray()
        except ImportError as error:
            require("crn-craft[xarray]" in str(error), "to_xarray hint")
        else:
            raise RuntimeError("to_xarray without xarray")
    else:
        for result in (framed, partial):
            data = result.to_xarray()
            cube = np.stack(
                [data[name].values for name in result.paths.columns], axis=-1
            )
            require(
                data[result.paths.columns[0]].dims == ("replicate", "time")
                and list(data["time"].values) == result.metadata["report_times"]
                and list(data["replicate"].values)
                == [p["id"] for p in result.metadata["paths"]]
                and data.attrs["sampler"] == result.metadata["sampler"],
                "to_xarray coordinates",
            )
            for kept in result.metadata["paths"]:
                rows = np.array(kept["values"])
                require(
                    np.array_equal(cube[kept["id"], : len(rows)], rows)
                    and np.isnan(cube[kept["id"], len(rows) :]).all(),
                    "to_xarray values",
                )
    # keep_paths: "all" and True are every replicate, False none; same request.
    shapes = {
        repr(keep): rules.ensemble(
            model, end=3.0, report_times=[0.0, 3.0], replicates=3, keep_paths=keep
        )
        .simulate()
        .metadata["run_sha256"]
        for keep in ("all", True, 3, False, 0)
    }
    require(
        shapes["'all'"] == shapes["True"] == shapes["3"]
        and shapes["False"] == shapes["0"] != shapes["3"],
        "keep_paths all, True and False",
    )
    for bad in ("some", -1, 1.5):
        try:
            rules.ensemble(model, end=3.0, report_times=[0.0, 3.0], keep_paths=bad)
        except (ValueError, TypeError):  # CrnError is a ValueError
            pass
        else:
            raise RuntimeError(f"keep_paths {bad!r} accepted")
    args.output.write_text(
        json.dumps(
            {
                "passed": True,
                "checks": [
                    "1/20 typed frame bytes",
                    "ordered statistics",
                    "frames match the native buffers and metadata",
                    "keep_paths all, True and False",
                    "partial paths and to_xarray",
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
                    "trackers without a method",
                ],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
