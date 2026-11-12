"""Exercise the installed wheel's rules API, buffers and callback contracts."""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
import xarray as xr

import crn_craft
from crn_craft import rules
from crn_craft._rules_loop import CompiledLoop, fingerprint_of
from check_rules import no_duplicates, require

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
            if key == "quantiles":
                value = [(p, frame_bytes(frame)) for p, frame in value.items()]
            elif not (value is None or isinstance(value, int)):
                value = frame_bytes(value)
            out.append((key, value))
    return out


ATTRS = ("sampler", "method", "root_seed", "model_sha256", "run_sha256")


def array_attrs(meta):
    """DataArray attrs: netCDF-storable, with the method as JSON text."""
    attrs = {key: meta.get(key) for key in ATTRS}
    attrs["method"] = json.dumps(attrs["method"], sort_keys=True)
    require(
        all(isinstance(value, str | int | float) for value in attrs.values()),
        "DataArray attrs are netCDF-storable",
    )
    return attrs


def check_path_array(result, meta):
    """Result.paths: replicates x report times x variables, NaN where unsaved."""
    observables = meta["observables"]
    times = meta["report_times"]
    kept = meta["paths"]
    paths = result.paths
    expected = np.full((len(kept), len(times), len(observables)), np.nan)
    for index, path in enumerate(kept):
        rows = np.array(path["values"], dtype=np.float64).reshape(-1, len(observables))
        expected[index, : len(rows), :] = rows
    require(
        isinstance(paths, xr.DataArray)
        and paths.dims == ("replicate", "time", "variable")
        and paths.name == "paths"
        and paths.dtype == np.float64
        and paths.shape == expected.shape
        and paths["time"].dtype == np.float64
        and paths["time"].values.tolist() == times
        and paths["variable"].values.tolist() == observables
        and paths["replicate"].dtype == np.int64
        and paths["replicate"].values.tolist() == [p["id"] for p in kept]
        and set(paths.coords) == {"time", "variable", "replicate"}
        and np.array_equal(paths.values, expected, equal_nan=True)
        and paths.attrs == array_attrs(meta),
        "path array against the raw paths",
    )
    # The 0.4.0 frame is the same data without the unsaved rows.
    frame = result.paths_frame()
    long = paths.values.reshape(-1, len(observables))
    saved = ~np.isnan(long).all(axis=1)
    require(
        np.array_equal(long[saved], frame.to_numpy())
        and list(frame.columns) == observables,
        "path array and paths_frame hold the same rows",
    )
    data = result.to_xarray()
    split = paths.to_dataset("variable")
    split.attrs = {key: meta[key] for key in ATTRS if key != "method"}
    require(
        isinstance(data, xr.Dataset) and data.identical(split),
        "to_xarray is paths as a Dataset over (replicate, time)",
    )


def named(values):
    """Rule names as nested lists, with "<missing>" for float NaN only."""
    if values.ndim > 1:
        return [named(row) for row in values]
    return [
        n if isinstance(n, str) else "<missing>" if n != n else repr(n) for n in values
    ]


def check_event_array(result, meta):
    """Result.event_columns(): replicates x entries x variables, NaN-padded."""
    observables = meta["observables"]
    rules_ = meta["rules"]
    observed = meta.get("observations", [])
    lengths = [len(o["log"]) for o in observed]
    width = max(lengths, default=0)
    shape = (len(observed), width)
    values = np.full((*shape, len(observables)), np.nan)
    time = np.full(shape, np.nan)
    rule = np.full(shape, -2, dtype=np.int64)
    name = np.full(shape, "<missing>", dtype=object)
    for row, observation in enumerate(observed):
        for entry, event in enumerate(observation["log"]):
            values[row, entry] = event["observables"]
            time[row, entry] = event["time"]
            if event["rule"] is None:
                rule[row, entry] = -1
            else:
                rule[row, entry] = event["rule"]
                name[row, entry] = rules_[event["rule"]]
    log = result.event_columns()
    require(
        isinstance(log, xr.DataArray)
        and log.dims == ("replicate", "entry", "variable")
        and log.name == "event_log"
        and log.dtype == np.float64
        and log.shape == values.shape
        and np.array_equal(log.values, values, equal_nan=True)
        and set(log.coords)
        == {"replicate", "entry", "variable", "length", "time", "rule", "rule_name"}
        and log["replicate"].dtype == np.int64
        and log["replicate"].values.tolist() == [o["id"] for o in observed]
        and log["entry"].dtype == np.int64
        and log["entry"].values.tolist() == list(range(width))
        and log["variable"].values.tolist() == observables
        and log["length"].dims == ("replicate",)
        and log["length"].dtype == np.int64
        and log["length"].values.tolist() == lengths
        and all(log[key].dims == ("replicate", "entry") for key in ("time", "rule"))
        and log["time"].dtype == np.float64
        and np.array_equal(log["time"].values, time, equal_nan=True)
        and log["rule"].dtype == np.int64
        and np.array_equal(log["rule"].values, rule)
        and log["rule_name"].dims == ("replicate", "entry")
        and log["rule_name"].dtype == object
        and named(log["rule_name"].values) == name.tolist()
        and log.attrs == array_attrs(meta),
        "event array against the raw logs",
    )
    for row, observation in enumerate(observed):
        one = result.event_columns(observation["id"])
        require(
            one.dims == ("entry", "variable")
            and one.shape == (lengths[row], len(observables))
            and one.identical(
                log.sel(replicate=observation["id"]).isel(entry=slice(0, lengths[row]))
            ),
            "one replicate's event array is its unpadded row",
        )
    try:
        result.event_columns(-1)
    except KeyError:
        pass
    else:
        raise RuntimeError("event_columns accepted an unobserved replicate")


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
    check_path_array(result, meta)
    check_event_array(result, meta)
    paths = result.paths_frame()
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
            frame = result.event_frame(observation["id"])
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
    everything = result.event_frame()
    require(
        list(everything.index.names) == ["replicate", "entry"]
        and len(everything) == sum(len(o["log"]) for o in observed)
        and all(
            everything.loc[o["id"]].equals(result.event_frame(o["id"]))
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


# Result.save and rules.load_report: the report layout of crn-craft.report/1.0.
REPORT_KEYS = ["schema", "crn_craft_version", "model", "request", "result"]
REPORT_RESULT_KEYS = [
    "model_sha256",
    "run_sha256",
    "sampler",
    "report_times",
    "observables",
    "terminations",
    "rule_events_mean",
    "statistics",
    "paths",
]
REPORT_STATISTICS = ("mean", "mean_standard_error", "sample_variance")
NON_FINITE = ("nan", "inf", "-inf")


def report_floats(value):
    """Report numbers as floats; non-finite ones are written as strings."""
    if isinstance(value, list):
        return [report_floats(v) for v in value]
    require(
        value in NON_FINITE
        or (
            isinstance(value, int | float)
            and not isinstance(value, bool)
            and math.isfinite(value)
        ),
        f"report number {value!r}",
    )
    return float(value)


def same_values(written, expected):
    """Written report numbers equal `expected` exactly, NaN where NaN."""
    return np.array_equal(
        np.array(report_floats(written), dtype=np.float64),
        np.asarray(expected, dtype=np.float64),
        equal_nan=True,
    )


def read_report(path):
    """A JSON report parsed strictly, and the same through load_report."""
    text = Path(path).read_text("utf-8")

    def constant(name):
        raise RuntimeError(f"non-standard JSON constant {name} in a report")

    data = json.loads(text, parse_constant=constant, object_pairs_hook=no_duplicates)
    require(
        text.endswith("}\n") and not text.endswith("\n\n") and "\r" not in text,
        "a JSON report ends with one newline",
    )
    require(rules.load_report(path) == data, "load_report returns the document")
    return data


def termination_counts(result):
    """Counts of the terminations that occurred, in TERMINATIONS order."""
    counts = result.records["termination"].value_counts()
    return {name: int(counts[name]) for name in rules.TERMINATIONS if counts[name]}


def mean_events(result):
    """Mean event count per rule, from the exact integer sums."""
    events = result.rule_events
    return {name: int(events[name].sum()) / len(events) for name in result.rule_names}


def report_statistics(result):
    """(key, title, DataFrame) of every statistic a report shows, in order."""
    stats = result.statistics
    out = [
        (key, f"### {key}", stats[key])
        for key in REPORT_STATISTICS
        if stats.get(key) is not None
    ]
    for p, frame in sorted((stats.get("quantiles") or {}).items()):
        # The percentage with the shortest exact decimal: 2.5, 97.5, 25.
        text = format(float(rules_decimal(p) * 100), "f")
        if "." in text:
            text = text.rstrip("0").rstrip(".")
        title = "### median" if p == 0.5 else f"### {text}% percentile"
        out.append((p, title, frame))
    return out


def rules_decimal(p):
    """The decimal that the shortest repr of p denotes."""
    from decimal import Decimal

    return Decimal(repr(float(p)))


def check_report_json(result, report, *, paths=True):
    """Every number of a JSON report against the Result it was saved from."""
    meta = result.metadata
    observables = meta["observables"]
    require(
        list(report) == REPORT_KEYS
        and report["schema"] == "crn-craft.report/1.0"
        and report["crn_craft_version"] == crn_craft.__version__
        and report["model"] == result.model
        and report["request"] == result.request,
        "report header, model and request",
    )
    body = report["result"]
    require(
        list(body) == REPORT_RESULT_KEYS
        and all(body[k] == meta[k] for k in ("model_sha256", "run_sha256", "sampler"))
        and same_values(body["report_times"], meta["report_times"])
        and body["observables"] == observables,
        "report provenance, grid and observables (result order)",
    )
    expected = termination_counts(result)
    events = mean_events(result)
    require(
        body["terminations"] == expected
        and list(body["terminations"]) == list(expected)
        and all(isinstance(v, int) for v in body["terminations"].values())
        and list(body["rule_events_mean"]) == list(result.rule_names)
        and all(
            same_values(body["rule_events_mean"][name], events[name])
            for name in result.rule_names
        ),
        "report terminations and mean rule events",
    )
    stats = result.statistics
    if stats is None:
        require(body["statistics"] is None, "withheld statistics are null")
    else:
        written = body["statistics"]
        quantiles = stats.get("quantiles") or {}
        require(
            list(written) == ["samples", *REPORT_STATISTICS, "quantiles"]
            and written["samples"] == stats["samples"]
            and [q["p"] for q in written["quantiles"]] == sorted(quantiles)
            and all(list(q) == ["p", "values"] for q in written["quantiles"]),
            "report statistics keys and quantile order",
        )
        tables = [(written[key], stats.get(key)) for key in REPORT_STATISTICS]
        tables += [(q["values"], quantiles[q["p"]]) for q in written["quantiles"]]
        for table, frame in tables:
            if frame is None:
                require(table is None, "a statistic not computed is null")
                continue
            require(
                list(table) == observables
                and all(
                    same_values(table[name], frame[name].to_numpy())
                    for name in observables
                ),
                "report statistics equal the DataFrames",
            )
    if paths and meta["paths"]:
        cube = result.paths
        kept = body["paths"]
        require(
            list(kept) == ["replicates", "values"]
            and kept["replicates"] == cube["replicate"].values.tolist()
            and list(kept["values"]) == observables
            and all(
                same_values(kept["values"][name], cube.sel(variable=name).values)
                for name in observables
            ),
            "report paths equal the kept paths",
        )
    else:
        require(body["paths"] is None, "a report without paths has paths null")


def markdown_sections(text):
    """(heading, data lines) of a markdown report, in file order.

    Checks the rendered layout on the way: a blank line before every heading
    but the first, after every heading, and around every fence; the data
    lines of a section in one ```text fence; no other text than the
    `Name:` paragraph of `# Model`.
    """
    lines = text.split("\n")
    require(lines[-1] == "" and lines[-2] != "", "one trailing newline")
    lines = lines[:-1]
    sections, fenced, index = [], set(), 0
    while index < len(lines):
        line = lines[index]
        if line.startswith("#"):
            require(
                (index == 0 or lines[index - 1] == "")
                and index + 1 < len(lines)
                and lines[index + 1] == "",
                f"blank lines around the heading {line!r}",
            )
            sections.append((line, []))
        elif line == "```text":
            heading = sections[-1][0]
            end = lines.index("```", index + 1)
            require(
                lines[index - 1] == ""
                and (end + 1 == len(lines) or lines[end + 1] == "")
                and heading not in fenced
                and end > index + 1
                and not any(x.startswith("```") for x in lines[index + 1 : end]),
                f"one fenced block of data lines under {heading!r}",
            )
            fenced.add(heading)
            sections[-1][1].extend(lines[index + 1 : end])
            index = end
        elif line:
            require(
                line.startswith("Name: ") and sections[-1] == ("# Model", []),
                f"markdown text outside a fence: {line[:40]!r}",
            )
            sections[-1][1].append(line)
        index += 1
    require(
        {heading for heading, data in sections if data and heading != "# Model"}
        == fenced,
        "every section with data lines is fenced",
    )
    return sections


def split_values(line, label):
    """The numbers of a `<label>: a, b, ...` line."""
    head, _, values = line.partition(": ")
    require(head == label, f"markdown line {line[:40]!r} is labeled {label}")
    return [v if v in NON_FINITE else float(v) for v in values.split(", ")]


def check_report_markdown(result, text, order, *, paths=True):
    """A markdown report's simulation and results sections against a Result.

    `order` is the observable order the report must use.
    """
    meta = result.metadata
    require(
        sorted(order) == sorted(meta["observables"])
        and text.startswith("# Model\n\nName: ")
        and text.endswith("\n")
        and not text.endswith("\n\n")
        and "\r" not in text,
        "markdown report frame",
    )
    sections = markdown_sections(text)
    headings = [heading for heading, _ in sections]
    found = dict(sections)
    require(
        headings.index("# Model")
        < headings.index("# Simulation")
        < headings.index("# Results")
        < headings.index("## Terminations")
        < headings.index("## rule_events (mean per replicate)")
        < headings.index("## Summary statistics"),
        "markdown report section order",
    )
    simulation = dict(line.split(" = ", 1) for line in found["# Simulation"])
    request = result.request
    require(
        list(simulation)[:5] == ["seed", "end", "replicates", "report_times", "method"]
        and list(simulation)[-4:]
        == ["sampler", "crn-craft", "model_sha256", "run_sha256"]
        and int(simulation["seed"]) == request["ensemble"]["root_seed"]
        and float(simulation["end"]) == request["time"]["end"]
        and int(simulation["replicates"]) == request["ensemble"]["replicates"]
        and [float(t) for t in simulation["report_times"].split(", ")]
        == meta["report_times"]
        and simulation["sampler"] == meta["sampler"]
        and simulation["crn-craft"] == crn_craft.__version__
        and simulation["model_sha256"] == meta["model_sha256"]
        and simulation["run_sha256"] == meta["run_sha256"],
        "markdown simulation section",
    )
    require(
        found["## Terminations"]
        == [f"{k} = {v}" for k, v in termination_counts(result).items()],
        "markdown terminations",
    )
    events = mean_events(result)
    lines = found["## rule_events (mean per replicate)"]
    require(
        [line.split(" = ")[0] for line in lines] == list(result.rule_names)
        and all(
            float(line.split(" = ")[1]) == events[name]
            for line, name in zip(lines, result.rule_names, strict=True)
        ),
        "markdown mean rule events",
    )
    heading = "## Single trajectories for " + ", ".join(order)
    if paths and meta["paths"]:
        cube = result.paths
        ids = cube["replicate"].values.tolist()
        width = max(3, len(str(max(ids))))
        labels = [f"{name}{i:0{width}d}" for i in ids for name in order]
        lines = found[heading]
        require(
            len(lines) == len(labels)
            and all(
                same_values(
                    split_values(line, label),
                    cube.sel(replicate=int(label[len(name) :]), variable=name).values,
                )
                for line, label, name in zip(
                    lines, labels, [n for _ in ids for n in order], strict=True
                )
            ),
            "markdown trajectories equal the kept paths",
        )
    else:
        require(
            not any(h.startswith("## Single trajectories") for h in headings),
            "no trajectories without kept paths",
        )
    position = headings.index("## Summary statistics")
    if result.statistics is None:
        reason = found["## Summary statistics"]
        require(
            len(reason) == 1
            and reason[0].startswith("statistics withheld (")
            and reason[0].endswith(")")
            and meta["termination"] in reason[0]
            and headings[position + 1 :] == [],
            "markdown withheld statistics line names the termination",
        )
        return
    shown = report_statistics(result)
    require(
        found["## Summary statistics"] == [f"samples = {result.statistics['samples']}"]
        and headings[position + 1 :] == [title for _, title, _ in shown],
        "markdown statistics subsections and their order",
    )
    for _, title, frame in shown:
        lines = found[title]
        require(
            len(lines) == len(order)
            and all(
                same_values(split_values(line, name), frame[name].to_numpy())
                for line, name in zip(lines, order, strict=True)
            ),
            f"markdown {title[4:]} equals the DataFrame",
        )


def saved(result, folder, name, **kwargs):
    """Save `result` under `folder` and return (returned path, text)."""
    path = result.save(Path(folder) / name, **kwargs)
    return path, path.read_text("utf-8")


def sir_model():
    """The SIR builder model of the report example."""
    sir = rules.Model("sir")
    sir.add_species({"S": 999, "I": 1, "R": 0})
    sir.add_parameter({"beta": 2e-4, "gamma": 0.01})
    sir.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
    sir.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    sir.set_stop("I == 0")
    return sir


def check_result_documents(model):
    """result.model and result.request are the documents that ran."""
    sir = sir_model()
    run = rules.ensemble(
        sir, end=50.0, report_times=[0.0, 50.0], replicates=4, seed=3, threads=2
    )
    result = run.simulate()
    require(
        result.model == sir.to_dict()
        and result.request == json.loads(run._documents[1])
        and result.request["model"]["sha256"] == result.metadata["model_sha256"]
        and result.model is not result.model,
        "result.model and result.request are new dictionaries of what ran",
    )
    # A hand-built Result knows neither document and cannot be saved.
    metadata, columns = run._native.simulate(None, None)
    bare = rules.Result(json.loads(metadata), columns)
    require(bare.model is None and bare.request is None, "a bare Result")
    with tempfile.TemporaryDirectory() as folder:
        try:
            bare.save(Path(folder) / "bare.md")
        except ValueError:
            pass
        else:
            raise RuntimeError("a Result without its documents was saved")
        require(not list(Path(folder).iterdir()), "a refused save writes nothing")
    # bind_parameters: the bound values and the candidate's model hash, and
    # the documents replay the candidate's numbers.
    fixed = {
        "end": 3.0,
        "report_times": [0.0, 1.0, 3.0],
        "replicates": 6,
        "seed": 1234,
        "method": "static_direct",
        "keep_paths": "all",
    }
    template = rules.ensemble(model, **fixed)
    candidate = template.bind_parameters({"mu": 0.3, "sm": 1.8})
    bound = candidate.simulate()
    again = candidate.bind_parameters({"mu": 0.2}).simulate()
    expected = json.loads(model)
    require(
        template.simulate().model == expected
        and bound.model["parameters"]
        == {
            k: float(v)
            for k, v in sorted({**expected["parameters"], "mu": 0.3, "sm": 1.8}.items())
        }
        and again.model["parameters"]["mu"] == 0.2
        and again.model["parameters"]["sm"] == 1.8
        and {k: v for k, v in bound.model.items() if k != "parameters"}
        == {k: v for k, v in expected.items() if k != "parameters"}
        and bound.request["model"]["sha256"] == bound.metadata["model_sha256"]
        and again.request["model"]["sha256"] == again.metadata["model_sha256"]
        and bound.metadata["model_sha256"] != again.metadata["model_sha256"],
        "a bound candidate reports its bound parameter values and model hash",
    )
    for result in (bound, again):
        request = result.request
        request["model"]["sha256"] = hashlib.sha256(
            rules.model_bytes(result.model)
        ).hexdigest()
        replay = rules.prepare(result.model, request).simulate()
        require(
            replay.metadata["sampler"] == result.metadata["sampler"]
            and result_bytes(replay) == result_bytes(result)
            and replay.paths.equals(result.paths),
            "a bound candidate's documents replay its numbers",
        )
    # The compiled path, with a stand-in loop that runs the interpreter: its
    # results know the documents of the run they came from.
    loop = CompiledLoop(
        fingerprint=fingerprint_of(template._native.loop_source()),
        model_sha256=hashlib.sha256(model).hexdigest(),
        version=crn_craft._native.__version__,
        directory=tempfile.gettempdir(),
        cached=True,
    )
    loop.execute = lambda request, cancel_file=None: template._native.simulate(
        None, None
    )
    compiled = template.simulate(compiled=loop)
    verified = loop.verify(template)
    require(
        compiled.model == verified.model == expected
        and compiled.request == verified.request == json.loads(template._documents[1]),
        "results of the compiled path know their documents",
    )


QUANTILES = (0.025, 0.25, 0.5, 0.75, 0.975)


def check_numpy_quantiles(result, probabilities, label):
    """statistics["quantiles"][p] is np.quantile of the kept paths, bit for
    bit, shaped like statistics["mean"]; "median" is quantiles[0.5]."""
    stats = result.statistics
    require(stats is not None, f"{label}: statistics present")
    mean = stats["mean"]
    quantiles = stats["quantiles"]
    require(
        list(quantiles) == list(probabilities)
        and all(isinstance(p, float) for p in quantiles),
        f"{label}: one float key per probability, in request order",
    )
    cube = result.paths.values
    require(
        cube.shape[0] == result.metadata["requested_replicates"],
        f"{label}: every path is kept",
    )
    for p, frame in quantiles.items():
        expected = np.quantile(cube, p, axis=0)
        require(
            isinstance(frame, pd.DataFrame)
            and frame.index.equals(mean.index)
            and list(frame.columns) == list(mean.columns)
            and all(str(t) == "float64" for t in frame.dtypes)
            and frame.to_numpy().shape == expected.shape
            and frame.to_numpy().tobytes() == expected.tobytes(),
            f"{label}: quantiles[{p}] equals numpy.quantile of the paths",
        )
    if 0.5 in quantiles:
        require(stats["median"] is quantiles[0.5], f"{label}: median")
    else:
        require("median" not in stats, f"{label}: no median without 0.5")


def without_quantiles(result):
    """result_bytes without the quantile entries."""
    return [
        entry
        for entry in result_bytes(result)
        if not (isinstance(entry, tuple) and entry[0] in ("quantiles", "median"))
    ]


def refused_with(action, expected, label):
    """`action` raises CrnError with the (code, path) `expected`; returns it."""
    try:
        action()
    except crn_craft.CrnError as error:
        first = error.diagnostics[0]
        require(
            (first["code"], first["path"]) == expected,
            f"{label}: expected {expected}, got {(first['code'], first['path'])}",
        )
        return first
    raise RuntimeError(f"{label}: accepted")


def check_fsp():
    """rules.fsp: the Poisson law of immigration-death and the refusals."""
    m = rules.Model("immigration_death")
    m.add_species({"X": 0})
    m.add_reaction(10.0, [], ["X"])
    m.add_reaction(1.0, ["X"], [])
    result = rules.fsp(m, 2.0, [0.0, 2.0], {"X": 80})
    require(
        result["schema"] == "crn-craft.fsp-result/1.0" and result["species"] == ["X"],
        "fsp result contract",
    )
    require(result["snapshots"][0]["probabilities"][0] == 1.0, "fsp initial state")
    mean = 10.0 * -math.expm1(-2.0)
    worst = max(
        abs(p - math.exp(-mean) * mean ** s[0] / math.factorial(s[0]))
        for s, p in zip(
            result["states"], result["snapshots"][1]["probabilities"], strict=True
        )
    )
    require(worst < 1e-12, f"fsp Poisson law, largest error {worst}")
    first = result["snapshots"][1]["projected_first_moment"][0]
    require(abs(first - mean) < 1e-9, "fsp mean")
    refused_with(
        lambda: rules.fsp(m, 2.0, [2.0], {"X": 80}, max_states=10),
        ("CRN-E070", "/fsp"),
        "fsp state cap",
    )
    refused_with(
        lambda: rules.fsp(m, 2.0, [2.0], {}),
        ("CRN-R020", "/method/bounds/X"),
        "fsp missing bound",
    )
    refused_with(
        lambda: rules.fsp(m, 2.0, [2.0], {"X": 80}, tolerance=1.0),
        ("CRN-R020", "/method/tolerance"),
        "fsp tolerance",
    )
    refused_with(
        lambda: rules.ensemble(m, end=2.0, report_times=[2.0], method="fsp"),
        ("CRN-R020", "/method"),
        "ensemble with the fsp method",
    )
    for bounds in ([80], {"X": 1.5}):
        try:
            rules.fsp(m, 2.0, [2.0], bounds)
        except TypeError:
            pass
        else:
            raise RuntimeError(f"fsp bounds {bounds!r} accepted")


def check_archives():
    """rules.archive, PreparedEnsemble.archive and rules.verify_archive."""
    sir = sir_model()
    run = rules.ensemble(
        sir,
        end=50.0,
        report_times=[0.0, 25.0, 50.0],
        replicates=6,
        seed=5,
        threads=2,
        keep_paths="all",
        quantiles=(0.5,),
    )
    result = run.simulate()
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        report = run.archive(folder / "run")
        require(
            report["schema"] == "crn-craft.archive-check/1.0"
            and report["valid"]
            and report["simulation_complete"]
            and report["termination"] == "completed"
            and report["result"] == "crn-craft.result/1.0"
            and report["model_sha256"] == result.metadata["model_sha256"]
            and report["run_sha256"] == result.metadata["run_sha256"],
            f"prepared run archive: {report}",
        )
        model, request = run._documents
        require(
            (folder / "run/model.json").read_bytes() == model
            and (folder / "run/request.json").read_bytes() == request,
            "archived input bytes",
        )
        archived = json.loads((folder / "run/result.json").read_bytes())
        require(
            archived["statistics"]["mean"] == result.metadata["statistics"]["mean"]
            and archived["statistics"]["quantiles"]
            == result.metadata["statistics"]["quantiles"]
            and [r["events"] for r in archived["replicates"]]
            == result.records["events"].tolist(),
            "archived result equals the simulated result",
        )
        verified = rules.verify_archive(folder / "run")
        require(
            verified["valid"] and verified["run_sha256"] == report["run_sha256"],
            "verify a prepared run archive",
        )
        # The same archive from the documents, and refusals as reports.
        pretty = json.dumps(json.loads(model), indent=1)
        same = rules.archive(pretty, request, folder / "documents")
        require(
            not same["valid"]
            and (same["diagnostics"][0]["code"], same["diagnostics"][0]["path"])
            == ("CRN-R020", "/model/sha256")
            and not (folder / "documents").exists(),
            "a re-serialized model has another hash",
        )
        same = rules.archive(model, request, folder / "documents")
        require(
            same == report
            and (folder / "documents/result.json").read_bytes()
            == (folder / "run/result.json").read_bytes(),
            "rules.archive equals PreparedEnsemble.archive",
        )
        again = run.archive(folder / "run")
        require(
            not again["valid"] and again["diagnostics"][0]["code"] == "CRN-E065",
            "an existing directory is refused",
        )
        # Tampering: a byte appended to the result.
        with (folder / "run/result.json").open("ab") as handle:
            handle.write(b" ")
        tampered = rules.verify_archive(folder / "run")
        require(
            not tampered["valid"] and tampered["diagnostics"][0]["code"] == "CRN-E065",
            "tampered archive accepted",
        )
        # A finite state projection request.
        decay = rules.Model("decay")
        decay.add_species({"X": 5})
        decay.add_reaction(1.0, ["X"], [])
        data = rules.model_bytes(decay.to_dict())
        fsp_request = {
            "schema": rules.RUN_1_0,
            "model": {"path": "model.json", "sha256": hashlib.sha256(data).hexdigest()},
            "time": {"end": 1.0, "report_times": [0.0, 1.0]},
            "method": {
                "kind": "fsp",
                "bounds": {"X": 5},
                "max_states": 10,
                "tolerance": 1e-12,
            },
            "budget": {"max_output_bytes": 100000},
        }
        projected = rules.archive(data, fsp_request, folder / "fsp")
        require(
            projected["valid"]
            and projected["result"] == "crn-craft.fsp-result/1.0"
            and projected["termination"] is None
            and rules.verify_archive(folder / "fsp")["valid"],
            "finite state projection archive",
        )
        # A cancelled run is a valid archive of an incomplete run.
        cancel = folder / "cancel"
        cancel.write_text("")
        cancelled = run.archive(folder / "cancelled", cancel_file=cancel)
        require(
            cancelled["valid"]
            and not cancelled["simulation_complete"]
            and cancelled["termination"] == "cancelled"
            and rules.verify_archive(folder / "cancelled")["valid"],
            "cancelled archive",
        )
        # The network engine's archive/0.1 is no longer read (M8b): its file
        # set, with a manifest and marker that match, is refused.
        old = folder / "old"
        old.mkdir()
        files = []
        for name in ("model.json", "run.json", "result.json", "provenance.json"):
            (old / name).write_bytes(b"{}\n")
            files.append(
                {
                    "path": name,
                    "bytes": 3,
                    "sha256": hashlib.sha256(b"{}\n").hexdigest(),
                }
            )
        manifest = json.dumps(
            {"schema": "crn-craft.archive/0.1", "files": files}
        ).encode()
        (old / "manifest.json").write_bytes(manifest)
        (old / "COMPLETE").write_text(hashlib.sha256(manifest).hexdigest() + "\n")
        checked = rules.verify_archive(old)
        require(
            not checked["valid"]
            and checked["schema"] == "crn-craft.archive-check/1.0"
            and checked["diagnostics"][0]["code"] == "CRN-E065",
            "a network archive/0.1 is refused",
        )
    try:
        run.bind_parameters({"beta": 3e-4}).archive("unused")
    except ValueError:
        pass
    else:
        raise RuntimeError("a bind_parameters candidate archived without its bytes")


def check_quantiles(model):
    """quantiles= on ensemble and prepare: engine quantiles equal
    numpy.quantile of the kept paths; the documents are the 1.0 contracts,
    rewriting only the schema strings; 0.6.0 documents without them."""
    sir = sir_model()
    settings = {
        "end": 250.0,
        "report_times": [0, 50, 100, 150, 200, 250],
        "replicates": 40,
        "seed": 1,
        "keep_paths": "all",
        "method": "static_direct",
    }
    plain_run = rules.ensemble(sir, **settings)
    run = rules.ensemble(sir, quantiles=QUANTILES, **settings)
    # Without quantiles the documents are those of 0.6.0, byte for byte.
    for none in (None, (), [], np.array([])):
        require(
            rules.ensemble(sir, quantiles=none, **settings)._documents
            == plain_run._documents,
            f"quantiles={none!r} writes the documents without quantiles",
        )
    old_model, old_request = plain_run._documents
    require(
        json.loads(old_request)["schema"] == "crn-craft.rules-run/0.6"
        and json.loads(old_model)["schema"] == "crn-craft.rules-model/0.3",
        "the SIR documents without quantiles",
    )
    data, request = run._documents
    expected_request = json.loads(old_request)
    expected_request["schema"] = rules.RUN_1_0
    expected_request["model"]["sha256"] = hashlib.sha256(data).hexdigest()
    expected_request["outputs"]["quantiles"] = list(QUANTILES)
    # run/1.0 defaults to "auto"; the linear scan of rules-run is written out.
    expected_request["method"]["selector"] = "linear"
    require(
        data
        == old_model.replace(
            b'"crn-craft.rules-model/0.3"', b'"crn-craft.model/1.0"', 1
        )
        and json.loads(request) == expected_request,
        "quantiles write run/1.0 and model/1.0, changing only schema strings "
        "and the selector",
    )
    plain, result = plain_run.simulate(), run.simulate()
    require(
        result.metadata["schema"] == "crn-craft.result/1.0"
        and result.model == {**sir.to_dict(), "schema": rules.MODEL_1_0}
        and result.request == json.loads(request)
        and result.metadata["sampler"] == plain.metadata["sampler"]
        and without_quantiles(result) == result_bytes(plain)
        and result.paths.equals(plain.paths),
        "quantiles change the schemas and hashes, not the numbers",
    )
    check_numpy_quantiles(result, QUANTILES, "SIR")
    # Widely spread values: x[j] + g (x[j+1] - x[j]) as written rounds
    # differently from numpy.quantile in some cells; the engine is NumPy's.
    spread = (0.01, 0.03, 0.07, 0.12, 0.2, 0.3, 0.38, 0.45, 0.55, 0.62)
    spread += (0.7, 0.8, 0.88, 0.93, 0.97, 0.99)
    epidemics = rules.ensemble(
        {**sir.to_dict(), "schema": rules.MODEL_1_0},
        end=100.0,
        report_times=[0.0, 20.0, 30.0, 40.0, 50.0, 100.0],
        replicates=10,
        seed=17,
        keep_paths="all",
        method="static_direct",
        quantiles=spread,
    ).simulate()
    check_numpy_quantiles(epidemics, spread, "SIR, ten epidemics")
    ordered = np.sort(epidemics.paths.values, axis=0)
    as_written = 0
    for p in spread:
        h = 9 * p
        j = int(np.floor(h))
        written = ordered[j] + (h - j) * (ordered[j + 1] - ordered[j])
        as_written += int(
            (written != np.quantile(epidemics.paths.values, p, axis=0)).sum()
        )
    require(as_written > 0, "the ten epidemics tell the two evaluations apart")
    # Every method, on models of its scope; bytes from rules.load keep every
    # byte but the schema string.
    data = model
    loaded = rules.ensemble(
        data, end=3.0, report_times=[0.0, 3.0], quantiles=(0.5,)
    )._documents[0]
    require(
        loaded
        == data.replace(b'"crn-craft.rules-model/0.1"', b'"crn-craft.model/1.0"'),
        "model bytes keep every byte but the schema string",
    )
    common = {
        "end": 3.0,
        "report_times": [0.0, 1.0, 3.0],
        "replicates": 30,
        "seed": 9,
        "keep_paths": "all",
        "quantiles": QUANTILES,
    }
    for method in (
        None,
        "static_direct",
        {"kind": "static_direct", "waiting_time": "inverse"},
        "next_reaction",
        {"kind": "tau", "epsilon": 0.03},
        "mtdna_direct",
    ):
        check_numpy_quantiles(
            rules.ensemble(model, method=method, **common).simulate(),
            QUANTILES,
            f"method {method}",
        )
    # crn_craft.rules keeps the linear scan of every rules-run version when
    # quantiles write run/1.0, so its paths equal those of the call without
    # quantiles even where "auto" would pick another selector (a ring of
    # eight first-order conversions with ten molecules each: rssa_cr); the
    # top-level crn_craft.ensemble leaves the selector to "auto".
    ring = {
        "schema": "crn-craft.rules-model/0.1",
        "name": "ring8",
        "parameters": {"k": 1.0},
        "types": [{"name": f"A{i}"} for i in range(8)],
        "initial": [{"type": f"A{i}", "count": 10} for i in range(8)],
        "observables": {f"A{i}": {"type": f"A{i}"} for i in range(8)},
        "rules": [
            {
                "name": f"r{i}",
                "reactants": [{"label": "x", "type": f"A{i}"}],
                "rate": "k",
                "products": [{"type": f"A{(i + 1) % 8}"}],
            }
            for i in range(8)
        ],
    }
    ring_common = {**common, "end": 5.0, "report_times": [0.0, 2.0, 5.0]}
    without_q = {k: v for k, v in ring_common.items() if k != "quantiles"}
    with_quantiles = rules.ensemble(ring, method="static_direct", **ring_common)
    plain = rules.ensemble(ring, method="static_direct", **without_q).simulate()
    top_level = crn_craft.ensemble(ring, method="static_direct", **ring_common)
    auto_pick = rules.ensemble(
        ring, method="static_direct", selector="auto", **ring_common
    ).simulate()
    require(
        json.loads(with_quantiles._documents[1])["method"].get("selector") == "linear"
        and with_quantiles.simulate().paths.equals(plain.paths)
        and "selector" not in json.loads(top_level._documents[1])["method"]
        and auto_pick.metadata["sampler"] != plain.metadata["sampler"],
        "rules.ensemble(quantiles=) keeps the linear scan and its paths; "
        "crn_craft.ensemble runs auto",
    )
    birth_death = {
        "schema": "crn-craft.rules-model/0.1",
        "name": "birth_death",
        "parameters": {"b": 0.9, "d": 1.0},
        "types": [{"name": "X"}],
        "initial": [{"type": "X", "count": 30}],
        "observables": {"X": {"type": "X"}},
        "rules": [
            {
                "name": "birth",
                "reactants": [{"label": "x", "type": "X"}],
                "rate": "b",
                "products": [{"type": "X"}, {"type": "X"}],
            },
            {
                "name": "death",
                "reactants": [{"label": "x", "type": "X"}],
                "rate": "d",
                "products": [],
            },
        ],
    }
    check_numpy_quantiles(
        rules.ensemble(
            birth_death, method={"kind": "binomial", "step": 0.01}, **common
        ).simulate(),
        QUANTILES,
        "binomial",
    )
    atp = rules.load(ROOT / "models/rules/mtdna-three-class-atp.rules.json")
    thinning = rules.ensemble(
        atp, method="thinning", **{**common, "end": 200.0, "report_times": [0, 200]}
    ).simulate()
    require("ATP" in thinning.statistics["mean"].columns, "thinning reports ATP")
    check_numpy_quantiles(thinning, QUANTILES, "thinning")
    ramp = rules.Model("ramp")
    ramp.add_species("A")
    ramp.add_ode("x", 0.0, "2", lower_bound=0, on_bound="clamp")
    ramp.add_reaction("x", [], ["A"], name="immigration")
    hybrid = rules.ensemble(ramp, method="hybrid", **common).simulate()
    require(
        list(hybrid.statistics["mean"].columns) == ["A", "x"]
        and hybrid.model["schema"] == rules.MODEL_1_0,
        "hybrid: the continuous column, model/1.0",
    )
    check_numpy_quantiles(hybrid, QUANTILES, "hybrid")
    # Thread counts do not change any frame.
    threads = {
        n: result_bytes(
            rules.ensemble(
                model, threads=n, method="static_direct", **common
            ).simulate()
        )
        for n in (1, 2, 4)
    }
    require(threads[1] == threads[2] == threads[4], "quantiles for 1, 2, 4 threads")
    # Only quantiles: the other statistics are None.
    only = rules.ensemble(sir, quantiles=(0.5,), **settings)
    document = json.loads(only._documents[1])
    document["outputs"]["statistics"] = []
    only = rules.prepare(only._documents[0], document).simulate()
    require(
        only.statistics["mean"] is None
        and only.statistics["samples"] == 40
        and list(only.statistics["quantiles"]) == [0.5],
        "quantiles without other statistics",
    )
    # Withheld with the other statistics.
    limited = rules.ensemble(
        model, max_events=20, method="static_direct", **common
    ).simulate()
    require(limited.statistics is None, "quantiles are withheld")
    # prepare(): a rules-run/0.4 Static Direct request keeps its meaning
    # (inverse waits, per-event trackers, and since M4c's selectors the
    # linear scan, which run/1.0 must name), as bytes, text or a dictionary.
    base = json.loads(
        rules.ensemble(model, method="static_direct", **common)._documents[1]
    )
    base["schema"] = "crn-craft.rules-run/0.4"
    base["model"]["sha256"] = hashlib.sha256(model).hexdigest()
    del base["outputs"]["quantiles"]
    del base["method"]["selector"]
    reference = rules.prepare(model, base).simulate()
    for given in (base, json.dumps(base), json.dumps(base).encode()):
        prepared = rules.prepare(model, given, quantiles=QUANTILES).simulate()
        require(
            prepared.request["schema"] == rules.RUN_1_0
            and prepared.request["method"]
            == {
                "kind": "static_direct",
                "waiting_time": "inverse",
                "sparse_trackers": False,
                "selector": "linear",
            }
            and prepared.request["model"]["sha256"] == prepared.metadata["model_sha256"]
            and prepared.model["schema"] == rules.MODEL_1_0
            and prepared.metadata["sampler"] == reference.metadata["sampler"]
            and without_quantiles(prepared) == result_bytes(reference)
            and prepared.paths.equals(reference.paths),
            f"prepare(quantiles=) keeps a rules-run/0.4 request's paths "
            f"({type(given).__name__})",
        )
        check_numpy_quantiles(prepared, QUANTILES, "prepare")
    require(base["schema"] == "crn-craft.rules-run/0.4", "the caller's request")
    # A run/1.0 request may carry the field itself, with a rules-model.
    unified = {**base, "schema": rules.RUN_1_0}
    unified["outputs"] = {**base["outputs"], "quantiles": [0.5]}
    carried = rules.prepare(model, unified)
    require(
        carried._documents[0] == model
        and list(carried.simulate().statistics["quantiles"]) == [0.5],
        "a run/1.0 request with outputs.quantiles and a rules-model",
    )
    try:
        rules.prepare(model, unified, quantiles=(0.5,))
    except ValueError as error:
        require("not both" in str(error), "quantiles given twice")
    else:
        raise RuntimeError("quantiles given twice were accepted")
    duplicate = json.dumps(base)[:-1] + ', "schema": "crn-craft.rules-run/0.4"}'
    refused_with(
        lambda: rules.prepare(model, duplicate, quantiles=(0.5,)),
        ("CRN-E002", ""),
        "a duplicate key in the request",
    )
    # bind_parameters candidates keep the quantiles and the 1.0 documents.
    template = rules.ensemble(model, method="static_direct", **common)
    candidate = template.bind_parameters({"mu": 0.3}).simulate()
    require(
        candidate.request["outputs"]["quantiles"] == list(QUANTILES)
        and candidate.model["schema"] == rules.MODEL_1_0
        and candidate.metadata["schema"] == "crn-craft.result/1.0",
        "a bound candidate keeps the quantiles",
    )
    check_numpy_quantiles(candidate, QUANTILES, "bound candidate")
    # Refusals: the engine's, with codes, and Python's type checks.
    for given, pointer in (
        ((0.5, 0.5), "/outputs/quantiles/1"),
        ((0.0,), "/outputs/quantiles/0"),
        ((0.2, 1.0), "/outputs/quantiles/1"),
        (tuple(i / 18 for i in range(1, 18)), "/outputs/quantiles"),
    ):
        refused_with(
            lambda q=given: rules.ensemble(sir, quantiles=q, **settings),
            ("CRN-R020", pointer),
            f"quantiles {given}",
        )
    for given in (0.5, "0.5", (True,), ("0.5",), {0.5: 1}):
        try:
            rules.ensemble(sir, quantiles=given, **settings)
        except TypeError:
            pass
        else:
            raise RuntimeError(f"quantiles={given!r} accepted")
    require(
        list(
            rules.ensemble(sir, quantiles=np.array([0.025, 0.975]), **settings)
            .simulate()
            .statistics["quantiles"]
        )
        == [0.025, 0.975],
        "a NumPy array of probabilities",
    )
    large = {
        "end": 3.0,
        "report_times": [i * 0.02 for i in range(150)],
        "replicates": 1_000_000,
        "max_output_bytes": 1 << 30,
    }
    hint = refused_with(
        lambda: rules.ensemble(model, quantiles=(0.5,), **large),
        ("CRN-R020", "/outputs/quantiles"),
        "the quantile buffer above the workspace limit",
    )["hint"]
    require(
        all(
            s in hint
            for s in ("fewer report times", "fewer replicates", "no quantiles")
        ),
        "the workspace refusal names its remedies",
    )
    rules.ensemble(model, **large)  # admitted without quantiles

    # Compiled loops serve exact bytes: one compiled for the model/1.0
    # document runs quantiles (here a stand-in that runs the interpreter);
    # one for the rules-model bytes is refused with a hint.
    def stand_in(document):
        loop = CompiledLoop(
            fingerprint=fingerprint_of(rules.loop_source(document)),
            model_sha256=hashlib.sha256(rules.model_bytes(document)).hexdigest(),
            version=crn_craft._native.__version__,
            directory=tempfile.gettempdir(),
            cached=True,
        )
        return loop

    old_loop = stand_in(sir.to_dict())
    hint = refused_with(
        lambda: rules.ensemble(sir, quantiles=QUANTILES, compiled=old_loop, **settings),
        ("CRN-R020", "/model"),
        "a loop for the rules-model bytes with quantiles",
    )["hint"]
    require("crn-craft.model/1.0" in hint, "the hint names the model/1.0 document")
    rules.ensemble(sir, compiled=old_loop, **settings)  # without quantiles: fine
    new_loop = stand_in({**sir.to_dict(), "schema": rules.MODEL_1_0})
    compiled_run = rules.ensemble(
        sir, quantiles=QUANTILES, compiled=new_loop, **settings
    )
    new_loop.execute = lambda request, cancel_file=None: compiled_run._native.simulate(
        None, None
    )
    verified = new_loop.verify(compiled_run)
    check_numpy_quantiles(verified, QUANTILES, "compiled stand-in")
    require(
        result_bytes(verified) == result_bytes(compiled_run.simulate()),
        "the compiled path reports the quantiles",
    )


def check_reports(model):
    """Result.save: markdown and JSON reports of builder, attributed and
    hybrid models, with and without paths, statistics and quantiles."""
    sir = rules.ensemble(
        sir_model(),
        end=250.0,
        report_times=[0, 50, 100, 150, 200, 250],
        replicates=10,
        seed=1,
        keep_paths="all",
        method="static_direct",
    ).simulate()
    before = (copy.deepcopy(sir.metadata), result_bytes(sir))
    with tempfile.TemporaryDirectory() as folder:
        path, markdown = saved(sir, folder, "sir.md")
        require(
            isinstance(path, Path) and path == Path(folder) / "sir.md",
            "save returns the path as a pathlib.Path",
        )
        check_report_markdown(sir, markdown, ["S", "I", "R"])
        require(
            markdown.startswith(
                "# Model\n\nName: sir\n\n## Reactions\n\n```text\n"
                "S + I -> I + I, beta * S * I\nI -> R, gamma * I\n```\n"
                "\n## Start values\n\n```text\nS = 999\nI = 1\nR = 0\n```\n"
                "\n## Parameter values\n\n```text\nbeta = 0.0002\ngamma = 0.01\n```\n"
                "\n## Stop condition\n\n```text\nI == 0\n```\n"
                "\n# Simulation\n\n```text\nseed = 1\nend = 250\n"
                "replicates = 10\nreport_times = 0, 50, 100, 150, 200, 250\n"
                "method = static_direct\n"
            )
            and "\n```\n\n# Results\n\n## Terminations\n\n```text\n" in markdown
            and markdown.endswith("\n```\n")
            and "## Species types" not in markdown,
            "markdown model section of the SIR example",
        )
        path, text = saved(sir, folder, "sir.json")
        check_report_json(sir, read_report(path))
        # Result numbers in the report form; the documents as they ran.
        require(
            '\n    "report_times": [0, 50, 100, 150, 200, 250],\n' in text
            and '\n      "report_times": [0.0, 50.0, 100.0, 150.0, 200.0, 250.0]\n'
            in text
            and '\n    "observables": ["I", "R", "S"],\n' in text,
            "JSON report number forms and observables in result order",
        )
        # The extension decides, case-insensitively; a file is replaced.
        (Path(folder) / "SIR.MARKDOWN").write_text("old")
        for name, same in (("SIR.MD", markdown), ("SIR.MARKDOWN", markdown)):
            require(saved(sir, folder, name)[1] == same, f"{name} is markdown")
        require(saved(sir, folder, "SIR.Json")[1] == text, "SIR.Json is JSON")
        # paths=False leaves the trajectories out.
        _, short = saved(sir, folder, "short.md", paths=False)
        check_report_markdown(sir, short, ["S", "I", "R"], paths=False)
        _, short = saved(sir, folder, "short.json", paths=False)
        check_report_json(sir, read_report(Path(folder) / "short.json"), paths=False)
        for bad in ("sir.txt", "sir", "sir.md.bak", "sir.jsonl"):
            try:
                sir.save(Path(folder) / bad)
            except ValueError as error:
                require(
                    all(s in str(error) for s in (".md", ".markdown", ".json")),
                    "the refusal names the supported extensions",
                )
            else:
                raise RuntimeError(f"unsupported report extension {bad} accepted")
            require(not (Path(folder) / bad).exists(), "a refused save writes nothing")
        try:
            sir.save(Path(folder) / "flag.md", paths=1)
        except TypeError:
            pass
        else:
            raise RuntimeError("paths=1 accepted")
        # A failed write raises OSError and leaves no file behind.
        names = sorted(p.name for p in Path(folder).iterdir())
        target = Path(folder) / "directory.md"
        target.mkdir()
        for path in (Path(folder) / "missing" / "sir.md", target):
            try:
                sir.save(path)
            except OSError:
                pass
            else:
                raise RuntimeError(f"saving to {path} did not fail")
        require(
            sorted(p.name for p in Path(folder).iterdir())
            == sorted([*names, "directory.md"])
            and not list(target.iterdir()),
            "a failed write leaves no temporary file",
        )
    require(
        (sir.metadata, result_bytes(sir)) == before,
        "saving changes nothing about the result",
    )
    # The engine's quantiles: after the sample variance, in ascending order,
    # 0.5 as the median.
    quantiled = rules.ensemble(
        sir_model(),
        end=250.0,
        report_times=[0, 50, 100, 150, 200, 250],
        replicates=10,
        seed=1,
        keep_paths="all",
        method="static_direct",
        quantiles=(0.025, 0.25, 0.5, 0.975),
    ).simulate()
    check_numpy_quantiles(quantiled, (0.025, 0.25, 0.5, 0.975), "report")
    with tempfile.TemporaryDirectory() as folder:
        _, markdown = saved(quantiled, folder, "q.md")
        check_report_markdown(quantiled, markdown, ["S", "I", "R"])
        require(
            [h for h in markdown.split("\n") if h.startswith("### ")]
            == [
                "### mean",
                "### mean_standard_error",
                "### sample_variance",
                "### 2.5% percentile",
                "### 25% percentile",
                "### median",
                "### 97.5% percentile",
            ],
            "quantile subsection titles",
        )
        path, _ = saved(quantiled, folder, "q.json")
        report = read_report(path)
        check_report_json(quantiled, report)
        require(
            [q["p"] for q in report["result"]["statistics"]["quantiles"]]
            == [0.025, 0.25, 0.5, 0.975],
            "JSON quantiles in ascending order",
        )
    # Withheld statistics and partial paths (NaN after a budget stop).
    partial = rules.ensemble(
        model,
        end=3.0,
        report_times=[0.0, 1.0, 2.0, 3.0],
        replicates=6,
        seed=5,
        max_events=20,
        keep_paths=True,
        method="static_direct",
    ).simulate()
    require(
        partial.statistics is None and bool(partial.paths.isnull().any()),
        "the partial run withholds statistics and has NaN paths",
    )
    order = partial.metadata["observables"]
    with tempfile.TemporaryDirectory() as folder:
        _, markdown = saved(partial, folder, "partial.md")
        check_report_markdown(partial, markdown, order)
        path, text = saved(partial, folder, "partial.json")
        check_report_json(partial, read_report(path))
        require(
            '"nan"' in text and "NaN" not in text and ", nan" in markdown,
            "non-finite values are the strings nan, inf and -inf",
        )
    # The three-class model: attributes, filtered observables, rules that
    # need JSON, attributed start values.
    attributed = rules.ensemble(
        model,
        end=3.0,
        report_times=[0.0, 1.0, 3.0],
        replicates=5,
        seed=2,
        keep_paths=2,
        method="static_direct",
    ).simulate()
    with tempfile.TemporaryDirectory() as folder:
        _, markdown = saved(attributed, folder, "mtdna.md")
        check_report_markdown(attributed, markdown, attributed.metadata["observables"])
        found = dict(markdown_sections(markdown))
        rules_json = [
            json.dumps(r, separators=(",", ":")) for r in attributed.model["rules"]
        ]
        require(
            found["## Species types"] == ["mtDNA: clone, s, wt"]
            and found["## Reactions"]
            == [
                f"{r['name']}: {j}"
                for r, j in zip(attributed.model["rules"], rules_json, strict=True)
            ]
            and found["## Start values"] == ["mtDNA {clone=0, s=sw, wt=1} = w0"]
            and found["## Parameter values"][0] == "K = 20"
            and "## Single trajectories for "
            + ", ".join(attributed.metadata["observables"])
            in found
            and len(
                found[
                    "## Single trajectories for "
                    + ", ".join(attributed.metadata["observables"])
                ]
            )
            == 2 * len(attributed.metadata["observables"]),
            "markdown sections of an attributed model",
        )
        path, _ = saved(attributed, folder, "mtdna.json")
        check_report_json(attributed, read_report(path))
    # A continuous variable: its section, and its column after the species.
    ramp = rules.Model("ramp")
    ramp.add_species("A")
    ramp.add_ode("x", 0.0, "2", lower_bound=0, on_bound="clamp")
    ramp.add_reaction("x", [], ["A"], name="immigration")
    hybrid = rules.ensemble(
        ramp,
        end=3.0,
        report_times=[0.0, 1.5, 3.0],
        replicates=8,
        seed=13,
        keep_paths=3,
        method="hybrid",
    ).simulate()
    with tempfile.TemporaryDirectory() as folder:
        _, markdown = saved(hybrid, folder, "ramp.md")
        check_report_markdown(hybrid, markdown, ["A", "x"])
        found = dict(markdown_sections(markdown))
        require(
            found["## Continuous variables"]
            == ["x: initial = 0, derivative = 2, lower_bound = 0, on_bound = clamp"]
            and found["## Reactions"] == ["0 -> A, x"]
            and found["## Start values"] == ["A = 0"]
            and "## Parameter values" not in found
            and "method = hybrid, rtol = 1e-06, atol = 1e-09" in found["# Simulation"],
            "markdown sections of a model with a continuous variable",
        )
        path, _ = saved(hybrid, folder, "ramp.json")
        check_report_json(hybrid, read_report(path))
    # load_report checks the schema string.
    with tempfile.TemporaryDirectory() as folder:
        other = Path(folder) / "other.json"
        for document in ({"schema": "crn-craft.report/0.9"}, [1], {"result": {}}):
            other.write_text(json.dumps(document))
            try:
                rules.load_report(other)
            except ValueError:
                pass
            else:
                raise RuntimeError(f"load_report accepted {document}")
    require(crn_craft.load_report is rules.load_report, "crn_craft.load_report")


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
        framed.paths.sizes["replicate"] == 20
        and framed.paths.sizes["time"] == 3
        and len(framed.paths_frame()) == 20 * 3
        and framed.paths_frame().index.get_level_values(0).nunique() == 20,
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
    saved = partial.paths_frame().groupby(level="replicate").size()
    filled = partial.paths.notnull().all("variable").sum("time")
    require(
        partial.statistics is None
        and set(partial.records["termination"]) == {"event_limit"}
        and saved.max() < 4
        and np.array_equal(saved.to_numpy(), partial.records["saved_points"].to_numpy())
        and np.array_equal(filled.values, saved.to_numpy())
        and bool(partial.paths.isnull().any()),
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
    for result in (framed, partial):
        data = result.to_xarray()
        names = result.metadata["observables"]
        cube = np.stack([data[name].values for name in names], axis=-1)
        require(
            list(data.data_vars) == names
            and data[names[0]].dims == ("replicate", "time")
            and list(data["time"].values) == result.metadata["report_times"]
            and list(data["replicate"].values)
            == [p["id"] for p in result.metadata["paths"]]
            and data.attrs
            == {
                key: result.metadata[key]
                for key in ("sampler", "model_sha256", "run_sha256", "root_seed")
            },
            "to_xarray coordinates",
        )
        for kept in result.metadata["paths"]:
            rows = np.array(kept["values"])
            require(
                np.array_equal(cube[kept["id"], : len(rows)], rows)
                and np.isnan(cube[kept["id"], len(rows) :]).all(),
                "to_xarray values",
            )
    # Event arrays: replicates without a logged event among logged ones, and
    # a run of one replicate.
    for replicates, seed in ((12, 3), (1, 5)):
        sparse = rules.ensemble(
            model,
            end=3.0,
            report_times=[0.0, 3.0],
            replicates=replicates,
            seed=seed,
            method="static_direct",
            event_log={"max_entries": 3, "rules": ["mutation_mtAdv"]},
        )
        metadata, columns = sparse._native.simulate(None, None)
        metadata = json.loads(metadata)
        logged = rules.Result(metadata, columns)
        check_frames(logged, metadata, columns)
        lengths = logged.event_columns()["length"].values
        require(
            logged.event_columns().sizes["replicate"] == replicates
            and lengths.max() > 0
            and (replicates == 1 or lengths.min() == 0),
            "event array fixtures: empty and nonempty logs",
        )
    # rules-run/0.6: "static_direct" defaults to Ziggurat waits and sparse
    # trackers; explicit options give the numbers of the 0.3.0 request.
    counts = [{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]
    base = {
        "end": 3.0,
        "report_times": [0.0, 1.0, 3.0],
        "replicates": 12,
        "seed": 9,
        "keep_paths": "all",
        "trackers": counts,
    }
    new = rules.ensemble(model, method="static_direct", **base).simulate()
    require(
        new.metadata["schema"] == "crn-craft.rules-result/0.3"
        and new.metadata["method"] == {"kind": "static_direct"}
        and new.metadata["sampler"]
        == "chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1"
        and not new.metadata["diagnostics"],
        "static_direct selects the rules-run/0.6 defaults",
    )
    inverse = {"kind": "static_direct", "waiting_time": "inverse"}
    restored = rules.ensemble(
        model, method=inverse | {"sparse_trackers": False}, **base
    ).simulate()
    released = rules.prepare(
        model,
        {
            "schema": "crn-craft.rules-run/0.2",
            "model": {
                "path": "model.rules.json",
                "sha256": hashlib.sha256(model).hexdigest(),
            },
            "time": {"end": 3.0, "report_times": [0.0, 1.0, 3.0]},
            "ensemble": {"replicates": 12, "root_seed": 9},
            "outputs": {
                "statistics": ["mean", "variance"],
                "keep_paths": 12,
                "keep_instances": 0,
            },
            "budget": {
                "threads": 1,
                "max_events_per_replicate": 1_000_000,
                "max_sampling_steps_per_replicate": 8_000_000,
                "max_instances": 100_000,
                "max_output_bytes": 64 << 20,
            },
            "trackers": counts,
            "method": {"kind": "static_direct"},
        },
    ).simulate()
    require(
        restored.metadata["sampler"]
        == released.metadata["sampler"]
        == "chacha8-u53/rules-static-direct/1"
        and result_bytes(restored) == result_bytes(released)
        and restored.paths.equals(released.paths)
        and restored.observations == released.observations,
        "explicit options reproduce the 0.3.0 static_direct numbers",
    )
    # Without an eligible tracker the default is silent.
    quiet = rules.ensemble(
        model,
        method="static_direct",
        **(
            base
            | {"trackers": [{"observable": "mutants", "first_passage": "mutants >= 2"}]}
        ),
    ).simulate()
    require(
        quiet.metadata["sampler"] == "chacha8-u53-ziggurat/rules-static-direct/1"
        and not quiet.metadata["diagnostics"],
        "static_direct without an eligible tracker",
    )

    # The documented mtdna_direct fallback: a stop outside the certified form
    # runs generic Static Direct as in 0.3.0 (inverse waits, per-event
    # trackers), with its numbers and a CRN-R051 warning.
    odd = json.loads(model)
    odd["stop"] = "mutants > 2"
    settings = {
        "end": 3.0,
        "report_times": [0.0, 1.0, 3.0],
        "replicates": 12,
        "seed": 9,
        "keep_paths": "all",
        "trackers": [{"observable": "mutants"}],
    }
    fell = rules.ensemble(odd, method="mtdna_direct", **settings).simulate()
    plain = rules.ensemble(
        odd, method=inverse | {"sparse_trackers": False}, **settings
    ).simulate()
    require(
        fell.metadata["sampler"]
        == plain.metadata["sampler"]
        == "chacha8-u53/rules-static-direct/1"
        and result_bytes(fell) == result_bytes(plain)
        and fell.paths.equals(plain.paths)
        and fell.observations == plain.observations
        and [d["code"] for d in fell.metadata["diagnostics"]] == ["CRN-R051"]
        and not plain.metadata["diagnostics"],
        "mtdna_direct fallback equals static_direct",
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
    check_result_documents(model)
    check_fsp()
    check_archives()
    check_quantiles(model)
    check_reports(model)
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
                    "path and event arrays match the raw metadata",
                    "event arrays with empty logs, leaps and one replicate",
                    "static_direct selects the rules-run/0.6 defaults",
                    "explicit options reproduce the 0.3.0 static_direct numbers",
                    "static_direct without an eligible tracker",
                    "mtdna_direct fallback equals static_direct",
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
                    "result model and request of every path",
                    "rules.fsp law and refusals",
                    "rules archives and their verification",
                    "quantiles equal numpy.quantile of the paths",
                    "markdown and JSON reports",
                ],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
