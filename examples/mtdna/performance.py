"""Timings behind docs/mtdna-performance-analysis.md.

Usage: python examples/mtdna/performance.py SECTION [CRN_CRAFT_BINARY]

SECTION is methods, threads, copy, observation, interfaces (needs the native
CLI), models, results or all. Writes performance-SECTION.json in the working
directory. Timing is host-specific.
"""

import hashlib
import json
import platform
import subprocess
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas
import xarray
from crn_craft import rules

OUT = {
    "platform": platform.platform(),
    "python": sys.version.split()[0],
    "numpy": np.__version__,
    "pandas": pandas.__version__,
    "xarray": xarray.__version__,
    "rows": [],
}
BASE = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
STOP = dict(BASE, stop="N > 0 && mutants >= 0.9 * N")
GRID = 24.0 * np.arange(0, 1081, 30)
TRACK = [{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]


def timed(label, ensemble, repeats=3, **kwargs):
    seconds = []
    for _ in range(repeats):
        start = time.perf_counter()
        result = ensemble.simulate(**kwargs)
        seconds.append(time.perf_counter() - start)
    meta = result.metadata
    stopped = int((result.records["termination"] == "stopped").sum())
    row = {
        "case": label,
        "seconds": sorted(seconds)[len(seconds) // 2],
        "all": seconds,
        "events": meta["events"],
        "cells": len(result.records),
        "stopped": stopped,
        "sampler": meta["sampler"],
    }
    row["ns_per_event"] = 1e9 * row["seconds"] / row["events"]
    OUT["rows"].append(row)
    print(
        f"{label:48s} {row['seconds']:7.2f} s {row['ns_per_event']:6.1f} ns/event "
        f"{row['events']:>12,} events stopped {stopped}",
        flush=True,
    )
    return result


def ensemble(model=STOP, cells=1000, threads=2, method="mtdna_direct", **kwargs):
    kwargs.setdefault("report_times", GRID)
    return rules.ensemble(
        model,
        end=25_920.0,
        replicates=cells,
        seed=17,
        threads=threads,
        method=method,
        **kwargs,
    )


section = sys.argv[1] if len(sys.argv) > 1 else "all"
if section in ("methods", "all"):
    for label, method, cells, repeats in [
        ("tree", None, 200, 1),
        # The rules-run/0.6 default: Ziggurat waits and sparse trackers.
        ("static_direct", "static_direct", 1000, 3),
        (
            "static_direct inverse, per-event trackers (0.3.0 default)",
            {
                "kind": "static_direct",
                "waiting_time": "inverse",
                "sparse_trackers": False,
            },
            1000,
            3,
        ),
        ("mtdna_direct", "mtdna_direct", 1000, 3),
        ("next_reaction", "next_reaction", 1000, 1),
        ("tau epsilon 0.03", {"kind": "tau", "epsilon": 0.03}, 200, 1),
        ("binomial 1 h", {"kind": "binomial", "step": 1.0}, 1000, 3),
    ]:
        timed(
            f"method {label} ({cells} cells, trackers)",
            ensemble(cells=cells, method=method, trackers=TRACK),
            repeats,
        )
if section in ("threads", "all"):
    for threads in (1, 2):
        timed(
            f"threads {threads} (2000 cells, mtdna_direct)",
            ensemble(cells=2000, threads=threads),
        )
if section in ("copy", "all"):
    for k in (250, 500, 1000, 2000):
        model = json.loads(json.dumps(STOP))
        model["parameters"].update({"K": k, "w0": k})
        timed(f"K {k} (500 cells, mtdna_direct)", ensemble(model, cells=500))
if section in ("observation", "all"):
    timed("observation none (1000 cells)", ensemble())
    timed("observation 2 count trackers", ensemble(trackers=TRACK))
    timed(
        "observation first passage tracker",
        ensemble(
            trackers=[
                {
                    "observable": "mutants",
                    "first_passage": "N > 0 && mutants >= 0.6 * N",
                }
            ]
        ),
    )
    timed(
        "observation mutation log (50 per cell)",
        ensemble(
            event_log={"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]},
            max_output_bytes=1 << 30,
        ),
    )
    small = ensemble(
        cells=20,
        threads=1,
        method={
            "kind": "static_direct",
            "waiting_time": "inverse",
            "sparse_trackers": False,
        },
    )
    timed("callback none (20 cells, static_direct, 1 thread)", small)
    timed(
        "callback Python (20 cells, static_direct, 1 thread)",
        small,
        callback=lambda replicate, event: True,
    )
if section in ("interfaces", "all"):
    work = Path(tempfile.mkdtemp())
    data = (json.dumps(STOP) + "\n").encode()
    (work / "p.rules.json").write_bytes(data)
    request = {
        "schema": "crn-craft.rules-run/0.3",
        "model": {"path": "p.rules.json", "sha256": hashlib.sha256(data).hexdigest()},
        "time": {"end": 25_920.0, "report_times": list(GRID)},
        "ensemble": {"replicates": 1000, "root_seed": 17},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": 0,
            "keep_instances": 0,
        },
        "budget": {
            "threads": 2,
            "max_events_per_replicate": 1_000_000,
            "max_sampling_steps_per_replicate": 8_000_000,
            "max_instances": 100_000,
            "max_output_bytes": 64 << 20,
        },
        "method": {"kind": "mtdna_direct"},
    }
    (work / "p.run.json").write_text(json.dumps(request))
    cli = sys.argv[2]
    seconds = []
    for _ in range(3):
        start = time.perf_counter()
        subprocess.run(
            [cli, "rules", "p.rules.json", "p.run.json"],
            cwd=work,
            check=True,
            capture_output=True,
        )
        seconds.append(time.perf_counter() - start)
    OUT["rows"].append(
        {
            "case": "interface CLI (1000 cells, mtdna_direct)",
            "seconds": sorted(seconds)[1],
            "all": seconds,
        }
    )
    print("interface CLI", sorted(seconds)[1], flush=True)
    timed("interface Python (1000 cells, mtdna_direct)", rules.prepare(data, request))
    template = ensemble(cells=1000)
    start = time.perf_counter()
    for i in range(20):
        template.bind_parameters({"mu": 3e-5 * (1 + i / 100)})
    bind = (time.perf_counter() - start) / 20
    OUT["rows"].append({"case": "bind_parameters per candidate", "seconds": bind})
    print(f"bind_parameters per candidate {bind * 1e3:.2f} ms", flush=True)
    start = time.perf_counter()
    for _ in range(20):
        ensemble(cells=1000)
    fresh = (time.perf_counter() - start) / 20
    OUT["rows"].append({"case": "fresh rules.ensemble per candidate", "seconds": fresh})
    print(f"fresh ensemble per candidate {fresh * 1e3:.2f} ms", flush=True)
if section in ("models", "all"):
    clones = json.loads(rules.load("models/rules/mtdna-fixed.rules.json"))
    clones["parameters"].update(
        {"K": 1000, "w0": 1000, "mu": 3e-5, "delta": BASE["parameters"]["delta"]}
    )
    timed(
        "clone model, tree (100 cells)",
        ensemble(clones, cells=100, method=None, report_times=[0.0, 25_920.0]),
        1,
    )
    timed(
        "three-class, tree (100 cells)",
        ensemble(BASE, cells=100, method=None, report_times=[0.0, 25_920.0]),
        1,
    )
    atp = json.loads(rules.load("models/rules/mtdna-three-class-atp.rules.json"))
    timed(
        "ATP model, thinning (512 cells)",
        ensemble(atp, cells=512, method="thinning"),
        1,
    )
if section in ("results", "all"):
    # What the Python result costs next to the simulation: 4000 cells with
    # every path kept. The native call returns the result document as JSON
    # and the per-cell records as buffers; each later step is timed alone.
    # paths (a DataArray) builds the path arrays that paths_frame() and
    # to_xarray() then reuse; "paths_frame() first" times the 0.4.0 frame
    # from a fresh Result, including those arrays.
    log = {"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]}
    full = {"max_entries": 50, "rules": []}  # every cell logs 50 events

    def array_bytes(array):
        """Bytes of a DataArray's values and coordinates (pointers for objects)."""
        return array.nbytes + sum(c.nbytes for c in array.coords.values())

    def frame_bytes(frame):
        return int(frame.memory_usage(deep=True, index=True).sum())

    for label, grid, extra in [
        ("monthly reports", GRID, {}),
        ("weekly reports", 24.0 * np.arange(0, 1081, 7), {}),
        ("monthly reports, mutation log", GRID, {"event_log": log}),
        ("monthly reports, full 50-entry log", GRID, {"event_log": full}),
    ]:
        prepared = ensemble(
            cells=4000,
            keep_paths="all",
            report_times=grid,
            max_output_bytes=1 << 30,
            **extra,
        )
        steps = {}
        for _ in range(3):
            # Free the previous repeat first, so that no step times its release.
            text = buffers = document = result = fresh = paths = None
            clock = time.perf_counter()
            text, buffers = prepared._native.simulate(None, None)
            lap = {"simulate (native)": time.perf_counter() - clock}
            clock = time.perf_counter()
            document = json.loads(text)
            lap["parse JSON"] = time.perf_counter() - clock
            clock = time.perf_counter()
            result = rules.Result(document, buffers)
            lap["records, rule_events, statistics"] = time.perf_counter() - clock
            clock = time.perf_counter()
            paths = result.paths
            lap["paths"] = time.perf_counter() - clock
            clock = time.perf_counter()
            result.paths_frame()
            lap["paths_frame()"] = time.perf_counter() - clock
            clock = time.perf_counter()
            result.to_xarray()
            lap["to_xarray"] = time.perf_counter() - clock
            fresh = rules.Result(document, buffers)
            clock = time.perf_counter()
            fresh.paths_frame()
            lap["paths_frame() first"] = time.perf_counter() - clock
            if result.observations:
                clock = time.perf_counter()
                result.event_columns()
                lap["event_columns()"] = time.perf_counter() - clock
                clock = time.perf_counter()
                result.event_frame()
                lap["event_frame()"] = time.perf_counter() - clock
            for key, value in lap.items():
                steps.setdefault(key, []).append(value)
        # Memory, untimed: the peak of traced allocations during each step
        # (NumPy reports its buffers to tracemalloc) and the size of what it
        # returns.
        memory = {}
        conversions = [
            ("paths", lambda r: r.paths, array_bytes),
            ("paths_frame() first", lambda r: r.paths_frame(), frame_bytes),
        ]
        if result.observations:
            conversions += [
                ("event_columns()", lambda r: r.event_columns(), array_bytes),
                ("event_frame()", lambda r: r.event_frame(), frame_bytes),
            ]
        for key, convert, size in conversions:
            fresh = rules.Result(document, buffers)
            tracemalloc.start()
            made = convert(fresh)
            peak = tracemalloc.get_traced_memory()[1]
            tracemalloc.stop()
            memory[key] = {"peak_traced_bytes": peak, "result_bytes": size(made)}
            made = None
        row = {
            "case": f"results 4000 cells, keep_paths all, {label}",
            "report_times": len(grid),
            "path_rows": int(paths.notnull().all("variable").sum()),
            "path_array_shape": dict(paths.sizes),
            "logged_events": sum(len(o["log"]) for o in result.observations),
            "events": result.metadata["events"],
            "seconds": {key: sorted(v)[1] for key, v in steps.items()},
            "memory": memory,
            "all": steps,
        }
        OUT["rows"].append(row)
        print(row["case"], json.dumps(row["seconds"]), flush=True)
        print("  memory", json.dumps(memory), flush=True)
Path(f"performance-{section}.json").write_text(json.dumps(OUT, indent=2) + "\n")
