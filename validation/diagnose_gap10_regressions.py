"""Independent follow-up for exposed regressions; never overwrites primary timings."""

import json
from pathlib import Path
import statistics
import subprocess
import time
from compare_julia_accelerated import rust_job
from compare_julia_gap10 import Baseline, ROOT

FOLDER = ROOT / "validation/julia-gap10"
CASES = [
    ("W3", "rssa_cr/inverse"),
    ("W4", "cached_direct/ziggurat"),
    ("W4", "direct_cr/inverse"),
    ("W4", "direct_cr/ziggurat"),
    ("W5", "cached_direct/inverse"),
    ("W7", "legacy_cached"),
]


def main():
    original = json.loads((FOLDER / "comparison.json").read_text())
    report = dict(
        note="Fresh diagnostic pairs; primary observations and acceptance decisions are unchanged.",
        python={},
        native={},
    )
    baseline = Baseline(Path("/tmp/crn-gap10-base-package"))
    try:
        for name, engine in CASES:
            entry = next(v for v in original["workloads"]["rows"] if v["name"] == name)
            end = entry["end"]
            paths = entry["timing_paths"]
            model = (
                ROOT / f"validation/ssa-acceleration/models/{name}.model.json"
            ).read_bytes()
            baseline.call(
                command="prepare",
                model_path=str(
                    ROOT / f"validation/ssa-acceleration/models/{name}.model.json"
                ),
                end=end,
                paths=paths,
            )
            job = rust_job(engine, model, end, paths, 4242, False, "prepared")
            job()
            baseline.call(command="run", engine=engine)
            rows = []
            for i in range(5):
                pair = {}
                for side in (
                    ["baseline", "final"] if i % 2 == 0 else ["final", "baseline"]
                ):
                    obs = []
                    for _ in range(2):
                        if side == "baseline":
                            obs.append(
                                baseline.call(command="run", engine=engine)["seconds"]
                            )
                        else:
                            start = time.perf_counter()
                            job()
                            obs.append(time.perf_counter() - start)
                    pair[side] = obs
                pair["speedup"] = statistics.median(
                    pair["baseline"]
                ) / statistics.median(pair["final"])
                rows.append(pair)
            key = name + "/" + engine
            before = (
                original["timings"][name]["Baseline/" + engine + "/prepared"][
                    "median_seconds"
                ]
                / original["timings"][name][engine + "/prepared"]["median_seconds"]
            )
            report["python"][key] = dict(
                original_speedup=before,
                diagnostic_speedup=statistics.median(v["speedup"] for v in rows),
                rows=rows,
            )
            (FOLDER / "regression-diagnostics.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
            print(key, report["python"][key]["diagnostic_speedup"], flush=True)
    finally:
        baseline.close()
    # Localize the largest apparent regression across the actual kept commits,
    # using the prepared native boundary and the same unchanged numerical paths.
    method = "direct_cr"
    for wait, section in [("inverse", "native"), ("ziggurat", "native_ziggurat")]:
        report[section] = {}
        for label in ["base", "C05", "C06", "C09"]:
            obs = []
            for _ in range(3):
                reply = json.loads(
                    subprocess.check_output(
                        [
                            f"/tmp/crn-gap10-binaries/{label}",
                            str(
                                ROOT
                                / "validation/ssa-acceleration/models/W4.model.json"
                            ),
                            "100",
                            "64",
                            method,
                            wait,
                            "1",
                        ]
                    )
                )
                obs.append(reply)
            report[section][label] = dict(
                observations=obs,
                median_seconds=statistics.median(v["seconds"][0] for v in obs),
            )
            print(section, label, report[section][label]["median_seconds"], flush=True)
            (FOLDER / "regression-diagnostics.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )


if __name__ == "__main__":
    main()
