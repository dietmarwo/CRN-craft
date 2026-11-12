"""Same-session speed comparison: the 0.6.0 network engine and the unified engine.

For every model of `validation/merge/baseline-0.6.0.json` (or `--models`), the
0.6.0 best method of the baseline's CLI summary is rerun now with the kept
0.6.0 CLI (`--old`), and the model converted by the engine under test
(`--new`, `crn-craft convert`) runs as a `crn-craft.run/1.0` Static Direct
request (its defaults: Ziggurat waits, sparse trackers where eligible) with
each selector of `--selectors`, at the sizes of the baseline row: same
horizon, report times, replicates, root seed, one thread, statistics mean and
variance, no kept paths. The time of a run is the wall time of one CLI process
(`crn-craft ssa` for 0.6.0, `crn-craft rules` for the unified engine), stdout
captured; one untimed warm-up, then `--repeats` timed runs, and the median.
All runs of a row must print the same bytes.

The output holds per model the 0.6.0 time and, per selector, the time, the
ratio to the 0.6.0 time, the sampler identity and the event count. Run it
under the CPU lock:

    flock /home/claude/craft5/cpu.lock python3 validation/merge/compare_speed.py \\
        --old /home/claude/m07/bin060/crn-craft --new target/release/crn-craft \\
        --models W1,W2,W3,W4,W5,W6,W7 --output speed.json

Only the standard library is needed. M4c records these timings as
information; the speed gate of M4d (declared in `validation/merge/RESULTS.md`)
compares the best selector against a same-session rerun of `baseline.py`.
"""

import argparse
import hashlib
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "validation" / "merge"
SELECTORS = ["linear", "composition_rejection", "sorting", "rssa_cr", "auto"]


def dump(value):
    return (json.dumps(value, separators=(",", ":")) + "\n").encode()


def timed(binary, args, cwd, repeats):
    """Median wall time of `repeats` runs after a warm-up, and the output."""
    runs = []
    output = None
    for i in range(repeats + 1):
        start = time.perf_counter()
        done = subprocess.run(
            [str(binary), *args],
            cwd=cwd,
            capture_output=True,
            check=False,
        )
        seconds = time.perf_counter() - start
        if done.returncode != 0:
            raise RuntimeError(f"{args}: {done.stderr.decode()[:400]}")
        if output is not None and done.stdout != output:
            raise RuntimeError(f"{args}: runs printed different bytes")
        output = done.stdout
        if i:
            runs.append(seconds)
    return statistics.median(runs), runs, json.loads(output)


def unified_request(row, model_bytes, converted, selector):
    request = row["request"]
    method = {"kind": "static_direct"}
    if selector != "default":
        method["selector"] = selector
    return {
        "schema": "crn-craft.run/1.0",
        "model": {
            "path": "conv.json",
            "sha256": hashlib.sha256(model_bytes).hexdigest(),
        },
        "time": {
            "end": request["time"]["end"],
            "report_times": request["time"]["save_times"],
        },
        "ensemble": {
            "replicates": request["ensemble"]["replicates"],
            "root_seed": request["ensemble"]["root_seed"],
        },
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": 0,
            "keep_instances": 0,
        },
        "budget": {
            "max_events_per_replicate": 10**9,
            "max_sampling_steps_per_replicate": 10**12,
            "max_instances": max(1, len(converted["types"])),
            "threads": 1,
            "max_output_bytes": 1 << 30,
        },
        "method": method,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--old", type=Path, help="the kept 0.6.0 CLI (omit to skip)")
    parser.add_argument("--new", type=Path, required=True, help="CLI under test")
    parser.add_argument("--models", help="comma-separated model names (default: all)")
    parser.add_argument(
        "--selectors",
        default=",".join(SELECTORS),
        help="comma-separated selectors; 'default' omits the field",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.new = args.new.resolve()
    if args.old is not None:
        args.old = args.old.resolve()
    baseline = json.loads((HERE / "baseline-0.6.0.json").read_text())
    summary = baseline["timing_summary"]["cli"]
    manifest = json.loads((HERE / "models.json").read_text())
    wanted = args.models.split(",") if args.models else None
    selectors = args.selectors.split(",")
    rows_by_model = {}
    for row in baseline["rows"]:
        rows_by_model.setdefault(row["model"], []).append(row)
    report = {
        "schema": "crn-craft.merge-speed-comparison/0.1",
        "host": {"machine": platform.machine(), "python": platform.python_version()},
        "old": None if args.old is None else str(args.old),
        "new": str(args.new),
        "repeats": args.repeats,
        "boundary": "wall time of one CLI process, median after one warm-up",
        "models": [],
    }
    for name, rows in rows_by_model.items():
        if wanted is not None and name not in wanted:
            continue
        best = summary[name]["best"]
        aggregator, waiting, fast = best.split("/")
        row = next(
            r
            for r in rows
            if r["aggregator"] == aggregator
            and r["waiting_time"] == waiting
            and r["fast"] == (fast == "fast")
        )
        entry = next(m for m in manifest["rows"] if m["name"] == name)
        network = (ROOT / entry["path"]).read_bytes()
        result = {
            "model": name,
            "old_method": best,
            "old_events": row["events"],
            "old_baseline_seconds": summary[name]["best_seconds"],
            "selectors": {},
        }
        with tempfile.TemporaryDirectory() as scratch:
            scratch = Path(scratch)
            (scratch / "model.json").write_bytes(network)
            if args.old is not None:
                (scratch / "request.json").write_bytes(dump(row["request"]))
                seconds, runs, _ = timed(
                    args.old,
                    ["ssa", "model.json", "request.json"],
                    scratch,
                    args.repeats,
                )
                result.update(old_seconds=seconds, old_runs=runs)
            done = subprocess.run(
                [str(args.new), "convert", "model.json"],
                cwd=scratch,
                capture_output=True,
                check=False,
            )
            if done.returncode != 0:
                raise RuntimeError(f"{name}: {done.stderr.decode()[:400]}")
            converted_bytes = done.stdout
            (scratch / "conv.json").write_bytes(converted_bytes)
            converted = json.loads(converted_bytes)
            for selector in selectors:
                request = unified_request(row, converted_bytes, converted, selector)
                (scratch / "run.json").write_bytes(dump(request))
                seconds, runs, out = timed(
                    args.new, ["rules", "conv.json", "run.json"], scratch, args.repeats
                )
                item = {
                    "seconds": seconds,
                    "runs": runs,
                    "sampler": out["sampler"],
                    "events": out["events"],
                    "termination": out["termination"],
                }
                if "old_seconds" in result:
                    item["ratio"] = seconds / result["old_seconds"]
                result["selectors"][selector] = item
                ratio = f"  ratio {item['ratio']:5.2f}" if "ratio" in item else ""
                print(
                    f"{name:18s} {selector:22s} {seconds:8.4f} s{ratio}  {out['sampler']}",
                    file=sys.stderr,
                    flush=True,
                )
        if "old_seconds" in result:
            print(
                f"{name:18s} {'0.6.0 ' + best:22s} {result['old_seconds']:8.4f} s",
                file=sys.stderr,
                flush=True,
            )
        report["models"].append(result)
        args.output.write_text(json.dumps(report, indent=1) + "\n")
    args.output.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
