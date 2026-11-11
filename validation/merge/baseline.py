"""Time every exact SSA method of the 0.6.0 engine on the merge models (task M0c).

For each model of `validation/merge/models.json` (its `timing` sizes) and each
of the 21 methods (`direct`, `cached_direct`, `direct_cr`, `sorting_direct`,
`rssa_cr` with `fast` false/true and inverse/Ziggurat waits, and
`next_reaction`), this runs one untimed warm-up and three timed runs at two
boundaries and writes `validation/merge/baseline-0.6.0.json`:

- `cli` (gated): wall time of one `crn-craft ssa model.json request.json`
  process, stdout written to a file;
- `python` (informative): the native `simulate_json` call of a prepared
  `crn_craft.prepare_ssa` handle, without preparation and JSON decoding.

Each row keeps the exact request, the event count, the work counters, the
sampler string and result hashes; timings sit under `timing` and
`timing_summary`, which also evaluates the selector rule of `RESULTS.md`.

This process needs only the standard library; the Python boundary runs in a
worker started with `--python`. Run it under the CPU lock:

    flock /home/claude/craft5/cpu.lock python validation/merge/baseline.py \
        --python /home/claude/m07/venv060/bin/python \
        --binary /home/claude/m07/bin060/crn-craft

`--check-against FILE` compares every non-timing field of the rows written
with those of FILE and exits 1 on a difference.
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_common import (  # noqa: E402
    MANIFEST,
    ROOT,
    binary_identity,
    encode,
    host_info,
    load_manifest,
    run_worker,
    sha256,
    ssa_request,
)

DEFAULT_OUTPUT = ROOT / "validation/merge/baseline-0.6.0.json"
REPETITIONS = 3
SELECTOR_FACTOR = 1.3
WORKLOADS = [f"W{i}" for i in range(1, 8)]
CONFIGS = [
    (aggregator, waiting, fast)
    for aggregator in [
        "direct",
        "cached_direct",
        "direct_cr",
        "sorting_direct",
        "rssa_cr",
    ]
    for fast in [False, True]
    for waiting in ["inverse", "ziggurat"]
] + [("next_reaction", "inverse", False)]
COMMAND = ["crn-craft", "ssa", "model.json", "request.json"]


def label(aggregator, waiting, fast):
    return f"{aggregator}/{waiting}/{'fast' if fast else 'checked'}"


def canonical(result):
    """Hash of a parsed result, independent of the JSON layout."""
    return sha256(json.dumps(result, sort_keys=True).encode())


# ---------------------------------------------------------------- worker side


def worker_timing(job):
    import crn_craft as crn

    model = (ROOT / job["path"]).read_bytes()
    answers = []
    for request in job["requests"]:
        start = time.perf_counter()
        prepared = crn.prepare_ssa(model, encode(request))
        preparation = time.perf_counter() - start
        outputs, seconds = [], []
        for _ in range(job["repetitions"] + 1):
            start = time.perf_counter()
            raw = prepared._prepared.simulate_json(None)
            seconds.append(time.perf_counter() - start)
            outputs.append(raw)
        if any(raw != outputs[0] for raw in outputs):
            raise RuntimeError(f"{job['path']}: prepared replays differ")
        answers.append(
            {
                "output_sha256": sha256(outputs[0].encode()),
                "result_sha256": canonical(json.loads(outputs[0])),
                "preparation_seconds": preparation,
                "warmup_seconds": seconds[0],
                "seconds": seconds[1:],
            }
        )
    return answers


def worker_identity():
    import crn_craft as crn

    return {"python": sys.version, "crn_craft": crn.__version__}


def worker(job_path, answer_path):
    job = json.loads(Path(job_path).read_text())
    answer = worker_identity() if job["kind"] == "identity" else worker_timing(job)
    Path(answer_path).write_text(json.dumps(answer))


# ----------------------------------------------------------------- outer side


def time_cli(binary, folder, repetitions):
    outputs, seconds = [], []
    for _ in range(repetitions + 1):
        with open(folder / "result.json", "wb") as sink:
            start = time.perf_counter()
            done = subprocess.run(
                [str(binary), *COMMAND[1:]],
                cwd=folder,
                stdout=sink,
                stderr=subprocess.PIPE,
            )
            seconds.append(time.perf_counter() - start)
        if done.returncode != 0:
            raise RuntimeError(f"{folder.name}: exit {done.returncode} {done.stderr!r}")
        outputs.append((folder / "result.json").read_bytes())
    if any(data != outputs[0] for data in outputs):
        raise RuntimeError(f"{folder.name}: CLI replays differ")
    return outputs[0], seconds[0], seconds[1:]


def timing_entry(warmup, seconds, **extra):
    return {
        **extra,
        "warmup_seconds": warmup,
        "seconds": seconds,
        "median_seconds": statistics.median(seconds),
    }


def summarize(rows, boundary):
    """Best method per model and the selector-rule ratios, for one boundary."""
    out = {}
    for model in dict.fromkeys(row["model"] for row in rows):
        mine = [row for row in rows if row["model"] == model]
        best = min(mine, key=lambda row: row["timing"][boundary]["median_seconds"])
        per_method = {}
        for row in mine:
            t = row["timing"][boundary]["median_seconds"]
            key = row["aggregator"]
            per_method[key] = min(per_method.get(key, t), t)
        out[model] = {
            "best": label(best["aggregator"], best["waiting_time"], best["fast"]),
            "best_seconds": best["timing"][boundary]["median_seconds"],
            "best_events": best["events"],
            "best_per_method_seconds": per_method,
            "direct_cr_over_sorting_direct": per_method["direct_cr"]
            / per_method["sorting_direct"],
            "direct_cr_over_rssa_cr": per_method["direct_cr"] / per_method["rssa_cr"],
        }
    return out


def selector_rule(summary):
    if not all(name in summary for name in WORKLOADS):
        return None
    rule = {}
    for method, key in [
        ("sorting_direct", "direct_cr_over_sorting_direct"),
        ("rssa_cr", "direct_cr_over_rssa_cr"),
    ]:
        ratios = {name: summary[name][key] for name in WORKLOADS}
        winner = max(ratios, key=ratios.get)
        rule[method] = {
            "ratios": ratios,
            "largest": ratios[winner],
            "workload": winner,
            "ported": ratios[winner] >= SELECTOR_FACTOR,
        }
    return rule


def strip(row):
    return {key: value for key, value in row.items() if key != "timing"}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--python", type=Path, help="Python with the 0.6.0 wheel")
    parser.add_argument("--binary", type=Path, help="0.6.0 crn-craft CLI")
    parser.add_argument("--models", help="comma-separated model names (default: all)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check-against", type=Path)
    parser.add_argument("--worker", nargs=2, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(*args.worker)
        return 0
    if not (args.python and args.binary):
        parser.error("--python and --binary are required")
    binary = args.binary.resolve()
    models = load_manifest(args.models.split(",") if args.models else None)
    load = os.getloadavg()
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-merge-baseline-") as scratch:
        scratch = Path(scratch)
        identity = run_worker(args.python, __file__, {"kind": "identity"}, scratch)
        for model in models:
            name = model["name"]
            size = model["timing"]
            model_bytes = (ROOT / model["path"]).read_bytes()
            folder = scratch / name
            folder.mkdir()
            (folder / "model.json").write_bytes(model_bytes)
            requests = []
            for aggregator, waiting, fast in CONFIGS:
                request = ssa_request(
                    model_bytes,
                    aggregator=aggregator,
                    waiting_time=waiting,
                    fast=fast,
                    end=size["end"],
                    save_times=size["save_times"],
                    replicates=size["replicates"],
                    seed=size["seed"],
                    threads=size["threads"],
                )
                (folder / "request.json").write_bytes(encode(request))
                output, warmup, seconds = time_cli(binary, folder, REPETITIONS)
                result = json.loads(output)
                if result["termination"] != "completed":
                    raise RuntimeError(f"{name} {aggregator}: {result['diagnostics']}")
                rows.append(
                    {
                        "model": name,
                        "aggregator": aggregator,
                        "waiting_time": waiting,
                        "fast": fast,
                        "request": request,
                        "request_sha256": sha256(encode(request)),
                        "termination": result["termination"],
                        "completed_replicates": result["completed_replicates"],
                        "sampler": result["sampler"],
                        "events": result["events"],
                        "work": result["work"],
                        "result_sha256": canonical(result),
                        "cli_output_sha256": sha256(output),
                        "timing": {"cli": timing_entry(warmup, seconds)},
                    }
                )
                requests.append(request)
                print(
                    f"{name:18} {label(aggregator, waiting, fast):34}"
                    f" cli {statistics.median(seconds):8.4f} s"
                    f" events {result['events']}",
                    file=sys.stderr,
                    flush=True,
                )
            job = {
                "kind": "timing",
                "path": model["path"],
                "requests": requests,
                "repetitions": REPETITIONS,
            }
            answers = run_worker(args.python, __file__, job, scratch)
            for row, answer in zip(rows[-len(requests) :], answers, strict=True):
                if answer["result_sha256"] != row["result_sha256"]:
                    raise RuntimeError(f"{name} {row['aggregator']}: boundaries differ")
                row["python_output_sha256"] = answer["output_sha256"]
                row["timing"]["python"] = timing_entry(
                    answer["warmup_seconds"],
                    answer["seconds"],
                    preparation_seconds=answer["preparation_seconds"],
                )
    summary = {boundary: summarize(rows, boundary) for boundary in ["cli", "python"]}
    report = {
        "schema": "crn-craft.merge-baseline/0.1",
        "declared_in": "validation/merge/RESULTS.md",
        "models_json_sha256": sha256(MANIFEST.read_bytes()),
        "host": host_info(),
        "tools": {"crn_craft_cli": binary_identity(binary), "worker": identity},
        "repetitions": REPETITIONS,
        "boundaries": {
            "cli": "wall time of `crn-craft ssa model.json request.json` in a "
            "directory holding the model bytes and the request bytes "
            "(compact JSON plus newline); stdout to a file; gated",
            "python": "native simulate_json of crn_craft.prepare_ssa(model, "
            "request bytes); preparation and JSON decoding excluded; "
            "informative",
        },
        "command": COMMAND,
        "rows": rows,
        "timing_summary": {
            "load_average_at_start": list(load),
            "cli": summary["cli"],
            "python": summary["python"],
            "selector_rule": selector_rule(summary["cli"]),
            "selector_rule_python_boundary": selector_rule(summary["python"]),
        },
    }
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    if args.check_against:
        other = json.loads(args.check_against.read_text())
        known = {
            (r["model"], r["aggregator"], r["waiting_time"], r["fast"]): strip(r)
            for r in other["rows"]
        }
        different = [
            label(row["aggregator"], row["waiting_time"], row["fast"])
            + f" on {row['model']}"
            for row in rows
            if known.get(
                (row["model"], row["aggregator"], row["waiting_time"], row["fast"])
            )
            != strip(row)
        ]
        for key in ["host", "tools", "models_json_sha256", "command"]:
            if other[key] != report[key]:
                different.append(key)
        print(
            json.dumps({"compared": len(rows), "different": different}, indent=2),
            flush=True,
        )
        return 1 if different else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
