"""Law gate of the 0.7.0 engine merge (M2, later M4a, M4b, M4c and after M8).

For every model of `validation/merge/models.json` this converts the network
model with the engine under test (`crn-craft convert`, `crn_craft.convert`),
runs the converted rules model with the given method on the declared grid
(`end`, `report_times`) with the declared replicates in chunks, chunk k with
root seed `gate_seed + k`, and accumulates exact integer power sums of every
species count at every report time. It then applies the law gate declared in
`validation/merge/RESULTS.md` with the formulas of `merge_common.py`:

1. mean against the 0.6.0 reference: |x - x_ref| <= 5 sqrt(se^2 + se_ref^2);
2. variance against the reference, with fourth-moment standard errors;
3. and 4. mean and variance against the FSP law where one is stored.

Species are matched to the converted model's observables by name. A model the
converter refuses is listed with its diagnostics; `--expect-refused` names the
models that must be refused (for M2: W4, W5, W7), and every other model must
convert. Every comparison row is written to the JSON report, failed rows
included; the report holds no timings, so reruns compare byte for byte.

The engine under test runs either through Python (`--python`, an interpreter
with the `crn_craft` wheel under test and numpy; the simulations run in a
worker process) or through the native CLI (`--binary`; this process then
needs only the standard library). The method is a JSON object, so other
methods and selectors reuse the script:

    flock /home/claude/craft5/cpu.lock python validation/merge/check_laws.py \\
        --python /home/claude/fb2/engine-venv/bin/python \\
        --method '{"kind": "static_direct"}' --expect-refused W4,W5,W7

`--replicates N` makes a smoke run with N replicates per model (chunks of at
most the declared size); its report says so and does not count as a gate run.
Exit status: 0 when every compared model passes and the refusals are as
expected, 1 otherwise.
"""

import argparse
import base64
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_common import (  # noqa: E402
    HERE,
    MANIFEST,
    ROOT,
    binary_identity,
    compare_to_law,
    compare_to_reference,
    encode,
    load_manifest,
    run_worker,
    sha256,
    statistics_from_power_sums,
    summarize,
)

REFERENCE = HERE / "reference"
DEFAULT_METHOD = {"kind": "static_direct"}
# The chunk requests: generous work limits, so a replicate never stops early.
MAX_EVENTS = 10**9
MAX_SAMPLING_STEPS = 10**12
MAX_OUTPUT_BYTES = 1 << 30


def dump(value):
    return (json.dumps(value, separators=(",", ":"), allow_nan=False) + "\n").encode()


def request_for(converted, row, *, chunk, seed, threads, method):
    """The run request of one chunk: every path kept, mean and variance."""
    schema = (
        "crn-craft.run/1.0"
        if converted["schema"] == "crn-craft.model/1.0"
        else "crn-craft.rules-run/0.6"
    )
    return {
        "schema": schema,
        "model": {"path": "model.json", "sha256": None},
        "time": {"end": row["end"], "report_times": row["report_times"]},
        "ensemble": {"replicates": chunk, "root_seed": seed},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": chunk,
            "keep_instances": 0,
        },
        "budget": {
            "max_events_per_replicate": MAX_EVENTS,
            "max_sampling_steps_per_replicate": MAX_SAMPLING_STEPS,
            "max_instances": max(1, len(converted["types"])),
            "threads": threads,
            "max_output_bytes": MAX_OUTPUT_BYTES,
        },
        "method": method,
    }


def checked_paths(result, chunk, times, name):
    """The kept series of a complete chunk as [replicate][time][observable]."""
    paths = result["paths"]
    if (
        result["termination"] != "completed"
        or result["completed_replicates"] != chunk
        or result["stopped_replicates"] != 0
        or [p["id"] for p in paths] != list(range(chunk))
        or any(len(p["values"]) != len(times) for p in paths)
    ):
        raise RuntimeError(f"{name}: incomplete chunk {result['diagnostics']}")
    return [p["values"] for p in paths]


def chunk_plan(row, replicates):
    chunk = min(row["chunk"], replicates)
    if replicates % chunk:
        raise SystemExit(f"{row['name']}: {replicates} replicates are not whole chunks")
    return chunk, replicates // chunk


# ---------------------------------------------------------------- worker side


def worker_model(job):
    """Convert and run one model through the Python API; return power sums."""
    import crn_craft as crn
    import numpy as np
    from crn_craft import rules
    from references import power_sums

    network = (ROOT / job["path"]).read_bytes()
    try:
        converted = crn.convert(network)
    except crn.CrnError as error:
        return {"refused": error.diagnostics}
    document = json.loads(converted)
    totals = None
    chunks = []
    sampler = None
    columns = None
    for k in range(job["chunks"]):
        seed = job["gate_seed"] + k
        request = request_for(
            document,
            job,
            chunk=job["chunk"],
            seed=seed,
            threads=job["threads"],
            method=job["method"],
        )
        request["model"]["sha256"] = hashlib.sha256(converted).hexdigest()
        result = rules.prepare(converted, request).simulate().metadata
        values = checked_paths(result, job["chunk"], job["report_times"], job["name"])
        if k == 0:
            sampler = result["sampler"]
            observables = result["observables"]
            missing = [s for s in job["species"] if s not in observables]
            if missing:
                raise RuntimeError(
                    f"{job['name']}: species without observable {missing}"
                )
            columns = [observables.index(s) for s in job["species"]]
        if result["sampler"] != sampler:
            raise RuntimeError(f"{job['name']}: sampler changed")
        x = np.asarray(values, dtype=np.float64)[:, :, columns]
        counts = x.astype(np.int64)
        if not np.array_equal(counts, x):
            raise RuntimeError(f"{job['name']}: non-integral counts")
        # Decoding check against the engine's own statistics.
        stats = result["statistics"]
        for mine, theirs in [
            (x.mean(axis=0), stats["mean"]),
            (x.var(axis=0, ddof=1), stats["sample_variance"]),
        ]:
            theirs = np.asarray(theirs)[:, columns]
            if not np.all(
                np.abs(mine - theirs) <= 1e-8 * np.maximum(1.0, np.abs(theirs))
            ):
                raise RuntimeError(f"{job['name']}: statistics disagree with the paths")
        sums = power_sums(counts, np)
        totals = (
            sums
            if totals is None
            else [a + b for a, b in zip(totals, sums, strict=True)]
        )
        chunks.append({"root_seed": seed, "events": result["events"]})
        print(
            f"  {job['name']} chunk {k} seed {seed} events {result['events']}",
            file=sys.stderr,
            flush=True,
        )
    return {
        "converted": base64.b64encode(converted).decode(),
        "request_schema": request["schema"],
        "sampler": sampler,
        "chunks": chunks,
        "power_sums": {
            "n": job["chunk"] * job["chunks"],
            **{f"s{p + 1}": totals[p].tolist() for p in range(4)},
        },
    }


def worker_identity():
    import crn_craft as crn
    import numpy as np
    from crn_craft import _native

    return {
        "python": sys.version,
        "crn_craft": crn.__version__,
        "native_sha256": sha256(Path(_native.__file__).read_bytes()),
        "numpy": np.__version__,
    }


def worker(job_path, answer_path):
    job = json.loads(Path(job_path).read_text())
    answer = worker_identity() if job["kind"] == "identity" else worker_model(job)
    Path(answer_path).write_text(json.dumps(answer))


# ------------------------------------------------------------ CLI backend


def cli_model(binary, job, scratch):
    """Convert and run one model through the native CLI; return power sums."""
    folder = Path(scratch) / job["name"]
    folder.mkdir()
    done = subprocess.run(
        [str(binary), "convert", str(ROOT / job["path"])], capture_output=True
    )
    if done.returncode != 0:
        return {"refused": json.loads(done.stdout)["diagnostics"]}
    converted = done.stdout
    document = json.loads(converted)
    (folder / "model.json").write_bytes(converted)
    times = job["report_times"]
    shape = (len(times), len(job["species"]))
    totals = [[[0] * shape[1] for _ in range(shape[0])] for _ in range(4)]
    chunks = []
    sampler = None
    columns = None
    for k in range(job["chunks"]):
        seed = job["gate_seed"] + k
        request = request_for(
            document,
            job,
            chunk=job["chunk"],
            seed=seed,
            threads=job["threads"],
            method=job["method"],
        )
        request["model"]["sha256"] = sha256(converted)
        (folder / "request.json").write_bytes(encode(request))
        done = subprocess.run(
            [str(binary), "rules", "model.json", "request.json"],
            capture_output=True,
            cwd=folder,
        )
        result = json.loads(done.stdout)
        if "termination" not in result:
            raise RuntimeError(f"{job['name']}: {result}")
        values = checked_paths(result, job["chunk"], times, job["name"])
        if k == 0:
            sampler = result["sampler"]
            observables = result["observables"]
            missing = [s for s in job["species"] if s not in observables]
            if missing:
                raise RuntimeError(
                    f"{job['name']}: species without observable {missing}"
                )
            columns = [observables.index(s) for s in job["species"]]
        if result["sampler"] != sampler:
            raise RuntimeError(f"{job['name']}: sampler changed")
        stats = result["statistics"]
        n = job["chunk"]
        for t in range(shape[0]):
            for s, column in enumerate(columns):
                counts = []
                for path in values:
                    value = path[t][column]
                    if value != int(value):
                        raise RuntimeError(f"{job['name']}: non-integral count")
                    counts.append(int(value))
                power = 1
                powers = [0, 0, 0, 0]
                for x in counts:
                    power = x
                    for p in range(4):
                        powers[p] += power
                        power *= x
                for p in range(4):
                    totals[p][t][s] += powers[p]
                # Decoding check against the engine's own statistics.
                mean = math.fsum(counts) / n
                variance = (
                    math.fsum((x - mean) ** 2 for x in counts) / (n - 1)
                    if n > 1
                    else 0.0
                )
                for mine, theirs in [
                    (mean, stats["mean"][t][column]),
                    (variance, stats["sample_variance"][t][column]),
                ]:
                    if abs(mine - theirs) > 1e-8 * max(1.0, abs(theirs)):
                        raise RuntimeError(
                            f"{job['name']}: statistics disagree with the paths"
                        )
        chunks.append({"root_seed": seed, "events": result["events"]})
        print(
            f"  {job['name']} chunk {k} seed {seed} events {result['events']}",
            file=sys.stderr,
            flush=True,
        )
    return {
        "converted": base64.b64encode(converted).decode(),
        "request_schema": request["schema"],
        "sampler": sampler,
        "chunks": chunks,
        "power_sums": {
            "n": job["chunk"] * job["chunks"],
            **{f"s{p + 1}": totals[p] for p in range(4)},
        },
    }


# ----------------------------------------------------------------- outer side


def compare(row, answer, n):
    """The law-gate rows of one model and their summaries."""
    reference = json.loads((REFERENCE / f"{row['name']}.moments.json").read_text())
    species = reference["species"]
    times = row["report_times"]
    if reference["report_times"] != times or reference["end"] != row["end"]:
        raise RuntimeError(f"{row['name']}: reference grid differs from models.json")
    stats = statistics_from_power_sums(answer["power_sums"])
    rows = compare_to_reference(stats, reference, times, species)
    out = {"reference": summarize(rows)}
    law_path = REFERENCE / f"{row['name']}.fsp-moments.json"
    if law_path.exists():
        law = json.loads(law_path.read_text())
        if law["species"] != species or law["report_times"] != times:
            raise RuntimeError(f"{row['name']}: FSP law grid differs")
        law_rows = compare_to_law(stats, law, n, times, species)
        out["fsp"] = summarize(law_rows)
        rows += law_rows
    out["all"] = summarize(rows)
    return out, rows, stats


def markdown(report):
    """Result rows for RESULTS.md."""
    lines = [
        "| Model | Cells | Comparisons (reference) | Failures | max \\|z\\| "
        "| Comparisons (FSP) | Failures | max \\|z\\| |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    def z(value):
        return "–" if value is None else f"{value:.2f}"

    for model in report["models"]:
        if model["status"] == "refused":
            codes = sorted({d["code"] for d in model["diagnostics"]})
            lines.append(
                f"| {model['name']} | refused ({', '.join(codes)}) | | | | | | |"
            )
            continue
        ref = model["summary"]["reference"]
        fsp = model["summary"].get("fsp")
        lines.append(
            f"| {model['name']} | {model['cells']} | {ref['comparisons']} "
            f"| {ref['failures']} | {z(ref['max_z'])} "
            f"| {fsp['comparisons'] if fsp else '–'} "
            f"| {fsp['failures'] if fsp else '–'} | {z(fsp['max_z']) if fsp else '–'} |"
        )
    total = report["summary"]
    lines.append(
        f"| **Total** | {total['cells']} | {total['comparisons']} comparisons "
        f"(N), {total['deterministic']} at t = 0 | **{total['failures']}** "
        f"| {z(total['max_z'])} | | | |"
    )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    engine = parser.add_mutually_exclusive_group()
    engine.add_argument("--python", type=Path, help="Python with the wheel under test")
    engine.add_argument("--binary", type=Path, help="crn-craft CLI under test")
    parser.add_argument(
        "--method",
        default=json.dumps(DEFAULT_METHOD),
        help='the run request\'s method object, default {"kind": "static_direct"}',
    )
    parser.add_argument("--label", help="report label (default: the method kind)")
    parser.add_argument("--models", help="comma-separated model names (default: all)")
    parser.add_argument(
        "--expect-refused", default="", help="models the converter must refuse"
    )
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--replicates", type=int, help="smoke run with N replicates")
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--manifest", type=Path, default=MANIFEST, help=argparse.SUPPRESS
    )
    parser.add_argument("--worker", nargs=2, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(*args.worker)
        return 0
    if not (args.python or args.binary):
        parser.error("one of --python and --binary is required")
    method = json.loads(args.method)
    if not isinstance(method, dict) or "kind" not in method:
        parser.error("--method must be a JSON object with a kind")
    label = args.label or method["kind"]
    expect_refused = {name for name in args.expect_refused.split(",") if name}
    rows = load_manifest(args.models.split(",") if args.models else None, args.manifest)
    unknown = expect_refused - {row["name"] for row in rows}
    if unknown:
        parser.error(f"--expect-refused names models not in the run: {sorted(unknown)}")
    output = args.output or (
        HERE / "laws" / f"{label}{'-smoke' if args.replicates else ''}.json"
    )
    report = {
        "schema": "crn-craft.merge-law-gate/0.1",
        "declared_in": "validation/merge/RESULTS.md",
        "models_json_sha256": sha256(args.manifest.read_bytes()),
        "label": label,
        "method": method,
        "smoke": args.replicates is not None,
        "backend": "python" if args.python else "cli",
        "engine": None,
        "models": [],
    }
    passed = True
    with tempfile.TemporaryDirectory(prefix="crn-merge-laws-") as scratch:
        if args.python:
            report["engine"] = run_worker(
                args.python, __file__, {"kind": "identity"}, Path(scratch)
            )
        else:
            report["engine"] = binary_identity(args.binary.resolve())
        for row in rows:
            name = row["name"]
            replicates = args.replicates or row["replicates"]
            chunk, chunks = chunk_plan(row, replicates)
            reference = json.loads((REFERENCE / f"{name}.moments.json").read_text())
            job = {
                "kind": "model",
                **row,
                "species": reference["species"],
                "chunk": chunk,
                "chunks": chunks,
                "threads": args.threads,
                "method": method,
            }
            print(f"{name}: {replicates} replicates", file=sys.stderr, flush=True)
            if args.python:
                answer = run_worker(args.python, __file__, job, Path(scratch))
            else:
                answer = cli_model(args.binary.resolve(), job, scratch)
            entry = {"name": name, "path": row["path"], "sha256": row["sha256"]}
            if "refused" in answer:
                entry["status"] = "refused"
                entry["diagnostics"] = answer["refused"]
                entry["expected"] = name in expect_refused
                passed &= entry["expected"]
                print(f"{name}: refused", file=sys.stderr, flush=True)
                report["models"].append(entry)
                continue
            converted = base64.b64decode(answer["converted"])
            summaries, comparisons, stats = compare(row, answer, replicates)
            entry.update(
                {
                    "status": "compared",
                    "expected": name not in expect_refused,
                    "converted_sha256": sha256(converted),
                    "request_schema": answer["request_schema"],
                    "sampler": answer["sampler"],
                    "replicates": replicates,
                    "chunk": chunk,
                    "chunks": answer["chunks"],
                    "events": sum(c["events"] for c in answer["chunks"]),
                    "species": reference["species"],
                    "cells": len(row["report_times"]) * len(reference["species"]),
                    "summary": summaries,
                    "statistics": {
                        key: stats[key]
                        for key in ["mean", "variance", "mean_se", "variance_se"]
                    },
                    "power_sums": answer["power_sums"],
                    "rows": comparisons,
                }
            )
            passed &= entry["expected"] and summaries["all"]["failures"] == 0
            print(f"{name}: {summaries}", file=sys.stderr, flush=True)
            report["models"].append(entry)
    compared = [m for m in report["models"] if m["status"] == "compared"]
    every_row = [r for m in compared for r in m["rows"]]
    report["summary"] = {
        **summarize(every_row),
        "cells": sum(m["cells"] for m in compared),
        "deterministic": sum(1 for r in every_row if r["time"] == 0.0),
        "models_compared": len(compared),
        "models_refused": len(report["models"]) - len(compared),
        "passed": passed,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(dump(report))
    print(markdown(report))
    print(json.dumps(report["summary"], indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
