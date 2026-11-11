"""Shared helpers of the 0.7.0 merge validation scripts (standard library only).

`references.py` and `baseline.py` (task M0c) use them to record the 0.6.0
engine; the law-gate and speed scripts of later tasks reuse the request
builder, the statistics and the comparisons, so every engine is measured with
the same formulas. The declarations they implement are in `RESULTS.md`.
"""

import gzip
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "models.json"
# Standard errors per comparison in the law gate (RESULTS.md, "Law gate").
SIGMAS = 5


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def encode(value):
    """Compact JSON bytes with a trailing newline, as `crn_craft.model_bytes`."""
    return (json.dumps(value, allow_nan=False, separators=(",", ":")) + "\n").encode()


def gzip_bytes(data):
    """Deterministic gzip: no file name and a zero time stamp."""
    return gzip.compress(data, compresslevel=9, mtime=0)


def load_manifest(names=None, path=MANIFEST):
    """Rows of `models.json` (or `path`), optionally only `names`, hashes checked."""
    manifest = json.loads(Path(path).read_text())
    rows = manifest["rows"]
    if names:
        known = {row["name"] for row in rows}
        unknown = sorted(set(names) - known)
        if unknown:
            raise SystemExit(f"unknown models: {unknown}")
        rows = [row for row in rows if row["name"] in names]
    for row in rows:
        data = (ROOT / row["path"]).read_bytes()
        if sha256(data) != row["sha256"]:
            raise SystemExit(f"{row['path']}: model bytes differ from models.json")
        times = row["report_times"]
        if times[0] != 0.0 or times[-1] != row["end"] or times != sorted(set(times)):
            raise SystemExit(f"{row['name']}: report times must rise from 0 to end")
        if row["replicates"] % row["chunk"]:
            raise SystemExit(f"{row['name']}: replicates must be a multiple of chunk")
    return rows


def ssa_request(
    model_bytes,
    *,
    aggregator,
    waiting_time,
    fast,
    end,
    save_times,
    replicates,
    seed,
    threads,
    keep=0,
    max_events=10**9,
    max_output_bytes=200_000_000,
    max_sampling_steps=10**9,
):
    """A 0.6.0 `ssa-run/0.3` request (0.4 for `next_reaction`), as `crn.ssa_spec`."""
    return {
        "schema": "crn-craft.ssa-run/0.4"
        if aggregator == "next_reaction"
        else "crn-craft.ssa-run/0.3",
        "model": {"path": "model.json", "sha256": sha256(model_bytes)},
        "method": {
            "kind": "ssa",
            "aggregator": aggregator,
            "waiting_time": waiting_time,
            "fast": fast,
        },
        "time": {"end": end, "save_times": list(save_times)},
        "ensemble": {"replicates": replicates, "root_seed": seed},
        "outputs": {"statistics": ["mean", "variance"], "keep_trajectories": keep},
        "budget": {
            "threads": threads,
            "max_events_per_replicate": max_events,
            "max_output_bytes": max_output_bytes,
            "max_sampling_steps_per_replicate": max_sampling_steps,
        },
    }


def host_info():
    """Static facts about the machine; equal for reruns on the same host."""

    def first(path, prefix):
        try:
            for line in Path(path).read_text().splitlines():
                if line.startswith(prefix):
                    return line.split(":", 1)[1].strip()
        except OSError:
            pass
        return None

    try:
        quota = Path("/sys/fs/cgroup/cpu.max").read_text().strip()
    except OSError:
        quota = None
    return {
        "cpu_model": first("/proc/cpuinfo", "model name"),
        "logical_cpus": os.cpu_count(),
        "cpu_affinity": len(os.sched_getaffinity(0)),
        "cpu_quota": quota,
        "memory_total": first("/proc/meminfo", "MemTotal"),
        "platform": platform.platform(),
    }


def binary_identity(binary):
    data = Path(binary).read_bytes()
    out = subprocess.run(
        [str(binary), "--version"], capture_output=True, check=True
    ).stdout
    return {"sha256": sha256(data), "version": json.loads(out)}


def worker_command(python, script, *arguments):
    """Start `script` in worker mode under `python`, passing on `-O`."""
    command = [str(python)]
    if sys.flags.optimize:
        command.append("-" + "O" * sys.flags.optimize)
    return command + [str(script), *map(str, arguments)]


def run_worker(python, script, job, scratch):
    """Run one worker job (a JSON value) and return its JSON answer."""
    job_path = Path(scratch) / "job.json"
    answer_path = Path(scratch) / "answer.json"
    job_path.write_text(json.dumps(job))
    subprocess.run(
        worker_command(python, script, "--worker", job_path, answer_path),
        check=True,
        cwd=ROOT,
    )
    return json.loads(answer_path.read_text())


def cell_statistics(n, s1, s2, s3, s4):
    """Mean, variance and standard errors of n integer samples from power sums.

    Exact rational arithmetic, rounded once to float (RESULTS.md, "Reference
    method and statistics"): unbiased variance v, fourth central moment
    m4 = sum((x - mean)^4)/n, se = sqrt(v/n) and the fourth-moment estimate
    sev = sqrt(max(0, (m4 - v^2 (n-3)/(n-1))/n)).
    """
    mean = Fraction(s1, n)
    variance = (s2 - mean * s1) / (n - 1)
    m4 = (s4 - 4 * mean * s3 + 6 * mean**2 * s2 - 3 * n * mean**4) / n
    var_of_variance = (m4 - variance**2 * Fraction(n - 3, n - 1)) / n
    return {
        "mean": float(mean),
        "variance": float(variance),
        "mean_se": math.sqrt(float(variance / n)),
        "variance_se": math.sqrt(float(max(var_of_variance, Fraction(0)))),
        "fourth_central_moment": float(m4),
    }


def statistics_from_power_sums(power_sums):
    """Per-cell statistics (lists indexed [time][species]) from stored power sums."""
    n = power_sums["n"]
    names = ["mean", "variance", "mean_se", "variance_se", "fourth_central_moment"]
    out = {name: [] for name in names}
    for row in zip(*(power_sums[k] for k in ["s1", "s2", "s3", "s4"]), strict=True):
        cells = [cell_statistics(n, *cell) for cell in zip(*row, strict=True)]
        for name in names:
            out[name].append([cell[name] for cell in cells])
    return out


def _row(kind, time, species, value, reference, scale):
    bound = SIGMAS * scale
    difference = abs(value - reference)
    return {
        "comparison": kind,
        "time": time,
        "species": species,
        "value": value,
        "reference": reference,
        "bound": bound,
        "z": difference / scale if scale > 0 else (0.0 if difference == 0 else None),
        "passed": difference <= bound,
    }


def compare_to_reference(new, reference, times, species):
    """Law-gate comparisons 1 and 2: means and variances of two samples."""
    rows = []
    for t, time in enumerate(times):
        for s, name in enumerate(species):
            for kind, value, se in [
                ("mean", "mean", "mean_se"),
                ("variance", "variance", "variance_se"),
            ]:
                scale = math.hypot(new[se][t][s], reference[se][t][s])
                rows.append(
                    _row(
                        kind,
                        time,
                        name,
                        new[value][t][s],
                        reference[value][t][s],
                        scale,
                    )
                )
    return rows


def compare_to_law(new, law, n, times, species):
    """Law-gate comparisons 3 and 4: sample means and variances against an exact law."""
    rows = []
    for t, time in enumerate(times):
        for s, name in enumerate(species):
            mu = law["mean"][t][s]
            sigma2 = law["variance"][t][s]
            mu4 = law["fourth_central_moment"][t][s]
            rows.append(
                _row(
                    "fsp_mean",
                    time,
                    name,
                    new["mean"][t][s],
                    mu,
                    math.sqrt(sigma2 / n),
                )
            )
            spread = max(0.0, mu4 - sigma2 * sigma2 * (n - 3) / (n - 1))
            rows.append(
                _row(
                    "fsp_variance",
                    time,
                    name,
                    new["variance"][t][s],
                    sigma2,
                    math.sqrt(spread / n),
                )
            )
    return rows


def summarize(rows):
    zs = [r["z"] for r in rows if r["z"] is not None]
    return {
        "comparisons": len(rows),
        "failures": sum(not r["passed"] for r in rows),
        "max_z": max(zs) if zs else None,
    }
