#!/usr/bin/env python3
"""P19 supplementary evidence: tau leaping with large leap means, two samplers.

Runs the W9 immigration-death model (stationary mean 10,000, epsilon 0.05,
horizon 5) with 4,000 paths under two binaries, the product-method sampler
(before) and PTRS (after), from different seeds. Both draw exact Poisson
variates for the same leap rule, so their final-count distributions must
agree; the two-sample Kolmogorov-Smirnov distance is compared with its
asymptotic critical value at alpha 0.001. The exact law of X(5) (binomial
survivors plus Poisson immigrants) is reported for context: tau leaping is
approximate, so its distance from the exact law is described, not gated.
This check was added after the declaration and is not one of P19's gates.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

from benchmark_iterations import tau_model

PATHS = 4000
HORIZON = 5.0


def run(binary, seed, folder):
    model = (json.dumps(tau_model()) + "\n").encode()
    (folder / "model.json").write_bytes(model)
    request = {
        "schema": "crn-craft.tau-request/0.1",
        "epsilon": 0.05,
        "max_leaps_per_replicate": 1000000,
        "reference_run": {
            "schema": "crn-craft.run/0.1",
            "model": {
                "path": "model.json",
                "sha256": hashlib.sha256(model).hexdigest(),
            },
            "method": {"kind": "ssa", "aggregator": "direct"},
            "ensemble": {"replicates": PATHS, "root_seed": seed},
            "time": {"end": HORIZON, "save_times": [0.0, HORIZON]},
            "outputs": {"statistics": ["mean", "variance"], "keep_trajectories": PATHS},
            "budget": {
                "threads": 1,
                "max_events_per_replicate": 100000000,
                "max_output_bytes": 100000000,
            },
        },
    }
    (folder / "tau.json").write_text(json.dumps(request))
    proc = subprocess.run(
        [str(binary), "tau", str(folder / "model.json"), str(folder / "tau.json")],
        capture_output=True,
        check=False,
        timeout=900,
    )
    if proc.returncode:
        raise RuntimeError(proc.stdout.decode() + proc.stderr.decode())
    result = json.loads(proc.stdout)
    finals = sorted(path["counts"][-1][0] for path in result["trajectories"])
    if len(finals) != PATHS:
        raise RuntimeError("missing trajectories")
    return result["algorithm"], finals, result.get("leaps")


def ks(a, b):
    i = j = 0
    distance = 0.0
    while i < len(a) and j < len(b):
        x = min(a[i], b[j])
        while i < len(a) and a[i] == x:
            i += 1
        while j < len(b) and b[j] == x:
            j += 1
        distance = max(distance, abs(i / len(a) - j / len(b)))
    return distance


def moments(values):
    mean = sum(values) / len(values)
    variance = sum((v - mean) ** 2 for v in values) / (len(values) - 1)
    return mean, variance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="crn-ptrs-") as temp:
        folder = Path(temp)
        old_algorithm, old, _ = run(args.before, 1101, folder)
        new_algorithm, new, _ = run(args.after, 2202, folder)
    distance = ks(old, new)
    critical = math.sqrt(-0.5 * math.log(0.001 / 2)) * math.sqrt(2 / PATHS)
    survive = math.exp(-HORIZON)
    exact_mean = 10000.0
    exact_variance = 10000.0 * survive * (1 - survive) + 10000.0 * (1 - survive)
    rows = {}
    for label, values, algorithm in [
        ("before", old, old_algorithm),
        ("after", new, new_algorithm),
    ]:
        mean, variance = moments(values)
        rows[label] = {
            "algorithm": algorithm,
            "mean": mean,
            "variance": variance,
            "mean_z_vs_exact": (mean - exact_mean) / math.sqrt(exact_variance / PATHS),
            "variance_ratio_vs_exact": variance / exact_variance,
        }
    report = {
        "schema": "crn-craft.tau-sampler-comparison/0.1",
        "paths_per_sampler": PATHS,
        "exact_mean": exact_mean,
        "exact_variance": exact_variance,
        "samplers": rows,
        "two_sample_ks": distance,
        "ks_critical_alpha_0_001": critical,
        "passed": distance <= critical,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps({"passed": report["passed"], "ks": distance, "critical": critical})
    )
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
