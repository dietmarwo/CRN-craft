"""Predeclared approximate bias envelopes against an analytic linear law."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.stats import binom, poisson


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    n = 10000
    p = math.exp(-1.0)
    immigration = 40 * (1 - p)
    law = np.convolve(
        binom.pmf(np.arange(31), 30, p), poisson.pmf(np.arange(200), immigration)
    )
    cdf = law.cumsum()
    mean = 30 * p + immigration
    variance = 30 * p * (1 - p) + immigration
    epsilon_stat = math.sqrt(math.log(6 / 0.001) / (2 * n))
    with tempfile.TemporaryDirectory(prefix="crn-tau-") as temp:
        folder = Path(temp)
        model = {
            "schema": "crn-craft.model/0.1",
            "name": "TauControl",
            "convention": "stochastic_counts",
            "omega": 1.0,
            "parameters": {},
            "species": [{"name": "X", "initial_count": 30}],
            "reactions": [
                {
                    "name": "birth",
                    "reactants": {},
                    "products": {"X": 1},
                    "rate": {"kind": "mass_action", "coefficient": 20.0},
                },
                {
                    "name": "death",
                    "reactants": {"X": 1},
                    "products": {},
                    "rate": {"kind": "mass_action", "coefficient": 0.5},
                },
            ],
        }
        m = json.dumps(model).encode()
        (folder / "model.json").write_bytes(m)
        for i, epsilon in enumerate([0.1, 0.05, 0.025]):
            run = {
                "schema": "crn-craft.run/0.1",
                "model": {
                    "path": "model.json",
                    "sha256": hashlib.sha256(m).hexdigest(),
                },
                "method": {"kind": "ssa", "aggregator": "direct"},
                "ensemble": {"replicates": n, "root_seed": 9101 + i},
                "time": {"end": 2.0, "save_times": [0.0, 2.0]},
                "outputs": {"statistics": ["mean", "variance"], "keep_trajectories": n},
                "budget": {
                    "threads": 1,
                    "max_events_per_replicate": 10000,
                    "max_output_bytes": 100000000,
                },
            }
            request = {
                "schema": "crn-craft.tau-request/0.1",
                "epsilon": epsilon,
                "max_leaps_per_replicate": 100000,
                "reference_run": run,
            }
            (folder / "tau.json").write_text(json.dumps(request))
            proc = subprocess.run(
                [
                    str(args.binary.resolve()),
                    "tau",
                    str(folder / "model.json"),
                    str(folder / "tau.json"),
                ],
                capture_output=True,
                check=False,
                timeout=120,
            )
            if proc.returncode:
                raise RuntimeError(proc.stdout.decode() + proc.stderr.decode())
            result = json.loads(proc.stdout)
            samples = np.array([t["counts"][-1][0] for t in result["trajectories"]])
            empirical = (
                np.bincount(samples, minlength=len(cdf))[: len(cdf)].cumsum() / n
            )
            error = float(np.max(np.abs(empirical - cdf)))
            bias = abs(float(samples.mean()) - mean)
            envelope = epsilon / 2 + epsilon_stat
            mean_envelope = epsilon * mean + 4.5 * math.sqrt(variance / n)
            rows.append(
                {
                    "epsilon": epsilon,
                    "paths": n,
                    "cdf_error": error,
                    "cdf_envelope": envelope,
                    "mean_error": bias,
                    "mean_envelope": mean_envelope,
                    "min_count": int(samples.min()),
                    "events": result["events"],
                    "leaps": result["leaps"],
                    "rejected_leaps": result["rejected_leaps"],
                    "passed": error <= envelope
                    and bias <= mean_envelope
                    and result["assessment"] == "approximate_ensemble"
                    and result["termination"] == "completed",
                    "result_sha256": hashlib.sha256(proc.stdout).hexdigest(),
                }
            )
    report = {
        "schema": "crn-craft.tau-validation/0.1",
        "passed": all(r["passed"] for r in rows),
        "epsilon_stat": epsilon_stat,
        "rows": rows,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["passed"]:
        raise RuntimeError("approximate gate failed; retained report")
    print(
        json.dumps(
            {
                "passed": True,
                "paths": 30000,
                "cdf_errors": [r["cdf_error"] for r in rows],
            }
        )
    )


if __name__ == "__main__":
    main()
