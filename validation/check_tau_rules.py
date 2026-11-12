"""The approximate bias envelopes of check_tau.py on the rules `tau` method.

The network model of check_tau.py is converted with `crn-craft convert` and
run as a crn-craft.run/1.0 request with `{"kind": "tau", "epsilon": ε}` by
`crn-craft rules`, with the declared sizes, seeds and tolerances of
check_tau.py (validation/merge/RESULTS.md, M5c). check_tau.py keeps testing
the network engine's `tau` command until it is removed.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.stats import binom, poisson

SAMPLER = "chacha8-u53/rules-approximate-cgp-poisson/2"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    binary = str(args.binary.resolve())
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
    with tempfile.TemporaryDirectory(prefix="crn-tau-rules-") as temp:
        folder = Path(temp)
        network = {
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
        (folder / "network.json").write_text(json.dumps(network))
        converted = subprocess.run(
            [binary, "convert", str(folder / "network.json")],
            capture_output=True,
            check=False,
            timeout=60,
        )
        if converted.returncode:
            raise RuntimeError(converted.stdout.decode() + converted.stderr.decode())
        m = converted.stdout
        if json.loads(m)["schema"] != "crn-craft.model/1.0":
            raise RuntimeError("the converter did not write crn-craft.model/1.0")
        (folder / "model.json").write_bytes(m)
        for i, epsilon in enumerate([0.1, 0.05, 0.025]):
            method = {"kind": "tau", "epsilon": epsilon}
            request = {
                "schema": "crn-craft.run/1.0",
                "model": {
                    "path": "model.json",
                    "sha256": hashlib.sha256(m).hexdigest(),
                },
                "method": method,
                "ensemble": {"replicates": n, "root_seed": 9101 + i},
                "time": {"end": 2.0, "report_times": [0.0, 2.0]},
                "outputs": {
                    "statistics": ["mean", "variance"],
                    "keep_paths": n,
                    "keep_instances": 0,
                },
                "budget": {
                    "threads": 1,
                    "max_events_per_replicate": 10000,
                    "max_sampling_steps_per_replicate": 100000000,
                    "max_instances": 1000,
                    "max_output_bytes": 100000000,
                },
            }
            (folder / "tau.json").write_text(json.dumps(request))
            proc = subprocess.run(
                [binary, "rules", str(folder / "model.json"), str(folder / "tau.json")],
                capture_output=True,
                check=False,
                timeout=300,
            )
            if proc.returncode:
                raise RuntimeError(proc.stdout.decode() + proc.stderr.decode())
            result = json.loads(proc.stdout)
            column = result["observables"].index("X")
            samples = np.array(
                [int(path["values"][-1][column]) for path in result["paths"]]
            )
            if len(samples) != n:
                raise RuntimeError("missing paths")
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
                    "rng_words": sum(r["rng_words"] for r in result["replicates"]),
                    "sampler": result["sampler"],
                    "passed": error <= envelope
                    and bias <= mean_envelope
                    and result["termination"] == "completed"
                    and result["sampler"] == SAMPLER
                    and result["method"] == method,
                    "result_sha256": hashlib.sha256(proc.stdout).hexdigest(),
                }
            )
    report = {
        "schema": "crn-craft.tau-rules-validation/0.1",
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
