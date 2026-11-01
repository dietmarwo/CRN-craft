"""Independent matrix-exponential, analytic and DSMTS checks for native FSP."""

import argparse
import hashlib
import json
import math
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply

from check_cme import FIXTURES, generator


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--dsmts", type=Path)
    args = ap.parse_args()
    np.random.seed(271828)  # SciPy sparse norm estimator has randomized probes.
    rows = []
    fixtures = dict(FIXTURES)
    for name, initial, birth, death, bound in [
        ("DSMTS_001_01", 100, 0.1, 0.11, 600),
        ("DSMTS_002_01", 0, 1.0, 0.1, 150),
    ]:
        fixtures[name] = {
            "species": [("X", initial)],
            "parameters": {},
            "cap": bound,
            "reactions": [
                {
                    "name": "birth",
                    "reactants": {"X": 1} if initial else {},
                    "products": {"X": 2} if initial else {"X": 1},
                    "rate": {"kind": "mass_action", "coefficient": birth},
                },
                {
                    "name": "death",
                    "reactants": {"X": 1},
                    "products": {},
                    "rate": {"kind": "mass_action", "coefficient": death},
                },
            ],
        }
    with tempfile.TemporaryDirectory(prefix="crn-fsp-") as temp:
        folder = Path(temp)
        for name, fixture in fixtures.items():
            states, edges, exits = generator(fixture)
            n = len(states)
            q = np.zeros((n + 1, n + 1))
            for i, j, a in edges:
                q[i, j] += a
            for i, a in enumerate(exits):
                q[i, n] = a - q[i].sum()
                q[i, i] -= a
            model = {
                "schema": "crn-craft.model/0.1",
                "name": name,
                "convention": "stochastic_counts",
                "omega": 1.0,
                "parameters": fixture["parameters"],
                "species": [
                    {"name": s, "initial_count": c} for s, c in fixture["species"]
                ],
                "reactions": fixture["reactions"],
            }
            m = (json.dumps(model) + "\n").encode()
            (folder / "model.json").write_bytes(m)
            times = [0.0, 0.5, 2.0, 10.0]
            run = {
                "schema": "crn-craft.run/0.1",
                "model": {
                    "path": "model.json",
                    "sha256": hashlib.sha256(m).hexdigest(),
                },
                "method": {
                    "kind": "fsp",
                    "bounds": {
                        s: max(x[j] for x in states)
                        for j, (s, _) in enumerate(fixture["species"])
                    },
                    "max_states": 10000,
                    "tolerance": 1e-14,
                },
                "time": {"end": times[-1], "save_times": times},
                "outputs": {"statistics": ["mean"], "keep_trajectories": 0},
                "budget": {"threads": 1, "max_output_bytes": 10000000},
            }
            (folder / "run.json").write_text(json.dumps(run))
            p = subprocess.run(
                [
                    str(args.binary.resolve()),
                    "fsp",
                    str(folder / "model.json"),
                    str(folder / "run.json"),
                ],
                capture_output=True,
                check=False,
                timeout=120,
            )
            require(p.returncode == 0, p.stdout.decode() + p.stderr.decode())
            result = json.loads(p.stdout)
            idx = {tuple(x): i for i, x in enumerate(states)}
            errors = []
            for snapshot in result["snapshots"]:
                initial = np.zeros(n + 1)
                initial[0] = 1
                expected = (
                    expm(q * snapshot["time"])[0]
                    if n < 150
                    else expm_multiply(csc_matrix(q.T) * snapshot["time"], initial)
                )
                actual = np.zeros(n + 1)
                for x, value in zip(
                    result["states"], snapshot["probabilities"], strict=True
                ):
                    actual[idx[tuple(x)]] = value
                actual[n] = snapshot["physical_leakage_estimate"]
                errors.append(float(np.max(np.abs(actual - expected))))
                require(
                    snapshot["series_tail_upper_bound"] <= 1.01e-14, "tail allocation"
                )
                require(
                    abs(float(actual.sum()) + snapshot["discarded_weight_estimate"] - 1)
                    < 2e-12,
                    "mass balance",
                )
            rows.append(
                {
                    "name": name,
                    "check": "independent_matrix_exponential",
                    "passed": max(errors) <= 2e-12,
                    "max_error": max(errors),
                    "states": n,
                    "result_sha256": hashlib.sha256(p.stdout).hexdigest(),
                }
            )
            if name.startswith("DSMTS"):
                for snapshot in result["snapshots"]:
                    t = snapshot["time"]
                    mean = snapshot["projected_first_moment"][0]
                    second = snapshot["projected_second_raw_moment"][0]
                    if name == "DSMTS_001_01":
                        expected_mean = 100 * math.exp(-0.01 * t)
                        expected_var = (
                            100
                            * 0.21
                            / (-0.01)
                            * math.exp(-0.01 * t)
                            * math.expm1(-0.01 * t)
                        )
                    else:
                        expected_mean = 10 * (-math.expm1(-0.1 * t))
                        expected_var = expected_mean
                    rows.append(
                        {
                            "name": name,
                            "time": t,
                            "check": "analytic_moments",
                            "passed": abs(mean - expected_mean) < 1e-9
                            and abs(second - mean * mean - expected_var) < 1e-7,
                            "mean_error": abs(mean - expected_mean),
                            "variance_error": abs(second - mean * mean - expected_var),
                        }
                    )
    provenance = None
    if args.dsmts:
        names = [
            "dsmts-001-01.mod",
            "dsmts-001-01-mean.csv",
            "dsmts-001-01-sd.csv",
            "dsmts-002-01.mod",
            "dsmts-002-01-mean.csv",
            "dsmts-002-01-sd.csv",
        ]
        provenance = {
            "repository": "https://github.com/darrenjw/dsmts",
            "revision": subprocess.check_output(
                ["git", "-C", str(args.dsmts), "rev-parse", "HEAD"], text=True
            ).strip(),
            "license": "LGPL; external validation inputs only, not copied into the MIT engine",
            "hashes": {
                n: hashlib.sha256((args.dsmts / n).read_bytes()).hexdigest()
                for n in names
            },
        }
        for prefix in ["001", "002"]:
            values = np.genfromtxt(
                args.dsmts / f"dsmts-{prefix}-01-mean.csv", delimiter=",", skip_header=1
            )
            expected = (
                100 * np.exp(-0.01 * values[:, 0])
                if prefix == "001"
                else 10 * (-np.expm1(-0.1 * values[:, 0]))
            )
            rows.append(
                {
                    "name": prefix,
                    "check": "upstream_rounded_mean",
                    "passed": bool(np.max(np.abs(expected - values[:, 1])) < 5.1e-6),
                }
            )
    report = {
        "schema": "crn-craft.fsp-validation/0.1",
        "passed": all(r["passed"] for r in rows),
        "rows": rows,
        "dsmts": provenance,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    require(report["passed"], "FSP checks failed; see retained report")
    print(
        json.dumps(
            {
                "passed": report["passed"],
                "checks": len(rows),
                "max_distribution_error": max(r.get("max_error", 0) for r in rows),
            }
        )
    )


if __name__ == "__main__":
    main()
