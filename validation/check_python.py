"""Binding/CLI equivalence, strict-input and held-out observable checks."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import crn_craft as crn
from crn_craft.observables import oscillation, sampled_first_passage


def require(value, message):
    if not value:
        raise RuntimeError(message)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    rows = []
    m = (root / "models/examples/immigration-death.model.json").read_bytes()
    run = crn.run_spec(m, end=2.0, replicates=64, seed=5001, threads=1, keep=64)
    expected = crn.simulate(m, run)
    for threads in [1, 2, 8]:
        run["budget"]["threads"] = threads
        actual = crn.simulate(m, run)
        normalized = dict(actual)
        normalized.pop("run_sha256")
        exp = dict(expected)
        exp.pop("run_sha256")
        require(normalized == exp, "binding thread identity")
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            (folder / "model.json").write_bytes(m)
            (folder / "run.json").write_bytes(crn.model_bytes(run))
            proc = subprocess.run(
                [
                    str(args.binary.resolve()),
                    "simulate",
                    str(folder / "model.json"),
                    str(folder / "run.json"),
                    str(folder / "out"),
                ],
                capture_output=True,
                check=False,
            )
            require(proc.returncode == 0, proc.stdout.decode())
            require(
                json.loads((folder / "out/result.json").read_bytes()) == actual,
                "CLI binding identity",
            )
        rows.append(
            {
                "check": f"CLI_binding_threads_{threads}",
                "passed": True,
                "numerical_sha256": hashlib.sha256(
                    json.dumps(normalized, sort_keys=True).encode()
                ).hexdigest(),
            }
        )
    for invalid in [b'{"name":"a","name":"b"}', b'{"schema":"crn-craft.model/99"}']:
        require(not crn.check(invalid)["valid"], "strict binding parser")
    rows.append({"check": "strict_inputs", "passed": True})
    t = np.linspace(0, 100, 501)
    for split, seeds in [
        ("development", [101, 102]),
        ("held_out", [201, 202, 203, 204]),
    ]:
        correct = 0
        for seed in seeds:
            rng = np.random.default_rng(seed)
            positive = np.sin(2 * np.pi * t / 10) + 0.1 * rng.normal(size=len(t))
            negative = rng.normal(size=len(t))
            correct += oscillation(t, positive)["classification"] == "oscillatory"
            correct += oscillation(t, negative)["classification"] == "nonoscillatory"
        require(correct == 2 * len(seeds), "synthetic held-out detector checks")
        rows.append(
            {
                "check": "synthetic_" + split,
                "passed": True,
                "correct": correct,
                "total": 2 * len(seeds),
                "scope": "labelled synthetic controls; no biological sensitivity/specificity claim",
            }
        )
    first = sampled_first_passage([0, 1, 2], [[0, 2, 3], [0, 0, 0], [3, 3, 3]], 2)
    require(
        first["censored"] == 1
        and first["restricted_mean_lower"] == 2 / 3
        and first["restricted_mean_upper"] == 1.0,
        "censor accounting",
    )
    rows.append({"check": "interval_censoring", "passed": True})
    circuit = (root / "models/examples/repressilator.model.json").read_bytes()
    study = []
    for split, seeds in [("development", [3101, 3102]), ("held_out", [3201, 3202])]:
        for seed in seeds:
            run = crn.run_spec(
                circuit,
                end=100.0,
                save_times=np.linspace(0, 100, 501).tolist(),
                replicates=4,
                seed=seed,
                keep=4,
            )
            result = crn.simulate(circuit, run)
            require(result["termination"] == "completed", "circuit incomplete")
            diagnostics = [
                oscillation(result["save_times"], np.asarray(path["counts"])[:, 0])
                for path in result["trajectories"]
            ]
            study.append(
                {
                    "split": split,
                    "seed": seed,
                    "events": result["events"],
                    "classification": [d["classification"] for d in diagnostics],
                    "scope": "synthetic fixed topology illustration; circuit-study source not supplied",
                }
            )
    args.output.write_text(
        json.dumps(
            {
                "schema": "crn-craft.python-validation/0.1",
                "passed": True,
                "rows": rows,
                "circuit_study": study,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"passed": True, "checks": len(rows), "circuit_paths": 16}))


if __name__ == "__main__":
    main()
