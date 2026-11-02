"""Binding/CLI equivalence, strict-input and held-out observable checks."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
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
    # Parameter search reuses topology while each candidate keeps full replay inputs.
    parameter_request = crn.ssa_spec(m, fast=True, end=2.0, replicates=8, keep=8)
    template = crn.prepare_parameters(m, parameter_request)
    values = [x * 1.2 for x in template.default_values]
    bound = template.bind(
        dict(zip(template.parameter_names, values, strict=True)), seed=5012
    )
    require(
        bound.simulate() == crn.simulate_ssa(bound.model_bytes, bound.request_bytes),
        "parameter binding replay",
    )
    require(
        template.simulate(values, seed=5012) == bound.simulate(),
        "repeated parameter binding",
    )
    try:
        template.bind([True] * len(values))
    except crn.CrnError:
        pass
    else:
        raise RuntimeError("boolean parameter was accepted")
    require(
        crn.mass_action_kernel_source(m).startswith("// Generated"),
        "AOT source interface",
    )
    rows.append({"check": "reusable_parameters_and_codegen", "passed": True})
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
        first["not_detected"] == 1
        and first["detected_fraction"] == 2 / 3
        and first["restricted_mean_upper"] == 1.0
        and "restricted_mean_lower" not in first,
        "saved-grid passage accounting",
    )
    # A crossing between saved times: the coarse grid misses it, so it must
    # neither report a passage nor a lower bound; the finer grid sees it.
    fine_t = [0.0, 0.5, 1.0, 1.5, 2.0]
    fine_x = [[0, 3, 1, 1, 0]]
    coarse = sampled_first_passage(fine_t[::2], [row[::2] for row in fine_x], 2)
    fine = sampled_first_passage(fine_t, fine_x, 2)
    require(
        coarse["detected"] == 0
        and fine["detected"] == 1
        and fine["restricted_mean_upper"] <= coarse["restricted_mean_upper"],
        "nested saved grids",
    )
    rows.append({"check": "saved_grid_first_passage", "passed": True})
    # A constant non-integer series is nonoscillatory, not an error.
    flat = oscillation(np.arange(100.0), np.full(100, 0.1))
    require(flat["classification"] == "nonoscillatory", "constant series")
    try:
        oscillation(np.arange(20.0), np.sin(np.arange(20.0)), burn_fraction=0.79)
    except ValueError:
        pass
    else:
        raise RuntimeError("too few samples after the burn-in were accepted")
    rows.append({"check": "oscillation_edge_cases", "passed": True})
    # The detector's three documented gates, each checked on its own:
    # removing the mean (an offset changes nothing), the autocorrelation gate
    # (a spectral peak whose lag-period autocorrelation is below 0.2 is not
    # called oscillatory) and the burn-in (samples before it are ignored).
    rng = np.random.default_rng(5)
    base = np.sin(2 * np.pi * t / 10) + 0.1 * rng.normal(size=len(t))
    plain, shifted = oscillation(t, base), oscillation(t, base + 1000.0)
    require(
        plain["classification"] == shifted["classification"] == "oscillatory"
        and abs(plain["autocorrelation_peak"] - shifted["autocorrelation_peak"]) < 1e-9,
        "oscillation is offset invariant",
    )
    noise = np.random.default_rng(12).normal(scale=np.sqrt(0.8), size=500)
    n = np.arange(500.0)
    aliased = oscillation(n, np.sin(2 * np.pi * 0.4 * n) + noise, burn_fraction=0)
    require(
        aliased["peak_fraction"] >= 0.3
        and aliased["cycles"] >= 3
        and aliased["autocorrelation_peak"] < 0.2
        and aliased["classification"] == "nonoscillatory",
        "autocorrelation gate",
    )
    transient = base.copy()
    transient[: int(len(t) * 0.2)] += 1000.0
    require(
        oscillation(t, transient) == plain
        and oscillation(t, transient, burn_fraction=0)["classification"]
        != "oscillatory",
        "burn-in excludes the transient",
    )
    rows.append({"check": "oscillation_gates", "passed": True})
    # Rejections carry their diagnostics; a cancel file stops a run; Ctrl-C
    # interrupts a long native call promptly instead of being ignored.
    model = crn.load_example("immigration-death")
    try:
        crn.simulate(model, dict(crn.run_spec(model), schema="crn-craft.run/9"))
    except crn.CrnError as error:
        require(
            isinstance(error, ValueError)
            and error.diagnostics[0]["code"] == "CRN-E003"
            and str(error).startswith("CRN-E003"),
            "CrnError diagnostics",
        )
    else:
        raise RuntimeError("an invalid run was accepted")
    with tempfile.TemporaryDirectory() as tmp:
        stop = Path(tmp) / "stop"
        stop.write_text("")
        long = crn.run_spec(model, end=1e6, replicates=8, max_events=10**9)
        cancelled = crn.simulate(model, long, cancel_file=stop)
        require(cancelled["termination"] == "cancelled", "binding cancel file")
    rows.append({"check": "binding_errors_and_cancel_file", "passed": True})
    if os.name == "posix":
        script = (
            "import crn_craft as c, time\n"
            "m = c.load_example('immigration-death')\n"
            "r = c.run_spec(m, end=1e9, replicates=4, max_events=10**12)\n"
            "print('started', flush=True)\n"
            "try:\n"
            "    c.simulate(m, r)\n"
            "except KeyboardInterrupt:\n"
            "    print('interrupted', flush=True)\n"
        )
        child = subprocess.Popen(
            [sys.executable, "-I", "-c", script],
            stdout=subprocess.PIPE,
            text=True,
        )
        require(child.stdout.readline().strip() == "started", "child started")
        time.sleep(0.5)
        sent = time.monotonic()
        child.send_signal(signal.SIGINT)
        out, _ = child.communicate(timeout=30)
        latency = time.monotonic() - sent
        require(
            out.strip() == "interrupted" and child.returncode == 0 and latency < 2.0,
            f"Ctrl-C interrupts a native run ({latency:.2f} s, {out!r})",
        )
        rows.append({"check": "keyboard_interrupt", "passed": True})
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
