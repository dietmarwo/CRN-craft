"""Exact-law check of the Static Direct selectors (M4c), safe under Python -O.

The structure of `validation/check_alternative_ssa.py`, on the unified engine:
each fixture is a network model converted with the CLI under test
(`crn-craft convert`) and run as a `crn-craft.run/1.0` Static Direct request
with every selector (`linear`, `composition_rejection`, `sorting`,
`rssa_cr`) and both waiting times (`inverse`, `ziggurat`), n = 100,000
replicates with every path kept, through `crn-craft rules`. At every positive
report time, the largest CDF distance between the replicates' marginal and
the exact law must be at most the DKW bound of the declared family:
epsilon = sqrt(ln(2 * TESTS / ALPHA) / (2 n)). The reported statistics must
equal the moments of the kept paths, and for the fixtures with a wrong law
(the falling-factorial convention) the check must reject that law at some
time (power).

Fixtures: the immigration-death laws of `validation/ssa-fixtures` (pure
decay, immigration from zero and from twelve; Poisson-binomial laws in closed
form) and the CME fixtures of `validation/check_cme.py` (dimerization,
heterodimer, Michaelis-Menten with Hill recycling, open annihilation; laws by
uniformization), with their root seeds. Declared in
`validation/merge/RESULTS.md` before the first run.
"""

import argparse
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation"))
from check_cme import (  # noqa: E402
    FIXTURES,
    TIMES,
    generator,
    marginal,
    max_cdf_error,
    uniformize,
)
from check_ssa import analytic_probabilities, cdf_distance  # noqa: E402

N = 100_000
SELECTORS = ["linear", "composition_rejection", "sorting", "rssa_cr"]
WAITS = ["inverse", "ziggurat"]
ALPHA = 0.001
TESTS = 168
EPSILON = math.sqrt(math.log(2 * TESTS / ALPHA) / (2 * N))


def fixtures():
    for name, initial, birth, death in [
        ("pure-decay", 20, 0.0, 0.7),
        ("immigration-zero", 0, 3.0, 0.8),
        ("immigration-nonzero", 12, 3.0, 0.8),
    ]:
        model = (ROOT / f"validation/ssa-fixtures/{name}.model.json").read_bytes()
        run = json.loads(
            (ROOT / f"validation/ssa-fixtures/{name}.run.json").read_bytes()
        )
        times = run["time"]["save_times"]
        laws = [analytic_probabilities(initial, birth, death, t) for t in times[1:]]
        species = json.loads(model)["species"][0]["name"]
        yield name, model, run["ensemble"]["root_seed"], times, species, laws, None
    for name, f in FIXTURES.items():
        model = {
            "schema": "crn-craft.model/0.1",
            "name": name,
            "convention": "stochastic_counts",
            "omega": 1.0,
            "parameters": f["parameters"],
            "species": [{"name": s, "initial_count": x} for s, x in f["species"]],
            "reactions": f["reactions"],
        }
        model_bytes = (json.dumps(model, indent=2) + "\n").encode()
        position = [s for s, _ in f["species"]].index(f["marginal"])
        states, edges, exits = generator(f)
        laws = [
            marginal(states, position, law)
            for law in uniformize(states, edges, exits)[1:]
        ]
        wrong = None
        if f.get("power"):
            ws, we, wx = generator(f, True)
            wrong = [marginal(ws, position, law) for law in uniformize(ws, we, wx)[1:]]
        yield name, model_bytes, f["seed"], TIMES, f["marginal"], laws, wrong


def cli(binary, folder, *args):
    done = subprocess.run(
        [str(binary), *args], cwd=folder, capture_output=True, check=False
    )
    if done.returncode != 0:
        raise RuntimeError(f"{args}: {done.stderr.decode()[:400]} {done.stdout[:400]}")
    return done.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--binary", type=Path, required=True, help="crn-craft CLI under test"
    )
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--replicates",
        type=int,
        help="smoke run with this many replicates; the report says so and is no gate run",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    global N, EPSILON
    if args.replicates is not None:
        N = args.replicates
        EPSILON = math.sqrt(math.log(2 * TESTS / ALPHA) / (2 * N))
    report = {
        "schema": "crn-craft.selector-law-validation/0.1",
        "smoke": args.replicates is not None,
        "replicates_per_fixture": N,
        "positive_time_cdf_tests": TESTS,
        "family_alpha": ALPHA,
        "dkw_epsilon": EPSILON,
        "checks": [],
        "power": [],
        "runs": [],
    }

    def save():
        report["passed"] = sum(r["passed"] for r in report["checks"])
        report["failed"] = sum(not r["passed"] for r in report["checks"])
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    try:
        for name, network, seed, times, species, laws, wrong in fixtures():
            with tempfile.TemporaryDirectory() as folder:
                folder = Path(folder)
                (folder / "network.json").write_bytes(network)
                converted = cli(binary, folder, "convert", "network.json")
                (folder / "model.json").write_bytes(converted)
                for selector in SELECTORS:
                    for wait in WAITS:
                        request = {
                            "schema": "crn-craft.run/1.0",
                            "model": {
                                "path": "model.json",
                                "sha256": hashlib.sha256(converted).hexdigest(),
                            },
                            "time": {"end": times[-1], "report_times": times},
                            "ensemble": {"replicates": N, "root_seed": seed},
                            "outputs": {
                                "statistics": ["mean", "variance"],
                                "keep_paths": N,
                                "keep_instances": 0,
                            },
                            "budget": {
                                "max_events_per_replicate": 1_000_000,
                                "max_sampling_steps_per_replicate": 100_000_000,
                                "max_instances": 16,
                                "threads": args.threads,
                                "max_output_bytes": 1 << 30,
                            },
                            "method": {
                                "kind": "static_direct",
                                "selector": selector,
                                "waiting_time": wait,
                            },
                        }
                        (folder / "run.json").write_bytes(
                            (json.dumps(request) + "\n").encode()
                        )
                        output = cli(binary, folder, "rules", "model.json", "run.json")
                        result = json.loads(output)
                        label = f"{name}/{selector}/{wait}"
                        if (
                            result["termination"] != "completed"
                            or result["completed_replicates"] != N
                            or len(result["paths"]) != N
                        ):
                            raise RuntimeError(f"{label}: incomplete")
                        column = result["observables"].index(species)
                        report["runs"].append(
                            {
                                "fixture": label,
                                "sampler": result["sampler"],
                                "events": result["events"],
                                "model_sha256": result["model_sha256"],
                                "run_sha256": result["run_sha256"],
                                "result_sha256": hashlib.sha256(output).hexdigest(),
                            }
                        )
                        power = []
                        for index, (t, law) in enumerate(
                            zip(times[1:], laws, strict=True), 1
                        ):
                            values = [
                                p["values"][index][column] for p in result["paths"]
                            ]
                            if any(v != int(v) for v in values):
                                raise RuntimeError(f"{label}: non-integer count")
                            values = [int(v) for v in values]
                            hist = Counter(values)
                            error = (
                                cdf_distance(hist, law, N)[0]
                                if isinstance(law, list)
                                else max_cdf_error(values, law)
                            )
                            mean = sum(values) / N
                            variance = (
                                sum(v * v for v in values) - sum(values) ** 2 / N
                            ) / (N - 1)
                            s = result["statistics"]
                            matches = abs(
                                s["mean"][index][column] - mean
                            ) <= 2e-12 * max(1, mean) and abs(
                                s["sample_variance"][index][column] - variance
                            ) <= 2e-12 * max(1, variance)
                            report["checks"].append(
                                {
                                    "fixture": label,
                                    "time": t,
                                    "max_cdf_error": error,
                                    "statistics_match": matches,
                                    "passed": error <= EPSILON and matches,
                                    "histogram": dict(sorted(hist.items())),
                                }
                            )
                            if wrong is not None:
                                power.append(max_cdf_error(values, wrong[index - 1]))
                        if wrong is not None:
                            rejected = any(v > EPSILON for v in power)
                            report["power"].append(
                                {
                                    "fixture": label,
                                    "wrong_law_errors": power,
                                    "rejected": rejected,
                                }
                            )
                            if not rejected:
                                report["checks"].append(
                                    {
                                        "fixture": label,
                                        "check": "wrong_law_power",
                                        "passed": False,
                                    }
                                )
                        save()
                        print(label, report["passed"], report["failed"], flush=True)
        if len([r for r in report["checks"] if "time" in r]) != TESTS:
            raise RuntimeError("test count changed")
    except Exception as error:
        report["checks"].append(
            {"check": "execution_error", "passed": False, "error": str(error)}
        )
        save()
        raise
    save()
    if report["failed"]:
        raise RuntimeError("declared statistical gate failed; rows retained")


if __name__ == "__main__":
    main()
