"""Independent analytic campaign, trajectory oracle and archive/CLI checks.

Targets are frozen in docs/decisions/iteration-04-validation.md. Timings are
written separately; the normal/-O acceptance report is fixed-work evidence.

Since 0.7.0 (M8a) the campaign runs on the unified engine: every network
fixture is converted with `crn-craft convert` and run as a `crn-craft.run/1.0`
request, Static Direct with the `linear` selector and `inverse` waits (the
checked reference variant, `chacha8-u53/rules-static-direct/1`), by `crn-craft
run MODEL REQUEST --archive DIR`, and the `crn-craft.archive/1.0` directory is
verified independently. The ChaCha8 stream is the one of `crn_craft::rules::
stream` (reference_rng.py). The sample sizes, seeds of the fixtures, the
family-wise error rate and the DKW bound are those of 0.6.0; the trajectories
differ from the 0.6.0 ones because the stream's key domain changed (decision
D4), which the statistical tests do not depend on.
"""

import argparse
from bisect import bisect_right
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
import time

from jsonschema import Draft7Validator
from check_cme import convert_network, network_model, run_request
from reference_rng import vectors, words

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 0.001
N = 100000
M = 9
EPSILON = math.sqrt(math.log(2 * M / ALPHA) / (2 * N))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate key {key}")
        result[key] = value
    return result


def load(path):
    return json.loads(path.read_bytes(), object_pairs_hook=no_duplicates)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def independent_archive(path, model, run):
    """Check a `crn-craft.archive/1.0` directory without the Rust verifier."""
    names = {"model.json", "request.json", "result.json", "manifest.json", "COMPLETE"}
    require({p.name for p in path.iterdir()} == names, "archive names differ")
    require(
        all(p.is_file() and not p.is_symlink() for p in path.iterdir()),
        "nonregular payload",
    )
    require((path / "model.json").read_bytes() == model, "original model bytes changed")
    require(
        (path / "request.json").read_bytes() == run, "original request bytes changed"
    )
    manifest_bytes = (path / "manifest.json").read_bytes()
    require(
        (path / "COMPLETE").read_text() == digest(manifest_bytes) + "\n",
        "completion hash",
    )
    manifest = load(path / "manifest.json")
    require(manifest["schema"] == "crn-craft.archive/1.0", "archive version")
    require(len(manifest["files"]) == 3, "manifest length")
    require(
        {f["path"] for f in manifest["files"]} == names - {"manifest.json", "COMPLETE"},
        "manifest duplicate/path",
    )
    for entry in manifest["files"]:
        data = (path / entry["path"]).read_bytes()
        require(
            entry["bytes"] == len(data) and entry["sha256"] == digest(data),
            "payload integrity",
        )
    # The output cap bounds the printed result (the model and request are the
    # caller's own bytes).
    require(
        (path / "result.json").stat().st_size
        <= json.loads(run)["budget"]["max_output_bytes"],
        "result exceeds output cap",
    )
    result = load(path / "result.json")
    require(result["schema"] == "crn-craft.result/1.0", "result contract")
    require(
        result["model_sha256"] == digest(model) and result["run_sha256"] == digest(run),
        "result identity",
    )
    require(
        [r["id"] for r in result["replicates"]]
        == list(range(result["requested_replicates"])),
        "replicate accounting",
    )
    require(
        sum(r["events"] for r in result["replicates"]) == result["events"],
        "event accounting",
    )
    require(
        sum(r["termination"] in ("completed", "stopped") for r in result["replicates"])
        == result["completed_replicates"],
        "completion accounting",
    )
    if result["termination"] != "completed":
        require(result["statistics"] is None, "conditional statistics mislabelled")
    return result


def analytic_probabilities(initial, birth, death, horizon):
    survival = math.exp(-death * horizon)
    binomial = [
        math.comb(initial, k) * survival**k * (1 - survival) ** (initial - k)
        for k in range(initial + 1)
    ]
    mu = birth / death * (1 - survival)
    poisson = [math.exp(-mu)]
    for k in range(1, 61):
        poisson.append(poisson[-1] * mu / k)
    # For these frozen means <= 3, omitted mass is far below 1e-40.
    require(
        mu <= 3 and poisson[-1] * mu / 61 / (1 - mu / 62) < 1e-40, "analytic tail bound"
    )
    probabilities = [0.0] * (len(binomial) + len(poisson) - 1)
    for i, first in enumerate(binomial):
        for j, second in enumerate(poisson):
            probabilities[i + j] += first * second
    return probabilities


def cdf_distance(histogram, probabilities, samples):
    empirical = 0
    exact = 0.0
    distance = 0.0
    location = 0
    for k in range(max(max(histogram), len(probabilities) - 1) + 1):
        empirical += histogram[k]
        if k < len(probabilities):
            exact += probabilities[k]
        difference = abs(empirical / samples - exact)
        if difference > distance:
            distance, location = difference, k
    return distance, location


def reference_paths(model, run):
    names = [s["name"] for s in model["species"]]
    output = []
    for replicate in range(run["ensemble"]["replicates"]):
        state = [s["initial_count"] for s in model["species"]]
        generator = words(run["ensemble"]["root_seed"], replicate)
        event_times = [0.0]
        states = [state.copy()]
        clock = 0.0
        while True:
            rates = []
            for reaction in model["reactions"]:
                rate = reaction["rate"]["coefficient"]
                for species, order in reaction["reactants"].items():
                    count = state[names.index(species)]
                    rate *= math.comb(count, order) if count >= order else 0
                rates.append(rate)
            total = sum(rates)
            if not total:
                break
            bits = next(generator) >> 11
            while not bits:
                bits = next(generator) >> 11
            # The engine's inverse wait is -log(1 - u) / total, u in (0, 1).
            clock += -math.log(1.0 - bits / 2**53) / total
            if clock > run["time"]["end"]:
                break
            threshold = (next(generator) >> 11) / 2**53 * total
            cumulative = 0.0
            for index, rate in enumerate(rates):
                cumulative += rate
                if threshold < cumulative:
                    selected = index
                    break
            else:
                raise RuntimeError("reference selection failed")
            reaction = model["reactions"][selected]
            for i, name in enumerate(names):
                state[i] += reaction["products"].get(name, 0) - reaction[
                    "reactants"
                ].get(name, 0)
            event_times.append(clock)
            states.append(state.copy())
        output.append(
            {
                "id": replicate,
                "events": len(states) - 1,
                "counts": [
                    states[bisect_right(event_times, t) - 1]
                    for t in run["time"]["save_times"]
                ],
            }
        )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--latency-output", type=Path)
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    rows = []
    latencies = {}

    def invoke(*args, expected=0, timeout=180):
        p = subprocess.run(
            [binary, *map(str, args)], capture_output=True, timeout=timeout, check=False
        )
        require(
            p.returncode == expected,
            f"{args[0]} exit {p.returncode}: {p.stdout[-2000:]!r} {p.stderr[-1000:]!r}",
        )
        require(not p.stderr, "unexpected CLI stderr")
        return json.loads(p.stdout)

    def execute(model_path, run_path, out, expected=0, extra=()):
        summary = invoke(
            "run", model_path, run_path, "--archive", out, *extra, expected=expected
        )
        require(
            summary["valid"] and summary["archive_complete"],
            "simulation archive not finalized",
        )
        verified = invoke("verify", out)
        require(verified["valid"], "Rust archive verification failed")
        return independent_archive(out, model_path.read_bytes(), run_path.read_bytes())

    def write_request(folder, model_bytes, name="run.json", **options):
        spec = run_request(model_bytes, **options)
        path = folder / name
        path.write_text(json.dumps(spec) + "\n")
        return spec, path

    try:
        require(
            vectors() == load(ROOT / "validation/rng-vectors-rules.json"),
            "independent RNG vectors changed",
        )
        rows.append({"check": "independent_rng_360_words", "passed": True})
        with tempfile.TemporaryDirectory(prefix="crn-ssa-acceptance-") as temporary:
            folder = Path(temporary)
            for name, initial, birth, death, seed in [
                ("pure-decay", 20, 0.0, 0.7, 2026110401),
                ("immigration-zero", 0, 3.0, 0.8, 2026110402),
                ("immigration-nonzero", 12, 3.0, 0.8, 2026110403),
            ]:
                base = ROOT / "validation/ssa-fixtures"
                network_path = base / (name + ".model.json")
                frozen = load(base / (name + ".run.json"))
                require(
                    frozen["ensemble"] == {"replicates": N, "root_seed": seed},
                    "frozen fixture changed",
                )
                work = folder / name
                work.mkdir()
                model_bytes = convert_network(binary, network_path.read_bytes(), work)
                model_path = work / "model.json"
                model_path.write_bytes(model_bytes)
                _, run_path = write_request(
                    work,
                    model_bytes,
                    times=frozen["time"]["save_times"],
                    replicates=N,
                    seed=seed,
                    threads=8,
                    keep=N,
                    max_events=frozen["budget"]["max_events_per_replicate"],
                    max_output_bytes=frozen["budget"]["max_output_bytes"],
                )
                out = work / "archive"
                result = execute(model_path, run_path, out)
                require(
                    result["completed_replicates"] == N
                    and result["termination"] == "completed",
                    "analytic campaign incomplete",
                )
                require(
                    len(result["paths"]) == N,
                    "all analytic paths must be retained",
                )
                column = result["observables"].index("X")
                require(
                    all(
                        p["values"][0][column] == initial and len(p["values"]) == 4
                        for p in result["paths"]
                    ),
                    "initial law/save grid",
                )
                stats = result["statistics"]
                for index, horizon in enumerate([0.25, 1.0, 2.0], start=1):
                    values = [p["values"][index][column] for p in result["paths"]]
                    require(
                        all(x == int(x) and x >= 0 for x in values),
                        "invalid counts",
                    )
                    values = [int(x) for x in values]
                    histogram = Counter(values)
                    probabilities = analytic_probabilities(
                        initial, birth, death, horizon
                    )
                    distance, location = cdf_distance(histogram, probabilities, N)
                    mean = sum(values) / N
                    variance = (sum(x * x for x in values) - sum(values) ** 2 / N) / (
                        N - 1
                    )
                    numerical = abs(stats["mean"][index][column] - mean) <= 2e-12 * max(
                        1, mean
                    )
                    numerical &= abs(
                        stats["sample_variance"][index][column] - variance
                    ) <= 2e-12 * max(1, variance)
                    numerical &= (
                        abs(
                            stats["mean_standard_error"][index][column]
                            - math.sqrt(variance / N)
                        )
                        <= 2e-12
                    )
                    rows.append(
                        {
                            "check": name,
                            "time": horizon,
                            "passed": distance <= EPSILON and numerical,
                            "cdf_distance": distance,
                            "cdf_argmax": location,
                            "epsilon": EPSILON,
                            "histogram": dict(sorted(histogram.items())),
                            "mean": mean,
                            "sample_variance": variance,
                            "statistics_match_independent_reduction": numerical,
                            "events": result["events"],
                            "model_sha256": result["model_sha256"],
                            "run_sha256": result["run_sha256"],
                            "result_sha256": digest((out / "result.json").read_bytes()),
                        }
                    )
                require(
                    cdf_distance(
                        Counter({0: N}),
                        analytic_probabilities(initial, birth, death, 1.0),
                        N,
                    )[0]
                    > EPSILON,
                    "CDF checker accepted a zeroed solver",
                )
            rows.append(
                {
                    "check": "cdf_checker_rejects_zeroed_solver_in_all_three_models",
                    "passed": True,
                }
            )

            network = network_model(
                "reference",
                {},
                [("A", 5), ("B", 0)],
                [
                    {
                        "name": name,
                        "reactants": left,
                        "products": right,
                        "rate": {"kind": "mass_action", "coefficient": coefficient},
                    }
                    for name, left, right, coefficient in [
                        ("birth", {}, {"A": 1}, 1.5),
                        ("death", {"A": 1}, {}, 0.125),
                        ("dimer", {"A": 2}, {"B": 1}, 0.5),
                        ("split", {"B": 1}, {"A": 2}, 1.0),
                    ]
                ],
                omega=137.0,
            )
            work = folder / "reference"
            work.mkdir()
            model_bytes = convert_network(
                binary, (json.dumps(network) + "\n").encode(), work
            )
            model_path = work / "model.json"
            model_path.write_bytes(model_bytes)
            names = [s["name"] for s in network["species"]]
            reference_times = [0.0, 0.25, 1.0, 2.0]
            options = {
                "times": reference_times,
                "replicates": 64,
                "seed": 123,
                "keep": 64,
                "max_events": 10000,
                "max_output_bytes": 4_000_000,
            }
            oracle_spec = {
                "ensemble": {"replicates": 64, "root_seed": 123},
                "time": {"end": 2.0, "save_times": reference_times},
            }
            oracle = reference_paths(network, oracle_spec)
            baseline = None
            validator = Draft7Validator(invoke("describe", "result-1.0"))
            order = [None] * len(names)
            for threads in [1, 2, 8]:
                _, run_path = write_request(
                    work, model_bytes, threads=threads, **options
                )
                result = execute(model_path, run_path, work / f"threads-{threads}")
                validator.validate(result)
                order = [result["observables"].index(n) for n in names]
                require(
                    [
                        {
                            "id": p["id"],
                            "counts": [[row[i] for i in order] for row in p["values"]],
                        }
                        for p in result["paths"]
                    ]
                    == [
                        {
                            "id": p["id"],
                            "counts": [[float(c) for c in row] for row in p["counts"]],
                        }
                        for p in oracle
                    ],
                    "independent trajectory oracle disagrees",
                )
                require(
                    [p["events"] for p in result["replicates"]]
                    == [p["events"] for p in oracle],
                    "independent event counts disagree",
                )
                del result["run_sha256"]
                if baseline is None:
                    baseline = result
                require(result == baseline, "numerical payload changes with threads")
                rows.append(
                    {
                        "check": f"reference_64_paths_threads_{threads}",
                        "passed": True,
                        "events": result["events"],
                    }
                )

            work = folder / "enzyme"
            work.mkdir()
            enzyme_network = (
                ROOT / "models/examples/enzyme-mechanism.model.json"
            ).read_bytes()
            enzyme_frozen = load(ROOT / "models/examples/enzyme-mechanism.run.json")
            enzyme_bytes = convert_network(binary, enzyme_network, work)
            enzyme_model = work / "model.json"
            enzyme_model.write_bytes(enzyme_bytes)
            enzyme_names = [s["name"] for s in json.loads(enzyme_network)["species"]]
            _, enzyme_run = write_request(
                work,
                enzyme_bytes,
                times=enzyme_frozen["time"]["save_times"],
                replicates=enzyme_frozen["ensemble"]["replicates"],
                seed=enzyme_frozen["ensemble"]["root_seed"],
                threads=8,
                keep=enzyme_frozen["ensemble"]["replicates"],
                max_events=enzyme_frozen["budget"]["max_events_per_replicate"],
                max_output_bytes=enzyme_frozen["budget"]["max_output_bytes"],
            )
            enzyme = execute(enzyme_model, enzyme_run, work / "archive")
            order = [enzyme["observables"].index(n) for n in enzyme_names]
            for path in enzyme["paths"]:
                for row in path["values"]:
                    substrate, catalyst, complex_count, product = (
                        row[i] for i in order
                    )
                    require(
                        catalyst + complex_count == 10
                        and substrate + complex_count + product == 50,
                        "enzyme conservation failed in SSA",
                    )
            rows.append(
                {
                    "check": "enzyme_conservation_paths",
                    "passed": True,
                    "replicates": enzyme["completed_replicates"],
                    "events": enzyme["events"],
                }
            )

            # Actual process cancellation, while events are being generated.
            network["reactions"] = [
                {
                    "name": "fast",
                    "reactants": {},
                    "products": {"A": 1},
                    "rate": {"kind": "mass_action", "coefficient": 1e9},
                }
            ]
            work = folder / "cancel"
            work.mkdir()
            model_bytes = convert_network(
                binary, (json.dumps(network) + "\n").encode(), work
            )
            model_path = work / "model.json"
            model_path.write_bytes(model_bytes)
            options = {
                "times": [0.0, 1e6],
                "replicates": 8,
                "seed": 123,
                "keep": 8,
                "max_events": 10_000_000_000,
                "max_output_bytes": 4_000_000,
            }
            _, run_path = write_request(work, model_bytes, threads=8, **options)
            cancel = work / "cancel"
            out = work / "cancelled"
            process = subprocess.Popen(
                [
                    binary,
                    "run",
                    str(model_path),
                    str(run_path),
                    "--archive",
                    str(out),
                    "--cancel-file",
                    str(cancel),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            ready_deadline = time.monotonic() + 10
            while (
                not (out / "model.json").exists()
                and process.poll() is None
                and time.monotonic() < ready_deadline
            ):
                time.sleep(0.005)
            time.sleep(0.03)
            start = time.monotonic()
            cancel.write_text("cancel\n")
            try:
                stdout, stderr = process.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.communicate()
                raise RuntimeError("cooperative cancellation timed out") from None
            latency = time.monotonic() - start
            latencies["active_cancel_seconds"] = latency
            require(
                process.returncode == 1 and not stderr,
                f"cancel CLI failure: {stdout!r} {stderr!r}",
            )
            require(invoke("verify", out)["valid"], "cancelled archive rejected")
            result = independent_archive(
                out, model_path.read_bytes(), run_path.read_bytes()
            )
            latencies["cancelled_events"] = result["events"]
            require(
                result["termination"] == "cancelled" and result["events"] > 0,
                "cancellation did not interrupt active work",
            )
            rows.append(
                {
                    "check": "active_cli_cancellation_under_one_second",
                    "passed": latency < 1.0,
                }
            )
            _, run_path = write_request(
                work, model_bytes, threads=8, wall_time_s=0.02, **options
            )
            result = execute(model_path, run_path, work / "deadline", expected=1)
            require(
                result["termination"] == "wall_time_limit", "deadline was not enforced"
            )
            rows.append({"check": "live_wall_time_stop", "passed": True})
    except Exception as error:
        rows.append({"check": "execution_error", "passed": False, "error": str(error)})
    report = {
        "schema": "crn-craft.ssa-validation/0.2",
        "replicates_per_analytic_fixture": N,
        "analytic_fixtures": 3,
        "positive_time_cdf_tests": M,
        "family_alpha": ALPHA,
        "dkw_epsilon": EPSILON,
        "passed": sum(row["passed"] for row in rows),
        "failed": sum(not row["passed"] for row in rows),
        "checks": rows,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    if args.latency_output:
        args.latency_output.write_text(json.dumps(latencies, indent=2) + "\n")
    print(text, end="")
    if report["failed"]:
        raise RuntimeError("SSA acceptance failed; failed rows retained")


if __name__ == "__main__":
    main()
