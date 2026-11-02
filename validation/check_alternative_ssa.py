"""Fixed 210-test gate for all exact SSA/wait variants; safe under Python -O."""

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import crn_craft as crn
from check_ssa import analytic_probabilities, cdf_distance
from check_cme import FIXTURES, TIMES, generator, uniformize, marginal, max_cdf_error

ROOT = Path(__file__).resolve().parents[1]
N = 100_000
TESTS = 210
ALPHA = 0.001
EPSILON = math.sqrt(math.log(2 * TESTS / ALPHA) / (2 * N))
VARIANTS = list(
    itertools.product(
        ["direct", "cached_direct", "direct_cr", "sorting_direct", "rssa_cr"],
        ["inverse", "ziggurat"],
    )
)


def fixtures():
    for name, initial, birth, death in [
        ("pure-decay", 20, 0.0, 0.7),
        ("immigration-zero", 0, 3.0, 0.8),
        ("immigration-nonzero", 12, 3.0, 0.8),
    ]:
        m = (ROOT / f"validation/ssa-fixtures/{name}.model.json").read_bytes()
        run = json.loads(
            (ROOT / f"validation/ssa-fixtures/{name}.run.json").read_bytes()
        )
        laws = [
            analytic_probabilities(initial, birth, death, t)
            for t in run["time"]["save_times"][1:]
        ]
        yield (
            name,
            m,
            run["ensemble"]["root_seed"],
            run["time"]["save_times"],
            0,
            laws,
            None,
        )
    for name, f in FIXTURES.items():
        model = dict(
            schema="crn-craft.model/0.1",
            name=name,
            convention="stochastic_counts",
            omega=1.0,
            parameters=f["parameters"],
            species=[dict(name=s, initial_count=x) for s, x in f["species"]],
            reactions=f["reactions"],
        )
        m = (json.dumps(model, indent=2) + "\n").encode()
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
        yield name, m, f["seed"], TIMES, position, laws, wrong


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    report = dict(
        schema="crn-craft.alternative-ssa-validation/0.1",
        replicates_per_fixture=N,
        positive_time_cdf_tests=TESTS,
        family_alpha=ALPHA,
        dkw_epsilon=EPSILON,
        checks=[],
        power=[],
        runs=[],
    )

    def save():
        report["passed"] = sum(r["passed"] for r in report["checks"])
        report["failed"] = sum(not r["passed"] for r in report["checks"])
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    try:
        for name, model, seed, times, position, laws, wrong in fixtures():
            for algorithm, wait in VARIANTS:
                request = crn.ssa_spec(
                    model,
                    aggregator=algorithm,
                    waiting_time=wait,
                    end=times[-1],
                    save_times=times,
                    replicates=N,
                    seed=seed,
                    threads=8,
                    keep=N,
                    max_events=1_000_000,
                    max_sampling_steps=100_000_000,
                    max_output_bytes=268_435_456,
                )
                result = crn.simulate_ssa(model, request)
                if (
                    result["termination"] != "completed"
                    or len(result["trajectories"]) != N
                ):
                    raise RuntimeError(f"{name}/{algorithm}/{wait}: incomplete")
                label = f"{name}/{algorithm}/{wait}"
                initial = json.loads(model)["species"][position]["initial_count"]
                if any(
                    t["counts"][0][position] != initial for t in result["trajectories"]
                ):
                    raise RuntimeError("initial state")
                report["runs"].append(
                    dict(
                        fixture=label,
                        events=result["events"],
                        work=result["work"],
                        model_sha256=result["model_sha256"],
                        run_sha256=result["run_sha256"],
                        result_sha256=hashlib.sha256(
                            json.dumps(
                                result, sort_keys=True, separators=(",", ":")
                            ).encode()
                        ).hexdigest(),
                    )
                )
                power = []
                for index, (t, law) in enumerate(zip(times[1:], laws, strict=True), 1):
                    values = [
                        v["counts"][index][position] for v in result["trajectories"]
                    ]
                    hist = Counter(values)
                    error = (
                        cdf_distance(hist, law, N)[0]
                        if isinstance(law, list)
                        else max_cdf_error(values, law)
                    )
                    mean = sum(values) / N
                    variance = (sum(v * v for v in values) - sum(values) ** 2 / N) / (
                        N - 1
                    )
                    s = result["statistics"]
                    matches = abs(s["mean"][index][position] - mean) <= 2e-12 * max(
                        1, mean
                    ) and abs(
                        s["sample_variance"][index][position] - variance
                    ) <= 2e-12 * max(1, variance)
                    report["checks"].append(
                        dict(
                            fixture=label,
                            time=t,
                            max_cdf_error=error,
                            statistics_match=matches,
                            passed=error <= EPSILON and matches,
                            histogram=dict(sorted(hist.items())),
                        )
                    )
                    if wrong is not None:
                        power.append(max_cdf_error(values, wrong[index - 1]))
                if wrong is not None:
                    rejected = any(v > EPSILON for v in power)
                    report["power"].append(
                        dict(fixture=label, wrong_law_errors=power, rejected=rejected)
                    )
                    if not rejected:
                        report["checks"].append(
                            dict(fixture=label, check="wrong_law_power", passed=False)
                        )
                save()
                print(label, report["passed"], report["failed"], flush=True)
        if len([r for r in report["checks"] if "time" in r]) != TESTS:
            raise RuntimeError("test count changed")
    except Exception as error:
        report["checks"].append(
            dict(check="execution_error", passed=False, error=str(error))
        )
        save()
        raise
    save()
    if report["failed"]:
        raise RuntimeError("declared statistical gate failed; rows retained")


if __name__ == "__main__":
    main()
