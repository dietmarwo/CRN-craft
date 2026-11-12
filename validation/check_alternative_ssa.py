"""Fixed 168-test gate for every exact Static Direct selector and wait variant.

The structure of the 0.6.0 gate on the unified engine and its Python API: each
fixture is a network model converted with `crn_craft.convert` and run with
`crn_craft.ensemble(...).simulate()` as Static Direct with every selector
(`linear`, `composition_rejection`, `sorting`, `rssa_cr`) and both waiting
times (`inverse`, `ziggurat`), n = 100,000 replicates with every path kept.
At every positive report time the largest CDF distance between the replicates'
marginal and the exact law (closed form for the immigration-death fixtures,
uniformization for the CME fixtures of check_cme.py) must be at most the DKW
bound of the declared family, and the statistics must equal the moments of the
kept paths. The fixtures with a wrong law (the falling-factorial convention)
must reject it at some time (power). Safe under Python -O.

Deviations from the 0.6.0 gate (docs/campaign/merge-0.7/M8a-log.md): the
variants `cached_direct` (removed algorithm) and `fast` (removed xoshiro256++
alternative) are gone, so the family has 21 positive-time tests in each of 8
variants, TESTS = 168 instead of 210, at the same alpha of 0.001 (the bound
is the DKW bound of the smaller family, epsilon = sqrt(ln(2 TESTS / alpha) /
(2 n)), slightly tighter). The seeds are those of 0.6.0; the trajectories
differ because the stream is the unified ChaCha8 stream (decision D4).
"""

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import crn_craft as crn
from check_ssa import analytic_probabilities, cdf_distance
from check_cme import (
    FIXTURES,
    TIMES,
    generator,
    marginal,
    max_cdf_error,
    network_model,
    uniformize,
)

ROOT = Path(__file__).resolve().parents[1]
N = 100_000
TESTS = 168
ALPHA = 0.001
EPSILON = math.sqrt(math.log(2 * TESTS / ALPHA) / (2 * N))
VARIANTS = list(
    itertools.product(
        ["linear", "composition_rejection", "sorting", "rssa_cr"],
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
            json.loads(m)["species"][0]["name"],
            laws,
            None,
        )
    for name, f in FIXTURES.items():
        model = network_model(name, f["parameters"], f["species"], f["reactions"])
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
        yield name, m, f["seed"], TIMES, f["marginal"], laws, wrong


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    report = dict(
        schema="crn-craft.alternative-ssa-validation/0.2",
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
        for name, network, seed, times, species, laws, wrong in fixtures():
            model = crn.convert(network)
            for algorithm, wait in VARIANTS:
                prepared = crn.ensemble(
                    model,
                    end=times[-1],
                    report_times=times,
                    replicates=N,
                    seed=seed,
                    threads=8,
                    keep_paths=N,
                    max_events=1_000_000,
                    max_sampling_steps=100_000_000,
                    max_output_bytes=268_435_456,
                    method={"kind": "static_direct", "waiting_time": wait},
                    selector=algorithm,
                    quantiles=(),
                )
                result = prepared.simulate()
                meta = result.metadata
                if meta["termination"] != "completed" or len(meta["paths"]) != N:
                    raise RuntimeError(f"{name}/{algorithm}/{wait}: incomplete")
                label = f"{name}/{algorithm}/{wait}"
                initial = next(
                    s for s in json.loads(network)["species"] if s["name"] == species
                )["initial_count"]
                series = result.paths.sel(variable=species).to_numpy()
                if (series[:, 0] != initial).any():
                    raise RuntimeError("initial state")
                report["runs"].append(
                    dict(
                        fixture=label,
                        sampler=meta["sampler"],
                        events=meta["events"],
                        rng_words=int(result.records["rng_words"].sum()),
                        model_sha256=meta["model_sha256"],
                        run_sha256=meta["run_sha256"],
                        result_sha256=hashlib.sha256(
                            json.dumps(
                                meta, sort_keys=True, separators=(",", ":")
                            ).encode()
                        ).hexdigest(),
                    )
                )
                power = []
                for index, (t, law) in enumerate(zip(times[1:], laws, strict=True), 1):
                    column = series[:, index]
                    if (column != column.astype(int)).any():
                        raise RuntimeError(f"{label}: non-integer count")
                    values = column.astype(int).tolist()
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
                    s = result.statistics
                    matches = bool(
                        abs(s["mean"][species].iloc[index] - mean)
                        <= 2e-12 * max(1, mean)
                        and abs(s["sample_variance"][species].iloc[index] - variance)
                        <= 2e-12 * max(1, variance)
                    )
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
