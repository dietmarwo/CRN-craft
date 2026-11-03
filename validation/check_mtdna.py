"""Frozen mtDNA A0–A6 law and event observation gates."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from check_rules import MTDNA, invoke, law_mtdna, require

ROOT = Path(__file__).resolve().parents[1]
N = 32_768
SEED = 20_261_103
ALPHA = 0.001
TIMES = [1.0, 3.0, 10.0]
MODEL_DIR = ROOT / "models/rules"
TRACKER = {
    "observable": "mutants",
    "first_passage": "N > 0 && mutants >= 0.5 * N",
    "keep_episodes": 0,
}


def read_model(name):
    return json.loads((MODEL_DIR / name).read_bytes())


def rules_run(
    binary,
    folder,
    model,
    *,
    method=None,
    times=TIMES,
    n=N,
    seed=SEED,
    threads=8,
    keep=True,
    trackers=(),
    log=None,
    max_events=1_000_000,
    max_sampling_steps=8_000_000,
    cancel=None,
):
    path = folder / "model.rules.json"
    path.write_text(json.dumps(model))
    request = {
        "schema": "crn-craft.rules-run/0.2"
        if method or trackers or log is not None
        else "crn-craft.rules-run/0.1",
        "model": {
            "path": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        },
        "time": {"end": times[-1], "report_times": times},
        "ensemble": {"replicates": n, "root_seed": seed},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": n if keep else 0,
            "keep_instances": 0,
        },
        "budget": {
            "max_events_per_replicate": max_events,
            "max_sampling_steps_per_replicate": max_sampling_steps,
            "max_instances": 100_000,
            "threads": threads,
            "max_output_bytes": 512 << 20,
        },
    }
    if method:
        request["method"] = {"kind": method} if isinstance(method, str) else method
    if trackers:
        request["trackers"] = list(trackers)
    if log is not None:
        request["event_log"] = log
    run = folder / "request.json"
    run.write_text(json.dumps(request))
    args = ["rules", path, run]
    if cancel:
        args.extend(["--cancel-file", cancel])
    return invoke(binary, *args)


def sample(binary, folder, model, **kwargs):
    code, result = rules_run(binary, folder, model, **kwargs)
    require(
        code == 0 and result["termination"] == "completed",
        str(result.get("diagnostics")),
    )
    paths = np.asarray([p["values"] for p in result["paths"]])
    return paths, result


def column(paths, result, name):
    return paths[:, :, result["observables"].index(name)]


def dkw(family):
    return math.sqrt(math.log(2 * family / ALPHA) / (2 * N))


def cdf_distance(x, y):
    points = np.sort(np.concatenate((x, y)))
    lhs = np.searchsorted(np.sort(x), points, side="right") / len(x)
    rhs = np.searchsorted(np.sort(y), points, side="right") / len(y)
    return float(np.max(np.abs(lhs - rhs)))


def row(rows, name, distance, bound, required=True):
    rows.append(
        {
            "test": name,
            "distance": float(distance),
            "bound": bound,
            "passed": bool(distance <= bound),
            "required": required,
        }
    )


def integer_cdf(xs, cdf):
    support = np.arange(len(cdf))
    require(
        np.all(xs >= 0) and np.all(xs == np.floor(xs)), "noninteger or negative count"
    )
    return np.max(
        np.abs(np.searchsorted(np.sort(xs), support, side="right") / len(xs) - cdf)
    )


def first_passage(result):
    values = (ob["trackers"][0]["first_passage"] for ob in result["observations"])
    return np.array([np.inf if x is None else x for x in values])


def check_mt(binary, folder, rows, summaries):
    oracle, _ = law_mtdna(MTDNA, TIMES)
    models = [
        read_model("mtdna-fixed.rules.json"),
        read_model("mtdna-three-class-cme.rules.json"),
    ]
    exact = None
    for i, (name, model, method) in enumerate(
        zip(
            ["clone", "tree", "static", "next"],
            [models[0], *[models[1]] * 3],
            [None, None, "static_direct", "next_reaction"],
            strict=True,
        )
    ):
        paths, result = sample(
            binary, folder, model, method=method, seed=SEED + i, trackers=[TRACKER]
        )
        columns = {key: column(paths, result, key) for key in oracle}
        for key, law in oracle.items():
            for j, t in enumerate(TIMES):
                row(
                    rows,
                    f"cme/{name}/{key}/{t}",
                    integer_cdf(columns[key][:, j], law[j]),
                    dkw(36),
                )
        if name == "tree":
            tree = columns
        if name in {"static", "next"}:
            for key in oracle:
                for j, t in enumerate(TIMES):
                    row(
                        rows,
                        f"two-sample/{name}/{key}/{t}",
                        cdf_distance(columns[key][:, j], tree[key][:, j]),
                        2 * dkw(36),
                    )
        if name == "static":
            exact = columns, first_passage(result)
        print("checked", name, flush=True)

    def outcomes(columns, passage):
        mut = columns["mutants"][:, -1]
        n = mut + columns["wild_type"][:, -1]
        fraction = np.divide(mut, n, out=np.zeros_like(mut), where=n > 0)
        return fraction, fraction >= 0.9, passage

    base, passage = exact
    ref = outcomes(base, passage)
    for i, epsilon in enumerate([0.03, 0.01]):
        paths, result = sample(
            binary,
            folder,
            models[1],
            method={"kind": "tau", "epsilon": epsilon},
            seed=SEED + 10 + i,
            trackers=[TRACKER],
        )
        columns = {key: column(paths, result, key) for key in oracle}
        for key in oracle:
            for j, t in enumerate(TIMES):
                require(np.min(columns[key][:, j]) >= 0, "negative tau sample")
                row(
                    rows,
                    f"tau/{epsilon}/{key}/{t}",
                    cdf_distance(columns[key][:, j], base[key][:, j]),
                    2 * dkw(48),
                    epsilon == 0.01,
                )
        approx = outcomes(columns, first_passage(result))
        for name, a, b in zip(
            ["fraction", "takeover", "passage"], approx, ref, strict=True
        ):
            row(
                rows,
                f"tau/{epsilon}/{name}",
                cdf_distance(a, b),
                2 * dkw(48),
                epsilon == 0.01,
            )
        summaries[str(epsilon)] = {
            "mean_fraction": float(approx[0].mean()),
            "reference_mean_fraction": float(ref[0].mean()),
            "takeover": float(approx[1].mean()),
            "reference_takeover": float(ref[1].mean()),
            "first_passage_fraction": float(np.isfinite(approx[2]).mean()),
            "reference_first_passage_fraction": float(np.isfinite(ref[2]).mean()),
        }
        print("checked tau", epsilon, flush=True)


def check_tracking(binary, folder, rows):
    model = read_model("mtdna-three-class-cme.rules.json")
    opts = dict(
        method="static_direct",
        n=32,
        trackers=[{**TRACKER, "keep_episodes": 128}],
        log={"max_entries": 1024},
        times=[0.0, 10.0],
    )
    _, one = sample(binary, folder, model, **opts, threads=1)
    _, many = sample(binary, folder, model, **opts, threads=20)
    one.pop("run_sha256")
    many.pop("run_sha256")
    require(one == many, "1/20 observed static runs differ")
    _, bare = sample(
        binary, folder, model, method="static_direct", n=32, times=[0.0, 10.0]
    )
    bare.pop("run_sha256")
    observed = one.pop("observations")
    require(one == bare, "tracking changed the numerical payload")
    idx = one["observables"].index("mutants")
    total_idx = one["observables"].index("N")
    for observation in observed:
        require(observation["omitted_log_entries"] == 0, "unexpected truncation")
        tracker = observation["trackers"][0]
        last_t, prev, integral = 0.0, 0.0, 0.0
        seen, gone = [], []
        episodes = []
        crossing = None
        start, peak, maximum = None, 0, 0
        counters = [0] * len(one["rules"])
        for event in observation["log"]:
            t = event["time"]
            value = event["observables"][idx]
            integral += prev * (t - last_t)
            if prev == 0 and value > 0:
                seen.append(t)
                start, peak = t, value
            if prev > 0 and value == 0:
                gone.append(t)
                episodes.append(
                    {"start": start, "end": t, "duration": t - start, "peak": peak}
                )
                start = None
            if start is not None:
                peak = max(peak, value)
            maximum = max(maximum, value)
            n = event["observables"][total_idx]
            if crossing is None and n > 0 and value >= 0.5 * n:
                crossing = t
            last_t, prev = t, value
            counters[event["rule"]] += 1
        integral += prev * (10.0 - last_t)
        if start is not None:
            episodes.append(
                {"start": start, "end": None, "duration": 10.0 - start, "peak": peak}
            )
        require(
            tracker["integral"] == integral and tracker["maximum"] == maximum,
            "integral or peak differs",
        )
        require(
            tracker["appearances"] == len(seen) and tracker["extinctions"] == len(gone),
            "crossing counts differ",
        )
        require(
            tracker["last_appearance"] == (seen[-1] if seen else None),
            "last appearance differs",
        )
        require(
            tracker["last_extinction"] == (gone[-1] if gone else None),
            "last extinction differs",
        )
        require(
            tracker["first_passage"] == crossing and tracker["episodes"] == episodes,
            "first passage or episodes differ",
        )
        require(
            counters == one["replicates"][observation["id"]]["rule_events"],
            "rule counters differ",
        )
    for method in [
        "tree",
        "next_reaction",
        {"kind": "tau", "epsilon": 0.01},
        {"kind": "binomial", "step": 1.0},
    ]:
        _, a = sample(binary, folder, model, method=method, n=32, threads=1)
        _, b = sample(binary, folder, model, method=method, n=32, threads=20)
        a.pop("run_sha256")
        b.pop("run_sha256")
        require(a == b, f"1/20 replay failed for {method}")
    flag = folder / "cancel"
    flag.touch()
    code, result = rules_run(
        binary, folder, model, n=1, method="static_direct", cancel=flag
    )
    require(
        code == 1
        and result["termination"] == "cancelled"
        and result["statistics"] is None,
        "cancel failed",
    )
    code, result = rules_run(
        binary, folder, model, n=1, method="static_direct", max_events=1
    )
    require(
        code == 1
        and result["termination"] == "event_limit"
        and result["statistics"] is None,
        "event budget failed",
    )
    _, capped = sample(
        binary,
        folder,
        model,
        n=1,
        method="static_direct",
        log={"max_entries": 2, "rules": ["degradation_wt"]},
        trackers=[TRACKER],
    )
    ob = capped["observations"][0]
    require(len(ob["log"]) == 2 and ob["omitted_log_entries"] > 0, "log cap failed")
    for reaction in model["rules"]:
        reaction["rate"] = "0"
    code, result = rules_run(
        binary,
        folder,
        model,
        n=1,
        method={"kind": "binomial", "step": 0.001},
        max_sampling_steps=2,
    )
    require(
        code == 1
        and result["termination"] == "sampling_limit"
        and result["replicates"][0]["rng_words"] == 0,
        "zero-event leap escaped its work budget",
    )
    row(rows, "events/exact-reconstruction-and-replay", 0.0, 0.0)
    print("checked tracking and replay", flush=True)


def poisson_cdf(mean, support=100):
    terms = [math.exp(-mean)]
    for k in range(1, support):
        terms.append(terms[-1] * mean / k)
    return np.cumsum(terms)


def check_atp(binary, folder, rows):
    model = read_model("immigration-death.rules.json")
    model["schema"] = "crn-craft.rules-model/0.2"
    model["types"].append({"name": "WT"})
    model["initial"].append({"type": "WT", "count": 1})
    model["observables"]["wt"] = {"type": "WT"}
    model["rules"] = model["rules"][:1]
    model["parameters"]["k_in"] = 4.0
    for initial in (0.25, 2.0):
        model["atp"] = {
            "wild_type": "wt",
            "total": "wt",
            "initial": initial,
            "production": 1.0,
            "maintenance": 0.0,
            "decay": 1.0,
            "scale": 1.0,
            "replication_rules": ["immigration"],
        }
        times = [0.5, 1.0, 3.0]
        paths, result = sample(binary, folder, model, times=times, method="thinning")
        for i, t in enumerate(times):
            rate_integral = 2.0 * math.log(
                (2.0 * math.exp(t) + initial - 1.0) / (1.0 + initial)
            )
            row(
                rows,
                f"atp/poisson/{initial}/{t}",
                integer_cdf(
                    column(paths, result, "X")[:, i], poisson_cdf(rate_integral)
                ),
                dkw(6),
            )
            require(
                np.max(
                    np.abs(
                        column(paths, result, "ATP")[:, i]
                        - (1 + (initial - 1) * math.exp(-t))
                    )
                )
                < 1e-12,
                "ATP flow grid mismatch",
            )
    model["parameters"]["k_in"] = 0.0
    model["atp"].update(initial=2.0, production=0.0, maintenance=1.0, decay=0.5)
    _, result = sample(
        binary, folder, model, method="thinning", n=1, times=[0.0, 0.5, 2.0]
    )
    require(
        result["replicates"][0]["termination"] == "stopped"
        and abs(result["replicates"][0]["stop_time"] - 2 * math.log(2)) < 1e-12,
        "ATP zero crossing",
    )
    row(rows, "atp/deterministic-stop", 0.0, 0.0)
    print("checked ATP", flush=True)


def check_leaps(binary, folder, rows):
    model = read_model("immigration-death.rules.json")
    model["parameters"]["k_in"] = 1000.0
    model["initial"] = [{"type": "X", "count": 1000}]
    times = [0.1, 0.5, 1.0]
    n = 8192
    method = {"kind": "tau", "epsilon": 0.01}
    paths, result = sample(binary, folder, model, method=method, n=n, times=times)
    for i, time in enumerate(times):
        p = math.exp(-time)
        binomial = np.array(
            [
                math.exp(
                    math.lgamma(1001)
                    - math.lgamma(k + 1)
                    - math.lgamma(1001 - k)
                    + k * math.log(p)
                    + (1000 - k) * math.log1p(-p)
                )
                for k in range(1001)
            ]
        )
        mean = 1000 * (1 - p)
        poisson = np.array(
            [
                math.exp(-mean + k * math.log(mean) - math.lgamma(k + 1))
                for k in range(1500)
            ]
        )
        cdf = np.cumsum(np.convolve(binomial, poisson))
        row(
            rows,
            f"tau/high-copy/{time}",
            integer_cdf(column(paths, result, "X")[:, i], cdf),
            math.sqrt(math.log(6 / ALPHA) / (2 * n)),
        )
    _, log = sample(
        binary,
        folder,
        model,
        method=method,
        n=1,
        times=times,
        log={"max_entries": 4096},
    )
    require(
        any(
            sum(k for _, k in event["firings"]) > 1
            for event in log["observations"][0]["log"]
        ),
        "adaptive run did not actually leap",
    )
    row(rows, "tau/multi-firing-leap", 0.0, 0.0)
    print("checked high-copy leaps", flush=True)


def network_sample(binary, folder, model, times):
    path = folder / "network.model.json"
    path.write_text(json.dumps(model))
    req = {
        "schema": "crn-craft.ssa-run/0.4",
        "model": {
            "path": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        },
        "method": {
            "kind": "ssa",
            "aggregator": "next_reaction",
            "waiting_time": "inverse",
        },
        "time": {"end": times[-1], "save_times": [0.0, *times]},
        "ensemble": {"replicates": N, "root_seed": SEED},
        "outputs": {"statistics": ["mean"], "keep_trajectories": N},
        "budget": {
            "max_events_per_replicate": 1_000_000,
            "max_sampling_steps_per_replicate": 8_000_000,
            "threads": 8,
            "max_output_bytes": 512 << 20,
        },
    }
    request = folder / "network.run.json"
    request.write_text(json.dumps(req))
    code, result = invoke(binary, "ssa", path, request)
    require(
        code == 0 and result["termination"] == "completed",
        str(result.get("diagnostics")),
    )
    return np.asarray([p["counts"] for p in result["trajectories"]])[:, 1:, 0]


def check_clocks(binary, folder, rows):
    for birth in [False, True]:
        rules = read_model("immigration-death.rules.json")
        k_in, k_local = (1.0, 0.5) if birth else (10.0, 1.0)
        rules["parameters"].update(k_in=k_in, k_out=k_local)
        products = [{"type": "X", "copy": "x"}] * 2 if birth else []
        rules["rules"][1]["products"] = products
        net = {
            "schema": "crn-craft.model/0.1",
            "name": "time_reuse",
            "convention": "stochastic_counts",
            "omega": 1.0,
            "parameters": {},
            "species": [{"name": "X", "initial_count": 0}],
            "reactions": [
                {
                    "name": "immigration",
                    "reactants": {},
                    "products": {"X": 1},
                    "rate": {"kind": "mass_action", "coefficient": k_in},
                },
                {
                    "name": "local",
                    "reactants": {"X": 1},
                    "products": {"X": 2} if birth else {},
                    "rate": {"kind": "mass_action", "coefficient": k_local},
                },
            ],
        }
        times = [0.5, 1.0, 3.0]
        paths, result = sample(
            binary, folder, rules, method="next_reaction", times=times
        )
        native = network_sample(binary, folder, net, times)
        for i, t in enumerate(times):
            if birth:
                p = math.exp(-k_local * t)
                cdf = np.cumsum([(k + 1) * p * p * (1 - p) ** k for k in range(200)])
            else:
                cdf = poisson_cdf(k_in / k_local * (1 - math.exp(-k_local * t)))
            for label, data in [
                ("rules", column(paths, result, "X")),
                ("network", native),
            ]:
                row(
                    rows,
                    f"clocks/{label}/{birth}/{t}",
                    integer_cdf(data[:, i], cdf),
                    dkw(12),
                )
    print("checked next reaction laws", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {"paths": N, "seed": SEED, "alpha": ALPHA, "checks": [], "tau": {}}
    try:
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            check_mt(args.binary, folder, report["checks"], report["tau"])
            check_tracking(args.binary, folder, report["checks"])
            check_atp(args.binary, folder, report["checks"])
            check_leaps(args.binary, folder, report["checks"])
            check_clocks(args.binary, folder, report["checks"])
        report["passed"] = all(
            r["passed"] or not r["required"] for r in report["checks"]
        )
    finally:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    require(report["passed"], "mtDNA gate failure")


if __name__ == "__main__":
    main()
