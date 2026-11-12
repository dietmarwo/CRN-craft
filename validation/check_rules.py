"""K1 law gates for the rules engine (dynamic, attributed species).

Targets are frozen in docs/decisions/rules-k0-k1-2026-11-02.md: five fixtures,
100,000 replicates each, one family of DKW one-sample bounds at family alpha
0.001, plus a wrong-law power check per fixture. Exact laws come from the
master equation (uniformization in numpy) or from the Poisson law; they share
no code with the Rust engine. The report holds fixed-work evidence only;
timings are printed, not written.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from jsonschema import Draft7Validator

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 0.001
N = 100_000
SEED = 20261102
# Marginal comparisons: 3 + (3 + 3) + 9 + 6 + 6.
K = 30
EPSILON = math.sqrt(math.log(2 * K / ALPHA) / (2 * N))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate key {key}")
        result[key] = value
    return result


def load(text):
    return json.loads(text, object_pairs_hook=no_duplicates)


# ---------------------------------------------------------------------------
# Exact laws


def uniformize(generator, p0, times):
    """Transient distributions p0 exp(Q t) at increasing times."""
    q = np.asarray(generator, dtype=float)
    rate = float(np.max(-np.diag(q)))
    if rate == 0.0:
        return [p0.copy() for _ in times]
    step = np.eye(len(p0)) + q / rate
    out = []
    for t in times:
        lam = rate * t
        # Poisson weights in log space; the truncated tail is below 1e-25.
        kmax = int(lam + 12 * math.sqrt(lam) + 50)
        v = p0.copy()
        acc = np.zeros_like(p0)
        logw = -lam
        total = 0.0
        for k in range(kmax + 1):
            if k:
                v = v @ step
                logw += math.log(lam) - math.log(k)
            w = math.exp(logw)
            acc += w * v
            total += w
        require(abs(1.0 - total) < 1e-9, f"uniformization weights sum to {total}")
        out.append(acc)
    return out


def chain(states, transitions):
    """Dense generator from a list of states and a function state -> [(state, rate)]."""
    index = {s: i for i, s in enumerate(states)}
    q = np.zeros((len(states), len(states)))
    for s in states:
        i = index[s]
        for target, rate in transitions(s):
            if rate > 0:
                q[i, index[target]] += rate
                q[i, i] -= rate
    return q, index


def marginal_cdfs(states, distributions, observables, support):
    """CDFs of integer-valued observables on 0..support."""
    out = {}
    for name, f in observables.items():
        values = np.array([f(s) for s in states])
        rows = []
        for p in distributions:
            pmf = np.bincount(values, weights=p, minlength=support + 1)[: support + 1]
            rows.append(np.cumsum(pmf))
        out[name] = rows
    return out


def law_immigration_death(k_in, k_out, times, support=60):
    out = []
    for t in times:
        mean = k_in / k_out * (1 - math.exp(-k_out * t))
        pmf = [
            math.exp(-mean + k * math.log(mean) - math.lgamma(k + 1))
            for k in range(support + 1)
        ]
        out.append(np.cumsum(pmf))
    return {"X": out}, 0.0


def law_dimerization(k_in, k2, times, nmax=120):
    # States 0..nmax and an absorbing leak state nmax + 1.
    states = list(range(nmax + 2))

    def moves(n):
        if n > nmax:
            return []
        return [(n + 1, k_in), (n - 2, k2 * n * (n - 1) / 2) if n >= 2 else (n, 0.0)]

    q, _ = chain(states, moves)
    p0 = np.zeros(len(states))
    p0[0] = 1
    dists = uniformize(q, p0, times)
    leak = max(float(p[-1]) for p in dists)
    cdfs = marginal_cdfs(states[:-1], [p[:-1] for p in dists], {"X": lambda n: n}, nmax)
    return cdfs, leak


def law_mtdna(p, times, kowald=True):
    """Fixed-population mtDNA chain on (wild types, sw mutants, sm mutants).

    Replication fills vacancies V = max(0, K - N) with weights s / D, where
    Kowald's competition term is D = w + sm (a + b); the alternative D = w +
    sw a + sm b weights every molecule by its own replication rate.
    """
    big_k = int(p["K"])
    states = [
        (w, a, b)
        for w in range(big_k + 1)
        for a in range(big_k + 1 - w)
        for b in range(big_k + 1 - w - a)
    ]

    def moves(s):
        w, a, b = s
        n = w + a + b
        out = [
            ((w - 1, a, b), p["delta"] * w),
            ((w, a - 1, b), p["delta"] * a),
            ((w, a, b - 1), p["delta"] * b),
        ]
        vacancy = max(0, big_k - n)
        if n and vacancy:
            d = w + p["sm"] * (a + b) if kowald else w + p["sw"] * a + p["sm"] * b
            g = vacancy / d
            out += [
                ((w + 1, a, b), g * p["sw"] * (1 - p["mu"]) * w),
                ((w, a, b + 1), g * p["sw"] * p["mu"] * w * p["fracAdv"]),
                ((w, a + 1, b), g * p["sw"] * p["mu"] * w * (1 - p["fracAdv"])),
                ((w, a + 1, b), g * p["sw"] * a),
                ((w, a, b + 1), g * p["sm"] * b),
            ]
        return out

    q, index = chain(states, moves)
    p0 = np.zeros(len(states))
    p0[index[(int(p["w0"]), 0, 0)]] = 1
    dists = uniformize(q, p0, times)
    observables = {
        "wild_type": lambda s: s[0],
        "mutants": lambda s: s[1] + s[2],
        "advantaged_mutants": lambda s: s[2],
    }
    return marginal_cdfs(states, dists, observables, big_k), 0.0


def partitions(n, largest=None):
    largest = n if largest is None else largest
    if n == 0:
        yield ()
        return
    for k in range(min(n, largest), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def law_aggregation(kernel, k, monomers, times):
    states = list(partitions(monomers))
    require(len(states) == {12: 77}.get(monomers, len(states)), "partition count")

    def moves(s):
        out = []
        counts = {}
        for size in s:
            counts[size] = counts.get(size, 0) + 1
        sizes = sorted(counts)
        for i, x in enumerate(sizes):
            for y in sizes[i:]:
                pairs = (
                    counts[x] * (counts[x] - 1) / 2 if x == y else counts[x] * counts[y]
                )
                if not pairs:
                    continue
                merged = list(s)
                merged.remove(x)
                merged.remove(y)
                merged.append(x + y)
                out.append(
                    (tuple(sorted(merged, reverse=True)), k * kernel(x, y) * pairs)
                )
        return out

    q, index = chain(states, moves)
    p0 = np.zeros(len(states))
    p0[index[(1,) * monomers]] = 1
    dists = uniformize(q, p0, times)
    observables = {"clusters": len, "largest": max}
    return marginal_cdfs(states, dists, observables, monomers), 0.0


def constant(_x, _y):
    return 1.0


def product(x, y):
    return float(x * y)


# ---------------------------------------------------------------------------
# Fixtures (frozen)

MTDNA = {
    "K": 20,
    "delta": 1.0,
    "sw": 1.0,
    "sm": 1.5,
    "mu": 0.1,
    "fracAdv": 0.5,
    "w0": 20,
}

FIXTURES = [
    {
        "name": "immigration_death",
        "model": "immigration-death.rules.json",
        "times": [0.5, 1.0, 3.0],
        "law": lambda times: law_immigration_death(10.0, 1.0, times),
        "wrong": ("k_in = 10.5", lambda times: law_immigration_death(10.5, 1.0, times)),
    },
    {
        "name": "immigration_dimerization",
        "model": "immigration-dimerization.rules.json",
        "times": [0.5, 1.0, 3.0],
        "law": lambda times: law_dimerization(20.0, 0.1, times),
        "wrong": ("k2 = 0.11", lambda times: law_dimerization(20.0, 0.11, times)),
        "network": True,
    },
    {
        "name": "mtdna_fixed_population",
        "model": "mtdna-fixed.rules.json",
        "times": [1.0, 3.0, 10.0],
        "law": lambda times: law_mtdna(MTDNA, times),
        "wrong": ("sm = 1.35", lambda times: law_mtdna({**MTDNA, "sm": 1.35}, times)),
        # Supplementary: each molecule weighted by its own rate in D instead of
        # Kowald's w + sm (a + b).
        "also_wrong": [
            ("D = w + sw a + sm b", lambda times: law_mtdna(MTDNA, times, kowald=False))
        ],
    },
    {
        "name": "aggregation_constant",
        "model": "aggregation-constant.rules.json",
        "times": [0.02, 0.1, 0.5],
        "law": lambda times: law_aggregation(constant, 1.0, 12, times),
        "wrong": ("k = 1.1", lambda times: law_aggregation(constant, 1.1, 12, times)),
    },
    {
        "name": "aggregation_multiplicative",
        "model": "aggregation-multiplicative.rules.json",
        "times": [0.02, 0.1, 0.5],
        "law": lambda times: law_aggregation(product, 1.0, 12, times),
        "wrong": ("k = 1.1", lambda times: law_aggregation(product, 1.1, 12, times)),
    },
]


# ---------------------------------------------------------------------------
# Engines


def invoke(binary, *args):
    process = subprocess.run(
        [str(binary), *map(str, args)], capture_output=True, text=True, check=False
    )
    return process.returncode, load(process.stdout)


def rules_run(binary, workdir, model_name, times, replicates, threads, validators):
    source = ROOT / "models/rules" / model_name
    model = workdir / model_name
    model.write_bytes(source.read_bytes())
    request = {
        "schema": "crn-craft.rules-run/0.1",
        "model": {
            "path": model_name,
            "sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
        },
        "time": {"end": times[-1], "report_times": times},
        "ensemble": {"replicates": replicates, "root_seed": SEED},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": replicates,
            "keep_instances": 0,
        },
        "budget": {
            "max_events_per_replicate": 10_000_000,
            "max_sampling_steps_per_replicate": 100_000_000,
            "max_instances": 100_000,
            "threads": threads,
            "max_output_bytes": 1 << 30,
        },
    }
    path = workdir / (model_name + f".{threads}.run.json")
    path.write_text(json.dumps(request))
    validators["rules-run"].validate(request)
    code, result = invoke(binary, "rules", model, path)
    validators["rules-result"].validate(result)
    require(
        code == 0 and result["termination"] == "completed",
        f"{model_name}: {result['termination']}",
    )
    return result


def network_dimerization(binary, workdir, times, replicates):
    model = {
        "schema": "crn-craft.model/0.1",
        "name": "immigration_dimerization",
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": {"k_in": 20.0, "k2": 0.1},
        "species": [{"name": "X", "initial_count": 0}],
        "reactions": [
            {
                "name": "immigration",
                "reactants": {},
                "products": {"X": 1},
                "rate": {"kind": "mass_action", "coefficient": "k_in"},
            },
            {
                "name": "annihilation",
                "reactants": {"X": 2},
                "products": {},
                "rate": {"kind": "mass_action", "coefficient": "k2"},
            },
        ],
    }
    path = workdir / "dimerization.network.json"
    path.write_text(json.dumps(model))
    # Since 0.7.0 (M8a): the network model is converted and run on the unified
    # engine as Static Direct with inverse waits and the linear selector, which
    # replaces the network engine's `ssa` command of 0.6.0.
    converted = subprocess.run(
        [str(binary), "convert", str(path)],
        capture_output=True,
        check=False,
    )
    require(converted.returncode == 0, "convert of the network dimerization")
    path = workdir / "dimerization.model.json"
    path.write_bytes(converted.stdout)
    request = {
        "schema": "crn-craft.run/1.0",
        "model": {
            "path": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        },
        "method": {
            "kind": "static_direct",
            "selector": "linear",
            "waiting_time": "inverse",
        },
        "time": {"end": times[-1], "report_times": [0.0, *times]},
        "ensemble": {"replicates": replicates, "root_seed": SEED},
        "outputs": {
            "statistics": ["mean"],
            "keep_paths": replicates,
            "keep_instances": 0,
        },
        "budget": {
            "max_events_per_replicate": 10_000_000,
            "max_sampling_steps_per_replicate": 100_000_000,
            "max_instances": 16,
            "threads": 4,
            "max_output_bytes": 1 << 30,
        },
    }
    rpath = workdir / "dimerization.run.json"
    rpath.write_text(json.dumps(request))
    code, result = invoke(binary, "run", path, rpath)
    require(
        code == 0 and result["termination"] == "completed", "network SSA incomplete"
    )
    column = result["observables"].index("X")
    counts = np.array(
        [[row[column] for row in p["values"]] for p in result["paths"]],
        dtype=np.int64,
    )
    # Drop the t = 0 row; one species.
    return {"X": counts[:, 1:]}


def ecdf_distance(samples, cdf):
    """sup_k |F_n(k) - F(k)| over the integer support of F and the samples."""
    top = max(len(cdf) - 1, int(samples.max()))
    counts = np.bincount(samples, minlength=top + 1)
    empirical = np.cumsum(counts) / len(samples)
    exact = np.concatenate([cdf, np.full(top + 1 - len(cdf), cdf[-1])])
    return float(np.max(np.abs(empirical - exact)))


def compare(label, columns, cdfs, times, leak, rows, epsilon=EPSILON):
    for name, rows_cdf in cdfs.items():
        for i, t in enumerate(times):
            samples = columns[name][:, i]
            require(
                np.all(samples == np.round(samples)) and samples.min() >= 0,
                "integer observables",
            )
            distance = ecdf_distance(samples.astype(np.int64), rows_cdf[i])
            rows.append(
                {
                    "check": f"{label}:{name}@{t}",
                    "distance": distance,
                    "bound": epsilon + leak,
                    "passed": distance <= epsilon + leak,
                }
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--replicates",
        type=int,
        default=N,
        help="for smoke runs only; gates use 100000",
    )
    args = parser.parse_args()
    require(
        args.replicates == N or args.output.name.startswith("smoke"),
        "gates need 100000 replicates",
    )

    # Smoke runs scale the bound to their size; the gate run uses EPSILON.
    epsilon = math.sqrt(math.log(2 * K / ALPHA) / (2 * args.replicates))
    report = {
        "alpha": ALPHA,
        "comparisons": K,
        "replicates": args.replicates,
        "epsilon": epsilon,
        "seed": SEED,
        "checks": [],
        "power": [],
        "laws": [],
    }
    validators = {}
    for kind in ("rules-model", "rules-check", "rules-run", "rules-result"):
        code, schema = invoke(args.binary, "describe", kind)
        Draft7Validator.check_schema(schema)
        validators[kind] = Draft7Validator(schema)
        snapshot = load((ROOT / f"schemas/{kind}.schema.json").read_text())
        report["checks"].append(
            {
                "check": f"schema_snapshot_{kind}",
                "passed": code == 0 and snapshot == schema,
            }
        )

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for fixture in FIXTURES:
            times = fixture["times"]
            model_path = ROOT / "models/rules" / fixture["model"]
            validators["rules-model"].validate(load(model_path.read_text()))
            code, check = invoke(args.binary, "rules-check", model_path)
            validators["rules-check"].validate(check)
            report["checks"].append(
                {
                    "check": f"rules_check_{fixture['name']}",
                    "passed": code == 0 and check["valid"],
                }
            )

            cdfs, leak = fixture["law"](times)
            report["laws"].append({"fixture": fixture["name"], "leakage": leak})
            result = rules_run(
                args.binary,
                work,
                fixture["model"],
                times,
                args.replicates,
                4,
                validators,
            )
            names = result["observables"]
            values = np.array([p["values"] for p in result["paths"]])
            columns = {name: values[:, :, names.index(name)] for name in cdfs}
            compare(
                f"rules:{fixture['name']}",
                columns,
                cdfs,
                times,
                leak,
                report["checks"],
                epsilon,
            )

            # Wrong law: the same samples must be rejected somewhere.
            for label, wrong in [fixture["wrong"], *fixture.get("also_wrong", [])]:
                wrong_cdfs, _ = wrong(times)
                rows = []
                compare("power", columns, wrong_cdfs, times, 0.0, rows, epsilon)
                rejected = any(not r["passed"] for r in rows)
                report["power"].append(
                    {
                        "fixture": fixture["name"],
                        "wrong_law": label,
                        "max_distance": max(r["distance"] for r in rows),
                        "rejected": rejected,
                    }
                )
                report["checks"].append(
                    {"check": f"power:{fixture['name']}:{label}", "passed": rejected}
                )

            if fixture.get("network"):
                columns = network_dimerization(
                    args.binary, work, times, args.replicates
                )
                compare(
                    "network:immigration_dimerization",
                    columns,
                    cdfs,
                    times,
                    leak,
                    report["checks"],
                    epsilon,
                )

            if fixture["name"] == "aggregation_multiplicative":
                # Ensembles replay identically with one and with four threads.
                small = 2000
                one = rules_run(
                    args.binary, work, fixture["model"], times, small, 1, validators
                )
                four = rules_run(
                    args.binary, work, fixture["model"], times, small, 4, validators
                )
                same = all(
                    one[k] == four[k]
                    for k in ("replicates", "paths", "statistics", "events")
                )
                report["checks"].append({"check": "replay_threads_1_4", "passed": same})
            print(
                fixture["name"],
                sum(c["passed"] for c in report["checks"]),
                "passed",
                flush=True,
            )

    comparisons = [c for c in report["checks"] if "@" in c["check"]]
    report["checks"].append(
        {"check": "declared_comparison_count", "passed": len(comparisons) == K}
    )
    report["passed"] = sum(c["passed"] for c in report["checks"])
    report["failed"] = sum(not c["passed"] for c in report["checks"])
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"passed": report["passed"], "failed": report["failed"]}))
    return 1 if report["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
