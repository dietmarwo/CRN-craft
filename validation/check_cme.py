"""Direct SSA marginals against exact CME laws for nonlinear and custom rates.

The rule is frozen in docs/decisions/iteration-04-review-cme.md. The exact
laws come from uniformization of the master-equation generator, with
propensities written here independently of the Rust code. Standard library
only; the report is fixed-work evidence and identical under ``python -O``.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
N = 100_000
TIMES = [0.0, 0.5, 2.0, 10.0]
ALPHA = 0.001
TESTS = 12
EPSILON = math.sqrt(math.log(2 * TESTS / ALPHA) / (2 * N))


def mass_action(c):
    return {"kind": "mass_action", "coefficient": c}


def reaction(name, reactants, products, rate):
    return {"name": name, "reactants": reactants, "products": products, "rate": rate}


FIXTURES = {
    "F1_dimerization": {
        "species": [("A", 20), ("B", 0)],
        "parameters": {},
        "reactions": [
            reaction("dimerize", {"A": 2}, {"B": 1}, mass_action(0.1)),
            reaction("split", {"B": 1}, {"A": 2}, mass_action(1.0)),
            reaction("null", {"B": 1}, {"B": 1}, mass_action(2.0)),
        ],
        "marginal": "B",
        "seed": 2026110451,
        "power": True,
    },
    "F2_heterodimer": {
        "species": [("A", 15), ("B", 10), ("C", 0)],
        "parameters": {},
        "reactions": [
            reaction("bind", {"A": 1, "B": 1}, {"C": 1}, mass_action(0.05)),
            reaction("unbind", {"C": 1}, {"A": 1, "B": 1}, mass_action(0.5)),
        ],
        "marginal": "C",
        "seed": 2026110452,
    },
    "F3_hill_mm": {
        "species": [("S", 30), ("P", 0)],
        "parameters": {
            "vmax": 5.0,
            "km": 20.0,
            "basal": 0.2,
            "hill_vmax": 3.0,
            "k": 10.0,
            "n": 2.5,
        },
        "reactions": [
            reaction(
                "convert",
                {"S": 1},
                {"P": 1},
                {
                    "kind": "michaelis_menten",
                    "substrate": "S",
                    "vmax": "vmax",
                    "km": "km",
                },
            ),
            reaction(
                "recycle",
                {"P": 1},
                {"S": 1},
                {
                    "kind": "hill_activation",
                    "regulator": "P",
                    "vmax": "hill_vmax",
                    "basal": "basal",
                    "k": "k",
                    "n": "n",
                },
            ),
        ],
        "marginal": "P",
        "seed": 2026110453,
    },
    "F4_open_annihilation": {
        "species": [("X", 0)],
        "parameters": {},
        "reactions": [
            reaction("birth", {}, {"X": 1}, mass_action(10.0)),
            reaction("annihilate", {"X": 2}, {}, mass_action(0.2)),
        ],
        "marginal": "X",
        "seed": 2026110454,
        "cap": 80,
        "power": True,
    },
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def propensity(rx, counts, parameters, falling):
    for species, order in rx["reactants"].items():
        if counts[species] < order:
            return 0.0
    rate = rx["rate"]
    if rate["kind"] == "mass_action":
        value = rate["coefficient"]
        for species, order in rx["reactants"].items():
            factor = math.perm if falling else math.comb
            value *= factor(counts[species], order)
        return value
    if rate["kind"] == "michaelis_menten":
        s = counts[rate["substrate"]]
        return parameters[rate["vmax"]] * s / (parameters[rate["km"]] + s)
    if rate["kind"] == "hill_activation":
        x = counts[rate["regulator"]]
        k, n = parameters[rate["k"]], parameters[rate["n"]]
        induced = parameters[rate["vmax"]] * x**n / (k**n + x**n)
        return parameters[rate["basal"]] + induced
    raise ValueError(rate["kind"])


def generator(fixture, falling=False):
    """Reachable states, transitions and exit rates; capped states leak mass."""
    names = [name for name, _ in fixture["species"]]
    start = tuple(count for _, count in fixture["species"])
    cap = fixture.get("cap")
    index, states, frontier = {start: 0}, [start], [start]
    edges, exits = [], {}
    while frontier:
        found = []
        for state in frontier:
            counts = dict(zip(names, state, strict=True))
            out = 0.0
            for rx in fixture["reactions"]:
                change = [
                    rx["products"].get(s, 0) - rx["reactants"].get(s, 0) for s in names
                ]
                a = propensity(rx, counts, fixture["parameters"], falling)
                if a == 0.0 or not any(change):
                    continue
                out += a
                target = tuple(x + d for x, d in zip(state, change, strict=True))
                if cap is not None and max(target) > cap:
                    continue
                if target not in index:
                    index[target] = len(states)
                    states.append(target)
                    found.append(target)
                edges.append((index[state], index[target], a))
            exits[index[state]] = out
        frontier = found
    return states, edges, [exits[i] for i in range(len(states))]


def uniformize(states, edges, exits):
    """p(t) at each save time; Poisson weights with a geometric tail bound."""
    rate = max(exits) * 1.02
    p = [1.0] + [0.0] * (len(states) - 1)
    laws, now = [], 0.0
    for t in TIMES:
        mu = rate * (t - now)
        now = t
        if mu == 0.0:
            laws.append(p[:])
            continue
        total = [0.0] * len(states)
        vector, k, log_weight = p[:], 0, -mu
        while True:
            weight = math.exp(log_weight)
            for i, value in enumerate(vector):
                total[i] += weight * value
            if k + 1 > mu and weight * (k + 1) / (k + 1 - mu) < 1e-17:
                break
            require(k < mu + 50 * math.sqrt(mu) + 1000, "uniformization diverged")
            step = [v * (1.0 - e / rate) for v, e in zip(vector, exits, strict=True)]
            for i, j, a in edges:
                step[j] += vector[i] * a / rate
            vector, k = step, k + 1
            log_weight += math.log(mu) - math.log(k)
        p = total
        laws.append(p[:])
    return laws


def marginal(states, position, law):
    distribution = {}
    for state, mass in zip(states, law, strict=True):
        distribution[state[position]] = distribution.get(state[position], 0.0) + mass
    return distribution


def max_cdf_error(samples, distribution):
    counts = {}
    for value in samples:
        counts[value] = counts.get(value, 0) + 1
    empirical = exact = worst = 0.0
    for value in sorted(set(distribution) | set(counts)):
        empirical += counts.get(value, 0) / len(samples)
        exact += distribution.get(value, 0.0)
        worst = max(worst, abs(empirical - exact))
    return worst


def invoke(binary, *args, expected=0):
    p = subprocess.run(
        [binary, *map(str, args)], capture_output=True, timeout=900, check=False
    )
    require(
        p.returncode == expected and not p.stderr,
        f"{args[0]} exit {p.returncode}: {p.stdout[-2000:]!r} {p.stderr[-1000:]!r}",
    )
    return json.loads(p.stdout)


def simulate(binary, name, fixture, folder, command="simulate"):
    model = {
        "schema": "crn-craft.model/0.1",
        "name": name,
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": fixture["parameters"],
        "species": [{"name": s, "initial_count": c} for s, c in fixture["species"]],
        "reactions": fixture["reactions"],
    }
    model_bytes = (json.dumps(model, indent=2) + "\n").encode()
    (folder / "model.json").write_bytes(model_bytes)
    run = {
        "schema": "crn-craft.run/0.1",
        "model": {
            "path": "model.json",
            "sha256": hashlib.sha256(model_bytes).hexdigest(),
        },
        "method": {"kind": "ssa", "aggregator": "direct"},
        "ensemble": {"replicates": N, "root_seed": fixture["seed"]},
        "time": {"end": TIMES[-1], "save_times": TIMES},
        "outputs": {"statistics": ["mean"], "keep_trajectories": N},
        "budget": {
            "max_events_per_replicate": 1_000_000,
            "threads": 8,
            "max_output_bytes": 268_435_456,
        },
    }
    (folder / "run.json").write_text(json.dumps(run, indent=2) + "\n")
    out = folder / "archive"
    report = invoke(binary, command, folder / "model.json", folder / "run.json", out)
    require(report["valid"] and report["simulation_complete"], "simulation incomplete")
    require(invoke(binary, "verify", out)["valid"], "archive verification failed")
    result_bytes = (out / "result.json").read_bytes()
    return json.loads(result_bytes), hashlib.sha256(result_bytes).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--command", choices=["simulate", "simulate-cached"], default="simulate"
    )
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    rows, power, fixtures = [], [], []
    try:
        with tempfile.TemporaryDirectory(prefix="crn-cme-") as tmp:
            for name, fixture in FIXTURES.items():
                folder = Path(tmp) / name
                folder.mkdir()
                result, result_hash = simulate(
                    binary, name, fixture, folder, args.command
                )
                names = [s for s, _ in fixture["species"]]
                position = names.index(fixture["marginal"])
                states, edges, exits = generator(fixture)
                laws = uniformize(states, edges, exits)
                wrong = None
                if fixture.get("power"):
                    wrong_states, wrong_edges, wrong_exits = generator(fixture, True)
                    wrong = (
                        wrong_states,
                        uniformize(wrong_states, wrong_edges, wrong_exits),
                    )
                fixtures.append(
                    {
                        "fixture": name,
                        "states": len(states),
                        "events": result["events"],
                        "result_sha256": result_hash,
                    }
                )
                for k, t in enumerate(TIMES):
                    samples = [
                        path["counts"][k][position] for path in result["trajectories"]
                    ]
                    if t == 0.0:
                        initial = fixture["species"][position][1]
                        require(all(s == initial for s in samples), "bad initial state")
                        continue
                    distribution = marginal(states, position, laws[k])
                    error = max_cdf_error(samples, distribution)
                    rows.append(
                        {
                            "check": f"{name}_t{t:g}",
                            "max_cdf_error": error,
                            "exact_mean": sum(v * m for v, m in distribution.items()),
                            "sample_mean": result["statistics"]["mean"][k][position],
                            "exact_mass": sum(laws[k]),
                            "passed": error <= EPSILON,
                        }
                    )
                    if wrong:
                        wrong_error = max_cdf_error(
                            samples, marginal(wrong[0], position, wrong[1][k])
                        )
                        power.append({"fixture": name, "time": t, "error": wrong_error})
        for name in ("F1_dimerization", "F4_open_annihilation"):
            rejected = any(p["error"] > EPSILON for p in power if p["fixture"] == name)
            rows.append(
                {"check": f"{name}_rejects_falling_factorial_law", "passed": rejected}
            )
    except Exception as error:
        rows.append({"check": "execution_error", "passed": False, "error": str(error)})
    report = {
        "schema": "crn-craft.cme-validation/0.1",
        "replicates_per_fixture": N,
        "positive_time_cdf_tests": TESTS,
        "family_alpha": ALPHA,
        "dkw_epsilon": EPSILON,
        "fixtures": fixtures,
        "power": power,
        "passed": sum(row["passed"] for row in rows),
        "failed": sum(not row["passed"] for row in rows),
        "checks": rows,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    if report["failed"]:
        raise RuntimeError("CME acceptance failed; failed rows retained")


if __name__ == "__main__":
    main()
