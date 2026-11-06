"""Python builder and whole-rule propensities: exact means and rebop.

The protocol is declared in validation/builder/RESULTS.md: 18 z tests on
means at fixed times, 4,000 replicates per arm, one Bonferroni family at
alpha 0.001. rebop 0.10.0 is optional (its rows are skipped without it) and
is not a dependency of the package. The report holds fixed-work evidence
only and is identical under python -O.
"""

import argparse
import json
import math
from pathlib import Path
from statistics import NormalDist
import tempfile

import numpy as np

from crn_craft import rules

ROOT = Path(__file__).resolve().parents[1]
N = 4000
ALPHA = 0.001
FAMILY = 18
Z_CRIT = NormalDist().inv_cdf(1 - ALPHA / (2 * FAMILY))

SIR_HAND = {
    "schema": "crn-craft.rules-model/0.3",
    "name": "sir",
    "parameters": {},
    "types": [{"name": "S"}, {"name": "I"}, {"name": "R"}],
    "initial": [{"type": "S", "count": 490}, {"type": "I", "count": 10}],
    "observables": {"S": {"type": "S"}, "I": {"type": "I"}, "R": {"type": "R"}},
    "rules": [
        {
            "name": "infection",
            "reactants": [{"label": "x", "type": "S"}, {"label": "y", "type": "I"}],
            "rate": "0.0005",
            "products": [{"type": "I"}, {"type": "I"}],
        },
        {
            "name": "recovery",
            "reactants": [{"label": "x", "type": "I"}],
            "rate": "0.1",
            "products": [{"type": "R"}],
        },
    ],
}

BIRTH_DEATH_HAND = {
    "schema": "crn-craft.rules-model/0.3",
    "name": "birth_death",
    "parameters": {"k": 10.0, "g": 1.0},
    "types": [{"name": "A"}],
    "initial": [],
    "observables": {"A": {"type": "A"}},
    "rules": [
        {
            "name": "birth",
            "reactants": [],
            "propensity": "k",
            "products": [{"type": "A"}],
        },
        {
            "name": "death",
            "reactants": [{"label": "x", "type": "A"}],
            "propensity": "g * A",
            "products": [],
        },
    ],
}


def sir_builder(strings):
    model = rules.Model("sir")
    model.add_species("S", 490).add_species("I", 10).add_species("R")
    if strings:
        model.add_parameter("beta", 5e-4).add_parameter("gamma", 0.1)
        model.add_reaction("beta * S * I", ["S", "I"], ["I", "I"], name="infection")
        model.add_reaction("gamma * I", ["I"], ["R"], name="recovery")
    else:
        model.add_reaction(5e-4, ["S", "I"], ["I", "I"], name="infection")
        model.add_reaction(0.1, ["I"], ["R"], name="recovery")
    return model


def birth_death_builder(strings):
    model = rules.Model("birth_death")
    if strings:
        model.add_parameter("k", 10.0).add_parameter("g", 1.0)
        model.add_reaction("k", [], ["A"], name="birth")
        model.add_reaction("g * A", ["A"], [], name="death")
    else:
        model.add_reaction(10.0, [], ["A"], name="birth")
        model.add_reaction(1.0, ["A"], [], name="death")
    return model


def simulate(model, end, times, seed):
    """Replicate × time × observable values, and the observable names."""
    run = rules.ensemble(
        model,
        end=end,
        report_times=times,
        replicates=N,
        seed=seed,
        threads=2,
        keep_paths=N,
        method="static_direct",
    )
    result = run.simulate()
    meta = result.metadata
    if meta["termination"] != "completed":
        raise RuntimeError(f"{model} ended with {meta['termination']}")
    paths = np.array([p["values"] for p in meta["paths"]], dtype=float)
    return paths, meta["observables"]


def summary(values):
    values = np.asarray(values, dtype=float)
    return float(values.mean()), float(values.var(ddof=1)), len(values)


def test_exact(name, values, mean, variance=None):
    m, v, n = summary(values)
    se = math.sqrt((v if variance is None else variance) / n)
    z = (m - mean) / se
    return {
        "test": name,
        "kind": "exact mean",
        "estimate": m,
        "expected": mean,
        "standard_error": se,
        "z": z,
        "passed": abs(z) <= Z_CRIT,
    }


def test_welch(name, first, second):
    m1, v1, n1 = summary(first)
    m2, v2, n2 = summary(second)
    se = math.sqrt(v1 / n1 + v2 / n2)
    z = (m1 - m2) / se if se > 0 else (0.0 if m1 == m2 else math.inf)
    return {
        "test": name,
        "kind": "two-sample mean (Welch)",
        "estimate": m1,
        "expected": m2,
        "standard_error": se,
        "z": z,
        "passed": abs(z) <= Z_CRIT,
    }


def skipped(name):
    return {"test": name, "kind": "skipped: rebop is not installed", "passed": None}


def rebop_runs(build, init, end, steps, seed):
    import rebop

    sim = rebop.Gillespie()
    build(sim)
    rng = np.random.default_rng(seed)
    runs = [sim.run(init, tmax=end, nb_steps=steps, rng=rng) for _ in range(N)]
    return runs, rebop.__version__


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks = []

    def check(name, condition):
        checks.append({"check": name, "passed": bool(condition)})

    check("sir_numeric_json", sir_builder(False).to_dict() == SIR_HAND)
    check(
        "birth_death_string_json",
        birth_death_builder(True).to_dict() == BIRTH_DEATH_HAND,
    )
    with tempfile.TemporaryDirectory(prefix="crn-builder-") as folder:
        path = Path(folder) / "sir.rules.json"
        path.write_text(sir_builder(True).to_json())
        # rules.load raises CrnError for an invalid model.
        check("builder_model_validates", rules.load(path) == path.read_bytes())

    tests = []
    # 1-2: birth-death against Poisson(10 (1 - e^-10)).
    mean = 10.0 * (1.0 - math.exp(-10.0))
    for i, strings in enumerate((False, True)):
        paths, names = simulate(birth_death_builder(strings), 10.0, [10.0], 201 + i)
        label = "string" if strings else "numeric"
        tests.append(
            test_exact(f"birth_death_{label}_mean_A_10", paths[:, 0, 0], mean, mean)
        )

    # 3: single-enzyme Michaelis-Menten with mass-action constants.
    kon, koff, kcat, enzymes, t = 2.0, 0.5, 1.0, 5, 20.0
    mm = rules.Model("enzyme")
    mm.add_species("E", enzymes)
    mm.add_reaction(kon, ["E"], ["C"], name="binding")
    mm.add_reaction(koff, ["C"], ["E"], name="unbinding")
    mm.add_reaction(kcat, ["C"], ["E", "P"], name="catalysis")
    paths, names = simulate(mm, t, [t], 203)
    lam = kon + koff + kcat
    exact = enzymes * kcat * (kon / lam) * (t - (1.0 - math.exp(-lam * t)) / lam)
    tests.append(
        test_exact(
            "michaelis_menten_mass_action_mean_P_20",
            paths[:, 0, names.index("P")],
            exact,
        )
    )

    # 4: Michaelis-Menten flux as a whole-rule propensity on S -> S + P.
    vmax, km, s, t = 4.0, 25.0, 50, 10.0
    flux = rules.Model("flux")
    flux.add_parameter("vmax", vmax).add_parameter("Km", km).add_species("S", s)
    flux.add_reaction("vmax * S / (Km + S)", ["S"], ["S", "P"], name="turnover")
    paths, names = simulate(flux, t, [t], 204)
    mean = t * vmax * s / (km + s)
    tests.append(
        test_exact(
            "michaelis_menten_propensity_mean_P_10",
            paths[:, 0, names.index("P")],
            mean,
            mean,
        )
    )

    # 5-16: SIR, hand-written mass-action JSON against the builder's string
    # propensities and against rebop.
    times = [10.0, 20.0, 40.0]
    hand, names = simulate(SIR_HAND, 40.0, times, 101)
    built, built_names = simulate(sir_builder(True), 40.0, times, 102)
    if built_names != names:
        raise RuntimeError("observable order differs")
    for k, t in enumerate(times):
        for species in ("S", "I"):
            j = names.index(species)
            tests.append(
                test_welch(
                    f"sir_json_vs_builder_{species}_{t:g}",
                    hand[:, k, j],
                    built[:, k, j],
                )
            )
    try:
        runs, version = rebop_runs(
            lambda g: (
                g.add_reaction(5e-4, ["S", "I"], ["I", "I"]),
                g.add_reaction(0.1, ["I"], ["R"]),
            ),
            {"S": 490, "I": 10},
            40.0,
            4,
            1001,
        )
    except ImportError:
        runs, version = None, None
    for k, t in enumerate(times):
        for species in ("S", "I"):
            name = f"sir_json_vs_rebop_{species}_{t:g}"
            if runs is None:
                tests.append(skipped(name))
                continue
            row = [float(r["time"].values[int(t // 10)]) for r in runs[:1]]
            if row != [t]:
                raise RuntimeError("rebop time grid mismatch")
            other = [float(r[species].values[int(t // 10)]) for r in runs]
            tests.append(test_welch(name, hand[:, k, names.index(species)], other))

    # 17-18: identical reactants. Builder rate 2k equals rebop's k.
    dimer = rules.Model("dimer")
    dimer.add_species("A", 100).add_reaction(0.02, ["A", "A"], [], name="pair")
    paths, names = simulate(dimer, 1.0, [0.5, 1.0], 301)
    try:
        runs, version = rebop_runs(
            lambda g: g.add_reaction(0.01, ["A", "A"], []), {"A": 100}, 1.0, 2, 1002
        )
    except ImportError:
        runs = None
    for k, t in enumerate((0.5, 1.0)):
        name = f"dimer_builder_2k_vs_rebop_k_A_{t:g}"
        if runs is None:
            tests.append(skipped(name))
            continue
        other = [float(r["A"].values[k + 1]) for r in runs]
        tests.append(test_welch(name, paths[:, k, 0], other))

    if len(tests) != FAMILY:
        raise RuntimeError(f"{len(tests)} tests, declared {FAMILY}")
    failed = [t["test"] for t in tests if t["passed"] is False]
    report = {
        "schema": "crn-craft.builder-validation/0.1",
        "protocol": "validation/builder/RESULTS.md",
        "replicates_per_arm": N,
        "family_alpha": ALPHA,
        "family_size": FAMILY,
        "z_critical": Z_CRIT,
        "rebop": version,
        "checks": checks,
        "tests": tests,
        "failed": failed + [c["check"] for c in checks if not c["passed"]],
        "skipped": [t["test"] for t in tests if t["passed"] is None],
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(text)
    print(text, end="")
    if report["failed"]:
        raise SystemExit("builder validation failed: " + ", ".join(report["failed"]))


if __name__ == "__main__":
    main()
