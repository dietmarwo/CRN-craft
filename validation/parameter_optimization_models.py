"""Frozen parameter-fit models and independent exact moment targets."""

import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "validation/parameter-optimization"


def reaction(name, reactants, products, coefficient):
    return dict(
        name=name,
        reactants=reactants,
        products=products,
        rate=dict(kind="mass_action", coefficient=coefficient),
    )


def model(name, species, parameters, reactions):
    return dict(
        schema="crn-craft.model/0.1",
        name=name,
        convention="stochastic_counts",
        omega=1.0,
        parameters=parameters,
        species=[dict(name=n, initial_count=x) for n, x in species],
        reactions=reactions,
    )


def fixtures():
    small = model(
        "ParameterBirthDeath",
        [("X", 0)],
        {"birth": 20.0, "death": 0.6},
        [
            reaction("birth", {}, {"X": 1}, "birth"),
            reaction("death", {"X": 1}, {}, "death"),
        ],
    )
    species = [(f"X{i}", 0) for i in range(33)]
    reactions = []
    for name, _ in species:
        reactions.extend(
            [
                reaction("birth_" + name, {}, {name: 1}, "birth"),
                reaction("death_" + name, {name: 1}, {}, "death"),
            ]
        )
    sparse = model(
        "ParameterSparse33", species, {"birth": 3.0, "death": 0.2}, reactions
    )
    binding = model(
        "ParameterBinding12",
        [("A", 24), ("B", 12), ("C", 0)],
        {"bind": 0.03, "unbind": 0.4},
        [
            reaction("bind", {"A": 1, "B": 1}, {"C": 1}, "bind"),
            reaction("unbind", {"C": 1}, {"A": 1, "B": 1}, "unbind"),
        ],
    )
    return [
        dict(case="birth_death", model=small, times=[10.0 * i / 5 for i in range(6)]),
        dict(case="sparse33", model=sparse, times=[20.0 * i / 5 for i in range(6)]),
        dict(case="binding12", model=binding, times=[2.0 * i / 5 for i in range(6)]),
    ]


def exact_moments(fixture, values):
    times = np.asarray(fixture["times"])
    if fixture["case"] != "binding12":
        birth, death = values
        mean = (birth / death * -np.expm1(-death * times))[:, None]
        mean = np.repeat(mean, len(fixture["model"]["species"]), axis=1)
        return dict(mean=mean.tolist(), variance=mean.tolist())
    bind, unbind = values
    q = np.zeros((13, 13))
    for n in range(13):
        if n < 12:
            q[n, n + 1] = bind * (24 - n) * (12 - n)
        if n:
            q[n, n - 1] = unbind * n
        q[n, n] = -q[n].sum()
    mean, variance = [], []
    x = np.arange(13.0)
    for t in times:
        p = expm(q * t)[0]
        if np.min(p) < -1e-12 or abs(p.sum() - 1.0) > 1e-11:
            raise RuntimeError("Independent finite-state oracle failed")
        m = float(p @ x)
        v = float(p @ (x * x) - m * m)
        mean.append([24.0 - m, 12.0 - m, m])
        variance.append([v] * 3)
    return dict(mean=mean, variance=variance)


def loss(actual, target):
    total = 0.0
    for key in ["mean", "variance"]:
        a, t = np.asarray(actual[key])[1:], np.asarray(target[key])[1:]
        if a.shape != t.shape or not np.all(np.isfinite(a)):
            raise RuntimeError("Invalid objective moments")
        total += float(np.mean(((a - t) / (1.0 + t)) ** 2)) / 2.0
    return total


def freeze():
    rows = []
    for f in fixtures():
        names = sorted(f["model"]["parameters"])
        truth = [f["model"]["parameters"][name] for name in names]
        target = exact_moments(f, truth)
        row = dict(
            case=f["case"],
            parameter_names=names,
            true_values=truth,
            times=f["times"],
            training_paths=32,
            heldout_paths=4096,
            training_seed=77201,
            heldout_seed=880431,
            guess=(np.asarray(truth) * [1.7, 0.6]).tolist(),
            lower=(np.asarray(truth) / 4.0).tolist(),
            upper=(np.asarray(truth) * 4.0).tolist(),
            target=target,
        )
        data = (json.dumps(f["model"], sort_keys=True, indent=2) + "\n").encode()
        path = FOLDER / "models" / f"{f['case']}.model.json"
        if path.exists() and path.read_bytes() != data:
            raise RuntimeError("Frozen model changed")
        path.write_bytes(data)
        row["model_sha256"] = hashlib.sha256(data).hexdigest()
        rows.append(row)
    payload = dict(
        schema="crn-craft.parameter-fit-fixtures/0.1",
        rows=rows,
        optimizer="fcmaes-rust 0.2.1 BiteOpt",
        evaluations=128,
        optimizer_seeds=[20261102, 20261103, 20261104],
    )
    path = FOLDER / "fixtures.json"
    data = json.dumps(payload, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != data:
        raise RuntimeError("Frozen targets changed")
    path.write_text(data)
    return payload


if __name__ == "__main__":
    print(json.dumps(freeze(), indent=2))


def guard_fixture():
    return model(
        "GeneratedArithmeticGuard",
        [("P", 8), ("Q", 5), ("D", 0)],
        {"rate": 0.1},
        [
            reaction("birth", {}, {"P": 1}, "rate"),
            reaction("first", {"P": 1}, {}, "rate"),
            reaction("dimer", {"P": 2}, {"D": 1}, "rate"),
            reaction("pair", {"P": 1, "Q": 1}, {"D": 1}, "rate"),
            reaction("higher", {"P": 3}, {"D": 1}, "rate"),
        ],
    )
