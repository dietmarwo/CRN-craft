"""Independent exact/Decimal oracle for M1a; no numerical Rust code is reused.

Since 0.7.0 (M5b) `crn-craft analyze` reads the network models through the
converter and reports `crn-craft.analysis/1.0`. Three expectations follow
the rules engine (validation/merge/RESULTS.md, M5b): a network without
reactions is refused with `CRN-R016`; a propensity that overflows f64 fails
with `CRN-R030`; a positive propensity that underflows is reported as 0 and
checked against the exact value with the same tolerance. A Hill law whose
induced term the overflow-free converted form evaluates to 0 (M2) is
accepted within the absolute floor vmax * 2^-1020 that M2 declared, and
listed in the report.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import tempfile

from jsonschema import Draft7Validator

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261101
RTOL = 2e-12
ATOL = 8 * math.ulp(0.0)
# The absolute floor M2 declared for the converted Hill laws.
HILL_FLOOR = 2.0**-1020


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rref(matrix, width):
    """Rational Gauss-Jordan oracle, independent of Rust's integer elimination."""
    rows = [[Fraction(value) for value in row] for row in matrix]
    pivots = []
    for column in range(width):
        pivot = next(
            (i for i in range(len(pivots), len(rows)) if rows[i][column]), None
        )
        if pivot is None:
            continue
        row = len(pivots)
        rows[row], rows[pivot] = rows[pivot], rows[row]
        divisor = rows[row][column]
        rows[row] = [value / divisor for value in rows[row]]
        for i in range(len(rows)):
            if i != row:
                multiplier = rows[i][column]
                rows[i] = [
                    a - multiplier * b for a, b in zip(rows[i], rows[row], strict=True)
                ]
        pivots.append(column)
    basis = []
    for free in range(width):
        if free not in pivots:
            vector = [Fraction(0)] * width
            vector[free] = Fraction(1)
            for row, pivot in enumerate(pivots):
                vector[pivot] = -rows[row][free]
            basis.append(vector)
    return len(pivots), basis


def scalar_value(model, value):
    """A rate argument as a Decimal: a literal or a parameter's value."""
    if isinstance(value, str):
        value = model["parameters"][value]
    return Decimal.from_float(float(value))


def exact_rate(model, reaction, counts):
    if any(counts[name] < order for name, order in reaction["reactants"].items()):
        return Decimal(0)

    def scalar(value):
        if isinstance(value, str):
            value = model["parameters"][value]
        return Decimal.from_float(float(value))

    rate = reaction["rate"]
    kind = rate["kind"]
    with localcontext() as context:
        context.prec = 100
        if kind == "mass_action":
            combinations = math.prod(
                math.comb(counts[name], order)
                for name, order in reaction["reactants"].items()
            )
            return scalar(rate["coefficient"]) * combinations
        if kind == "michaelis_menten":
            x = Decimal(counts[rate["substrate"]])
            return scalar(rate["vmax"]) * x / (scalar(rate["km"]) + x)
        x = Decimal(counts[rate["regulator"]])
        k, n = scalar(rate["k"]), scalar(rate["n"])
        xn, kn = x**n, k**n
        numerator = kn if kind == "hill_repression" else xn
        return scalar(rate["basal"]) + scalar(rate["vmax"]) * numerator / (kn + xn)


def make_model(species, reactions, parameters=None):
    return {
        "schema": "crn-craft.model/0.1",
        "name": "oracle_fixture",
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": parameters or {},
        "species": [{"name": name, "initial_count": count} for name, count in species],
        "reactions": reactions,
    }


def reaction(name, left, right, rate):
    return {"name": name, "reactants": left, "products": right, "rate": rate}


def mass_action(coefficient):
    return {"kind": "mass_action", "coefficient": coefficient}


def random_model(rng, index):
    names = [f"X{i}" for i in range(rng.randint(2, 7))]
    rng.shuffle(names)
    species = [(name, rng.randint(0, 50)) for name in names]
    reactions = []
    for i in range(rng.randint(0, 10)):
        if index % 2 == 0:
            # Closed conversions exercise nontrivial conservation spaces.
            a, b = rng.sample(names, 2)
            left, right = {a: rng.randint(1, 2)}, {b: rng.randint(1, 2)}
        else:
            left = {name: rng.randint(1, 3) for name in names if rng.random() < 0.3}
            right = {name: rng.randint(1, 3) for name in names if rng.random() < 0.3}
        kind = rng.choice(
            ["mass_action", "hill_activation", "hill_repression", "michaelis_menten"]
        )
        if kind == "mass_action":
            rate = mass_action(rng.choice([0, 0.5, 1, "c"]))
        elif kind == "michaelis_menten":
            rate = {
                "kind": kind,
                "substrate": rng.choice(names),
                "vmax": "v",
                "km": "k",
            }
        else:
            rate = {
                "kind": kind,
                "regulator": rng.choice(names),
                "vmax": "v",
                "basal": 0.25,
                "k": "k",
                "n": "n",
            }
        reactions.append(reaction(f"r{i}", left, right, rate))
    return make_model(species, reactions, {"c": 0.125, "v": 3.0, "k": 5.0, "n": 1.5})


def compare_structure(model, data):
    names = [s["name"] for s in model["species"]]
    require(data["species"] == names, "species index ordering changed")
    indices = {name: i for i, name in enumerate(names)}
    columns, reads = [], []
    for expected, actual in zip(model["reactions"], data["reactions"], strict=True):
        column = [
            expected["products"].get(n, 0) - expected["reactants"].get(n, 0)
            for n in names
        ]
        columns.append(column)
        for label in ("reactants", "products"):
            entries = [
                {"species": indices[name], "coefficient": count}
                for name, count in expected[label].items()
            ]
            require(actual[label] == sorted(entries, key=lambda t: t["species"]), label)
        require(
            actual["net_changes"]
            == [{"species": i, "coefficient": c} for i, c in enumerate(column) if c],
            "net stoichiometric column mismatch",
        )
        names_read = set(expected["reactants"])
        rate = expected["rate"]
        if "regulator" in rate:
            names_read.add(rate["regulator"])
        if "substrate" in rate:
            names_read.add(rate["substrate"])
        read = {indices[name] for name in names_read}
        reads.append(read)
        require(actual["reads"] == sorted(read), "rate or availability read missing")
    graph = [
        [j for j, read in enumerate(reads) if any(column[i] for i in read)]
        for column in columns
    ]
    require(data["dependencies"] == graph, "dependency graph mismatch")
    rank, oracle_basis = rref(columns, len(names))
    conservation = data["conservation"]
    require(conservation["status"] == "exact", "small-network conservation unavailable")
    require(conservation["rank"] == rank, "stoichiometric rank mismatch")
    laws = conservation["basis"]
    basis = [[int(v) for v in law["coefficients"]] for law in laws]
    nullity = len(names) - rank
    require(len(basis) == nullity, "conservation dimension mismatch")
    require(rref(basis, len(names))[0] == nullity, "dependent conservation basis")
    require(
        rref(basis + oracle_basis, len(names))[0] == nullity, "null-space span mismatch"
    )
    for row, law in zip(basis, laws, strict=True):
        for column in columns:
            require(
                sum(w * c for w, c in zip(row, column, strict=True)) == 0,
                "nonzero exact conservation residual",
            )
        total = sum(
            w * s["initial_count"] for w, s in zip(row, model["species"], strict=True)
        )
        require(str(total) == law["initial_total"], "initial conserved total mismatch")
    for i, bound in enumerate(conservation["initial_class_upper_bounds"]):
        candidates = [
            int(law["initial_total"]) // row[i]
            for row, law in zip(basis, laws, strict=True)
            if row[i] > 0 and all(w >= 0 for w in row)
        ]
        require(
            bound == (str(min(candidates)) if candidates else None),
            "unsupported count bound",
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    schema = json.loads(
        subprocess.check_output([str(binary), "describe", "analysis-1.0"])
    )
    Draft7Validator.check_schema(schema)
    validator = Draft7Validator(schema)
    checks = []
    maximum_relative_error = 0.0
    evaluated_rates = 0
    converter_refusals, overflow_failures, underflow_zero, hill_floor = [], [], [], []
    fixture_hash = hashlib.sha256()
    with tempfile.TemporaryDirectory(prefix="crn-network-") as folder:
        path = Path(folder) / "model.json"

        def case(name, model, structure=True):
            nonlocal maximum_relative_error, evaluated_rates
            raw = (json.dumps(model, sort_keys=True) + "\n").encode()
            fixture_hash.update(raw)
            path.write_bytes(raw)
            try:
                result = subprocess.run(
                    [str(binary), "analyze", str(path)],
                    capture_output=True,
                    timeout=15,
                    check=False,
                )
                require(not result.stderr, "unexpected stderr")
                report = json.loads(result.stdout)
                validator.validate(report)
                counts = {s["name"]: s["initial_count"] for s in model["species"]}
                references = [exact_rate(model, r, counts) for r in model["reactions"]]
                overflow = any(not math.isfinite(float(v)) for v in references)
                if not model["reactions"]:
                    # A rules model needs a rule; the converter refuses.
                    require(
                        result.returncode == 1
                        and not report["success"]
                        and report["analysis"] is None,
                        "network without reactions accepted",
                    )
                    require(
                        any(
                            d["code"] == "CRN-R016" and d["path"] == "/reactions"
                            for d in report["diagnostics"]
                        ),
                        "wrong refusal of a network without reactions",
                    )
                    converter_refusals.append(name)
                elif overflow:
                    require(
                        result.returncode == 1
                        and not report["success"]
                        and report["analysis"] is None,
                        "unrepresentable propensity accepted",
                    )
                    require(
                        any(d["code"] == "CRN-R030" for d in report["diagnostics"]),
                        "wrong numeric failure code",
                    )
                    overflow_failures.append(name)
                else:
                    require(
                        result.returncode == 0 and report["success"],
                        f"analysis failed: {report['diagnostics']}",
                    )
                    require(
                        report["simulation_available"] is False,
                        "analysis claims simulation",
                    )
                    require(
                        report["model_sha256"] == hashlib.sha256(raw).hexdigest(),
                        "model hash mismatch",
                    )
                    data = report["analysis"]
                    if structure:
                        compare_structure(model, data)
                    for reaction, actual, exact in zip(
                        model["reactions"],
                        data["propensities"],
                        references,
                        strict=True,
                    ):
                        expected = float(exact)
                        error = abs(actual - expected)
                        if exact > 0 and expected == 0:
                            underflow_zero.append(name)
                        rate = reaction["rate"]
                        if (
                            error > RTOL * abs(expected) + ATOL
                            and rate["kind"] in ("hill_activation", "hill_repression")
                            and actual == float(scalar_value(model, rate["basal"]))
                        ):
                            # The converted induced term is 0 (M2's floor).
                            vmax = float(scalar_value(model, rate["vmax"]))
                            require(
                                error <= vmax * HILL_FLOOR,
                                f"propensity {actual!r} below {expected!r} beyond the Hill floor",
                            )
                            hill_floor.append(name)
                            evaluated_rates += 1
                            continue
                        require(
                            error <= RTOL * abs(expected) + ATOL,
                            f"propensity {actual!r} differs from {expected!r}",
                        )
                        if expected != 0:
                            maximum_relative_error = max(
                                maximum_relative_error, error / abs(expected)
                            )
                        evaluated_rates += 1
                checks.append({"case": name, "passed": True})
            except Exception as error:
                checks.append({"case": name, "passed": False, "error": str(error)})

        rng = random.Random(SEED)
        for index in range(200):
            case(f"random_{index:03d}", random_model(rng, index))
        for model_path in sorted((ROOT / "models/legacy").glob("*.model.json")):
            case(model_path.name, json.loads(model_path.read_text()))
        for n, k, coefficient in [
            (50, 2, 0.5),
            (2000, 1000, 1e-300),
            (2000, 1000, 1),
            (1000, 1000, math.ulp(0.0)),
            (49, 50, 1),
        ]:
            case(
                f"mass_{n}_{k}_{coefficient}",
                make_model(
                    [("A", n)], [reaction("r", {"A": k}, {}, mass_action(coefficient))]
                ),
            )
        for kind, x, k, n, vmax in [
            ("hill_activation", 1, 1e300, 2, 1e300),
            ("hill_repression", 1, math.ulp(0.0), 2, 1),
            ("hill_activation", 1, math.nextafter(1.0, 2.0), 1e15, 1),
            ("hill_repression", 50, 1e-300, 0.01, 7),
            ("hill_activation", 0, 3, 0.5, 5),
        ]:
            case(
                f"{kind}_{x}_{k}",
                make_model(
                    [("X", x)],
                    [
                        reaction(
                            "r",
                            {},
                            {"X": 1},
                            {
                                "kind": kind,
                                "regulator": "X",
                                "vmax": vmax,
                                "basal": 0,
                                "k": k,
                                "n": n,
                            },
                        )
                    ],
                ),
            )
        case(
            "mm_underflow",
            make_model(
                [("S", 1)],
                [
                    reaction(
                        "r",
                        {"S": 1},
                        {},
                        {
                            "kind": "michaelis_menten",
                            "substrate": "S",
                            "vmax": math.ulp(0.0),
                            "km": 1e308,
                        },
                    )
                ],
            ),
        )
        # Deliberately bad independent outputs must be rejected by the checker.
        sample = make_model(
            [("A", 3), ("B", 0)], [reaction("r", {"A": 1}, {"B": 1}, mass_action(1))]
        )
        path.write_text(json.dumps(sample))
        clean = json.loads(
            subprocess.check_output([str(binary), "analyze", str(path)])
        )["analysis"]
        for name, mutate in [
            (
                "bad_stoichiometry",
                lambda data: data["reactions"][0]["net_changes"][0].update(
                    coefficient=1
                ),
            ),
            ("bad_dependency", lambda data: data["dependencies"][0].clear()),
            (
                "bad_conservation",
                lambda data: data["conservation"]["basis"][0][
                    "coefficients"
                ].__setitem__(0, "2"),
            ),
        ]:
            altered = json.loads(json.dumps(clean))
            mutate(altered)
            try:
                compare_structure(sample, altered)
                checks.append(
                    {
                        "case": name,
                        "passed": False,
                        "error": "checker accepted corrupted output",
                    }
                )
            except RuntimeError:
                checks.append({"case": name, "passed": True})
    report = {
        "schema": "crn-craft.network-validation/0.2",
        "scope": "network_analysis_and_propensities_no_simulation",
        "seed": SEED,
        "random_networks": 200,
        "fixture_stream_sha256": fixture_hash.hexdigest(),
        "relative_tolerance": RTOL,
        "absolute_tolerance": ATOL,
        "hill_floor": HILL_FLOOR,
        "evaluated_rates": evaluated_rates,
        "maximum_relative_error": maximum_relative_error,
        "refused_without_reactions": converter_refusals,
        "overflow_failures": overflow_failures,
        "underflow_to_zero": underflow_zero,
        "hill_floor_cases": hill_floor,
        "passed": sum(row["passed"] for row in checks),
        "failed": sum(not row["passed"] for row in checks),
        "checks": checks,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    require(
        report["failed"] == 0,
        "independent network checks failed; inspect retained rows",
    )


if __name__ == "__main__":
    main()
