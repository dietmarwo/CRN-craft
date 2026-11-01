"""Independent support-enumeration cone oracle and reachability graph checker."""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import tempfile

from jsonschema import Draft7Validator
from check_network import make_model, mass_action, reaction, random_model, require, rref

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261103


def primitive(vector):
    common = math.lcm(*(v.denominator for v in vector))
    values = [int(v * common) for v in vector]
    divisor = math.gcd(*values)
    return [v // divisor for v in values]


def support_rays(columns, n):
    """All supports; a positive 1D nullspace is exactly an extreme ray."""
    rays = []
    for size in range(1, n + 1):
        for support in itertools.combinations(range(n), size):
            restricted = [[column[i] for i in support] for column in columns]
            _, basis = rref(restricted, size)
            if len(basis) != 1:
                continue
            row = basis[0]
            if all(v < 0 for v in row):
                row = [-v for v in row]
            if not all(v > 0 for v in row):
                continue
            small = primitive(row)
            full = [0] * n
            for i, value in zip(support, small, strict=True):
                full[i] = value
            rays.append(full)
    return sorted(rays)


def graph_reference(model, names, rank):
    complexes, edges = [], []
    for r in model["reactions"]:
        edge = []
        for side in ("reactants", "products"):
            complex_ = [r[side].get(name, 0) for name in names]
            if complex_ not in complexes:
                complexes.append(complex_)
            edge.append(complexes.index(complex_))
        edges.append(edge)
    n = len(complexes)
    reach = [{i} for i in range(n)]
    weak = [{i} for i in range(n)]
    for a, b in edges:
        reach[a].add(b)
        weak[a].add(b)
        weak[b].add(a)
    # Set-valued transitive closure, independent of the Rust DFS/Kosaraju walk.
    for graph in (reach, weak):
        for k in range(n):
            for i in range(n):
                if k in graph[i]:
                    graph[i] |= graph[k]
    classes = sorted({tuple(sorted(c)) for c in weak})
    strong = sorted(
        {
            tuple(j for j in range(n) if j in reach[i] and i in reach[j])
            for i in range(n)
        }
    )
    terminal = [c for c in strong if not any(a in c and b not in c for a, b in edges)]
    return {
        "complexes": [
            [{"species": i, "coefficient": c} for i, c in enumerate(col) if c]
            for col in complexes
        ],
        "reaction_edges": [{"source": a, "target": b} for a, b in edges],
        "linkage_classes": [list(c) for c in classes],
        "strong_linkage_classes": [list(c) for c in strong],
        "terminal_strong_linkage_classes": [list(c) for c in terminal],
        "weakly_reversible": classes == strong,
        "deficiency": n - len(classes) - rank,
    }


def verify(model, data):
    names = [s["name"] for s in model["species"]]
    columns = [
        [r["products"].get(n, 0) - r["reactants"].get(n, 0) for n in names]
        for r in model["reactions"]
    ]
    rays = support_rays(columns, len(names))
    cone = data["nonnegative_conservation"]
    require(cone["status"] == "exact", "small cone unavailable")
    actual = [[int(c) for c in law["coefficients"]] for law in cone["extreme_rays"]]
    require(actual == rays, "complete primitive extreme ray set differs")
    bounds = [None] * len(names)
    for row, law in zip(rays, cone["extreme_rays"], strict=True):
        total = sum(
            v * s["initial_count"] for v, s in zip(row, model["species"], strict=True)
        )
        require(law["initial_total"] == str(total), "wrong exact total")
        for col in columns:
            require(
                sum(a * b for a, b in zip(row, col, strict=True)) == 0,
                "nonzero residual",
            )
        for i, weight in enumerate(row):
            if weight:
                candidate = total // weight
                bounds[i] = (
                    candidate if bounds[i] is None else min(bounds[i], candidate)
                )
    bounds = [str(v) if v is not None else None for v in bounds]
    require(cone["initial_class_upper_bounds"] == bounds, "wrong cone bounds")
    require(
        data["initial_class_upper_bounds"] == bounds,
        "aggregate bounds did not select exact cone",
    )
    rank, _ = rref(columns, len(names))
    require(
        data["linkage"] == graph_reference(model, names, rank),
        "complex/linkage graph mismatch",
    )
    old = data["conservation"]["initial_class_upper_bounds"]
    for previous, current in zip(old, bounds, strict=True):
        if previous is not None:
            require(
                current is not None and int(current) <= int(previous),
                "cone lost or weakened a sound bound",
            )
    return (
        sum(v is not None for v in old),
        sum(v is not None for v in bounds),
        len(rays),
    )


def special_models():
    review = json.loads((ROOT / "models/analysis/cone-gap.model.json").read_text())
    adjacency = make_model(
        [(n, i + 1) for i, n in enumerate("ABCD")],
        [
            reaction("r0", {"C": 1, "D": 1}, {"A": 1, "B": 1}, mass_action(1)),
            reaction("r1", {"B": 1, "D": 1}, {"A": 1, "C": 1}, mass_action(1)),
        ],
    )
    weighted = make_model(
        [("A", 3), ("B", 2)], [reaction("r", {"A": 2}, {"B": 1}, mass_action(1))]
    )
    enzyme = json.loads(
        (ROOT / "models/examples/enzyme-mechanism.model.json").read_text()
    )
    return [
        ("review_counterexample", review),
        ("adjacency", adjacency),
        ("weighted", weighted),
        ("enzyme", enzyme),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    schema = json.loads(subprocess.check_output([binary, "describe", "analysis"]))
    validator = Draft7Validator(schema)
    checks, old_bounds, new_bounds, total_rays = [], 0, 0, 0
    fixtures_hash = hashlib.sha256()
    rng = random.Random(SEED)
    fixtures = [
        (f"random_{i:03}", random_model(rng, i)) for i in range(300)
    ] + special_models()
    sample = None
    with tempfile.TemporaryDirectory(prefix="crn-structure-") as folder:
        path = Path(folder) / "model.json"
        for label, model in fixtures:
            raw = (json.dumps(model, sort_keys=True) + "\n").encode()
            fixtures_hash.update(raw)
            path.write_bytes(raw)
            try:
                result = subprocess.run(
                    [binary, "analyze", str(path)],
                    capture_output=True,
                    timeout=15,
                    check=False,
                )
                require(result.returncode == 0 and not result.stderr, "analyze failed")
                report = json.loads(result.stdout)
                validator.validate(report)
                require(
                    report["schema"] == "crn-craft.analysis/0.2" and report["success"],
                    "wrong analysis contract",
                )
                before, after, rays = verify(model, report["analysis"])
                old_bounds += before
                new_bounds += after
                total_rays += rays
                checks.append(
                    {
                        "case": label,
                        "passed": True,
                        "basis_bounds": before,
                        "cone_bounds": after,
                        "rays": rays,
                    }
                )
                if label == "review_counterexample":
                    sample = (model, report["analysis"])
                    require(
                        before == 0 and after == 4, "review regression not resolved"
                    )
            except Exception as error:
                checks.append({"case": label, "passed": False, "error": str(error)})
        require(sample is not None, "missing mutation sample")
        for label, mutate in [
            (
                "missing_ray",
                lambda d: d["nonnegative_conservation"]["extreme_rays"].pop(),
            ),
            (
                "wrong_bound",
                lambda d: d["nonnegative_conservation"][
                    "initial_class_upper_bounds"
                ].__setitem__(0, "999"),
            ),
            (
                "reversed_edge",
                lambda d: d["linkage"]["reaction_edges"][0].update(source=1, target=0),
            ),
            ("wrong_scc", lambda d: d["linkage"]["strong_linkage_classes"].clear()),
        ]:
            altered = json.loads(json.dumps(sample[1]))
            mutate(altered)
            try:
                verify(sample[0], altered)
                checks.append(
                    {
                        "case": label,
                        "passed": False,
                        "error": "checker accepted corrupted output",
                    }
                )
            except RuntimeError:
                checks.append({"case": label, "passed": True})
    report = {
        "schema": "crn-craft.structure-validation/0.1",
        "seed": SEED,
        "random_networks": 300,
        "fixture_stream_sha256": fixtures_hash.hexdigest(),
        "basis_bounded_species": old_bounds,
        "cone_bounded_species": new_bounds,
        "additional_bounded_species": new_bounds - old_bounds,
        "extreme_rays_checked": total_rays,
        "passed": sum(c["passed"] for c in checks),
        "failed": sum(not c["passed"] for c in checks),
        "checks": checks,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    require(
        report["failed"] == 0, "structural validation failed; inspect retained rows"
    )


if __name__ == "__main__":
    main()
