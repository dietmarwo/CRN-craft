"""Frozen connected workloads for the declared five-step SSA campaign."""

import hashlib
import json
from pathlib import Path
from benchmark_campaign import models

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "validation/ssa-acceleration/models"


def encode(model):
    return (json.dumps(model) + "\n").encode()


def base(name, species, reactions):
    return dict(
        schema="crn-craft.model/0.1",
        name=name,
        convention="stochastic_counts",
        omega=1.0,
        parameters={},
        species=[dict(name=s, initial_count=x) for s, x in species],
        reactions=reactions,
    )


def rx(name, a, b, k):
    return dict(
        name=name, reactants=a, products=b, rate=dict(kind="mass_action", coefficient=k)
    )


def workloads():
    for i, (_, model) in enumerate(models(), 1):
        yield f"W{i}", model, 100.0, 64
    species = [(f"X{i}", 4) for i in range(256)]
    reactions = []
    for i in range(255):
        reactions.extend(
            [
                rx(f"right{i}", {f"X{i}": 1}, {f"X{i + 1}": 1}, 0.5),
                rx(f"left{i}", {f"X{i + 1}": 1}, {f"X{i}": 1}, 0.5),
            ]
        )
    yield "W5", base("DiffusionChain256", species, reactions), 20.0, 16
    species = [(f"X{i}", 32) for i in range(32)]
    reactions = [
        rx(
            f"convert{i}_{j}",
            {f"X{i}": 1},
            {f"X{j}": 1},
            0.005 * (1 + (7 * i + 11 * j) % 8),
        )
        for i in range(32)
        for j in range(32)
        if i != j
    ]
    yield "W6", base("DenseConversion32", species, reactions), 20.0, 16
    species = [
        (f"{s}{i}", x) for i in range(32) for s, x in [("A", 24), ("B", 12), ("C", 4)]
    ]
    reactions = []
    for i in range(32):
        reactions.extend(
            [
                rx(f"bind{i}", {f"A{i}": 1, f"B{i}": 1}, {f"C{i}": 1}, 0.015),
                rx(f"unbind{i}", {f"C{i}": 1}, {f"A{i}": 1, f"B{i}": 1}, 0.5),
            ]
        )
    for i in range(31):
        reactions.extend(
            [
                rx(f"right{i}", {f"A{i}": 1}, {f"A{i + 1}": 1}, 0.1),
                rx(f"left{i}", {f"A{i + 1}": 1}, {f"A{i}": 1}, 0.1),
            ]
        )
    yield "W7", base("CompartmentBinding32", species, reactions), 20.0, 16


def freeze():
    FOLDER.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, model, end, paths in workloads():
        data = encode(model)
        target = FOLDER / f"{name}.model.json"
        if target.exists() and target.read_bytes() != data:
            raise RuntimeError("frozen model changed")
        target.write_bytes(data)
        rows.append(
            dict(
                name=name,
                path=target.name,
                sha256=hashlib.sha256(data).hexdigest(),
                species=len(model["species"]),
                reactions=len(model["reactions"]),
                end=end,
                timing_paths=paths,
                validation_paths=512,
            )
        )
    manifest = dict(
        schema="crn-craft.acceleration-workloads/0.1",
        validation_engines=17,
        validation_comparisons=6946,
        family_alpha=0.01,
        rows=rows,
    )
    (FOLDER / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    print(json.dumps(freeze(), indent=2))
