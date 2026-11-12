"""Build actual GillesPy2 exports; independently verify count conventions and smoke SSA.

Since 0.7.0 (M5f) `crn-craft export` reads network models through the
converter, and reactions without a constant rate of order at most two
export as custom propensity strings. The Michaelis-Menten and repressilator
examples, refused by 0.6.0, are exported and checked like the others: their
custom propensities (evaluated with `pow`, which GillesPy2 writes for `**`)
against the exact rate law and the engine's rates, and a reference SSA smoke
run of the Michaelis-Menten model. The refusal checks use a rate with a
function and a stop condition, which GillesPy2 cannot evaluate.
"""

import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import subprocess
import tempfile

import gillespy2 as gp
import numpy as np
from jsonschema import Draft7Validator
from check_network import (
    ATOL,
    RTOL,
    exact_rate,
    make_model,
    mass_action,
    reaction,
    require,
)

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(gp.__version__ == "1.8.3", "requires pinned GillesPy2 1.8.3")
    binary = str(args.binary.resolve())
    validator = Draft7Validator(
        json.loads(subprocess.check_output([binary, "describe", "export"]))
    )
    checks, rate_checks, max_error = [], 0, 0.0
    source_hashes = {}
    fixture = make_model(
        [("vol", 3), ("K0", 2), ("S1", 1), ("__import__", 0), ("pow", 1)],
        [
            reaction("source", {}, {"vol": 1}, mass_action("c0")),
            reaction("first", {"vol": 1}, {"K0": 1}, mass_action(0.25)),
            reaction("homo", {"vol": 2}, {"S1": 1}, mass_action("c2")),
            reaction("hetero", {"vol": 1, "K0": 1}, {"pow": 1}, mass_action(0.125)),
            reaction(
                "catalyst", {"vol": 1, "K0": 1}, {"vol": 1, "S1": 1}, mass_action(0.5)
            ),
            reaction("zero", {"K0": 1}, {}, mass_action(0)),
        ],
        {"c0": 1.5, "c2": 2.5},
    )
    fixture["omega"] = 137.0
    models = [("grid_all_elementary", fixture)] + [
        (
            name,
            json.loads((ROOT / "models/legacy" / (name + ".model.json")).read_text()),
        )
        for name in (
            "immigration-death",
            "enzyme-mechanism",
            "michaelis-menten",
            "repressilator",
        )
    ]
    with tempfile.TemporaryDirectory(prefix="crn-export-") as folder:
        path, state_path = Path(folder) / "model.json", Path(folder) / "state.json"
        for label, model in models:
            raw = (json.dumps(model, sort_keys=True) + "\n").encode()
            path.write_bytes(raw)
            try:
                exported = subprocess.run(
                    [binary, "export", str(path), "gillespy2"],
                    capture_output=True,
                    check=False,
                    timeout=15,
                )
                require(
                    exported.returncode == 0 and not exported.stderr,
                    "export CLI failed",
                )
                report = json.loads(exported.stdout)
                validator.validate(report)
                require(
                    report["success"]
                    and report["model_sha256"] == hashlib.sha256(raw).hexdigest(),
                    "export identity failed",
                )
                artifact = report["artifact"]
                source_hashes[label] = hashlib.sha256(
                    artifact["source"].encode()
                ).hexdigest()
                namespace = {}
                exec(
                    compile(artifact["source"], "generated_reference.py", "exec"),
                    namespace,
                )
                target = namespace["build_model"]()
                target.compile_prep()
                require(
                    target.volume == 1.0,
                    "target volume rescales stochastic coefficients",
                )
                names = [s["name"] for s in model["species"]]
                aliases = artifact["species_aliases"]
                for source_species in model["species"]:
                    require(
                        target.listOfSpecies[
                            aliases[source_species["name"]]
                        ].initial_value
                        == source_species["initial_count"],
                        "initial population changed",
                    )
                target_reactions = [
                    target.listOfReactions[artifact["reaction_aliases"][r["name"]]]
                    for r in model["reactions"]
                ]
                for original, translated in zip(
                    model["reactions"], target_reactions, strict=True
                ):
                    for side in ("reactants", "products"):
                        actual = {
                            s.name: int(n) for s, n in getattr(translated, side).items()
                        }
                        require(
                            actual
                            == {aliases[k]: v for k, v in original[side].items()},
                            "export stoichiometry changed",
                        )
                if label.startswith("grid"):
                    states = [
                        {**dict.fromkeys(names, 1), "vol": a, "K0": b}
                        for a in range(51)
                        for b in range(51)
                    ]
                elif label == "immigration-death":
                    states = [{"A": a} for a in range(51)]
                else:
                    states = [
                        {s["name"]: s["initial_count"] for s in model["species"]},
                        dict.fromkeys(names, 0),
                        dict.fromkeys(names, 7),
                    ]
                for counts in states:
                    state_path.write_text(
                        json.dumps({"schema": "crn-craft.state/0.1", "counts": counts})
                    )
                    analyzed = json.loads(
                        subprocess.check_output(
                            [binary, "analyze", str(path), str(state_path)], timeout=15
                        )
                    )
                    require(analyzed["success"], "Rust rate evaluation failed")
                    values = {aliases[name]: count for name, count in counts.items()}
                    values.update(
                        {name: p.value for name, p in target.listOfParameters.items()}
                    )
                    values["vol"] = target.volume
                    for original, translated, rust_rate in zip(
                        model["reactions"],
                        target_reactions,
                        analyzed["analysis"]["propensities"],
                        strict=True,
                    ):
                        # Custom propensities use pow, as GillesPy2 writes `**`.
                        actual = eval(
                            translated.propensity_function,
                            {"__builtins__": {}, "pow": pow},
                            values,
                        )
                        expected = float(exact_rate(model, original, counts))
                        require(
                            math.isfinite(actual) and actual >= 0,
                            "invalid target propensity",
                        )
                        if label.startswith("grid"):
                            require(
                                actual == rust_rate == expected,
                                f"exact elementary rate mismatch at {counts}: {actual}, {rust_rate}, {expected}",
                            )
                        else:
                            require(
                                abs(actual - expected) <= RTOL * abs(expected) + ATOL,
                                "target rate mismatch",
                            )
                            require(
                                abs(rust_rate - expected)
                                <= RTOL * abs(expected) + ATOL,
                                "Rust reference mismatch",
                            )
                        if expected:
                            max_error = max(
                                max_error, abs(actual - expected) / abs(expected)
                            )
                        rate_checks += 1
                checks.append({"case": label, "passed": True, "states": len(states)})
                if label == "michaelis-menten":
                    # The custom propensity runs in the reference solver.
                    target.timespan(np.array([0.0, 1.0, 2.0]))
                    result = target.run(
                        solver=gp.NumPySSASolver(model=target),
                        number_of_trajectories=2,
                        seed=20261111,
                    )
                    require(len(result) == 2, "wrong custom smoke replicate count")
                    for trajectory in result:
                        total = trajectory["S0"] + trajectory["S1"]
                        require(
                            all(float(v) == 50 for v in total)
                            and all(float(v) >= 0 for v in trajectory["S0"]),
                            "custom propensity smoke run lost molecules",
                        )
                    checks.append(
                        {
                            "case": "reference_ssa_smoke_custom_propensity",
                            "passed": True,
                        }
                    )
                if label == "immigration-death":
                    target.timespan(np.array([0.0, 1.0, 2.0]))
                    result = target.run(
                        solver=gp.NumPySSASolver(model=target),
                        number_of_trajectories=2,
                        seed=20261103,
                    )
                    require(len(result) == 2, "wrong smoke replicate count")
                    payload = []
                    for trajectory in result:
                        require(
                            list(trajectory["time"]) == [0, 1, 2],
                            "wrong reference time grid",
                        )
                        counts = [float(x) for x in trajectory["S0"]]
                        require(
                            all(
                                math.isfinite(v) and v >= 0 and v.is_integer()
                                for v in counts
                            ),
                            "invalid reference populations",
                        )
                        payload.append(counts)
                    checks.append(
                        {
                            "case": "reference_ssa_smoke_only",
                            "passed": True,
                            "trajectories": payload,
                        }
                    )
            except Exception as error:
                checks.append({"case": label, "passed": False, "error": str(error)})
        # What GillesPy2 cannot evaluate is refused (0.6.0 refused these two
        # examples, which now export custom propensities).
        rules = json.loads(
            (ROOT / "models/rules/immigration-death.rules.json").read_text()
        )
        function = json.loads(json.dumps(rules))
        # A rate that reads the observable X through a function.
        function["rules"][-1]["rate"] = (
            function["rules"][-1]["rate"] + " * exp(X / 1000)"
        )
        stopped = json.loads(json.dumps(rules))
        stopped["stop"] = "1 > 2"
        for name, document, path in [
            ("function_rate", function, f"/rules/{len(rules['rules']) - 1}/rate"),
            ("stop_condition", stopped, "/stop"),
        ]:
            raw = (json.dumps(document, sort_keys=True) + "\n").encode()
            path_file = Path(folder) / "unsupported.json"
            path_file.write_bytes(raw)
            result = subprocess.run(
                [binary, "export", str(path_file), "gillespy2"],
                capture_output=True,
                check=False,
                timeout=15,
            )
            report = json.loads(result.stdout)
            validator.validate(report)
            checks.append(
                {
                    "case": "unsupported_" + name,
                    "passed": result.returncode == 1
                    and not report["success"]
                    and report["artifact"] is None
                    and any(
                        d["code"] == "CRN-E050" and d["path"] == path
                        for d in report["diagnostics"]
                    ),
                }
            )
    report = {
        "schema": "crn-craft.reference-export-validation/0.2",
        "gillespy2": gp.__version__,
        "numpy": importlib.metadata.version("numpy"),
        "scope": "exported_propensities_and_reference_smoke_no_crn_simulation",
        "source_sha256": source_hashes,
        "rate_comparisons": rate_checks,
        "maximum_relative_error": max_error,
        "relative_tolerance": RTOL,
        "absolute_tolerance": ATOL,
        "passed": sum(c["passed"] for c in checks),
        "failed": sum(not c["passed"] for c in checks),
        "checks": checks,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    require(
        report["failed"] == 0, "reference export checks failed; inspect retained rows"
    )


if __name__ == "__main__":
    main()
