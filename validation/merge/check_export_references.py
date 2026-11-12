"""Export equivalence check of the 0.7.0 engine merge (M5f).

For the twelve models of `models.json` and the fixture of
`validation/check_export.py` this runs `crn-craft export MODEL gillespy2` of
the engine under test and of the 0.6.0 binary, as declared in
`validation/merge/RESULTS.md` ("M5f ... declaration"):

1. a model 0.6.0 exported must export again with a byte-identical
   `artifact` (source, aliases, notes);
2. a model 0.6.0 refused with `CRN-E050` must now export; each custom
   propensity string, evaluated with Python floats (`+ - * /` and `pow`,
   the operators GillesPy2 1.8.3 keeps), must equal the engine's propensity
   (`crn-craft analyze MODEL STATE`) within relative 2e-12 plus absolute
   8 * 2^-1074 at four states: the initial counts, all 0, all 7, and
   initial count + i + 1 for species i;
3. any other 0.6.0 outcome must be repeated with the same code.

Every failing comparison is kept in the JSON report, which holds no timings,
so reruns compare byte for byte. Standard library only:

    flock /home/claude/craft5/cpu.lock python validation/merge/check_export_references.py \\
        --binary target/release/crn-craft --reference-binary BINARY_0_6_0

Exit status: 0 when every check passes, 1 otherwise.
"""

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_common import HERE, ROOT, binary_identity, load_manifest  # noqa: E402

OUTPUT = HERE / "export" / "references.json"
RTOL = 2e-12
ATOL = 8 * math.ulp(0.0)
CUSTOM = re.compile(r"name='R(\d+)', .*propensity_function=(\".*\")\)\)$")


def reaction(name, left, right, rate):
    return {"name": name, "reactants": left, "products": right, "rate": rate}


def mass_action(coefficient):
    return {"kind": "mass_action", "coefficient": coefficient}


def fixture():
    """The model of `validation/check_export.py`, as that script writes it
    (its helpers of `check_network.py`, repeated to stay standard library)."""
    species = [("vol", 3), ("K0", 2), ("S1", 1), ("__import__", 0), ("pow", 1)]
    model = {
        "schema": "crn-craft.model/0.1",
        "name": "oracle_fixture",
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": {"c0": 1.5, "c2": 2.5},
        "species": [{"name": n, "initial_count": c} for n, c in species],
        "reactions": [
            reaction("source", {}, {"vol": 1}, mass_action("c0")),
            reaction("first", {"vol": 1}, {"K0": 1}, mass_action(0.25)),
            reaction("homo", {"vol": 2}, {"S1": 1}, mass_action("c2")),
            reaction("hetero", {"vol": 1, "K0": 1}, {"pow": 1}, mass_action(0.125)),
            reaction(
                "catalyst", {"vol": 1, "K0": 1}, {"vol": 1, "S1": 1}, mass_action(0.5)
            ),
            reaction("zero", {"K0": 1}, {}, mass_action(0)),
        ],
    }
    model["omega"] = 137.0
    return (json.dumps(model, sort_keys=True) + "\n").encode()


def export(binary, path):
    done = subprocess.run(
        [str(binary), "export", str(path), "gillespy2"],
        capture_output=True,
        check=False,
    )
    return done.returncode, json.loads(done.stdout)


def propensities(binary, path, state_path, counts):
    state_path.write_text(
        json.dumps({"schema": "crn-craft.state/0.1", "counts": counts})
    )
    done = subprocess.run(
        [str(binary), "analyze", str(path), str(state_path)],
        capture_output=True,
        check=False,
    )
    report = json.loads(done.stdout)
    if done.returncode != 0:
        raise RuntimeError(f"analyze failed: {report['diagnostics']}")
    return report["analysis"]["propensities"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--reference-binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    binary, reference = args.binary.resolve(), args.reference_binary.resolve()
    rows, artifacts, compared, failures = [], 0, 0, 0
    with tempfile.TemporaryDirectory(prefix="crn-export-gate-") as temp:
        folder = Path(temp)
        cases = [(row["name"], ROOT / row["path"]) for row in load_manifest()]
        (folder / "fixture.model.json").write_bytes(fixture())
        cases.append(("check_export_fixture", folder / "fixture.model.json"))
        for name, path in cases:
            model = json.loads(path.read_text())
            old_code, old = export(reference, path)
            new_code, new = export(binary, path)
            row = {"model": name, "failures": []}
            fail = row["failures"].append
            if old["success"]:
                row["outcome_0_6_0"] = "exported"
                if new_code != 0 or not new["success"]:
                    fail({"check": "exports", "diagnostics": new["diagnostics"]})
                elif new["artifact"] != old["artifact"]:
                    fail({"check": "artifact_identical"})
                artifacts += 1
            elif any(d["code"] == "CRN-E050" for d in old["diagnostics"]):
                row["outcome_0_6_0"] = "CRN-E050"
                if new_code != 0 or not new["success"]:
                    fail({"check": "exports", "diagnostics": new["diagnostics"]})
                else:
                    aliases = new["artifact"]["species_aliases"]
                    custom = {}
                    for line in new["artifact"]["source"].splitlines():
                        match = CUSTOM.search(line)
                        if match:
                            custom[int(match.group(1))] = json.loads(match.group(2))
                    row["custom_propensities"] = len(custom)
                    names = [s["name"] for s in model["species"]]
                    initial = {s["name"]: s["initial_count"] for s in model["species"]}
                    states = [
                        initial,
                        dict.fromkeys(names, 0),
                        dict.fromkeys(names, 7),
                        {n: initial[n] + i + 1 for i, n in enumerate(names)},
                    ]
                    largest = 0.0
                    for counts in states:
                        engine = propensities(
                            binary, path, folder / "state.json", counts
                        )
                        values = {aliases[n]: float(c) for n, c in counts.items()}
                        for channel, text in sorted(custom.items()):
                            actual = eval(  # noqa: S307
                                text, {"__builtins__": {}, "pow": pow}, values
                            )
                            expected = engine[channel]
                            error = abs(actual - expected)
                            compared += 1
                            if not error <= RTOL * abs(expected) + ATOL:
                                fail(
                                    {
                                        "check": "custom_propensity",
                                        "reaction": channel,
                                        "counts": counts,
                                        "python": actual,
                                        "engine": expected,
                                    }
                                )
                            if expected:
                                largest = max(largest, error / abs(expected))
                    row["largest_relative_difference"] = largest
            else:
                codes = sorted({d["code"] for d in old["diagnostics"]})
                row["outcome_0_6_0"] = codes
                new_codes = sorted({d["code"] for d in new["diagnostics"]})
                if new["success"] or new_codes != codes:
                    fail({"check": "same_refusal", "codes": new_codes})
            failures += len(row["failures"])
            rows.append(row)
            print(
                f"{name}: {row['outcome_0_6_0']}, {len(row['failures'])} failures",
                file=sys.stderr,
            )
    report = {
        "schema": "crn-craft.merge-export-references/0.1",
        "declared_in": "validation/merge/RESULTS.md (M5f declaration)",
        "binary": binary_identity(binary),
        "reference_binary": binary_identity(reference),
        "relative_tolerance": RTOL,
        "absolute_tolerance": ATOL,
        "identical_artifacts_compared": artifacts,
        "custom_propensities_compared": compared,
        "failures": failures,
        "passed": failures == 0,
        "models": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {"artifacts": artifacts, "propensities": compared, "failures": failures}
        )
    )
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
