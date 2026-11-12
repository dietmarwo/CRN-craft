"""FSP reference gate of the 0.7.0 engine merge (M5a).

For every stored 0.6.0 FSP reference (M0c: W1, W4, enzyme-mechanism,
immigration-death, michaelis-menten; M4b: third-order-decay, trimolecular)
this converts the network model with the engine under test (`crn-craft
convert`), solves a `crn-craft.run/1.0` request with method `fsp` through
`crn-craft rules MODEL REQUEST` with the reference's bounds, `max_states`,
tolerance, end and report times, and compares with the reference by species
name and state, as declared in `validation/merge/RESULTS.md` ("M5a ...
declaration"):

1. the projected state set equals the reference's;
2. every probability within 1e-10 absolute;
3. the sink probability at every report time within 1e-10 absolute.

The uniformization rate, the edge count and the largest moment difference
are reported as information. Every failing comparison is kept in the JSON
report, which holds no timings, so reruns compare byte for byte. Standard
library only:

    flock /home/claude/craft5/cpu.lock python validation/merge/check_fsp_references.py \\
        --binary target/release/crn-craft

Exit status: 0 when every check passes, 1 otherwise.
"""

import argparse
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_common import (  # noqa: E402
    HERE,
    ROOT,
    binary_identity,
    encode,
    load_manifest,
    sha256,
)

REFERENCE = HERE / "reference"
HIGHER_ORDER = HERE / "higher-order" / "models.json"
TOLERANCE = 1e-10
OUTPUT = HERE / "fsp" / "references.json"


def cases():
    """(name, network model path, fsp settings, reference file) per model."""
    out = []
    for row in load_manifest():
        if row["fsp"] is not None:
            out.append((row["name"], row["path"], row["fsp"], f"{row['name']}.fsp"))
    for row in load_manifest(path=HIGHER_ORDER):
        out.append(
            (row["name"], row["path"], row["fsp"], f"higher-order/{row['name']}.fsp")
        )
    return out


def solve(binary, path, fsp, reference, folder):
    """Convert the model and solve the 1.0 FSP request; the parsed result."""
    done = subprocess.run(
        [str(binary), "convert", str(ROOT / path)], capture_output=True, check=False
    )
    if done.returncode != 0:
        raise RuntimeError(f"{path}: conversion refused: {done.stdout.decode()}")
    converted = done.stdout
    (folder / "model.json").write_bytes(converted)
    times = [s["time"] for s in reference["snapshots"]]
    request = {
        "schema": "crn-craft.run/1.0",
        "model": {"path": "model.json", "sha256": sha256(converted)},
        "time": {"end": times[-1], "report_times": times},
        "method": {
            "kind": "fsp",
            "bounds": fsp["bounds"],
            "max_states": fsp["max_states"],
            "tolerance": fsp["tolerance"],
        },
        "budget": {"max_output_bytes": 1 << 30},
    }
    (folder / "request.json").write_bytes(encode(request))
    done = subprocess.run(
        [str(binary), "rules", "model.json", "request.json"],
        capture_output=True,
        cwd=folder,
        check=False,
    )
    result = json.loads(done.stdout)
    if done.returncode != 0 or result.get("schema") != "crn-craft.fsp-result/1.0":
        raise RuntimeError(f"{path}: {result}")
    return result


def compare(name, result, reference):
    """The gated checks of one model and its informative differences."""
    species = reference["species"]
    if sorted(result["species"]) != sorted(species):
        return {
            "model": name,
            "passed": False,
            "failures": [{"check": "species", "labels": result["species"]}],
        }
    order = [result["species"].index(s) for s in species]
    states = [tuple(state[i] for i in order) for state in result["states"]]
    index = {state: i for i, state in enumerate(states)}
    reference_states = [tuple(state) for state in reference["states"]]
    same_states = sorted(states) == sorted(reference_states)
    failures = []
    if not same_states:
        failures.append({"check": "state_set", "states": len(states)})
    largest_p = 0.0
    largest_leak = 0.0
    largest_moment = 0.0
    comparisons = 0
    for k, (new, old) in enumerate(
        zip(result["snapshots"], reference["snapshots"], strict=True)
    ):
        if new["time"] != old["time"]:
            failures.append({"check": "time", "snapshot": k})
        for j, state in enumerate(reference_states):
            i = index.get(state)
            p = new["probabilities"][i] if i is not None else 0.0
            difference = abs(p - old["probabilities"][j])
            largest_p = max(largest_p, difference)
            comparisons += 1
            if not difference <= TOLERANCE:
                failures.append(
                    {
                        "check": "probability",
                        "time": old["time"],
                        "state": list(state),
                        "difference": difference,
                    }
                )
        difference = abs(
            new["physical_leakage_estimate"] - old["physical_leakage_estimate"]
        )
        largest_leak = max(largest_leak, difference)
        comparisons += 1
        if not difference <= TOLERANCE:
            failures.append(
                {"check": "leakage", "time": old["time"], "difference": difference}
            )
        for field in ("projected_first_moment", "projected_second_raw_moment"):
            for s, i in enumerate(order):
                largest_moment = max(largest_moment, abs(new[field][i] - old[field][s]))
    return {
        "model": name,
        "species": species,
        "states": len(states),
        "report_times": len(reference["snapshots"]),
        "state_set_equal": same_states,
        "comparisons": comparisons + 1,
        "max_probability_difference": largest_p,
        "max_leakage_difference": largest_leak,
        "information": {
            "uniformization_rate": [
                result["uniformization_rate"],
                reference["uniformization_rate"],
            ],
            "sparse_edges": [result["sparse_edges"], reference["sparse_edges"]],
            "max_raw_moment_difference": largest_moment,
        },
        "failures": failures,
        "passed": not failures,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    binary = args.binary.resolve()
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-fsp-references-") as scratch:
        for name, path, fsp, stem in cases():
            reference = json.loads(
                gzip.decompress((REFERENCE / f"{stem}.json.gz").read_bytes())
            )
            folder = Path(scratch) / name
            folder.mkdir()
            result = solve(binary, path, fsp, reference, folder)
            rows.append(compare(name, result, reference))
    report = {
        "schema": "crn-craft.merge-fsp-references/0.1",
        "declared_in": "validation/merge/RESULTS.md",
        "tolerance": TOLERANCE,
        "binary": binary_identity(binary),
        "rows": rows,
        "checks": sum(r.get("comparisons", 0) for r in rows),
        "failures": sum(len(r["failures"]) for r in rows),
        "passed": all(r["passed"] for r in rows),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "passed": report["passed"],
                "checks": report["checks"],
                "failures": report["failures"],
                "max_probability_difference": max(
                    r.get("max_probability_difference", 0.0) for r in rows
                ),
            }
        )
    )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
