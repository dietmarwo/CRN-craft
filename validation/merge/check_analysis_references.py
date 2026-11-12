"""Analysis reference gate of the 0.7.0 engine merge (M5b).

For each of the twelve models of `models.json` this runs `crn-craft analyze`
of the engine under test twice, on the network model file (read through the
converter) and on the model converted by `crn-craft convert`, and compares
each `crn-craft.analysis/1.0` report with the stored 0.6.0 report
(`reference/<name>.analysis.json.gz`, `crn-craft.analysis/0.2`) leaf by leaf,
as declared in `validation/merge/RESULTS.md` ("M5b ... declaration"):

- every object key, list length and value must be equal;
- `schema` is ignored (checked to be `crn-craft.analysis/1.0`), and in the
  run on the converted model `model_sha256` (checked to be the hash of the
  converted bytes);
- an element of `propensities` may differ for a Hill law (relative
  difference <= 1e-14 * max(1, n)) or mass action over two molecules
  (relative difference <= 4 eps); each such difference is listed.

Every failing comparison is kept in the JSON report, which holds no timings,
so reruns compare byte for byte. Standard library only:

    flock /home/claude/craft5/cpu.lock python validation/merge/check_analysis_references.py \\
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
from merge_common import HERE, ROOT, binary_identity, load_manifest, sha256  # noqa: E402

REFERENCE = HERE / "reference"
OUTPUT = HERE / "analysis" / "references.json"
EPSILON = 2.0**-52
SCHEMA = "crn-craft.analysis/1.0"


def tolerance(model, channel):
    """The declared relative tolerance of a propensity, or None (exact)."""
    reaction = model["reactions"][channel]
    rate = reaction["rate"]
    if rate["kind"] in ("hill_activation", "hill_repression"):
        n = rate["n"]
        if isinstance(n, str):
            n = model["parameters"][n]
        return 1e-14 * max(1.0, float(n))
    if rate["kind"] == "mass_action" and sum(reaction["reactants"].values()) == 2:
        return 4 * EPSILON
    return None


def compare(reference, actual, path, model, out):
    """Appends (path, reference, actual, note) for every differing leaf;
    returns the number of leaves compared."""
    if isinstance(reference, dict) and isinstance(actual, dict):
        count = 0
        for key in sorted(set(reference) | set(actual)):
            where = f"{path}/{key}"
            if key not in reference or key not in actual:
                out.append((where, reference.get(key), actual.get(key), "missing key"))
                count += 1
                continue
            count += compare(reference[key], actual[key], where, model, out)
        return count
    if isinstance(reference, list) and isinstance(actual, list):
        if len(reference) != len(actual):
            out.append((path, len(reference), len(actual), "list length"))
            return 1
        return sum(
            compare(r, a, f"{path}/{i}", model, out)
            for i, (r, a) in enumerate(zip(reference, actual, strict=True))
        )
    if type(reference) is not type(actual) or reference != actual:
        out.append((path, reference, actual, "value"))
    return 1


def check(name, model, reference, report, code, converted_hash):
    """The rows of one run: failures and listed propensity differences."""
    rows, listed = [], []
    if code != 0:
        rows.append({"path": "exit_code", "reference": 0, "actual": code})
    if report.get("schema") != SCHEMA:
        rows.append(
            {"path": "/schema", "reference": SCHEMA, "actual": report.get("schema")}
        )
    if converted_hash is not None and report.get("model_sha256") != converted_hash:
        rows.append(
            {
                "path": "/model_sha256",
                "reference": converted_hash,
                "actual": report.get("model_sha256"),
            }
        )
    ignored = {"schema"} | ({"model_sha256"} if converted_hash is not None else set())
    left = {k: v for k, v in reference.items() if k not in ignored}
    right = {k: v for k, v in report.items() if k not in ignored}
    differences = []
    leaves = compare(left, right, "", model, differences)
    for path, old, new, note in differences:
        parts = path.split("/")
        if (
            note == "value"
            and len(parts) == 4
            and parts[1:3] == ["analysis", "propensities"]
            and isinstance(old, float)
            and isinstance(new, float)
        ):
            channel = int(parts[3])
            bound = tolerance(model, channel)
            relative = abs(new - old) / abs(old) if old else float("inf")
            row = {
                "path": path,
                "reference": old,
                "actual": new,
                "relative": relative,
                "tolerance": bound,
            }
            if bound is not None and relative <= bound:
                listed.append(row)
                continue
            rows.append(row)
            continue
        rows.append({"path": path, "reference": old, "actual": new, "note": note})
    return leaves, rows, listed


def run(binary, path):
    done = subprocess.run(
        [str(binary), "analyze", str(path)], capture_output=True, check=False
    )
    if done.stderr:
        raise RuntimeError(f"{path}: stderr {done.stderr.decode()!r}")
    return done.returncode, json.loads(done.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/crn-craft"
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    binary = args.binary.resolve()
    results, total, failures = [], 0, 0
    with tempfile.TemporaryDirectory(prefix="crn-analysis-gate-") as temp:
        folder = Path(temp)
        for row in load_manifest():
            name = row["name"]
            network = ROOT / row["path"]
            model = json.loads(network.read_text())
            reference = json.loads(
                gzip.decompress((REFERENCE / f"{name}.analysis.json.gz").read_bytes())
            )
            done = subprocess.run(
                [str(binary), "convert", str(network)], capture_output=True, check=False
            )
            if done.returncode != 0:
                raise RuntimeError(
                    f"{name}: conversion refused: {done.stdout.decode()}"
                )
            converted = folder / f"{name}.model.json"
            converted.write_bytes(done.stdout)
            for label, path, converted_hash in [
                ("network", network, None),
                ("converted", converted, sha256(done.stdout)),
            ]:
                code, report = run(binary, path)
                leaves, rows, listed = check(
                    name, model, reference, report, code, converted_hash
                )
                total += leaves
                failures += len(rows)
                largest = max((r["relative"] for r in listed), default=0.0)
                results.append(
                    {
                        "model": name,
                        "input": label,
                        "leaves": leaves,
                        "failures": rows,
                        "propensity_differences": listed,
                        "largest_relative_propensity_difference": largest,
                    }
                )
                print(
                    f"{name} ({label}): {leaves} leaves, {len(rows)} failures, "
                    f"{len(listed)} listed propensity differences",
                    file=sys.stderr,
                )
    report = {
        "schema": "crn-craft.merge-analysis-references/0.1",
        "declared_in": "validation/merge/RESULTS.md (M5b declaration)",
        "binary": binary_identity(binary),
        "leaves": total,
        "failures": failures,
        "passed": failures == 0,
        "models": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"leaves": total, "failures": failures}))
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
