"""Exercise four source mutations in an isolated copy; never edits the checkout."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MUTATIONS = [
    (
        "missing_binomial_divisor",
        "crates/crn-craft/src/rates.rs",
        "(count - order + i) as f64 / i as f64",
        "(count - order + i) as f64 / 1.0",
    ),
    (
        "wrong_product_sign",
        "crates/crn-craft/src/network.rs",
        "*changes.entry(term.species).or_insert(0) += term.coefficient;",
        "*changes.entry(term.species).or_insert(0) -= term.coefficient;",
    ),
    (
        "missing_custom_rate_read",
        "crates/crn-craft/src/network.rs",
        "reads.insert(species);",
        "let _ = species;",
    ),
    (
        "missing_reactant_gate",
        "crates/crn-craft/src/rates.rs",
        "state[term.species] < term.coefficient",
        "state[term.species] < 0",
    ),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cargo", default="cargo")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-mutations-") as folder:
        clone = Path(folder) / "project"
        shutil.copytree(
            ROOT,
            clone,
            ignore=shutil.ignore_patterns(
                ".git", "target", "__pycache__", ".ruff_cache", "iteration-02"
            ),
        )
        for name, relative, old, new in MUTATIONS:
            path = clone / relative
            original = path.read_text()
            if original.count(old) != 1:
                raise RuntimeError(f"{name}: mutation anchor not unique")
            path.write_text(original.replace(old, new))
            try:
                result = subprocess.run(
                    [
                        args.cargo,
                        "test",
                        "--release",
                        "--locked",
                        "--offline",
                        "-p",
                        "crn-craft",
                        "--test",
                        "network",
                    ],
                    cwd=clone,
                    capture_output=True,
                    text=True,
                    timeout=120,
                    check=False,
                )
                combined = result.stdout + result.stderr
                failed_tests = [
                    line
                    for line in result.stdout.splitlines()
                    if line.startswith("test ") and line.endswith(" ... FAILED")
                ]
                caught = (
                    result.returncode == 101
                    and bool(failed_tests)
                    and "could not compile" not in combined
                )
                rows.append(
                    {
                        "mutation": name,
                        "file": relative,
                        "caught_by_test_failure": caught,
                        "exit_code": result.returncode,
                        "failed_tests": failed_tests,
                    }
                )
            finally:
                path.write_text(original)
    report = {
        "schema": "crn-craft.mutation-check/0.1",
        "mutations": rows,
        "caught": sum(row["caught_by_test_failure"] for row in rows),
        "total": len(MUTATIONS),
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    if report["caught"] != report["total"]:
        raise RuntimeError("a mutation survived or did not produce a test failure")


if __name__ == "__main__":
    main()
