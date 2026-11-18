#!/usr/bin/env python3
"""Check the script and the tables of guide chapter 16 (Kowald & Kirkwood 2013).

    python validation/check_kk13.py --output kk13.json
    python -O validation/check_kk13.py --output kk13-optimized.json

1. Runs `examples/kk13/kowald_kirkwood_2013.py --size quick` into a fresh
   directory and compares every result CSV, and the cell counts, seeds and
   event counts of its JSON files, with the stored results in
   `examples/kk13/results/quick/`, exactly.
2. Checks that the 0.8 tables of the chapter's pages equal what
   `examples/kk13/make_tables.py` makes from `examples/kk13/results/draft/`.

The report has no timings, so a normal and an optimized run write equal
reports. Requires the installed wheel and numpy.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "examples" / "kk13" / "kowald_kirkwood_2013.py"
TABLES = ROOT / "examples" / "kk13" / "make_tables.py"
STORED = ROOT / "examples" / "kk13" / "results" / "quick"
PAGES = sorted((ROOT / "docs" / "guide").glob("16*.md"))
# Simulation outputs compare exactly; the logistic fit's estimates (numpy
# linear algebra, whose last bits depend on the CPU) to 6 significant digits.
JSON_FIELDS = ("cells", "seeds", "events", "cox_negative")


def _rounded(value):
    if isinstance(value, float):
        return float(f"{value:.6g}")
    if isinstance(value, dict):
        return {k: _rounded(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_rounded(v) for v in value]
    return value


def _python(*args):
    flags = ["-O"] if sys.flags.optimize else []
    return subprocess.run(
        [sys.executable, *flags, *map(str, args)],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


def check_quick(report):
    with tempfile.TemporaryDirectory() as tmp:
        run = _python(SCRIPT, "--size", "quick", "--threads", "2", "--results", tmp)
        if run.returncode != 0:
            report["failures"].append(f"quick run failed: {run.stderr[-2000:]}")
            return
        fresh = Path(tmp)
        stored = sorted(p.name for p in STORED.glob("*.csv"))
        made = sorted(p.name for p in fresh.glob("*.csv"))
        if stored != made:
            report["failures"].append(
                f"result files differ: stored {stored}, made {made}"
            )
        for name in stored:
            if name not in made:
                continue
            differs = []
            if (STORED / name).read_bytes() != (fresh / name).read_bytes():
                differs.append("the CSV")
            info = STORED / name.replace(".csv", ".json")
            if info.exists():
                a = json.loads(info.read_text())
                b = json.loads((fresh / info.name).read_text())
                differs += [k for k in JSON_FIELDS if a.get(k) != b.get(k)]
                if _rounded(a.get("estimates")) != _rounded(b.get("estimates")):
                    differs.append("estimates (6 significant digits)")
            report["checks"].append({"file": name, "equal": not differs})
            if differs:
                report["failures"].append(
                    f"{name}: differs from the stored quick result in {', '.join(differs)}"
                )


def check_tables(report):
    run = _python(TABLES, "--check", *PAGES)
    report["checks"].append(
        {"tables": [p.name for p in PAGES], "equal": run.returncode == 0}
    )
    if run.returncode != 0:
        report["failures"].append(f"tables: {run.stderr.strip()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {"schema": "crn-craft.validation.kk13/0.1", "checks": [], "failures": []}
    check_quick(report)
    check_tables(report)
    report["passed"] = not report["failures"]
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    for failure in report["failures"]:
        print(failure, file=sys.stderr)
    print(f"kk13: {len(report['checks'])} checks, {len(report['failures'])} failures")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
