"""P02 secondary control: pre-existing dimer topology, 64 streams, horizon 100."""

import argparse
import json
from pathlib import Path
import statistics
from benchmark_campaign import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", type=Path, required=True)
    ap.add_argument("--after", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    model = {
        "schema": "crn-craft.model/0.1",
        "name": "DimerCost",
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": {},
        "species": [
            {"name": "A", "initial_count": 200},
            {"name": "B", "initial_count": 0},
        ],
        "reactions": [
            {
                "name": "dimer",
                "reactants": {"A": 2},
                "products": {"B": 1},
                "rate": {"kind": "mass_action", "coefficient": 0.01},
            },
            {
                "name": "split",
                "reactants": {"B": 1},
                "products": {"A": 2},
                "rate": {"kind": "mass_action", "coefficient": 1.0},
            },
        ],
    }
    before = []
    after = []
    for _ in range(7):
        before += run_case(args.before.resolve(), "run", model, 1)["rows"]
        after += run_case(args.after.resolve(), "run", model, 1)["rows"]
    identical = len({r["numerical_sha256"] for r in before + after}) == 1
    result = {
        "secondary_control": True,
        "topology": "M3 dimer, larger initial conservation class for timing",
        "identical": identical,
        "before": before,
        "after": after,
        "paired_median_speedup": statistics.median(
            b["process_seconds"] / a["process_seconds"]
            for b, a in zip(before, after, strict=True)
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if not identical:
        raise RuntimeError("dimer numerical payload changed")
    print(json.dumps({"identical": True, "speedup": result["paired_median_speedup"]}))


if __name__ == "__main__":
    main()
