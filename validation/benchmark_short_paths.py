"""Secondary startup/queue control: same two-channel model, horizon 0.1."""

import argparse
import json
from pathlib import Path
import statistics
from benchmark_campaign import models, run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", type=Path, required=True)
    ap.add_argument("--after", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    model = next(models())[1]
    before = []
    after = []
    for _ in range(21):
        before += run_case(args.before.resolve(), "run", model, 1, args.threads, 0.1)[
            "rows"
        ]
        after += run_case(args.after.resolve(), "run", model, 1, args.threads, 0.1)[
            "rows"
        ]
    identical = len({r["numerical_sha256"] for r in before + after}) == 1
    result = {
        "secondary_control": True,
        "horizon": 0.1,
        "threads": args.threads,
        "replicates": 64,
        "before": before,
        "after": after,
        "identical": identical,
        "paired_median_speedup": statistics.median(
            b["process_seconds"] / a["process_seconds"]
            for b, a in zip(before, after, strict=True)
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if not identical:
        raise RuntimeError("short-path payload changed")
    print(json.dumps({"identical": True, "speedup": result["paired_median_speedup"]}))


if __name__ == "__main__":
    main()
