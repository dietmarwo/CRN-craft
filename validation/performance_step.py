"""One serial performance iteration: correctness, paired timings, retained evidence."""

import argparse
import json
from pathlib import Path
import subprocess


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("iteration")
    ap.add_argument("--before", type=Path, required=True)
    ap.add_argument("--repetitions", type=int, default=7)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = root / "validation/campaign"
    target = Path("/tmp/crn-campaign-target")

    def command(argv, suffix):
        with (out / f"{args.iteration}-{suffix}.log").open("w") as file:
            proc = subprocess.run(
                argv,
                cwd=root,
                stdout=file,
                stderr=subprocess.STDOUT,
                check=False,
                env={**__import__("os").environ, "CARGO_TARGET_DIR": str(target)},
            )
        if proc.returncode:
            raise RuntimeError(f"{args.iteration} {suffix} failed; see retained log")

    command(["cargo", "fmt", "--all", "--", "--check"], "fmt")
    command(["cargo", "test", "--release", "--workspace", "--locked"], "tests")
    command(
        [
            "cargo",
            "clippy",
            "--workspace",
            "--all-targets",
            "--locked",
            "--",
            "-D",
            "warnings",
        ],
        "clippy",
    )
    command(["cargo", "build", "--release", "--locked"], "build")
    command(
        [
            "python3",
            "validation/benchmark_campaign.py",
            "--before",
            str(args.before.resolve()),
            "--after",
            str(target / "release/crn-craft"),
            "--before-command",
            "simulate-cached",
            "--after-command",
            "simulate-cached",
            "--repetitions",
            str(args.repetitions),
            "--output",
            str(out / f"{args.iteration}-timings.json"),
        ],
        "benchmark",
    )
    report = json.loads((out / f"{args.iteration}-timings.json").read_text())
    print(
        json.dumps(
            {
                "iteration": args.iteration,
                "speedups": {
                    k: v["paired_speedup_median"]
                    for k, v in report["workloads"].items()
                },
                "identical": all(v["identical"] for v in report["workloads"].values()),
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
