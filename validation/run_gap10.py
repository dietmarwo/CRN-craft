"""Per-iteration numerical identity and predeclared paired performance gate."""

import argparse
import json
import math
from pathlib import Path
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
SPECS = {
    "W1": (100.0, 64),
    "W2": (100.0, 64),
    "W3": (100.0, 64),
    "W4": (100.0, 64),
    "W5": (20.0, 16),
    "W6": (20.0, 16),
    "W7": (20.0, 16),
}
ROWS = [
    ("W1", "legacy_cached"),
    ("W2", "direct_cr"),
    ("W2", "sorting_direct"),
    ("W3", "legacy_cached"),
    ("W3", "sorting_direct"),
    ("W4", "cached_direct"),
    ("W4", "rssa_cr"),
    ("W5", "direct_cr"),
    ("W5", "rssa_cr"),
    ("W6", "rssa_cr"),
    ("W7", "rssa_cr"),
]
TREE = {"direct_cr", "sorting_direct", "rssa_cr"}


def target(iteration, w, m):
    return {
        "C01": True,
        "C02": not m.startswith("legacy_"),
        "C03": m in TREE,
        "C04": m in TREE,
        "C05": m == "rssa_cr" and w in ["W5", "W6", "W7"],
        "C06": (w, m) == ("W4", "rssa_cr"),
        "C07": m in ["direct_cr", "rssa_cr"],
        "C08": m == "sorting_direct",
        "C09": (w, m)
        in [
            ("W1", "legacy_cached"),
            ("W2", "direct_cr"),
            ("W5", "direct_cr"),
            ("W5", "rssa_cr"),
            ("W6", "rssa_cr"),
            ("W7", "rssa_cr"),
        ],
        "C10": True,
    }[iteration]


def geometric(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


def bench(binary, w, m, paths, repeats):
    end, default = SPECS[w]
    wait = "inverse" if m.startswith("legacy_") else "ziggurat"
    command = [
        str(binary),
        str(ROOT / f"validation/ssa-acceleration/models/{w}.model.json"),
        str(end),
        str(default if paths is None else paths),
        m,
        wait,
        str(repeats),
    ]
    value = json.loads(subprocess.check_output(command, timeout=180))
    return value


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parent", type=Path, required=True)
    p.add_argument("--candidate", type=Path, required=True)
    p.add_argument("--iteration", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--identity-only", action="store_true")
    args = p.parse_args()
    report = dict(
        protocol="docs/decisions/julia-gap10-2026-11-02.md",
        iteration=args.iteration,
        identity=[],
        rows={},
        started=time.time(),
    )

    def save():
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    for w in SPECS:
        for m in [
            "legacy_cached",
            "cached_direct",
            "direct_cr",
            "sorting_direct",
            "rssa_cr",
        ]:
            a = bench(args.parent, w, m, 4, 1)
            b = bench(args.candidate, w, m, 4, 1)
            identical = a["numerical_digest"] == b["numerical_digest"] and (
                not m.startswith("legacy_") or a["digest"] == b["digest"]
            )
            report["identity"].append(
                dict(
                    workload=w,
                    method=m,
                    passed=identical,
                    parent=a["numerical_digest"],
                    candidate=b["numerical_digest"],
                )
            )
            save()
            if not identical:
                raise RuntimeError(f"identity failed {w}/{m}")
    report["identity_passed"] = True
    save()
    if args.identity_only:
        return
    for w, m in ROWS:
        report["rows"][w + "/" + m] = dict(
            target=target(args.iteration, w, m), parent=[], candidate=[], ratios=[]
        )
    for round_index in range(5):
        offset = (round_index * 3) % len(ROWS)
        order = ROWS[offset:] + ROWS[:offset]
        if round_index % 2:
            order = list(reversed(order))
        for j, (w, m) in enumerate(order):
            row = report["rows"][w + "/" + m]
            pair = {}
            for kind in (
                ["parent", "candidate"]
                if (round_index + j) % 2 == 0
                else ["candidate", "parent"]
            ):
                value = bench(getattr(args, kind), w, m, None, 2)
                row[kind].append(value)
                pair[kind] = statistics.median(value["seconds"])
            if (
                row["parent"][-1]["numerical_digest"]
                != row["candidate"][-1]["numerical_digest"]
            ):
                raise RuntimeError("full-path identity failure")
            row["ratios"].append(pair["parent"] / pair["candidate"])
        save()
        print(args.iteration, "round", round_index + 1, flush=True)
    for key, row in report["rows"].items():
        if statistics.median(row["ratios"]) < 0.95:
            w, m = key.split("/")
            for i in range(5):
                pair = {}
                for kind in (
                    ["candidate", "parent"] if i % 2 else ["parent", "candidate"]
                ):
                    value = bench(getattr(args, kind), w, m, None, 2)
                    row[kind].append(value)
                    pair[kind] = statistics.median(value["seconds"])
                row["ratios"].append(pair["parent"] / pair["candidate"])
            save()
        row["speedup"] = statistics.median(row["ratios"])
    targets = [r for r in report["rows"].values() if r["target"]]
    target_gain = geometric(r["speedup"] for r in targets)
    rounds = [geometric(r["ratios"][i] for r in targets) for i in range(5)]
    guard_gain = geometric(r["speedup"] for r in report["rows"].values())
    report.update(
        target_speedup=target_gain,
        target_rounds=rounds,
        guard_speedup=guard_gain,
        accepted=target_gain >= 1.02
        and sum(v > 1 for v in rounds) >= 4
        and guard_gain >= 0.995
        and all(r["speedup"] >= 0.95 for r in report["rows"].values()),
    )
    save()
    print(
        json.dumps(
            {
                k: report[k]
                for k in [
                    "iteration",
                    "accepted",
                    "target_speedup",
                    "guard_speedup",
                    "target_rounds",
                ]
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
