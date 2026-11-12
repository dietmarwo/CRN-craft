"""Speed ratios quoted in the mtDNA user manual, on its section 3 workload.

Declared in validation/merge/RESULTS.md ("Review: manual speed ratios")
before the first timing. Run under the CPU lock, with the wheel and the
release build of the checkout, `CRN_CRAFT_SOURCE` naming the checkout:

    python validation/merge/manual-ratios/bench.py --output timings.json

Arms: A `mtdna_direct`; B `static_direct` with its defaults; C B with the
generated loop of `crn.compile_loop`; D Static Direct with inverse waits and
per-event trackers. One untimed warm-up run per arm, then three rounds in
the orders A B C D, D C B A, A B C D. Times are the wall time of
`simulate()`.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import time

import numpy as np

import crn_craft as crn

ROOT = Path(__file__).resolve().parents[3]
ORDERS = ["ABCD", "DCBA", "ABCD"]
METHODS = {
    "A": "mtdna_direct",
    "B": "static_direct",
    "C": "static_direct",
    "D": {"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False},
}
LABELS = {
    "A": "mtdna_direct",
    "B": "static_direct (defaults)",
    "C": "static_direct, compiled loop",
    "D": "static_direct, inverse waits and per-event trackers",
}


def workload():
    """The model and ensemble arguments of the manual's section 3 block."""
    model = json.loads(
        crn.load(str(ROOT / "models/rules/mtdna-three-class.rules.json"))
    )
    model["stop"] = "N > 0 && mutants >= 0.9 * N"
    days = np.arange(0, 1081, 90)
    arguments = dict(
        end=25_920.0, report_times=24.0 * days, replicates=200, seed=1, threads=2
    )
    return model, arguments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    model, arguments = workload()
    runs = {
        arm: crn.ensemble(model, method=METHODS[arm], **arguments) for arm in "ABCD"
    }
    start = time.perf_counter()
    loop = crn.compile_loop(model, progress=False)
    build_seconds = time.perf_counter() - start

    def simulate(arm):
        compiled = loop if arm == "C" else None
        start = time.perf_counter()
        result = runs[arm].simulate(compiled=compiled)
        seconds = time.perf_counter() - start
        meta = result.metadata
        return {
            "arm": arm,
            "seconds": seconds,
            "events": meta["events"],
            "sampler": meta["sampler"],
            "termination": meta["termination"],
            "mean_sha256": hashlib.sha256(
                result.statistics["mean"].to_numpy().tobytes()
            ).hexdigest(),
        }

    warmup = [simulate(arm) for arm in "ABCD"]
    rows = []
    for round_index, order in enumerate(ORDERS):
        for arm in order:
            row = simulate(arm)
            row["round"] = round_index + 1
            rows.append(row)
            print(json.dumps(row), flush=True)
    loop.verify(runs["B"])
    medians = {
        arm: statistics.median(r["seconds"] for r in rows if r["arm"] == arm)
        for arm in "ABCD"
    }
    events = {
        arm: rows[[r["arm"] for r in rows].index(arm)]["events"] for arm in "ABCD"
    }
    per_event = {arm: medians[arm] / events[arm] for arm in "ABCD"}
    checks = {
        "a_sampler_is_mtdna_direct": all(
            r["sampler"] == "chacha8-u53-ziggurat/rules-mtdna-direct/1"
            for r in rows + warmup
            if r["arm"] == "A"
        ),
        "b_and_c_equal_events_and_means": all(
            r["events"] == events["B"]
            and r["mean_sha256"]
            == next(x for x in rows if x["arm"] == "B")["mean_sha256"]
            for r in rows + warmup
            if r["arm"] in "BC"
        ),
        "loop_verify_passed": True,
        "every_run_completed": all(
            r["termination"] == "completed" for r in rows + warmup
        ),
    }
    report = {
        "schema": "crn-craft.manual-ratios/0.1",
        "workload": "docs/mtdna-user-manual.md section 3: three-class model, takeover stop, 25,920 h, reports every 90 days, 200 cells, seed 1, two threads",
        "labels": LABELS,
        "python": platform.python_version(),
        "crn_craft": crn.__version__,
        "loop_build_seconds": build_seconds,
        "warmup": warmup,
        "rows": rows,
        "median_seconds": medians,
        "events": events,
        "median_seconds_per_event": per_event,
        "ratios": {
            "B/A wall": medians["B"] / medians["A"],
            "B/A per event": per_event["B"] / per_event["A"],
            "D/A wall": medians["D"] / medians["A"],
            "D/A per event": per_event["D"] / per_event["A"],
            "B/C wall": medians["B"] / medians["C"],
        },
        "checks": checks,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["ratios"], indent=2))
    print(json.dumps(checks))


if __name__ == "__main__":
    main()
