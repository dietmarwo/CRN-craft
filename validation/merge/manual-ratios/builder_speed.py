"""The builder form of the three-class mtDNA model after task E1.

Declared in validation/merge/RESULTS.md ("E1 (builder-built models as fast
as model files): declaration") before the first timing. Run under the CPU
lock, with the wheel of the E1 branch and `CRN_CRAFT_SOURCE` naming the
checkout:

    python validation/merge/manual-ratios/builder_speed.py --output builder-speed.json

Arms A, B, E and F as in `builder_form.py`; G the builder model of arm E
under `mtdna_direct`. One untimed warm-up run per arm, then three rounds in
the orders A B E F G, G F E B A, A B E F G. Times are the wall time of
`simulate()`. The event counts and mean digests of A, B and E are checked
against the 0.7.0 run recorded in `builder-form.json`.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import sys
import time

import numpy as np

import crn_craft as crn

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from builder_form import three_class  # noqa: E402

ROOT = HERE.parents[2]
ORDERS = ["ABEFG", "GFEBA", "ABEFG"]
ARMS = "ABEFG"
LABELS = {
    "A": "model file, mtdna_direct",
    "B": "model file, static_direct",
    "E": "builder model, static_direct",
    "F": "builder model, static_direct, compiled loop",
    "G": "builder model, mtdna_direct",
}
FAST = "chacha8-u53-ziggurat/rules-mtdna-direct/1"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    before = json.loads((HERE / "builder-form.json").read_text())
    document = json.loads(
        crn.load(str(ROOT / "models/rules/mtdna-three-class.rules.json"))
    )
    document["stop"] = "N > 0 && mutants >= 0.9 * N"
    builder = three_class(takeover=0.9)
    arguments = dict(
        end=25_920.0,
        report_times=24.0 * np.arange(0, 1081, 90),
        replicates=200,
        seed=1,
        threads=2,
        method="static_direct",
    )
    runs = {
        "A": crn.ensemble(document, **(arguments | {"method": "mtdna_direct"})),
        "B": crn.ensemble(document, **arguments),
        "E": crn.ensemble(builder, **arguments),
        "G": crn.ensemble(builder, **(arguments | {"method": "mtdna_direct"})),
    }
    runs["F"] = runs["E"]
    start = time.perf_counter()
    loop = crn.compile_loop(builder, progress=False)
    build_seconds = time.perf_counter() - start

    def simulate(arm):
        compiled = loop if arm == "F" else None
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
            "warnings": sorted({d["code"] for d in meta["diagnostics"]}),
            "mean_sha256": hashlib.sha256(
                result.statistics["mean"].to_numpy().tobytes()
            ).hexdigest(),
        }

    warmup = [simulate(arm) for arm in ARMS]
    rows = []
    for round_index, order in enumerate(ORDERS):
        for arm in order:
            row = simulate(arm)
            row["round"] = round_index + 1
            rows.append(row)
            print(json.dumps(row), flush=True)
    loop.verify(runs["E"])
    medians = {
        arm: statistics.median(r["seconds"] for r in rows if r["arm"] == arm)
        for arm in ARMS
    }
    events = {arm: next(r["events"] for r in rows if r["arm"] == arm) for arm in ARMS}
    per_event = {arm: medians[arm] / events[arm] for arm in ARMS}
    first_e = next(r for r in rows if r["arm"] == "E")
    recorded = {r["arm"]: r for r in before["warmup"]}
    checks = {
        "a_sampler_is_mtdna_direct": all(
            r["sampler"] == FAST for r in rows + warmup if r["arm"] == "A"
        ),
        "g_sampler_is_mtdna_direct_without_r051": all(
            r["sampler"] == FAST and "CRN-R051" not in r["warnings"]
            for r in rows + warmup
            if r["arm"] == "G"
        ),
        "e_and_f_equal_events_and_means": all(
            r["events"] == first_e["events"]
            and r["mean_sha256"] == first_e["mean_sha256"]
            for r in rows + warmup
            if r["arm"] in "EF"
        ),
        "loop_verify_passed": True,
        "every_run_completed": all(
            r["termination"] == "completed" for r in rows + warmup
        ),
    }
    identity = {
        "abe_events_and_means_as_on_13_november": all(
            r["events"] == recorded[r["arm"]]["events"]
            and r["mean_sha256"] == recorded[r["arm"]]["mean_sha256"]
            for r in rows + warmup
            if r["arm"] in "ABE"
        ),
    }
    target = per_event["E"] <= 1.5 * per_event["B"]
    report = {
        "schema": "crn-craft.manual-ratios/0.1",
        "workload": "guide chapter 12a, A first ensemble: three-class model, takeover stop at 90%, 25,920 h, reports every 90 days, 200 cells, seed 1, two threads",
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
            "E/B per event": per_event["E"] / per_event["B"],
            "E/A wall": medians["E"] / medians["A"],
            "E/A per event": per_event["E"] / per_event["A"],
            "E/B wall": medians["E"] / medians["B"],
            "B/A wall": medians["B"] / medians["A"],
            "E/F wall": medians["E"] / medians["F"],
            "F/A wall": medians["F"] / medians["A"],
            "G/A wall": medians["G"] / medians["A"],
            "G/E wall": medians["G"] / medians["E"],
        },
        "target_e1_t1_met": target,
        "checks": checks,
        "identity": identity,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["ratios"], indent=2))
    print(json.dumps({"target": target, **checks, **identity}))


if __name__ == "__main__":
    main()
