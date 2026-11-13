"""The builder form of the three-class mtDNA model against its file form.

Declared in validation/merge/RESULTS.md ("Docs 0.7.1: the builder form of
the mtDNA model") before the first timing. Run under the CPU lock, with the
wheel of the checkout and `CRN_CRAFT_SOURCE` naming the checkout:

    python validation/merge/manual-ratios/builder_form.py --output builder-form.json

Arms: A the model file under `mtdna_direct`; B the model file under
`static_direct`; E the builder model of guide chapter 12 under
`static_direct`; F E with the loop of `crn.compile_loop`. One untimed warm-up
run per arm, then three rounds in the orders A B E F, F E B A, A B E F.
Times are the wall time of `simulate()`.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import time

import numpy as np

import crn_craft as crn

ROOT = Path(__file__).resolve().parents[3]
ORDERS = ["ABEF", "FEBA", "ABEF"]
LABELS = {
    "A": "model file, mtdna_direct",
    "B": "model file, static_direct",
    "E": "builder model, static_direct",
    "F": "builder model, static_direct, compiled loop",
}


def three_class(
    K=1000, w0=1000, m0=0, mu=3e-5, sm=1.5, half_life_days=10, takeover=None
):
    """The three-class mtDNA model of guide chapter 12."""
    model = crn.Model("mtdna_three_class")
    model.add_parameter(
        {
            "K": K,
            "delta": math.log(2) / (24 * half_life_days),
            "sw": 1.0,
            "sm": sm,
            "mu": mu,
            "fracAdv": 0.2,
        }
    )
    model.add_species({"W": w0, "M": m0, "A": 0})
    free = "max(0, K - W - M - A) / (W + sm * (M + A))"
    model.add_reaction(
        f"{free} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt"
    )
    model.add_reaction(
        f"{free} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt"
    )
    model.add_reaction(
        f"{free} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv"
    )
    model.add_reaction(f"{free} * sw * M", ["M"], ["M", "M"], "replication_mt")
    model.add_reaction(f"{free} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
    for label, species in (("wt", "W"), ("mt", "M"), ("mtAdv", "A")):
        model.add_reaction(f"delta * {species}", [species], [], f"degradation_{label}")
    if takeover is not None:
        model.set_stop(f"W + M + A > 0 && M + A >= {takeover} * (W + M + A)")
    return model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
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
            "mean_sha256": hashlib.sha256(
                result.statistics["mean"].to_numpy().tobytes()
            ).hexdigest(),
        }

    warmup = [simulate(arm) for arm in "ABEF"]
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
        for arm in "ABEF"
    }
    events = {arm: next(r["events"] for r in rows if r["arm"] == arm) for arm in "ABEF"}
    per_event = {arm: medians[arm] / events[arm] for arm in "ABEF"}
    first_e = next(r for r in rows if r["arm"] == "E")
    checks = {
        "a_sampler_is_mtdna_direct": all(
            r["sampler"] == "chacha8-u53-ziggurat/rules-mtdna-direct/1"
            for r in rows + warmup
            if r["arm"] == "A"
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
            "E/A wall": medians["E"] / medians["A"],
            "E/A per event": per_event["E"] / per_event["A"],
            "E/B wall": medians["E"] / medians["B"],
            "B/A wall": medians["B"] / medians["A"],
            "E/F wall": medians["E"] / medians["F"],
            "F/A wall": medians["F"] / medians["A"],
        },
        "checks": checks,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["ratios"], indent=2))
    print(json.dumps(checks))


if __name__ == "__main__":
    main()
