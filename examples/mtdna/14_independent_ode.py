"""A daily rhythm as independent ODEs (rules-model/0.4, hybrid method).

The three-class model with degradation that follows a 24-hour rhythm. The
rhythm is a harmonic oscillator, two continuous variables whose derivatives
read only parameters and each other. `"independent": true` lets `hybrid`
integrate them once for the whole run instead of in every replicate. The
model is written as JSON: rules-model/0.4 adds the field.
"""

import json
import math
import time

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["schema"] = "crn-craft.rules-model/0.4"
model["name"] = "mtdna_three_class_daily"
model["parameters"].update({"w": 2.0 * math.pi / 24.0, "amplitude": 0.3})


def with_rhythm(independent):
    document = json.loads(json.dumps(model))
    document["continuous"] = [
        {
            "name": "cosine",
            "initial": 1,
            "derivative": "-w * sine",
            "independent": independent,
        },
        {
            "name": "sine",
            "initial": 0,
            "derivative": "w * cosine",
            "independent": independent,
        },
    ]
    for rule in document["rules"]:
        if rule["name"].startswith("degradation"):
            rule["rate"] = "delta * (1 + amplitude * cosine)"
    return document


times = [6.0, 2_886.0, 5_766.0, 8_646.0]  # a quarter day past midnight
for independent in (True, False):
    start = time.perf_counter()
    result = rules.ensemble(
        with_rhythm(independent),
        end=times[-1],
        report_times=times,
        replicates=32,
        seed=14,
        threads=2,
        method="hybrid",
    ).simulate()
    seconds = time.perf_counter() - start
    mean = result.statistics["mean"]
    print(
        f"independent={independent!s:5}  {seconds:5.2f} s  "
        f"{result.metadata['events']:,} events"
    )
    print(mean[["N", "mutants"]].round(1))
    # The exact rhythm is cos(w t) = 0 and sin(w t) = 1 at these times. The
    # shared solution takes long steps over the whole year; per replicate,
    # the oscillator restarts at every event and its error accumulates.
    drift = max(
        abs(mean["cosine"] - [math.cos(model["parameters"]["w"] * t) for t in times])
    )
    print(f"largest error of the rhythm: {drift:.1e}")
