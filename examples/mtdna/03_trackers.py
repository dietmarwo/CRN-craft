"""Event-level bookkeeping with trackers."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=500,
    seed=3,
    threads=2,
    method="mtdna_direct",
    trackers=[
        {"observable": "neutral_mutants", "keep_episodes": 3},
        {"observable": "advantaged_mutants", "keep_episodes": 3},
        # Time at which mutants first reach 60% of the molecules.
        {"observable": "mutants", "first_passage": "N > 0 && mutants >= 0.6 * N"},
    ],
)
result = run.simulate()
print("sampler:", result.metadata["sampler"])
# One DataFrame per tracker, one row per cell.
neutral = result.tracker_columns(0)
advantaged = result.tracker_columns(1)
passage = result.tracker_columns(2)
print(advantaged[["appearances", "extinctions", "maximum", "integral"]].head(3))
print("neutral lineages arising per cell:", neutral["appearances"].mean())
print("neutral lineages lost per cell:   ", neutral["extinctions"].mean())
print(
    "advantaged lineages arising/lost: ",
    advantaged["appearances"].mean(),
    advantaged["extinctions"].mean(),
)
print("largest advantaged class seen:    ", advantaged["maximum"].max())
# The integral runs until the replicate's stop or the horizon.
load = advantaged["integral"] / result.records["stop_time"]
print("time-averaged advantaged count, mean over cells:", load.mean())
reached = passage["first_passage"].dropna()  # NaN: never reached 60%
print(
    f"cells reaching 60% mutants: {len(reached)}; median day "
    f"{reached.median() / 24:.0f}"
)
# Episodes: each excursion of a class above zero, with start, end, peak.
first = next(o for o in result.observations if o["trackers"][1]["episodes"])
for episode in first["trackers"][1]["episodes"]:
    print("cell", first["id"], "advantaged episode:", episode)
