"""Seeds, threads and reproducibility."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"


def run(seed, threads):
    ensemble = rules.ensemble(
        model,
        end=25_920.0,
        report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
        replicates=64,
        seed=seed,
        threads=threads,
        method="mtdna_direct",
        keep_paths="all",  # every cell's series, as result.paths
    )
    return ensemble.simulate()


one, two, other = run(7, 1), run(7, 2), run(8, 2)
same = one.records.equals(two.records) and one.rule_events.equals(two.rule_events)
print("1 vs 2 threads, same seed, identical records:", same)
print("same statistics:", one.statistics["mean"].equals(two.statistics["mean"]))
print("same paths:", one.paths.equals(two.paths))
print(
    "seed 7 vs 8, identical stop times:",
    one.records["stop_time"].equals(other.records["stop_time"]),
)
# Replicate i always uses the same stream, whatever the ensemble size.
small = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
    replicates=8,
    seed=7,
    threads=1,
    method="mtdna_direct",
).simulate()
print(
    "first 8 of 64 equal an 8-cell run:",
    small.records["events"].equals(one.records["events"].head(8)),
)
paths = one.paths  # cells x report times x observables, an xarray DataArray
print(dict(paths.sizes))
print(paths.sel(replicate=0, variable=["N", "mutants"]).to_pandas())
