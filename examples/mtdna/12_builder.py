"""The three-class model built with rules.Model instead of JSON."""

import json

import numpy as np
from crn_craft import rules

# Species are counts: wild type W, neutral mutants M, advantaged mutants A.
# Numeric rates are mass-action constants; strings are whole propensities
# over parameters and species counts.
model = rules.Model("mtdna_three_species")
model.add_parameter(
    {
        "K": 1000,
        "delta": 0.0028881132523331052,
        "sw": 1.0,
        "sm": 1.5,
        "mu": 3e-05,
        "fracAdv": 0.2,
    }
)
model.add_species({"W": 1000, "M": 0, "A": 0})
# Replication fills the vacancy K - N, shared by the competition term.
vacancy = "max(0, K - W - M - A) / (W + sm * (M + A))"
model.add_reaction(
    f"{vacancy} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt"
)
model.add_reaction(
    f"{vacancy} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt"
)
model.add_reaction(
    f"{vacancy} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv"
)
model.add_reaction(f"{vacancy} * sw * M", ["M"], ["M", "M"], "replication_mt")
model.add_reaction(f"{vacancy} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
for species in ("W", "M", "A"):
    model.add_reaction(
        model.parameters["delta"], [species], [], f"degradation_{species}"
    )
print(model)
print("first rule:", json.dumps(model.to_dict()["rules"][0]))

days = np.arange(0, 1081, 180)
settings = dict(
    end=25_920.0,
    report_times=24.0 * days,
    replicates=200,
    threads=2,
    keep_paths="all",
    method="static_direct",
)
built = rules.ensemble(model, seed=1, **settings).simulate()
# The JSON file of the same model, with independent random streams.
reference = rules.ensemble(
    rules.load("models/rules/mtdna-three-class.rules.json"), seed=2, **settings
).simulate()


def mutant_fraction(result, mutants, total):
    """Mean and standard error over cells of each cell's mutant fraction."""
    paths = result.paths  # cells x report times x observables, a DataArray
    count = paths.sel(variable=mutants).sum("variable")
    fraction = count / paths.sel(variable=total).sum("variable")  # cells x times
    cells = fraction.sizes["replicate"]
    mean = fraction.mean("replicate")
    return mean.values, fraction.std("replicate", ddof=1).values / np.sqrt(cells)


ours = mutant_fraction(built, ["M", "A"], ["W", "M", "A"])
theirs = mutant_fraction(reference, ["mutants"], ["N"])
print("sampler:", built.metadata["sampler"])
print(
    f"events: builder {built.metadata['events']:,}, JSON {reference.metadata['events']:,}"
)
print("mean mutant fraction per cell (standard error)")
for day, a, sa, b, sb in zip(days, *ours, *theirs, strict=True):
    print(f"day {day:4d}: builder {a:.4f} ({sa:.4f})   JSON {b:.4f} ({sb:.4f})")
