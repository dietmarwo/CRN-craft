"""Model variants: copy number, half-life, initial heteroplasmy, threshold."""

import json
import math

from crn_craft import CrnError, rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["parameters"].update({"K": 500, "w0": 400, "delta": math.log(2) / (20 * 24)})
# Start heteroplasmic: 400 wild type and 100 neutral mutants (class 1).
model["initial"].append(
    {"type": "mtDNA", "attributes": {"clone": 1, "s": "sw", "wt": 0}, "count": 100}
)
model["stop"] = "N > 0 && mutants >= 0.6 * N"  # 60% threshold
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
    replicates=300,
    seed=9,
    threads=2,
    method="mtdna_direct",
)
result = run.simulate()
stopped = result.records["termination"] == rules.TERMINATIONS.index("stopped")
print(result.metadata["sampler"])
print(f"cells past 60% within three years: {stopped.mean():.2%}")

# mtdna_direct certifies the model at preparation. A model it does not
# recognize runs Static Direct and says why in a CRN-R051 warning.
odd = json.loads(json.dumps(model))
odd["stop"] = "mutants > 250"
result = rules.ensemble(
    odd, end=100.0, report_times=[0.0, 100.0], replicates=2, method="mtdna_direct"
).simulate()
print(result.metadata["sampler"])
print([d["message"] for d in result.metadata["diagnostics"]])

# Mistakes are reported with a code and a JSON path before anything runs.
broken = json.loads(json.dumps(model))
broken["rules"][1]["rate"] = "max(0, K - N) / DD * sw"
try:
    rules.ensemble(broken, end=1.0, report_times=[0.0, 1.0])
except CrnError as error:
    print("rejected:", error)
