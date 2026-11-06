"""mtdna_direct recognizes a model structure, not a file or names."""

import copy
import time

from crn_craft import rules

# Three strains share `room` places. The resident strain (0) gives rise to
# the two variants; every strain dies at one rate and grows into free room.
grow = "max(0, room - total) / crowding"
model = {
    "schema": "crn-craft.rules-model/0.1",
    "name": "three_strains",
    "parameters": {"room": 200, "loss": 0.01, "g": 0.5, "eps": 0.0002},
    "types": [{"name": "Cell", "attributes": ["strain"]}],
    "initial": [{"type": "Cell", "attributes": {"strain": 0}, "count": 200}],
    "observables": {
        "total": {"type": "Cell"},
        "crowding": {"type": "Cell", "weight": "1 + (x.strain == 2)"},
        "variants": {"type": "Cell", "where": "x.strain > 0"},
    },
    "stop": "total > 0 && variants >= 0.5 * total",
    "rules": [],
}
for strain, rate in ((0, "g"), (1, "g"), (2, "2 * g")):
    cell = [{"label": "c", "type": "Cell", "where": f"c.strain == {strain}"}]
    model["rules"] += [
        {"name": f"die_{strain}", "reactants": cell, "rate": "loss", "products": []},
        {
            "name": f"grow_{strain}",
            "reactants": cell,
            "rate": f"{grow} * {rate}",
            "products": [{"type": "Cell", "copy": "c"}, {"type": "Cell", "copy": "c"}],
        },
    ]
for strain in (1, 2):
    model["rules"].append(
        {
            "name": f"switch_to_{strain}",
            "reactants": [{"label": "c", "type": "Cell", "where": "c.strain == 0"}],
            "rate": f"{grow} * eps",
            "products": [
                {"type": "Cell", "copy": "c"},
                {"type": "Cell", "attributes": {"strain": strain}},
            ],
        }
    )


def run(model, method):
    start = time.perf_counter()
    result = rules.ensemble(
        model,
        end=2000.0,
        report_times=[0.0, 2000.0],
        replicates=4000,
        seed=1,
        threads=2,
        method=method,
    ).simulate()
    seconds = time.perf_counter() - start
    taken = (result.records["termination"] == "stopped").mean()
    print(
        f"{method:14s} {seconds:5.2f} s  variants took over: {taken:.3f}  "
        f"{result.metadata['sampler']}"
    )
    return result


run(model, "static_direct")
run(model, "mtdna_direct")
# A strain that dies faster breaks the shared degradation rate, so the run
# falls back to Static Direct and says why.
faster = copy.deepcopy(model)
next(r for r in faster["rules"] if r["name"] == "die_1")["rate"] = "2 * loss"
result = run(faster, "mtdna_direct")
print([d["message"] for d in result.metadata["diagnostics"]])
