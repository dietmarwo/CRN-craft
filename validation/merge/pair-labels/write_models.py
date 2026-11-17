"""Write the models of the pair-label law check (validation/merge/RESULTS.md).

Each row is a crn-craft.model/1.0 model whose pair rule ("x replaces y")
has products that depend on the labels, simulated with the tree method, and
an attribute-free network with the same law for the observed counts, which
gives the sample reference (static_direct) and the exact law (fsp):

- `pair-neutral`: 10 + 10 molecules of clones 0 and 1, one rule "x
  replaces y" at rate g per unordered pair. With random labels the counts
  C0, C1 follow C0 + C1 -> 2 C1 and C0 + C1 -> 2 C0, each at g/2.
- `pair-mutation`: 20 wild-type molecules (clone 0); "x replaces y" at
  g(1 - p) and "x replaces y with a new mutant" at g p (`new_id()`). The
  counts W (clone 0) and M (clone > 0) follow W + M -> 2 M at g/2,
  W + M -> 2 W at g(1 - p)/2 and 2 W -> W + M at g p.

Standard library only:

    python validation/merge/pair-labels/write_models.py
"""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
G, P = 0.1, 0.2
TIMES = [0.0, 0.25, 0.5, 1.0, 2.0]
PAIR = [{"label": "x", "type": "m"}, {"label": "y", "type": "m"}]
KEEP_X = [{"type": "m", "copy": "x"}, {"type": "m", "copy": "x"}]


def neutral():
    return {
        "schema": "crn-craft.model/1.0",
        "name": "pair_neutral",
        "parameters": {"g": G},
        "types": [{"name": "m", "attributes": ["clone"]}],
        "initial": [
            {"type": "m", "attributes": {"clone": 0}, "count": 10},
            {"type": "m", "attributes": {"clone": 1}, "count": 10},
        ],
        "observables": {
            "C0": {"type": "m", "where": "x.clone == 0"},
            "C1": {"type": "m", "where": "x.clone == 1"},
        },
        "rules": [
            {"name": "x_replaces_y", "reactants": PAIR, "rate": "g", "products": KEEP_X}
        ],
    }


def mutation():
    return {
        "schema": "crn-craft.model/1.0",
        "name": "pair_mutation",
        "parameters": {"g": G, "p": P},
        "types": [{"name": "m", "attributes": ["clone"]}],
        "initial": [{"type": "m", "attributes": {"clone": 0}, "count": 20}],
        "observables": {
            "W": {"type": "m", "where": "x.clone == 0"},
            "M": {"type": "m", "where": "x.clone > 0"},
        },
        "rules": [
            {"name": "x_replaces_y", "reactants": PAIR, "rate": "g * (1 - p)",
             "products": KEEP_X},
            {"name": "x_mutant_replaces_y", "reactants": PAIR, "rate": "g * p",
             "products": [{"type": "m", "copy": "x"},
                          {"type": "m", "attributes": {"clone": "new_id()"}}]},
        ],
    }


def network(name, species, initial, rules, parameters):
    return {
        "schema": "crn-craft.model/1.0",
        "name": name,
        "parameters": parameters,
        "types": [{"name": s} for s in species],
        "initial": [{"type": s, "count": initial[s]} for s in species],
        "observables": {s: {"type": s} for s in species},
        "rules": [
            {"name": rule,
             "reactants": [{"label": label, "type": t}
                           for label, t in zip(["x", "y"], reactants)],
             "rate": rate,
             "products": [{"type": t} for t in products]}
            for rule, reactants, rate, products in rules
        ],
    }


def neutral_network():
    return network(
        "pair_neutral_counts", ["C0", "C1"], {"C0": 10, "C1": 10},
        [("c1_replaces_c0", ["C0", "C1"], "g / 2", ["C1", "C1"]),
         ("c0_replaces_c1", ["C0", "C1"], "g / 2", ["C0", "C0"])],
        {"g": G},
    )


def mutation_network():
    return network(
        "pair_mutation_counts", ["W", "M"], {"W": 20, "M": 0},
        [("m_replaces_w", ["W", "M"], "g / 2", ["M", "M"]),
         ("w_replaces_m", ["W", "M"], "g * (1 - p) / 2", ["W", "W"]),
         ("w_mutates", ["W", "W"], "g * p", ["W", "M"])],
        {"g": G, "p": P},
    )


def write(path, model):
    data = (json.dumps(model, indent=1) + "\n").encode()
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def main():
    rows = []
    for name, model, counts, species in [
        ("pair-neutral", neutral(), neutral_network(), ["C0", "C1"]),
        ("pair-mutation", mutation(), mutation_network(), ["W", "M"]),
    ]:
        path = HERE / f"{name}.model.json"
        reference = HERE / f"{name}.counts.model.json"
        rows.append({
            "name": name,
            "path": str(path.relative_to(ROOT)),
            "sha256": write(path, model),
            "reference_path": str(reference.relative_to(ROOT)),
            "reference_sha256": write(reference, counts),
            "form": "rules",
            "stop": False,
            "species_names": species,
            "end": TIMES[-1],
            "report_times": TIMES,
            "replicates": 200000,
            "chunk": 10000,
            "reference_seed": 1000000,
            "gate_seed": 2000000,
            "max_instances": 1000,
            "fsp": {"bounds": {s: 20 for s in species}, "max_states": 10000,
                    "tolerance": 1e-12},
        })
    manifest = {
        "schema": "crn-craft.merge-models/0.1",
        "declared_in": "validation/merge/RESULTS.md (0.8 pair labels: law check declaration)",
        "seeds": "chunk k (0-based) of a model uses root seed reference_seed + k "
                 "(counts network, static_direct reference) or gate_seed + k "
                 "(engine under test, tree method)",
        "rows": rows,
    }
    (HERE / "models.json").write_text(json.dumps(manifest, indent=1) + "\n")
    for row in rows:
        print(row["name"])


if __name__ == "__main__":
    main()
