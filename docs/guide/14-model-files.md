# 14. Model files and species with attributes

Every model so far was built in Python. This chapter is for the two cases
that need a model as a JSON document: a model that exists as a file, and
species that carry attributes, which the builder cannot write. The setup
block builds the immigration–death model of [chapter 2](02-first-model.md).
The blocks of this chapter continue from it:

```python
import copy
import json
import pathlib
import time

import crn_craft as crn

model = crn.Model("immigration_death")
model.add_parameter({"birth": 4.0, "death": 0.2})
model.add_species("A", 0)
model.add_reaction("birth", [], ["A"], name="birth")
model.add_reaction("death * A", ["A"], [], name="death")
```

## What the builder writes

The builder writes a `crn-craft.model/1.0` document. Every species is a
*type* without attributes, and every reaction is a *rule* over the types:

```python
document = model.to_dict()
print(document["schema"])
print(json.dumps(document["rules"][1]))
```

```text
crn-craft.model/1.0
{"name": "death", "reactants": [{"label": "x", "type": "A"}], "propensity": "death * A", "products": []}
```

A string rate becomes a `propensity`, the whole expression. A number would
become a `rate`, the constant that multiplies the count of each reactant.

## Save and load a model file

`model.to_json()` is the text of a model file. `crn.load(path)` reads a
file, checks it, and returns its exact bytes. `crn.ensemble` takes those
bytes, a dictionary or a `Model`:

```python
path = pathlib.Path("immigration-death.model.json")
path.write_text(model.to_json())
data = crn.load(path)  # checked; raises crn.CrnError for an invalid file
loaded = crn.Model.from_dict(json.loads(data))
print(type(data).__name__, loaded)

settings = dict(end=30.0, report_times=[0.0, 30.0], replicates=400, seed=42,
                method="static_direct")
built = crn.ensemble(model, **settings).simulate()
from_file = crn.ensemble(data, **settings).simulate()
print("same numbers:", built.statistics["mean"].equals(from_file.statistics["mean"]))
print("same hash:", built.metadata["model_sha256"] == from_file.metadata["model_sha256"])

document = json.loads(data)  # a dictionary can be edited like any JSON
document["parameters"]["birth"] = 8.0
print(crn.ensemble(document, **settings).simulate().statistics["mean"]["A"].iloc[-1])
```

```text
bytes Model('immigration_death': 1 species, 2 reactions)
same numbers: True
same hash: False
39.61500000000001
```

The file and the built model give the same numbers. Their hashes differ,
because a hash names exact bytes and the file is indented. `Model.from_dict`
wraps a document as a `Model`, with `fsp`, `analyze` and `to_gillespy2`; the
`add_*` calls work on it only if the builder could have written it.

## Packaged examples

The package holds five example models: `immigration-death`,
`michaelis-menten`, `enzyme-mechanism`, `repressilator` and `vilar-noise`.
`crn.load_example` returns one as a `Model`. The packaged immigration–death
model runs the same paths as the model built above:

```python
example = crn.load_example("immigration-death")
for candidate in (model, example):
    result = crn.ensemble(
        candidate, end=30.0, report_times=[0.0, 30.0], replicates=400, seed=42,
        method="static_direct",
    ).simulate()
    print(result.metadata["events"], result.statistics["mean"]["A"].iloc[-1])
```

```text
88037 20.092500000000012
88037 20.092500000000012
```

## When species are not a fixed list

Some models cannot list their species in advance: every new mutant is its
own clone, a polymer has a length, a cell carries a lineage. For such
systems a species *type* carries numeric *attributes*. An *instance* is a
type with one value per attribute and a count. Instances appear when a rule
produces a new combination of values and disappear when their count reaches
zero. A rule is then a pattern: it stands for the reactions of every
instance that matches. The builder has no call for attributes, so such a
model is written as a document, in Python as a dictionary
([Model language](../reference/model-language.md)).

## A model document for three strains

Three strains of cells share `room` places. The resident strain 0 gives rise
to two variants, every strain dies at one rate, and strain 2 grows twice as
fast into the free room. Each cell is an instance of the type `Cell` with
the attribute `strain`:

```python
grow = "max(0, room - total) / crowding"
twin = {"type": "Cell", "copy": "c"}
strains = {
    "schema": "crn-craft.model/1.0",
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
    strains["rules"] += [
        {"name": f"die_{strain}", "reactants": cell, "rate": "loss", "products": []},
        {"name": f"grow_{strain}", "reactants": cell, "rate": f"{grow} * {rate}",
         "products": [twin, twin]},
    ]
for strain in (1, 2):
    resident = [{"label": "c", "type": "Cell", "where": "c.strain == 0"}]
    variant = {"type": "Cell", "attributes": {"strain": strain}}
    strains["rules"].append({"name": f"switch_to_{strain}", "reactants": resident,
                             "rate": f"{grow} * eps", "products": [twin, variant]})
print(crn.check(strains)["valid"], len(strains["rules"]), "rules")
pathlib.Path("three-strains.model.json").write_text(json.dumps(strains, indent=2) + "\n")
```

```text
True 8 rules
```

The parts of the document:

- **`initial`** lists instances: 200 cells with `strain` 0.
- **`observables`** are sums over the instances of one type. `where` selects
  instances (`variants`), and `weight` weighs them (`crowding` counts a cell
  of strain 2 twice). Results have one column per observable, and rates and
  the `stop` read observables by name.
- **A reactant** is a pattern: a `label`, a `type` and an optional `where`
  over the attributes. The `rate` is per matching cell, as a mass-action
  constant is.
- **A product** copies a reactant (`copy`), gives all attributes
  (`attributes`), or copies and overrides some.

## Run it

The three strains are a closed set of attribute values, so the model runs
under the static methods like a reaction network. This block loads the file
of the section before:

```python
strains = json.loads(crn.load("three-strains.model.json"))


def run(model, method):
    start = time.perf_counter()
    result = crn.ensemble(model, end=2000.0, report_times=[0.0, 2000.0], replicates=4000,
                          seed=1, threads=2, method=method).simulate()
    seconds = time.perf_counter() - start
    taken = (result.records["termination"] == "stopped").mean()
    print(f"{method:14s} {seconds:5.2f} s  variants took over: {taken:.3f}  "
          f"{result.metadata['sampler']}")
    return result


run(strains, "static_direct")
run(strains, "mtdna_direct")
# A strain that dies faster breaks the shared death rate, so the run
# falls back to Static Direct and says why.
faster = copy.deepcopy(strains)
next(r for r in faster["rules"] if r["name"] == "die_1")["rate"] = "2 * loss"
result = run(faster, "mtdna_direct")
print([d["message"] for d in result.metadata["diagnostics"]])
```

```text
static_direct   3.22 s  variants took over: 0.421  chacha8-u53-ziggurat/rules-static-direct/1
mtdna_direct    0.97 s  variants took over: 0.427  chacha8-u53-ziggurat/rules-mtdna-direct/1
mtdna_direct    3.00 s  variants took over: 0.434  chacha8-u53/rules-static-direct/1
['mtdna_direct ran Static Direct: degradation is not one constant rate shared by all classes']
```

`mtdna_direct` is an exact fast path for one structure: three classes in a
fixed room, one shared death rate, and growth into the free room. It reads
the structure from the rules, whatever the names, and certifies each model
before it runs. The first two runs agree within sampling error (with 4,000
replicates each, the standard error of their difference is about 0.011). A
model that does not fit runs Static Direct instead, and a `CRN-R051` warning
names the failed check: read `metadata["sampler"]` to see what ran. The
certificate also reads whole-rule propensities, as the builder writes them,
when the request is `crn-craft.run/1.0`, which `crn.ensemble` writes; so
the builder form of the three-class mtDNA model runs the fast path too
([page 12d](12d-mtdna-methods.md#the-fast-path)).
[Page 14a](14a-mtdna-model-files.md#the-fast-path) uses the fast path on the
mtDNA model files.

## New identities and random attributes

A product attribute may be a random draw: `new_id()` returns an identity
that no instance has had, and `choice(p, a, b)` is `a` with probability `p`,
else `b`. Here a growing cell sometimes founds a new strain with its own
identity and a drawn weight `w`. The observable `strains` counts distinct
instances:

```python
grow = "max(0, room - total) / crowding * g"
cell = {"label": "c", "type": "Cell"}
weighted = dict(cell, factor="c.w")  # the rate counts a cell w times
twin = {"type": "Cell", "copy": "c"}
founder = {"type": "Cell", "attributes": {"strain": "new_id()", "w": "choice(0.2, 2, 1)"}}
founders = {
    "schema": "crn-craft.model/1.0",
    "name": "founders",
    "parameters": {"room": 200, "loss": 0.01, "g": 0.5, "eps": 0.002},
    "types": [{"name": "Cell", "attributes": ["strain", "w"]}],
    "initial": [{"type": "Cell", "attributes": {"strain": 0, "w": 1}, "count": 200}],
    "observables": {
        "total": {"type": "Cell"},
        "crowding": {"type": "Cell", "weight": "x.w"},
        "strains": {"type": "Cell", "reduce": "distinct"},
    },
    "rules": [
        {"name": "die", "reactants": [cell], "rate": "loss", "products": []},
        {"name": "grow", "reactants": [weighted], "rate": grow, "products": [twin, twin]},
        {"name": "found", "reactants": [weighted], "rate": f"{grow} * eps",
         "products": [twin, founder]},
    ],
}
result = crn.ensemble(founders, end=2000.0, report_times=[0.0, 1000.0, 2000.0],
                      replicates=200, seed=1, threads=2).simulate()
print(result.metadata["sampler"])
print(result.statistics["mean"].round(2))
print("strains founded per replicate:", result.rule_events["found"].mean())
print("most strains alive at once:", result.records["peak_instances"].max())
try:
    crn.ensemble(founders, end=2000.0, report_times=[0.0, 2000.0], method="static_direct")
except crn.CrnError as error:
    print(error)
```

```text
chacha8-u53/rules-direct/2
        crowding  strains   total
time                             
0.0       200.00     1.00  200.00
1000.0    224.12     2.07  196.12
2000.0    270.94     2.27  196.27
strains founded per replicate: 7.435
most strains alive at once: 8
CRN-R020: random product attributes require the dynamic engine at /method
```

One `grow` rule now serves every strain. The set of instances changes
during the run, so the static methods refuse the model, and the exact tree
method, which runs when `method` is absent, takes it.
[Page 14a](14a-mtdna-model-files.md#every-lineage-with-the-clone-model)
follows every mutant lineage of the mtDNA model in this way.

## Limits and versions

A `crn-craft.model/1.0` document takes up to 4,096 types, 16,384 parameters
and 65,536 rules ([limits](../reference/model-language.md#limits)). Model
files of earlier releases (`crn-craft.rules-model/0.x`) are still read
([Contracts and versions](../reference/contracts.md)).

See also: [Model language](../reference/model-language.md),
[The mtDNA model files](14a-mtdna-model-files.md).
