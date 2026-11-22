# 12. Model files and species with attributes

Every model so far had species that are plain counts. Some models need
more: each molecule carries attributes, such as the clone a mutant belongs
to, and new kinds of molecules appear during the run. This chapter builds
such a model with `crn.Model`, then turns to the model as a JSON document:
what the builder writes, model files, and documents written by hand, which
are the special case. The setup block builds the immigration–death model of
[chapter 2](02-first-model.md). The blocks of this chapter continue from it:

```python
import copy
import json
import math
import pathlib
import time

import crn_craft as crn

model = crn.Model("immigration_death")
model.add_parameter({"birth": 4.0, "death": 0.2})
model.add_species("A", 0)
model.add_reaction("birth", [], ["A"], name="birth")
model.add_reaction("death * A", ["A"], [], name="death")
```

## When species are not a fixed list

Some models cannot list their species in advance: every new mutant is its
own clone, a polymer has a length, a cell carries a lineage. For such
systems a species *type* carries numeric *attributes*. An *instance* is a
type with one value per attribute and a count. Instances appear when a rule
produces a new combination of values and disappear when their count reaches
zero. A reaction is then a pattern: it stands for the reactions of every
instance that matches. The builder declares such a species with
`add_species(name, attributes=[...])`; the
[model language](../reference/model-language.md) describes the document it
writes.

## Species with attributes in the builder

A cell holds 1000 mtDNA molecules. Each molecule degrades with a half-life
of ten days, free places are refilled by copies of a random survivor, and a
copy is a new mutant with probability `pmut`. Every mutation founds its own
clone, so the number of clones is not known in advance. The molecules carry
the attribute `clone`, with 0 for the wild type, and time is in hours:

```python
clones = crn.Model("mtdna_clones")
clones.add_parameter({"K": 1000, "delta": math.log(2) / 240, "pmut": 1e-3, "r": 10.0})
clones.add_species("mtDNA", attributes=["clone"])  # molecules carry a clone id
clones.add_initial("mtDNA", 1000, clone=0)  # 1000 molecules of clone 0, the wild type
clones.add_observable("mutants", "mtDNA", where="x.clone > 0")
clones.add_observable("types", "mtDNA", reduce="distinct")
clones.add_observable("mutant_types", "mtDNA", where="x.clone > 0", reduce="distinct")
clones.add_observable("largest", "mtDNA", where="x.clone > 0", reduce="max_count")
refill = "r * max(0, K - mtDNA) / mtDNA"  # per survivor; mtDNA is the total count
clones.add_reaction(crn.mass_action("delta"), ["mtDNA"], [], name="degradation")
clones.add_reaction(crn.mass_action(f"{refill} * (1 - pmut)"), ["mtDNA"],
                    [crn.copy("x"), crn.copy("x")], name="replication")
clones.add_reaction(crn.mass_action(f"{refill} * pmut"), ["mtDNA"],
                    [crn.copy("x"), crn.new("mtDNA", clone="new_id()")], name="mutation")
print(clones, crn.check(clones)["valid"])

year = 24 * 360.0
result = crn.ensemble(clones, end=year, report_times=[0.0, year / 2, year],
                      replicates=100, seed=1, threads=2).simulate()
print(result.metadata["sampler"])
print(result.statistics["mean"].round(2))
```

```text
Model('mtdna_clones': 1 species, 3 reactions) True
chacha8-u53/rules-direct/2
        largest    mtDNA  mutant_types  mutants  types
time                                                  
0.0        0.00  1000.00          0.00     0.00   1.00
4320.0     7.47   999.54          2.35    11.25   3.35
8640.0    12.82   999.67          3.15    20.26   4.15
```

After a year a cell holds about 20 mutant molecules in three clones, and
its largest clone has about 13 molecules. The calls:

- **`add_species("mtDNA", attributes=["clone"])`** declares the species and
  its attribute. Its molecules at time 0 come from
  **`add_initial("mtDNA", 1000, clone=0)`**, one call per combination of
  attribute values.
- **A reactant** given as a species name is any molecule of it, labeled `x`
  (the second reactant is `y`). A replication consumes `x` and gives two
  copies of it, **`crn.copy("x")`**, so the clone keeps its identity. A
  mutation gives one copy and a molecule with every attribute set,
  **`crn.new("mtDNA", clone="new_id()")`**; `new_id()` is an identity no
  molecule has had in this run.
- **`crn.mass_action(expression)`** is a rate per reactant molecule, as a
  number is: the degradation fires at `delta` times the number of molecules,
  and each molecule is equally likely to be the one degraded. The refill
  divides by the count of survivors, so all replications together fill the
  free places at rate `r * (K - mtDNA)`.

The exact tree method runs the model, which happens when `method` is left
out, because the set of clones changes during the run. Each cell had about
50,000 events in the year.

### Observables

Every species has an observable under its own name, its count: here
`mtDNA`. `add_observable(name, species, where=None, weight=None,
reduce="sum")` adds others. `where` selects instances by their attributes,
`weight` gives each instance a number (both read the instance as `x`), and
`reduce` says how the matching instances become one number:

| `reduce` | Value |
| --- | --- |
| `"sum"` (default) | the sum of weight × count: with no weight, the number of matching molecules (`mutants`) |
| `"distinct"` | the sum of the weights, each instance once: with no weight, the number of matching instances (`types`, `mutant_types`) |
| `"max"`, `"min"` | the largest or smallest weight; 0 when no instance matches |
| `"max_count"` | the largest count among the matching instances, times the weight; 0 when none matches (`largest`, the size of the largest mutant clone) |

Observables are columns of the results, and rates, the stop condition and
trackers read them by name, as the refill rate reads `mtDNA`.

### Reactants, products and rates

| In `add_reaction` | Meaning |
| --- | --- |
| `"mtDNA"` as a reactant | any molecule of the species; labeled `x`, the second `y` |
| `crn.pattern("mtDNA", where="x.clone > 0", factor="x.size")` | the molecules that satisfy `where`, each counted `factor` times |
| `crn.copy("x")`, `crn.copy("x", size="x.size + 1")` | a product with the attributes of reactant `x`, some of them overridden |
| `crn.new("mtDNA", clone="new_id()")` | a product with every attribute given |
| a number, or `crn.mass_action("expression")` | a rate per combination of reactant molecules |
| a plain string | the reaction's whole propensity |

Attribute values are numbers or expressions over parameters, the reactants'
attributes (`x.clone`) and random draws such as `new_id()` or
`choice(p, a, b)`. A product of a species with attributes must say where
its attributes come from, so a bare species name is refused:

```python
try:
    clones.add_reaction(0.01, ["mtDNA"], ["mtDNA", "mtDNA"])
except ValueError as error:
    print(error)
```

```text
species 'mtDNA' has attributes, so a product of it says where they come from: crn_craft.copy(label) copies a reactant, crn_craft.new('mtDNA', clone=...) gives every attribute
```

A plain string rate keeps its meaning from [chapter 2](02-first-model.md):
the whole propensity. On a species with attributes the reaction then fires
at that rate and takes a molecule uniformly, whatever its clone, so
`"delta * mtDNA"` and `crn.mass_action("delta")` give the same degradation.
Only `mass_action` (or a number) combines with a pattern's `where` or
`factor`, and the static methods accept a whole propensity only over a
species with a single instance. Use `mass_action` for rates per molecule.

## What the builder writes

The builder writes a `crn-craft.model/1.0` document. Every species is a
*type*, here one without attributes, and every reaction is a *rule* over
the types:

```python
document = model.to_dict()
print(document["schema"])
print(json.dumps(document["rules"][1]))
```

```text
crn-craft.model/1.0
{"name": "death", "reactants": [{"label": "x", "type": "A"}], "propensity": "death * A", "products": []}
```

A string rate becomes a `propensity`, the whole expression. A number or
`crn.mass_action` becomes a `rate`, the factor that multiplies the count of
each reactant. A species with attributes becomes a type with `attributes`,
`add_initial` an entry of `initial`, `add_observable` an entry of
`observables`, `crn.copy("x")` the product `{"type": "mtDNA", "copy": "x"}`
and `crn.new` a product with `attributes`.

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
wraps a document as a `Model`, with `fsp`, `analyze` and `to_gillespy2`. The
`add_*` calls work on it when the builder could have written it, which
includes documents with attributes, filtered observables and the rules of
any of these sections.

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

## A model document for three strains

A model can also be written as a document by hand, in Python as a
dictionary; model files from other tools and earlier releases look like
this. Three strains of cells share `room` places. The resident strain 0
gives rise to two variants, every strain dies at one rate, and strain 2
grows twice as fast into the free room. Each cell is an instance of the type
`Cell` with the attribute `strain`:

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

The parts of the document, and the builder calls that write them:

- **`initial`** lists instances: 200 cells with `strain` 0
  (`add_initial("Cell", 200, strain=0)`).
- **`observables`** reduce the instances of one type to a number. `where`
  selects instances (`variants`), and `weight` weighs them (`crowding`
  counts a cell of strain 2 twice); `add_observable` takes the same fields.
  Results have one column per observable, and rates and the `stop` read
  observables by name. This document has no count observable named after
  the type, which the builder would add.
- **A reactant** is a pattern: a `label`, a `type` and an optional `where`
  over the attributes (`crn.pattern("Cell", where="c.strain == 0",
  label="c")`). The `rate` is per matching cell, as a mass-action constant
  is (`crn.mass_action`).
- **A product** copies a reactant (`copy`, `crn.copy("c")`), gives all
  attributes (`attributes`, `crn.new`), or copies and overrides some.

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
([page S1d](../tutorials/simulation/S1d-mtdna-methods.md#the-fast-path)).
[Page 12a](12a-mtdna-model-files.md#the-fast-path) uses the fast path on the
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
method, which runs when `method` is absent, takes it. With the builder the
weighted reactant is `crn.pattern("Cell", factor="c.w", label="c")` and the
founder `crn.new("Cell", strain="new_id()", w="choice(0.2, 2, 1)")`.
[Page 12a](12a-mtdna-model-files.md#every-lineage-with-the-clone-model)
follows every mutant lineage of the mtDNA model in this way.

## Limits and versions

A `crn-craft.model/1.0` document takes up to 4,096 types, 16,384 parameters
and 65,536 rules ([limits](../reference/model-language.md#limits)). Model
files of earlier releases (`crn-craft.rules-model/0.x`) are still read
([Contracts and versions](../reference/contracts.md)).

See also: [Model language](../reference/model-language.md),
[Fig 1C/1D of a published study, built with the builder](../tutorials/simulation/S2c-fig1cd.md),
[The mtDNA model files](12a-mtdna-model-files.md).
