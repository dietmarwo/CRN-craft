# The mtDNA model files: clones, the fast path and stored models

<!-- fixtures: models/rules -->

The case study built its model in Python. The repository also has the mtDNA
models as JSON files, in which every molecule is an instance of the type
`mtDNA` with attributes ([chapter 14](14-model-files.md)). Three things need
those files: following every mutant lineage, the fast path `mtdna_direct`,
and the stored form of the ATP model. To run this page, copy the folder
`models/rules/` of the repository into your working folder. The setup block
is the model of the [index page](12-mtdna.md#the-model), for comparison with
the files. The blocks of this page continue from it:

```python
import math

import numpy as np
import crn_craft as crn

END = 25_920.0  # three years (1,080 days), in hours


def three_class(K=1000, w0=1000, m0=0, mu=3e-5, sm=1.5, half_life_days=10, takeover=None):
    """The three-class mtDNA model: wild type W, neutral M, advantaged A."""
    model = crn.Model("mtdna_three_class")
    model.add_parameter({"K": K, "delta": math.log(2) / (24 * half_life_days),
                         "sw": 1.0, "sm": sm, "mu": mu, "fracAdv": 0.2})
    model.add_species({"W": w0, "M": m0, "A": 0})
    free = "max(0, K - W - M - A) / (W + sm * (M + A))"
    model.add_reaction(f"{free} * sw * (1 - mu) * W", ["W"], ["W", "W"], "replication_wt")
    model.add_reaction(f"{free} * sw * mu * (1 - fracAdv) * W", ["W"], ["W", "M"], "mutation_mt")
    model.add_reaction(f"{free} * sw * mu * fracAdv * W", ["W"], ["W", "A"], "mutation_mtAdv")
    model.add_reaction(f"{free} * sw * M", ["M"], ["M", "M"], "replication_mt")
    model.add_reaction(f"{free} * sm * A", ["A"], ["A", "A"], "replication_mtAdv")
    for label, species in (("wt", "W"), ("mt", "M"), ("mtAdv", "A")):
        model.add_reaction(f"delta * {species}", [species], [], f"degradation_{label}")
    if takeover is not None:
        model.set_stop(f"W + M + A > 0 && M + A >= {takeover} * (W + M + A)")
    return model
```

## The five model files

Applies: [chapter 14, When species are not a fixed list](14-model-files.md#when-species-are-not-a-fixed-list).

| File in `models/rules/` | What it models | Methods |
| --- | --- | --- |
| `mtdna-three-class.rules.json` | the model of the index page, as three classes of one type | every generic method, and the fast path `mtdna_direct` |
| `mtdna-three-class-cme.rules.json` | the same with `K = 20` and fast rates, for exact checks against the master equation | as above |
| `mtdna-three-class-atp-continuous.rules.json` | the ATP model of [page 12c](12c-mtdna-atp.md), with ATP as a continuous variable | `hybrid` |
| `mtdna-three-class-atp.rules.json` | the same in the stored form, with a fixed `atp` block | `thinning` |
| `mtdna-fixed.rules.json` | every new mutant is its own clone with an identity | the tree method |

Each molecule has three attributes: `clone` (0 for wild type, 1 and 2 for
the two mutant classes, or an identity in the clone model), `s` (its
replication weight) and `wt` (1 for wild type, 0 for mutants). The files
declare observables that the builder model writes as sums: `N` is
`W + M + A`, `wild_type` is `W`, `mutants` is `M + A`, `neutral_mutants` and
`advantaged_mutants` are `M` and `A`, and `D` is the competition
denominator. The parameters are those of the index page, with `w0` as a
parameter. The files have no stop; the blocks below add the takeover stop
`N > 0 && mutants >= 0.9 * N`.

## The fast path

Applies: [chapter 14, Run it](14-model-files.md#run-it), and [page 12d](12d-mtdna-methods.md).

`mtdna_direct` is an exact method fused for one structure: three classes in
a fixed population. At preparation it certifies the model:

1. three classes and eight channels, with the roles of degradation,
   replication and mutation read from the channels;
2. one constant degradation rate shared by all classes;
3. every replication and mutation rate is `max(0, K - N) / D` times
   constant factors, with one integer `K` of at most 1,000,000;
4. every observable is a weighted sum over the three classes;
5. no stop, or the guarded fraction `N > 0 && mutants >= f * N`;
6. every positive hazard stays finite and resolvable in every reachable
   state.

The builder model writes whole-rule propensities, which the certificate
does not read, so the fast path needs the file. The same 64 cells under the
fast path and under Static Direct, and the builder model beside them:

```python
import json
import time

document = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
document["stop"] = "N > 0 && mutants >= 0.9 * N"
arms = (("file", document, "mtdna_direct"), ("file", document, "static_direct"),
        ("builder", three_class(takeover=0.9), "static_direct"))
for label, model, method in arms:
    run = crn.ensemble(model, end=END, report_times=[0.0, 12_960.0, END], replicates=64,
                       seed=5, threads=2, method=method)
    start = time.perf_counter()
    result = run.simulate()
    seconds = time.perf_counter() - start
    stopped = (result.records["termination"] == "stopped").sum()
    print(f"{label:8s} {method:14s} {seconds:5.2f} s  {result.metadata['events']:>10,} events  "
          f"taken over {stopped:2d}  {result.metadata['sampler']}")
```

```text
file     mtdna_direct    0.17 s   8,969,870 events  taken over  8  chacha8-u53-ziggurat/rules-mtdna-direct/1
file     static_direct   0.55 s   8,861,769 events  taken over 10  chacha8-u53-ziggurat/rules-static-direct/1
builder  static_direct   1.66 s   9,231,158 events  taken over  7  chacha8-u53-ziggurat/rules-static-direct/1
```

The fast path samples the same process as Static Direct from other random
numbers, so the results agree in distribution, not cell by cell. On 12
November 2026, on a two-CPU host, it took 2.8 times less time than Static
Direct on the file. The builder model under Static Direct took 9.9 times as
long as the fast path, and its compiled loop 1.4 times as long
([measurements](https://github.com/dietmarwo/CRN-craft/blob/main/validation/merge/RESULTS.md)).

A model the certificate does not accept still runs: as Static Direct with
inverse waits, with a `CRN-R051` warning that names the failed check. Read
`metadata["sampler"]` to see what ran:

```python
odd = dict(document, stop="mutants > 250")  # not the guarded fraction
for model in (odd, three_class(takeover=0.9)):
    result = crn.ensemble(model, end=100.0, report_times=[0.0, 100.0], replicates=2,
                          method="mtdna_direct").simulate()
    print(result.metadata["sampler"])
    print([d["message"] for d in result.metadata["diagnostics"]])
```

```text
chacha8-u53/rules-static-direct/1
['mtdna_direct ran Static Direct: stop is neither absent nor the guarded mutant fraction']
chacha8-u53/rules-static-direct/1
['mtdna_direct ran Static Direct: the model has whole-rule propensities']
```

To record another threshold at the full speed of the fast path, keep the
guarded stop and add a first-passage tracker
([page 12b](12b-mtdna-events.md#class-appearances-and-extinctions)).
[Methods and their options](../reference/methods.md) has the full
certificate.

## Every lineage with the clone model

Applies: [chapter 14, New identities and random attributes](14-model-files.md#new-identities-and-random-attributes).

In `mtdna-fixed.rules.json` every new mutant gets its own `clone` identity
(`new_id()`) and a drawn weight (`choice(fracAdv, sm, sw)`). The set of
instances changes during a run, so only the tree method runs it. The
observable `mutant_clones` counts distinct mutant clones, and a complete
event log gives every lineage's birth, death and largest size. A small cell
(`K = 100`) with a high mutation probability keeps every event in the log:

```python
import json

model = json.loads(crn.load("models/rules/mtdna-fixed.rules.json"))
model["parameters"].update({"K": 100, "w0": 100, "mu": 0.01, "delta": 0.05})
result = crn.ensemble(
    model, end=400.0, report_times=[0.0, 100.0, 200.0, 400.0], replicates=20,
    seed=6, threads=2, keep_instances=3, trackers=[{"observable": "mutant_clones"}],
    event_log={"max_entries": 5000, "rules": []}, max_output_bytes=1 << 30,
).simulate()
print("sampler:", result.metadata["sampler"])
print("mean distinct mutant clones:", result.statistics["mean"]["mutant_clones"].to_numpy())
print("peak instances, max over cells:", result.records["peak_instances"].max())
print("instances created per cell:", result.records["created_instances"].mean())
print("log entries omitted:", sum(o["omitted_log_entries"] for o in result.observations))

# Lineages from the full log: replay each cell's net changes per clone id.
born, died, peak = {}, {}, {}
for observation in result.observations:
    counts = {}
    for entry in observation["log"]:
        change = {}
        for sign, records in ((-1, entry["reactants"]), (1, entry["products"])):
            for record in records:
                clone, _, wt = record["attributes"]
                if wt == 0:  # mutant clones only
                    key = (observation["id"], clone)
                    change[key] = change.get(key, 0) + sign * record["count"]
        for key, delta in change.items():
            counts[key] = counts.get(key, 0) + delta
            born.setdefault(key, entry["time"])
            peak[key] = max(peak.get(key, 0), counts[key])
            if counts[key] == 0:
                died[key] = entry["time"]
life = np.array([died[k] - born[k] for k in born if k in died])
print(f"{len(born)} clones born, {len(life)} lost within the horizon")
print("median lifetime of lost clones (h):", round(float(np.median(life)), 2))
print("largest clone ever (molecules):", max(peak.values()))
```

```text
sampler: chacha8-u53/rules-direct/2
mean distinct mutant clones: [0.   1.9  2.6  1.85]
peak instances, max over cells: 7
instances created per cell: 10.1
log entries omitted: 0
202 clones born, 165 lost within the horizon
median lifetime of lost clones (h): 18.43
largest clone ever (molecules): 100
```

Most lineages are lost within a day of their birth; a few grow. For
`K = 1000` over three years, a cell has about 150,000 events, too many to
log whole. Follow lineages there with `mutant_clones`, `peak_instances` and
`created_instances`, and log only the mutation and degradation rules of a
few cells.

## The stored ATP block and thinning

Applies: [chapter 7, Other methods](07-methods.md#other-methods), and [page 12c](12c-mtdna-atp.md).

`mtdna-three-class-atp.rules.json` writes the ATP model of page 12c with a
fixed `atp` block, a form of an earlier release. `method="thinning"`
samples it exactly, against a bound on each interval. Stored files and
requests of this form keep working and replay byte for byte. For new
models, use a continuous variable (page 12c):

```python
import json

data = crn.load("models/rules/mtdna-three-class-atp.rules.json")
print("ATP block:", json.loads(data)["atp"])
result = crn.ensemble(data, end=END, report_times=[0.0, END], replicates=64, seed=7,
                      threads=2, method="thinning").simulate()
print("sampler:", result.metadata["sampler"])
print("stopped (ATP exhausted):", (result.records["termination"] == "stopped").sum())
print(result.statistics["mean"][["N", "mutants", "ATP"]].round(2))
```

```text
ATP block: {'wild_type': 'wild_type', 'total': 'N', 'initial': 450.0, 'production': 0.1, 'decay': 0.2, 'maintenance': 0.01, 'scale': 1.3034153747670567, 'replication_rules': ['replication_wt', 'mutation_mt', 'mutation_mtAdv', 'replication_mt', 'replication_mtAdv']}
sampler: chacha8-u53/rules-atp-thinning/1
stopped (ATP exhausted): 8
               N  mutants     ATP
time                             
0.0      1000.00     0.00  450.00
25920.0  5919.92  4587.14  394.85
```

The final means of this run and of page 12c differ. Mutant expansion is
rare and lets the population grow, so the final population is
heavy-tailed, and 64 cells are too few to compare methods.
[`validation/check_hybrid.py`](https://github.com/dietmarwo/CRN-craft/blob/main/validation/check_hybrid.py)
compares `hybrid` with thinning on 400 cells each, with declared tests.

See also: [Model files and species with attributes](14-model-files.md),
[Methods and their options](../reference/methods.md).
