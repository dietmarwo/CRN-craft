# Large studies: sweeps, budgets and the command line

A study runs the model many times: for many parameter values, for long
cells, or as a batch of files on a cluster. The setup block is the model of
the [index page](12-mtdna.md#the-model). The blocks of this page continue
from it:

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

## Parameter sweeps

Applies: [chapter 5, Parameter sweeps](05-speed.md#parameter-sweeps).

How does the takeover fraction depend on the mutation probability `mu` and
the advantage `sm`? One prepared run serves as the template, and
`bind_parameters` gives each candidate its values:

```python
import time

template = crn.ensemble(three_class(takeover=0.9), end=END, report_times=[0.0, END],
                        replicates=200, seed=8, threads=2, method="mtdna_direct")
start = time.perf_counter()
print("takeover fraction after three years")
print("   mu \\ sm " + "".join(f"{sm:8.2f}" for sm in (1.0, 1.25, 1.5, 2.0)))
for mu in (1e-5, 3e-5, 1e-4):
    row = []
    for sm in (1.0, 1.25, 1.5, 2.0):
        result = template.bind_parameters({"mu": mu, "sm": sm}).simulate()
        row.append((result.records["termination"] == "stopped").mean())
    print(f"{mu:9.0e}  " + "".join(f"{v:8.3f}" for v in row))
print(f"12 candidates in {time.perf_counter() - start:.1f} s")
```

```text
takeover fraction after three years
   mu \ sm     1.00    1.25    1.50    2.00
    1e-05     0.000   0.025   0.040   0.075
    3e-05     0.000   0.040   0.110   0.190
    1e-04     0.000   0.155   0.325   0.510
12 candidates in 6.1 s
```

With equal weights (`sm = 1`), drift alone is too slow for a takeover
within three years at `K = 1000`. The takeover fraction grows roughly in
proportion to `mu` and steeply with `sm`. The same seed gives every
candidate the same random streams
([common random numbers](05-speed.md#common-random-numbers)), so the
differences between the cells of the table come from the parameters, not
from different seeds. With 200 cells, a fraction has a standard error of up
to 0.035; use more cells for a table to publish. Every candidate runs on
the fast path. Under Static Direct the sweep takes about 3.6 times as long
([page 12d](12d-mtdna-methods.md#the-fast-path)). A compiled loop does not
help here: it serves one parameter set. [Chapter 13](13-fitting.md) turns
the sweep into a fit.

## Budgets for long cells

Applies: [chapter 3, Budgets](03-running.md#budgets).

The default event budget is 1,000,000 events per cell. A cell at
`K = 1000` has about 150,000 events in three years, and two to three times
more with a doubled `K` or horizon. Raise `max_events` before a large study
rather than average over truncated cells: an ensemble in which a cell ended
on a budget withholds its statistics. The
[performance analysis](../evidence/mtdna-performance-analysis.md#8-sizing-a-study)
sizes studies by events, cells and threads.

## From the command line

Applies: [chapter 10, A request file](10-command-line.md#a-request-file), and [Run, archive and verify](10-command-line.md#run-archive-and-verify).

On a cluster, a study is often a set of files: a model file, a request per
job and an archive per result. This block writes the model with `to_json()`
and a request that binds it by path and SHA-256:

```python
import hashlib
import json
import pathlib

data = three_class(takeover=0.9).to_json().encode()
pathlib.Path("takeover.model.json").write_bytes(data)
request = {
    "schema": "crn-craft.run/1.0",
    "model": {"path": "takeover.model.json", "sha256": hashlib.sha256(data).hexdigest()},
    "method": {"kind": "mtdna_direct"},
    "time": {"end": END, "report_times": [0.0, 8640.0, 17_280.0, END]},
    "ensemble": {"replicates": 200, "root_seed": 11},
    "outputs": {"statistics": ["mean", "variance"], "quantiles": [0.5],
                "keep_paths": 0, "keep_instances": 0},
    "budget": {"threads": 2, "max_events_per_replicate": 10_000_000,
               "max_sampling_steps_per_replicate": 100_000_000,
               "max_instances": 10, "max_output_bytes": 64 << 20},
    "trackers": [{"observable": "A"}],
}
pathlib.Path("takeover.run.json").write_text(json.dumps(request, indent=2) + "\n")
print("model sha256:", request["model"]["sha256"])
```

```text
model sha256: 6bfa75e9f57beb6c5489e2d81c889c1100e9b10091d246a9ba1da037ee3376cc
```

The `crn-craft` command of the Python package, or the native tool, checks
the two files, runs them into a sealed archive, and verifies the archive:

```bash
crn-craft check takeover.model.json takeover.run.json > check.json
crn-craft run takeover.model.json takeover.run.json --archive takeover-archive > run.json
crn-craft verify takeover-archive > verify.json
```

Each command writes a JSON document. With `--archive`, `run` writes the
result into the archive and prints a short check of it:

```python
def read(name):
    return json.loads(pathlib.Path(name).read_text())


check, archived, verify = read("check.json"), read("run.json"), read("verify.json")
result = read("takeover-archive/result.json")
print("check:", check["valid"], " archive:", archived["archive_complete"],
      " verify:", verify["valid"])
print(result["termination"], result["sampler"])
print("cells:", result["completed_replicates"], "stopped:", result["stopped_replicates"])
print("events:", result["events"])
```

```text
check: True  archive: True  verify: True
completed chacha8-u53-ziggurat/rules-mtdna-direct/1
cells: 200 stopped: 20
events: 28757292
```

A result names its model and request by hash, and the archive holds the
exact bytes of both, so a job can be checked and replayed later on any host
([chapter 10, Replay](10-command-line.md#replay)).

## Checks behind these results

Applies: [chapter 8, What a run does not show](08-reproducibility.md#what-a-run-does-not-show).

- Every exact method, the fast path included, agrees with an exact solution
  of the master equation for the model at `K = 20`
  ([`validation/check_mtdna.py`](https://github.com/dietmarwo/CRN-craft/blob/main/validation/check_mtdna.py)).
- The selectors and the quantiles are checked against exact laws and NumPy
  ([`validation/merge/RESULTS.md`](https://github.com/dietmarwo/CRN-craft/blob/main/validation/merge/RESULTS.md)).
- Trackers equal their recomputation from complete event logs, and results
  with and without observers are byte-identical.
- On 10,000 cells, the fast path agreed with an independent Java
  implementation of the model
  ([`validation/mtdna/java-gap.json`](https://github.com/dietmarwo/CRN-craft/blob/main/validation/mtdna/java-gap.json)).

A numerical check confirms that the simulator solves the model. It does not
validate the model against biology.

See also: [Fitting parameters to data](13-fitting.md),
[The command line, request files and archives](10-command-line.md).
