# Migrating from earlier releases

This page is for code, models and results of earlier releases. A new user
can skip it. The setup block builds the strain model of
[chapter 7](07-methods.md#static-direct-and-its-selectors). The blocks of
this page continue from it:

```python
import crn_craft as crn


def strains(count, builder=crn.Model):
    """Strains S0, S1, ... share room K; strain i gives rise to strain i + 1."""
    names = [f"S{i}" for i in range(count)]
    model = builder(f"strains_{count}")
    model.add_parameter({"K": 500.0, "g": 0.02, "mu": 0.001})
    model.add_species({name: 200 if i == 0 else 0 for i, name in enumerate(names)})
    free = f"max(0, K - ({' + '.join(names)})) / K"
    for i, name in enumerate(names):
        model.add_reaction(f"g * {name} * {free}", [name], [name, name], f"grow_{i}")
        model.add_reaction(0.01, [name], [], f"die_{i}")
        if i + 1 < count:
            model.add_reaction(f"mu * {name}", [name], [name, names[i + 1]], f"mutate_{i}")
    return model
```

## One engine since 0.7.0

Releases up to 0.6.0 had two engines: a network engine for fixed reaction
networks, at the top level of the Python package, and a rules engine for
models whose species carry attributes, in `crn_craft.rules`. Release 0.7.0
merged them into one engine, which runs every model of both. Network models
of 0.6.0 convert with `crn-craft convert`, and code written for
`crn_craft.rules` keeps running: that module keeps the API and the defaults
of release 0.6.0 under its old name. The
[changelog](https://github.com/dietmarwo/CRN-craft/blob/main/CHANGELOG.md#migration)
has the full account.

## From the 0.6.0 module to the top level

Code written for 0.6.0 with `from crn_craft import rules` runs unchanged:
the same names, defaults, documents, hashes and numbers. Code that moves to
`import crn_craft as crn` changes as follows:

| `crn.rules` (0.6.0) | `crn` (0.7.0) |
| --- | --- |
| `rules.Model` writes `rules-model/0.3` | `crn.Model` writes `crn-craft.model/1.0`: larger limits, three or more reactant molecules |
| `rules.ensemble` writes the lowest `rules-run` version, and Static Direct scans | `crn.ensemble` writes `crn-craft.run/1.0`; Static Direct uses the selector `auto`, so its paths can differ for models with many reactions (`selector="linear"` keeps them) |
| `statistics` has `mean`, `sample_variance` and `mean_standard_error` | also `median` and `quantiles` by default |
| `result.metadata["schema"]` is `crn-craft.rules-result/0.x` | `crn-craft.result/1.0`, with new hashes |
| `crn_craft.rules.describe()` | `crn.describe()`: `crn-craft.describe/1.0`; the 0.6.0 fields are under `legacy.rules` |

The old module always scans, whatever the size of the model, so adding
`quantiles=` to a 0.6.0 script never changes its paths. Its models keep the
limits of `rules-model`, which a hundred strains exceed:

```python
settings = dict(end=300.0, report_times=[0.0, 300.0], replicates=200, seed=1, threads=2)
old = crn.rules.ensemble(strains(20), method="static_direct", **settings).simulate()
print("crn.rules:", old.metadata["sampler"], sorted(old.statistics))
new = crn.ensemble(strains(20), method="static_direct", **settings).simulate()
print("crn:      ", new.metadata["sampler"], sorted(new.statistics))

big = strains(100, crn.rules.Model)
print("rules.Model schema:", big.to_dict()["schema"])
try:
    crn.rules.ensemble(big, method="static_direct", **settings)
except crn.CrnError as error:
    print("refused:", error.diagnostics[0]["code"], error.diagnostics[0]["path"])
result = crn.ensemble(strains(100), method="static_direct", **settings).simulate()
print("Model:", result.metadata["sampler"], len(result.rule_names), "reactions")
```

```text
crn.rules: chacha8-u53-ziggurat/rules-static-direct/1 ['mean', 'mean_standard_error', 'sample_variance', 'samples']
crn:       chacha8-u53-ziggurat/rules-static-direct/1 ['mean', 'mean_standard_error', 'median', 'quantiles', 'sample_variance', 'samples']
rules.Model schema: crn-craft.rules-model/0.3
refused: CRN-R015 /types
Model: chacha8-u53-ziggurat/rules-static-direct-cr/1 299 reactions
```

On twenty strains `auto` also picks the linear scan, so the two modules
give the same paths; the top level adds `median` and `quantiles`. On a
hundred strains the old module refuses the model, and the top level runs it
with composition–rejection. To lift the limits of an old model file, change
its `schema` string to
`crn-craft.model/1.0`: the two mean the same thing.

## From the network engine of 0.6.0, in Python

The network functions of the top level were removed. Each has a successor:

| 0.6.0 | 0.7.0 |
| --- | --- |
| `crn.load_example(name)` (network bytes) | `crn.load_example(name)` returns a `crn.Model` of the converted model |
| `crn.run_spec(...)` + `crn.simulate(model, run)` | `crn.ensemble(model, end=..., report_times=..., replicates=..., method="static_direct").simulate()` |
| `crn.prepare(model, run, cached=...)` | `crn.prepare(model, request)` with a rules or `run/1.0` request |
| `crn.ssa_spec`, `crn.simulate_ssa`, `crn.prepare_ssa` | `crn.ensemble(..., method="static_direct", selector=...)` or `method="next_reaction"` |
| `aggregator: "direct"`, `"cached_direct"`, `"direct_cr"`, `"sorting_direct"`, `"rssa_cr"` | `selector`: `"linear"`, `"linear"`, `"composition_rejection"`, `"sorting"`, `"rssa_cr"` |
| `aggregator: "next_reaction"` | `method="next_reaction"` |
| `fast: true` or `false` | no such switch |
| `crn.prepare_parameters(model, request)`, `ParameterTemplate.bind(values)` | `crn.ensemble(...)`, then `run.bind_parameters({name: value})` ([chapter 5](05-speed.md#parameter-sweeps)) |
| `template.parameter_names` | the model's `parameters` |
| `crn.mass_action_kernel_source(model)` | `crn.loop_source(model)`, `crn.compile_loop(model)` ([chapter 5](05-speed.md#compiled-models)) |
| `crn.fsp(model, run)` | `model.fsp(end=..., report_times=..., bounds={...})` |
| `crn.tau(model, request)` | `crn.ensemble(..., method={"kind": "tau", "epsilon": ...})` |
| `crn.analyze(model)` | `crn.Model(...).analyze()` or `crn.load_example(name).analyze()` (`analysis/1.0`) |
| `crn.check(model, run)` (`check/0.1`) | `crn.check(model, request)` (`check/1.0`) |
| `crn.model_bytes(document)` | `crn.rules.model_bytes(document)` |

## From the network engine, on the command line

| 0.6.0 | 0.7.0 |
| --- | --- |
| `crn-craft simulate MODEL RUN DIR`, `simulate-cached` | `crn-craft convert MODEL > model.json`, then `crn-craft run model.json request.json --archive DIR` |
| `crn-craft ssa MODEL REQUEST` | `crn-craft run` with `"method": {"kind": "static_direct", "selector": ...}` or `next_reaction` |
| `crn-craft tau MODEL REQUEST` | `crn-craft run` with `"method": {"kind": "tau", "epsilon": ...}` |
| `crn-craft fsp MODEL RUN` | `crn-craft run MODEL REQUEST` with the method `fsp` |
| `crn-craft check MODEL [RUN]` (`check/0.1`) | `crn-craft check MODEL [REQUEST]` (`check/1.0`) |
| `crn-craft verify DIR` of an `archive/0.1` | rerun with `crn-craft run ... --archive DIR`; `verify` reads `archive/1.0` |
| `crn-craft rules`, `rules-check` | `crn-craft run`, `check` (the old names are kept as aliases) |

[Chapter 11](11-analysis-import-export.md#convert-a-network-model-of-release-060)
shows `crn-craft convert` and `crn.convert`.

## What replays and what does not

- **Rules documents and `crn_craft.rules` code of 0.6.0** run unchanged
  with the same results, byte for byte.
- **Network-engine results of 0.6.0 do not replay byte for byte.** The
  network loops and their random stream (xoshiro256++) are gone. The one
  engine samples the same laws: a law gate checked it over 4,814
  comparisons on twelve models, including references from finite state
  projection, with no failure, for every selector.
- **The top level and the old module differ in paths** for models on which
  `auto` picks a selector other than the linear scan, and every result of
  the top level has new hashes, because its request and result documents
  are new versions.

## Contracts

`model/0.1` → `model/1.0` (by `crn-craft convert`); `run/0.1`,
`ssa-run/0.1`–`0.4`, `tau-request/0.1` → `run/1.0`; `result/0.1` and `0.2`,
`ssa-result/0.1`–`0.4`, `tau-result/0.1` → `result/1.0`; `fsp-result/0.1`
→ `fsp-result/1.0`; `check/0.1` → `check/1.0`; `analysis/0.1` and `0.2` →
`analysis/1.0`; `archive/0.1`, `archive-check/0.1` → `archive/1.0`,
`archive-check/1.0`; `describe/0.11` → `describe/1.0`. The `rules-*`
contracts are unchanged, and `model/1.0` means what `rules-model/0.4` means.
[Contracts and versions](../reference/contracts.md) lists every version
that is still read.

## Older releases

Python packages: pandas is required since 0.3.0 and xarray since 0.5.0. The
extra `crn-craft[xarray]` still installs.

Code written for release 0.5.0:

| 0.5.0 | Now |
| --- | --- |
| `result.paths` with dimensions (`time`, `variable`, `replicate`) | dimensions (`replicate`, `time`, `variable`); `result.paths.transpose("time", "variable", "replicate")` gives the old order |
| `result.paths.values[t, v, r]` | `result.paths.values[r, t, v]`; code that selects by name (`sel`, `mean("replicate")`) is unchanged |

R code written for 0.5.0 that indexes `paths` by position must use the new
order; code that reads the dimension names, as in
[chapter 15](15-from-r.md#example-sir-model), works with both.

Code written for release 0.4.0:

| 0.4.0 | Now |
| --- | --- |
| `result.paths`, a DataFrame | `result.paths_frame()`; `result.paths` is a DataArray |
| `result.paths.loc[3]` | `result.paths.sel(replicate=3)`, or `result.paths_frame().loc[3]` |
| `result.paths["mutants"].unstack("time")` | `result.paths.sel(variable="mutants").to_pandas()` |
| `result.event_columns()`, `result.event_columns(id)`, DataFrames | `result.event_frame()`, `result.event_frame(id)`; `event_columns` returns DataArrays |
| `pip install "crn-craft[xarray]"` | `pip install crn-craft`: xarray is required; the extra still installs |
| `Model.add_continuous(...)` | `Model.add_ode(...)` with the same arguments; `add_continuous` and the alias `add_ODE` still work, and `add_continuous` emits a `DeprecationWarning` |

Code written for the NumPy results of release 0.2.0:

| 0.2.0 | Now |
| --- | --- |
| `records["termination"] == crn.rules.TERMINATIONS.index("stopped")` | `records["termination"] == "stopped"`; the old form is silently all `False` |
| `records["id"]` | `records.index` |
| `records["rule_events"]` | `result.rule_events` |
| `records[name]` (an array) | `records[name].to_numpy()`; integer termination codes are `records["termination"].cat.codes` |
| `statistics["mean"][:, names.index("N")]` | `statistics["mean"]["N"]` |
| `statistics["mean"]` (an array) | `statistics["mean"].to_numpy()` |
| `metadata["paths"]` rows | `result.paths_frame()`, or `result.paths_frame().to_numpy()` |
| `tracker_columns(i)["id"]` and fields | `tracker_columns(i).index` and its columns |
| `event_columns(id)["observables"]` | `event_columns(id).to_numpy()` |
| `keep_paths=True` (one cell) | `keep_paths=1`; `True` now keeps every cell |

## Rust

The Rust changes of 0.7.0 (the crate path `crn_craft::rules`, new struct
fields, `Kinetics::from_bytes`, parameter binding) are in
[Rust API](../developer/rust-api.md).

See also: [Contracts and versions](../reference/contracts.md),
[Convert a network model](11-analysis-import-export.md#convert-a-network-model-of-release-060).
