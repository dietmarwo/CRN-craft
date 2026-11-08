# CRN-craft

![CRN-craft: chemical reaction networks in native Rust with stochastic simulation](https://raw.githubusercontent.com/dietmarwo/CRN-craft/main/docs/images/crn-craft-banner.jpg)

Reproducible simulation and analysis of well-mixed chemical reaction networks,
with a native Rust core. Models, requests and results are versioned JSON
documents with strict input checks, explicit budgets and replayable seeds.

```bash
pip install crn-craft
```

Wheels for CPython 3.10 and newer (stable ABI) on Linux x86_64, Windows x64
and macOS (x86_64 and arm64); other platforms build from the source
distribution with a Rust toolchain.

## Example

Immigration–death (birth rate 4, death rate 0.2, starting from zero) has a
Poisson law with mean `20 (1 - exp(-0.2 t))`: 12.64 at t = 5 and 17.29 at
t = 10.

```python
import crn_craft as crn

model = crn.load_example("immigration-death")
request = crn.ssa_spec(
    model, end=10.0, save_times=[0.0, 5.0, 10.0], replicates=1000, seed=1
)
result = crn.simulate_ssa(model, request)
print(result["sampler"])
print(result["statistics"]["mean"], result["statistics"]["mean_standard_error"])
```

```text
xoshiro256pp-u53/direct_cr/inverse/fast-1
[[0.0], [12.618], [17.22500000000001]] [[0.0], [0.10959967193421745], [0.13107843665743252]]
```

Both means lie within one standard error of the exact values. The same seed
gives the same result on every run and for every thread count.

## What is included

- Exact stochastic simulation: `direct`, `cached_direct`, `direct_cr`,
  `sorting_direct` and `rssa_cr`, with inverse or Ziggurat waiting times.
  `fast` (the default) emulates JumpProcesses.jl; `fast=False` selects the
  checked policy, which rejects lost increments and channels below the 53-bit
  selection resolution.
- Prepared simulations and reusable parameter templates for optimization
  loops (`crn.prepare_ssa`, `crn.prepare_parameters`).
- Sparse finite state projection (`crn.fsp`), explicitly approximate Poisson
  leaping (`crn.tau`), network analysis and conservation checks.
- Strict SBML and BioNetGen import subsets (`crn-craft-import`) and the
  `crn-craft` command for archives, verification and export.
- `pip install "crn-craft[fit]"` adds the pinned SciPy, fcmaes-rust and
  SQP-craft versions used by the fitting examples.

The rules engine for dynamic, attributed species is in the wheel as
`crn_craft.rules`: prepared parallel ensembles, trackers, event logs and
callbacks, and the mtDNA methods; see the
[mtDNA user manual](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-user-manual.md).
`rules.Model` builds models from Python calls, in the style of rebop: a
numeric rate is a mass-action constant, a string the reaction's whole
propensity. `add_species` and `add_parameter` take one name and value or a
mapping such as `{"S": 999, "I": 1}`. Results are pandas DataFrames:
per-replicate records and rule counts, statistics by report time, trackers
and event logs. Kept paths are an xarray DataArray: `result.paths` has the
dimensions (time, variable, replicate). pandas and xarray are installed with
the package. `paths_frame()` returns the DataFrame of release 0.4.0, and
`to_xarray()` the paths as a Dataset. The same immigration–death model:

```python
from crn_craft import rules

model = rules.Model("immigration_death")
model.add_parameter("k_in", 4.0).add_species("X")
model.add_reaction("k_in", [], ["X"], name="immigration")  # string: whole propensity
model.add_reaction(0.2, ["X"], [], name="death")  # number: mass action, 0.2 X
result = rules.ensemble(
    model,
    end=10.0,
    report_times=[0.0, 5.0, 10.0],
    replicates=1000,
    seed=1,
    threads=2,
    method="static_direct",
    keep_paths="all",
).simulate()
print(result.statistics["mean"])
print(result.rule_events.sum().to_dict())
print(result.paths.sel(variable="X").T.to_pandas().head(3))  # cells x times
```

```text
           X
time
0.0    0.000
5.0   12.623
10.0  17.305
{'immigration': 39800, 'death': 22495}
time       0.0   5.0   10.0
replicate
0           0.0  14.0  19.0
1           0.0  19.0  20.0
2           0.0  17.0  17.0
```

The means lie within one standard error of the exact 12.64 and 17.29.
`method="static_direct"` draws Ziggurat waiting times and integrates count
trackers sparsely (rules-run/0.6). `method={"kind": "static_direct",
"waiting_time": "inverse", "sparse_trackers": False}` gives the numbers of
release 0.3.0, which printed 12.550 and 17.351 here.
Models may also declare continuous variables, ODEs between events that
rates read (for example an ATP pool), which `method="hybrid"` simulates at
declared tolerances.

New in 0.4.0, opt-in, with a Rust toolchain (1.85 or newer): `loop = rules.compile_loop(model)`
compiles a model's Static Direct event loop into a cached native program,
and `rules.ensemble(..., method="static_direct", compiled=loop)` runs it
with results byte-identical to the interpreter's, about 2.2× faster on
the three-class mtDNA model. The first build takes a minute or more; see
the [user manual](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-user-manual.md#optional-compiled-event-loop)
for requirements and limitations. Release 0.3.0 changed the rules results from NumPy
arrays to DataFrames; the
[changelog](https://github.com/dietmarwo/CRN-craft/blob/main/CHANGELOG.md)
lists what to migrate.
Its command-line tools (`rules-check`, `rules`) are in the native CLI
(`cargo install crn-craft-cli --locked`, crates
[crn-craft](https://crates.io/crates/crn-craft) and
[crn-craft-rules](https://crates.io/crates/crn-craft-rules)), which installs a
`crn-craft` command of the same name as the wheel's. The optional diffsol BDF
integrator with forward sensitivities and the job worker for
[optimization-jobs](https://github.com/dietmarwo/optimization-jobs) build from
the repository and are not part of the wheel.

## Documentation

- [Repository and README](https://github.com/dietmarwo/CRN-craft)
- [Tutorial](https://github.com/dietmarwo/CRN-craft/blob/main/docs/tutorial.md)
- [Alternative exact SSA methods and the `fast` switch](https://github.com/dietmarwo/CRN-craft/blob/main/docs/alternative-ssa.md)
- [Rules language](https://github.com/dietmarwo/CRN-craft/blob/main/docs/rules-language.md) and [mtDNA user manual](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-user-manual.md)
- [Optimization tutorials C1–C4](https://github.com/dietmarwo/CRN-craft/blob/main/tutorials/optimization/README.md)
- [Agent context](https://github.com/dietmarwo/CRN-craft/blob/main/docs/ai-context.md)
- [Changelog](https://github.com/dietmarwo/CRN-craft/blob/main/CHANGELOG.md)

A successful numerical run is not biological validation. MIT licensed; see
[THIRD_PARTY_NOTICES.md](https://github.com/dietmarwo/CRN-craft/blob/main/THIRD_PARTY_NOTICES.md)
for the Ziggurat tables adapted from Julia.
