# CRN-craft

<!-- check: standalone -->

CRN-craft simulates chemical reaction networks stochastically, from Python,
with a native engine.

```sh
pip install crn-craft
```

This example builds a model, runs 1,000 trajectories and prints their means:

```python
import crn_craft as crn

model = crn.Model("immigration_death")
model.add_species("A", 0)
model.add_reaction(4.0, [], ["A"], name="birth")  # molecules arrive at rate 4
model.add_reaction(0.2, ["A"], [], name="death")  # each one decays at rate 0.2
run = crn.ensemble(model, end=10.0, report_times=[0.0, 5.0, 10.0],
                   replicates=1000, seed=1, threads=2, method="static_direct")
result = run.simulate()
print(result.statistics["mean"])  # exact: 20 * (1 - exp(-0.2 t)) = 0, 12.64, 17.29
```

```text
           A
time        
0.0    0.000
5.0   12.623
10.0  17.305
```

The 1,000 trajectories ran on two threads in one call, and the same seed
gives the same numbers on any number of threads.

## Why CRN-craft

- **Many trajectories in one call.** An ensemble runs on native threads, and
  the same seed gives the same result on any number of threads
  ([chapter 3](guide/03-running.md#many-trajectories-in-one-call)).
- **Compiled models from Python.** `crn.compile_loop(model)` compiles a
  model's event loop, and the results stay byte-identical. On 12 November
  2026 the compiled mtDNA model ran 2.1 times as fast as the interpreted one
  ([chapter 5](guide/05-speed.md#compiled-models)).
- **Exact methods and exact distributions.** Static Direct with five
  selectors, next reaction and finite state projection; approximate methods
  are labeled as approximate ([chapter 7](guide/07-methods.md)).
- **More than reaction networks.** Continuous variables
  ([chapter 9](guide/09-continuous-variables.md)) and species with attributes,
  such as clones, sizes and lineages ([chapter 14](guide/14-model-files.md)).
- **Runs that can be checked.** Strict input checks, explicit budgets, and
  sealed archives with a verifier ([chapter 10](guide/10-command-line.md)).

Coming from GillesPy2 or rebop? See [how CRN-craft compares](comparison.md).

## Where to start

- **New here:** read [chapters 1 to 4](guide/01-install.md) in order.
- **A question:** the [task index](howto.md) links fifty tasks to their answers.
- **mtDNA models:** the [case study](guide/12-mtdna.md).
- **Code of an earlier release:** [Migrating from earlier releases](guide/migrating.md).

## What is included

| Part | What it does |
| --- | --- |
| Model input | the builder `crn.Model`; JSON model files; strict SBML and BioNetGen subsets ([chapter 11](guide/11-analysis-import-export.md)) |
| Exact methods | Static Direct with five selectors, next reaction, the tree method for species that appear and disappear |
| Approximate methods | tau leaping, binomial splitting, `hybrid` for continuous variables |
| Distributions | finite state projection with separate error terms |
| Results | pandas and xarray; mean, variance, standard error, quantiles; trackers, event logs and callbacks; reports in markdown or JSON |
| Reproducibility | one random stream per replicate, results independent of the thread count, model and request named by SHA-256 |
| Interfaces | Python, the `crn-craft` command line, Rust; R through reticulate |
| Analysis | conservation laws and structure, GillesPy2 export; ODE right-hand sides and Jacobians in Rust |

## The rest of the site

- [Reference](reference/model-language.md): the model language, methods, the Python API, contracts, results, archives and diagnostics.
- [Developer](developer/architecture.md): the architecture, the Rust API, extending and checking the project.
- [Evidence](evidence/technical-report.md): the technical report, measurements and comparisons.
- [History](history/index.md): plans, decisions and campaign logs.

Release 0.7.0 merged two earlier engines into one; the
[changelog](../CHANGELOG.md)
says what changed and how to migrate. CRN-craft is one of 18 related
repositories; the
[repository overview](https://github.com/dietmarwo/fcmaes-rust/blob/main/docs/repositories.md)
shows how they connect.
