# CRN-craft

<!-- check: standalone -->

![CRN-craft: chemical reaction networks in native Rust with stochastic simulation. Models, budgets, evidence: define, simulate, verify.](docs/images/crn-craft-banner.jpg)

> **One of 18 related repositories.** The [repository overview](https://github.com/dietmarwo/fcmaes-rust/blob/main/docs/repositories.md) shows what each one does and how they connect.

CRN-craft simulates chemical reaction networks stochastically, from Python,
with a native engine.

```sh
pip install crn-craft
```

## Example

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
  the same seed gives the same result on any number of threads.
- **Compiled models from Python.** `crn.compile_loop(model)` compiles a
  model's event loop, and the results stay byte-identical. On 12 November
  2026 the compiled mtDNA model ran 2.1 times as fast as the interpreted one.
- **Exact methods and exact distributions.** Static Direct with five
  selectors, next reaction and finite state projection; approximate methods
  are labeled as approximate.
- **More than reaction networks.** Species with attributes, such as clones,
  sizes and lineages, and continuous variables that rates read.
- **Runs that can be checked.** Strict input checks, explicit budgets, and
  sealed archives with a verifier.

Coming from GillesPy2 or rebop? See [how CRN-craft compares](docs/comparison.md).

## Documentation

- [Documentation site](https://dietmarwo.github.io/CRN-craft/): everything below, with search
- [User guide](docs/guide/01-install.md): read chapters 1 to 4 first
- [How do I ...?](docs/howto.md): fifty tasks, each with a link to the answer
- [Case study: mtDNA mutants](docs/guide/12-mtdna.md)
- [Reference](docs/reference/model-language.md): model language, methods, Python API, contracts, diagnostics
- [Changelog](CHANGELOG.md)

## Other ways to use it

The wheel also installs the `crn-craft` command line, which runs model and
request files and seals results in archives. The engine is the Rust crate
`crn-craft` (`cargo add crn-craft`), and R uses the Python package through
reticulate.

## Status and license

Version 0.7.0. Release 0.7.0 merged two earlier engines into one; the
changelog says what changed and how to migrate. A successful numerical run
is not biological validation. MIT licensed; see `LICENSE` and
`THIRD_PARTY_NOTICES.md`.
