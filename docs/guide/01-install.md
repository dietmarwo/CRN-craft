# 1. Install

CRN-craft simulates chemical reaction networks stochastically, from Python,
with a native engine. This chapter installs it and says how the guide is
organized.

## Install the Python package

```sh
python -m pip install crn-craft
```

The package contains the engine, so nothing has to be compiled. It installs
the module `crn_craft` and two commands, `crn-craft` and `crn-craft-import`.
pip also installs NumPy, pandas and xarray, because results are pandas data
frames and xarray arrays. Wheels exist for CPython 3.10 and newer on Linux
x86_64, Windows x64 and macOS (x86_64 and arm64). Other platforms build from
the source distribution with a Rust toolchain.

## Check the installation

One call shows the version, and a second one shows that the engine answers:

```python
import crn_craft as crn

print(crn.__version__)
print(crn.describe()["schema"])
```

```text
0.7.1
crn-craft.describe/1.0
```

## Coming from another tool

Two things are worth trying first. Many trajectories run as one ensemble in
one call, on native threads ([chapter 3](03-running.md)). A model can be
compiled from Python, with the same results in less time
([chapter 5](05-speed.md)). [How CRN-craft compares](../comparison.md) puts
CRN-craft beside GillesPy2 and rebop, with features and measured speed.

## The command-line tool is optional

The wheel already has the command `crn-craft`, so this guide needs no other
tool. A native tool with the same commands and reports installs with Rust
1.85 or newer: `cargo install crn-craft-cli --locked`. If both are installed,
the one that comes first on your `PATH` runs, and `crn-craft --version` tells
them apart: the Python command prints `"interface": "python"`.
[Chapter 10](10-command-line.md) uses either.

## Build from source

To use the current source instead of a release, clone the repository, build
the command-line tool with cargo and the wheel with maturin:

```sh
git clone https://github.com/dietmarwo/CRN-craft
cd CRN-craft
cargo build --release --workspace --locked
./target/release/crn-craft --help
python -m pip install maturin==1.9.6
maturin build --release --locked --out dist
python -m pip install --force-reinstall dist/*.whl
```

## Optional extras

- **Plots.** The plotting example of [chapter 4](04-results.md) needs
  matplotlib: `python -m pip install matplotlib`.
- **Fitting.** `python -m pip install "crn-craft[fit]"` adds the pinned
  SciPy, fcmaes-rust and SQP-craft versions of the fitting examples
  ([chapter 13](13-fitting.md)).
- **Compiled models.** `crn.compile_loop` needs a Rust toolchain, 1.85 or
  newer, from [rustup](https://rustup.rs) ([chapter 5](05-speed.md)).
- **R.** The Python package works from R through reticulate
  ([chapter 15](15-from-r.md)).

## How to read this guide

Read chapters 1 to 4 in order. After them you can build a model, run an
ensemble and inspect, plot and save the result. The later chapters can be
read in any order, and each assumes only chapters 1 to 4.

- Every chapter starts with a setup block that imports the package as
  `import crn_craft as crn` and builds the chapter's model. The blocks of a
  section continue from the setup block and from earlier blocks of the same
  section, so the setup block and one section are a complete program.
- Work in an empty folder. A chapter writes every file it reads, or names
  the repository file to copy.
- A checker runs every block and compares each printed output with the block
  before it. For a single question, start from the [task index](../howto.md).

See also: [A first model](02-first-model.md),
[Troubleshooting](troubleshooting.md),
[Migrating from earlier releases](migrating.md).
