# Optimization tutorials

Ten executable lessons put an optimizer around a CRN-craft simulation. Each
separates what the optimizer returned from the numerical checks and from
the scientific reading of the result. They are the optimization group of
the [tutorials](../../docs/tutorials.md); the simulation group, which
rebuilds published stochastic models without an optimizer, is in the
[guide](../../docs/tutorials.md#simulation-tutorials).

## Method lessons C1–C4

These use synthetic data and teach the route. Start with
[C1](C1-calibration.md), then choose the lesson that matches your model.
They connect CRN-craft to fcmaes-rust and SQP-craft.

| Lesson | Scientific question | Native route |
|---|---|---|
| [C1 — calibration](C1-calibration.md) | Does reusable preparation preserve a complete fit? | `crn.ensemble` + `bind_parameters` + BiteOpt |
| [C2 — sensitivities](C2-sensitivities.md) | Are derivatives correct, rates identifiable, and ODE fits transferable? | Optional BDF integration + SQP/BiteOpt + SSA |
| [C3 — oscillators](C3-oscillators.md) | Which trade-offs and behaviors survive independent paths? | SSA + MODE + MAP-Elites + ReBop compatibility control |
| [C4 — rules](C4-rules.md) | Can finite-path fits of exact aggregation laws distinguish kernels? | Bounded native `crn-craft run` + BiteOpt |

## Publication studies C5–C10

Each of these reconstructs a published optimization study, checks the
objective against the paper or its archive before any search, runs
[fcmaes](https://github.com/dietmarwo/fast-cma-es) and then verifies the
finalists. A lesson that misses the published value says by how much and
why.

| Lesson | Study | Route | Outcome |
|---|---|---|---|
| [C5 — RNA count calibration](C5-rna-counts.md) | Poovathingal & Gunawan 2010: three rates from count distributions | SSA, prepared FSP and DE | With 10,000 trajectories the three rates are recovered within 3%. With 100 trajectories two of four data sets miss the second rate by a factor 1.7: the score has a flat valley |
| [C6 — stimulus model selection](C6-stimulus-fit.md) | Gómez-Schiavon et al. 2017: seven hypotheses on 376 measured cells | Stationary and transient FSP, joint likelihood, AIC and BIC | Both criteria select the published model; every fitted rate is within one posterior standard deviation of the paper's table |
| [C7 — virtual reactor design](C7-reactor.md) | Felton et al. 2021: yield against waste in 50 evaluations | Concentration BDF and MODE | At 50 evaluations MODE cannot be told from random search and is far below the paper's best strategy; with 500 evaluations it passes that strategy's median |
| [C8 — sampling design](C8-experiment-design.md) | Fox & Munsky 2019: sampling times and a stimulus pulse by Fisher information | FSP forward sensitivities | The sampling optimum and the ranking of the two published designs in 97 of 100 parameter samples are reproduced, the archived matrices entry by entry; no continuous design search is claimed |
| [C9 — STAT5 calibration](C9-stat5.md) | Boehm et al. 2014: nine parameters of a signaling model, as a PEtab problem | PEtab v1 tables, concentration BDF, DE and CMA-ES | The archived likelihood is reproduced to 1e-11 and the native model to 2e-6 nats; plain retries stop 9 nats above it, a search with profiled noise scales reaches it |
| [C10 — pathway calibration](C10-pathway.md) | Moles, Mendes & Banga 2003: 36 parameters over 16 experiments | Batched concentration BDF, DE and CMA-ES | The objective is verified, but searches of up to 20,000 evaluations reach costs of 62 to 667 where the paper reports 0.0013 after 2.8 million evaluations |

Use the shared [installation and command-job guide](https://github.com/dietmarwo/optimization-jobs/blob/main/tutorials/simulators/README.md).
It pins the Python dependencies, builds the optional binaries and lists all
ten lessons. The CRN Python package must be built from this checkout; the
`crn-craft` command and the optional diffsol binary are separate builds. C3's reference
binary belongs to fcmaes-rust's existing oscillator tutorial workspace.

From this repository, with the environment activated:

```bash
python tutorials/optimization/run.py C1 --output /tmp/crn-c1-smoke
python tutorials/optimization/run.py C1 --preset study --include-next-reaction --output /tmp/crn-c1-study
python -m unittest discover -s tutorials/optimization -p 'test_*.py'
python -O -m unittest discover -s tutorials/optimization -p 'test_*.py'
```

Choose a new output directory for every run. The driver preserves a
`summary.json` with `passed: false` and an error when a numerical check fails;
`state` is `running` until the chapter ends, so a summary read mid-run is not
mistaken for a result. Chapter-specific options a chapter does not use are
rejected.
`passed` means the chapter's declared checks passed; it does not mean the
optimizer identified the true parameters or discovered a robust oscillator.
Checks use explicit exceptions and remain active under `python -O`.

SSA requests are `crn-craft.run/1.0` ensembles of the Static Direct method
(ChaCha8 streams). C1's `--include-next-reaction` adds an arm with the
next-reaction method; the other chapters, and the BDF integrator (C2), have no
such switch.

The smoke preset exercises the complete workflow with small search budgets.
The study preset increases calls, paths or optimizer seeds; neither is a
statistical power guarantee. Numerical policies and validation thresholds are
fixed in the scripts before execution. Build/import costs are outside a
chapter's search timer; use [X1](https://github.com/dietmarwo/experiment-craft/blob/main/examples/campaigns/simulators/README.md)
for complete process and campaign accounting. Per-arm `*.run.json` files use
the existing tutorial version-1 convention alongside simulator-owned evidence.

These are library/command examples. They do not register new stock CRN job
models or replace the synthetic `app.crn.circuit-score`. Model files, requests, actual work and held-out results are the
primary evidence. See [validation notes](VALIDATION.md) for the recorded checks.


## Setup for C5–C10

Build the CRN-craft wheel from this checkout
(`maturin build --release --locked --out dist`, then
`pip install dist/*.whl`) and install the other packages:

```bash
pip install "crn-craft[fit]" fcmaes==2.2.0 matplotlib PyYAML python-libsbml
```

fcmaes brings loguru and scikit-learn. PyYAML and python-libsbml are needed
by C9 only; C1 to C4 need the `fit` extra. C7, C9 and C10 use the optional
concentration solver, which `--ode-binary` finds at its default path after

```bash
cargo build --release --locked --manifest-path integrations/diffsol/Cargo.toml
```

C5 and C10 generate their own data. C6 to C9 read pinned files that are
checked against the manifest
[`examples/publication_tutorials/sources.json`](../../examples/publication_tutorials/sources.json);
fetch them once, then pass `--data-dir DIR` or set `CRN_PUBLICATION_DATA`:

```bash
python tutorials/optimization/publications/sources.py --data-dir DIR C6 C7 C8 C9
```

Single files are pinned by SHA-256. The C8 archive (62 MB) is pinned by a
digest over its member files, so a GitHub tarball and a `git archive` of
the pinned revision both pass. Later runs read the files locally.

```bash
python tutorials/optimization/run.py C6 --preset smoke --data-dir DIR --output /tmp/crn-c6-smoke
python tutorials/optimization/run.py C7 --preset study --data-dir DIR --output /tmp/crn-c7-study
python tutorials/optimization/figures.py
```

Each result directory keeps the model, the numerical checks and every
proposal in evaluation order. A search can continue from a finished run
with `--resume-from`; it refuses a checkpoint of another objective, seed or
optimizer. Use `--threads` for one layer of concurrent candidate
evaluation and set `OPENBLAS_NUM_THREADS=1` with it, so the numerical
libraries do not start their own worker pools. Smoke presets test the full
route with small budgets; study options are on each lesson's page.

The tables behind the figures of each lesson are saved under
[`results/C#/`](results/), with a README that names the run, the seeds and
whether the numbers come from the delivered eight-thread study or from a
rerun. `figures.py` redraws every figure of C5 to C10 from those tables;
`--tables DIR` draws one lesson from a run of your own.
