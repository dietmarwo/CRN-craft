# Optimization with stochastic reaction networks

Four executable lessons connect CRN-craft to fcmaes-rust and SQP-craft. Each
separates an optimizer's result from numerical acceptance and scientific
interpretation. Start with [C1](C1-calibration.md), then choose the lesson that
matches your model.

| Chapter | Scientific question | Native route |
|---|---|---|
| [C1 — calibration](C1-calibration.md) | Does reusable preparation preserve a complete fit? | SSA parameter templates + BiteOpt |
| [C2 — sensitivities](C2-sensitivities.md) | Are derivatives correct, rates identifiable, and ODE fits transferable? | Optional BDF integration + SQP/BiteOpt + SSA |
| [C3 — oscillators](C3-oscillators.md) | Which trade-offs and behaviors survive independent paths? | SSA + MODE + MAP-Elites + ReBop compatibility control |
| [C4 — rules](C4-rules.md) | Can finite-path fits of exact aggregation laws distinguish kernels? | Bounded native rules CLI + BiteOpt |

Use the shared [installation and command-job guide](https://github.com/dietmarwo/optimization-jobs/blob/main/tutorials/simulators/README.md).
It pins the Python dependencies, builds the optional binaries and lists all
eight lessons. The CRN Python package must be built from this checkout; the
rules CLI and the optional diffsol binary are separate builds. C3's reference
binary belongs to fcmaes-rust's existing oscillator tutorial workspace.

From this repository, with the environment activated:

```bash
python tutorials/optimization/run.py C1 --output /tmp/crn-c1-smoke
python tutorials/optimization/run.py C1 --preset study --include-checked --output /tmp/crn-c1-study
python -m unittest discover -s tutorials/optimization -p 'test_*.py'
python -O -m unittest discover -s tutorials/optimization -p 'test_*.py'
```

Choose a new output directory for every run. The driver preserves a
`summary.json` with `passed: false` and an error when a numerical check fails;
`state` is `running` until the chapter ends, so a summary read mid-run is not
mistaken for a result. Options a chapter does not use are rejected.
`passed` means the chapter's declared checks passed; it does not mean the
optimizer identified the true parameters or discovered a robust oscillator.
Checks use explicit exceptions and remain active under `python -O`.

SSA requests use CRN-craft's default `fast` policy (`crn-craft.ssa-run/0.3`);
C1's `--include-checked` adds a checked-policy arm. The rules engine (C4) and
the BDF integrator (C2) have no such switch.

The smoke preset exercises the complete workflow with small search budgets.
The study preset increases calls, paths or optimizer seeds; neither is a
statistical power guarantee. Numerical policies and validation thresholds are
fixed in the scripts before execution. Build/import costs are outside a
chapter's search timer; use [X1](https://github.com/dietmarwo/experiment-craft/blob/main/examples/campaigns/simulators/README.md)
for complete process and campaign accounting. Per-arm `*.run.json` files use
the existing tutorial version-1 convention alongside simulator-owned evidence.

These are library/command examples. They do not register new stock CRN job
models, replace the synthetic `app.crn.circuit-score`, or add a rules parameter
template API. Model files, requests, actual work and held-out results are the
primary evidence. See [validation notes](VALIDATION.md) for the recorded checks.
