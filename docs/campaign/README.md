# CRN-craft campaign guide

Baseline: GitHub main 5a0733444aa46c8cea2a5b61bd8601daf5890d79. Read
[acceptance.md](acceptance.md) for the frozen plan, and M3.md through M9.md for
implementation, executed evidence and limitations. Performance has ten separate
experiments, including retained/rejected candidates, in performance.md.
Ten more (P11–P20), with their declaration and results, are in
performance-p11-p20.md. For an introduction read the [tutorial](../tutorial.md),
the [architecture overview](../architecture.md) and the
[technical report](../technical-report.md).

| Use | Command / interface |
|---|---|
| Exact SSA and archives | `crn-craft simulate MODEL RUN NEW_DIRECTORY` |
| Alternative exact SSA | `crn-craft ssa MODEL REQUEST`; see [guide](../alternative-ssa.md) |
| Exact dependency-cached SSA | `crn-craft simulate-cached MODEL RUN NEW_DIRECTORY` |
| Verify exact archive | `crn-craft verify DIRECTORY` |
| Sparse FSP | `crn-craft fsp MODEL RUN` |
| Experimental Poisson leaps | `crn-craft tau MODEL TAU_REQUEST` |
| Python | `pip install wheel.whl`; `crn_craft.load_example`, `run_spec`, `simulate`, `fsp`, `tau` |
| Python native-backed CLI | installed `crn-craft`, or `python -m crn_craft` |
| Trajectory diagnostics | `crn_craft.observables` |
| Explicit import subsets | `python -m crn_craft.importers --help` |
| fcmaes worker | build `integrations/jobs/Cargo.toml`; existing retry/surrogate/evaluate profiles |
| BDF and forward sensitivities | build `integrations/diffsol/Cargo.toml` with Rust 1.89+ |
| Global to local fit | `examples/circuit/fit_ode.py` with `crn-diffsol` |
| Conservative topology screen | `examples/circuit/screen.py` |

The original run/model/result/archive schemas are unchanged. Discovery advances
to 0.5, listing additional commands; fsp-result, tau-request and tau-result are new
0.1 schemas generated from Rust. Approximate output cannot be verified by the
legacy exact-archive verifier or labelled an exact ensemble. FSP moments are
unnormalized projections; leakage/series/roundoff have separate meanings.

Default library and CLI have no new dependency over reviewed M2. Optional pins:
PyO3 0.29.0, diffsol 0.18.0, fcmaes-rust 583f33a and its actual SDK revision
602d31a. The tested Python fit uses fcmaes-rust 0.2.1 and SQP-craft 0.17.0.

Build the package with maturin 1.9.6. The delivered wheel is Linux x86_64 abi3
(Python 3.10+ ABI); execution is locally checked on Python 3.12. Other-platform
wheels are a configured CI job, not locally validated. No package is published.
Julia/reference coverage and broad biological detector validation remain explicit
limits, not implied passes. None of this replaces biological model validation.

Both exact implementations run the full analytic SSA and nonlinear CME acceptance
campaigns; their fixed-work reports are byte-identical. Optional worker/ODE JSON
schemas are generated from their Rust types by `describe`. The base ABI3 package
is tagged Python 3.10+; the pinned fitting extras require newer Python according
to their package metadata (the executed study uses 3.12).
