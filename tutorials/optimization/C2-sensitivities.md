# C2 — Fit kinetics with BDF sensitivities and test stochastic transfer

This lesson extends the closed-form consecutive-reaction example in
[SQP6](https://github.com/dietmarwo/SQP-craft/blob/main/tutorials/reaction-kinetics/README.md)
to nonlinear enzyme kinetics. It uses the native optional diffsol integrator,
not a replacement ODE implementation inside the optimizer. See
[sensitivity.py](sensitivity.py).

## Separate four questions

The model is `S + E ⇌ ES → E + P`, initially 50 S and 10 E molecules. Fit `kon`
and `kcat`, keeping `koff=4` fixed. Synthetic observations come from an
independently written DOP853 integration at nine times from 0 to 4. Its tight
tolerances and explicit reaction equations provide an integration reference;
they are separate from the native BDF sensitivities.

Ask whether the integrator is accurate, whether the derivatives are consistent,
whether the observations identify the rates, and whether a deterministic fit
transfers to stochastic paths. Passing one question does not answer the others.

## Run the experiment

Build the optional workspace with Rust 1.89 as described in [setup](README.md).
From CRN-craft:

```bash
python tutorials/optimization/run.py C2 --output /tmp/c2-smoke
python tutorials/optimization/run.py C2 --preset study --threads 8 --output /tmp/c2-study
```

Use `--ode-binary /absolute/path/crn-diffsol` for a nondefault build location.
Each native subprocess has a 30-second outer timeout because this optional
integrator has no internal step/wall budget. The smoke preset permits 48 native
solves per strategy; the study permits 160. Audits and final SSA validation are
additional work and are counted separately.

First, compare native trajectories with DOP853 at a point away from truth.
Require normalized state error below `2e-6`; tighten BDF tolerances and require
error below `2e-7`. Then compare the analytic objective gradient with central
differences in **log-rate coordinates**, requiring absolute error below `2e-5`.
The chain rule multiplies rate sensitivities by the corresponding rate. The
code explicitly reorders native sensitivity parameters before using them.

Next, inspect the scaled observation Jacobian. All observed species and times
have rank two in this fixture; only terminal product concentration has rank
one. This local rank demonstration shows an ambiguity. It does not establish
global uniqueness or practical parameter confidence.

## Compare solvers with honest work accounting

Three strategies use the same bounds and observations: BiteOpt alone, two
SQP starts, and BiteOpt followed by SQP. The hybrid reserves one third of the
native-solve cap for global exploration. SQP uses its active-set QP backend and
an analytic gradient. Repeated function/derivative requests at the same point
share a cached native result. A budget-exhausted evaluation is explicitly
undefined; the best finite point actually evaluated is retained, even if a
local solve stops before stationarity.

`summary.json` reports each native solve count, solver termination, native and
independent losses, relative rate error, and held-out SSA events/paths. Every
successful native call retains its model and result in the strategy directory.
Inspect a stopped SQP status before making a stationarity claim. A finite
incumbent does not imply convergence.

At the selected rates, 256 smoke or 4,096 study SSA paths use a disjoint seed.
Substrate and enzyme conservation are checked. The report measures normalized
SSA-mean versus ODE discrepancy rather than demanding equality: nonlinear
count kinetics do not generally close under expectation.

## Prompt and follow-up

> Fit the nonlinear enzyme model with independently audited sensitivities.
> Compare global search, multistart SQP and global-to-local refinement under a
> shared native-solve cap. Demonstrate a rank ambiguity and explain stochastic
> transfer separately from deterministic fit quality and solver termination.

Remove intermediate observations, add declared measurement noise, or move to a
new initial population. Re-audit gradients and integration before comparing
fits. Do not tune integration tolerances after inspecting solver rankings.
The [existing ODE fitting example](../../examples/circuit/fit_ode.py) and
[optional workspace](../../integrations/diffsol/Cargo.toml) provide further context.
