# C1 — Calibrate a stochastic network without rebuilding topology

A reaction network often keeps its species and reactions while an optimizer
changes only rate constants. This lesson measures what reusable parameter
preparation saves over a **complete optimization**, and checks that it preserves
the result. The implementation is [calibration.py](calibration.py).

## Question and reference models

Fit immigration–death and reversible conversion `A ⇌ B`. Both have independent
finite-time moments, making simulator error distinguishable from fitting error.
For immigration rate `a`, death rate `b` and zero initial count, the distribution
is Poisson with mean and variance `a/b * (1-exp(-b*t))`. With 30 initial A
molecules, the B population is binomial with success probability
`a/(a+b) * (1-exp(-(a+b)*t))`; A+B stays 30.

The synthetic observations are these exact means and variances at seven times
from 0 to 6. The objective sums mean-square discrepancies of the simulated mean
and sample variance, each scaled by one plus its exact target. It fits log
rates within 0.35 to 2.5 times the generating rates. This is a controlled
identifiability exercise with known truth, not an inference model for arbitrary
experimental measurements.

## Run and inspect

Follow [setup](README.md), then run from CRN-craft:

```bash
python tutorials/optimization/run.py C1 --output /tmp/c1-checked
python tutorials/optimization/run.py C1 --include-fast --output /tmp/c1-fast
python tutorials/optimization/run.py C1 --preset study --threads 8 --output /tmp/c1-study
```

The smoke run requests 32 objective calls per route, 16 training paths and 1,024
held-out paths. The study uses 128 calls, 32 training paths, three optimizer
seeds and 4,096 held-out paths. `--include-fast` adds a separate explicit fast
arm. It does not silently replace the checked policy.

For each model, policy and optimizer seed, the fresh arm reconstructs the model
and request for every callback. The reusable arm prepares once and calls
`template.simulate(rates)`. Both use the same optimizer initialization, rate
bounds, numerical policy and root training seed. Inspect the `*.trace.json`
files: candidate coordinates and objective values must agree exactly between
these two routes. Parameter coordinate order is checked against
`template.parameter_names`.

At the fitted point, ordinary execution of the binding's serialized model and
request must equal `bound.simulate()`. This replay statement concerns the two
CRN routes under their documented policy; it says nothing about another
simulator's random-number stream.

## Read the evidence

Each arm in `summary.json` records requested and actual calls, preparation,
search and validation time, training and held-out losses, rate error, and
simulation paths/events. The setup timer is the common template preparation
cost and must be charged once to the reusable route. The fresh route's repeated
preparation is inside its callback time. Validation is deliberately separate
from search and includes ordinary/bound replay.

A disjoint root seed supplies the final ensembles. Fitted-point sample means
must lie within six estimated standard errors of the independently known
moments, with a declared small denominator floor. This is a conservative
engineering check, not a confidence statement about fitted parameters.
The true-rate held-out loss supplies a finite-sample reference. A low training
loss alone is insufficient: fixed training seeds produce a deterministic but
finite-path objective, and an optimizer can fit its simulation error.

## Prompt and follow-up

> Fit both networks with fresh and reusable preparation. Keep candidate streams,
> numerical policies and training seeds matched. Account for preparation,
> optimization and final validation separately; explain rate error and held-out
> moment error before interpreting any timing difference.

Repeat the study over additional optimizer and validation seeds. Increase the
training paths while holding the requested search calls fixed. Does prediction
improve before rate recovery? Then run the [X1 CPU layouts](https://github.com/dietmarwo/experiment-craft/blob/main/examples/campaigns/simulators/README.md):
compare complete fits at a fixed total CPU allowance. Do not infer full-fit
speedup from one SSA timing or from one smoke measurement.

Background: [reusable parameters](../../docs/reusable-parameters.md) and
[the existing parameter comparison](../../docs/comparison-parameter-optimization.md).
