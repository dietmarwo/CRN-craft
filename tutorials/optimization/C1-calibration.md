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
rates within 0.35 to 2.5 times the generating rates, starting from 1.25 times
them. Start and box depend on the known truth, which a real calibration would
not have; `best_is_start` flags a fit that never improved on the start. This is
a controlled identifiability exercise with known truth, not an inference model
for arbitrary experimental measurements.

## Run and inspect

Follow [setup](README.md), then run from CRN-craft:

```bash
python tutorials/optimization/run.py C1 --output /tmp/c1-smoke
python tutorials/optimization/run.py C1 --include-next-reaction --output /tmp/c1-both
python tutorials/optimization/run.py C1 --preset study --threads 8 --output /tmp/c1-study
```

The smoke run requests 32 objective calls per route, 16 training paths and 1,024
held-out paths. The study uses 128 calls, 32 training paths, three optimizer
seeds and 4,096 held-out paths.

The models are `crn-craft.model/1.0` documents (the packaged
immigration–death example and a two-rule binding model written with the
builder) and the requests are `crn-craft.run/1.0` ensembles of the Static Direct
method (sampler `chacha8-u53-ziggurat/rules-static-direct/1`), written with
`crn.ensemble`. `--include-next-reaction` adds a separate arm with the
`next_reaction` method (sampler `chacha8-u53/rules-next-reaction/1`) and its own
check; both are exact, but the two methods draw different paths and are never
pooled. Each row records its policy (`direct` or `next_reaction`), sampler and
request schema.

For each model, policy and optimizer seed, the fresh arm reconstructs the model
and request for every callback. The reusable arm prepares once with
`crn.ensemble(...)` and calls `template.bind_parameters(rates).simulate()`. Both use the same optimizer initialization, rate
bounds, numerical policy and root training seed. Inspect the `*.trace.json`
files: candidate coordinates and objective values must agree exactly between
these two routes. The parameter names are sorted, and the objective maps the
coordinates to them by that order.

At the fitted point, a model prepared afresh from the fitted values must give
the same result as the binding: the same statistics, events and sampler. Only
the hashes of the serialized model and request may differ, because a binding
serializes the bound model itself. This replay statement concerns the two CRN
routes under their documented method; it says nothing about another simulator's
random-number stream.

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
finite-path objective, and an optimizer can fit its simulation error. Every
row therefore also records `truth_training_loss`, the training loss of the
true rates on the same sample. A fitted loss far below it means the search
fitted that sample's noise. One training sample (seed 77201) serves every
optimizer seed, so the three study seeds replicate the search, not the data:
a bias they share, such as a low binding reverse rate, comes from that sample.

## Prompt and follow-up

> Fit both networks with fresh and reusable preparation. Keep candidate streams,
> numerical policies and training seeds matched. Account for preparation,
> optimization and final validation separately; explain rate error and held-out
> moment error before interpreting any timing difference.

Repeat the study over additional training, optimizer and validation seeds.
Increase the training paths while holding the requested search calls fixed. Does prediction
improve before rate recovery? Then run the [X1 CPU layouts](https://github.com/dietmarwo/experiment-craft/blob/main/examples/campaigns/simulators/README.md):
compare complete fits at a fixed total CPU allowance. Do not infer full-fit
speedup from one SSA timing or from one smoke measurement.

Background: [reusable parameters](../../docs/reusable-parameters.md) and
[the existing parameter comparison](../../docs/comparison-parameter-optimization.md).
