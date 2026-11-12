# C3 — Oscillator Pareto fronts and behavior archives

Use one trace evaluator for multi-objective search, quality diversity and final
validation. The executable [oscillators.py](oscillators.py) builds an additive
Hill network as a `crn-craft.model/1.0` document, runs native MODE and MAP-Elites through their Python APIs, and
checks a second simulator backend before interpreting the archive.

## Define the model before comparing it

The search model has three species. Each species has basal production and
linear degradation; each active signed edge adds a separate Hill production
channel (a propensity expression such as `vmax / (1 + (X2 / k)^n)`). The default topology is a three-edge repression ring. Splitting a
summed production rate into channels changes random-number consumption but
preserves the continuous-time Markov generator when propensities and state
changes agree.

The generator audit includes zero populations, half-saturation, random states,
multiple incoming edges and both signs. It compares CRN's aggregated production
and death rates, read from `Model.analyze()` at each state, with the explicit
additive formula. A second native reference,
[crn_reference.rs](https://github.com/dietmarwo/fcmaes-rust/blob/main/tutorials/oscillator-topology-search/examples/crn_reference.rs),
uses the existing topology tutorial's **local ReBop 0.9.7 compatibility backend**
and six aggregate channels. This is not a claim to benchmark unmodified
upstream ReBop. The small reference accepts only bounded parameters, counts,
paths and time grids and is invoked with a 30-second timeout.

The canonical nine-species, sixteen-reaction Vilar model is kept separate. A
short control checks its two promoter conservation laws. That control is not a
Vilar oscillation optimization, and the three-species teaching model is not a
replacement for it.

## Run and inspect

Follow [setup](README.md), including the fcmaes-rust reference example, then:

```bash
python tutorials/optimization/run.py C3 --output /tmp/c3-smoke
python tutorials/optimization/run.py C3 --preset study --threads 8 --output /tmp/c3-study
```

Override the reference location with `--reference-binary /absolute/crn_reference`.
MODE and MAP-Elites each receive 32 smoke or 256 study requests. Each candidate
uses 4 or 16 training paths, sampled at 201 equally spaced times over 125 time
units, with a 200,000-event cap per path. The three fitted coordinates are log
Hill strength, half-saturation and degradation rate. Topology is held fixed.
All SSA runs are Static Direct ensembles (`crn.ensemble`, request
`crn-craft.run/1.0`) whose candidates are bound with `bind_parameters`; the
summary records the request schema and every evaluation's sampler.

## State the scoring rule

After a 20% burn-in, the evaluator uses CRN's oscillation diagnostic on X0.
Its FFT resolves only periods `W/k`, where `W` (about 100.6 time units) is the
retained window, and a path counts as oscillatory only with `k >= 3` cycles,
so oscillatory periods are at most about 33.5. The horizon is chosen so that
the 20-unit target is one of them (`k = 5`, about 20.1); its neighbors are
about 25.2 and 16.8. Quality combines relative period error against 20 time
units and one minus the spectral peak fraction. Cost is the mean total molecule count divided by 100.
Fragility combines the nonoscillatory-path fraction and period coefficient of
variation. Amplitude is the X0 90th-minus-10th percentile. Missing periods have
a declared fallback of 160; the original failure classification is retained.
These are teaching metrics, shared exactly across both backends. They differ
from the earlier fcmaes-rust oscillator tutorials' scoring and parameterization.

MODE minimizes the three objectives. MAP-Elites minimizes their sum and indexes
period/amplitude on a fixed 4-by-4 grid over `[0,80] × [0,80]`. This transparent
grid allows independent niche-index checks and held-out niche-drift reporting;
it is not a CVT archive. Out-of-range descriptors clip to an edge cell and are
counted in the summary; a fallback-period cell is not evidence of oscillation.
The period descriptor averages over all paths, nonoscillatory ones included, so
period columns beyond about 33.5 hold averages of slow or ambiguous paths, not
slow oscillators. Every evaluation, archive elite and validation row therefore
records its `oscillatory_fraction`. Noncompleted evaluations are retained as
failures and cannot become valid archive entries.

## Validate the selections

Re-evaluate three Pareto extremes and every occupied archive cell on the
validation seed (9101), with 16 smoke or 128 study paths. `pareto.json`, `archive.json` and
`evaluations.json` retain candidates and trace metrics. The summary reports
changes in descriptors and niche assignments, including nonoscillatory paths.
Numerical completion alone is not evidence of a reliable oscillator.

At the default parameters and one selected design, run both backends through
the same score. CRN uses its own control seed (9201) and the reference another
(9301), so neither side reuses the training paths that selected the design.
Compare path-averaged molecule counts with a fixed six-standard-error gate
using independent paths as the sampling units. This avoids treating
correlated time samples as independent observations. Both backend metrics are
reported; paths are never required to match. The reference backend does not
export an event count, so its work field is explicitly unavailable.

## Prompt and follow-up

> Verify additive generator equivalence, then build a robust-oscillator Pareto
> set and a behavior archive. Use one scoring implementation, seeds that never
> reuse selection paths, and explicit failed paths. Explain niche drift and preserve the distinct
> Vilar reference rather than merging different oscillator models.

Increase validation paths, inspect trajectories behind extreme scores, and only
then compare topology proposers using the same inner protocol. Adding an
agent-guided topology search or stock MO/QD job application is a separate
extension. See the original [topology tutorial](https://github.com/dietmarwo/fcmaes-rust/blob/main/tutorials/oscillator-topology-search/README.md)
and [Vilar tutorial](https://github.com/dietmarwo/fcmaes-rust/blob/main/tutorials/rebop-oscillator/README.md).
