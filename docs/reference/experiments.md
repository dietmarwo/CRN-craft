# Reproducible stochastic experiments

A protocol gives the time and conditions of an experiment. A count model
describes the reactions. Timed interventions and reaction completions use
the native count scheduler. Prepared projections and likelihoods are
available in Python with the optional `crn-craft[projection]` dependency.

## Scheduled stochastic simulation

```python
import crn_craft as crn

model = crn.Model("timed_birth")
model.add_species("X", 0)
model.add_parameter("birth", 2.0)
model.add_reaction(crn.mass_action("birth"), [], ["X"], name="arrival")

protocol = crn.Protocol(time_unit="hour")
protocol.at(3.0, parameters={"birth": 4.0}, partition={"X": 0.5})
run = crn.ensemble(
    model, end=5, report_times=[0, 3, 5], protocol=protocol,
    seed=123, replicates=20_000, threads=4,
)
result = run.simulate()
```

At hour three the daughter count is Poisson with mean 3. At hour five
it is Poisson with mean 11. Reports are right-continuous. Completions due
at a boundary precede external action groups. Equal-time action groups
follow their order in the protocol. Each group reads its pre-group counts.
A stopped trajectory retains its counts.

The methods are linear `static_direct` and `next_reaction` for closed
attribute-free count models. At each boundary the engine updates the
counts, recompiles numeric rates during preparation when parameters
change, and restarts the clock from the conditional law. Each replicate
keeps one random stream. Count actions include signed increments and
absolute assignments. A partition independently retains each molecule
with a specified probability. The `time_unit` field documents the unit
shared by the inputs.

`crn-craft.run/1.1` is a versioned request with required `protocol`. The
native CLI supports its usual `check`, `run`, `--archive` and `verify`
commands. `crn-craft.result/1.1` logs calendar events with `kind`,
`pending` and `cancelled`. An unfiltered event log stores changed counts
before and after a deterministic event in its participant fields. Report
`events` and `rule_events` count stochastic initiations. Interventions
and completions consume work from the same replicate event budget. To
migrate a 1.0 count request, change its version and add `protocol`, then
choose linear Static Direct or Next Reaction.

## Fixed-delay completion

`Protocol.delay(reaction, duration, products=..., owner=...)` adds a
fixed-time completion to a reaction. The named reaction still applies
its immediate stoichiometry at initiation. An optional owner is a busy
count species created once by that initiation rule. Degradation,
partition and deletion of a busy copy cancel its individual completion.
A successful completion consumes the owner and creates the specified
products. The pending queue is bounded by `max_pending`. Count
projections need Markov state and therefore refuse this delay history.

The delayed-reaction implementation follows the event-queue semantics
of the [modified next-reaction method](https://arxiv.org/abs/0708.0370).
The [delayed feedback example](../../examples/delayed_feedback/study.py)
compares a power spectrum and an autocorrelation against an independent
scheduler. For a paper specifying hourly degradation trials, model that
discrete sequence separately from the continuous-time degradation hazard.
[Fixed replication-time example](../guide/17a-mtdna-mechanisms-models.md)

## Prepared probability laws

```python
projection = crn.prepare_fsp(model, bounds={"X": 48}, protocol=protocol)
law = projection.solve(5, [0, 3, 5], tolerance=1e-12)
daughter = law.marginal("X", at=1)
continued = projection.solve(5, [5], initial=law.snapshots[1],
                             tolerance=1e-12)
```

Prepare with `initial={(0,): 0.25, (1,): 0.75}` to seed a distribution.
`solve(parameters=...)` reuses the reachable graph while evaluating
numerical rates. Topology, initial-support or bound changes require a
fresh graph when they would introduce states or edges. An action maps
the retained law through its count kernel. Escaped probability and
previously stopped probability carry forward.

Pass `signals={"arrival": crn.Signal.table(times, rates)}` to replace
the count-independent rate of a named rule by a prescribed time series.
Tables use linear interpolation or an explicit right-continuous step
curve. `Signal(function, change_times=...)` also accepts a Python rate
function. Python callables are supplied by the caller when rerunning the
experiment. The two integrators for smoothly forced segments are sparse
BDF and adaptive exponential midpoint. Constant segments use positive
uniformization. Each snapshot exposes the absorbing projection sink,
Poisson series-tail estimate, mass defect and minimum probability.
Check tolerance refinement on the selected observable for a smooth
signal. See the [RNA fit replay](../../examples/yeast_rna_fit/native_fit.py).

Use `stationary_fsp(model, bounds)` for one finite closed communicating
class or specify `class_state`. With a `return_state`, escaped rates
redirect to a retained state. Successive enlargements compare marginal
CDFs and outward flux. These diagnostics are empirical on unbounded
networks. Certified error bounds require additional recurrence evidence.
[Stationary FSP method](https://arxiv.org/abs/1704.07259)

## Spatial models

`SpatialModel` expands local reaction templates and directed first-order
hops to an ordinary count model. For equal-volume lattice neighbors,
`diffusion(..., diffusivity=D, spacing=h)` adds both directions with rate
`D / h**2` per copy. Specify transport rates directly for unequal
volumes and other graph edges. An explicit reservoir maintains its
initial concentration of copies at the boundary.

`spatial_turnover` is a grouped exact sampler for the two-genotype chain
in the [spatial spread study](../guide/19-mtdna-spread.md). It updates
the affected regions after each event and can carry neutral clone labels
when a mutation creates a new identity. Opt in through
`study.py --source-chain --wave-sampler local`. Fiber IDs, checkpoints and
bootstraps all operate on independent fibers.

## Observation models

`crn_craft.likelihood` maps states to counts, sums, joint measurements
and intervals. It implements binomial detection and binned Gaussian
fluorescence. `histogram_score` returns grouped log scores with the
effect of missing probability bounded separately. `mean_score` uses a
sample covariance. The RNA example binds candidates to the prepared
projection and supplies the score to fcmaes. The historical `legacy_tail`
option reproduces the guide's coarse overflow convention.

## Published experiment files

Use `crn.read_experiment(path, convention="stochastic_count_propensity")`
and then `preflight()` or `run()`. The wheel also installs
`crn-craft-experiment EXPERIMENT --convention stochastic_count_propensity`,
with `--check` and `--archive DIRECTORY`. The archive contains the sealed
native task runs and a report of the source members and interpreted
conditions.

| Format | Supported subset |
|---|---|
| SBML L3V1 and L3V2 Core | Fixed compartments, integer item amounts, concentration references with explicit compartment size, nonnegative binomial mass-action hazards, acyclic scalar parameter assignments, and instantaneous time-triggered assignments with resolved priorities |
| SED-ML L1V3 and L1V4 | Local SBML documents and model references, parameter `changeAttribute`, Gillespie direct uniform time courses, reset-only repeated tasks with a single uniform range, time/species report datasets and plot series |
| COMBINE/OMEX | One resolved master SED-ML document, local model members, bounded archive reads and member hashes |

The explicit `split_signed_mass_action` convention decomposes a signed
law into a forward jump and a reverse jump only when exact polynomial
factorization gives two nonnegative hazards. It handles the two signed
promoter-binding laws in the [published oscillator archive](https://github.com/biosimulators/Biosimulators_test_suite/blob/ae6afe13c37aa4b896eb5f839ab0e26ba3bc2f5b/examples/sbml-core/Vilar-PNAS-2002-minimal-circardian-clock-discrete-SSA.omex).
Use the strict count convention for ordinary nonnegative irreversible laws.
Ambiguous rate conventions, unsupported constructs and unresolved units
produce a diagnostic during translation. [SBML](https://sbml.org/documents/specifications/level-3/version-2/core/) · [SED-ML](https://sed-ml.org/specifications.html)
