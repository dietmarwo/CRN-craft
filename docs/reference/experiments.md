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
Groups at time zero apply before the initial report. The stop predicate
is tested after each group, and a stopped trajectory retains its counts.
With no reaction possible, the trajectory waits for the next group.

The methods are linear `static_direct` and `next_reaction` for closed
attribute-free count models. Other methods, and a Static Direct selector
other than `linear`, are refused with `CRN-R020` at `/method`. At each
boundary the engine updates the counts and switches to rate constants
that were compiled during preparation for every group that changes
parameters. Static Direct discards a wait that passes the boundary and
draws a new one under the new rates. Next Reaction draws fresh clocks for
all channels at every boundary. Both are exact because the exponential
waits are memoryless. A reaction drawn at exactly the time of a boundary
fires after it. Reports between boundaries draw no random numbers. Each
replicate keeps one random stream across all phases, so results do not
depend on `threads` or on the report grid.

Count actions include signed increments (`change`), absolute assignments
(`set`) and partitions. A partition independently retains each molecule
with the given probability. A group may hold at most one count operation
per species, and every operation reads the counts before the group.
Parameter values replace the current values until a later group changes
them. A parameter that sets an initial count cannot be changed. The
`time_unit` field is a label. The engine does not convert units: action
times, durations, the horizon and the model's rates must use the same
unit.

A request is refused with `CRN-R020` at `/protocol/...` for an unknown
species or parameter, conflicting operations on one species, a
probability outside [0, 1], or action times that decrease, are negative,
or lie after `end`. A protocol holds at most 4,096 action groups and
4,096 delays. A count change that would leave a species negative, or
above the count limit, fails the replicate with `CRN-R030` at
`/protocol/actions/i`. The message names the group, its time, the species
and its count.

`crn-craft.run/1.1` is a versioned request with required `protocol`. A
`crn-craft.run/1.0` request with `protocol` is refused with `CRN-R020`, and
1.0 requests without it run as before. The native CLI supports its usual
`check`, `run`, `--archive` and `verify` commands. `crn-craft.result/1.1`
logs calendar events with `kind`, `pending` and `cancelled`. An unfiltered
event log stores changed counts before and after a deterministic event in
its participant fields. Report `events` and `rule_events` count
stochastic reaction events. Interventions and completions consume work
from the same replicate event budget. To migrate a 1.0 count request,
change its version and add `protocol`, then choose linear Static Direct
or Next Reaction.

Protocol runs have their own sampler identities:
`chacha8-u53-ziggurat/rules-protocol-direct/1` (Static Direct with the
default Ziggurat waits), `chacha8-u53/rules-protocol-direct/1` (inverse
waits) and `chacha8-u53/rules-protocol-next-reaction/1`. The Rust tests in
`crates/crn-craft/tests/protocol.rs` check exact laws for all three: a
rate change and a channel that is switched on, the moments across a
binomial division, and the fixed-delay laws below.
`validation/check_paper_capabilities.py` checks the example above against
its two Poisson laws.

## Fixed-delay completion

`Protocol.delay(reaction, duration, products=..., owner=...)` adds a
fixed-time completion to a reaction. The named reaction still applies
its immediate stoichiometry at initiation. The duration is a positive
number in the run's time unit. An optional owner is a busy count species
that starts at zero, and only that initiation reaction may create it,
one copy per firing. Busy copies are exchangeable. When a reaction, a
partition or a count action removes k busy copies, the engine cancels k
pending completions of that owner, chosen uniformly at random from the
replicate's stream. An action that would increase a busy count fails the
replicate with `CRN-R030`. A completion that comes due consumes one busy
copy and creates the specified products. Pending completions are ordered
by due time, then by creation. At an exact tie, due completions run
first, then action groups in their declared order, then the report. The
queue holds at most `max_pending` entries (1 to 4,096, default 4,096). A
replicate that needs more fails with `CRN-R030`. Pending completions
survive parameter changes. `prepare_fsp` raises `ValueError` for a
protocol with delays, because a count state cannot represent the pending
history.

The delayed-reaction implementation follows the event-queue semantics
of the [modified next-reaction method](https://arxiv.org/abs/0708.0370).
`crates/crn-craft/tests/protocol.rs` checks three laws for each sampler.
Products completed by time t from Poisson starts at rate r with delay τ
are Poisson with mean r max(0, t − τ). With owner loss at rate δ that
mean is multiplied by exp(−δτ). When a partition keeps one of two busy
owners, each owner is kept with probability 1/2.

The [delayed feedback example](../../examples/delayed_feedback/study.py)
compares the first power-spectrum peak and the autocorrelation at twice
the delay with an independent delayed direct method written in Python.
With 12 paths per sampler, the peaks were 0.15390 (native) and 0.12877
(independent) rad per time unit, within the 0.03 tolerance. Over six
batches of 12 paths, the peak varied with a standard deviation of 0.0051
(native) and 0.0097 (independent). The tolerance is about ten frequency
bins, so the check detects only large changes of the period. The batch
means, 0.1487 and 0.1392, differ by 2.1 standard errors. This is
inconclusive, and more paths would be needed to resolve it. The example
is not run in CI.

For a paper specifying hourly degradation trials, model that discrete
sequence separately from the continuous-time degradation hazard.
[Fixed replication-time example](../guide/17a-mtdna-mechanisms-models.md)

## Prepared probability laws

```python
projection = crn.prepare_fsp(model, bounds={"X": 48}, protocol=protocol)
law = projection.solve(5, [0, 3, 5], tolerance=1e-12)
daughter = law.marginal("X", at=1)
continued = projection.solve(5, [5], initial=law.snapshots[1],
                             tolerance=1e-12)
```

Install the `projection` extra (`pip install "crn-craft[projection]"`)
for `prepare_fsp`, `stationary_fsp` and the SciPy-based readouts in
`crn_craft.likelihood`. `import crn_craft` works without SciPy; these
functions then raise `ModuleNotFoundError` naming the extra.

The native core lowers the model and evaluates the rates, including
count-dependent rates and stop conditions. Python enumerates the states
and applies interventions. SciPy propagates the law: sparse matrices,
the uniformization series, `solve_ivp` BDF and `expm_multiply`. The
existing `Model.fsp` is unchanged and fully native. Use it for one fixed
model from one initial state. Use `prepare_fsp` for protocols, prescribed
signals, initial mixtures, checkpoints or many parameter vectors on one
graph.

Prepare with `initial={(0,): 0.25, (1,): 0.75}` to seed a distribution.
The labels are count tuples in the order of `projection.species`. The
probabilities must be finite, nonnegative and sum to one within 1e-12.
Preparation enumerates every transition whose reactants are present,
whatever its rate, and every state an action can reach. So a rule with
zero rate at preparation can be switched on later. `solve(parameters=...)`
reuses this graph and only re-evaluates the rates. A binding that changes
the topology or the initial counts raises `ValueError`; prepare a new
projection for new bounds or initial support. A prepared projection
keeps no per-solve state, so threads can share it. It cannot be pickled;
prepare one per process.

An action maps the retained law through its count kernel, one species at
a time; a partition is binomial. Escaped probability stays in the
absorbing sink; an action can add to it but never returns it. Previously stopped
probability carries forward. A snapshot from `law.snapshots` continues a
solve under the same parameters and signals. To change parameters at a
time point, use a protocol action.

Pass `signals={"arrival": crn.Signal.table(times, rates)}` to replace
the count-independent rate constant of a named rule by a prescribed time
series in the rule's rate units and the model's time unit. Tables use
linear interpolation or an explicit right-continuous step curve and hold
their end values outside the table. Step times split the integration;
the knots of a linear table do not. `Signal(function, change_times=...)`
accepts a Python rate function; list its discontinuities in
`change_times`. Python callables are supplied by the caller when
rerunning the experiment.

Constant segments use uniformization with a truncated Poisson series.
Smoothly forced segments use SciPy's sparse BDF (`solver="bdf"`) or an
adaptive exponential midpoint rule (`solver="midpoint"`). `result.solver`
names the method used. Each snapshot reports the absorbing sink, the
Poisson series-tail bound, the mass defect and the minimum probability.
The series-tail bound covers uniformization only. It excludes rounding
and adds nothing for BDF or midpoint segments, which report no
integration error estimate. BDF can return slightly negative
probabilities; check `minimum_probability`. The midpoint rule is much
slower and stops with an error when its tolerance is below the rounding
level of its step estimate. Check tolerance refinement on the selected
observable for a smooth signal. See the
[RNA fit replay](../../examples/yeast_rna_fit/native_fit.py).

Use `stationary_fsp(model, bounds)` for one finite closed communicating
class or specify `class_state`. With a `return_state`, escaped rates
redirect to a retained state. SciPy's sparse LU solves the linear
system, and the result reports its L1 residual. Successive enlargements
compare marginal CDFs and outward flux. These diagnostics are empirical
on unbounded networks. Certified error bounds require additional
recurrence evidence. `result.method` separates the exact finite law
(`finite-closed-class`) from the redirected approximation
(`stationary-fsp-return-state`).
[Stationary FSP method](https://arxiv.org/abs/1704.07259)

## Spatial models

`SpatialModel(name, volumes)` expands named, well-mixed compartments to an
ordinary count model. `volumes` maps region IDs (letters, digits and `_`)
to positive, fixed volumes. Each species gets one count per region,
labeled `<species>_<region>`; `label(species, region)` returns that name
and `profile(run.paths, species)` selects the region columns in region
order. Ensemble means of a profile are the means of those columns. There
is no separate helper for conserved totals; sum the region columns.

- `add_species(name, initial, reservoirs=())` takes one count for every
  region or a mapping by region.
- `local(name, rate, reactants, products, combinatorial=False,
  regions=None)` adds one reaction per region. A string rate may name the
  local species as `{X}` and the region volume as `{volume}`. Volumes are
  not applied automatically: write a second-order rate constant as, for
  example, `"k / {volume}"`. A local reaction may not change a reservoir
  species, so pass `regions` to leave reservoir regions out.
- `transport(species, origin, target, rate)` adds a directed hop at
  `rate` per copy and day. Volumes do not enter the rate. For detailed
  balance between unequal volumes choose rates with
  `V_i * k_ij = V_j * k_ji`.
- `diffusion(species, first, second, diffusivity, spacing)` adds both
  directions at `diffusivity / spacing**2` per copy. It requires equal
  volumes.

A region without an edge is a reflecting end. A reservoir is a species
whose region count never changes. A hop out of it fires at
`rate * initial count` and adds one copy to the target. A hop into it
removes the copy from the origin; the copy is absorbed.

`spatial_turnover(initial, end=, target=, feedback=, mutant_weight=,
death=, hop=, seed=0, mutation=0.0, fixed_left=False, labels=False,
max_events=3_000_000_000)` samples the endpoint of one exact path of the
two-genotype chain of the [spatial spread study](../guide/19-mtdna-spread.md).
`initial` lists `(wild, mutant)` counts by region, at most 4,096 regions.
Every copy is born at `max(0, death + feedback * (target - W - mutant_weight * M))`,
dies at `death` and hops to each nearest neighbor at `hop` per day. With
`mutation > 0` a birth makes a new mutant clone with that probability.
With `labels=True` the result lists the clone ID of every mutant copy by
region: 0 for initial mutants, then 1, 2, … in order of creation. Labels
do not change the hazards. With `fixed_left=True` region 0 is a reservoir
with the same semantics as above. The result holds `wild`, `mutant`,
`labels`, `events`, `mutations` and the sampler identity
`chacha8-u53-ziggurat/spatial-turnover/2`. Seed `s` uses the same ChaCha8
stream as replicate 0 of a generic run with root seed `s`; do not use
equal seeds for a comparison between the two. The sampler raises `CrnError`
with `CRN-R030` for invalid input, `CRN-R031` for a positive weight below
16 × 2⁻⁵³ of its total (the resolution of the selection draw), `CRN-R040`
when the event budget is reached, and `CRN-R044` when the run is cancelled.

The sampler keeps the region weights in a binary sum tree that recomputes
each parent from its children, so no running total drifts. It updates the
event's region and the destination of a hop. Against the generic engine on
a four-region chain with the clamp and the reservoir, five count CDFs over
20,000 paths each differed by at most 0.0099, against a simultaneous 99%
bound of 0.0276. On the 500-region chain it took about 150 ns per event,
4 to 7 times less than Next Reaction. Use it through
`study.py --source-chain --wave-sampler local`; fiber `i` uses seed
`seed + 1000 + i` and its own checkpoint.

## Observation models

`crn_craft.likelihood` maps states to counts, sums, joint measurements
and intervals. It implements binomial detection and binned Gaussian
fluorescence. Fluorescence bins integrate the Gaussian over the bin
edges. `histogram_score` returns the negative log likelihood in nats per
cell and the contribution of each group. Empty bins contribute zero. An
occupied bin with zero retained probability makes the score infinite
and is counted in `impossible_bins`. The optimistic score adds the
missing probability of a group to each occupied bin. `mean_score` takes
a symmetric positive definite covariance. The RNA example binds
candidates to the prepared projection and supplies the score to fcmaes.
The example's `legacy_tail` argument reproduces the guide's coarse
overflow convention.

## Published experiment files

Use `crn.read_experiment(path, convention="stochastic_count_propensity")`
and then `preflight()` or `run()`. `path` is a COMBINE/OMEX archive
(`.omex` or `.zip`) or a SED-ML file whose models lie in its folder. The
wheel also installs
`crn-craft-experiment EXPERIMENT --convention stochastic_count_propensity`,
with `--check` (preflight only), `--seed`, `--replicates`, `--threads` and
`--archive DIRECTORY`. `DIRECTORY` must not exist; it receives one sealed
native run per task, named by the task ID, and `experiment.json` with the
source member hashes, parameter changes, effective rates, events and
output mappings; `preflight()` returns the same without simulating. Task `i`
in document order uses seed `seed + i`. A refused input raises
`ImportFailure` (`CRN-E090`); the command prints `CRN-E090: message` and
exits with status 2.

| Format | Supported subset |
|---|---|
| SBML L3V1 and L3V2 Core | Fixed compartments of any size and dimension; species with item amount units, as amounts or as concentrations in a sized compartment; an explicit time unit; acyclic scalar parameter assignment rules, recalculated after each edit; time triggers `time >= t` with a literal `t`, instantaneous assignments to non-constant independent parameters or species counts, literal priorities |
| SED-ML L1V3 and L1V4 | Local SBML models and `#model` references; `changeAttribute` of a parameter `value` given as the absolute path `/sbml/model/listOfParameters/parameter[@id='…']/@value` (namespace prefixes allowed and not resolved); uniform time courses with Gillespie direct (`KISAO:0000029`) and no algorithm parameters; repeated tasks with one subtask, one uniform range and `resetModel="true"`, which add replicates; data generators that select one time or species variable; reports and 2D plots |
| COMBINE/OMEX | One master SED-ML document from the manifest; members read from the archive, never extracted |

**Rate convention.** A kinetic law is read as the jump hazard in items per
model time unit, with species amounts in items and concentrations as items
per compartment size. A law must equal `c × Π C(n, k)` for the reactants of
its reaction, the convention of every native rate, so `k*A*(A - 1)` for
`2 A → …` is accepted with `c = 2k`, and `k*A*A` is refused. The explicit
`split_signed_mass_action` convention also accepts a law
`u × Π C(n_reactant, k) − v × Π C(n_product, k)` with `u, v ≥ 0`, found by
exact rational arithmetic; it becomes a forward and a reverse jump, both
nonnegative in every state. This handles the two signed promoter-binding
laws in the [published oscillator archive](https://github.com/biosimulators/Biosimulators_test_suite/blob/ae6afe13c37aa4b896eb5f839ab0e26ba3bc2f5b/examples/sbml-core/Vilar-PNAS-2002-minimal-circardian-clock-discrete-SSA.omex).
Every law of that model is multiplied by the size of compartment `cell`
(4.18879), so its rates are those of the paper times 4.18879 per archive
time unit. Use the strict count convention for ordinary irreversible laws.

**Grid and events.** `numberOfPoints` (L1V3) or `numberOfSteps` (L1V4)
counts intervals: 400 gives 401 report times from `outputStartTime` to
`outputEndTime`. The model starts at `initialTime`. Reports at an event
time show the state after the event. An event whose trigger is true at
time 0 with `initialValue="true"` never fires; with `initialValue="false"`
it fires at time 0. Events at the same time need distinct priorities and
run from the highest; their assignments read the state before that time
unless `useValuesFromTriggerTime="false"`. State triggers, delays,
assignments that read species, and events before `initialTime` are
refused.

**Input safety and limits.** XML must be UTF-8 and may not contain a
document type or entities. Model sources must be relative paths inside the
archive or the SED-ML folder; `..`, absolute paths, URLs, URNs and links
are refused. Limits: 16 MiB per XML document, 64 MiB per archive and for
its declared uncompressed members, 256 members, 64 species, 4,096
parameters and compartments, 1,024 reactions, 256 events, and 16,777,216
returned values (paths × report times × species) per run.

Not supported: SBML Level 2 (the official `L1V3_repeated-stochastic-runs.omex`
is refused), local functions, algebraic and rate rules, other KiSAO
methods, seeds or other algorithm parameters, data generator expressions,
`setValue` in repeated tasks, and report labels (the series are keyed by
data set ID). Unsupported constructs and unresolved units produce a
diagnostic during translation. [SBML](https://sbml.org/documents/specifications/level-3/version-2/core/) · [SED-ML](https://sed-ml.org/specifications.html)

**Verification.** `validation/check_paper_capabilities.py` runs offline in
CI: a two-compartment pulse with a non-unit volume against the same native
protocol (32 paths, identical counts), the pair hazard (survival 0.3535 at
4,000 paths against e^−1 = 0.3679), derived rates after an edit, event
timing, and 33 refused inputs. CI also downloads the oscillator archive and
SBML Test Suite stochastic case 00001 by commit and checks their SHA-256.
Both oscillator tasks give 401 report times; the decay rate of `R` is
0.837758 and 0.209440 per archive time unit. `validation/check_public_experiments.py`
compares the counts of `A` and `R` at times 0.5, 1 and 2 with a separately
written direct SSA of the same reactions: at 512 paths per task the largest
of 12 two-sample CDF gaps was 0.066 against a simultaneous bound of 0.123,
which detects an edit that is not applied but not a 10% change of one rate.
With `--vilar-paths 8704` (about 3 minutes on one core) the largest gap was
0.023 against 0.030. In a review run with a reference SSA changed on
purpose, 512 paths missed a 10% change of the `R` decay rate or of the
promoter unbinding rate (largest gaps 0.086 and 0.057), and 8,704 paths
caught both (0.062 and 0.038). Case 00001 at 5,000 paths and 50
report times gave |Z| ≤ 1.71 and |Y| ≤ 2.73, within the suite's ranges of 3
and 5. These runs exercise the two 2-unit tasks of the deposited
experiment, not the paper's long oscillation.
