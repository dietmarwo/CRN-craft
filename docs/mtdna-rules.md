# Three-class mtDNA ensembles

For a task-oriented introduction, read the [user manual](mtdna-user-manual.md),
the [optimization tutorial](mtdna-optimization-tutorial.md) and the
[performance analysis](mtdna-performance-analysis.md). This page is the
reference behind them.

CRN-craft represents wild type, neutral mutants and advantaged mutants
with three fixed attribute tuples. The model
`models/rules/mtdna-three-class.rules.json` starts at 1000 wild types. It
uses degradation ln(2)/240 per hour and replication proportional to
`max(0, K - N) / (wt + sm * (mt + mtAdv))`. Mutation is split into two
explicit channels. The horizon used in validation is 25,920 hours.
Set `stop` to `N > 0 && mutants >= 0.9 * N` for the Java fixed-model
takeover convention. The clone fixture `mtdna-fixed.rules.json` retains
lineage identity through its `clone` attribute.

## Python and event observation

```python
import crn_craft as crn

model = crn.load("models/rules/mtdna-three-class.rules.json")
run = crn.ensemble(
    model, end=25920.0, report_times=range(0, 25921, 720),
    replicates=10_000, threads=20, seed=17, method="static_direct",
    trackers=[{
        "observable": "mutants",
        "first_passage": "N > 0 && mutants >= 0.9 * N",
        "keep_episodes": 8,
    }],
)
result = run.simulate()
print(result.metadata["termination"])
print(result.rule_events.sum())  # events per rule, summed over replicates
passage = result.tracker_columns(0)["first_passage"]  # NaN: never reached
print(passage.dropna().describe())
```

The blocks use the top-level names, which write `crn-craft.model/1.0` and
`crn-craft.run/1.0` documents; `crn_craft.rules` is the same API under its
0.6.0 name and writes the `rules-*` versions
([user manual](mtdna-user-manual.md#18-reference)). For this model both
run the same Static Direct sampler: the default selector `auto` picks the
linear scan here.

Since 0.3.0 Python results are pandas objects (pandas is a required
dependency): `result.records` is a DataFrame indexed by replicate id with a
categorical `termination` (compare with a name such as `"stopped"`),
`result.rule_events` has one column per rule name, `result.statistics` holds
one DataFrame per statistic indexed by report time with one column per
observable. Since 0.5.0 `result.paths` holds the kept series as an xarray
DataArray, with dimensions (`replicate`, `time`, `variable`) since 0.6.0,
and xarray is a required dependency; `result.paths_frame()` returns the DataFrame indexed
by (replicate, time) of earlier releases, and `result.to_xarray()` a
Dataset. `keep_paths="all"` keeps every replicate. Every frame and array has
`.to_numpy()`; `result.metadata` keeps the full result document, and
`run.rule_names` and `result.rule_names` name the rules in declaration
order. The [user manual](mtdna-user-manual.md#18-reference) lists the
changes from the NumPy results of 0.2.0. Ctrl-C and cancel files use
cooperative cancellation. Interrupted ensembles withhold statistics. A
prepared ensemble can be executed again with independent execution state.

Trackers capture zero crossings, last crossing times, extrema, time
integrals and first passage. A tracked interval runs until the replicate's
stop or horizon. A predicate that evaluates to NaN does not count as a
hit. Episodes keep the configured prefix and summarize omitted duration
and peak. `event_log={"max_entries": 100, "rules": ["mutation_mt"]}`
retains up to 100 selected events per path and counts omissions. Logs
include rule index, participant type and attributes, and post-event
observables. Retention counts enter admission checks, so use small cohorts
for long logs. The three-class model observes class appearances and
extinctions. Per-clone fates can be reconstructed from the clone fixture
and its complete event log for a small cohort.

`simulate(callback=fn)` calls `fn(replicate_id, event)` on native
workers. Returning `False` stops that replicate. Each path preserves
event order. Callbacks acquire the GIL, and interleaving across replicates
is scheduler-dependent. Exceptions cancel the run and propagate. A
`batch_callback` receives chunks of the retained log after simulation.
The Rust `execute_with_observer` entry point exposes the same committed
events and supports a per-path stop. The default observer is compiled
away.

## Algorithms and contract migration

`crn-craft.model/1.0` and `crn-craft.run/1.0` mean what rules-model/0.4 and
rules-run/0.6 mean, plus the features that arrived with them: Static Direct
selectors, quantile statistics, the `fsp` method, higher-order rules and
larger limits ([model reference](rules-language.md)). The history of the
older versions follows. Rules model and run version 0.1 still select the
tree method and the original sampler. Use rules-run/0.2 for a `method` object, `trackers`
or `event_log`. It returns rules-result/0.2 with a sampler identity and
optional ordered `observations`. Rules-run/0.3 adds `mtdna_direct` and 0.4
the Static Direct options (below); both return rules-result/0.2.
Rules-run/0.5 adds `hybrid` with its tolerances and returns
rules-result/0.3, which adds the request's `method` and the names of the
`continuous` columns. Rules-run/0.6 adds no field: Static Direct defaults
to Ziggurat waits and sparse trackers (below), and the result is
rules-result/0.3. Rules-model/0.2 adds the legacy `atp` block, and
rules-model/0.3 whole-rule `propensity` and `continuous` variables; from
0.3 a rule gives exactly one of `rate` and `propensity`. Rules-model/0.4
adds `"independent": true` for continuous variables whose derivative reads
only parameters and other independent variables; `hybrid` integrates them
once per prepared run and shares the solution by every replicate
([rules language](rules-language.md#independent-variables-rules-model04)).
A 0.4 document without the field means what it means as 0.3. Stored
rules-model 0.1–0.3 files and rules-run 0.1–0.6 requests keep their
meaning and replay byte-identically. Each backend has its own replay identity. Worker count
preserves a completed numerical payload apart from its request hash.

| Method | Scope |
| --- | --- |
| `tree` | Exact Direct over dynamic attributed instances |
| `static_direct` | Exact Direct over a closed set of tuples and fixed channels |
| `mtdna_direct` | Exact grouped Direct for the fixed three-class vacancy model, with Static Direct fallback |
| `next_reaction` | Gibson–Bruck indexed clocks over static channels |
| `tau` with `epsilon` | Approximate adaptive Poisson leaps with competing critical events and exact fallback |
| `binomial` with `step` | Approximate death-then-birth splitting for single-parent unit reactions |
| `thinning` | Exact stochastic thinning with analytic ATP flow and a deterministic zero crossing |
| `hybrid` | Approximate piecewise-deterministic simulation of continuous variables (rules-model/0.3; independent variables once per run in 0.4) at declared tolerances; rules-run/0.5 |

Higher-order rules of `model/1.0` (a pattern `count` above one, or more than
two patterns) run on `static_direct`, `next_reaction`, `tau` and `hybrid`;
`tree`, `mtdna_direct` and `binomial` refuse them with `CRN-R020` at
`/method`.

Rules-model/0.3 adds whole-rule propensities
([rules language](rules-language.md#whole-rule-propensities-rules-model03)),
which `crn_craft.rules.Model` writes for string rates. Every method above
samples them. The static methods need each reactant type of such a rule to
lower to one tuple, as for attribute-free species, and refuse other models
with `CRN-R020` at `/method`; `tree` has no such limit. `mtdna_direct` does
not certify them and runs Static Direct with a `CRN-R051` warning, and the
generated rates kernels refuse them. Generated event loops emit them, as
they emit every rate expression and observable that Static Direct runs
([checklist](specialized-kernels.md#generated-rules-event-loops)).

Static lowering evaluates local predicates and factors once and compiles
constant and vacancy-rate kernels. It accepts up to 256 tuples and 4,096
channels for the rules-model versions, and 4,096 tuples and 65,536 channels
for `model/1.0`. Static Direct and leap fallbacks select among all channels in one
level, so every positive channel must be at least 2⁻⁵³ of the total
(`CRN-R031`); the tree engine resolves rules and instances separately and can
run such models. Random product attributes use the tree method. Identical pairs
in the static path require a single matching tuple. Under `run/1.0` Static
Direct also has the selectors of [exact methods](direct-ssa.md#selectors).

### Specialized fixed-population Direct

Use `{"schema": "crn-craft.rules-run/0.3", "method": {"kind": "mtdna_direct"}}`
in the corresponding request fields (or any later request version,
`crn-craft.run/1.0` included), or pass `method="mtdna_direct"` to
`crn.ensemble`. The result uses the existing rules-result/0.2 shape.
The sampler `chacha8-u53-ziggurat/rules-mtdna-direct/1` confirms that the
specialization executed. On a fallback the result carries the Static Direct
sampler identity and a `CRN-R051` warning that names the reason.
The same seed gives a different trajectory under the two samplers.
Each sampler preserves its own completed numerical results across worker counts.

Preparation verifies the lowered channels, unit state updates, constant
degradation, common vacancy law, integer capacity up to 1,000,000 and a
positive competition denominator. It certifies that in every reachable state
each positive hazard is at least 16 × 2⁻⁵³ of the total, evaluated at the
worst state (one vacancy, one molecule of the class). It supports no stop or
the guarded mutant fraction predicate shown above. Unsupported static models or parameter
ranges use Static Direct. Models with changing instance identities use
`tree`. The legacy ATP block requires `thinning`, and continuous variables
require `hybrid`.

The kernel tracks the total population directly and computes the common
vacancy factor once per event. It groups degradations and wild-type birth
outcomes into cumulative intervals with each original reaction's hazard.
Report observables are evaluated on demand. Simple single-class count trackers
integrate between changes in that class. Last-bit rounding in integrals
can differ from integration across every event. First-passage expressions,
event logs and live callbacks retain the general observation path.
The unit exponential uses the existing Julia-derived Ziggurat sampler.
ChaCha8 streams, accounting and work budgets remain in use.

A prepared template can bind an independent candidate without reconstructing
the run request in Python. Each candidate validates and lowers the new
parameter values once, since rules parameters can change attributes and
topology. The resulting prepared kernel is reusable for subsequent calls.
A candidate's `model_sha256` and `run_sha256` hash the re-serialized bound
model and request, not a file on disk; even an empty override changes them.

```python
import json
import crn_craft as crn

model = json.loads(crn.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
template = crn.ensemble(
    model, end=25920.0, report_times=range(0, 25921, 720),
    replicates=10_000, threads=8, seed=17, method="mtdna_direct",
    keep_paths=10_000, max_output_bytes=512 << 20,
    trackers=[{"observable": name} for name in
              ("neutral_mutants", "advantaged_mutants")],
)
candidate = template.bind_parameters({"mu": 4e-5, "sm": 1.55})
result = candidate.simulate()
print(result.metadata["sampler"], result.metadata["stopped_replicates"])
```

Tau follows the step selection of Cao, Gillespie and Petzold. The species
drift and variance bound uses g=2 for expression rates and falls back to
exact steps for a submolecule allowance. A propensity condition bounds the
mean and spread of every rate change by ε·a_j, using first differences
a_j(x + ν_k) − a_j(x), so rates that read global observables such as the
vacancy K − N are controlled directly. If one firing moves a rate by more
than ε·a_j, the step is exact. Low populations, creation into a rare class
and channels that ten firings could switch off are critical and fire at
most once per leap. Unacceptable proposals halve the step before an exact
fallback. The cost per step grows with the square of the channel count.
In the fixed-population model at K = 1000 the vacancy holds only a few
places, so this rule runs exact events and tau gives no speedup there; the
review gate at K = 1000 found no multi-firing leaps. Leaps occur when every
rate changes slowly, as in the high-copy immigration/death gate.
Binomial splitting recomputes the replication rates after degradation
and groups competing replication outcomes per parent. Each parent's
outcome probabilities are the competing-risk values (a/Λ)(1 − e^(−Λ·dt/n))
with rates frozen over the step. Steps lie on the fixed lattice k·step;
reports between lattice points hold the last step's state, so the report
grid does not change the path. Frozen rates can put N a few molecules
above K. Java's hourly
implementation groups simultaneous new mutants differently and keeps
clone identities. Both approximate methods report aggregate firings at
leap boundaries, with their own sampler identities. DOI
[10.1063/1.2159468](https://doi.org/10.1063/1.2159468).

Gibson–Bruck updates the fired channel and the channels affected by its
count changes. Active clocks are rescaled when their rates change.
Sparse networks can benefit when few clocks change at each event.
DOI [10.1021/jp993732q](https://doi.org/10.1021/jp993732q).

`models/rules/mtdna-three-class-atp-continuous.rules.json` writes the same
model with a continuous variable `ATP` (rules-model/0.3): its derivative is
the expression `f * wild_type - v1 * ATP - v2 * N`, each coupled rate ends in
`* c / (c + ATP)`, and `lower_bound: 0` with `on_bound: "stop"` keeps the
depletion rule below. The `hybrid` method integrates it with Dormand–Prince
5(4) at the request's tolerances
([rules language](rules-language.md#continuous-variables-rules-model03));
`validation/hybrid/RESULTS.md` compares it with thinning on the legacy file.
The legacy file and `thinning` are unchanged.

The ATP fixture uses rules-model/0.2 with an `atp` block. Between events
the solver integrates dA/dt = f × wt − v1 × A − v2 × N analytically.
Selected replication rules are multiplied by c/(c+A). The minimum ATP
over the fixed-count interval gives the thinning envelope. The solver
reports ATP on the grid and stops at its zero crossing. ATP at zero counts as
depleted only if the forcing f × wt − v2 × N is not positive; with positive
forcing it rises again. An ATP stop and a model `stop` both report `stopped`;
ATP is zero on every report after an ATP stop. The additional
ATP column is continuous. Event trackers are restricted to discrete
observables and predicates that read those observables. The supplied Java
SSA uses the same ATP flow but freezes rates over its sampled waiting
interval.

## Reproduction

```bash
cargo build --release --workspace --locked
python validation/check_mtdna.py --binary target/release/crn-craft --output mtdna-gates.json
python validation/check_rules_python.py --output mtdna-python.json
python validation/benchmark_mtdna.py --binary target/release/crn-craft --output mtdna-timing.json
```

Build and install the wheel before the Python checks and benchmark.
The optional Java comparison accepts `--java`, `--java-classpath` and
`--java-source`. Use JDK 25 and Apache Commons RNG 1.6, compile the supplied
source locally, and supply a classpath containing its compiled classes
and the four Commons RNG jars. The benchmark writes only summaries to the
patch. The original Java material and per-trajectory output stay local.

For leap methods, the sampling-step budget bounds both random words
and leap iterations, including iterations with zero firings.

`result.event_columns(replicate)` returns the bounded log of one replicate
as an xarray DataArray with dimensions (`entry`, `variable`), the
observables after each event, and the coordinates `time`, `rule` (the
index, -1 for a leap) and `rule_name`. Without an argument it returns
every log with dimensions (`replicate`, `entry`, `variable`), padded with
NaN (and `rule` -2) past each replicate's `length`.
`result.event_frame()` returns the DataFrames of earlier releases: one row
per logged event, indexed by `entry` or by (`replicate`, `entry`).
`result.tracker_columns(index)` returns the scalar tracker summaries as a
DataFrame indexed by replicate id, with NaN for missing times or peaks.
Ragged participant identities, leap multiplicities and episodes remain
in `result.observations`.

The ATP envelope limits each interval using the current ATP derivative
and half the replication denominator. The derivative magnitude shrinks
while counts stay fixed. This keeps the envelope tight near rapid changes
and advances the analytic flow at each interval boundary.

For the eight-worker comparison, build the optional phase-timing example
with `cargo build --release --locked -p crn-craft-rules --example benchmark_rules`.
Run `python validation/benchmark_mtdna_gap.py --help` for the baseline,
candidate and Java inputs. The saved evidence is
`validation/mtdna/java-gap.json`. The Java phase probe was built locally
by timing from before pool creation to after all futures were collected.
The end-to-end comparison used the original Java classes.

### Generic Static Direct options

Static Direct has two options in the method object, from rules-run/0.4:

```json
{
  "kind": "static_direct",
  "waiting_time": "ziggurat",
  "sparse_trackers": true
}
```

- `waiting_time`: `"inverse"` inverts one uniform per wait; `"ziggurat"`
  draws it with the Ziggurat sampler, which consumes the seeded stream
  differently, so the numbers differ cell by cell but not in distribution.
- `sparse_trackers`: `true` integrates each eligible tracker, an
  unweighted single-class count without a first-passage predicate, over
  its constant spans. The integral can differ in the last bits from
  per-event integration. First-passage predicates and other projection
  families keep per-event observation. Without an eligible tracker,
  `true` changes nothing and reports the per-event sampler identity with a
  `CRN-R051` warning.

What an absent option means depends on the request version:

| Request | Absent `waiting_time` | Absent `sparse_trackers` |
| --- | --- | --- |
| rules-run/0.2 to 0.5 | `"inverse"` | `false` |
| rules-run/0.6 and `crn-craft.run/1.0` | `"ziggurat"` | sparse where a tracker is eligible, otherwise per-event integration without a warning |

Under `run/1.0` the method object also takes `selector`; absent means
`auto`, which picks the linear scan for this model.

So `{"kind": "static_direct"}` in rules-run/0.6 runs
`chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1` with eligible
trackers and `chacha8-u53-ziggurat/rules-static-direct/1` without, and
`{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers":
false}` runs the earlier default with its numbers. The Python rules helper
selects rules-run/0.6 for every `static_direct` method. Rust callers
constructing `Method::StaticDirect` with `waiting_time: None,
sparse_trackers: None, selector: None` get the defaults of their request's
version.
Attaching a Continue observer or an event log leaves the tracker policy
consistent. `mtdna_direct`'s fallback keeps inverse waits and per-event
trackers in every version.

Measured on releases 0.5.0 and 0.6.0, on a two-CPU host, Ziggurat waits
with sparse trackers took about 100–110 ns per event on the 128-cell probe,
against about 35 ns for `mtdna_direct`; on the 1000-cell three-class case they were 1.17× faster
than inverse waits with per-event trackers
([performance analysis](mtdna-performance-analysis.md#2-methods)). Inline
factors, a fused resolution check, exact incremental projections and a
grouped evaluator did not change that by 5%, so none was adopted
([decision record](decisions/generic-rules-structural-2026-11-03.md)).
From Rust, a generated event loop runs the same requests with
byte-identical results at about 48 ns per event, whether the stop
is written in the fused fraction form or not
([generated event loops](specialized-kernels.md#generated-rules-event-loops)).
From Python, `crn.compile_loop` builds that loop with a local Rust
toolchain and runs it with `compiled=`, about 2.2× faster than interpreted
Static Direct with its defaults, end to end, as measured on 0.6.0. Release
0.7.0 made interpreted Static Direct faster (0.72 times the time on the
[three-class workload](mtdna-performance-analysis.md#release-070)). `mtdna_direct` remains the fastest route for the
fixed three-class topology and the only fast one from the CLI.

Cached rates and dependency-aware projection refresh apply to the default
path. Potentially overflowing report-only projections retain eager checks.
See [the rules kernel review checklist](specialized-kernels.md) for the
conditional Rust generator and its workload-specific adoption gate.
