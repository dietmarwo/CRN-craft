# Three-class mtDNA ensembles

The rules engine represents wild type, neutral mutants and advantaged mutants
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
from crn_craft import rules

model = rules.load("models/rules/mtdna-three-class.rules.json")
run = rules.ensemble(
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
print(result.records["rule_events"].sum(axis=0))
```

Python exposes read-only NumPy record columns, an ordered rule-count
matrix, report-time statistics and bounded observation results. Termination
codes index `rules.TERMINATIONS`. Ctrl-C and cancel files use cooperative
cancellation. Interrupted ensembles withhold statistics. A prepared
ensemble can be executed again with independent execution state.

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

Rules model and run version 0.1 still select the tree engine and the
original sampler. Use rules-run/0.2 for a `method` object, `trackers`
or `event_log`. It returns rules-result/0.2 with a sampler identity and
optional ordered `observations`. Each backend has its own replay identity.
Worker count preserves a completed numerical payload apart from its
request hash.

| Method | Scope |
| --- | --- |
| `tree` | Exact Direct over dynamic attributed instances |
| `static_direct` | Exact Direct over a closed set of tuples and fixed channels |
| `mtdna_direct` | Exact grouped Direct for the fixed three-class vacancy model, with Static Direct fallback |
| `next_reaction` | Gibson–Bruck indexed clocks over static channels |
| `tau` with `epsilon` | Approximate adaptive Poisson leaps with competing critical events and exact fallback |
| `binomial` with `step` | Approximate death-then-birth splitting for single-parent unit reactions |
| `thinning` | Exact stochastic thinning with analytic ATP flow and a deterministic zero crossing |

Static lowering evaluates local predicates and factors once and compiles
constant and vacancy-rate kernels. It accepts up to 256 tuples and 4096
channels. Static Direct and leap fallbacks select among all channels in one
level, so every positive channel must be at least 2⁻⁵³ of the total
(`CRN-R031`); the tree engine resolves rules and instances separately and can
run such models. Random product attributes use the tree engine. Identical pairs
in the static path require a single matching tuple. The network engine
also exposes Gibson–Bruck under ssa-run/0.4, inverse waiting time and
checked policy. Other SSA schemas retain their existing defaults.

### Specialized fixed-population Direct

Use `{"schema": "crn-craft.rules-run/0.3", "method": {"kind": "mtdna_direct"}}`
in the corresponding request fields, or pass `method="mtdna_direct"` to
`rules.ensemble`. The result uses the existing rules-result/0.2 shape.
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
`tree`. ATP requires `thinning`.

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
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
template = rules.ensemble(
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

`result.event_columns(replicate)` exposes time, rule index and observable
arrays for the bounded log. `result.tracker_columns(index)` exposes the
scalar tracker summaries as arrays, with NaN for missing times or peaks.
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

Existing `static_direct` requests keep inverse-transform waits, the same
sampler identity and per-event tracker integration. Rules-run/0.4 accepts
two independent opt-ins in the method object.

```json
{
  "kind": "static_direct",
  "waiting_time": "ziggurat",
  "sparse_trackers": true
}
```

The Python rules helper selects rules-run/0.4 for this method dictionary.
Rust callers constructing `Method::StaticDirect` use the struct variant
with `waiting_time: None, sparse_trackers: None` for the previous defaults.
The result still uses rules-result/0.2. Ziggurat changes the seeded stream
consumption. Sparse tracking integrates each eligible unweighted single-class
count over its constant spans. The integral can differ in the last bits
from per-event integration. First-passage predicates and other projection
families retain per-event observation. If no tracker is eligible, the run
reports the default sampler identity and a `CRN-R051` warning. Attaching a Continue observer or an
event log leaves the opt-in tracker policy consistent.

On a two-CPU host the opt-in path took about 100–110 ns per event on the
128-cell probe, against about 35 ns for `mtdna_direct`. Inline factors, a
fused resolution check, exact incremental projections and a grouped
evaluator did not change that by 5%, so none was adopted
([decision record](decisions/generic-rules-structural-2026-11-03.md)).
From Rust, a generated event loop runs the same requests with
byte-identical results at about 49 ns per event
([generated event loops](specialized-kernels.md#generated-rules-event-loops)).
`mtdna_direct` remains the fastest route for the fixed three-class topology
and the only fast one from Python or the CLI.

Cached rates and dependency-aware projection refresh apply to the default
path. Potentially overflowing report-only projections retain eager checks.
See [the rules kernel review checklist](specialized-kernels.md) for the
conditional Rust generator and its workload-specific adoption gate.
