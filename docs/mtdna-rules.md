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
| `next_reaction` | Gibson–Bruck indexed clocks over static channels |
| `tau` with `epsilon` | Approximate adaptive Poisson leaps with competing critical events and exact fallback |
| `binomial` with `step` | Approximate death-then-birth splitting for single-parent unit reactions |
| `thinning` | Exact stochastic thinning with analytic ATP flow and a deterministic zero crossing |

Static lowering evaluates local predicates and factors once and compiles
constant and vacancy-rate kernels. It accepts up to 256 tuples and 4096
channels. Random product attributes use the tree engine. Identical pairs
in the static path require a single matching tuple. The network engine
also exposes Gibson–Bruck under ssa-run/0.4, inverse waiting time and
checked policy. Other SSA schemas retain their existing defaults.

Tau follows the species drift and variance step selection of Cao,
Gillespie and Petzold. The expression-rate heuristic uses g=2 and falls
back to exact steps for a submolecule allowance. Low populations,
including creation into a rare class, are treated as critical.
Unacceptable proposals halve the step before an exact fallback.
Binomial splitting recomputes the replication rates after degradation
and groups competing replication outcomes per parent. Java's hourly
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
reports ATP on the grid and stops at its zero crossing. The additional
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
