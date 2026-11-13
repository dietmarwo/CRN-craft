# Methods and their options

This page describes every simulation method of the one engine: the exact
methods and their event selectors, the approximate methods, the fast path
`mtdna_direct`, the ATP block, the hybrid method and the compiled event
loops. It also has the options of Static Direct, the random stream and the
numerical policy. The model language is in the [model reference](model-language.md),
the request and model versions are in the [contracts](contracts.md), and the
cost of each method is in the
[performance analysis](../evidence/mtdna-performance-analysis.md). The guide
teaches the choice of a method in
[chapter 7](../guide/07-methods.md). An exact method samples the model's
process with finite-precision arithmetic and Monte Carlo uncertainty. It is
not a biological validation claim.

## The methods

| `method` | Exact | Scope | Sampler identity (prefix `chacha8-u53/` or `chacha8-u53-ziggurat/`) |
| --- | --- | --- | --- |
| `static_direct` | yes | closed models: a fixed set of attribute vectors, any number of reactant molecules in `model/1.0` | `rules-static-direct`, with a suffix for the selector and the trackers |
| `next_reaction` | yes | the same models | `rules-next-reaction/1` |
| `tree` (the default) | yes | every model without an ATP block or continuous variables; instances may change identity | `rules-direct/2` |
| `mtdna_direct` | yes | the fixed three-class vacancy model, with a fallback to Static Direct | `rules-mtdna-direct/1` ([the fast path](#the-fast-path-mtdna_direct)) |
| `thinning` | yes | models with an ATP block, and only those | its own identity ([ATP and thinning](#the-atp-block-and-thinning)) |
| `tau` with `epsilon` | no | closed models: adaptive Poisson leaps with competing critical events and exact fallback | its own identity ([leaping](#approximate-methods-tau-and-binomial)) |
| `binomial` with `step` | no | closed models of single-parent unit reactions: death-then-birth splitting | its own identity |
| `hybrid` | no | models with continuous variables, and models without; piecewise-deterministic at declared tolerances | `rules-approximate-hybrid-dopri54/1` ([hybrid](#the-hybrid-method)) |

Static Direct lowers the model to channels once, at preparation: species,
channels, the changes each channel makes and the channels whose rates each
change affects. Each event then recomputes only the affected rates, draws a
waiting time and selects a channel. A zero total rate is absorbing: the state
is held to the horizon without consuming random numbers. Report values are
right-continuous, so an event exactly at a report time is included.

Three rules hold for every selector:

- A positive channel must be at least 2⁻⁵³ of the total rate. Below that, a
  53-bit uniform cannot select it exactly, and the replicate fails with
  `CRN-R031` rather than dropping the channel.
- A non-finite or negative rate fails the replicate (`CRN-R030`). There is no
  fallback to the last channel.
- Under `run/1.0` a positive wait too small to move the clock is accepted;
  1,000,000 such steps in a row fail with `CRN-R030`. Every `rules-run`
  version fails at the first one.

## Selectors

Under `crn-craft.run/1.0`, `method.selector` of `static_direct` chooses how
the next channel is selected. Every selector samples the same law; each uses
the random stream differently, so the same seed gives different paths, and the
sampler identity names the selector.

| `selector` | How it selects | Sampler suffix |
| --- | --- | --- |
| `linear` | a scan of the channels in order, as in every `rules-run` version | none |
| `composition_rejection` | channels grouped by the binary exponent of their rate; a group by its mass, a member by rejection | `-cr` |
| `sorting` | a scan in an order that moves each fired channel one place forward | `-sorting` |
| `rssa_cr` | candidates by composition–rejection from upper bounds of the rates over population intervals, accepted against lower bounds or the exact rate | `-rssa-cr` |
| `auto` (the default in `run/1.0`) | one of the others, chosen at preparation | the chosen one |

`auto` chooses RSSA-CR from 8 channels when no rate reads two or more
species, and from 128 channels when the species that start with molecules
average fewer than 8 molecules. Otherwise it chooses composition–rejection
from 128 channels, and otherwise the linear scan. It never chooses `sorting`.
It decides from the initial counts: a model whose populations later grow or
shrink far from them can end up with a slower selector. The thresholds are in
`describe()["selectors"]["auto"]`.

The `rules-run` versions run the linear scan and refuse the field
(`CRN-R020` at `/schema`). A request that the Python helpers write as
`run/1.0` for a `crn_craft.rules` call names `"selector": "linear"`, so
adding quantiles never changes the paths. The top-level `crn_craft.ensemble`
leaves the choice to `auto`. Generated event loops and rates kernels run the
linear scan only, and refuse any other selector, also when `auto` would pick
one (`CRN-R020` at `/method/selector`).

Selectors are ports of published algorithms: composition–rejection (Slepoy,
Thompson and Plimpton, J. Chem. Phys. 128, 205101 (2008)), sorting direct
(McCollum et al., Comput. Biol. Chem. 30 (2006), 39–49) and RSSA-CR
(Thanh, Priami and Zunino, J. Chem. Phys. 141, 134116 (2014), with
composition–rejection for the candidate draw).

## Waiting times and the options of Static Direct

Static Direct has two options in the method object, from rules-run/0.4:

```json
{
  "kind": "static_direct",
  "waiting_time": "ziggurat",
  "sparse_trackers": true
}
```

- `waiting_time` is `"inverse"` or `"ziggurat"`. The inverse transform inverts
  one uniform per wait. Ziggurat draws the wait with the Ziggurat sampler,
  which consumes the seeded stream differently, so the numbers differ cell by
  cell but not in distribution.
- `sparse_trackers` is `true` or `false`. `true` integrates each eligible
  tracker, an unweighted single-class count without a first-passage
  predicate, over its constant spans. The integral can differ in the last
  bits from per-event integration. First-passage predicates and other
  projection families keep per-event observation. Without an eligible
  tracker, an explicit `true` changes nothing and reports the per-event
  sampler identity with a `CRN-R051` warning.

What an absent option means depends on the request version:

| Request | Absent `waiting_time` | Absent `sparse_trackers` |
| --- | --- | --- |
| rules-run/0.2 to 0.5 | `"inverse"` | `false` |
| rules-run/0.6 and `crn-craft.run/1.0` | `"ziggurat"` | sparse where a tracker is eligible, otherwise per-event integration without a warning |

Under `run/1.0` the method object also takes `selector` (above); absent means
`auto`. So `{"kind": "static_direct"}` in rules-run/0.6 runs
`chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1` with eligible
trackers and `chacha8-u53-ziggurat/rules-static-direct/1` without, and
`{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers":
false}` runs the earlier default with its numbers. The Python rules helper
selects rules-run/0.6 for every `static_direct` method. Rust callers that
construct `Method::StaticDirect` with `waiting_time: None,
sparse_trackers: None, selector: None` get the defaults of their request's
version. Attaching a Continue observer or an event log leaves the tracker
policy consistent. The fallback of `mtdna_direct` keeps inverse waits and
per-event trackers in every version.

Ziggurat draws exponential waits from tables adapted from the Julia ecosystem
(see [THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md)). The inverse
transform is `-ln(1 - u) / total`, with `u` a 53-bit uniform in [0, 1) that is
drawn again when it is 0.

## The random stream

Every replicate draws from its own ChaCha8 stream. The key is the SHA-256 of
the domain `crn-craft.rules/0.1` and the root seed; the replicate
index is the stream id. A 53-bit uniform is the high 53 bits of a 64-bit word
times 2⁻⁵³. Results therefore depend on the seed and the method, never on the
thread count, and a replicate's path is the same in a larger ensemble.

Statistical equivalence of the selectors is checked against exact
distributions: master-equation solutions and the finite state projection, with
predeclared tests and a recorded family error rate (`validation/merge/`).

## Next reaction

Gibson–Bruck next reaction updates the fired channel and the channels that
its count changes affect. Active clocks are rescaled when their rates change.
Sparse networks can benefit when few clocks change at each event (DOI
[10.1021/jp993732q](https://doi.org/10.1021/jp993732q)).

## Approximate methods: tau and binomial

Tau follows the step selection of Cao, Gillespie and Petzold (DOI
[10.1063/1.2159468](https://doi.org/10.1063/1.2159468)). The species drift and
variance bound uses g=2 for expression rates and falls back to exact steps for
a submolecule allowance. A propensity condition bounds the mean and spread of
every rate change by ε·a_j, using first differences a_j(x + ν_k) − a_j(x), so
rates that read global observables such as the vacancy K − N are controlled
directly. If one firing moves a rate by more than ε·a_j, the step is exact.
Low populations, creation into a rare class and channels that ten firings
could switch off are critical and fire at most once per leap. Unacceptable
proposals halve the step before an exact fallback. The cost per step grows
with the square of the channel count. In the fixed-population mtDNA model at
K = 1000 the vacancy holds only a few places, so this rule runs exact events
and tau gives no speedup there; the review gate at K = 1000 found no
multi-firing leaps. Leaps occur when every rate changes slowly, as in the
high-copy immigration/death gate.

Binomial splitting recomputes the replication rates after degradation and
groups competing replication outcomes per parent. Each parent's outcome
probabilities are the competing-risk values (a/Λ)(1 − e^(−Λ·dt/n)) with rates
frozen over the step. Steps lie on the fixed lattice k·step; reports between
lattice points hold the last step's state, so the report grid does not change
the path. Frozen rates can put N a few molecules above K. The hourly
implementation of the Java reference groups simultaneous new mutants
differently and keeps clone identities.

Both methods report aggregate firings at leap boundaries, each with its own
sampler identity. For both, the sampling-step budget bounds the random words
and the leap iterations, including iterations with zero firings.

## Higher-order rules and static lowering

`crn-craft.model/1.0` and `crn-craft.run/1.0` mean what rules-model/0.4 and
rules-run/0.6 mean, plus the features that arrived with them: Static Direct
selectors, quantile statistics, the `fsp` method, higher-order rules and
larger limits ([model reference](model-language.md)). The versions are in
[contracts](contracts.md#versions-of-the-request-and-the-model).

Higher-order rules of `model/1.0` (a pattern `count` above one, or more than
two patterns) run on `static_direct`, `next_reaction`, `tau` and `hybrid`.
`tree`, `mtdna_direct` and `binomial` refuse them with `CRN-R020` at
`/method`.

Whole-rule propensities (rules-model/0.3, which `crn_craft.Model` writes for
string rates; see the [model reference](model-language.md#whole-rule-propensities-rules-model03))
are sampled by every method. The static methods need each reactant type of
such a rule to lower to one tuple, as for attribute-free species, and refuse
other models with `CRN-R020` at `/method`; `tree` has no such limit.
`mtdna_direct` does not certify them and runs Static Direct with a `CRN-R051`
warning, and the generated rates kernels refuse them. Generated event loops
emit them, as they emit every rate expression and observable that Static
Direct runs ([checklist](../developer/specialized-kernels.md#generated-rules-event-loops)).

Static lowering evaluates local predicates and factors once and compiles
constant and vacancy-rate kernels. It accepts up to 256 tuples and 4,096
channels for the rules-model versions, and 4,096 tuples and 65,536 channels
for `model/1.0`. Static Direct and the leap fallbacks select among all
channels in one level, so every positive channel must be at least 2⁻⁵³ of the
total (`CRN-R031`); the tree method resolves rules and instances separately
and can run such models. Random product attributes use the tree method.
Identical pairs in the static path require a single matching tuple.

Models with changing instance identities use `tree`. The legacy ATP block
requires `thinning`, and continuous variables require `hybrid`.

## The tree method and its numerical policy

- The tree method is Gillespie's direct method. Rule selection scans the rule propensities in
  order; instance selection descends a pairwise sum tree per distinct pattern.
  Tree totals are recomputed from their children on each update (no
  accumulated drift); their rounding depends on the slot layout, which is a
  deterministic function of the replicate's history.
- Identical-pattern pairs use a pair tree: each node stores its pattern
  weight S and the weight P of pairs of distinct molecules inside it,
  P = P_left + P_right + S_left·S_right (a leaf holds h²·n(n−1)/2). The
  propensity is g × P at the root, a sum of nonnegative terms without the
  cancellation of ½[(Σ h n)² − Σ h² n]. Selection descends by P into a child or
  stops at the cross term, then picks one molecule on each side by S; no
  rejection (sampler `chacha8-u53/rules-direct/2`; version 1, used only for the
  first K1 gate run, drew two molecules and rejected the same one).
- A positive rule propensity or instance weight below 2⁻⁵³ of its total cannot
  be selected exactly by a 53-bit uniform and fails with `CRN-R031`. For pair
  rules the check is conservative: the smallest factor squared must reach
  2⁻⁵³ of the pair total.
- Overflowing weights, pair totals and observable totals fail the replicate
  (`CRN-R030`) rather than becoming 0 or `null`.
- Each replicate draws from its own ChaCha8 stream keyed by SHA-256 of the
  domain `crn-craft.rules/0.1` and the root seed, with the replicate index as
  stream id. Results do not depend on the thread count.
- Counts stay at most 2⁵³ (`CRN-R042` otherwise) and identities below 2⁵³.

## The fast path `mtdna_direct`

`mtdna_direct` is an exact direct method fused for one model structure: three
classes in a fixed population. It reads that structure from the lowered rules;
the model file and its names do not matter. Select it with `method` set to
`{"kind": "mtdna_direct"}` in a request of version `rules-run/0.3` or later
(`run/1.0` included), or with `method="mtdna_direct"` in `crn.ensemble`. The
guide shows it in [chapter 14a](../guide/14a-mtdna-model-files.md).

At preparation it lowers the model to a fixed network, as Static Direct does,
and certifies six points:

1. **Three classes, eight channels.** The network has exactly three attribute
   vectors (classes) and eight channels. Each channel consumes one molecule
   and changes one class by one molecule: every class degrades and
   replicates, and one class, the wild type, also yields one molecule of each
   other class (the two mutation channels). The roles are read from the
   channels.
2. **Degradation.** All three classes degrade at one constant rate.
3. **Replication.** Every replication and mutation rate is `max(0, K - N) / D`
   times constant factors. All five share one `K`, an expression of
   parameters whose value is an integer from 1 to 1,000,000 and at least the
   initial population. `N` is an observable that counts every molecule, and
   `D` an observable that sums the classes with finite positive weights.
4. **Observables.** Every observable is a finite weighted sum over the three
   classes.
5. **Stop.** The stop is absent, or `N > 0 && mutants >= f * N` with a number
   `f` in [0, 1], where `mutants` counts the two classes other than the wild
   type.
6. **Resolution.** Bounds from `K`, the weights and the rate constants show
   that every hazard stays finite, and every positive hazard at least
   16 · 2⁻⁵³ of the total, in every reachable state (evaluated at the worst
   state: one vacancy, one molecule of the class).

The certificate holds for one model and one parameter set: `bind_parameters`
certifies each candidate again. If a check fails, the run uses generic Static
Direct with inverse waits and per-event trackers. `metadata["sampler"]` is
then `chacha8-u53/rules-static-direct/1`, the numbers equal those of
`{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers":
false}`, and a `CRN-R051` warning in `metadata["diagnostics"]` names the
failed check. A model that cannot be lowered at all, such as the clone model,
is rejected with `CRN-R020`, as by every method except the tree. A model
written with the builder is not certified: its rules carry whole-rule
propensities, so the run falls back with `CRN-R051`.

`mtdna_direct` samples the same process as Static Direct with a different
random stream: each cell's ChaCha8 words go to Ziggurat waiting times and to
a selection over grouped hazards. Its sampler identity is
`chacha8-u53-ziggurat/rules-mtdna-direct/1`, which confirms that the
specialization ran. Results agree with Static Direct in distribution, not
cell by cell, and the same seed gives a different trajectory under the two
samplers. Each sampler keeps its own completed numerical results across
worker counts.

How it works: the kernel tracks the total population directly and computes
the common vacancy factor once per event. It groups degradations and
wild-type birth outcomes into cumulative intervals with each original
reaction's hazard, so it computes the eight hazards from three counts with
one division and reads its stop from the counts. Report observables are
evaluated on demand. Simple single-class count trackers integrate between
changes in that class; the integral can differ in the last bits from
integration across every event. First-passage expressions, event logs and
live callbacks keep the general observation path. The unit exponential uses
the Julia-derived Ziggurat sampler. ChaCha8 streams, accounting and work
budgets stay in use.

What it gains, measured on 12 November 2026 on the file form of the model on
the two-CPU host (medians of three runs, `validation/merge/manual-ratios/`):
about 2.8 times less time per event than Static Direct with its default
Ziggurat waits and sparse trackers, and 3.1 times less than Static Direct with
inverse waits and per-event trackers. The [performance
analysis](../evidence/mtdna-performance-analysis.md#2-methods) has older
measurements. It is the fastest route for the fixed three-class topology and
the only fast one from the command line.

A prepared run can bind an independent candidate without writing the request
again: each candidate validates and lowers the new parameter values once,
since the parameters of an attributed model can change attributes and
topology. A candidate's `model_sha256` and `run_sha256` hash the re-serialized
bound model and request, not a file on disk; even an empty override changes
them.

## The ATP block and thinning

The legacy ATP fixture is a rules-model/0.2 document with an `atp` block, and
`thinning` is the exact method that runs it. Between events the solver
integrates dA/dt = f × wt − v1 × A − v2 × N analytically. Selected
replication rules are multiplied by c/(c+A). The minimum ATP over the
fixed-count interval gives the thinning envelope: it limits each interval
using the current ATP derivative and half the replication denominator. The
derivative magnitude shrinks while counts stay fixed, which keeps the
envelope tight near rapid changes and advances the analytic flow at each
interval boundary. The solver reports ATP on the grid and stops at its zero
crossing. ATP at zero counts as depleted only if the forcing f × wt − v2 × N
is not positive; with positive forcing it rises again. An ATP stop and a
model `stop` both report `stopped`; ATP is zero on every report after an ATP
stop. The additional ATP column is continuous. Event trackers are restricted
to discrete observables and predicates that read those observables. The
supplied Java SSA uses the same ATP flow but freezes rates over its sampled
waiting interval.

`models/rules/mtdna-three-class-atp-continuous.rules.json` writes the same
model with a continuous variable `ATP` (rules-model/0.3): its derivative is
the expression `f * wild_type - v1 * ATP - v2 * N`, each coupled rate ends in
`* c / (c + ATP)`, and `lower_bound: 0` with `on_bound: "stop"` keeps the
depletion rule above. The `hybrid` method integrates it with Dormand–Prince
5(4) at the request's tolerances
([continuous variables](model-language.md#continuous-variables-rules-model03));
`validation/hybrid/RESULTS.md` compares it with thinning on the legacy file.
The legacy file and `thinning` are unchanged.

## The hybrid method

Only `hybrid` runs a model with continuous variables; every other method
refuses one with `CRN-R020` at `/method`. The request gives its tolerances in the method object:

```json
{"kind": "hybrid", "rtol": 1e-6, "atol": 1e-9}
```

with `rtol` in [1e-12, 1e-2] and `atol` in (0, 1e-2] (`CRN-R020` at
`/method/rtol` or `/method/atol`); `hybrid` needs rules-run/0.5 (`CRN-R020`
at `/schema`). The Python helper writes these values for `method="hybrid"`.

The model is a piecewise-deterministic Markov process (Davis 1984, DOI
[10.1111/j.2517-6161.1984.tb01308.x](https://doi.org/10.1111/j.2517-6161.1984.tb01308.x)),
simulated by integrated propensities (Alfonsi et al. 2005, DOI
[10.1051/proc:2005001](https://doi.org/10.1051/proc:2005001)):

1. After each event draw a unit exponential threshold E from one open
   uniform, E = −ln u.
2. Integrate the variables and the integrated total propensity
   Λ(t) = ∫ a₀ ds, with Λ = 0 at the event, by the Dormand–Prince 5(4)
   pair (DOI [10.1016/0771-050X(80)90013-3](https://doi.org/10.1016/0771-050X(80)90013-3)).
   Steps follow the standard error control of Hairer, Nørsett and Wanner
   (*Solving Ordinary Differential Equations I*, 2nd ed., 1993, II.4): the
   root-mean-square of the error estimate over all components, each scaled
   by `atol + rtol·max(|y₀|, |y₁|)`, must not exceed 1; the next step is
   h·min(10, max(0.2, 0.9·err^(−1/5))). The first step uses their starting
   rule, later ones the previous prediction, and each step is at most four
   times (E − Λ)/a₀ at its start (but not below 16 rounding units of the
   time).
3. The event happens where Λ reaches E, located by the Illinois method on
   the fourth-order dense output (II.6, after Shampine 1986, DOI
   [10.1090/S0025-5718-1986-0815836-3](https://doi.org/10.1090/S0025-5718-1986-0815836-3)).
   It selects channel j with probability a_j/a₀ evaluated at that point,
   from one uniform, with the resolution rule of Static Direct
   (`CRN-R031`).
4. Report values and bound crossings come from the same dense output; a
   bound crossing that comes first ends the step there.

Rates that do not read a variable are computed once per interval.

Independent variables (rules-model/0.4) are integrated once per prepared
run, by the same Dormand–Prince pair, step control and bound handling,
over [0, end] at the request's `rtol` and `atol`, with their own error
norm. The accepted steps and their fourth-order dense output are kept, and
every replicate and worker reads the variables from them, so results stay
independent of `threads`. Per replicate, the method then integrates only
the dependent variables and Λ, reading the independent variables from the
dense output at each stage time; the error norm covers those components
only. A stop bound of an independent variable ends the shared solution, and
every replicate that reaches its time stops there (a tie with the
threshold goes to the stop, as above). If the shared solve fails (a
non-finite derivative, a step that underflows, or more attempted steps than
`max_sampling_steps_per_replicate`), it ends at the start of the failing
step, and every replicate that reaches that time fails with the solve's
diagnostic (`CRN-R030`, or `CRN-R041` for the step budget). The solve runs
in `prepare` (and again for every `bind_parameters` candidate); the step
budget bounds it.

The method lowers the model like the static methods, so it refuses models with
dynamic instances; it refuses ATP models (use `thinning`); and it also runs
models without continuous variables, as an exact jump process up to
rounding. It is approximate at its tolerances. Its sampler is
`chacha8-u53/rules-approximate-hybrid-dopri54/1`; rules-result/0.3 records
the request's `method`, tolerances included. Integration steps, accepted or
rejected, count against `max_sampling_steps_per_replicate` like leap
iterations (`CRN-R041`), and a step below four rounding units of the time
fails the replicate (`CRN-R030`). Validation:
[`validation/hybrid/RESULTS.md`](../../validation/hybrid/RESULTS.md).

## Compiled event loops

`crn.compile_loop` builds the Static Direct event loop of one model and
parameter set into a separate program that `compiled=` then runs. The loop is
opt-in, and it needs a Rust toolchain. It gives byte-identical results to the
interpreter. The guide teaches when to use it in
[chapter 5](../guide/05-speed.md#compiled-models); the Rust side is in the
[Rust API](../developer/rust-api.md#generated-event-loops).

**Where `compiled=` goes.** A loop can be given in two places, and both run
the same program with the same results.

- `crn.ensemble(..., compiled=loop)`, or `crn.prepare(model, request,
  compiled=loop)`, attaches the loop to the prepared run. Preparation checks
  once that the loop was compiled for these model bytes and that the request
  selects `static_direct`, and raises `crn.CrnError` `CRN-R020` (at `/model`
  or `/method`) otherwise, before anything runs. Every `run.simulate()` then
  uses the loop.
- `run.simulate(compiled=loop)` uses the loop for that call only, on a run
  prepared without one; the first such call makes the same checks. The same
  prepared run can then also run interpreted.

A run with an attached loop always runs compiled: `simulate(callback=...)`
raises `ValueError`, so prepare the request without `compiled` to watch
events. `run.bind_parameters(...)` returns a candidate without the loop,
which runs interpreted, since a loop serves only the parameter set it was
compiled for.

**Requirements.** Rust 1.85 or newer, the MSRV of the crates (install it with
[rustup](https://rustup.rs)), and network access to crates.io: the program
depends on `crn-craft` at exactly the installed version. Offline, or with a
version that is not on crates.io (a wheel built from a clone), pass `source=`
a CRN-craft checkout of that version, or set `CRN_CRAFT_SOURCE`. `cargo=`
names another cargo executable.

**Build time and cache.** The first build compiles `crn-craft` and its
dependencies: about 40 s on the two-CPU host after the downloads, several
minutes on slower machines. `compile_loop` prints cargo's progress and the
path of the build log to stderr (`progress=False` silences it). Builds are
cached in the user cache directory (`~/.cache/crn-craft/rules-loops` on
Linux; `cache_dir=` or `CRN_CRAFT_CACHE_DIR` choose another), keyed by the
loop's fingerprint, the installed version, `rustc -vV`, the target and the
Rust flags. A cached call took 0.06 s. All builds share one Cargo target
directory (`CARGO_TARGET_DIR` if set), so a second model built in 3 s. The
cache takes about 160 MB; delete the directory to reclaim it. Loops cached by
release 0.6.0 are built again once on first use, with the same results, and
their entries stay in the cache until you delete it.

**Limitations.**

- Only `method="static_direct"`, with any options, and only the linear scan.
  `crn.ensemble` leaves the selector to `auto`, so for a model whose `auto`
  choice is another selector, name `selector="linear"` together with
  `compiled=` (`CRN-R020` at `/method/selector` otherwise). Other methods are
  refused with `CRN-R020` at `/method`.
- The generator emits every model that `method="static_direct"` runs:
  constant, vacancy and expression rates, string rates (`propensity` rules)
  and `sum`, `distinct`, `max` and `min` observables. It refuses, with
  `CRN-R020` at `/method` and before cargo runs, models that only another
  method runs (continuous variables, ATP blocks, instances that change
  identity), and loops that would exceed 1 MiB of source, which densely
  coupled networks with many channels reach. `crn.loop_source(model)` returns
  the source, or that error, without cargo.
- The loop is compiled for the exact model bytes. `crn.compile_loop` builds
  for the document `crn.ensemble` runs (`model/1.0`), so the same dictionary
  or `crn.Model` works for both calls. `crn.rules.compile_loop` builds for the
  bytes as given, which `crn.rules.ensemble` runs; a loop compiled for other
  bytes is refused with `CRN-R020` at `/model`.
- Parameters are compiled in, so another parameter set needs another loop. A
  loop passed for other model bytes, or for a `bind_parameters` candidate, is
  refused with `CRN-R020` at `/model`; candidates run interpreted.
- `callback` and `batch_callback` are refused with `ValueError`: the loop runs
  in its own process. Trackers and event logs work.
- Threads, report times, stops, budgets, `wall_time_s` and `cancel_file` work
  as interpreted. Ctrl-C stops the program and raises `KeyboardInterrupt`.
- On networks with many constant-rate channels the loop can be slower than the
  interpreter
  ([details](../developer/specialized-kernels.md#generated-rules-event-loops)).
- Missing cargo, a failed build (`error.log` is the build log) and a result
  that differs raise `crn.CompiledLoopError`, whose `reason` is
  `"cargo_missing"`, `"build_failed"` or `"mismatch"`.

**Confirming identity.** Every `compile_loop` call, cached or not, runs short
requests through the program and the interpreter, with both waiting times,
both tracker policies and an event log, and returns the loop only if the
metadata and every record column are byte-identical. `loop.verify(run)` does
the same for your own request and returns the result. On a mismatch both
results are saved next to the build for inspection.

**Speed.** Measured from Python including the subprocess, on 12 November
2026, a compiled loop ran the three-class model about 2.1 times as fast as
interpreted Static Direct with its defaults. On a susceptible-infected-removed
model with 2,000 replicates (37.6 million events) the loop took 1.64 s
against 3.35 s with inverse waits (2.00 times) and 1.36 s against 2.97 s with
Ziggurat waits (2.18 times), medians of three on the two-CPU host
(`validation/rules/compiled-loop-string-rates.json`). Release 0.7.0 made
interpreted Static Direct faster (0.72 times the time on the [three-class
workload](../evidence/mtdna-performance-analysis.md#release-070)). Cached
rates and dependency-aware projection refresh apply to the default path;
potentially overflowing report-only projections keep their eager checks. The
[rules kernel review checklist](../developer/specialized-kernels.md) has the
conditional Rust generator and its adoption gate.
