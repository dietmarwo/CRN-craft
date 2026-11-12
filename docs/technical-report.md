# CRN-craft technical report

**Version:** 0.7.0, source of 12 November 2026: the release that merged
the network engine and the rules engine into one engine.
**Scope:** design, numerical methods, validation evidence and measured
performance of the one engine. Companion documents: the
[architecture overview](architecture.md), the [contracts](contracts.md), the
[exact methods](direct-ssa.md), the [model reference](rules-language.md) and
the [merge implementation log](campaign/merge-0.7/README.md).
**Evidence:** every gate below was declared in
[`validation/merge/RESULTS.md`](../validation/merge/RESULTS.md) before it
ran, and its report is committed under `validation/merge/`. Failed rows are
kept there.
**History:** the report of release 0.6.0 described the network engine of
releases 0.2.0 to 0.6.0. Its numbers are summarized, and marked, in
section 9.

## 1. Summary

CRN-craft is a Rust library and command-line tool, with a Python package,
for stochastic simulation and analysis of well-mixed chemical reaction
networks, including networks whose species carry attributes and appear
during a run. One engine provides exact simulation (Direct with four event
selectors, next reaction, a dynamic tree method, thinning and a certified
specialization), approximate tau leaping, binomial and hybrid methods, a
finite state projection solver for the chemical master equation, exact
conservation and graph analysis, ODE right-hand sides with Jacobians, and
strict importers and exporters.

- **Reproducible to the bit.** A result depends only on the exact input
  bytes, the sampler identity and the seed. Replicate *i* draws from its own
  keyed ChaCha8 stream, sums are formed in a fixed order with explicit
  resolution checks, and ensembles reduce in replicate order, so one, two or
  eight threads give identical results. Every stored request of the former
  rules engine (`rules-run` 0.1 to 0.6) replays byte for byte, checked on
  every build.
- **Honest outputs.** Incomplete ensembles withhold statistics; FSP reports
  leakage, series truncation and roundoff separately and never renormalizes;
  approximate methods carry their own sampler identities.
- **One engine, the same laws.** Results of the former network engine do not
  replay: its loops and its xoshiro256++ stream are gone. Its models convert,
  and the unified engine samples the same laws: the law gate of 4,814
  comparisons over twelve models passed with 0 failures for every selector
  and for the default (section 6.2). The features the network engine had
  alone (FSP, analysis, kinetics, export) were ported to a view of the
  lowered model and checked against references recorded from release 0.6.0
  (section 6.3).
- **Speed.** The declared speed gate allowed the best unified method 1.25
  times the time of the best of 21 methods of release 0.6.0 on each of
  twelve models. Its rerun passed, with ratios from 0.66 to 1.20; the
  margins of five models lie within this host's run-to-run noise
  (section 7.1). On the mtDNA workload, 0.7.0 takes 0.72 to 1.00 times the
  time of 0.6.0, with identical results (section 7.2).

## 2. Design goals

1. **Exactness is a contract.** A method called exact samples the stochastic
   process up to floating-point arithmetic whose order is fixed and whose
   failure modes (lost increments, unresolvable rates, stalled clocks) are
   detected, not hidden. Two implementations that claim one sampler
   identity must agree bit for bit.
2. **Traceability.** Every result names the SHA-256 of the exact model and
   request bytes and its sampler identity; archives seal a run with a
   manifest and a completion marker.
3. **Explicit limits.** Every algorithm has published size, work, memory and
   output limits. A run is admitted only if its memory reservation fits.
   Exceeding a limit produces a diagnostic, never a truncated result
   presented as complete.
4. **Small trusted core.** The engine has no Python, job or ODE dependency,
   forbids unsafe code, pins every dependency, and uses no GPL code.
   Optional integrations live in separate workspaces.

## 3. Architecture

The [architecture overview](architecture.md) describes the layers. In brief:
versioned JSON documents are parsed strictly into Rust types that also
generate the JSON Schemas; the model is compiled (patterns into classes,
predicates into bytecode); a request is validated and admitted against its
memory reservation; a closed model is lowered to a flat network of species
and channels whose rates are compiled into programs with folded parameters;
and one of the methods runs the ensemble. The flat view of the lowered model serves
FSP, analysis, kinetics and export. The command-line tool and the PyO3
binding are thin layers over the same crate, `crn-craft`. After the merge,
`crates/*/src` holds 29,372 lines of Rust: 28,492 in `crn-craft` (48
files), 788 in the command-line crate and 92 in the facade, whose crate
documentation the final review shortened. Before the network engine was
deleted it held 37,157 (M8b log of the
[implementation log](campaign/merge-0.7/README.md)).

## 4. Numerical methods

### 4.1 Models, lowering and propensities

A model declares types with numeric attributes, rules whose reactant
patterns match instances, observables (weighted sums, counts of distinct
instances, maxima and minima) and an optional stop expression. A rule's
`rate` is per combination of reactant molecules: a rule with patterns of
counts *k₁, …, k_m* over species with counts *n₁, …, n_m* has propensity
*rate · factor · ∏ C(nᵢ, kᵢ)*, the combinatorial mass action of the
stochastic convention (`A + A` gives *c·n(n − 1)/2*). A whole-rule
`propensity` is used as written while every pattern can take its molecules.
Higher-order products (a pattern `count` above one, or more than two
patterns) are evaluated as a scaled mantissa–exponent product that does not
overflow before the result does (`crn_craft::rules::combination_rate`).

A closed model is lowered once, at preparation: species (attribute
vectors), channels (one per choice of matching instances), the count
changes of each channel, the species each rate reads, and for each channel
the sorted list of channels whose rates its firing can change. Lowering is
index-based and accepts 4,096 instance tuples and 65,536 channels in
`model/1.0`. Rates are compiled into flat programs; the converted Hill and
Michaelis–Menten forms have dedicated programs. A randomized test compares
9,000 compiled programs, and the two law shapes at 2,000 points each, with
the expression tree evaluator: the values are equal.

### 4.2 Static Direct

For state *x* with channel propensities *a_j(x)* and total *a₀*, the Direct
method (Gillespie 1977) draws a waiting time and a channel with probability
*a_j/a₀*.

- **Random words.** Replicate *i* uses ChaCha8 (`rand_chacha`) keyed by
  SHA-256 of the domain `crn-craft.rules/0.1` and the root seed as eight
  little-endian bytes, with stream number *i*. A uniform is the high 53 bits
  of a 64-bit word times 2⁻⁵³.
- **Waiting time.** The inverse transform *τ = −ln(1 − u)/a₀*, with *u*
  redrawn when it is 0, or the Ziggurat exponential (tables adapted from the
  Julia ecosystem, see [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md)).
  Ziggurat waits are the default from `rules-run/0.6` and in `run/1.0`.
- **Selection with the linear scan.** The total is the sequential sum of the
  rates in channel order. After an event only the affected channels change,
  so the sum resumes at the first of them with the same additions. The
  channel is found by subtracting rates from *u·a₀* in channel order.
- **Resolution and failures.** A positive rate below 2⁻⁵³ of the total
  cannot be selected by a 53-bit uniform and fails the replicate
  (`CRN-R031`); a lower bound on the smallest positive rate avoids the exact
  check on most events. A negative or non-finite rate fails it
  (`CRN-R030`). A zero total is absorbing: the state is held to the horizon.
- **Time.** Reports are filled before an event that would pass them, with
  the pre-event state; an event exactly at a report time is included. Under
  `run/1.0` a positive wait too small to move the clock is accepted, and
  1,000,000 consecutive such steps fail with `CRN-R030`; a step then leaves
  the clock unchanged with probability about *1 − exp(−a₀·ulp(t)/2)*, so the
  limit is reached only when about 22 or more events are expected per
  representable step of the clock. Every `rules-run` version fails at the
  first stalled step.
- **Updates.** Count updates are checked and applied in one pass that undoes
  itself on a range or instance error. Observables are refreshed along
  dependencies; report-only observables that cannot overflow are computed
  at report times.

### 4.3 Selectors

Under `run/1.0` Static Direct takes a `selector`. Every selector samples the
same law; each consumes the random stream differently and has its own
sampler identity ([exact methods](direct-ssa.md)).

- **Composition–rejection** (Slepoy, Thompson and Plimpton 2008): channels
  grouped by the binary exponent of their rate, a group drawn by its mass
  from pairwise sums with the strictness of the linear scan, a member by
  rejection with Lemire's joint member-and-acceptance draw.
- **Sorting** (McCollum et al. 2006): the linear scan in an order that
  moves each fired channel one place forward; a threshold at or above the
  sequential sum is redrawn.
- **RSSA-CR** (Thanh, Priami and Zunino 2014): each species has an interval
  around its count; each channel has lower and upper bounds of its rate over
  the box. A rule whose rate reads one species is evaluated at every count of
  that species' interval with the engine's arithmetic, so its range is exact
  for any expression; a rule that reads two or more species pins them, and
  its bounds are its exact value. A candidate drawn by composition–rejection
  from the upper bounds is accepted below its lower bound, otherwise against
  its exact rate, reusing the composition's accepting threshold. So RSSA-CR
  runs every model that Static Direct runs.
- **`auto`**, the default of `run/1.0`: RSSA-CR from 8 channels when no rate
  reads two or more species (from 128 channels when the species that start
  with molecules average fewer than 8), composition–rejection from 128
  channels, otherwise the linear scan. The rule was fitted by declared
  criteria to a benchmark of eight synthetic families and the twelve gate
  models; every pick on the gate models is the fastest selector except
  michaelis-menten (linear, 1.08 times sorting).

### 4.4 Other exact methods

- **Next reaction** (Gibson and Bruck 2000) keeps one clock per channel in an
  indexed heap and rescales the clocks of the affected channels.
- **Tree** recomputes rule propensities over sum trees of instances. It runs
  models whose instances change identity or draw random product attributes,
  and refuses higher-order rules.
- **`mtdna_direct`** is a grouped Direct method for the fixed three-class
  mtDNA structure. Preparation certifies the model and each parameter set:
  in every reachable state each positive hazard is at least 16 × 2⁻⁵³ of the
  total. Otherwise Static Direct runs with a `CRN-R051` warning.
- **Thinning** integrates the analytic ATP flow of models with the ATP block
  between events and thins against a checked intensity envelope.

### 4.5 Approximate methods

- **Tau leaping** follows the step selection of Cao, Gillespie and Petzold
  (2006): the species bound uses *g = 2* for expression rates, and a
  propensity condition bounds the mean and spread of every rate change by
  *ε·a_j* through first differences, so rates that read global observables
  are controlled directly. Channels that ten firings could switch off are
  critical and fire at most once per leap; unacceptable proposals halve the
  step before an exact fallback. Poisson variates with mean at most 16 use
  the multiplication method, larger means Hörmann's PTRS (1993). The leap
  condition costs *O(channels²)* per attempted leap (section 6.6).
- **Binomial splitting** resolves fixed steps by competing-risk
  probabilities *(a/Λ)(1 − e^(−Λ·dt/n))* per parent for single-parent unit
  births and deaths.
- **Hybrid** integrates continuous variables with Dormand–Prince 5(4) at the
  request's tolerances between stochastic events.

Their results are approximate at their declared parameters, and their
sampler identities say so.

### 4.6 Ensembles, statistics and quantiles

Workers claim replicate ids from a shared counter; a window of four
replicates per worker bounds how far they run ahead, and results reduce in
ascending id order, so statistics do not depend on the thread count.
Quantiles keep one value per replicate, report time and observable, and
use NumPy's default definition (Hyndman and Fan type 7) evaluated as NumPy
evaluates it, so the values equal `numpy.quantile` bit for bit. Admission
counts that buffer against the 2 GiB workspace limit. Statistics,
quantiles included, are withheld unless every replicate completed or
stopped.

### 4.7 Finite state projection

Following Munsky and Khammash (2006), an `fsp` request enumerates the
reachable states of the flat view within per-species bounds by breadth-first
search (up to `max_states`), with every transition leaving the box sent to
an absorbing sink. A state where the model's stop holds is absorbing as
well. The probability vector is propagated by uniformization, each interval
split so that *ΛΔ ≤ 16*, and each series truncated when its geometric tail
majorant falls below the step's share of the tolerance. Leakage, the series
tail bounds, the omitted series weight, the total mass defect and the
roundoff residual are reported separately per snapshot, and probabilities
are never renormalized. The generator is stored by column and each series
term gathers in ascending state order.

### 4.8 Analysis, kinetics and export

- **Conservation** computes an exact rational basis of the left nullspace of
  the stoichiometric matrix with checked i128 arithmetic; the **cone**
  enumerates the extreme rays of the nonnegative conservation cone by a
  bounded double description method (after Burton) and bounds every species
  of the initial class; **linkage** computes complexes, linkage classes,
  deficiency and weak reversibility. Each stops at explicit limits.
- **Kinetics** gives the right-hand side *rate · factor · ∏ yᵢ^kᵢ/kᵢ!* (or a
  whole-rule propensity without its availability gate, accepted only where
  the expression vanishes structurally when a taken species is zero) and its
  state and parameter Jacobians, from dual-number evaluation of the rule
  expressions. At most 64 species, 32 parameters and 256 channels.
- **Export** writes GillesPy2 source: mass action for constant rates of at
  most two molecules, custom propensity strings otherwise, and refusals
  (`CRN-E050`) for what GillesPy2 cannot evaluate in both of its solvers.

### 4.9 Conversion of network models

The converter writes each species as an attribute-free type, mass action
as a `rate`, Hill laws as `basal + vmax / (1 + (k / R)^n)` and
`basal + vmax / (1 + (R / k)^n)`, and Michaelis–Menten laws as
`vmax * (S / (km + S))`. Against the network engine's own evaluation
(`crates/crn-craft/src/rules/convert.rs` tests): mass action of order 0
and 1 and Michaelis–Menten are bit-identical, pairs within 4 ε, Hill laws
within 1e-14 · max(1, n) relative, and 116 of 5,880 Hill cases where the
induced term is below *vmax · 2⁻¹⁰²⁰* give the basal rate.

## 5. Interfaces and contracts

The command line is `crn-craft describe`, `check`, `run` (with `--archive`
and `--cancel-file`), `verify`, `analyze`, `convert` and `export`; all
output is JSON with documented exit codes. The contracts are
`crn-craft.model/1.0`, `run/1.0`, `result/1.0`, `fsp-result/1.0`,
`check/1.0`, `analysis/1.0` and `archive/1.0`, with the `rules-*` versions
still read ([contracts](contracts.md)). Contract types generate JSON
Schemas (`crn-craft describe KIND`), and any shape change gets a new
version.

The Python package (ABI3 wheel for Python 3.10 and later) exposes the
builder, prepared ensembles, pandas and xarray results, `Result.save`, FSP
and analysis on models, importers and compiled event loops. Native calls run
off the GIL, so Ctrl-C cancels them within about 50 ms; rejected inputs raise
`CrnError` with parsed diagnostics. `crn_craft.rules` keeps the 0.6.0 API
with its defaults.

## 6. Validation

### 6.1 Test suites and invariants

| Gate | Content |
|---|---|
| Rust tests | `cargo test --release --workspace --locked`: 397 tests (29 of them doctests) on the merged tree (final review, 12 November 2026; 538 before M8b removed the network engine's tests); contracts, conversion, limits, higher-order rules, selectors, quantiles, the clock rule, FSP, analysis, kinetics, export, archives, generated code and the CLI |
| I1: stored requests | `crates/crn-craft/tests/stored_requests.rs` and `static_direct_defaults.rs`: every stored `rules-run` request replays byte for byte; the same requests as `run/1.0` give the `rules-run/0.6` result apart from schema and hashes |
| I2, I2b: manuals | `check_manual.py` runs every block of the user manual (26) and of the frozen 0.6.0 manual (22), normally and under `python -O`, with identical reports |
| I3: generated code | 20 checked-in generated sources (19 event loops, 1 rates kernel) equal the generators' output and replay the interpreter byte for byte |
| Mutation gate | `mutation_check.py`: 47 source mutants, each turning a passing suite into named failures; 47 of 47 caught on the unified code (M8a and M8b logs, and again in the final review) |
| Exact laws | `check_ssa.py` (analytic campaigns, an independent Python oracle that reproduces the engine's paths, thread identity, cancellation latency), `check_cme.py` (exact master-equation checks), `check_alternative_ssa.py` (168 DKW tests: 4 selectors × 2 waiting times, N = 100,000, family α = 0.001) |
| Features | `check_fsp.py` (including the DSMTS suite), `check_tau.py`, `check_network.py`, `check_structure.py`, `check_export.py`, `check_importers.py`, `check_contracts.py` |
| Python and documents | `check_python.py`, `check_rules_python.py`, `check_builder.py`, `check_hybrid.py`, `check_installed.py`, `check_tutorial.py`, `check_compiled_loop.py`, `check_links.py` |
| Dependency gate | the default build's dependency set and licenses |

Every Python checker runs normally and under `python -O`, and CI diffs the
reports.

### 6.2 The law gate (N = 4,814)

Declared before the first run ("Law gate" in `RESULTS.md`): twelve models
(the five examples and the workloads W1 to W7, up to 256 species, 992
reactions and order 8), each simulated with 200,000 replicates in chunks
with declared seeds, and compared per species and report time with
references recorded from release 0.6.0 (200,000 replicates of its cached
Direct method) and, for five models, with FSP laws. A comparison fails when
a mean or a variance differs by more than five standard errors of the
difference. Of the 4,814 comparisons, 940 are at time 0 and must be equal;
about 0.002 failures are expected per run of a correct engine, so 0 are
expected.

| Engine configuration | Report | Comparisons | Failures | Largest \|z\| |
|---|---|---:|---:|---:|
| `linear` | `laws/m4c-linear.json` | 4,814 | 0 | 3.88 |
| `composition_rejection` | `laws/m4c-composition_rejection.json` | 4,814 | 0 | 3.83 |
| `sorting` | `laws/m4c-sorting.json` | 4,814 | 0 | 3.51 |
| `rssa_cr`, first run | `laws/m4c-rssa_cr.json`, `laws/m4c-rssa_cr-vilar-noise-failure.json` | 4,724 (11 models) | **vilar-noise failed** | 3.94 |
| `rssa_cr`, after the clock rule | `laws/m4f-rssa_cr.json` | 4,814 | 0 | 3.88 |
| `auto`, the `run/1.0` default | `laws/m4d-auto.json` | 4,814 | 0 | 3.88 |

Report paths are under `validation/merge/`. The first RSSA-CR run failed on
vilar-noise: a replicate stopped at t = 7.16 because a trial's wait was
below half an ulp of t (`CRN-R030`, clock failed to advance), as the 0.6.0
engine also did. The clock rule of section 4.2 accepts such steps under
`run/1.0`; the rerun with the declared seeds passed. The `auto` gate gave
event counts and power sums equal to the gates of the selectors it picked.
Earlier runs on subsets passed as well (`laws/static_direct.json`: 1,102
comparisons before the higher limits; `laws/m4-static_direct.json`: all
twelve models after the higher limits and higher-order rules, 4,814, 0
failures).

### 6.3 References of the ported features

| Feature | Report (`validation/merge/`) | Result |
|---|---|---|
| FSP against the 0.6.0 solver | `fsp/references.json` | 11,278 checks, 0 failures; largest probability difference 1.4e-17, six of seven models bit-identical |
| Analysis of the twelve models against 0.6.0 | `analysis/references.json` | 194,388 leaves, 0 failures |
| Kinetics against 0.6.0 | `kinetics/references.json` | 6,252 values and 12 acceptance checks, 0 failures; 6,192 bit-identical, largest relative difference 2.6e-16 |
| GillesPy2 export against 0.6.0 | `export/references.json` | 9 artifacts byte-identical, 160 custom-propensity evaluations equal, 0 failures |
| Higher-order rules against FSP laws | `laws/higher-order.json` | 272 comparisons, 0 failures |
| Tau envelopes | `tau/check_tau_rules.json` | 3 rows (ε = 0.1, 0.05, 0.025; 10,000 paths each), 0 failures |
| Quantiles against the Poisson law | `laws/quantiles-immigration-death.json` | 9 comparisons, 0 failures |
| Selector laws | `laws/selectors-exact-m4f.json` | 168 of 168 DKW tests, 16 of 16 power checks |

### 6.4 Limitations

- Validation uses analytic and master-equation references for small
  networks and synthetic circuits; it does not establish biological
  validity of any model.
- The law gate compares marginal means and variances, not covariances or
  whole distributions; FSP references exist for five of the twelve models.
- Tau leaping is approximate, its envelopes cover immigration–death
  controls only, and on the gate model most of its steps are exact
  fallbacks.
- Timings come from one shared two-CPU host; absolute numbers vary by
  machine.
- The optional integrations (`crn-diffsol`, the job worker) are built and
  checked locally, not in CI. `check_ode.py` is not a CI step.
- macOS and Windows wheels are built by a manual workflow and are not
  executed here.

## 7. Performance

### 7.1 Speed gate

Declared before any measurement ("Speed gate (M4d)" in `RESULTS.md`):
twelve models (W1 to W7 at the workload sizes, the five examples), the wall
time of one CLI process, median of three after a warm-up, one thread. The
0.6.0 best is the smallest median over 21 methods of release 0.6.0 (five
aggregators with the `fast` policy on and off and both waiting times, and
next reaction), rerun in the same session with the kept 0.6.0 build. The
unified best is the smallest median over `linear`, `composition_rejection`,
`sorting`, `rssa_cr`, `auto` and `next_reaction`. Gate: unified best at most
1.25 times the 0.6.0 best on each model.

The first run failed on W5 (1.28; `validation/merge/speed-m4d-run1.json`).
After the one allowed fix round, the rerun passed
(`validation/merge/speed-m4d.json`):

| Model | 0.6.0 best | s | Unified best | s | Ratio |
|---|---|---:|---|---:|---:|
| W1 | `direct/ziggurat/fast` | 0.0482 | `auto` (linear) | 0.0511 | 1.06 |
| W2 | `direct_cr/ziggurat/fast` | 0.0906 | `rssa_cr` | 0.1031 | 1.14 |
| W3 | `cached_direct/ziggurat/fast` | 0.2248 | `auto` (linear) | 0.2707 | 1.20 |
| W4 | `rssa_cr/ziggurat/fast` | 0.1119 | `auto` (rssa_cr) | 0.1084 | 0.97 |
| W5 | `direct_cr/ziggurat/fast` | 0.0614 | `auto` (rssa_cr) | 0.0685 | 1.12 |
| W6 | `rssa_cr/inverse/fast` | 0.0639 | `auto` (rssa_cr) | 0.0599 | 0.94 |
| W7 | `cached_direct/ziggurat/fast` | 0.0185 | `rssa_cr` | 0.0216 | 1.17 |
| enzyme-mechanism | `cached_direct/ziggurat/fast` | 0.0731 | `sorting` | 0.0622 | 0.85 |
| immigration-death | `cached_direct/ziggurat/fast` | 0.0687 | `auto` (linear) | 0.0474 | 0.69 |
| michaelis-menten | `direct/ziggurat/fast` | 0.0561 | `auto` (linear) | 0.0368 | 0.66 |
| repressilator | `cached_direct/ziggurat/fast` | 0.2809 | `linear` | 0.3178 | 1.13 |
| vilar-noise | `cached_direct/ziggurat/fast` | 1.0514 | `auto` (linear) | 0.9451 | 0.90 |

**Noise.** `auto` runs the same computation as the selector it picks, so
each model has two timings of one computation; in the rerun they differ by
up to 1.55 times (W6: 0.0599 against 0.0930 s). With the slower of the two,
the ratios of W1 to W7 would be 1.21 to 1.46. Pooled over both runs (each
side's smaller median), the ratios are W1 1.05, W2 1.15, W3 1.18, W4 1.04,
W5 1.14, W6 0.94, W7 1.05, repressilator 1.18 and 0.66 to 0.90 for the other
examples. The 0.6.0 side, a minimum over 21 close variants, profits from the
same effect. So the gate passes by its declared rule, but for W2, W3, W5, W7
and the repressilator (1.12 to 1.20) the margin lies within this host's
run-to-run noise, and W5's move from 1.28 to 1.12 is mostly noise.

Where the time goes (callgrind of the rerun's CLI, M4d log): on W1 the two
ChaCha8 words per event (about 3.7 ns per word against 1.5 ns for 0.6.0's
xoshiro256++) and the Ziggurat wait; on W3 and the repressilator one
mispredicted scan exit per event and the Hill programs; on the RSSA-CR
models the trials per event and the exact evaluation of candidates. The
unified engine runs fewer instructions than 0.6.0's best on every profiled
model (0.63 to 0.99 times, M4f log), but not always less wall time.

Compiled event loops, timed at the Python boundary (informative, same
section of `RESULTS.md`), took 0.61 to 0.90 times the interpreted linear
scan on eight models and 1.15 times on W7; W4 to W6 exceed the 1 MiB
source limit.

### 7.2 mtDNA workload

Measured on 12 November 2026 with **no declared target**
(`validation/merge/mtdna-workload/`, results section in `RESULTS.md`): the
three-class model `models/rules/mtdna-three-class-cme.rules.json`, a
37-point grid from 0 to 25,920 h, seed 1, wall time of
`crn_craft.rules.ensemble(...).simulate()`, medians of five runs alternating
0.6.0 and 0.7.0 under a CPU lock. Statistics and records hash identically
between the versions.

| Method | Cells | Threads | 0.6.0 s | 0.7.0 s | Ratio |
|---|---:|---:|---:|---:|---:|
| `mtdna_direct` | 1000 | 1 | 10.67 | 10.27 | 0.96 |
| `mtdna_direct` | 1000 | 2 | 5.42 | 5.26 | 0.97 |
| `static_direct` | 400 | 1 | 13.20 | 9.49 | 0.72 |
| `static_direct` | 400 | 2 | 6.70 | 5.63 | 0.84 |
| `tree` | 100 | 1 | 10.04 | 9.94 | 0.99 |
| `tree` | 100 | 2 | 5.46 | 5.45 | 1.00 |

Static Direct's linear scan runs a tuned loop in 0.7.0 (compiled rate
programs, running totals, precomputed count factors) with byte-identical
results. The [performance analysis](mtdna-performance-analysis.md) sizes
mtDNA studies.

### 7.3 Other measurements

- **Parameter workers.** Binding parameters on a prepared run, against the
  0.6.0 templates in one session (`validation/merge/parameter-workers-m4d.json`,
  declared target 1.25): five of six ratios between 1.04 and 1.15;
  `binding12` with reuse **failed** at 1.55, because `bind_parameters`
  prepares the whole model again on each call (about 0.1 ms on that
  12-species model).
- **Tau leaping on large networks** (`validation/merge/tau/leap-condition-cost.json`):
  on W5 and W6 (510 and 992 channels) `tau` took 14.1 and 12.9 times as long
  as Static Direct; neither model leaps at the declared ε, so this is the
  cost of its exact fallback. Each attempted leap evaluates every channel
  after one firing of each enabled channel. Use Static Direct there.
- **GillesPy2 1.8.3** ([comparison](comparison-gillespy2.md),
  `validation/merge/gillespy2-comparison-0.7.json`, protocol declared in
  `RESULTS.md`): on the four frozen workloads (64 paths, one thread) its
  compiled C++ solver took 2.9, 10.6, 3.6 and 9.3 times as long as
  CRN-craft's default request (`auto`, Ziggurat waits), and 2.8, 5.6, 3.7
  and 1.3 times as long as the linear scan with inverse waits; it needs 2.1
  to 2.9 s to compile each model. Its NumPy solver took 134 to 358 times as
  long as the default. All 28 agreement tests passed the declared rule.
- **Julia (Catalyst and JumpProcesses).** Not repeated for 0.7.0: Julia is
  not available on this host. The [comparisons](comparison-catalyst-jumpprocesses.md)
  measured the 0.6.0 network engine and carry a note saying so.

## 8. Reproducing

```sh
cargo build --release --workspace --locked
cargo test --release --workspace --locked
python validation/merge/check_laws.py --help          # the law gate
python validation/merge/compare_speed.py --help       # the speed gate
python validation/compare_gillespy2.py --binary target/release/crn-craft \
    --output gillespy2.json --figure gillespy2.png
python validation/merge/mtdna-workload/bench.py \
    models/rules/mtdna-three-class-cme.rules.json static_direct 1 400 timings.jsonl
```

The speed gate needs the kept 0.6.0 build; `validation/merge/RESULTS.md`
names the commands, sizes and seeds of every gate.

## 9. History: the network engine of releases 0.2.0 to 0.6.0

The report of release 0.6.0 (`git show v0.6.0:docs/technical-report.md`)
describes the network engine and its performance campaign; the numbers below
come from it and from the campaign files it cites, and none of them was
measured on 0.7.0.

- The network engine had two Direct implementations, `simulate` and the
  dependency-cached `simulate-cached`, with byte-identical results, a
  separate stream domain (`crn-craft.chacha8/0.1`), FSP, tau leaping with
  its own request and result contracts, and sealed archives with
  provenance. The `ssa` command added DirectCR, SortingDirect, RSSA-CR and
  next reaction, and the `fast` policy with xoshiro256++.
- Its performance campaign P01–P20 ([P01–P10](campaign/performance.md),
  [P11–P20](campaign/performance-p11-p20.md)) brought `simulate-cached` to
  48, 152, 151 and 553 ns per event on the four frozen workloads (two
  channels, sparse 66 reactions, a 12-reaction circuit, dense 66 reactions),
  on the same kind of shared two-CPU host.
- Against GillesPy2 1.8.3 on those workloads, its compiled C++ solver took
  1.1 to 3.3 times as long as `simulate-cached`.

## References

- D. T. Gillespie (1977). Exact stochastic simulation of coupled chemical
  reactions. *J. Phys. Chem.* 81, 2340–2361.
- M. A. Gibson and J. Bruck (2000). Efficient exact stochastic simulation of
  chemical systems with many species and many channels. *J. Phys. Chem. A*
  104, 1876–1889.
- A. Slepoy, A. P. Thompson and S. J. Plimpton (2008). A constant-time
  kinetic Monte Carlo algorithm for simulation of large biochemical reaction
  networks. *J. Chem. Phys.* 128, 205101.
- J. M. McCollum, G. D. Peterson, C. D. Cox, M. L. Simpson and
  N. F. Samatova (2006). The sorting direct method for stochastic simulation
  of biochemical systems with varying reaction execution behavior.
  *Comput. Biol. Chem.* 30, 39–49.
- V. H. Thanh, C. Priami and R. Zunino (2014). Efficient rejection-based
  simulation of biochemical reactions with stochastic noise and delays.
  *J. Chem. Phys.* 141, 134116.
- B. Munsky and M. Khammash (2006). The finite state projection algorithm for
  the solution of the chemical master equation. *J. Chem. Phys.* 124, 044104.
- Y. Cao, D. T. Gillespie and L. R. Petzold (2006). Efficient step size
  selection for the tau-leaping simulation method. *J. Chem. Phys.* 124,
  044109.
- W. Hörmann (1993). The transformed rejection method for generating Poisson
  random variables. *Insurance: Mathematics and Economics* 12, 39–45.
- R. J. Hyndman and Y. Fan (1996). Sample quantiles in statistical packages.
  *The American Statistician* 50, 361–365.
- B. A. Burton (2010). Optimizing the double description method for normal
  surface enumeration. *Mathematics of Computation* 79, 453–484;
  arXiv:0808.4050.
- D. J. Bernstein (2008). ChaCha, a variant of Salsa20.
- A. Dvoretzky, J. Kiefer and J. Wolfowitz (1956), with Massart's constant
  (1990): the DKW band used by the statistical checks.
