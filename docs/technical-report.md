# CRN-craft technical report

**Version:** `0.2.0-dev.1`, source of 1 November 2026 (campaign M3–M10,
performance iterations P01–P20 and the review of this date).
**Scope:** design, numerical methods, validation evidence and measured
performance of the exact and approximate solvers. Companion documents: the
[architecture overview](architecture.md), the [tutorial](tutorial.md) and the
performance declarations [P01–P10](campaign/performance.md) and
[P11–P20](campaign/performance-p11-p20.md).
**Later additions** are documented separately: the `ssa` command with six
exact aggregators and the `fast` policy ([alternative SSA](alternative-ssa.md)),
discovery 0.10 (native CLI: 0.11 in 0.2.0, 0.12 in 0.3.0 and 0.13 in
0.4.0 with rules-run/0.6), and the experimental rules engine
([rules language](rules-language.md), [mtDNA ensembles](mtdna-rules.md)).
All of them first shipped in release 0.2.0 (4 November 2026). The comparisons cited in
section 7 were measured on development builds before that release.

## 1. Summary

CRN-craft is a Rust library and command-line tool, with an optional Python
package, for stochastic simulation and analysis of well-mixed chemical
reaction networks. It provides two exact implementations of Gillespie's Direct
method, a finite state projection solver for the chemical master equation,
experimental tau leaping, exact conservation and graph analysis, strict
importers and exporters, and integrations with the fcmaes optimization worker
and the diffsol ODE solver.

Main properties and results:

- **Reproducible to the bit.** Results depend only on the exact input bytes:
  replicate *i* draws from its own keyed ChaCha8 stream, sums are formed in a
  fixed order with explicit resolution checks, and ensembles are reduced in
  replicate order. One, two and eight threads give identical results; the
  plain and dependency-cached implementations give byte-identical results; and
  the final source reproduces the numerical digests of the original M2
  baseline (5a07334) on every benchmark workload.
- **Honest outputs.** Incomplete ensembles withhold statistics; FSP reports
  leakage, series truncation and roundoff separately and never renormalizes;
  approximate methods use their own schema and assessment.
- **Validated.** 166 Rust tests in the core and CLI and 9 in the
  integrations; a 45-mutant source mutation gate; analytic
  SSA campaigns of 300,000 paths and exact master-equation checks of
  4 × 100,000 paths with second-order, Hill and Michaelis–Menten channels;
  FSP checks including the DSMTS suite; tau envelopes over 30,000 paths;
  Python, importer, installed-package and tutorial checkers, each run normally
  and under `python -O` with identical reports.
- **Fast where it matters.** One core of a shared 2-vCPU host runs the exact
  SSA at 48 ns per event on a two-channel network and 151 ns per event on
  12- and sparse 66-reaction networks with `simulate-cached`. Against the
  original baseline, `simulate` is 1.4–2.0 times faster on the frozen
  workloads, and `simulate-cached` is a further 2.1–3.3 times faster on the
  12- and sparse 66-reaction networks. FSP is 1.4–2.3 times faster than at P10,
  and tau leaping with large leap means about 300 times faster.

## 2. Design goals

1. **Exactness is a contract.** A method called exact reproduces the
   stochastic process up to floating-point arithmetic whose order is fixed and
   whose failure modes (lost increments, unresolvable rates, stalled clocks)
   are detected, not hidden. Two implementations of the same method must agree
   bit for bit.
2. **Traceability.** Every result names the SHA-256 of the exact model and run
   bytes; archives seal results with a manifest and a completion marker; the
   provenance records toolchain, lockfile and executable hashes.
3. **Explicit limits.** Every algorithm has published size, work, memory and
   output limits. A run is admitted only if its memory reservation fits.
   Exceeding a limit produces a diagnostic, never a truncated result presented
   as complete.
4. **Small trusted core.** The engine has no Python, job or ODE dependency,
   forbids unsafe code, pins every dependency, and uses no GPL code.
   Optional integrations live in separate workspaces.

## 3. Architecture

The [architecture overview](architecture.md) describes the layers in detail.
In brief: versioned JSON documents (generated JSON Schemas) are parsed
strictly, validated into diagnostics with stable codes, compiled into an
immutable indexed network with sparse stoichiometry and a reaction-dependency
graph, and executed by one of four engines (Direct, cached Direct, FSP, tau).
SSA results can be sealed into verifiable archives. The command-line tool and
the PyO3 binding are thin layers over the same core.

| Component | Lines (Rust/Python) | Role |
|---|---|---|
| `simulation.rs` | 900 | Direct and cached Direct, ensembles, budgets, polling |
| `tau.rs` | 640 | tau leaping, CGP step size, PTRS sampling |
| `fsp.rs` | 480 | reachable-state projection, uniformization |
| `archive.rs` | 560 | sealed result directories and the verifier |
| `validate.rs` | 490 | structural and domain validation |
| `sampling.rs` | 380 | random streams, strict sums, selection |
| `network.rs` | 360 | compiled network, checked updates |
| `cone.rs`, `conservation.rs`, `linkage.rs` | 760 | exact analysis |
| `python/crn_craft` | 1,010 | dictionaries API, observables, importers, CLI |

## 4. Numerical methods

### 4.1 Direct SSA

For state *x* with propensities *a_j(x)*, the Direct method (Gillespie 1977)
draws the waiting time *τ = −ln(u₁)/a₀* with *a₀ = Σ a_j* and the channel *j*
with *Σ_{k<j} a_k ≤ u₂ a₀ < Σ_{k≤j} a_k*.

Implementation details that make this exact and reproducible:

- **Random words.** Replicate *i* uses ChaCha8 (`rand_chacha` 0.3.1) with a
  256-bit key equal to SHA-256 of `"crn-craft.chacha8/0.1\0"` and the root
  seed, and stream number *i*. Uniforms are mapped from 53-bit words by a
  pinned mapping to (0, 1) or [0, 1). The sampler identity
  (`crn-craft.chacha8-direct/0.1`) is recorded in every result.
- **Strict ordered sums.** *a₀* is accumulated in declaration order. A
  positive rate whose addition does not change the sum (lost increment), a
  non-finite sum, or a positive minimum rate below 2^-53 of the total (a
  channel that `u₂ a₀` could never select) is a resolution failure
  (`CRN-E061`). The prefix sums of this pass are kept, and the channel is
  found by binary search over them (`partition_point(|s| s <= u₂ a₀)`), which
  returns exactly the channel of a linear scan, including at dyadic ties (P11).
- **Propensities.** Mass action uses falling factorials
  (*c ∏ x(x−1)…/order!*), with exact fast paths for orders 0–2 and a scaled
  fallback for extreme values (P01, P02). Hill and Michaelis–Menten laws are
  complete propensities with reactant gates.
- **Updates.** Firing is transactional with checked i64 arithmetic; a
  negative count (a scheduler defect) is `CRN-E043`, a species population
  limit stops the replicate (`CRN-E044`, `population_limit`).
- **Time.** Save points are filled before an event that would pass them; the
  clock must advance (`CRN-E062` otherwise); the final save is at the end time.

### 4.2 Dependency-cached Direct

`simulate-cached` keeps the propensities and, after firing reaction *r*,
recomputes only the dependents of *r* (reactions whose propensity reads a
species that *r* changes, including regulators of custom laws and reactant
gates). Because dependents are sorted, all prefix sums below the lowest
dependent are unchanged; the strict pass resumes from there, reusing the
recorded running minimum (P12). Every quantity that influences a decision is
computed with the same operations in the same order as in `simulate`, so the
two implementations draw the same words and choose the same channels. CI
runs the full SSA and CME validation campaigns with both and requires
identical reports; new tests compare them on every termination path
(failure, population limit, overflow, event limit) and with partial
resumption and resolution failures.

### 4.3 Ensembles and concurrency

Replicates are assigned round-robin to up to 64 workers. Each worker reuses
its state, rate, prefix and minimum buffers across replicates (P16) and sends
paths through a bounded channel (capacity 4 when the larger reservation fits,
else 1); the caller consumes replicates strictly in order and updates Welford
mean and variance sums. With one worker no thread is spawned (P06). Saved
states use one flat buffer per path (P15), which makes the memory reservation
formula an actual upper bound: before P15 a run with many save points
allocated about 375 MiB against a 241 MiB reservation.

Stopping is cooperative. A per-worker poller reads the shared cancellation
flag at every poll and the wall clock, cancel file and an optional external
callback every 64th poll, the file at most once per millisecond. The Direct
loop polls once per event and every 64 channels inside large recomputations
(P14).

### 4.4 Finite state projection

Following Munsky and Khammash (2006), `fsp` builds the reachable states within
per-species bounds by breadth-first search (up to `max_states`), with every
transition leaving the box sent to an absorbing sink. The probability vector
*p(t)* is propagated by uniformization: with *Λ ≥* the largest exit rate and
*P = I + Q/Λ*, *p(t+Δ) = Σ_k e^{−ΛΔ}(ΛΔ)^k/k! · p P^k*. Each save interval is
split so that *ΛΔ ≤ 16*; each series is truncated when the geometric tail
majorant falls below the per-step share of the requested tolerance. Reported
separately per snapshot: probability in the sink (physical leakage), the sum
of the series tail bounds, the omitted series weight (summed explicitly from
the first omitted term, not as one minus the kept weight), the total mass
defect and the roundoff residual. Probabilities are never renormalized.

Since P18 the generator is stored by column (CSC with the diagonal at its own
position and the sink last), and each series term is one pass that gathers
*(pP)_j = Σ_i p_i P_ij* in ascending *i* and adds it to the running series sum.
This is the same summation order as the earlier row-oriented scatter, so the
results are bitwise unchanged.

### 4.5 Tau leaping

`tau` implements explicit Poisson leaping. The leap *τ* is the largest step
for which the expected change and the standard deviation of every reactant
species stay below *max(ε x_i / g_i, 1)* (Cao, Gillespie and Petzold 2006,
with their factors *g_i* (eq. 27) for orders up to three), capped so that no
leap mean exceeds 10^7. Proposals that would make a count negative are
rejected and the leap halved; after 20 rejections, or when channels are
close to exhausting a reactant, the step falls back to exact Direct events.
Tau uses its own stream domain (`crn-craft.chacha8-tau/0.1`). Poisson variates
with mean at most 16 use the multiplication method; larger means use
Hörmann's PTRS (1993), an exact transformed-rejection sampler with about
1.1 pairs of uniforms per variate (P19). Results carry the algorithm
identifier `adaptive_explicit_poisson_cgp_reject_halve_direct_fallback/0.3`
and `assessment: approximate_ensemble`.

Accuracy note (measured, see §7.3): when ε exceeds the relative fluctuation
size of a species (about 1/√x near equilibrium), the leaps move the state by more
than its own fluctuations and the stationary variance is inflated: by 18% at
ε = 0.03 and 68% at ε = 0.1 for a mean of 1,000, and by a factor of about
1.7 at ε = 0.05 for a mean of 10,000. Means remain accurate.

### 4.6 Analysis

`conservation` computes an exact rational basis of the left nullspace of the
stoichiometric matrix with checked i128 arithmetic; `cone` enumerates the
extreme rays of the nonnegative conservation cone by a bounded double
description method (after Burton, arXiv:0808.4050) and derives upper bounds for
every species in the initial conservation class; `linkage` computes complexes,
linkage and strong linkage classes, terminal classes, deficiency and weak
reversibility. Every algorithm stops at explicit limits with an `unavailable`
status.

## 5. Interfaces and contracts

The command-line tool exposes `describe`, `check`, `analyze` (and
`analyze-v1`), `export`, `simulate`, `simulate-cached`, `verify`, `fsp` and
`tau`; all output is JSON, and exit codes are documented. Contract types
derive JSON Schemas (`crn-craft describe KIND`), and any shape change gets a
new schema version (result 0.2, analysis 0.2, discovery 0.5, FSP and tau 0.1).

The Python package (`crn-craft` 0.2.0.dev1, ABI3 wheel for Python ≥ 3.10)
mirrors the commands and adds observables, importers and a Python CLI. Since
this review, native calls run on a helper thread so that Ctrl-C cancels them
within about 50 ms, `simulate`, `fsp` and `tau` accept a `cancel_file`,
rejected requests raise `CrnError(ValueError)` with parsed diagnostics, and
`__version__` reports the installed version.

## 6. Validation

### 6.1 Test suites and gates

| Gate | Content |
|---|---|
| Rust tests | 166 tests in the core and CLI (plus 9 in the integrations) for contracts, validation, network, sampling, SSA, archives, FSP, tau, analysis, export and the CLI |
| Mutation gate | 45 source mutants, each required to turn a passing suite into named failures (35 reviewed, 8 campaign, 2 added here for the prefix-sum and resumed-sum checks) |
| `check_contracts.py`, `check_network.py`, `check_structure.py` | schemas, CLI behavior, exact/Decimal oracle on 200 seeded networks, 300-network cone/graph oracle |
| `check_ssa.py` | three predeclared 100,000-path analytic campaigns, an independent scalar RNG, thread identity, cancellation latency |
| `check_cme.py` | four 100,000-path exact master-equation checks (second-order, Hill, Michaelis–Menten), 12 DKW tests at family α = 0.001 |
| `check_fsp.py` | 16 checks including the pinned DSMTS stochastic test suite |
| `check_tau.py` | 30,000 analytic control paths at ε = 0.1, 0.05, 0.025 within predeclared envelopes |
| `check_python.py`, `check_importers.py`, `check_installed.py`, `check_tutorial.py` | binding/CLI identity, observables with synthetic controls and their gates, structured errors, cancel files and Ctrl-C, 48 importer checks, installed scripts, the tutorial's 17 runnable blocks |
| Dependency gate | the default build's dependency set and licenses |

Every Python checker runs normally and under `python -O`, and CI diffs the
reports. `simulate` and `simulate-cached` must give byte-identical SSA and CME
campaign reports.

### 6.2 Review of 1 November 2026

The supplied campaign patch was reviewed by independent reviewers for the
core, FSP/tau, the Python layer and the integrations. Fixed in this round:

- **Core:** a negative count was reported as a population-limit stop (now
  `CRN-E043`); the memory reservation was exceeded for many save points
  (fixed by P15, now tested); the resolution check used the first rather than
  the smallest positive rate in one path (tested); new tests cover every
  termination path of the cached implementation.
- **FSP:** the output budget was estimated, not enforced (the serialized
  size is now checked exactly, `CRN-E070`); the omitted series weight was
  computed as one minus the kept weight and lost to cancellation (now summed
  explicitly; tests check it stays below the tail bound).
- **Tau:** large populations exceeded the sampler's range (leap means are
  capped at 10^7); tau shared the Direct streams (now its own domain); the
  *g_i* factors are tested against the published cases; `--cancel-file` and
  help texts were added.
- **Python:** first-passage "lower bounds" were invalid when a path crosses
  between saved times (removed; only grid-supported bounds remain); the
  oscillation detector failed on constant series and accepted too few
  samples; importers accepted fixed species, missing stoichiometry, implicit
  units, `fast` reactions, colliding identifiers and pathological literals,
  and could hang or crash on adversarial input (all rejected with
  `CRN-E090` within one second); Ctrl-C was ignored; errors were unstructured;
  the installed-package check could test a different `crn-craft` on PATH.
- **ODE (`crn-diffsol`):** BDF trial states that overshoot slightly below
  zero aborted plain decays for horizons of 40 or more; rates are now
  evaluated at max(y, 0), and accepted or reported values must stay within
  one solver error weight of zero. Population limits are enforced
  (`CRN-E044`), kinetics diagnostics are no longer lost, and times near the
  largest double give a diagnostic instead of a panic. New dimerization, Hill
  and decay fixtures and Jacobian checks at distinct states fail under seven
  mutants that the old checks missed; `check_ode.py` grows from 16 to 50 rows.
- **Fitting and screening examples:** the fit is gated on its derivative
  report, its parameter order and a gradient check; the topology screen's
  gate, formerly constant, now has five declared criteria and is shown to
  fail on perturbed results. Its horizon was raised from 200 to 400 after
  observing an ambiguous classification at 200 (documented as post hoc).
- **Job worker:** the circuit-score penalty is a mean over paths (it grew
  with the replicate count), the amplitude species is named, save times up to
  a burn-in are excluded, and `seed_set` is required (problem version 2).
  `check_jobs.py` now checks every objective value against an independent
  recomputation and fails under five worker mutants.

### 6.3 Limitations

- Validation uses analytic and master-equation references for small networks
  and synthetic circuits; it does not establish biological validity of any
  model.
- Tau leaping is approximate (§4.5), and its envelopes cover three
  immigration-death controls only.
- The oscillation detector is calibrated on synthetic controls; its error
  rates on real data are unknown.
- Timings come from one shared 2-vCPU cloud host; absolute numbers vary by
  machine. The comparisons are paired and alternate binaries to control
  drift, but are not a cross-machine ranking.
- macOS and Windows wheels are built by a manual workflow and are not locally
  executed.
- The optional integrations (`crn-diffsol`, the job worker) are built and
  checked locally, not in CI: they need newer Rust toolchains, pinned git
  dependencies and the optimization-jobs runner. The ODE integration has no
  step or wall-time budget.

## 7. Performance

### 7.1 Method

Performance claims follow a predeclared protocol. Before any code or
measurement, a declaration fixes workloads, candidate hypotheses and the
acceptance rule. Each candidate is implemented with a test, built in release
mode, and timed in seven alternating before/after pairs (21 for W7, five for
W8b and W9) against the previous accepted binary, one worker, serial runs.
Exact candidates must reproduce every numerical digest. A candidate is
accepted if the median paired speedup on a target workload is at least 1.03
with the lower end of its 95% bootstrap interval (20,000 resamples, seed
20261101) above 1.0, and no frozen workload W1–W4 has an interval entirely
below 0.97. Rejected candidates are reverted, and their diffs and timings are
kept. Core time is the engine's own timer (simulation and output-budget
check); process time includes startup and archive writing.

| ID | Workload | Events or work |
|---|---|---|
| W1 | two channels (immigration–death), 64 paths, horizon 100 | 1,268,165 events |
| W2 | sparse network, 66 reactions | 1,234,576 events |
| W3 | three-gene circuit, 12 reactions | 3,663,350 events |
| W4 | dense network, 66 reactions | 1,225,923 events |
| W5, W6 | W3 with 2,001 save points; W6 keeps all 64 paths | 3,663,350 events |
| W7 | two channels, horizon 0.1 (latency) | 675 events |
| W8a | FSP birth–death, 401 states, horizon 40 | 87.1 M work units |
| W8b | FSP two species, 10,000 states, horizon 5 | 120.6 M work units |
| W9 | tau, immigration–death with mean 10^4, ε = 0.05, 256 paths | 25.6 M firings |

### 7.2 Iterations

P01–P10 (campaign, against the M2 baseline on W1–W4 and FSP):

| Iteration | Change | Result |
|---|---|---|
| P01 | fast paths for order 0–1 mass action | accepted, 1.01–1.15 |
| P02 | fast paths for dimer factors | accepted, 1.18 on a dimer control |
| P03 | skip a redundant state scan on trusted firings | accepted, 1.01–1.05 |
| P04 | one checked update for the selected channel | accepted |
| P05 | preallocate saved-path vectors | rejected (0.98–1.02) |
| P06 | no thread for one worker | accepted, 11.3 on short paths |
| P07 | bounded queues of four paths | accepted, 1.17 on eight-worker short paths |
| P08 | one pass for sum and minimum, one resolution check | accepted, 1.04–1.31 |
| P09 | cached Hill logarithms | rejected (layout regressions) |
| P10 | FSP rows normalized once, buffers reused | accepted, 1.19 |

P11–P20 (this review, against the previous accepted binary):

| Iteration | Change | Target result [95% interval] | Decision |
|---|---|---|---|
| P11 | binary search over recorded prefix sums | W2 cached 1.073 [1.070, 1.107] | accepted |
| P12 | resume strict sums at the lowest changed channel | W2 cached 1.589 [1.584, 1.678] | accepted |
| P13 | skip bitwise-unchanged rates | W1 0.927 [0.916, 0.945] | rejected |
| P14 | one scheduler poll per event | W1 1.095 [1.082, 1.132] | accepted |
| P15 | flat save buffer per path | W5 1.005, W6 1.003 | no speedup; kept as memory fix |
| P16 | reuse worker buffers across replicates | W7 1.041 [1.003, 1.075] | accepted |
| P17 | FSP: skip exact-zero rows | W8a 1.245, W8b 1.140 | accepted |
| P18 | FSP: column gather with fused accumulation | W8a 1.218, W8b 1.723 | accepted |
| P19 | tau: PTRS for leap means above 16 | W9 312 [296, 330] | accepted |
| P20 | archive: serialize once | W6 process 0.989 [0.987, 1.001] | rejected |

The full table, including the guard workloads, is in
[performance-p11-p20.md](campaign/performance-p11-p20.md).

### 7.3 Final numbers

Cost per event, median core time over the event count, one worker:

| Workload | 5a07334 `simulate` | final `simulate` | final `simulate-cached` |
|---|---|---|---|
| W1 two channels | 80.6 ns | 49.1 ns | 48.3 ns |
| W2 sparse 66 | 991.6 ns | 502.3 ns | 151.7 ns |
| W3 circuit 12 | 458.2 ns | 323.2 ns | 151.3 ns |
| W4 dense 66 | 902.7 ns | 659.5 ns | 553.3 ns |

Paired speedups of the final source:

| Workload | vs 5a07334 (`simulate`) | vs P10 (`simulate`) | vs P10 (`simulate-cached`) | cached vs plain (final) |
|---|---|---|---|---|
| W1 | 1.647 [1.429, 1.663] | 1.146 [0.952, 1.387] | 1.123 [1.085, 1.154] | 1.014 [0.999, 1.029] |
| W2 | 1.974 [1.956, 1.989] | 1.058 [1.054, 1.078] | 1.714 [1.679, 1.735] | 3.256 [3.241, 3.264] |
| W3 | 1.418 [1.404, 1.426] | 1.023 [0.994, 1.033] | 1.043 [1.043, 1.050] | 2.141 [2.047, 2.197] |
| W4 | 1.369 [1.320, 1.392] | 1.063 [1.054, 1.136] | 1.037 [1.027, 1.057] | 1.176 [1.164, 1.275] |
| W5 | 1.424 [1.400, 1.438] | 1.030 [1.018, 1.063] | 1.052 [1.046, 1.053] | 2.136 [2.121, 2.206] |
| W6 | 1.401 [1.396, 1.414] | 1.025 [1.013, 1.047] | 1.059 [1.051, 1.061] | 2.127 [2.115, 2.223] |
| W7 | 19.86 [18.57, 20.82] | 1.038 [0.998, 1.083] | 0.988 [0.955, 1.024] | 0.974 [0.949, 1.022] |

Combining the paired comparisons, the final `simulate-cached` runs the sparse
66-reaction network about 6.5 times and the circuit about 3.0 times faster
than the original `simulate` (ratios of median ns per event above).

FSP and tau against P10: W8a 1.424 [1.413, 1.444], W8b 2.340 [2.220, 2.396]
(probabilities identical; the review's FSP fixes changed only the `work`
count and the omitted-weight and roundoff fields), W9 335 [314, 343].

Where the time goes: for W1 (two channels) the cost is dominated by the two
random words, the logarithm of the waiting time and the bookkeeping per event;
caching cannot help. For W2 and W3 a reaction changes few propensities, so
`simulate-cached` saves most rate evaluations and, since P12, most of the
strict summation. W4 is dense: most reactions change most propensities, so
caching saves little. FSP is memory-bound in its sparse product; the column
form reads the generator sequentially once per series term. Tau's cost before
P19 was proportional to the leap means; after, it is proportional to the
number of leaps.

External reference (final source, M4a harness `compare_references.py`, same
host, five repeats, 64 paths to t = 100): rebop 0.10.0 through its Python API
took 39 ms on the two-channel workload against 62 ms of CRN-craft core time
(`simulate`) and 68 ms (`simulate-cached`), so rebop is about 1.6 times
faster there; on the sparse 66-reaction workload it took 2.21 s against
0.61 s (`simulate`, 3.6 times slower) and 0.19 s (`simulate-cached`, 11.8
times slower). The rebop times include its Python and result-array overhead,
CRN-craft core times include the ordered statistics, and the two use
independent random streams, so this is a cost comparison, not a path
comparison. Reports: `validation/campaign/p11-p20/cumulative/rebop-*.json`.

GillesPy2 1.8.3 on the same four workloads (details in the
[comparison](comparison-gillespy2.md)): its compiled C++ solver takes 87, 513,
170 and 617 ns per event on W1–W4, against 48, 154, 151 and 568 ns for
`simulate-cached`, so it is 1.1 to 3.3 times slower. It also needs 1.5–2.1 s
to compile each model. GillesPy2's NumPy solver is 43 to 126 times slower than
`simulate-cached`.

Catalyst 17 / JumpProcesses 10.0.1 on seven workloads (details in the
[comparison](comparison-catalyst-jumpprocesses.md)): after the alternative
exact methods were added, Julia's fastest aggregator was 1.54–2.39 times
faster than the fastest CRN-craft option (geometric mean 2.02). The follow-up
[Julia-gap campaign](campaign/julia-gap-implementation-log-2026-11-02.md)
sped DirectCR, SortingDirect and RSSACR up 1.07–2.4 times with unchanged
numbers or statistically gated paths; projected from those paired speedups, the remaining gap is
1.11–1.83 times (geometric mean 1.46). That projection is not a Julia
measurement. A later direct same-host measurement with the `fast` switch
([fast mode versus Julia](comparison-julia-8-threads.md)) gives a
geometric mean of best fast over best Julia of 1.15 at one worker and 0.69 at
eight workers on the original batches (1.12 and 1.04 on 256-path batches).
`fast` was opt-in then; it is the default of `ssa` requests since
`crn-craft.ssa-run/0.3`.

### 7.4 Memory

Each SSA run reserves, before executing, *(8 c + 24 t)(b w + k) + 24 c +
1024 r + w(2 MiB + 8 R) + 32 MiB* bytes for *c* = species × save times, *t*
save times, *b* buffers per worker (3, or 6 with four-path queues), *w*
workers, *k* kept paths, *r* replicates and *R* reactions, and rejects runs
above 256 MiB. P15 made this an upper bound for the path buffers, and a test
checks the formula for one and eight workers and the fallback to three
buffers. FSP reserves per state, edge and save time, and checks the exact
serialized size of its result against `max_output_bytes`.

## 8. Reproducing

```sh
cargo build --release --workspace --locked
cargo test --release --workspace --locked
python validation/benchmark_iterations.py --iteration NAME \
    --before OLD_BINARY --after NEW_BINARY --command simulate-cached \
    --workloads W1 W2 W3 W4 --targets W2 --output report.json
python validation/compare_tau_samplers.py --before OLD --after NEW --output tau.json
python examples/tutorial/make_assets.py --binary target/release/crn-craft \
    --output docs/images/tutorial
```

Raw timing reports, digests and logs: `validation/campaign/` (P01–P10) and
`validation/campaign/p11-p20/` (P11–P20 and the cumulative comparison).

## References

- D. T. Gillespie (1977). Exact stochastic simulation of coupled chemical
  reactions. *J. Phys. Chem.* 81, 2340–2361.
- M. A. Gibson and J. Bruck (2000). Efficient exact stochastic simulation of
  chemical systems with many species and many channels. *J. Phys. Chem. A*
  104, 1876–1889 (dependency graphs).
- B. Munsky and M. Khammash (2006). The finite state projection algorithm for
  the solution of the chemical master equation. *J. Chem. Phys.* 124, 044104.
- Y. Cao, D. T. Gillespie and L. R. Petzold (2006). Efficient step size
  selection for the tau-leaping simulation method. *J. Chem. Phys.* 124,
  044109.
- W. Hörmann (1993). The transformed rejection method for generating Poisson
  random variables. *Insurance: Mathematics and Economics* 12, 39–45.
- B. A. Burton (2010). Optimizing the double description method for normal
  surface enumeration. *Mathematics of Computation* 79, 453–484;
  arXiv:0808.4050 (double description details used for the cone
  enumeration).
- D. J. Bernstein (2008). ChaCha, a variant of Salsa20 (the stream cipher
  behind the random streams).
- A. Dvoretzky, J. Kiefer and J. Wolfowitz (1956), with Massart's constant
  (1990): the DKW band used by the statistical checks.
