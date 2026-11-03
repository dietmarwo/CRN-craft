# CRN-craft mtDNA implementation log

> Contributed with the patch and kept as written, apart from this note and
> neutral wording for the external comparator. The
> [review](mtdna-a0-a6-review-log-2026-11-03.md) changed several results
> described here. Tau step control is new: the K = 20 tau rows below never
> leaped, and at K = 1000 tau runs exact events. Binomial steps use
> competing-risk probabilities on a fixed lattice; both approximate
> samplers are now `/2`. The `mtdna_direct` certification is per state and
> fallbacks carry `CRN-R051`. `ssa-run/0.4` keeps the 0.3 defaults, and
> discovery is `describe/0.10` (CLI `0.11`). Commit hashes are the
> contributor's local ones.

## Fixed-model performance follow-up

The opt-in `mtdna_direct` kernel brought the 10,000-cell eight-worker
median from **43.40 s to 8.62 s**, a **5.04×** speedup over
Static Direct. The supplied Java SSA took **14.86 s** on the same host and
worker count. The specialized run was **1.72×** as fast as Java
for this pooled-class workload.

| Exact SSA implementation | Two elapsed runs, seconds | Median, seconds |
| --- | ---: | ---: |
| CRN Static Direct | 45.53 / 41.26 | 43.40 |
| Supplied Java fixed SSA | 14.46 / 15.26 | 14.86 |
| CRN specialized Direct | 8.59 / 8.65 | 8.62 |

All rows used 25,920 hours, 720-hour reports, a 90% takeover stop,
retained per-cell paths and an eight-CPU quota. CRN retained two class
trackers and serialized about 20.9 MB of JSON. Java wrote about 6.9 MB
of clone-resolved reports and separate summary counts. The timings
include each program's native output format, and the CRN measurements
also include Python JSON decoding. The specialized run processed
**1,434,907,957** events. Java's birth count implies about
1.44 billion births plus degradations.

### Four implementation and measurement iterations

1. **Profile and establish the baseline.** Single-worker native clock samples
   attributed about 19% of Static Direct samples to channel rate evaluation,
   14% to expression evaluation, 12% to the logarithm and 9% to observable
   refresh. Java Flight Recorder found 392 of 407 execution samples in the
   fixed SSA loop's stacks, with many around its wait and reaction draw.
   Sample counts indicate where to investigate and have sampling error.
2. **Fuse the exact event loop.** Recognize the lowered three-class topology,
   common vacancy law and unit count changes at preparation. Maintain total
   and live count directly, compute the shared factor once, group deaths and
   wild-type birth outcomes, and record the original rule for each event.
   Require conservative floating-point range and resolution bounds. The
   opt-in method falls back to Static Direct when it cannot certify a closed
   model. It reports which sampler actually ran.
3. **Defer observation work.** Evaluate reporting projections only when needed.
   Track eligible single-class counts on changes to the class. Arbitrary
   tracker predicates, event logging and live callbacks use full event
   observation. The event stream is identical between the fused and deferred
   builds. Scalar integrals can differ at the last bits because the
   summation order changes. Python and Rust can bind a reusable candidate
   to a prepared template. Candidate binding still validates and lowers
   parameter-dependent rules once per candidate.
4. **Replace the remaining log hotspot.** The inverse log accounted for about
   41% of the lazy loop's native profile samples. Reused the existing
   Ziggurat sampler with a shared random-word interface. Kept ChaCha8 and
   explicit event, sampling, instance and cancellation controls. The new
   sampler has its own replay identity and uses different random words from
   Static Direct. In the final profile random-word sampling accounted for
   about 13% of samples and the rare transcendental paths under 1%. The
   remaining event loop mostly appeared as inlined code.

These 512-cell intermediate timings are single observations. The one-worker
column includes Python decoding for CRN. The eight-worker column reports
CRN child-process times to separate a Python decoding outlier. Java columns
are process wall times in both cases.

| Iteration | One worker, seconds | Eight workers, process seconds |
| --- | ---: | ---: |
| Static Direct | 13.174 | 2.317 |
| Fused rates and state | 5.954 | 0.989 |
| On-demand reports and count trackers | 3.264 | 0.540 |
| Shared Ziggurat | 2.615 | 0.448 |
| Java SSA | 4.359 | 1.442 |

The separate Rust phase probe took 7.47 s in
prepared execution, 0.34 ms in preparation
and 0.10 s writing output. Python decoding took
0.35 s. The instrumented Java copy spent
13.59 s from pool creation to all futures
collected and 0.81 s elsewhere in
the process. Worker-side report construction is included in both solver
phases. The main Java rows used the original compiled classes.
A smaller Static Direct phase probe with 512 cells measured
2.26 s in prepared execution and
0.005 s in output writing.

### Acceptance and limits

- The independent K=20 CME gate passed with 32,768 paths. The specialized
  eight-worker production distribution comparison against Static Direct
  had a largest CDF distance of 0.0119 against
  the predeclared 0.0449 bound across five final counts and first passage.
- Java and the specialized kernel recorded 1,016 and
  1,008 takeovers. Their report-binned CDF distance was
  0.0029 against 0.0407.
  Java records individual clone extinctions. The specialized kernel tracks
  extinction of a pooled class.
- Rust fmt, workspace Clippy with warnings denied, release tests, CLI and
  schema checks, dependency gate, Python wheel and Python rules checks passed.
  The validation suites produced identical evidence under normal Python and
  Python `-O`. Specialized checks covered one/eight-worker replay, scalar
  and full observer agreement, endpoint mutation probabilities, supported
  stop behavior, budget stops, cancellation and Static Direct fallback.
  The bound Python candidate matched a fresh preparation numerically.
- The packaged build replayed the benchmark output byte for byte on a
  512-cell case. A full packaged build recheck took
  8.73 s and reproduced the event and takeover totals.
- The first benchmark attempt recorded three valid timings before a NumPy
  boolean failed JSON serialization. The fixed harness repeated the gates.
  Timing and sampling evidence is in `validation/mtdna/java-gap.json`.
  The run used the checked-out GitHub base `a34484b`. A fresh remote fetch
  could not be verified from this environment.

The release suite passed 224 Rust tests. The cumulative patch, follow-up patch and ordered mail series independently reproduced tree `861d295e20fadf692f759c5e74c3d9a86ede62b7` at local head `6edc100c24ce1996ad386d2de25dd6c07b18945e`.

## Earlier A0–A6 campaign

3 November 2026. The patch implements A0–A6 from the attached mtDNA plan on GitHub main `a34484b2cb2ce28a4e6015962478217e7a7d8462`. The verified patch reconstructs tree `e383b6ed2cc6e0af8ab51358828845d5cecec141` at local head `850914cdedab20b1761379d73cb7d857d1b75213`.

Static Direct processes **2.99×** as many events per second as tree Direct with one worker and **3.01×** with twenty workers on the 200-cell fixture. The full Python run processed **1,492,964,721 events in 35.73 s** with twenty workers. Its one-worker numerical digest matched exactly after excluding the request hash.

## Iterations

| Iteration | Implementation and evidence |
| --- | --- |
| A0 | Added the three-class rules model at the production defaults and a K=20 reference fixture. The independent CME marginals passed for both the three-class and original clone models. |
| A1 | Added prepared Python rules ensembles over the existing native workers, cancellation, budgets and read-only NumPy record buffers. The 10,000-cell one- and twenty-worker runs completed and matched. |
| A2 | Added Rust observers, live Python callbacks and batch delivery from the bounded log. Trackers record crossings, episodes, extrema, integrals and first passage. A full event log reconstructed every tested tracker and rule count. |
| A3 | Lowered closed attributed instances to flat channels. Precomputed local factors and compiled constant and vacancy-rate kernels. The measured throughput exceeds the 2× target. |
| A4 | Added indexed Gibson–Bruck clocks with time reuse in rules and network SSA. The exact-law gates passed. On the 48-species sparse rules case, next reaction improved event throughput by 6.24× over static Direct. |
| A5 | Added CGP-style adaptive Poisson leaps, critical events, rejection, exact fallback and split-binomial steps. The low-copy and high-copy law gates passed. The high-copy log contains actual multi-firing leaps. |
| A6 | Added the ATP-coupled three-class model, analytic flow, continuous thinning and a deterministic ATP-zero stop. The time-varying Poisson and deterministic zero-crossing gates passed. |

Rules-run/result 0.2 adds algorithms and observations. Rules-model 0.2 adds the ATP block. Network next reaction uses ssa-run/result 0.4 with inverse waiting times and the checked policy. The Rust-generated schemas and migration guide are included.

## Same-host timings

The container had an **eight-CPU quota** (`800000 100000`). Twenty workers on this host share that capacity. Each main row is a single elapsed measurement. The CRN CLI rows include request preparation, the solver and JSON transport. The full Python rows time the prepared API call, including the returned arrays. Java rows include JVM startup, simulation and its report file. The fixed takeover and ATP comparisons keep per-trajectory reports on both sides. Java retains clone detail and CRN uses pooled classes.

| Workload | Method | Cells | Workers | Seconds | CRN events |
| --- | --- | ---: | ---: | ---: | ---: |
| fixed_no_stop | tree | 200 | 1 | 13.843 | 29,862,599 |
| fixed_no_stop | static_direct | 200 | 1 | 4.634 | 29,862,187 |
| fixed_no_stop | next_reaction | 200 | 1 | 5.551 | 29,862,755 |
| fixed_no_stop | tree | 200 | 20 | 2.347 | 29,862,599 |
| fixed_no_stop | static_direct | 200 | 20 | 0.780 | 29,862,187 |
| fixed_no_stop | next_reaction | 200 | 20 | 0.961 | 29,862,755 |
| sparse_48 | static_direct | 64 | 1 | 1.182 | 1,225,977 |
| sparse_48 | next_reaction | 64 | 1 | 0.189 | 1,224,253 |
| fixed_no_stop_python | static_direct | 10,000 | 20 | 35.731 | 1,492,964,721 |
| fixed_no_stop_python | static_direct | 10,000 | 1 | 219.254 | 1,492,964,721 |
| fixed_cox_stop | static_direct | 10,000 | 20 | 54.224 | 1,431,587,716 |
| fixedPop | java_ssa | 10,000 | 20 | 14.795 | — |
| fixedPop | java_hourly | 10,000 | 20 | 46.647 | — |
| fixed_cox_stop_retuned | static_direct | 10,000 | 20 | 48.448 | 1,431,587,716 |
| variable_atp | thinning | 512 | 20 | 2.815 | 77,829,452 |
| variablePop | java_ssa | 512 | 20 | 3.896 | — |

The observed twenty-worker exact Python run takes 35.73 s for the full no-stop workload. The previously reported 11 s remains a different-host figure with an unspecified solver. The supplied Java source defaults to hourly binomial steps. Here the explicit Java SSA and hourly modes took 14.79 s and 46.65 s respectively. Rebop’s 480 s remains a reported figure because the supplied archive contains the Java implementation.

The scalar tracker retime reduced the same retained-report case from
54.22 s to 48.45 s by skipping participant
allocation. It reproduced event, takeover, mutation and extinction counts.
The retimed case is 3.27×
the supplied Java SSA wall time on this host.

### Observation cost

Each row below is the median of three calls on the same 64-cell small fixture and one worker. The live callback decodes a committed event and performs no user work.

| Configuration | Events | Median ms | Individual ms |
| --- | ---: | ---: | --- |
| bare | 12,505 | 2.212 | 2.908, 2.020, 2.212 |
| tracker_log | 12,505 | 156.483 | 213.750, 156.483, 148.789 |
| live_callback | 12,505 | 110.496 | 118.187, 105.468, 110.496 |

## Java cross-checks

The fixed exact comparison used 10,000 cells, a 25,920-hour horizon, 720-hour reports and the 90% mutant-fraction stop. CRN recorded **1058** takeovers and Java **1016**. The report-binned first-passage CDF distance, including never-takeovers, was **0.00650** against the declared **0.04073** bound. The gate **passed**.

| Count, neutral / advantaged | CRN | Java |
| --- | --- | --- |
| Mutation events | 16,936 / 4,166 | 16,975 / 4,277 |
| Extinctions | 14,696 / 2,683 classes | 16,015 / 2,849 clones |

The extinction rows have different event definitions. The pooled model preserves class observables, and the clone model and complete event logs provide the identity needed for lineage reconstruction. Direct equality of lineage extinction counts remains untested.

For ATP, CRN stopped **60** of 512 cells at zero ATP and Java stopped **63**. CRN used continuous thinning and Java froze propensities over each sampled waiting interval. These runs provide a descriptive cross-check. The analytic Poisson law supplies the exactness gate for thinning.

Java ran from the supplied source on openjdk version "25.0.5.0.3" 2026-10-24 LTS with Apache Commons RNG 1.6. Source SHA-256 is `5f31ab5dc7d1e592e4db6d48fcdd318c03fe617d24f7e57e9751b45020c76c11`. The ZIP contains the CRN patch, the benchmark harness and summary evidence. The supplied source remains separate.

## Validation and retained failures

- **223 release tests** passed. Rust formatting, Clippy with warnings denied, and the Python wheel checks passed.
- **102 A0–A6 checks** and **47 legacy rules checks** passed. Normal and optimized Python reports matched byte for byte. CLI/schema, alternative SSA and dependency checks also passed with matching reports.
- The first tau attempts exceeded the unchanged CDF bound on the K=20 wild-type marginal at one hour. The solver now uses exact events for submolecule allowances. A separate high-copy test exercises the approximate path. The failed rows remain in `validation/mtdna/gates-attempt1.json` and `gates-attempt2.json`.
- The first production ATP run had one numerical failure in 512 paths.
  The shorter monotone envelope intervals passed the repeated full-horizon
  cohort and the analytic gates. The original diagnostic remains in
  `validation/mtdna/atp-production-attempt.json`.
- The observation timing request needed a larger explicit reservation
  than the default. The timing resumed from the completed replay and
  retained all 12,505 events in the corrected log benchmark.
- Raw results, timing and the 10,000-cell replay digest are in `validation/mtdna/`. Reproduction commands and supported model shapes are in `docs/mtdna-rules.md`.

## Patch series

```text
2de1e8c A0: freeze three-class mtDNA fixtures and validation targets
f4aabe1 A1: expose prepared parallel rules ensembles and NumPy record columns
fe20b6e A2: observe committed rule events with bounded trackers and callbacks
17e79dd A3: lower closed instance sets to static Direct channels
c990b73 A4: share Gibson–Bruck indexed clocks between rules and network SSA
51ff127 A5: add adaptive Poisson and split-binomial rules methods
32623c5 A6: add analytic ATP flow and certified thinning for static rules
755d622 Gate A0–A6 laws and Python callbacks, tighten leap budgets and document migration
1989bfa Preserve undefined vacancy arithmetic and record final statistical gates
850914c Tighten ATP envelopes and skip tracker-only payload allocation; record Java comparisons
```
