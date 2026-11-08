# mtDNA simulations: performance analysis

This analysis answers four questions for studies of mtDNA mutant
accumulation with the CRN-craft rules engine:

1. How long does a study take?
2. Which method should it use?
3. What do observation and the interfaces cost?
4. How does the engine compare with other implementations?

It combines new measurements on a two-CPU host with the earlier recorded
measurements in `validation/mtdna/`. Read the
[user manual](mtdna-user-manual.md) for the models and methods.

## Summary

- **Cost model.** A three-year cell at K = 1000 has about 150,000 events
  (2 · δ · K · T). Run time is close to cells × events per cell × cost per
  event ÷ cores.
- **The generic methods.** Static Direct, the generic exact method for
  closed models, took about 134 ns of CPU time per event on the two-CPU
  host with its default since rules-run/0.6, Ziggurat waits and sparse
  trackers, and 156 ns with inverse waits and per-event trackers, its
  default up to release 0.3.0. Next reaction took 214 ns and the tree
  method, which also runs models whose instances change identity, 561 ns.
- **The opt-in fast path.** `mtdna_direct`, certified at preparation for
  the three-class fixed-population structure only, took about 40 ns per
  event: 3.3× less than Static Direct (3.9× less than with inverse waits).
  Other models fall back to Static Direct with inverse waits and
  `CRN-R051`. The numbers marked `mtdna_direct` below hold only for that
  structure; multiply them by about 3.3 for Static Direct.
- **A full study is short.** 10,000 cells for three years with the takeover
  stop take about 8.6 s on eight CPUs with `mtdna_direct` and 43 s with
  Static Direct. A supplied Java implementation took 14.9 s on the same
  host.
- **Threads scale almost linearly.** Results do not depend on the thread
  count.
- **Most bookkeeping is cheap.** Measured under `mtdna_direct`:
  - count trackers cost nothing measurable;
  - a first-passage tracker doubles the run time;
  - a bounded event log adds 70%.
  - A Python callback costs about 14 µs per event, a hundred times the
    simulation itself, so it is for a few cells only.
- **The approximate methods are not faster here.**
  - Tau leaping cannot take multi-event leaps at K = 1000, because the
    vacancy `K − N` holds only a few places. It ran 6.6× slower than
    Static Direct.
  - Hourly binomial steps cost about as much as exact Static Direct.
- **Interfaces add little.** The command line and Python took the same
  time, binding new parameters costs 0.3 ms per candidate, and the pandas
  results add at most 3% to a simulation.

## Measurement setup

New measurements:

- **Host.** Linux, two CPUs, Python 3.13, `crn-craft` built from `main` at
  `ac91628`; its Python API and engine are those of release 0.2.0.
- **Script.** `examples/mtdna/performance.py`. Raw timings are in
  `validation/mtdna/user-scenarios-timing.json`.
- **Workload.** Unless stated, the three-class model with
  `N > 0 && mutants >= 0.9 * N`, three years (25,920 h), monthly reports,
  seed 17 and two threads.
- **Repeats.** Medians of three runs, or a single run for the slow
  methods.
- **Fast path.** Most tables time the opt-in `mtdna_direct`, which applies
  only to the three-class structure
  ([user manual, section 8](mtdna-user-manual.md#8-choosing-a-method)).
  Their headers say so. Section 2 gives the factor for every other method.
- **Metrics.** *Wall ns per event* divides wall time by all events of the
  ensemble. With two busy threads, the CPU time per event is twice that.
- **Static Direct defaults.** Rows labeled *0.3.0 default* ran inverse
  waits and per-event trackers, which `method="static_direct"` selected up
  to release 0.3.0. Rows labeled *0.6 default* ran Ziggurat waits and
  sparse trackers, which it selects now. Section 2 gives a separate,
  alternating measurement of the two
  (`validation/mtdna/static-direct-defaults-timing.json`).

Recorded measurements, cited where used:

- `validation/mtdna/java-gap.json` and `timing.json`: an eight-CPU quota, the
  supplied Java implementation, 10,000 cells.
- `validation/mtdna/generated-loop.json`: the Rust generated loop.

Timing is host-specific. Ratios between methods carry over better than
absolute times.

## 1. The cost model

Events per cell are set by the turnover. In the fixed-population model
degradations and replications balance, each at about δ · K per hour, so a
cell sees about 2 · δ · K · T events. For K = 1000, half-life 10 days and
three years, that is 2 × 0.00289 × 1000 × 25,920 ≈ 150,000 events. A
takeover stop ends about a tenth of the cells early. That brought 1000
cells to 143.8 million events instead of 150 million.

| K (500 cells, `mtdna_direct`) | Events | Wall s | Wall ns per event |
| ---: | ---: | ---: | ---: |
| 250 | 18.4 M | 0.33 | 17.8 |
| 500 | 36.7 M | 0.84 | 22.8 |
| 1000 | 72.4 M | 1.61 | 22.2 |
| 2000 | 139.4 M | 3.15 | 22.6 |

The time grows linearly with `K`, and the cost per event stays flat. The
same holds for the horizon and the half-life: the cost is in the number of
events, not in the size of the state.

So, on this host:

> time ≈ cells × 2 · δ · K · T × cost per event ÷ cores

with a CPU cost per event of about 40 ns for `mtdna_direct` and 134 ns for
Static Direct (156 ns with inverse waits; section 2).

## 2. Methods

1000 cells (200 for the tree method and tau), two count trackers, two
threads:

| Method | Kind | Wall s per 1000 cells | Wall ns per event | Relative to Static Direct, 0.3.0 default |
| --- | --- | ---: | ---: | ---: |
| Static Direct, 0.3.0 default (inverse waits, per-event trackers) | exact | 11.2 | 78.2 | 1.00 |
| Static Direct, 0.6 default (Ziggurat waits, sparse trackers) | exact | 9.6 | 66.9 | 0.86 |
| next reaction | exact | 15.4 | 107.2 | 1.37 |
| tree | exact | 40 (8.1 for 200) | 280.5 | 3.59 |
| binomial, 1 h steps | approximate | 11.3 | 79.0 | 1.01 |
| tau, ε = 0.03 | approximate | 75 (15.1 for 200) | 519.3 | 6.64 |
| `mtdna_direct`, opt-in fast path | exact | 2.9 | 20.3 | 0.26 (3.9× faster) |

**The 0.6 default.** rules-run/0.6 made the faster row the default of
`static_direct`; Python's `method="static_direct"` selects it. Measured on
6 November 2026 on another two-CPU host (Python 3.12.3, commit `099c08d`,
`validation/benchmark_static_direct_defaults.py`), the same 1000 cells
with three runs of each request in alternating order took:

| Static Direct request | Median wall s | Median wall ns per event | Speedup |
| --- | ---: | ---: | ---: |
| rules-run/0.5 default (0.3.0 default) | 7.00 | 48.8 | 1.00 |
| rules-run/0.6 default | 5.99 | 41.6 | 1.17 |

The declared target was 1.10×. The table above gives the same ratio,
1.17×, on its slower host. The new default samples the same process with
another use of the random stream, so its numbers differ cell by cell; it
passed the master-equation gate of the Ziggurat sampler
([log](campaign/static-direct-defaults-log-2026-11-06.md)).
`{"kind": "static_direct", "waiting_time": "inverse", "sparse_trackers": False}`
reproduces the numbers of release 0.3.0.

Why the methods differ:

- **Static Direct** evaluates the eight channels from a data-driven plan,
  checks resolution and refreshes the observables at every event. A profile
  of 18.4 million events put:
  - 33% in rate evaluation;
  - 19% in the event loop;
  - 18% in selection with the waiting time;
  - 13% in random words;
  - the rest in state updates, observable refresh and the stop predicate.

  No part dominates. That is why interpreted rewrites could not close the
  gap
  ([record](decisions/generic-rules-structural-2026-11-03.md)) and a
  generated loop could (section 6).
- **Next reaction** keeps one clock per channel. Every event changes `N` and
  `D`, so every replication clock must be rescaled. That is more work than
  Direct's fresh draw.
- **The tree method** recomputes every rule over a sum tree of instances.
  It is needed only for models whose instances change identity, such as the
  clone model.
- **Tau** controls every rate's change within a leap. The vacancy term
  changes by a large relative amount with each event when only a few places
  are free, so tau runs exact events at extra cost. It pays off only where
  every rate changes slowly.
- **Binomial steps** resolve each hour by competing-risk probabilities per
  parent. That is about as much work as exact simulation at this event
  density, and the result is approximate.
- **`mtdna_direct`**, the opt-in fast path, fuses the whole event into
  fixed code for the three-class structure.
  - It computes the eight hazards from three counts with one division.
  - It certifies resolution once at preparation instead of at every event.
  - Its stop and count trackers read the counts directly.

  It is the only route that refreshes no observables between reports. It
  samples the same process as Static Direct with a different use of the
  random stream, so results agree in distribution, not cell by cell.

For a model with the three-class structure, opt in with `mtdna_direct` and
check that `metadata["sampler"]` names it. Preparation certifies the model
and each parameter set; if a check fails, the run is Static Direct and a
`CRN-R051` warning says why. The
[user manual](mtdna-user-manual.md#the-opt-in-fast-path-mtdna_direct) lists
the checks. Every other closed model runs on the generic methods above.

## 3. Threads

| Threads (2000 cells, `mtdna_direct`) | Wall s | Speedup |
| ---: | ---: | ---: |
| 1 | 12.47 | 1.0 |
| 2 | 5.87 | 2.1 |

Recorded with an eight-CPU quota and 20 threads, against one thread
(Static Direct with its 0.3.0 default):

| Method (200 cells, no stop) | 1 thread | 20 threads on 8 CPUs | Speedup |
| --- | ---: | ---: | ---: |
| Static Direct | 4.63 s | 0.78 s | 5.9 |
| next reaction | 5.55 s | 0.96 s | 5.8 |
| tree | 13.84 s | 2.35 s | 5.9 |
| Static Direct from Python, 10,000 cells | 219.3 s | 35.7 s | 6.1 |

Cells are independent, so threads scale until the CPUs are busy. Each cell
has its own random stream, so results are identical for any thread count.
Use as many threads as CPUs.

## 4. Observation

1000 cells, `mtdna_direct`:

| Observation | Wall s | Relative |
| --- | ---: | ---: |
| none | 3.00 | 1.00 |
| two count trackers (`neutral_mutants`, `advantaged_mutants`) | 2.90 | 0.97 |
| a first-passage tracker (`N > 0 && mutants >= 0.6 * N`) | 6.01 | 2.0 |
| an event log of mutation events, 50 per cell | 5.19 | 1.7 |

Static Direct, 0.3.0 default, 20 cells, one thread:

| Observation | Wall s | Wall ns per event |
| --- | ---: | ---: |
| none | 0.37 | 127 |
| a Python callback that always continues | 39.5 | 13,550 |

Under `mtdna_direct`, plain count trackers integrate only when their class
changes, so they are free; Static Direct does the same by default since
rules-run/0.6 (`sparse_trackers`). A first-passage predicate or an event log needs the general
observation path, which refreshes the observables and evaluates the
predicate after every event. That costs about the same again as the
simulation.

A Python callback receives every event as a parsed JSON document and holds
the interpreter lock, which serializes the threads. At 14 µs per event, one
three-year cell costs about 2 s. It suits exploring a handful of cells.
For whole ensembles, use trackers, event logs or a Rust observer.

## 5. Interfaces and parameter binding

| Operation | Time |
| --- | ---: |
| CLI `crn-craft rules`, 1000 cells, `mtdna_direct` | 2.90 s |
| Python `rules.prepare(...).simulate()`, the same request | 2.81 s |
| `template.bind_parameters({...})` per candidate | 0.28 ms |
| a fresh `rules.ensemble(...)` per candidate | 0.29 ms |

Both front ends run the same native engine. Preparing a candidate is
negligible next to its simulation either way. `bind_parameters` is
preferred because it keeps the settings identical, not because it is
faster.

Python parses the result document and returns pandas DataFrames. For 4000
cells with every path kept (`keep_paths="all"`, `mtdna_direct`, two
threads, pandas 3.1.1; `validation/mtdna/python-results-timing.json`,
medians of three):

| Step | 37 monthly reports | 155 weekly reports |
| --- | ---: | ---: |
| simulation (native) | 7.56 s | 7.92 s |
| parsing the result JSON | 0.083 s | 0.35 s |
| records, rule counts and statistics as DataFrames | 0.0012 s | 0.0013 s |
| `result.paths` (148,000 or 620,000 rows) | 0.039 s | 0.19 s |
| `result.to_xarray()` | 0.004 s | 0.015 s |

The DataFrames add at most 3% to the simulation, almost all of it for
`paths`, which is built on first use. That costs about as much as the
`np.array` over the path lists that code for release 0.2.0 needed (0.034 s
and 0.17 s). The first `to_xarray()` call also imports xarray, about
0.06 s. With a mutation log, `result.event_columns()` gave the 8626 logged
events in 0.012 s.

## 6. Other models and routes

| Workload | Wall s | Events | Wall ns per event |
| --- | ---: | ---: | ---: |
| clone model, tree, K = 1000, 100 cells, no stop | 2.65 | 14.9 M | 177.6 |
| three-class model, tree, the same | 4.63 | 14.9 M | 310.2 |
| ATP model, thinning, 512 cells | 13.18 | 77.8 M | 169.3 |

- **Clone model.** At K = 1000 the clone model is cheap: it has three rules
  instead of eight, and at this mutation probability few clones are alive
  at once. The
  tree method recomputes every rule per event, so fewer rules mean less
  work. Exact per-lineage bookkeeping is therefore affordable for
  thousands of cells. Event logs for it are limited by memory, not time.
- **ATP model.** Thinning integrates the ATP flow between events and
  rejects a share of candidate events. Recorded on eight CPUs, 512 cells
  took 2.8 s, against 3.9 s for the Java implementation of the same model.
- **Generated loop.** *Opt-in.* From Rust, the generated
  event loop runs Static Direct requests with byte-identical results. It
  needs constant or vacancy rates and sum observables, and is compiled for
  one model and parameter set; other requests are refused with `CRN-R020`.
  On
  the 128-cell probe, with one thread, Ziggurat waits and sparse trackers
  (`validation/mtdna/generated-loop.json`), the costs per event were:

  | Route | ns per event |
  | --- | ---: |
  | interpreted Static Direct | 112 |
  | Static Direct with the generated rates kernel | 94 |
  | generated event loop | 48 |
  | `mtdna_direct` | 38 |

  It brings Static Direct's identity close to the specialized sampler's
  speed. It also works for models `mtdna_direct` does not certify, such
  as other stop conditions. From Python, `rules.compile_loop` builds the
  same loop with a local Rust toolchain (new in 0.4.0; user manual,
  section 8). On the 1000-cell workload of section 2, end to end from
  Python, it took 18.1 against 40.1 ns per event under the Static Direct
  default (Ziggurat waits, sparse trackers; 2.21×) and 26.4 against
  49.1 ns with inverse waits and per-event trackers (1.86×), on the
  two-CPU host
  (`validation/mtdna/compiled-loop-python.json`). Both routes were timed
  together on 6 November 2026, at a later commit than the table of
  section 2, so compare the ratio, not the absolute times.

## 7. Comparison with other implementations

Recorded on the same eight-CPU host, 10,000 cells, three years, takeover
stop, eight workers, two runs each
(`validation/mtdna/java-gap.json`):

| Implementation | Median wall s |
| --- | ---: |
| CRN-craft Static Direct, 0.3.0 default | 43.40 |
| supplied Java SSA | 14.86 |
| CRN-craft `mtdna_direct` | 8.62 |

The specialized sampler ran 1.72× faster than the Java SSA. The timings
include each program's own output: CRN-craft serialized about 20.9 MB of
JSON with two class trackers, and Java wrote about 6.9 MB of reports.

The cross-checks agreed:

- **Takeover counts:** 1008 of 10,000 cells for CRN-craft, 1016 for Java.
- **Takeover-time distributions:** a distance of 0.003, within the declared
  bound of 0.041.
- **Mutation counts:** within 1% for both classes (17,030 and 4,236
  against 16,975 and 4,277).
- **Against Static Direct:** every observable and the first passage of
  `mtdna_direct` agreed in distribution within the family bound.

The checks are in the same file.

Earlier recorded rows with 20 threads on the eight-CPU quota show the
other methods:

- the Java hourly binomial scheme took 46.6 s;
- Static Direct with the takeover stop took 54.2 s before tuning and 48.4 s
  after.

A rebop implementation of the three-class model was reported to take 480 s
for 10,000 cells, because it computes trajectories one after another. It
was not rerun here. At its reported 11 s with 20 threads, the Java version
was 44× faster than that. `mtdna_direct` on eight CPUs is faster than both.

## 8. Sizing a study

Every column is for the fast path `mtdna_direct`. The two-CPU column uses
40 ns of CPU time per event, as measured here. The eight-CPU column uses
48 ns, as recorded for the 10,000-cell run on that host.

| Study | Events | `mtdna_direct`, 2 CPUs | `mtdna_direct`, 8 CPUs |
| --- | ---: | ---: | ---: |
| 1,000 cells, 3 years, K = 1000 | 0.14 G | 3 s | 1 s |
| 10,000 cells, 3 years | 1.4 G | 29 s | 8.6 s (measured) |
| 10,000 cells, 30 years | at most 15 G | 5 min | 1.5 min |
| 100,000 cells, 3 years, K = 2000 | 29 G | 10 min | 2.9 min |
| the tutorial's two-stage fit (100 × 800 and 33 × 2400 cells) | 23 G | 7.5 min (measured) | 2.3 min |

Static Direct, for any other closed model, takes about 3.3× these times
(3.9× with inverse waits), and the tree method about 14×.
Add the observation costs of section 4. Long horizons need a larger event
budget: `max_events` is per cell and must exceed 2 · δ · K · T.

Output and memory grow with what is kept, not with the events:

- **Kept paths:** cells × report times × observables × 8 bytes. For example,
  10,000 cells × 37 reports × 6 observables is 18 MB.
- **Event logs:** up to `max_entries` per cell. A request is admitted only
  if its worst-case reservation fits `max_output_bytes` (default 64 MiB).
  For large logs, raise it or log fewer cells.
- **Statistics and records:** small. Statistics accumulate while the
  ensemble runs.

## 9. Recommendations

1. **Pick the method by model.**
   - Any closed model: Static Direct, `method="static_direct"`, whose
     default uses Ziggurat waits and sparse trackers. Add
     `"waiting_time": "inverse", "sparse_trackers": False` to a method
     dictionary only to reproduce results of release 0.3.0.
   - Lineage identities: the clone model with the tree method.
   - Variable population: `thinning`.
   - The three-class structure: opt in with `mtdna_direct`, with the stop
     absent or written as `N > 0 && mutants >= f * N`, and check the
     sampler. Record other thresholds with first-passage trackers, or kept
     paths, so that the stop keeps the certified form.
2. **Use as many threads as CPUs.** Results do not change.
3. **Use exact methods at K = 1000.** The approximate methods do not save
   time there.
4. **Choose observation by cost.** Count trackers are free; a first-passage
   tracker or an event log doubles the cost; a Python callback is for a few
   cells.
5. **Optimize with prepared templates.** Use evaluations of a few hundred
   to a few thousand cells, with `mtdna_direct` where the model has its
   structure ([optimization tutorial](mtdna-optimization-tutorial.md)).
6. **For the largest studies, use Rust.** Call the engine directly, or use
   the opt-in generated loop for closed models with constant or vacancy
   rates that `mtdna_direct` does not certify.

## Reproducing

```sh
python examples/mtdna/performance.py all target/release/crn-craft
python validation/benchmark_mtdna.py --binary target/release/crn-craft --output mtdna-timing.json
python validation/check_generated_loop.py target/release/crn-craft target/release/examples --output generated-loop.json
```

The Java comparison needs the supplied Java sources, which are not part
of the repository. See [mtDNA ensembles](mtdna-rules.md#reproduction).
