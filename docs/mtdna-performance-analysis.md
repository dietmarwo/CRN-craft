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
- **`mtdna_direct` is the fast route.** On the two-CPU host it took about
  40 ns of CPU time per event. That is 3.9× less than Static Direct, 3.3×
  less than Static Direct with Ziggurat waits and sparse trackers, 5.3× less
  than next reaction and 14× less than the tree method.
- **A full study is short.** 10,000 cells for three years with the takeover
  stop take about 8.6 s on eight CPUs. A supplied Java implementation took
  14.9 s on the same host, and Static Direct 43 s.
- **Threads scale almost linearly.** Results do not depend on the thread
  count.
- **Most bookkeeping is cheap.**
  - Count trackers cost nothing measurable under `mtdna_direct`.
  - A first-passage tracker doubles the run time.
  - A bounded event log adds 70%.
  - A Python callback costs about 14 µs per event, a hundred times the
    simulation itself, so it is for a few cells only.
- **The approximate methods are not faster here.**
  - Tau leaping cannot take multi-event leaps at K = 1000, because the
    vacancy `K − N` holds only a few places. It ran 26× slower than
    `mtdna_direct`.
  - Hourly binomial steps cost about as much as exact Static Direct.
- **Interfaces add nothing.** The command line and Python took the same
  time, and binding new parameters costs 0.3 ms per candidate.

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
- **Metrics.** *Wall ns per event* divides wall time by all events of the
  ensemble. With two busy threads, the CPU time per event is twice that.

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

So, for `mtdna_direct` on this host:

> time ≈ cells × 2 · δ · K · T × 40 ns ÷ cores

## 2. Methods

1000 cells (200 for the tree method and tau), two count trackers, two
threads:

| Method | Kind | Wall s per 1000 cells | Wall ns per event | Relative |
| --- | --- | ---: | ---: | ---: |
| `mtdna_direct` | exact | 2.9 | 20.3 | 1.0 |
| Static Direct, Ziggurat + sparse trackers | exact | 9.6 | 66.9 | 3.3 |
| Static Direct | exact | 11.2 | 78.2 | 3.9 |
| binomial, 1 h steps | approximate | 11.3 | 79.0 | 3.9 |
| next reaction | exact | 15.4 | 107.2 | 5.3 |
| tree | exact | 40 (8.1 for 200) | 280.5 | 13.8 |
| tau, ε = 0.03 | approximate | 75 (15.1 for 200) | 519.3 | 25.6 |

Why the methods differ:

- **`mtdna_direct`** fuses the whole event into fixed code for the
  three-class structure.
  - It computes the eight hazards from three counts with one division.
  - It certifies resolution once at preparation instead of at every event.
  - Its stop and count trackers read the counts directly.

  It is the only route that refreshes no observables between reports.
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

Choose `mtdna_direct` for the three-class model, and check that
`metadata["sampler"]` names it. If it falls back, the result says why
(`CRN-R051`).

## 3. Threads

| Threads (2000 cells, `mtdna_direct`) | Wall s | Speedup |
| ---: | ---: | ---: |
| 1 | 12.47 | 1.0 |
| 2 | 5.87 | 2.1 |

Recorded with an eight-CPU quota and 20 threads, against one thread:

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

Static Direct, 20 cells, one thread:

| Observation | Wall s | Wall ns per event |
| --- | ---: | ---: |
| none | 0.37 | 127 |
| a Python callback that always continues | 39.5 | 13,550 |

Plain count trackers integrate between changes of their class, so they are
free. A first-passage predicate or an event log needs the general
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

Both front ends run the same native engine. Python returns NumPy columns
without copying the records. Preparing a candidate is negligible next to its
simulation either way. `bind_parameters` is preferred because it keeps the
settings identical, not because it is faster.

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
- **Rust generated loop.** *Rust only.* From Rust, the generated
  event loop runs Static Direct requests with byte-identical results. On
  the 128-cell probe, with one thread, Ziggurat waits and sparse trackers
  (`validation/mtdna/generated-loop.json`), the costs per event were:

  | Route | ns per event |
  | --- | ---: |
  | interpreted Static Direct | 112 |
  | Static Direct with the generated rates kernel | 94 |
  | generated event loop | 48 |
  | `mtdna_direct` | 38 |

  It brings Static Direct's identity close to the specialized sampler's
  speed. It also works for models `mtdna_direct` does not recognize, such
  as other stop conditions.

## 7. Comparison with other implementations

Recorded on the same eight-CPU host, 10,000 cells, three years, takeover
stop, eight workers, two runs each
(`validation/mtdna/java-gap.json`):

| Implementation | Median wall s |
| --- | ---: |
| CRN-craft Static Direct | 43.40 |
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

The two-CPU column uses 40 ns of CPU time per event, as measured here. The
eight-CPU column uses 48 ns, as recorded for the 10,000-cell run on that host.

| Study | Events | `mtdna_direct`, 2 CPUs | 8 CPUs |
| --- | ---: | ---: | ---: |
| 1,000 cells, 3 years, K = 1000 | 0.14 G | 3 s | 1 s |
| 10,000 cells, 3 years | 1.4 G | 29 s | 8.6 s (measured) |
| 10,000 cells, 30 years | at most 15 G | 5 min | 1.5 min |
| 100,000 cells, 3 years, K = 2000 | 29 G | 10 min | 2.9 min |
| the tutorial's two-stage fit (100 × 800 and 33 × 2400 cells) | 23 G | 7.5 min (measured) | 2.3 min |

Static Direct takes about 3.9× these times, and the tree method about 14×.
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
   - The three-class model: `mtdna_direct`, with the stop absent or written
     as `N > 0 && mutants >= f * N`.
   - Other thresholds: first-passage trackers, or kept paths.
   - Any other closed model: Static Direct with
     `{"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": True}`
     when its own replay identity is acceptable, plain Static Direct
     otherwise.
   - Lineage identities: the clone model with the tree method.
   - Variable population: `thinning`.
2. **Use as many threads as CPUs.** Results do not change.
3. **Use exact methods at K = 1000.** The approximate methods do not save
   time there.
4. **Choose observation by cost.** Count trackers are free; a first-passage
   tracker or an event log doubles the cost; a Python callback is for a few
   cells.
5. **Optimize with prepared templates.** Use `mtdna_direct` evaluations of
   a few hundred to a few thousand cells
   ([optimization tutorial](mtdna-optimization-tutorial.md)).
6. **For the largest studies, use Rust.** Call the engine directly, or use a
   generated loop for models outside `mtdna_direct`.

## Reproducing

```sh
python examples/mtdna/performance.py all target/release/crn-craft
python validation/benchmark_mtdna.py --binary target/release/crn-craft --output mtdna-timing.json
python validation/check_generated_loop.py target/release/crn-craft target/release/examples --output generated-loop.json
```

The Java comparison needs the supplied Java sources, which are not part
of the repository. See [mtDNA ensembles](mtdna-rules.md#reproduction).
