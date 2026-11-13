# The `fast` switch: implementation log (2 November 2026)

Branch `julia-gap10-fast`, from main `add6cca`. Host: 2 vCPU Intel Xeon
(AVX2/AVX-512), Linux, Rust 1.97 stable, Python 3.13, one numerical worker.
Evidence: `validation/fast/`.

## 1. Import of Astra's C01–C10 campaign

The uploaded ZIP was truncated: it had no central directory and ended inside
`W5.samples.json.gz`. I recovered its entries from their local headers with
CRC checks. The root aggregate patch (`CRN-craft-julia-gap10.patch`, SHA-256
`6d8d88b2…2eb1`) was complete and matched `SHA256SUMS`. Only the `repo/`
snapshot copies were lost, and the patch carries those files too.

It applied cleanly to `add6cca` and reproduced the declared tree `135f330`.
It is imported unchanged as commit `08c0e9c`, keeping C05, C06 and C09. The
188 release tests pass. Astra's log is
[julia-gap10-implementation-log-2026-11-02.md](julia-gap10-implementation-log-2026-11-02.md).

## 2. `method.fast`

**Purpose.** An opt-in policy that emulates what JumpProcesses.jl does, so
that benchmarks compare the same work. Off by default.

| Aspect | Default | `fast` |
|---|---|---|
| Random words | ChaCha8 per replicate | Xoshiro256++ per replicate, seeded from a keyed ChaCha8 stream in its own domain |
| Totals (DirectCR, SortingDirect, RSSACR) | Pairwise tree, lazily checked | Leaves plus a running total (`+= new - old`); recomputed from the leaves after 65,536 updates or below 2^-30 of the recent peak |
| Prefix sums (Direct, cached Direct) | Checked ordered sums | Additions only |
| Loss and 53-bit resolution checks | Yes | No |
| Stop polls | Before every random word | Every 256 loop iterations and before rejection retries |

In both modes the selection rule itself is unchanged, and so are sampling
budgets, cancellation, checked state updates, output limits and replay
independent of the thread count.

**Contract.** `fast` is part of the request. Because that changes the
request's shape, the request becomes `crn-craft.ssa-run/0.2`:

- `fast` is rejected under 0.1 with CRN-E020.
- A 0.2 request is answered with `crn-craft.ssa-result/0.2`.
- 0.1 requests keep byte-identical results.
- Schemas are regenerated from the Rust types; discovery moves to
  `describe/0.7`.
- Migration notes are in [alternative-ssa.md](../guide/migrating.md).
- Python gains `ssa_spec(fast=...)`.
- `bench_ssa` and `measure_julia_gap.py` accept a fast flag
  (`NAME=PATH@fast`).

**Why off by default.** The checked policy guarantees things a result cannot
show after the fact:

- no positive channel silently loses its probability;
- totals carry no drift;
- a cancellation request is seen before the next random word.

This is explained in the [README](../../README.md#exact-ssa-methods-and-the-fast-switch).

## 3. Correctness evidence

**Default mode unchanged.** The default 210-comparison law report is
identical, check for check, to the earlier A08 gate. All workload results
under the default are unchanged.

**Fast law gate.** `check_alternative_ssa.py --fast` runs the same fixtures,
seeds, 100,000 paths, family alpha and DKW limit as the default gate, as a
separate family:

- 210 of 210 comparisons pass;
- the largest CDF error is 0.005255 (F2 heterodimer, RSSACR/inverse, t=2)
  against a limit of 0.008046;
- all 20 wrong-law power checks reject;
- the normal and `-O` reports are byte-identical;
- CI now runs this gate.

**Rust tests:**

- Xoshiro256++ reference vectors;
- the additions-only prefix;
- fast rejected under 0.1;
- fast identity, and replay across 1 and 3 threads for all five methods and
  both waits;
- `fast: false` under 0.2;
- sampling budget and cancellation in fast mode.

**Checkers.** `check_contracts.py` (69 checks) and
`check_alternative_contracts.py` pass normally and under `-O` with identical
reports. The latter now covers all twenty method, wait and fast combinations.

**fmt and clippy.** Both are clean.

## 4. Measurements

There is no acceptance rule here: `fast` is accepted on correctness, and the
timings below are measurements, not decisions.

**fast versus default, W1–W7.** Paired, five rotated rounds, three repeats
each (`fast-vs-default.json`).

| Method | Speedup range |
|---|---|
| Cached Direct | 1.39–4.26× (W5 4.26×, W6 3.61×, W7 3.19×) |
| DirectCR | 1.36–2.05× |
| SortingDirect | 1.09–1.62× |
| RSSACR | 1.22–1.37× |

Legacy cached rows are unaffected (0.99–1.00×). All rows were faster in 5/5
rounds, except W7 RSSACR at 4/5.

**GillesPy2 1.8.3 rerun with `fast`** (`gillespy2-comparison-fast.json`):

| | W1 | W2 | W3 | W4 |
|---|---|---|---|---|
| fast cached Direct, faster than SSACSolver by | 2.92× | 8.18× | 3.61× | 1.53× |
| best fast method, faster by | 2.92× | 9.05× | 3.61× | 9.19× (RSSACR) |
| checked default `simulate-cached`, faster by | 1.83× | 3.59× | 2.27× | 1.12× |

Validation against GillesPy2:

- The default samples repeat the 1 November run exactly.
- fast cached Direct against new GillesPy2 paths: all ten KS p-values are at
  least 0.18, and every |z| is at most 1.47.
- Against the exact Poisson laws, pooled chi-square p is 0.59 (W1) and 0.30
  (W2).

**Julia, projected.** Julia's hosts are still blocked in this session. Scaling
Astra's same-host CRN medians by the fast speedups measured here gives a
CRN/Julia ratio of 0.80–1.40× (geometric mean 1.02×, against 1.44× measured
for the default). CRN is ahead on W1, W2, W3 and W7, and behind on W4 (1.06×),
W5 (1.20×) and W6 (1.40×).

## Limitations

- **Hosts.** Timings come from one shared 2-vCPU host. The Julia projection
  combines two hosts and two API boundaries; a direct `fast` rerun of the
  Julia comparison is still needed.
- **Precision.** `fast` gives up the resolution guarantee. Models whose rates
  span more than about 16 orders of magnitude can silently lose their
  smallest channels.
- **SortingDirect order.** Julia keeps the learned SortingDirect order across
  a batch. `fast` does not, because that would make results depend on the
  thread schedule.
- **Archives.** Legacy `simulate`/`simulate-cached` archives have no `fast`
  mode; their contract is frozen.
