# E1: builder-built models as fast as model files

Task E1 of the 0.7.1 campaign, branch `engine/E1` (from `docs-0.7.1`,
`3d80c6d`). Goals, in order: profile arm E (the builder form of the guide's
three-class mtDNA model under `static_direct`); make it fast without
changing any result; let `mtdna_direct` certify the builder form under
`run/1.0`; gates; repeat the measurement with a new arm G.

## Status

- 2026-11-13: worktree, venv and log set up.

## 1. Profile of arm E (before any code change)

Setup: the documents that `crn.ensemble` writes for arms A, B and E
(`validation/merge/manual-ratios/builder_form.py`: 200 cells, seed 1, two
threads, takeover stop at 90%) run through the release CLI of this
worktree (0.7.0 code, built with `CARGO_PROFILE_RELEASE_DEBUG=line-tables-only`).
The CLI reproduces the guide's picture: A 0.58 s, B 1.44 s, E 5.23 s, the
same event counts as the 13 November measurement (28,536,814, 28,847,586,
28,779,334). Profiler: `perf record -F 2999 --call-graph dwarf` (perf
6.8.12, task-clock), under the CPU lock. E: 33,487 samples; B: 9,088
(ratio 3.68, as the wall times).

Arm E, self time by function (inlined code counted to its caller):

| Where | Share of E |
| --- | ---: |
| rate programs: `expr::run` 41.2%, `deterministic_call` (the `max`) 4.0%, `RateKernel::eval` 3.1%, the call site in `RateCache::value` 4.6% | 52.9% |
| Static Direct linear loop (`simulate_linear`, inlined: affected-channel loop, checks, running sums, selection, waits, count update) | 27.6% |
| stop predicate: `Predicate::run` 13.3%, `stop_values` 1.1%, `Predicate::eval` 1.0% | 15.4% |
| projections (`Network::project`) | 2.1% |
| ChaCha8 refill, `exp`, kernel | 2.0% |

Arm B for comparison: linear loop 68.4%, projections 11.2%, the shared
vacancy factor (`RateCache::global`) 9.4%, ChaCha8 4.3%, stop 5.5% (the
fused fraction form `N > 0 && mutants >= 0.9 * N`).

Where the time goes in E:

- **Expression evaluation dominates.** Every event changes one species,
  and every replication rule reads all three, so 6 of the 8 channels are
  recomputed per event (5 replication rules and the degradation of the
  changed class). Each replication program is about 25 interpreted steps
  (`max(0, K - W - M - A) / (W + sm * (M + A))` alone is 21), so about 127
  steps and 5 calls of `max` per event. The vacancy subexpression is the
  same in the 5 rules and is evaluated 5 times. The file form shares it:
  its rates lower to `RateKernel::Vacancy`, one cached base per event.
- **Dependency fan-out is already tight.** `affected` is built from the
  species each rule's rate reads; the 6 recomputed channels all read a
  changed species. A tighter graph saves nothing here.
- **Constant folding is already done:** `Program::compile` folds every
  parameter-only subtree (`1 - mu`, `sm`, `K`). `free * sw * (1 - mu) * W`
  parses left-associated, so `sw` and `1 - mu` cannot be merged without
  changing the rounding.
- **The stop predicate costs 15%.** The builder writes
  `W + M + A > 0 && M + A >= 0.9 * (W + M + A)`. It reads three
  observables, so the fused fraction form does not apply; the register
  bytecode runs 14 instructions per event, zeroes a 32-slot register file
  per call and evaluates `W + M + A` twice. Every channel affects the stop,
  so it runs after every event.
- **The loop itself** costs more in E than in B (27.6% of E is about 1.4 s
  of CPU against B's 1.0 s), mostly latency after the program calls
  returns (the hottest instruction is the multiply right after the call).
  It shrinks with cheaper rates.

Consequences for the plan: share the vacancy subexpression across the
programs of one model (computed once per event, by the same steps, so its
bits are unchanged); shorten the programs with fused steps that keep the
operation order (`acc op observable`, `f(c, acc)`); evaluate a stop that
the rate compiler accepts with the accumulator program and its shared
subexpressions instead of the register bytecode. No change of what is
computed, only of how often and with how many dispatches.
