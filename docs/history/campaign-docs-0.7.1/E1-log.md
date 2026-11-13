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

## 2. Byte-identical speedups

What changed (all in `crates/crn-craft/src/rules/expr.rs` and `flat.rs`;
no contract, request, result or sampler identity changes):

- **Shared subexpressions** (`SharedPrograms`, `RateCache::expression`). At
  lowering, every subexpression of the expression rates of one model is
  keyed by its own compiled steps (constants by their bits). One that occurs
  at least twice, has at least four steps and contains no subexpression that
  occurs more often gets a slot; the outermost such subexpressions are
  replaced by a `Shared` step. The rate cache computes a slot the first time
  a rate needs it in a state (NaN until then, like the cache's other
  entries) and resets it with them. A subexpression's steps set the
  accumulator first and use the stack only above its top, so the slot's
  value is the value the inline steps compute. On the builder model the
  vacancy term `max(0, K - W - M - A) / (W + sm * (M + A))` becomes one slot,
  computed once per event instead of five times. A first version shared the
  largest repeated term, `free * sw` (four rules), and `free` again for the
  fifth rule (`free * sm * A`), so the term was computed twice per event; the
  rule "no slot for a term that contains a more frequent one" fixed that
  (E 2.25 s to 1.87 s on the CLI).
- **Fused steps** that keep the operation and its operand order: `acc op
  observable` (`AddLoad`, `SubLoad`, `MulLoad`, `DivLoad`, `BinLoad`), `f(c,
  acc)` and `f(acc, c)` with dedicated `max`/`min` steps, `PushLoad`, and
  the four arithmetic `Bin` steps without the logic checks of `apply`. The
  stack depth that decides between a program and the tree is still counted
  for the unfused steps, so the set of expressions that keep the tree is
  unchanged.
- **Products**: a program that is a load or shared value times constants,
  then times observables (`delta * W`, `free * sw * (1 - mu) * W`), is
  evaluated by two loops without dispatch.
- **Stop predicates** compile to the rate program when they can (the value
  of the tree; only a NaN payload may differ, and a NaN stop is an error
  whatever its payload), and `T > 0 && P >= f * T` with `T`, `P` sums of
  observables has a fused form like the existing one for single
  observables. `Predicate::fraction` is unchanged, so generated loops are
  unchanged.
- Tried and removed: compiling programs into nested closures instead of the
  step interpreter. Five alternating runs: 1.896 s against 1.920 s median
  (1.3%), not worth 200 lines.

Byte identity, checked on the way: the CLI result documents of the arm-B
and arm-E requests equal those of the 0.7.0 code byte for byte after every
step. Tests added: fused steps are covered by the existing random-program
test (3,000 expressions, 9,000 programs); `predicate_bytecode_matches_tree`
also checks the stop's program form; new
`sum_fraction_stops_have_the_value_of_the_tree` (10 forms, 5 refused, 512
states each with NaN, infinities, fractions),
`shared_subexpressions_keep_the_bits_of_every_program` (400 random model
rate sets, 3 parameter sets, 5 states), `the_builder_vacancy_term_is_shared_once`,
and in `flat.rs` `shared_subexpressions_leave_every_rate_unchanged` (every
closed model of the repository and the builder model, 64 random states:
the cache's value has the bits of the rule's own program).

Informal CLI timings on this host during development (not the declared
measurement): B 1.32–1.57 s; E 5.23 s (0.7.0), 4.0 s (sharing), 2.65 s
(products, sum stop), 2.37 s (one dispatch level less), 2.25 s (closures),
1.87 s (no dominated slot, closures), 1.92 s (closures removed; median of
five).

What the profile says after the change: the vacancy slot is a chain of
dependent operations ending in a division, so its cost is mostly latency;
the closure experiment did not shorten it. The rest of the gap to the file
form is the per-channel product and the stop.

## 3. `mtdna_direct` for the builder form

Code: `crates/crn-craft/src/rules/mtdna.rs`, gate in `mod.rs` (`builder`
is `request.schema == RunVersion::V10`).

- `Kernel::compile` runs the 0.7.0 recognition unchanged first (now
  `file_form`); its results, reasons and their order are those of 0.7.0.
  Only if it refuses and the request is `run/1.0` does `builder_form` read
  the model again. A model without whole-rule propensities that both refuse
  keeps the file form's reason, so existing `run/1.0` fallbacks report what
  they reported.
- `builder_form` reads each channel's law: a per-instance channel through
  its rate kernel (as before), a whole-rule channel through `Reader::whole`,
  which splits the propensity into factors over `*` and `/` in any
  association and requires constants, exactly one factor equal to the
  parent's count, and for a birth exactly one `max(0, L)` (or `max(L, 0)`)
  in the numerator and one linear `D` in the denominator, with `L = K - W -
  M - A` and `D` a positive weighted sum. Linear forms are computed with an
  exactness check on every coefficient (TwoSum for sums, a fused residual
  for products), so `0.1 * M + 0.2 * M` is refused rather than rounded. The
  laws must agree (one `delta`, one `K`, one weight vector to the bit), and
  the shared checks of 0.7.0 (capacity, projections, range, resolution)
  then run unchanged (`certify`, the 0.7.0 tail).
- The stop: the file form's `N > 0 && mutants >= f * N`, or `T > 0 && P >=
  f * T` (also `T * f`) with `T`, `P` sums of integer-weighted observables
  equal to the total and the mutant count. All partial sums are integers
  below 2^53 (checked from `K` and the weights), so the stop is true in
  exactly the states where the fast path's threshold test is.
- Sampler identity: unchanged, `chacha8-u53-ziggurat/rules-mtdna-direct/1`.
  The builder form of the guide's model gets the file form's kernel to the
  bit (unit test), and a cell makes the same draws: same events, and `M + A`
  equal to `mutants` at every report time
  (`mtdna_direct_runs_the_builder_form_under_run_1_0`). On the guide's
  workload G has A's 28,536,814 events and 23 stopped cells, and its means
  of `W`, `M`, `A` equal A's `wild_type`, `neutral_mutants`,
  `advantaged_mutants` exactly (the mean of `M + A` against `mutants`
  differs by at most 4.3e-14, the rounding of a sum of means).
- Tests (written to fail on the 0.7.0 certificate: with `builder_form`
  stubbed to refuse, the three new certificate tests fail): the guide's
  model accepted under `run/1.0` and refused (`WHOLE`) otherwise, with the
  file form's kernel; five variants (other names, inlined parameters,
  reordered sums and `max(K - ..., 0)`, other species order, factors
  reordered and reassociated, `W * delta`), each also checked hazard by
  hazard against Static Direct's propensities at seven states; thirteen
  near-misses and six stops, each with its reason (changed weight in one rule, vacancy
  without a class, another capacity in one rule, `max(1, ...)`, missing or
  wrong count, a squared count, another death rate, a death without its
  count, an inexact weight sum, a fourth reaction, a capacity below the
  initial population, six stops). At the request level: the builder form
  runs the fast path under `run/1.0` with the file form's draws, and a
  `rules-model/0.4` copy under `rules-run/0.6` keeps the fallback with
  `CRN-R051` ("whole-rule").
- Exact-law check (declared, then amended before any comparison because
  the FSP of the first sizes exceeded its work limit): two builder models
  of capacity 20 (with and without the takeover stop), 200,000 cells each,
  against the 0.7.0 Static Direct sample and the FSP law: 120 comparisons,
  0 failures, max |z| 2.09; every chunk on the fast path. `check_laws.py`
  gained `--references`, `--expect-sampler` and manifest rows of form
  `rules` with `stop`; existing runs are unaffected.

## 4. Measurement

`validation/merge/manual-ratios/builder_speed.py` (arms A, B, E, F, G), as
declared: E 68.4 ns per event against B 48.7, ratio **1.40** (target 1.5,
met; 0.7.0: 186.7 and 3.60). E/A 3.6, E/F 2.6, F/A 1.4, G/A 0.97. A, B and E
reproduce the 13 November event counts and mean digests; the CLI result
documents of B and E, also with 200 kept paths, equal those of `3d80c6d`
byte for byte.
