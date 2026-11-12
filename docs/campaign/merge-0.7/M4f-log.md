# M4f log: per-event cost, the RSSA-CR clock and `auto`

Task M4f, added by the coordinator before the speed gate M4d of
[the plan](../../plan-0.7.0-engine-merge.md): bring the unified engine's
per-event cost within the 1.25 × gate of M4d, fix the clock failure that the
M4c RSSA-CR law gate met on vilar-noise, and revisit `auto`. The gate itself
is in the [M4d log](M4d-log.md).

## Summary

| Change | Results | Commits |
| --- | --- | --- |
| A positive wait below half an ulp of the clock is accepted under `crn-craft.run/1.0` (Static Direct, every selector, generated kernels and loops); 1,000,000 consecutive such steps fail with `CRN-R030` | changed only where M4c failed | `5d99986` |
| Rate expressions compiled to flat programs with folded parameters; dedicated programs for the converted Hill and Michaelis–Menten laws | byte-identical | `91d06f3`, `6669e5b`, `9562fe2` |
| `crn-craft rules` writes the result with sorted keys without a `serde_json::Value` tree; other commands through a `BufWriter` | byte-identical output | `efee994`, `cfafa58` (sorted writer, second form) |
| Linear Static Direct in its own loop, monomorphized on the waiting time: running totals resumed at the first changed channel, a lower bound for the resolution check, a branch-free scan of up to four rates | byte-identical | `cfafa58`, `6669e5b`, `9562fe2` |
| Count factors precomputed per channel; error paths out of line; one-pass count updates with undo | byte-identical | `39601ed`, `778fefe` |
| RSSA-CR: range caches shared by rules with one rate program; branch-free interval test with an out-of-line refresh; corners as an enum | byte-identical | `f4b894f`, `778fefe` |
| Lemire draws compute their division only when needed | byte-identical | `9562fe2` |
| RSSA-CR reuses the composition's accepting threshold instead of drawing another uniform | RSSA-CR changes (law gate rerun) | `dccc6da` |
| Selector tables reused between the replicates of a worker | byte-identical | `bbbef17` |
| `auto` refitted with population as a feature | default selector changes for many run/1.0 requests | `9c4b5f9` |

"Byte-identical" was checked on 156 CLI requests after every step
(`/home/claude/m07/m4f/harness.py`): the twelve speed-gate models with every
selector and `auto`, both waiting times, `next_reaction`, and a variant with
kept paths, kept instances and nine report times; plus FSP, a refused
request, an event-limit run with quantiles and trackers, and an output-limit
run for the CLI writer. The I1 stored requests and the generated-loop
replays (interpreter against checked-in compiled loops) ran after each step
too.

## 1. The RSSA-CR clock (coordinator's decision)

M4c's RSSA-CR gate failed on vilar-noise: a positive, finite wait was below
half an ulp of t = 7.16, so t + wait == t (`CRN-R030` "clock failed to
advance").

- `engine::Clock` (new): `advance(time, wait)` returns `time + wait`. When
  the sum does not exceed `time`, every `rules-run` version fails as before
  ("clock failed to advance"); under `crn-craft.run/1.0` a positive finite
  wait is accepted, the step happens at the same representable time, and
  only the 1,000,001st consecutive such step fails ("clock failed to advance
  in 1000000 consecutive steps"). A zero wait always fails.
- Used by the linear loop (inverse and Ziggurat waits), generated rates
  kernels (same loop), generated event loops (`LoopHost::next_time`, so
  `loop.verify` stays byte-identical) and composition–rejection, sorting
  and RSSA-CR (a step is an event or an RSSA-CR trial). `next_reaction`,
  tau, binomial, hybrid, thinning and `mtdna_direct` keep their strict
  checks (they are not selectors of Static Direct).
- The wait is now computed first and then added (`time + wait`); for the
  inverse wait `time + (-log(1-u)/total)` is exactly the former
  `time - log(1-u)/total`, so results are unchanged where the clock
  advances.

**Why 1,000,000.** A step leaves the clock unchanged with probability
p ≈ 1 − exp(−a·ulp(t)/2) for total rate a. At the vilar-noise failure
a ≈ 9,729 (RSSA-CR's upper total) and ulp(7.16)/2 = 4.4e-16, so p ≈ 4e-12:
two stalls in a row have probability ≈ 2e-23, and a million are impossible.
A run reaches the limit with probability above 1e-9 only if
p^1,000,000 ≥ 1e-9, i.e. a·ulp(t)/2 ≥ 10.8: more than about 22 events
are expected per representable step of the clock. There the time axis no
longer resolves the dynamics, which is the runaway the check is for. The
error of an accepted stall is below half an ulp of t (up to half an ulp at
an exact tie), the same as the rounding of any other addition to the
clock. Detecting a true runaway costs 10⁶ steps, 0.05–0.5 s, and the event
and sampling budgets bound the run as before.

Tests written first (`crates/crn-craft/tests/clock.rs`; 2 of 3 failed on
the unchanged engine, the third guards I1): a burst of 50 molecules that
decay at 1e20 each completes under run/1.0 for every selector and both
waits (51 events) and fails under rules-run/0.2, 0.4 and 0.6 with the 0.6.0
message after one event; an A ⇄ B flip at 1e20 fails under run/1.0 after
exactly 1,000,001 events for every selector. A generated-loop replay
(`generated_loops_keep_the_clock_rule_of_the_request_version`, two new
fixtures `loop_clock_stall` and `loop_clock_runaway`) failed first as well:
the compiled loop must agree with the interpreter under both rules.

## 2. Per-event cost

Profiles: callgrind (valgrind 3.22) on a CLI built with line tables in a
separate target directory, one process at the speed-gate sizes; `perf` is
not installed. Instruction counts (Ir) are deterministic, so they decided
each step; wall time was checked with alternating runs of the 0.6.0 best
and the unified best (`abtime.py`, 5–9 pairs), because single timings on
this host varied by up to 50 % while another agent was active.

Instructions of the 0.6.0 best (fast variant, xoshiro256++) and of the
unified engine's best method, before and after M4f:

| Model | 0.6.0 best (M0c method) | Ir | Unified before M4f (best profiled) | Ir | Unified after M4f | Ir | After / 0.6.0 |
| --- | --- | ---: | --- | ---: | --- | ---: | ---: |
| W1 | `cached_direct` | 703 M | `linear` | 815 M | `linear` | 523 M | 0.74 |
| W2 | `direct_cr` | 1,054 M | `rssa_cr` | 1,230 M | `rssa_cr` | 861 M | 0.82 |
| W3 | `cached_direct` | 2,575 M | `linear` | 5,585 M | `linear` | 2,433 M | 0.94 |
| W4 | `rssa_cr` | 1,429 M | `rssa_cr` | 3,550 M | `rssa_cr` | 1,332 M | 0.93 |
| W5 | `direct_cr` | 660 M | `rssa_cr` | 712 M | `rssa_cr` | 503 M | 0.76 |
| W6 | `rssa_cr` | 701 M | `rssa_cr` | 753 M | `rssa_cr` | 693 M | 0.99 |
| W7 | `rssa_cr` | 171 M | `rssa_cr` | 189 M | `rssa_cr` | 143 M | 0.84 |
| enzyme-mechanism | `cached_direct` | 766 M | `linear` | 837 M | `linear` | 506 M | 0.66 |
| immigration-death | `cached_direct` | 732 M | not profiled | — | `linear` | 461 M | 0.63 |
| michaelis-menten | `direct` | 574 M | `sorting` | 745 M | `linear` | 368 M | 0.64 |
| repressilator | `sorting_direct` | 3,413 M | `linear` | 6,244 M | `linear` | 2,746 M | 0.80 |
| vilar-noise | `cached_direct` | 12,494 M | `linear` | 16,557 M | not profiled | — | — |

Ir counts the whole CLI process (start, run, result output). "Before M4f" is
the M4c engine (`head-ir.txt`), "after" the engine of `57cd574` (`ir-final.txt`);
each column takes the method with the fewest instructions among those
profiled. Fewer instructions did not always mean less wall time: W2, W3, W5
and W7 run fewer instructions than 0.6.0 but are 12–20 % slower in the M4d
rerun. The unified engine's ChaCha8 words, the remaining data-dependent
branch of the linear scan and memory traffic are the likely causes; they were
not separated. vilar-noise (16.6 G instructions before M4f) was not
profiled again.

Where the time went and what was done:

- **CLI output (examples, 10,000 replicate records).** `to_value` plus the
  drop of the `Value` tree and an unbuffered `StdoutLock` (a `memrchr` per
  token) took about 37 % of immigration-death. `sorted.rs` (CLI) serializes
  once into a buffer per object and copies members in key order: the bytes
  of `to_writer(to_value(x))` (tested against it on structs, enums, maps
  with numeric and repeated keys, flatten, bytes, non-finite numbers).
- **Expressions (W3, repressilator, W4 bounds).** The tree evaluator
  recursed with a generic draw argument and allocated a vector per function
  call. `expr::Program` compiles once per lowered network: parameters are
  constants, parameter-only subtrees are evaluated by the tree evaluator at
  compile time, every other operation keeps its operands and order (`&&`,
  `||` and `if` evaluate both sides, which is value-identical without side
  effects); `^2` is `x * x` because `libm::pow` returns exactly that for the
  exponent 2. The converted Hill and Michaelis–Menten forms are recognized
  as step sequences and evaluated without dispatch. A randomized test
  compares 9,000 programs (and the two law shapes at 2,000 points each) with
  the tree; NaN payloads may differ (Rust does not specify them; a NaN rate
  is an error whatever its payload).
- **Linear Static Direct loop.** Running totals: the total is the
  sequential sum from −0.0 (`Iterator::sum`), and after an event only the
  sorted `affected` channels change, so the pass resumes at the first of
  them with the same additions. The resolution check (some positive rate
  below 2⁻⁵³ of the total) uses the smallest positive rate since the last
  exact pass as a lower bound and computes the exact minimum only when the
  bound fails. Scans of up to four rates run without a data-dependent
  branch (the target stays nonnegative, so zero rates never select).
- **Rates, counts and updates.** Channel count factors are an enum decided
  at lowering; first-order channels skip the jump table (callgrind showed
  5.5 M indirect mispredictions on W3 from it); expression rates of rules
  with one channel are evaluated without the cache lookup; diagnostics are
  built in cold functions; count updates apply in one pass and undo
  themselves on a range or instance error (tested: counts and the live
  count unchanged).
- **RSSA-CR.** W4 has 33 copies of one Michaelis–Menten law; rules whose
  rates compile to the same program now share one range cache (tested:
  shared bounds equal unshared ones). The interval test after an event is
  branch-free and the refresh runs out of line only when a species left
  its interval (W4: 7 % of events). Box corners are an enum instead of
  `&dyn Fn` closures.
- **Random words.** ChaCha8 (D4) costs about 3.7 ns per word against
  1.5 ns for xoshiro256++ (measured in a loop); doubling the words in situ
  slowed W1 and W5 by 4–16 %. No extra draws per event were found against
  0.6.0's fast path (linear: wait and selection; RSSA-CR: wait, group,
  joint member and acceptance). RSSA-CR now reuses the composition's
  accepting threshold — given the candidate, uniform below its upper bound
  — for its own test, one word less per trial (17–20 % of the words on W2,
  W5 and vilar-noise). Lemire's rejection threshold is a 64-bit division;
  it is now computed only when the low product is below the range (same
  accepted words; tested at the edges and on 100,000 random cases).
- **Allocations.** Selector tables (2,098 composition groups among them)
  are reused between the replicates of a worker thread; a nested run on the
  same thread gets new ones.

Not done: changing the linear scan's arithmetic (prefix sums with a binary
search would change selections, so `linear` would no longer equal
rules-run/0.6); a different RNG (D4); build-profile changes (LTO, one
codegen unit), which would also change the 0.6.0 comparison build.

## 3. Law gates and identity checks

Declared in `validation/merge/RESULTS.md` ("M4f ... declaration", commit
`c98dcfe`) before any run.

| Check | Result |
| --- | --- |
| RSSA-CR law gate, twelve models, gate seeds (`laws/m4f-rssa_cr.json`) | 4,814 comparisons, **0 failures**, max \|z\| 3.88; vilar-noise passes (chunk 1, seed 2,000,001, completes) |
| `linear`, `composition_rejection`, `sorting` unchanged: chunks 0 and 1 of every model (W5: 0–3) against the M4c chunks | 78 of 78 chunks identical (power sums, events, samplers) per selector run, 3 × 26 |
| Exact-law check of every selector (`laws/selectors-exact-m4f.json`) | 168 of 168 DKW tests, 16 of 16 power checks; byte-identical under `python -O` |

The unchanged check was declared as 72 chunks (two per model); W5's chunks
have 5,000 replicates, so `--replicates 20000` ran four of them (78 in
all). Every M4c chunk of these three selectors had completed under the
strict clock, so no stall occurred in them, and the other M4f changes are
byte-identical by construction; the identical chunks confirm it.

## 4. `auto`

M4c fixed `auto` = composition–rejection from 512 channels, else linear,
because its rule required RSSA-CR to win in every family, and RSSA-CR lost
narrowly with four molecules per species. M4f declared a benchmark with
population as a feature (eight families, 4–256 molecules per species) and
a rule family fitted by declared criteria to the twelve target models
(`selector_rule.py`):

- **Rule** (`select::auto`): RSSA-CR when no rule's rate reads two or more
  species and n ≥ 8 (n ≥ 128 when the mean initial count of the species that
  start with molecules is below 8; a network that starts empty counts as
  large); else composition–rejection from 128 channels; else linear.
- **Targets**: W1, W3 and the five examples pick `linear`, W2 and W4–W7
  `rssa_cr`; every pick is the fastest selector except michaelis-menten
  (linear, 1.08 × sorting). Largest ratio on the synthetic grid 1.58 (8
  channels with 32 molecules, where the linear scan is fastest), forced by
  the `michaelis` family, where RSSA-CR is 1.8 × faster at 8–16 channels.
- `describe()["contracts_1_0"]["selectors"]["auto"]` lists the thresholds;
  the run-1.0 and result-1.0 schemas carry the new description (regenerated
  with `crn-craft describe`; the `rules-*` schemas are unchanged).

Test written first: `auto_follows_the_declared_threshold` (rings at the
thresholds with 7 and 8 molecules, two-read rings, W1–W7 and the examples)
failed against the M4c rule.

## Checks

| Check | Result |
| --- | --- |
| Byte identity after every step: 156 CLI requests (`/home/claude/m07/m4f/harness.py`), the I1 stored requests, generated-loop replays | equal, except the two declared changes (RSSA-CR threshold reuse, `auto`) |
| Randomized expression test: 9,000 programs and the Hill and Michaelis–Menten shapes at 2,000 points against the tree evaluator | equal values |
| Clock rule tests (`crates/crn-craft/tests/clock.rs`, two generated-loop fixtures) | failed first, pass |
| Law gate of `rssa_cr` after the threshold reuse (`laws/m4f-rssa_cr.json`) | 4,814 comparisons, 0 failures (vilar-noise included, which failed in M4c) |
| Exact selector laws (`laws/selectors-exact-m4f.json`) | pass |
| `auto_follows_the_declared_threshold` | failed first against the M4c rule, passes |
| Speed gate (M4d) | passed on the rerun; see the [M4d log](M4d-log.md) |
| Rust on the integrated branch (`5c30f0a`): fmt, clippy, rustdoc, 532 release tests | pass |

## Limitations

- RSSA-CR results differ from M4c's (threshold reuse); run/1.0 is not
  released, so the sampler identity stays `rules-static-direct-rssa-cr/1`.
- The clock rule's limit is a count of consecutive stalls, not a measure of
  the accumulated error; a run with occasional stalls is accepted.
- `auto` decides from initial counts; a model whose populations grow or
  shrink far from them can get a slower selector. Its thresholds are fitted
  on one host.
- The CLI writer copies each nested object once per level.
- The linear scan still ends with a data-dependent branch for more than four
  channels (one misprediction per event on W3), and ChaCha8 costs about
  twice the words of the 0.6.0 fast path.

## Deviations from the plan and the task

- Task-level: the unchanged check ran 78 chunks instead of the declared 72
  (W5's chunk size), all identical.
- `next_reaction` keeps its strict clock (the decision names selectors).
- The selector timings for the `auto` fit (`speed-m4f-selectors.json`) ran
  after the benchmark, as declared, on the commit before the `auto` change
  (the selectors themselves did not change).

## Migration notes

- `crn-craft.run/1.0` Static Direct: a wait below half an ulp of the clock
  no longer fails a replicate; 1,000,000 consecutive ones fail with
  `CRN-R030` "clock failed to advance in 1000000 consecutive steps".
  `rules-run` requests are unchanged.
- `auto` (the run/1.0 default) now picks RSSA-CR for most networks of 8 or
  more channels, so run/1.0 results without a selector differ from M4c's;
  `"selector": "linear"` still reproduces rules-run/0.6. Generated loops
  and kernels refuse a request whose `auto` picks RSSA-CR (`CRN-R020` at
  `/method/selector`): name `"selector": "linear"` with `compiled=`.
- RSSA-CR results under run/1.0 differ from M4c's (one uniform less per
  trial).
- `describe()["contracts_1_0"]["selectors"]["auto"]` has a new shape:
  `rssa_cr_from_channels` {`population_below_cut`, `otherwise`},
  `population_cut`, `composition_rejection_from_channels`.
- The CLI depends on `serde` directly (already in the lockfile and the
  dependency policy).
