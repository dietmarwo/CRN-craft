# Generated rules event loop: implementation log

3 November 2026, from `main` at `cae67c1`. The plan, its rules and the
outcome are in the
[decision record](../decisions/generated-rules-loop-2026-11-03.md). This
log describes the design as built and how it was checked.

## Design as built

**Public surface.** `crn-craft-rules` gains:

- `RulesLoop`: an object-safe trait with `fingerprint` and `run`;
- `LoopHost` and the opaque `LoopHalt`;
- `generate_rules_loop` and `execute_with_rules_loop`;
- the example `generate_rules_loop`.

The rates-only `RulesKernel` route is unchanged.

**Dispatch.** `simulate_flat` hands the replicate to the loop after the
same preparation as the interpreter: the initial refresh, the monitor start
and the initial stop check. The replicate's `Words` holds a reference to the
ChaCha8 stream. A run with a loop rebuilds `Words` for the concrete stream
type through an `Any` downcast, so `Path::run` and `simulate_flat` now
require `R: 'static`; production always passes `ChaCha8Rng`. The word count
is copied back after the run, so `rng_words` matches the interpreter.

**Emitted code.** The arrays have fixed sizes: counts, values and rates.

- *Rates.* A `match` on the previous channel updates the rates that channel
  affects. The first event updates all of them. Each rate uses the
  interpreter's expression, `global × (weight × factor)`, with the vacancy
  base as `(limit − count).max(0) / denominator` and the same weight-zero
  shortcut. Invalid rates are reported in the interpreter's order.
- *Selection.* The total is summed from −0.0 in channel order. The
  resolution test compares every positive rate with `total × 2⁻⁵³`, and the
  scan falls back to the last positive channel.
- *State update.* Each channel's arm checks every count change against
  `0..=2⁵³`, checks the live-instance count with the host and then writes
  the counts. It refreshes the projections the channel changes, in order,
  with the overflow check. When the stop predicate reads a projection that
  changed, the arm evaluates it: a guarded fraction is emitted, and the
  host evaluates any other predicate.

**Host.** `LoopHost` mirrors each interpreter step:

- polling every 1024 events;
- the waiting-time policy and its clock check;
- report saving with the cold refresh, and the horizon;
- `check_events` and `check_instances`;
- the counters and the time update;
- `observe_flat` with the monitor;
- the stop sequence, including the diagnostic path that refreshes cold
  projections before the error.

## Iterations

L1 produced the first loop. On the generic opt-in probe request it took
49 ns per event against 107 ns for the interpreter, with byte-identical
output. The profile then showed 11% in the stop bytecode. L2 emitted the
guarded fraction and reached 46 ns, a 7% gain, so it was kept. In L3 the
rule index of each channel became a constant for the event counters. The
generic request ran at 0.98× of L2 and the default one at 1.06×, so L3 was
reverted.

## Review

An independent review compared the emitted code and the host with the
interpreter step by step. It covered the order of every step, the error
paths, the refresh of deferred projections and the floating-point
expressions, and found no difference. Its own replays of six edge models
were byte-identical. Five of those models became fixtures. The review also
found the following, and each was addressed:

- **Quadratic source.** Rate code is repeated in every dependency row, so a
  network whose channels share one vacancy base grows with the square of
  the channel count. At 128 channels the emitted 3.1 MB exhausted a 7 GB
  host's compiler. The generator now refuses loops above 1 MiB.
- **Unused label.** A network without channels left `'pick` unused, which
  breaks builds with warnings denied. `unused_labels` is now allowed in
  the emitted function.
- **Cancellation test.** The original test cancelled before the run, so
  the loop was never entered. A new test cancels while replicates run and
  checks that the interrupted replicate stopped at a multiple of 1024
  events. With the emitted `host.tick()` removed, it fails.
- **Further gaps.** Missing coverage of event logs, the count range,
  projection overflow, an undefined stop, shared vacancy rates, empty
  networks and fingerprint coverage was added.
- **Dead limit.** The 1024-term projection limit could not trigger,
  because lowering allows 256 species. It was removed.
- **Refusal reasons.** `execute_with_rules_loop` now reports why a model is
  unsupported instead of a fingerprint mismatch.

The review noted that the fingerprint guards against stale code, not
against a hand-written `RulesLoop`. It also noted that a change to a host
call that leaves the emitted text alone needs a manual bump of the
`rules-loop/1` tag. The code comment says so.

## Verification

- `tests/generated_loop.rs` checks the nine checked-in loops against the
  generator. It replays them against `execute` under all four Static
  Direct policies: with and without trackers, first-passage trackers and
  event logs, with one to three workers and with tight budgets. It asserts
  that each of these occurs: the event-limit (`CRN-R040`), sampling-limit
  (`CRN-R041`) and instance-limit (`CRN-R042`) outcomes, `CRN-R031`, the
  invalid-rate, count-range, overflow and undefined-stop failures, a stop
  and a completed path. It also checks that cancellation is polled inside
  the loop, and that other model bytes (`CRN-R020`), other methods,
  expression rates and oversized networks are refused. A unit test checks
  that the fingerprint changes with the emitted arithmetic.
- With the total summed in reverse, six of the first ten tests failed. The
  correct generator passes them all.
- `validation/check_generated_loop.py` built loops for eleven cases and
  replayed them against the CLI byte for byte with one and eight workers,
  under default Static Direct and the Ziggurat and sparse-tracker options.
  The cases include an event log, a projection that overflows, an undefined
  stop and a resolution failure. It confirmed that the max-projection
  aggregation model is refused and that a loop for other parameters is
  rejected, then timed the probe.
EOF
echo ok