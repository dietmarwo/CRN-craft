# Closing the generic rules gap: iterative plan

3 November 2026. Declared before any change was measured.

## Starting point

On a two-CPU host the 128-cell mtDNA probe (25,920 h, 720 h reports, 90%
takeover stop, two class trackers, one worker) takes:

- 108 ns per event on the generic opt-in path (Static Direct with Ziggurat
  waits and sparse trackers);
- 134 ns per event on default Static Direct;
- 35 ns per event on `mtdna_direct`.

The generic path is thus 3.1× slower than the specialized sampler. Its profile
attributes 33% to rate evaluation, 18% to selection with the per-event
resolution check, 16% to the event loop and 11% to observable refresh. The
harness is `validation/probe_structural_rules.py` (five workloads, median of
alternating runs against a reference binary, byte comparison of outputs).
Its results are recorded per iteration in
`validation/mtdna/structural-iterations.json`.

## Rules for every iteration

- **Keep identity.** Iterations T1–T4 must leave every existing sampler
  byte-identical to `main` (`b4b5f21`): default Static Direct, Ziggurat and
  sparse-tracker options, next reaction, thinning, tau, binomial,
  `mtdna_direct`, generated kernels and the tree. The gate is the 96-case
  rules replay plus the probe's five workloads.
- **Pay for complexity.** An identity-preserving iteration is kept only if it
  improves the generic opt-in path by at least 5% without slowing any other
  workload by more than 5%. Otherwise it is reverted, and the attempt is
  recorded.
- **Test first.** Each change starts with a failing test where behavior
  could change, and runs fmt, Clippy with warnings denied, release tests and
  the Python suites.

## Iterations

1. **T1 — rate evaluation.** Store each channel's rate factors inline in its
   plan, so evaluation no longer reaches through `net.rates`. Compute each
   needed vacancy base once per event instead of checking a sentinel per
   channel. Keep the arithmetic order: ((base × f₁) × f₂ …) × (weight ×
   channel factor).
2. **T2 — selection.** Fuse the total, the resolution floor (computed once,
   as `total × 2⁻⁵³`) and the scan into one pass, with the same summation
   order and the same error conditions.
3. **T3 — observables.** Update sum observables incrementally when every
   weight is a nonnegative dyadic rational (an integer multiple of 2⁻³⁰ at
   most) and the value stays below 2²³. Both the incremental update and the
   recomputation are then exact, so the bits match. Values outside the guard
   are recomputed as now.
4. **T4 — event loop.** Profile-guided removal of per-event work that does
   not change results: the report check, live-instance accounting, the stop
   and tracker dispatch, bounds checks.
5. **T5 — decision.** Measure the gap. If the generic opt-in path is within
   1.5× of `mtdna_direct`, stop and document. Otherwise do T6.
6. **T6 — opt-in grouped evaluation (conditional).** A new request version
   adds `evaluation: "grouped"` on `static_direct`. It precombines each
   channel's constant factors, groups channels that share a vacancy base, and
   selects a group and then a channel. Resolution is required within each
   group and between groups, as in the tree engine. It changes rounding and
   random-word use, so it gets its own sampler identity and gates:
   - the K = 20 CME gate with 32,768 paths, family α = 0.001 and the existing
     bounds;
   - a production distribution gate against Static Direct with 2048 paths per
     side and nine CDFs (N, wild type and mutants at 250, 500 and 1000 h),
     with bound 2·sqrt(ln(2·9/0.001)/(2·2048));
   - replay across worker counts.

   Target: within 1.5× of `mtdna_direct`. Added after T1–T4 and before T6
   was measured: T6 is merged only if it passes every gate and makes the
   generic opt-in path at least 25% faster, enough to pay for a new
   contract version. Otherwise it is recorded and left out.
7. **T7 — documentation and PR.** Update every affected document, write the
   implementation log, open a PR and merge it when every check is green.

## Targets

- The generic opt-in path is within 1.5× of `mtdna_direct` on the probe,
  through T1–T4 if possible and otherwise through T6.
- No workload gets more than 5% slower than on `main`.
- Every existing sampler stays byte-identical, and failed attempts are kept
  in the log.

## Outcome

Recorded 3 November 2026, after the loop. No change was kept, and the
contract, the samplers and their identities are unchanged. Every attempt,
including the null runs that bound the noise, is in
`validation/mtdna/structural-iterations.json`. The patches for T1, T3 and
the T6 prototype are in `validation/mtdna/structural-attempts/` and apply
to `b4b5f21`. The probe is `validation/probe_structural_rules.py`.

Nanoseconds per event on the generic opt-in path, main against candidate,
median of alternating runs on the same two-CPU host:

| Iteration | Change | Generic opt-in | Other workloads | Decision |
| --- | --- | ---: | --- | --- |
| H0 | none (main against main) | 105 → 108 (0.97×) | 0.96–1.00× | noise floor |
| T1 | inline factors, eager vacancy bases | 102 → 101 (1.01×) | sparse Static Direct 0.92× | reverted |
| T2a | total and smallest rate in one pass | 108 → 110 (0.97×) | sparse Static Direct 1.31× | reverted |
| T2b | as T2a, branchless minimum | 111 → 118 (0.94×) | sparse Static Direct 0.44× | reverted |
| T3 | exact incremental sum projections | 104 → 100 (1.04×) | sparse paths 0.91× | reverted |
| T4 | event-loop trimming | not attempted | — | no candidate |
| T6b | grouped evaluation, one uniform | 107 → 103 (1.04×) | default Static Direct 1.12× | not merged |
| final | branch head | 108 → 109 (0.99×) | 0.93–1.02× | — |

The identical-code runs (H0 and final) moved by up to 7% on the short
sparse workloads and by 3% on the mtDNA workloads. T2a's sparse gain is
real, at 96 channels, but the declared rule targets the generic opt-in
path, where it lost 3%. A separate plan for many-channel networks can take
it up. T4 found no single per-event item that could clear 5% on its own.
The loop's own time sits in inlined code without line attribution.

T6 grouped channels by rate base and selected a group and then a channel.
The first prototype drew a second uniform inside a group. The tuned one
continued with the same uniform and checked every positive channel rate,
base × term, against the total. Both took about 100 ns per event, as the
per-channel path does. The work moved from `RateCache::channel` to the
group update without shrinking: terms for the changed species, two group
sums with their smallest terms, the vacancy division and the scans. It
stayed far below the 25% criterion, so no rules-run/0.5 was added.

### Where the time goes

A profile of the generic opt-in path on the T3 build (18.4 million events)
puts 33% in rate evaluation, 19% in the event loop, 18% in selection with
the Ziggurat wait, about 13% in random words, 6% in state updates, 4% each
in projection refresh and the stop predicate. Before T3, projection
refresh took 11–13%. No single part is large enough to close a 3× gap.
The generated rules kernel (S5 in the fast generic rules campaign) replaced
rate evaluation with emitted code and gained 1.11× on default Static
Direct. T6 gained no more: 1.12× there and 1.04× on the opt-in path.

`mtdna_direct` costs about 35 ns per event because it fuses the whole
event into fixed-arity code. It certifies resolution once at preparation,
and it computes its eight hazards from three counts with one division. Its
takeover stop and class trackers read counts directly, so without an
observer it refreshes no projections between reports. The generic path
must maintain every hot projection, evaluate the stop bytecode and run the
dependency lists after each event.

### What would close the gap

Interpreted structure stays within about 10% of today's generic path.
Closing the gap needs generated code for the whole event, not only for the
rates. The next step is to extend the ahead-of-time rules generator to
emit the event loop for a closed model with constant and vacancy rates:
fixed-order rates and selection, count updates, exact projection updates,
the stop predicate and count trackers as straight-line code. Its gates
would be those of the current generator (byte identity with Static Direct
where it keeps the arithmetic order) or, where it groups, the T6 gates
above. Until then, use `mtdna_direct` for the fixed three-class topology.
Other models use generic Static Direct with the 0.4 opt-ins.

Follow-up: the generated event loop is planned, measured and adopted in
[generated-rules-loop-2026-11-03.md](generated-rules-loop-2026-11-03.md).
