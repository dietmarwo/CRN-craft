# Generated rules event loop: iterative plan

3 November 2026. Declared before any change was measured. It follows the
[structural iterations](generic-rules-structural-2026-11-03.md), which
showed that interpreted changes stay within about 10% of the generic path
and proposed generating code for the whole event.

## Starting point

On the two-CPU host, `validation/probe_structural_rules.py` measured about
108 ns per event for Static Direct with Ziggurat waits and sparse trackers
on the 128-cell three-class mtDNA probe, against about 36 ns for
`mtdna_direct` (3.0×). The existing generated rules kernel replaces rate
evaluation only and gained 1.11× on default Static Direct.

## Design

A new opt-in Rust route, next to the existing rates-only kernel:

- `generate_rules_loop(&prepared)` emits a module that defines
  `GeneratedRulesLoop`, which implements a new `RulesLoop` trait.
  `execute_with_rules_loop(&prepared, &control, &kernel)` runs it. Like the
  rates kernel, its fingerprint is the SHA-256 of the emitted source keyed
  by the model bytes, and a stale or foreign kernel is refused with
  `CRN-R020`.
- The emitted code owns the event: rates in the interpreter's arithmetic
  order, the dependency updates after each fired channel, the fixed-order
  total, the per-event resolution check, the scan, the state update with
  its range and instance checks, and the projection refresh. Everything is
  straight-line code over fixed-size arrays.
- A host type in the crate supplies everything that depends on the request
  or on run control: waiting times and uniforms, polling, reports and the
  horizon, budgets and counters, tracker observation and the stop
  predicate. One kernel therefore serves every Static Direct request for
  its model bytes, with inverse or Ziggurat waits and either tracker
  policy.
- Results must be byte-identical to interpreted Static Direct for the same
  request. The emitted code performs the same floating-point operations in
  the same order. It refreshes every projection that the fired channel
  changes, a superset of the interpreter's hot set. The extra values are
  pure functions of the counts, cold projections cannot overflow by
  construction, and they are recomputed before every report, so no output
  changes.

Observers stay unsupported on this path, as for the rates kernel.

## Rules for every iteration

- **Keep identity.** Interpreted samplers are not changed apart from the
  plumbing that selects the new route; the probe's five workloads must stay
  byte-identical to `main`. The generated loop must reproduce interpreted
  Static Direct byte for byte on the probe requests and on a replay matrix
  in the Rust tests.
- **Pay for complexity.** An iteration after the first is kept only if it
  makes the generated loop at least 5% faster on the generic opt-in
  request. Failed attempts are recorded in
  `validation/mtdna/generated-loop-iterations.json`.
- **Test first.** Each behavior starts with a failing test. Every iteration
  runs fmt, Clippy with warnings denied, release tests and the Python
  suites.

## Iterations

1. **L0 — baseline.** Build a runner crate for the probe requests and
   measure the interpreter, the rates kernel and `mtdna_direct` on the same
   host.
2. **L1 — generated loop.** Add the host type, the trait, the generator and
   the execution route. Gates: byte identity on the probe requests (generic
   opt-in and default) and a Rust replay matrix covering inverse and
   Ziggurat waits, both tracker policies, a model stop, reports at event
   times, the event, sampling-step and instance limits, a horizon before the
   first event, extinction and an invalid rate.
3. **L2 — stop predicate.** Emit the stop predicate as Rust when it uses
   only projections, parameters, constants, arithmetic, comparisons and
   logical operators, with the interpreter's NaN diagnostic. Other
   predicates stay with the host.
4. **L3 — host trimming.** Profile-guided reduction of per-event host
   calls without changing results.
5. **L4 — decision.** Adopt the route if every identity gate passes and the
   generated loop is at least 1.5× faster than the interpreter on the
   generic opt-in request. Otherwise record it and leave it out.
6. **L5 — documentation and PR.** Check in a generated fixture, update every
   affected document, write the implementation log, open a PR and merge it
   when every check is green.

## Targets

- Within 1.5× of `mtdna_direct` on the generic opt-in probe request.
- Byte-identical results to interpreted Static Direct.
- No interpreted workload more than 5% slower than on `main`.

## Outcome

Recorded 3 November 2026, after the loop. The generated event loop is
adopted. Every attempt is in
`validation/mtdna/generated-loop-iterations.json`, and the final validation
run is in `validation/mtdna/generated-loop.json`.

Nanoseconds per event on the 128-cell probe, two-CPU host, medians of
alternating runs:

| Iteration | Change | Generic opt-in | Default | Decision |
| --- | --- | ---: | ---: | --- |
| L0 | interpreter | 107 | 135 | baseline |
| L0 | `mtdna_direct` | 36 | — | reference |
| L1 | generated loop | 49 | 80 | kept |
| L2 | emitted guarded fraction stop | 46 | 75 | kept (7% over L1) |
| L3 | constant rule index for counters | 0.98× of L2 | 1.06× | reverted |
| L4 | validation run on the branch head | 49 (2.31×) | 72 (1.85×) | adopted |

In the L4 run, after the review fixes, the interpreter took 114 and 134 ns,
the rates-only kernel 99 and 127 ns and `mtdna_direct` 36 ns, so the
generated loop is within 1.36× of the specialized sampler. The interpreted
samplers stayed byte-identical to `main` on the probe's five workloads,
and none was more than 3% slower. Both criteria hold: at least 1.5× faster
than the interpreter, and within 1.5× of `mtdna_direct`.

The identity gates held on every case:

- **Probe.** The generic opt-in and default probe requests are
  byte-identical.
- **Rust replay matrix.** The nine checked-in fixtures replay every Static
  Direct policy, and the event, sampling-step and instance limits. They
  also cover a horizon before the first event, a model stop, a general and
  an undefined stop predicate, extinction, loss of resolution, an invalid
  rate, a count range error, projection overflow, shared vacancy rates,
  event logs and an empty network.
- **Validation script.** Eleven generated replay cases, with one and eight
  workers, include an event log, an overflowing projection and an
  undefined stop.

A deliberately wrong generator, one that summed the total in reverse order,
failed six of the first ten Rust tests.

An independent review found no difference from the interpreter in code
reading or in its own replays of six further edge models, which were then
adopted as fixtures. It did find three problems, all fixed before merging:

- Rate code repeated per dependency row made densely coupled networks
  emit source that exhausts the compiler. Loops above 1 MiB of source are
  now refused.
- An empty network left an unused label, which fails builds that deny
  warnings.
- The cancellation test returned before entering the loop. A new test
  cancels a running replicate and checks that it ends at a polling point.
  It fails when the emitted poll is removed.

The review also showed that a projection limit of 1024 terms could not
trigger (lowering allows 256 species), so it was dropped. It also made
unsupported models report their own reason instead of a fingerprint
mismatch.

The remaining profile is the division in the vacancy base and the wait, the
rate validity checks that wait on it, random words, tracker observation
and the per-event budget calls. `mtdna_direct` still saves the per-event
resolution check and projection refresh by certifying resolution once and
reading counts directly. Those savings change replay identity, so they are
out of scope for a byte-identical loop.

General stop predicates, other than the guarded fraction, are still
evaluated by the host's bytecode. Expression rates and max, min or distinct
projections are refused. The route is Rust only.
