# Rules kernel review checklist

Use the generic Static Direct path first. Measure the default cached path,
then the opt-in Ziggurat and sparse count trackers when the changed replay
identity is appropriate. Compare the generic engine with the candidate on
the same model, output contract, workload and hardware. Try the ahead-of-time
generator when its supported rate and projection families cover the model.

Use the following checklist when the remaining measured gap warrants a
hand-written kernel.

1. **Profile the generic path.** Record the profile and the bottleneck the
   candidate addresses beyond shared rates, dependency updates, deferred
   reports, bytecode, incremental counts and optional sampling policies.
2. **Recognize the lowered network at preparation.** Check channel shapes,
   rate kernels and observables. Model and rule names are user labels.
3. **Certify exactness for reachable states.** Bound each positive hazard
   against the total in the same state. Require at least
   16 × 2⁻⁵³ of that total for a 53-bit grouped selector. Account for
   arithmetic order, zero hazards, threshold rounding and representability.
4. **Expose fallback explicitly.** Return Static Direct and include
   `CRN-R051` with the reason for an unsupported specialization. Exercise
   parameter sweeps through `bind_parameters` and check sampler and warning.
5. **Version replay.** Give the kernel its own sampler identity. Verify
   completed outputs at one and several workers, excluding request identity
   fields when worker count is part of the request.
6. **Preserve the scheduler contract.** An infinite wait reaches the horizon.
   Fill reports strictly before the next event with the pre-event state.
   At an exact tie, the existing right-continuous contract includes the event.
   Check the stop after each committed event and hold the stopped state.
   Match budgets, cancellation and polling.
7. **Keep observers out of the numbers.** A Continue observer must leave
   serialized simulation results identical. An observer or retained log
   must use the same tracker integration policy as a run without either.
8. **Run law gates.** Use an independent chemical master equation where a
   small fixture exists. Compare production distributions with Static
   Direct under a predeclared bound and keep any failed rows.
9. **Add Rust regressions.** Cover reports exactly at an event, fallback,
   threads, Continue observer identity and overflow.
10. **Request independent review before merging.** Include the recognition
    proof, numerical gate and timing evidence in the review.

## Ahead-of-time generated rules kernels

A Rust caller can generate a kernel from a prepared Static Direct request.
The generator covers constant and recognized vacancy rates and weighted sum
projections. It emits shared vacancy expressions, per-channel arithmetic,
static dependency sets, projection sums and a fingerprint: the SHA-256 of the
emitted source, keyed by the model bytes. Floating-point products and sums keep the
interpreter's order. The Static Direct scheduler handles state changes,
selection, waits, budgets, reports and observation.

```sh
cargo run --release --locked -p crn-craft-rules --example generate_rules_kernel -- MODEL REQUEST > kernel.rs
```

Compile the resulting module into the Rust application and use
`execute_with_rules_kernel(&prepared, &control, &GeneratedRulesKernel)`.
`execute_with_rules_kernel` regenerates the source and refuses a kernel
whose fingerprint differs. That covers other model bytes, including parameter
rebinding, and kernels compiled from an older generator. A request may select inverse or Ziggurat waits using
the same compiled kernel. Build and benchmark each generated kernel against
the interpreter under the same sampler. Python use would require a wheel
built with the generated module and a binding to dispatch it.

`crates/crn-craft-rules/tests/generated_kernel.rs` runs a checked-in kernel
for the three-class fixture under `cargo test`. It checks the generator's
output against that source, byte-identical replay with both waiting-time
policies and rejection of a kernel for other model bytes. After a generator
change, regenerate the fixture with the example above. Run
`validation/check_generated_rules.py` for a fresh build and timed
comparisons.

On a two-CPU host, profiling the generic path with Ziggurat waits and sparse
trackers (128 cells, one worker, 18.4 million events) attributed 33% of the
time to rate evaluation, 18% to selection with the per-event resolution check
and 11% to observable refresh. That path took 2.07 s, against 0.76 s for
`mtdna_direct`, whose resolution is certified once at preparation and whose
selection is grouped.

A later loop of structural changes to the interpreter tried inline rate
factors, a fused total and resolution check, exact incremental sum
projections and a grouped evaluator that selects a group and then a
channel. None gained 5% on that path, and the grouped evaluator gained
1.04× there and 1.12× on default Static Direct. The cost is spread over
rate evaluation, selection, random words, projection maintenance and the
stop predicate. `mtdna_direct` fuses all of them and refreshes no
projections between reports. Closing the gap therefore needs generated
code for the whole event, not only for the rates. The
[decision record](decisions/generic-rules-structural-2026-11-03.md) has
the measurements, the profile and that next step.
