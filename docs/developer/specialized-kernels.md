# Rules kernel review checklist

Use the generic Static Direct path first. Measure its default, which since
rules-run/0.6 uses Ziggurat waits and sparse count trackers, and inverse
waits with per-event trackers where the earlier replay identity is needed.
Under `crn-craft.run/1.0` also measure the selectors (`linear`,
`composition_rejection`, `sorting`, `rssa_cr`) and the `auto` default
([exact methods](direct-ssa.md#selectors)): on large networks a selector
often closes the gap that a kernel would. Compare the generic engine with
the candidate on the same model, output contract, workload and hardware. Try
the ahead-of-time generated event loop when its supported rate and
projection families cover the model.

The timings on this page were measured on releases 0.4.0 to 0.6.0, before
release 0.7.0 tuned the interpreted Static Direct loop (compiled rate
programs, running totals resumed at the first changed channel). The
interpreted rows are therefore slower than in 0.7.0; on the three-class
workload of `validation/merge/mtdna-workload/summary.json`, one thread,
0.7.0's Static Direct took 0.72 times the time of 0.6.0's. Measure again
before relying on a ratio below.

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

The generator examples (`generate_rules_kernel`, `generate_rules_loop`,
`benchmark_rules`) live in the deprecated facade crate `crn-craft-rules`;
the code they emit needs only `crn-craft`.

Compile the resulting module into the Rust application and use
`execute_with_rules_kernel(&prepared, &control, &GeneratedRulesKernel)`.
`execute_with_rules_kernel` regenerates the source and refuses a kernel
whose fingerprint differs. That covers other model bytes, including parameter
rebinding, and kernels compiled from an older generator. A request may select inverse or Ziggurat waits using
the same compiled kernel. Build and benchmark each generated kernel against
the interpreter under the same sampler. Python use would require a wheel
built with the generated module and a binding to dispatch it.

`crates/crn-craft/tests/generated_kernel.rs` runs a checked-in kernel
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
1.04× there and 1.12× with inverse waits and per-event trackers (the
Static Direct default before rules-run/0.6). The cost is spread over
rate evaluation, selection, random words, projection maintenance and the
stop predicate. `mtdna_direct` fuses all of them and refreshes no
projections between reports. Closing the gap therefore needs generated
code for the whole event, not only for the rates. The
[decision record](decisions/generic-rules-structural-2026-11-03.md) has
the measurements, the profile and that next step. The generated event loop
below implements it.

## Generated rules event loops

The loop generator emits the whole Static Direct event for every model
that Static Direct runs: constant, vacancy and expression rates, whole-rule
`propensity` expressions (the string rates of `crn_craft.Model`),
higher-order rules (through `crn_craft::rules::combination_rate`, the
interpreter's combination count), and `sum`, `distinct`, `max` and `min`
projections. The loop runs the linear scan: a request whose `selector`, or
the selector that `auto` picks, is another one is refused with `CRN-R020`
at `/method/selector`, so name `"selector": "linear"` in a `run/1.0`
request for a loop. The emitted
`GeneratedRulesLoop` holds the rates in a fixed-size array and, after each
event, updates the channels that the fired channel affects. It sums the
total in channel order, checks resolution, scans for the next channel,
validates and applies the count changes and refreshes the projections the
channel changes. It also evaluates the model's stop predicate. A guarded
fraction such as `N > 0 && mutants >= 0.9 * N` keeps the interpreter's
fused form. Any other predicate is emitted from its expression with the
bytecode's semantics: NaN propagation, short-circuit `&&` and `||`, and
the same `libm` functions. An undefined value is passed back to the host
for the interpreter's diagnostic. Rate and propensity expressions are
emitted the same way, with parameters compiled in as exact bit patterns and
observables read from the current values. The global factor is multiplied
by the reactant weight times the lowered channel factor, as in the
interpreter; for a whole-rule propensity the weight is 1 while every
pattern can take a distinct molecule and 0 otherwise, and the expression is
not evaluated at weight 0. A rate expression shared by several channels of
one rule is evaluated once per dependency row, as the interpreter's
per-event cache does. Projections use the interpreter's reductions over
the same terms in the same order. All of it is straight-line code with the
interpreter's operation order. `LoopHost` supplies the rest at run time:
polling, inverse or Ziggurat waits, uniforms, reports and the horizon, the
event, instance and sampling budgets, tracker observation and event logs.

```sh
cargo run --release --locked -p crn-craft-rules --example generate_rules_loop -- MODEL REQUEST > loop.rs
```

Compile the module into the Rust application and run it with
`execute_with_rules_loop(&prepared, &control, &GeneratedRulesLoop)`.
Results are byte-identical to `execute` for the same request. Request
settings are not compiled in, so one loop serves every Static Direct
request for its model bytes: either waiting time, either tracker policy,
any report grid, budget or worker count. As for the rates kernel, the
fingerprint is the SHA-256 of the emitted source keyed by the model bytes,
and a loop for other model bytes or from another generator version is
refused with `CRN-R020`. The generator refuses, with `CRN-R020` at
`/method`, exactly two kinds of request. One selects a method other than
`static_direct`; models with ATP, continuous variables or instances that
change identity need another method, so they never reach the generator.
The other is a loop whose source would exceed 1 MiB: the rate code is repeated in every dependency
row, so densely coupled networks grow with the square of their channel
count. A 64-channel network with one shared vacancy base emits about
0.8 MB, and 128 channels would exhaust a 7 GB build host. Such networks use
the rates kernel or the interpreter. A table-driven compact layout was tried
for them, but it gained only 1.15× over the interpreter on a 256-channel
dense network
([record](decisions/generated-loop-stops-and-compact-2026-11-04.md)).
Observers are not supported on this path; trackers and event logs are.

The loop refreshes every projection the fired channel changes. The
interpreter, without an observer, refreshes only the projections that rates,
the stop predicate or trackers read, and recomputes the rest before each
report. The extra values depend only on the counts and are recomputed
before every report, so outputs do not change. Projections that could
overflow are refreshed by both. A deferred projection cannot overflow: the
interpreter defers a projection only when each weight is at most
f64::MAX / (max(1024, terms) · 2⁵³), with 2⁵³ the largest count, so a sum
over its terms stays finite.

`crates/crn-craft/tests/generated_loop.rs` runs checked-in loops for
nineteen fixtures under `cargo test` (`crates/crn-craft/tests/generated/`):

- the three-class model with its takeover stop;
- identical and distinct pairs with extinction and a general stop
  predicate;
- the same pairs with a stop that uses every operator and function a stop
  may use;
- a channel below 53-bit resolution, with an influx;
- an invalid vacancy rate;
- vacancy rates shared by several channels, with a fraction stop or with
  a stop that becomes undefined;
- a count that leaves its range;
- a projection that overflows;
- a network without channels;
- the SIR model as `rules.Model` writes it, with the string rates
  `"beta * S * I / N"` and `"gamma * I"`, and the same model with a stop;
- Michaelis–Menten conversion `V * A / (Km + A)` beside constant rates;
- the three-class model with degradation rates `delta * (1 + mutants / K)`
  and an `if` in a replication rate;
- rates and propensities with every operator and function a rate may use,
  `distinct`, `max` and `min` observables and a whole-rule identical pair;
- expression rates that become negative or undefined;
- higher-order rules of a `crn-craft.model/1.0` model;
- a wait too small to move the clock, once and a million times in a row,
  under `run/1.0` and the `rules-run` versions.

The tests replay every Static Direct policy and the event, sampling-step
and instance limits. They also cover trackers, first-passage trackers,
event logs, reports close together, a horizon before the first event, the
stops and the errors. They check that cancellation is polled inside the
loop, and that other model bytes, other methods and oversized networks are
refused. `choice` and the other random functions never reach a rate: a rate
that uses one is invalid (`CRN-R012`), and a model that draws product
attributes is open, so Static Direct refuses it before a loop exists. The
string-rate fixtures run every Static Direct policy with one and two
threads, trackers, first-passage trackers, an event log and the model's
stop, and the event, sampling and instance limits.
`validation/check_generated_loop.py` generates loops for thirteen replay
cases, among them a stop outside the fused form and a closed aggregation
model with a max projection. It replays them against the CLI with one and
eight workers. It checks that an oversized dense network is refused, and
it times the probes.

On the two-CPU host, the 128-cell probe took about 48 ns per event with the
generated loop under Ziggurat waits and sparse trackers. That is 2.4× faster
than the interpreter (112 ns) and 1.25–1.34× the cost of `mtdna_direct`
(36–38 ns across runs). With inverse waits and per-event trackers (the
default before rules-run/0.6) the loop took 69 ns, against 135 ns. The rates-only kernel took 94 ns and 122 ns. With the stop
written as `N > 0 && mutants / N >= 0.9`, the loop took 48 ns against
138 ns.

Benchmark before adopting the loop for networks with many channels,
constant rates and short dependency rows. There, every event's work is the
total, the resolution check and the scan over all channels, and both paths
do that work in the same order. On the 96-channel sparse network the
unrolled loop was no faster than the interpreter (0.85–1.03× across runs).
On a 512-channel network it was clearly slower: 950 against 630 ns per
event. See the
[decision record](decisions/generated-rules-loop-2026-11-03.md),
[the follow-ups](decisions/generated-loop-stops-and-compact-2026-11-04.md)
and `validation/mtdna/generated-loop.json`.

### From Python: `compile_loop`

Opt-in, off by default, for users with a Rust toolchain at least at the
MSRV (1.85). `crn_craft.compile_loop(model)` (also
`crn_craft.rules.compile_loop`) takes the loop source from the extension
(`loop_source(model)`, which needs no cargo and raises the generator's
`CRN-R020` for unsupported models) and writes a small Cargo project: the
generated loop, the exact model bytes, a driver
(`bindings/python/src/rules_loop_main.rs`) and the result encoder that the
extension compiles itself (`bindings/python/src/rules_payload.rs`). The
project depends on `crn-craft` at exactly the installed version, from
crates.io or from a checkout (`source=` or `CRN_CRAFT_SOURCE`, a tree with
`crates/crn-craft/Cargo.toml`), and builds a standalone program. Generated
code names `crn_craft::rules::` paths, so a Rust project that compiles it
depends on `crn-craft` directly; code generated by 0.6.0 named
`crn_craft_rules::` and is refused by its fingerprint until it is
regenerated. `simulate(compiled=loop)` sends the request to it on
standard input and reads the metadata and record columns in the
extension's own form. No Rust value crosses a library boundary, so the
program's toolchain need not match the wheel's.

Builds are cached in the user cache directory, keyed by the loop
fingerprint, the installed version, `rustc -vV`, the target and Rust
flags, the driver and, for a checkout, its engine sources; one Cargo target
directory serves every model. Loops cached by release 0.6.0 are rebuilt
once, into new cache entries, with the same results. Every `compile_loop` call runs short requests
through the program and the interpreter and returns the loop only if the
result bytes are identical; `CompiledLoop.verify(run)` checks a user's
request the same way. Methods other than `static_direct`, other model
bytes, bound candidates and callbacks are refused.

On the 1000-cell three-class workload with two threads, end to end from
Python, the compiled loop took 2.6 s against 5.8 s for interpreted Static
Direct with its rules-run/0.6 defaults, Ziggurat waits and sparse trackers
(2.21×), and 3.8 s against 7.0 s with inverse waits and per-event trackers
(1.86×; declared target 1.5×). The first build took 42 s on the two-CPU
host, a second model 3 s and a cached call 0.06 s
([record](decisions/compiled-rules-loop-python-2026-11-06.md),
`validation/mtdna/compiled-loop-python.json`).
For string rates the declared target was the same 1.5×, on the SIR model
as `rules.Model` writes it, scaled to N = 10,000 (2,000 replicates to time
400, two threads, 37.6 million events). The compiled loop took 1.64 s
against 3.35 s with inverse waits (2.00×) and 1.36 s against 2.97 s with
Ziggurat waits (2.18×), medians of three alternating runs end to end from
Python (`validation/benchmark_compiled_string_rates.py`,
`validation/rules/compiled-loop-string-rates.json`).
`validation/check_compiled_loop.py` checks identity over threads, waiting
times, tracker policies, event logs, stops, budgets and cancellation, on
the three-class model and on four string-rate models (the builder SIR,
Michaelis–Menten, expression rates that read observables, and `max`, `min`
and `if` over `distinct`, `max` and `min` observables), and the refusals;
it runs in CI against the checkout. See the
[user manual](mtdna-user-manual.md#optional-compiled-event-loop) for usage.
