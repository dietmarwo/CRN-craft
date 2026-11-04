# Generated-loop follow-ups: implementation log

4 November 2026, from `main` at `a78a3a5`. Three next steps from the
[generated loop record](../decisions/generated-rules-loop-2026-11-03.md):

- emit every stop predicate;
- add a compact layout for networks above the 1 MiB limit;
- plan and run the many-channel selection change found in the structural
  loop.

The plans are
[generated-loop-stops-and-compact-2026-11-04.md](../decisions/generated-loop-stops-and-compact-2026-11-04.md)
and
[many-channel-selection-2026-11-04.md](../decisions/many-channel-selection-2026-11-04.md).
Both were committed before any measurement. Every run is in
`validation/mtdna/generated-loop-followups.json`.

## Many-channel selection (interpreter)

`direct_with` now calls `fused_total` for networks of at least 32 channels.
It sums the rates from −0.0 in channel order, which is the order
`Iterator::sum` uses, and records the smallest positive rate. The
resolution test then compares that minimum with `total × 2⁻⁵³`. Smaller
networks keep the original sum and `any` test.

The probe gained a 512-channel workload of 256 attributed instances that
replicate and die. Plain species could not be used: the model limit is 64
types.

Against `main`, sparse Static Direct ran 1.30× faster and the wide workload
1.39× faster. Every probe output and all 96 rules-replay cases were
byte-identical.

## Emitted stop predicates

`Expr::rust_predicate` writes the stop as Rust from the expression tree,
following `Predicate::run`:

- Literals and parameters become bit patterns, and observables are read
  from `values`.
- `!`, `&&`, `||` and `if` pass a NaN test or condition through.
  Short-circuit results are 0 or 1.
- Comparisons produce 0 or 1, and `min` and `max` give NaN for a NaN
  argument.
- `^`, `exp`, `log`, `sqrt`, `floor` and `ceil` call `libm`, which the
  crate now re-exports as a hidden item so that emitted code needs no
  dependency of its own.

The loop treats a NaN value as undefined and calls `LoopHost::stop`, which
repeats the interpreter's sequence and raises its diagnostic. Any other
value decides the stop on the spot. The fused fraction form is unchanged.

A new fixture's stop uses every operator and function a stop may use. When
the generator emitted `||` as returning 0, the general-stop replay failed;
the correct generator passes it. With the stop
`N > 0 && mutants / N >= 0.9`, the probe's generated loop took 53 ns per
event instead of 90 ns.

## Compact layout (not adopted)

The compact layout emitted:

- one rate function with a branch per channel;
- vacancy bases computed once per event;
- the dependency rows, count changes, changed projections and stop flags as
  constant tables;
- one function per projection;
- loops for selection.

A dependency row covering more than half of the channels updated every
channel in channel order. That is identical, because the other channels
have unchanged inputs and were valid when last computed. Automatic
selection used the compact layout for loops above 1 MiB, and an explicit
`LoopLayout` let callers choose. `execute_with_rules_loop` accepted either
layout.

It replayed byte for byte on four fixtures in compact form and on a dense
256-channel fixture. When the full-row update was dropped, the replays
failed. On the dense probe the first version reached 0.85× of the
interpreter. The second, which wrote the full update out as straight-line
code, reached 1.15×, below the declared 1.25×. The layout, its tests and
the extended validation script are kept as
`validation/mtdna/followup-attempts/compact-layout.patch`. The patch also
contains the stop emission and applies to `a78a3a5`.

## Fused selection in generated loops (not adopted)

In unrolled loops, tracking the minimum inline made the sparse loop 1.39×
slower and the wide loop 1.86× slower, so unrolled loops keep the per-rate
test. The compact layout benefited, but it was not adopted.

The same measurements showed that the unrolled loop is no faster than the
interpreter on many-channel constant-rate networks:

- 0.85–1.03× on 96 channels, across runs;
- 0.66× on 512 channels.

The kernel checklist and architecture overview now say to benchmark first
in that case.

## Verification

- Ten checked-in loops and 17 integration tests, plus the fingerprint unit
  test, pass.
- The final `validation/check_generated_loop.py` run replayed every
  supported case byte for byte with one and eight workers, under default
  Static Direct and the Ziggurat and sparse-tracker options. It also
  confirmed that the aggregation model and the dense network are refused.
- Against `main`, the probe's six workloads were byte-identical and the
  96-case rules replay was identical.
- The Rust and Python checks passed, the Python suites normally and with
  `-O` with identical reports.
