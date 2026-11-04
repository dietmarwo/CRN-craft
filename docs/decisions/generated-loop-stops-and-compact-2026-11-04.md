# Generated loops: emitted stop predicates and a compact layout

4 November 2026. Declared before any change was measured. It carries out
two of the next steps listed in the
[generated loop record](generated-rules-loop-2026-11-03.md).

## Starting point

The generated event loop emits only the guarded fraction stop
(`N > 0 && part >= f * N`). Every other stop predicate runs through the
host's bytecode after each event that changes a projection it reads. The
generator refuses loops above 1 MiB of source. Rate code is repeated in
every dependency row, so a network with one shared vacancy base and 128
instances (256 channels) is refused.

## S — emitted stop predicates

Emit every stop predicate from the compiled expression tree, with the
bytecode's semantics:

- Literals and parameters are compiled in, and projections are read from
  the values.
- `!`, `&&`, `||` and `if` return a NaN test or condition unchanged.
  `&&` and `||` short-circuit and return 0 or 1.
- Comparisons return 0 or 1. `min` and `max` return NaN for a NaN
  argument.
- `^`, `exp`, `log`, `sqrt`, `floor` and `ceil` call the same `libm`
  functions, which the crate re-exports for emitted code.

A NaN result defers to the host, which evaluates the bytecode and raises
the interpreter's diagnostic. Any other result decides the stop directly.
The fused fraction form stays as it is.

Gate: byte identity on the Rust replay matrix and the validation cases.
Keep only if the generated loop gains at least 5% on a probe request whose
stop is not in the fused form (`N > 0 && mutants / N >= 0.9`).

## C — compact layout

Add a table-driven layout:

- one rate function with a branch per channel;
- vacancy bases computed once per event;
- dependency rows as constant tables;
- count changes, changed projections and stop flags as constant tables;
- one function per projection;
- loops for the total, the resolution check and the scan.

A dependency row that covers more than half of the channels is replaced by
an update of every channel, in channel order. Channels outside the row have
unchanged inputs and were valid when last computed. Recomputing them
therefore gives the same values and cannot report an earlier error.

`generate_rules_loop` uses the unrolled layout up to 1 MiB and the compact
layout beyond it. Compact loops are refused above 8 MiB.
`generate_rules_loop_with` selects a layout explicitly, and
`execute_with_rules_loop` accepts a loop of either layout.

Gate: byte identity on the replay matrix, run with compact loops of
several fixtures, and on a dense fixture that only the compact layout can
serve. Adopt the layout if it is at least 1.25× faster than the
interpreter on the dense probe (128 instances sharing one vacancy base).
Record the cost of forcing the compact layout on the mtDNA probe for
information.

## Rules

- Interpreted samplers stay byte-identical to `main`, and no probe workload
  becomes more than 5% slower.
- Failing-first tests, the full Rust and Python checks, and every attempt,
  including failed ones, are recorded in
  `validation/mtdna/generated-loop-followups.json`.

## Outcome

Recorded 4 November 2026. S is adopted and C is not. Every run is in
`validation/mtdna/generated-loop-followups.json`. The final validation run
is in `validation/mtdna/generated-loop.json`.

**S.** On the probe request whose stop is `N > 0 && mutants / N >= 0.9`,
the generated loop took 90 ns per event with the host's bytecode and 53 ns
with the emitted predicate (1.69×). The interpreter took 144 ns. Under
default Static Direct the figures were 115 and 76 ns. Outputs were
byte-identical in every case:

- the probe;
- the Rust replay matrix, including a new fixture whose stop uses every
  operator and function a stop may use;
- an undefined stop, which still reaches the host's diagnostic;
- the validation cases.

An emitted `||` that returned 0 instead of 1 failed the replay test of the
general stop. In the final run the loop took 48 ns per event on that
request, against 138 ns for the interpreter.

**C.** The compact layout replayed byte for byte, but it missed the
criterion on the dense probe:

- *First version.* Full-row updates went through the rate function. It
  took 1621 ns per event against 1372 ns for the interpreter (0.85×).
- *Second version.* It wrote the update of every channel out once as
  straight-line code. It took 1218 ns against 1398 ns (1.15×).

Both are below the declared 1.25×. Each dense event must sum 256 rates, scan
them and recompute two projections of 128 terms, in the interpreter's
order. Those serial floating-point chains bound any byte-identical layout,
and the interpreter already comes close to that bound. The layout is left
out, and densely coupled networks above 1 MiB are still refused. The
attempt, with its tests and fixtures, is in
`validation/mtdna/followup-attempts/compact-layout.patch`.

The same runs exposed a limit of the unrolled layout. For networks of many
channels with constant rates and short dependency rows, the unrolled loop
is no faster than the interpreter and can be slower:

- on the 96-channel sparse network, 0.85× in one run and 1.03× in the
  final one, so no faster;
- on the 512-channel wide network, 950 against 630 ns per event (0.66×).

There the cost is dominated by work that both do in the same order: the
total, the resolution check and the scan over every channel. The generated
loop's gains come from models whose events are dominated by vacancy rates,
projections and the stop predicate, as in the mtDNA probe. The
documentation now says to benchmark before adopting the loop for
many-channel networks.
