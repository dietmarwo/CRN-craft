# Generic rules structural iterations: implementation log

3 November 2026, from `main` at `b4b5f21`. The plan, its rules and the
outcome are in the
[decision record](../decisions/generic-rules-structural-2026-11-03.md). This
log describes how each iteration was run and what it showed.

## Procedure

`validation/probe_structural_rules.py` times five workloads on a reference
and a candidate CLI binary, alternating the order between repeats. It
records every timing, the median, the event counts of both sides,
nanoseconds per event and whether the outputs are byte-identical. The
workloads are the generic opt-in path (Static Direct with Ziggurat waits
and sparse trackers), default Static Direct and `mtdna_direct`, all on the
three-class mtDNA model with 128 cells, plus default Static Direct and next
reaction on a 48-species birth–death network with 96 channels. The
reference was a copy of the `main` release binary. Profiles used gprofng on
a release build with line tables.

Each attempt was built, probed and then kept or reverted by the declared
rule. A null run of identical code (H0) came first, and another (final)
came last. Together they bound the noise on this host: up to 3% on the
mtDNA workloads and up to 7% on the sparse ones, which take under 0.3 s.
After T2, the rule's decisions were based on nine-repeat runs.

## Iterations

**T1 — rate evaluation.** Each plan carried its vacancy factors inline,
and `RateCache::begin` computed the needed vacancy bases once per event, in
the same arithmetic order. The outputs were identical. The generic path
gained 1% and sparse Static Direct lost 8%, so T1 was reverted.

**T2 — selection.** `direct_with` summed the rates and tracked the
smallest positive rate in one pass. The resolution check then compares
that minimum with `total × 2⁻⁵³`, which is equivalent to checking every
rate. The error order and the random-word use did not change. Sparse
Static Direct gained 1.31×, because the separate check pass over 96 rates
went away. The generic path, with eight channels, lost 3%. A branchless
minimum made sparse Static Direct 2.3× slower, through the loop-carried
dependency. Both variants were reverted.

**T3 — projections.** Sum projections with nonnegative weights that are
integer multiples of 2⁻³⁰ were updated by an exact per-channel delta while
the stored and updated values stayed in [0, 2²³). In that range every
partial sum is exact, so the update has the same bits as a recomputation.
Other projections, including the weighted denominator `D` with weight
1.5, were recomputed. Zero deltas were skipped. The profile showed refresh
falling from 11–13% to 4.5% of the time. The wall-clock gain on the generic
path was 4%, and both sparse workloads, which have no projections, were
9% slower. T3 was reverted.

**T4 — event loop.** No candidate was large enough to measure on its own.
The loop's own 19% sits in inlined code that gprofng could not attribute
to lines. The remaining items (polling, report checks, live-instance
accounting, rule-event counters and the stop dispatch) are a few
nanoseconds each.

**T5 — decision.** After T1–T4 the gap was still about 3×, so the loop
moved on to T6. The decision record was amended first with the T6
criterion: merge only with every gate passing and at least a 25% gain on
the generic opt-in path.

**T6 — grouped evaluation, as a prototype.** A temporary environment
switch selected a grouped evaluator in the Static Direct loop, so the
contract stayed untouched until the gain was known. Channels were grouped
by rate base: one group for constant rates and one per vacancy base. Each
channel's constant factors and combinatorial factor were folded into one
coefficient at preparation. Per event, the evaluator recomputed the terms
of channels whose reactants changed and the sums and smallest terms of the
touched groups. It then formed each group rate as base × sum.

The first variant drew a second uniform inside groups with more than one
positive channel. The second continued with the first uniform and checked
every positive channel rate, base × term, against the total, as flat
selection does. On the sparse model, with one constant group, the second
variant reproduced Static Direct's output byte for byte. On mtDNA it took
103 ns per event against 107 ns on `main` (1.04×), and default Static
Direct gained 1.12×. Profiles showed the time moving from per-channel
evaluation into the group update without shrinking. The prototype was
reverted and saved as a patch. No rules-run/0.5 was added, and the T6
gates were not run.

## Verification

The branch changes no crate source. It adds the probe script, the
iteration record, the attempt patches and documentation. Every reverted
attempt's identity was checked with the probe's byte comparison. The
96-case rules replay was not run for reverted attempts, because none was
kept. `cargo fmt`, Clippy with warnings denied, the release workspace tests,
ruff format and check, and the Python suites (normally and with `-O`, with
identical reports) passed on the branch head.
