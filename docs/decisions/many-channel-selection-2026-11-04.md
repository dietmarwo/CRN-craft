# Selection for many-channel networks

4 November 2026. Declared before any change was measured. It follows the
T2a result of the
[structural iterations](generic-rules-structural-2026-11-03.md). Summing the
rates and tracking the smallest positive rate in one pass made Static
Direct on the 96-channel sparse network 1.31× faster. The same change made
the 8-channel mtDNA probe 3% slower, so the rule for that loop rejected it.

## Change

`direct_with` computes the total and the smallest positive rate in one
pass when a network has at least 32 channels. The threshold was fixed
before measuring, between the two networks measured. The resolution check
then compares the smallest positive rate with `total × 2⁻⁵³`, which is
equivalent to testing every rate. Summation order, error order and random
words are unchanged, so results stay byte-identical. Networks with fewer
channels keep the current two passes.

The generated loops use the same rule. Unrolled and compact loops for
networks of at least 32 channels track the smallest positive rate while
summing, and smaller ones keep the per-rate test.

## Measurements

`validation/probe_structural_rules.py` gains a wider workload: 256
attributed instances that replicate and die, 512 channels. (Plain species
were declared first, but the model limit is 64 types.) The interpreter is measured with the
probe against `main`. The generated loops are measured with
`validation/check_generated_loop.py` on the sparse network, before and
after.

## Rules

- **Byte identity.** Identity on the probe's workloads, the 96-case rules
  replay and the generated-loop replay matrix.
- **Interpreter.** Kept only if sparse Static Direct gains at least 5% and
  no probe workload is more than 5% slower.
- **Generated loops.** Kept only if the generated loop on the sparse
  network gains at least 5% and the mtDNA loop is not more than 5% slower.
- **Records.** Failed attempts are recorded in
  `validation/mtdna/generated-loop-followups.json`.

## Outcome

Recorded 4 November 2026. The interpreter change is adopted, and the change
to generated loops is not. Runs are in
`validation/mtdna/generated-loop-followups.json`.

**Interpreter.** Against `main` (`a78a3a5`), with nine repeats, sparse Static
Direct went from 254 to 196 ns per event (1.30×) and the wide workload from
924 to 666 ns (1.39×). The mtDNA workloads and next reaction stayed within
2%. The probe's outputs and the 96-case rules replay were byte-identical. In
the final seven-repeat run the figures were 1.25× and 1.34×, and no
workload was more than 3% slower.

**Generated loops.** Tracking the minimum in unrolled code made the sparse
loop 1.39× slower (222 to 309 ns per event) and the wide loop 1.86× slower
(967 to 1799 ns), so unrolled loops keep testing every rate. In the compact
layout the fused pass helped (1.08× sparse, 1.32× wide), but that layout
was not adopted
([record](generated-loop-stops-and-compact-2026-11-04.md)).
