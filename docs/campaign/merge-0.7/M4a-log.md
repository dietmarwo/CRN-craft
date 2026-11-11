# M4a log: limits of `crn-craft.model/1.0`

Task M4a of [the plan](../../plan-0.7.0-engine-merge.md).

## Declared before measuring

- Compile-time target: `prepare` of a 4,096-species, 8,190-reaction
  reversible chain (`crn-craft.model/1.0`, one observable per species,
  `static_direct` and `next_reaction`, `crn-craft.run/1.0`) takes under
  **2 s** on this host (Intel Xeon @ 2.10 GHz, 2 logical CPUs), release
  build, under the CPU lock. The target is the constant `CHAIN_PREPARE` of
  `crates/crn-craft/tests/limits_v1.rs`, written with the test before the
  implementation. The test ran once (pass) after the implementation and
  before this paragraph was written, with its output captured, so no time
  was seen; the measurement below is the first one read.
- Law gate for W5 and W7: declared in `validation/merge/RESULTS.md`
  ("M4a ... declaration"), committed before the run.
