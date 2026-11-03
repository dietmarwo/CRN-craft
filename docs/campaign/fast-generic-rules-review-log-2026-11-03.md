# Review of the fast generic rules patch (S1–S6)

3 November 2026. Contributed series of three commits, which implements the
plan "fast generic rules simulation instead of hand-written kernels". The
contributed log is
[fast-generic-rules-implementation-log-2026-11-03.md](fast-generic-rules-implementation-log-2026-11-03.md).

## Integration

- The series was built on the pre-review head `6edc100`, whose tree equals
  `6501319` in this repository. It applied there cleanly and reproduced the
  stated tree `abf845c`.
- `main` (`f255ee9`) was then merged in. Conflicts were in `Method`, the wait
  draw, the observation spec, version checks and sampler identities. Both
  sides were kept:
  - the empty struct variants from `main`, with presence-checked
    `waiting_time` and `sparse_trackers` options;
  - the open uniform for waits;
  - the `/event_log` error path;
  - the `/2` approximate samplers;
  - the explained `mtdna_direct` fallback.
- The tau propensity step now uses its own rate cache, since
  `Network::rate` became test-only.

## Verification

- **Replay identity.** 96 rules requests give byte-identical output to
  `main`. They cover:
  - methods: tree, Static Direct, next reaction, `mtdna_direct`, tau,
    binomial and ATP thinning;
  - with and without the takeover stop, trackers, first-passage predicates
    and event logs;
  - the CME fixture and the five base rules models;
  - one and three workers.
- **Independent reviews** of the caching and observation changes and of the
  predicate bytecode and generated kernels.
  - Adversarial cases: more than 60 against `main`, all identical.
  - Fuzzing: 4,000 random predicates over grids with NaN, ±∞, ±0 and
    subnormals, all bit-identical to the tree evaluator.
  - Generated kernels: bit-identical on two models with four request
    variants.
- **Speed on this host** (two CPUs; 128 cells, 25,920 h, takeover stop, two
  class trackers, one worker):
  - Default Static Direct: 3.65 s on `main` against 2.37 s here (median of
    three, alternating), a 1.54× speedup with identical output.
  - Ziggurat waits plus sparse trackers: 2.07 s.
  - `mtdna_direct`: 0.76 s.
  - Evidence: `validation/mtdna/review-fast-generic-probe.json`.

## Findings and changes

Each fix has a test that failed before the change.

| Severity | Finding | Change |
| --- | --- | --- |
| Medium | The generated-kernel fingerprint covered the model bytes and channel shapes but not the emitted arithmetic; a kernel from an older generator could run silently after a change in lowering | Fingerprint = hash of the emitted source keyed by the model bytes; checked by regeneration before each run |
| Low | No Rust test covered generated kernels; the script that builds one is not in CI | Checked-in kernel for the three-class fixture with golden, replay and rejection tests under `cargo test` |
| Low | An observer declaring `ENABLED = false` was still called and, by stopping, could hold stale report-only observables | Such observers are not called |
| Low | `sparse_trackers: true` with no eligible tracker reported the sparse sampler identity for a default run | Default identity and a `CRN-R051` warning |
| Low | The public kernel API was undocumented | Doc comments on the trait, generator and runner |
| Merge | A series test expected CRN-R030 for a sub-resolution channel | Now CRN-R031, as on `main` |

## Remaining gap

The generic path still trails the specialized sampler: 2.7× here and 3.3× on
the contributor's host. A profile of the generic opt-in path attributes:

- 33% to rate evaluation;
- 18% to selection, with its per-event resolution check over every
  channel;
- 11% to observable refresh.

`mtdna_direct` certifies resolution once at preparation and selects among
grouped outcomes. Closing the gap would take a certified resolution bound or
a grouped selector in the generic path. The generated kernel gave 1.1× on
the contributor's host.

## Limitations

- The contributor's eight-worker timings were not repeated; this host has two
  CPUs.
- `validation/check_generated_rules.py` (fresh build and timing) stays outside
  CI; the Rust tests cover correctness.
