# Review of the mtDNA A0–A6 patch

3 November 2026. Review of the contributed A0–A6 patch and its fixed-model
performance follow-up (12 commits on base `a34484b`, tree `861d295`). The
contributed implementation log is
[mtdna-a0-a6-implementation-log-2026-11-03.md](mtdna-a0-a6-implementation-log-2026-11-03.md).

## Method

- Applied the ordered mail series; the tree matched the stated hash.
- Ran fmt, Clippy with warnings denied, 224 release tests, the dependency
  gate and every Python suite of the CI jobs, each under normal Python and
  `-O`. All passed before any change.
- Compared a build of the base with a build of the patch on 416 default
  invocations: `describe` for every schema, `check`, `analyze`,
  `analyze-v1`, `export`, `simulate`, `simulate-cached` and `fsp` on the
  example models, `ssa` for all five aggregators × two waiting times ×
  ssa-run 0.1, 0.2 and 0.3 with every `fast` setting, `tau`, and rules-run
  0.1 on every base rules model at one and three workers. Output was
  identical except the five intentionally versioned schemas and the
  archive's executable hash.
- Four independent read-only reviews, split by area: the specialized
  kernel; static lowering, next reaction and observers; tau, binomial and
  ATP thinning; the network engine and the Python binding.
- Checked provenance: no Java source, class names or private paths are in
  the patch. The harness only invokes a separately built comparator.

## Findings and changes

Each fix has a test that failed before the change.

| Area | Finding | Change |
| --- | --- | --- |
| Contracts | Unit variants of the rules `method` ignored extra keys (`{"kind": "static_direct", "epsilon": 0.5}`) | Empty struct variants; the JSON form is unchanged |
| Contracts | rules-run/0.1 accepted `trackers: []`, `method: null`, `event_log: null`; rules-model/0.1 accepted `atp: null` | Presence is recorded and rejected under 0.1; a present field must hold a value |
| Contracts | ssa-run/0.4 defaulted `fast` to false for every aggregator, a third default | 0.4 keeps the 0.3 defaults; next reaction defaults to false and rejects `fast: true` |
| Contracts | Next-reaction rejections returned CRN-E060 with an unrelated hint | CRN-E020 at `/method/aggregator`, `/method/fast` or `/method/waiting_time` |
| Contracts | Discovery gained keys without a version bump | describe/0.10 (core, Python) and 0.11 (CLI) |
| Python | `ssa_spec(aggregator="next_reaction")` was rejected with its own defaults | The default policy resolves per aggregator |
| Kernel | The resolution bound paired two extremes that cannot occur together; K = 1e5, fracAdv = 1e-6 or mu = 1e-7 fell back silently | Per-state bound with the same 16 × 2⁻⁵³ margin |
| Kernel | A fallback was visible only in the sampler string | CRN-R051 warning with the reason |
| Kernel | An overflowing wait failed the path; Static Direct completes it | Treated as beyond the horizon |
| Kernel | A Continue observer changed tracker integrals in the last bits | Integration independent of observers |
| Static | `-log(1 - u)` with u = 0 gave a zero wait and a failed path | Open uniform; other words unchanged |
| Static | Sub-resolution channels reported CRN-R030 | CRN-R031, pointing to the tree engine |
| Observers | Paths failing before tracking emitted empty tracker lists | No observation for such paths |
| Messages | Thinning on a model without ATP said the reverse; unknown log rules reported at `/trackers` | Corrected message and path |
| ATP | ATP at zero stopped at once even under positive forcing | Depletion only when the forcing is not positive |
| ATP | The zero-crossing time overflowed to infinity for a huge ratio | Logarithm term by term |
| ATP | The span cap fell below clock resolution near a crossing | Cap dropped when it cannot advance; the endpoint bound stays valid |
| Tau | Species-only step control overfilled the vacancy (N up to 1013 > K) | CGP propensity condition with first differences, and rate criticality |
| Tau | The K = 20 tau rows never leaped and so tested exact SSA | Documented; new production gate at K = 1000 |
| Binomial | p = min(rate · dt, 1) saturated: pure decay at rate 1 left zero instead of 1000/e | Competing-risk probabilities |
| Binomial | Steps were cut at report times, so the grid changed the path | Fixed lattice k · step; reports hold the last step |
| Python | NumPy arguments raised a raw TypeError | Converted; counts must be integral |
| Python | `tracker_columns` returned empty arrays for a bad index | IndexError |
| Python | The live callback checked the cancel file on every event | Flag only; the engine polls the file |

The two approximate samplers moved to `/2` identities. Exact samplers keep
their identities: the Static Direct and next-reaction change alters only
paths that previously failed.

## Pre-registered tau gate

The gate was declared in
[the decision record](../decisions/mtdna-a0-a6-2026-11-03.md) before it
was run: K = 1000, μ = 3 × 10⁻³, 2048 tau paths (ε = 0.01) against 2048
Static Direct paths, nine CDF comparisons with the two-sample DKW triangle
bound 0.0978 (family α = 0.001), and no report with N above K.

- Attempt 1, with the CGP floor max(ε·a_j, c_j), kept N ≤ K but failed N at
  250 h: 0.0991 against 0.0978 (`validation/mtdna/gates-review-attempt1.json`).
- With a few vacancies, one degradation moves the replication rate by about
  a third. A single firing that exceeds ε·a_j now forces an exact step, as
  a submolecule species allowance already did.
- The unchanged gate then passed all ten rows, largest distance 0.0264
  (`validation/mtdna/gates-review.json`, 122 checks).
- A 16-path logged cohort had no multi-firing leaps. At this copy number tau
  runs exact events and gives no speedup. The high-copy immigration/death
  gate still requires and finds real leaps.

## Validation after the changes

- 245 Rust release tests (224 before, 21 new), Clippy with warnings denied
  and rustfmt pass. The new tests are in
  `crates/crn-craft-rules/src/tests_review.rs` and
  `crates/crn-craft/tests/ssa_alternatives.rs`.
- The default-output comparison with the base gives the same result as
  before the changes.
- The Python suites pass under normal Python and `-O` with identical
  reports, as does the installed-package check.
- `mtdna_direct` on 512 cells (25,920 h, 90% stop, two class trackers,
  one worker) gives byte-identical output to the unreviewed head at the same
  speed: median 2.75 s against 2.81 s
  (`validation/mtdna/review-kernel-probe.json`). The eight-worker Java
  comparison was not repeated; this host has two CPUs and no Java.

## Limitations and open items

- Rust 1.85 and 1.88 could not be installed here; CI covers them. A manual
  scan found no newer APIs.
- The binomial sampler draws O(n·min(p, 1 − p)) uniforms per binomial, with
  cancellation polled once per step. That is fine at mtDNA copy numbers,
  but slow for counts near 2⁵³.
- Frozen rates over a binomial step can put N a few molecules above K.
- The tau step's cost grows with the square of the channel count.
- Static paths report a count above 2⁵³ as CRN-R030, the tree as CRN-R042.
- ATP stops and model stops share the `stopped` termination; ATP is zero on
  every report after an ATP stop.
- The heap behind next reaction has no randomized unit test; the exact-law
  gates and the new replay test cover it.
