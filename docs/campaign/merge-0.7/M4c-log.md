# M4c log: Static Direct selectors

Task M4c of [the plan](../../plan-0.7.0-engine-merge.md): event selectors
for Static Direct in `crn-craft.run/1.0`, ported from the 0.6.0 network
engine. The work was interrupted twice by server restarts; "Run history"
says what that affected.

## Step 0: the M4e merge

On `206fd0f` (M4e merged by the coordinator), wheel and CLI built from the
root, under the CPU lock: `check_rules_python.py` (25 checks),
`check_installed.py` (`python -I` from `/tmp`, 17 checks),
`check_contracts.py` (90 of 90), I2 `check_manual.py --require-cargo` with
`CRN_CRAFT_SOURCE` set (22 blocks) and `check_compiled_loop.py
--require-cargo` (57 rows) all pass, each byte-identical under `python -O`.
No merge problem was found, so there is no fix commit.

## What changed

### Contract (`crates/crn-craft/src/rules/mod.rs`)

- `Method::StaticDirect` gains `selector: Option<Selector>`, with the public
  enum `Selector` (`linear`, `composition_rejection`, `sorting`, `rssa_cr`,
  `auto`). Absent means `auto` in `crn-craft.run/1.0`.
- Every `rules-run` version refuses the field: `CRN-R020` at `/schema`,
  "static_direct.selector requires crn-craft.run/1.0" (the convention of the
  quantile gate of M4e), and always runs the linear scan. `null` and unknown
  values are `CRN-E003`; other methods have no such field (`CRN-E003`).
- Schemas: `ONLY_1_0` now also removes properties of tagged-enum variants
  (`oneOf` entries), and lists `("Method", "selector")`;
  `ONLY_1_0_DEFINITIONS` gains `Selector`. The `rules-*` schemas are
  byte-identical; `run-1.0` and `result-1.0` were regenerated with
  `crn-craft describe`. The run-1.0 description names the selector and its
  default.
- Sampler identity: the linear scan keeps the 0.6.0 identities; another
  selector adds its name before the sparse-tracker suffix:
  `chacha8-u53[-ziggurat]/rules-static-direct{,-cr,-sorting,-rssa-cr}[-sparse-trackers]/1`.
  `auto` reports the selector it picked. The result's `method` echoes the
  request (with or without `selector`).
- Discovery: `describe()["contracts_1_0"]["selectors"]` (values, default,
  `loops: ["linear"]`, and the `auto` thresholds).
- `docs/diagnostics.md`: the new `CRN-R020`, `CRN-R030` and `CRN-R031`
  cases.

### Engine

- `rules/select.rs` (new), ported from `ssa/selectors.rs`:
  - `RateTree`: pairwise sums of the channel rates with the lazily checked
    strictness of 0.6.0's checked tree. Leaves are staged and `commit`
    recomputes their ancestors one level at a time, each once (the same
    sums as per-leaf updates; tested against the fully checked reference
    tree of 0.6.0). Its errors follow the linear scan: a non-finite total
    is `CRN-R030`, a positive rate lost in the sums or below 2⁻⁵³ of the
    total `CRN-R031`.
  - `Composition`: groups by binary exponent, envelope sums rebuilt after a
    membership change, Lemire's joint member-and-acceptance draw for groups
    of up to 1,024 members.
  - `sorting_sample`: the linear scan in a learned order; a threshold at or
    above the sequential sum is redrawn.
  - `resolve` and `auto` (below).
- `rules/rssa.rs` (new), ported from `ssa/bounds.rs`: species intervals
  (± clamp(n/10, 4, 64), points for zero counts and pinned species), dirty
  readers after an interval exit, lower and upper bounds and the smallest
  positive value over the box. The network engine had four rate laws; here
  rates are expressions, so the bounds are built as follows:
  - count factors of channels of at most two molecules are products of
    correctly rounded operations on nonnegative values, monotone in each
    count, so the box corners are exact bounds;
  - a scaled combination count (`combination_rate`, `Multi` channels) is
    within 2m + 1 roundings of the exact count for m ratios; its corners
    are widened by twice that bound, rounded outward;
  - a rule whose rate reads one species is evaluated at every count of that
    species' interval (at most 129 values) with the engine's own
    arithmetic (`Network::global_fresh`, tested bit-equal to the cached
    path), so its range over the box is exact for any expression; ranges
    are cached per rule and interval;
  - a rule whose rate reads two or more species pins them (their
    intervals are their counts), so its bounds are the exact value at the
    current state, recomputed whenever a read species changes;
  - a box where a value is non-finite, negative or below 2⁻⁵³ of the upper
    bound pins the channel's reads for the rest of the path (0.6.0's point
    fallback), and the point bound is the engine's value at the current
    state, so an invalid current state fails exactly as in the linear scan.
  So RSSA-CR runs every model Static Direct runs; nothing is refused.
- `rules/flat_path.rs`: `simulate_selected`, the Static Direct loop for the
  three new selectors; waits (inverse and Ziggurat), the clock check,
  saves, the horizon, events, observation, trackers and stops are those of
  the linear loop. RSSA-CR thins: every trial draws a wait against the
  upper total and is accepted below the lower bound, otherwise against the
  exact propensity (a value outside its bounds fails with `CRN-R030`); a
  zero-propensity state with positive upper bounds is detected as in
  0.6.0. The linear loop is unchanged.
- `rules/flat.rs`: `Network` keeps `readers`, `rule_reads` and
  `rule_projections`; `project` (shared with `refresh_indices`),
  `needs`, `constant_global`, `scaled_combination`, `global_fresh` and
  `value_with`; the count factor of `RateCache::channel` moved into
  `ChannelRate::count_weight` with the same operations.
- Generated event loops and rates kernels run only the linear scan:
  `execute_with_rules_loop` and `execute_with_rules_kernel` refuse another
  selector, also one that `auto` picked, with `CRN-R020` at
  `/method/selector` before the fingerprint check. Generation does not
  depend on the selector.
- After profiling (below), `Bounds::bound` reuses a buffer instead of a
  vector per call (about 12% of the instructions on W6) and
  `RateTree::commit` skips sorting sorted lists; results byte-identical on
  the nine profiled requests.

### `auto`

Fixed from the declared benchmark (`validation/merge/selector-threshold.json`):
`AUTO_COMPOSITION_REJECTION = 512` channels, `AUTO_RSSA_CR = None`. So
`auto` picks composition–rejection from 512 lowered channels on and the
linear scan below; sorting and RSSA-CR are selected explicitly. On the
workloads `auto` picks the linear scan for W1–W5 and W7 and
composition–rejection for W6.

### Python

`python/crn_craft/_rules_loop.py`: the compiled-loop self-check of a
`crn-craft.model/1.0` model names `"selector": "linear"`, so
`compile_loop` also works for models where `auto` would pick another
selector. Nothing else in Python changed; a `method` dictionary with a
`selector` reaches the engine as written.

### Validation scripts (`validation/merge/`)

- `selector_threshold.py`: the declared benchmark; `--resume` keeps
  finished points.
- `check_selectors.py`: the exact-law check per selector and waiting time.
- `compare_speed.py`: same-session timings against the 0.6.0 best, with
  `--selectors` (M4d can use it).
- `check_laws.py --state DIR`: stores each finished chunk and loads it on a
  rerun with the same engine, model, method and sizes; reports are
  byte-identical with and without stored chunks (tested on a smoke run).

## Tests written first, and how they failed first

- `crates/crn-craft/tests/selectors.rs`, written against the unchanged
  engine: 10 of 10 failed (`CRN-E003` "unknown field `selector`"). They
  cover the field and its refusals, the sampler names (four selectors, two
  waits, sparse trackers), the equality of `linear` under `run/1.0` with
  `rules-run/0.6` apart from schema fields (below), a quick
  immigration–death law per selector and wait, thread independence (1, 2,
  3 threads on W5, W4 and repressilator), every channel shape (sources,
  single, pair, cross with factors, `Multi` with a count and three
  patterns, whole propensities reading one and three species) with
  trackers, logs and invariants, absorbing states, `CRN-R031`, the
  generated-loop refusal and the schemas.
  `auto_follows_the_declared_threshold` was added after the benchmark and
  failed first against the provisional constants (a 511-channel ring got
  composition–rejection).
- The ported unit tests (`select.rs`: the categorical DKW test after
  membership changes, the sorting redraw, exponent groups, the lazily
  checked tree against the fully checked one (now also comparing every
  committed node), drift; new: large groups, sorting order, a single
  channel) were written in the file before it was wired in, so they first
  ran with the implementation, as in M4b.
- Brute-force bound tests (`rssa.rs`): random boxes over models of every
  shape, rate form and observable reduction; every integer state of each
  box is evaluated with the engine's channel arithmetic and must lie in the
  bounds, positive values above the stated minimum; counts and observable
  values are restored. Also: states near the count limit for scaled
  combination counts (k = 2 … 8), the cached enumeration, pinning of
  rules reading two species, the point fallback and the error at an invalid
  current state, and the fresh global factor against the cached one.
- Mutation check of the bound tests (each mutation applied alone, then
  reverted): upper bound at the low corner, enumeration missing its last
  count, minimum positive at the high corner, and the global minimum
  replaced by the maximum were caught. Two survived: removing the widening
  of scaled combination counts (no nonmonotone step of `combination_rate`
  exists to observe: 30,000 random windows of 64 counts between 2⁵⁰ and
  2⁵³, k = 2 … 16, found none, so the widening is a proven margin the brute
  force cannot exercise), and ignoring the availability gate for the lower
  bound (an equivalent mutation: `value_with` is 0 at a corner below a
  reactant's need).
- The M3 equivalence test (`contracts_v1.rs`) now names `"selector":
  "linear"` in its run/1.0 Static Direct requests and removes the echoed
  `selector` before comparing; it holds for all 17 stored requests.

## Declared and measured

Declarations in `validation/merge/RESULTS.md` ("M4c (selectors):
declaration", commit `2d4cde6`), before any run; results in the same file.

| Check | Result |
| --- | --- |
| Selector-threshold benchmark (5 families × 10 sizes × 4 selectors) | N_cr = 512, N_rssa none; `auto` fixed as declared |
| Exact-law check (168 DKW tests, 16 power checks) | 168 of 168 pass, largest distance 0.00499 (ε 0.00798); power 16 of 16; `-O` identical |
| Law gate `linear` | 4,814 comparisons, **0 failures**, max \|z\| 3.88; power sums and events equal the M4b regression run for all twelve models |
| Law gate `composition_rejection` | 4,814 comparisons, **0 failures**, max \|z\| 3.83 |
| Law gate `sorting` | 4,814 comparisons, **0 failures**, max \|z\| 3.51 |
| Law gate `rssa_cr` | eleven models: 4,724 comparisons, 0 failures, max \|z\| 3.94; **vilar-noise failed** (below) |

**The benchmark rule excluded RSSA-CR.** RSSA-CR was the fastest selector
at every point of the `ring` and `michaelis` families and at every `dense`
point from 12 channels on (for example 414 against 3,347 ns per event for
composition–rejection at 4,032 dense channels), but slower than
composition–rejection at every `ring_small` point (four molecules per
species), by 2–4% at the two largest. The declared rule required RSSA-CR to
win in every family, so N_rssa is none. A rule that also reads the initial
populations would pick RSSA-CR for the large-count and dense families; it
needs a new declaration and benchmark run (not done here).

**RSSA-CR on vilar-noise.** Chunk 1 (root seed 2,000,001) failed:
replicate 6298 stopped at t = 7.1624 with `CRN-R030` "clock failed to
advance" (`/time`). One investigation (a temporary print, not committed):
the upper total was 9,729.05 (exact total 3,353.06), and the Ziggurat wait
unit / 9,729 was below half an ulp of t, so t + wait == t. The linear scan
has the same check; RSSA-CR meets it more often because it draws a wait for
every trial against a larger total. The 0.6.0 engine does the same: its
RSSA-CR trial loop (`ssa.rs`, `run_path`, mode 4) computes
`next = sampling::event_time(time, wait)?` before the composition draw and
the acceptance test, and `event_time` fails the replicate (`CRN-E062`,
"event time cannot advance in finite f64 arithmetic") when the clock does
not advance, also for a trial that would be rejected. Decision: keep these
semantics; the engine is unchanged, the failed row is kept
(`validation/merge/laws/m4c-rssa_cr-vilar-noise-failure.json`), the other
18 chunks of vilar-noise were not run, and nothing was rerun with other
seeds. **The declared RSSA-CR law gate is therefore not passed.**

Options for the owner: (a) accept it as an inherent floating-point limit
of RSSA-CR on fast models at late times, document it, and exclude
vilar-noise from RSSA-CR's gate; (b) let a rejected trial whose clock does
not advance be redrawn without failing (a departure from 0.6.0; accepted
events would still fail as in the linear scan), then rerun the vilar-noise
RSSA-CR gate with the declared seeds; (c) carry the time as a compensated
sum. `auto` never picks RSSA-CR, so the default is not affected.

### Timings (informative; table in RESULTS.md)

Same-session ratios to the 0.6.0 best (CLI, one thread, speed-gate sizes),
best selector: W1 1.32 (`linear`), W2 1.21 (`rssa_cr`), W3 1.92 (`linear`
through `auto`), W4 2.36 (`rssa_cr`), W5 1.27 (`rssa_cr`), W6 1.10
(`rssa_cr`), W7 1.06 (`rssa_cr`). Before M4c (coordinator, linear only):
W1 1.85, W2 1.98, W3 1.89, W5 4.15, W6 5.41, W7 1.64. Row timings varied
by up to about 30% between sessions on this host. Against the 1.25 gate of
M4d, W1, W3, W4 and W5 are over and W2 is close.

## Where the time goes (input for M4f)

Callgrind (valgrind 3.22, release CLI without debug info, so inlined code
counts in its caller), one `crn-craft rules` run at the speed-gate sizes,
self instructions. Instructions per event: W1 `linear` 644, W2 `rssa_cr`
969, W3 `linear` 1,497, W4 `rssa_cr` 2,879 (2,855 after the buffer
change), W5 `rssa_cr` 2,132, W5 `composition_rejection` 3,095, W6
`rssa_cr` 3,686 (3,232 after the buffer change), W7 `rssa_cr` 1,736.

| Workload | Largest self costs |
| --- | --- |
| W1 `linear` (2 channels) | `direct_with` (sum, scan and inlined Ziggurat wait) 25%, the inlined linear loop in `Path::run` (environment, rate-cache reset, polling, saves) 24%, `apply_flat` (checked count update, live-instance and limit bookkeeping, rule counters) 18%, `RateCache::channel` + `count_weight` + `global` 16%, ChaCha8 words (`Words::word` + refill) 11%, `refresh_channel` 5% |
| W2 `rssa_cr` | inlined selector loop (wait, composition draw, acceptance) 42%, ChaCha8 22%, `apply_flat` 12%, `Bounds::refresh` 11%, `RateTree::total` 6% |
| W3 `linear` (Hill rates) | expression interpreter `eval_node` 28% plus `libm::pow` 3%, `direct_with` 23%, linear loop 12%, `RateCache::channel` 9%, `apply_flat` 8% |
| W4 `rssa_cr` (Michaelis–Menten, order 8) | `eval_node` 30% (the bound enumeration of each rule's single read species), `Bounds::refresh` 24%, `apply_flat` 15% (eight species per event), selector loop 12%, `project` 5%, `global_fresh` 3% |
| W5 `rssa_cr` | selector loop 37%, ChaCha8 21%, `apply_flat` 8%, `Bounds::refresh` 7%, exact evaluations at acceptance 6%, `RateTree::total` 5%; model hashing and allocation about 4% (once per run) |
| W5 `composition_rejection` | `RateTree::total` with the level-wise commit 24%, the rate-update closure 17%, the loop 13%, `Composition::update` 11%, `RateCache::channel` 9% |
| W6 `rssa_cr` | `Bounds::refresh` 36% (31 readers per exiting species), allocation in `bound` about 15% (removed), `RateTree::total` 5%, `Composition::update` 5%, `count_weight` 4% |
| W7 `rssa_cr` | selector loop 31%, ChaCha8 15%, `Bounds::refresh` 13%, `apply_flat` 11%, `RateTree::total` 6% |

Pointers for M4f, none done here:

- Per-event overhead common to every selector: `apply_flat` (8–18%), the
  random words (11–22%; ChaCha8 is fixed by D4, so fewer words per event
  is the lever, for example a joint wait-and-selection draw), the rate-cache
  reset and environment construction per event, and polling.
- Linear scan: cached prefix sums resumed at the first affected channel
  (0.6.0's `cached_direct`, the 0.6.0 best on W1 and W3) would change the
  scan's arithmetic, so it needs its own sampler identity or a proof of
  equal selection.
- Expressions (W3, W4): `eval_node` dominates; specialized kernels for the
  converted Hill and Michaelis–Menten forms, or a cheaper interpreter.
- `RateTree::commit` for composition–rejection: per-leaf ancestor updates
  for small staged sets, or a cheaper exact total.
- RSSA-CR on W4: enumerate the bound range once per rule and interval
  (it is cached, but intervals of the eight reactant species move often).

## Checks

On `4fe800d` (wheel `_native` SHA-256 `4f14a81d…`), under the CPU lock.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass, 28 test binaries: `crn-craft` lib 137 (16 new in `select.rs` and `rssa.rs`), `selectors` 11 (new), `contracts_v1` 7, `higher_order` 8, `quantiles` 13, `limits_v1` 5, `convert` 10, CLI 14, I1 (`stored_requests` 2, `static_direct_defaults` 9), I3 (`generated_loop` 26, `generated_kernel` 3), all others |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (the known output-filename collision warning) |
| Schemas against `crn-craft describe` | `rules-model`, `rules-run`, `rules-result`, `rules-check`, `model-1.0`, `check-1.0` byte-identical; `run-1.0` and `result-1.0` regenerated |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| `check_dependencies.py` | pass; `-O` identical |
| `check_contracts.py` | 90 of 90; `-O` identical |
| `check_rules.py` | 47 of 47; `-O` identical |
| `check_rules_python.py` | 25 checks pass; `-O` identical |
| `check_installed.py` (`python -I`, from `/tmp`) | 17 checks pass; `-O` identical |
| `check_mtdna.py --binary target/release/crn-craft` | pass; `-O` identical |
| `check_builder.py` | 18 tests, 0 failed; `-O` identical |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` (`CRN_CRAFT_SOURCE`) | 22 blocks pass; `-O` identical |
| `check_compiled_loop.py --require-cargo` (`CRN_CRAFT_SOURCE`, `CARGO_TARGET_DIR`) | 57 rows pass; `-O` identical (I3 by the Rust fixtures) |
| `check_hybrid.py` | 0 failed; `-O` identical |
| `validation/merge/check_selectors.py` | 168 of 168; `-O` identical |

## Run history

- First restart: during the threshold benchmark, after 27 of 50 points
  (each written when finished). The script gained `--resume`; the other 23
  points ran afterwards, one family per command. The first 27 points ran on
  a CLI built before `cargo fmt` and an MSRV-compatible `next_up` helper
  (only higher-order widening differs; no such channel in the benchmark).
- The law gates ran in pieces of under ten minutes with `check_laws.py
  --state`; the second restart stopped the `sorting` gate during
  vilar-noise, which resumed from its stored chunks.
- Working files (profiles, stored chunks, logs) are in
  `/home/claude/m07/m4c/`.

## Limitations

- RSSA-CR fails a replicate when a trial's wait is below half an ulp of the
  time (0.6.0 semantics); on fast models at late times this is far more
  likely than for the linear scan (vilar-noise above).
- RSSA-CR pins every species read by a rule that reads two or more
  species; on models made of such rules it behaves like
  composition–rejection with exact recomputation, slower than
  composition–rejection itself (`two_reads` family).
- The widening of scaled combination counts is a proven error bound that
  the brute-force tests cannot exercise (no nonmonotone step found).
- `auto` never picks RSSA-CR or sorting, although RSSA-CR is the fastest
  selector on most measured families and on W2, W4–W7.
- Sorting is never the best selector on the workloads in this engine.
- Generated event loops and rates kernels run only the linear scan.
- Selector workspaces are allocated per replicate (a composition table of
  2,098 groups among them); negligible at the measured sizes.

## Deviations from the plan

- Sorting and RSSA-CR are ported, as the M0c selector rule decided.
- `auto` is decided by a declared benchmark over five synthetic families,
  with a rule that could also pick RSSA-CR; the plan's form ("linear below
  a channel-count threshold, composition–rejection at or above it") is what
  the measurement produced.
- RSSA-CR runs every Static Direct model (point bounds where an interval
  bound is not available) instead of refusing channels.
- The old-version refusal is `CRN-R020` at `/schema` (M4e's convention),
  not 0.6.0's `CRN-E003` unknown field.
- The M3 equivalence test names `"selector": "linear"` and ignores the
  echoed selector.
- File layout: `rules/select.rs` and `rules/rssa.rs` (the plan names
  `select.rs` only).
- `check_laws.py` gained `--state`, `selector_threshold.py` `--resume`
  (restarts).

## Migration notes

- New in `crn-craft.run/1.0`: `method.selector` for `static_direct`
  (`linear`, `composition_rejection`, `sorting`, `rssa_cr`, `auto`;
  absent means `auto`). A run/1.0 Static Direct request without a selector
  on a model of 512 or more lowered channels now runs composition–rejection
  (sampler `rules-static-direct-cr`), so its numbers differ from
  `rules-run/0.6`; `"selector": "linear"` reproduces them. `rules-run`
  requests are unchanged and refuse the field (`CRN-R020` at `/schema`).
- Rust: `Method::StaticDirect` has a third field, `selector`; struct
  literals must add `selector: None`. New public enum
  `crn_craft::rules::Selector`.
- Generated event loops and rates kernels refuse a request whose selector,
  or `auto`'s choice, is not `linear` (`CRN-R020` at `/method/selector`).
- `describe()["contracts_1_0"]["selectors"]` lists the values and the
  `auto` thresholds.
