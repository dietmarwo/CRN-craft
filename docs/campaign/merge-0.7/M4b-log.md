# M4b log: higher-order reactions

Task M4b of [the plan](../../plan-0.7.0-engine-merge.md): reactant
multiplicity (`count`) and any number of patterns in
`crn-craft.model/1.0`, the `Multi` shape, the overflow-safe combination
count, method support and refusals, and the converter for any
stoichiometry.

## Declared before running

In `validation/merge/RESULTS.md` ("M4b (higher-order reactions):
declaration", commit `39aecef`), before the FSP references existed and
before any run: the two exact-law models and their FSP bounds, the exact-law
check (methods, n = 200,000 per model and method, seeds 3,000,000 + k,
moment, state and remainder comparisons, factor 5, 0 failures expected), the
W4 law gate and the regression run over all twelve models. Unit-test
tolerances are stated in the test code before the implementation
(`combination.rs`: a binomial of m ratios within m + 1 roundings).

## What changed

- **Contract** (`rules/model.rs`): `Pattern` gains the optional integer
  `count` (1 to 2,147,483,647, `MAX_PATTERN_COUNT`; `null` rejected), and a
  `crn-craft.model/1.0` rule may have any number of patterns. Conditions,
  each a `CRN-R013` at the pattern's path: a pattern with `count` > 1, and
  every pattern of a rule with more than two patterns, is over an
  attribute-free type (`/type`) without `where` (`/where`) or `factor`
  (`/factor`); a count outside the range is refused at `/count`; patterns
  over one type must be identical (the existing "two different patterns over
  the same type" at `/reactants`); labels are unique over all patterns
  (`CRN-R010` at the repeated label).
- **0.x documents keep their diagnostics** (I1): a `count` in a
  `rules-model/0.1`–`0.4` document gives exactly the 0.6.0 strict-parser
  error (`CRN-E003` "unknown field `count`, expected one of `label`,
  `type`, `where`, `factor`" at `/rules/i/reactants/j/count`, same hint),
  for the first such pattern, from `check_model` and `prepare`; more than
  two patterns keep "more than two reactant patterns" (`CRN-R013`). The
  `rules-model` schema is byte-identical (the generator removes `count` from
  its `Pattern`); the `model-1.0` schema has `count` with minimum 1 and
  maximum 2147483647 and 1.0 texts for `reactants`, `label` and `type`.
- **Shapes**: a rule of at most two molecules keeps the 0.6.0 shapes and
  arithmetic: one pattern with `count: 2` is compiled as the pair of two
  identical patterns (`Shape::Pair`, event logs with two reactant entries),
  so it runs on every method, including `tree`, with the results of the
  two-pattern form. Every other rule with a count above 1 or more than two
  patterns is `Shape::Multi(Vec<(class, count)>)`, one entry per pattern.
  The plan's `u32` count is `u64` here, because patterns over one type add
  their counts.
- **Combination count** (new `rules/combination.rs`,
  `crn_craft::rules::combination_rate`): rate × factor × Π C(nᵢ, kᵢ) as the
  scaled mantissa/exponent product of the network engine's general
  mass-action law (`rates.rs`, "binomial(n,k) may overflow"), with the
  channel factor folded in as a second scaled factor; +inf once the exponent
  passes 1024 (bounded work: every ratio is at least 2, so at most about
  3,200 ratios); NaN, infinite and negative inputs carry through and are
  refused by the channel (`CRN-R030`). With factor 1 it gives the network
  engine's bits wherever the network engine gives a value.
- **Lowering** (`rules/flat.rs`): a `Multi` rule lowers to one channel per
  choice of a member for each pattern (patterns of one class are
  attribute-free, so at most one member); `Channel::multi` holds the merged
  `(species, count)` list sorted by species, `Channel::reactants` one
  species per pattern (the event participants), changes subtract the merged
  counts. The rate plan marks such channels (`HIGHER_ORDER`); the channel
  rate is 0 unless every species has its count, then the whole propensity
  (`global * (1.0 * factor)`, the 0.6.0 whole-rule arithmetic) or
  `combination_rate(global, factor, counts, list)`. The global factor
  computation moved into a helper with the same operations; every existing
  channel keeps its code path (I1, I3).
- **Methods**: `static_direct` (both waiting times, sparse trackers),
  `next_reaction`, `tau` and `hybrid` evaluate channels through the same
  rate cache and apply the merged changes, so they run `Multi` (tested
  below). `thinning` (ATP models) uses the same channel evaluation; no test
  combines ATP with a higher-order rule.
  `tau`'s criticality uses count × 10 for a higher-order reactant (a channel
  is critical within ten firings of losing a reactant), the old test for
  every other channel. `tree` (also as the default method) and
  `mtdna_direct` refuse a model with a `Multi` rule: `CRN-R020` at `/method`
  naming the rule. `binomial` refuses it through its existing scope check
  (`CRN-R020` at `/method`, "binomial splitting supports single-parent unit
  births and deaths"), since a `Multi` channel takes at least three
  molecules.
- **Generated event loops**: a higher-order channel's weight is emitted as
  its availability test and its rate as a call of
  `crn_craft_rules::combination_rate` with the global expression, the factor
  bits and the literal `(species, count)` list, inside the existing
  `if weight == 0.0` structure, so the loop runs the interpreter's
  arithmetic without restructuring the generator; whole propensities use the
  existing `global * (weight * factor)` line. The emitted text of every
  other channel is unchanged (I3). **Generated rates kernels** refuse
  higher-order rules (`CRN-R020` at `/method`), as they refuse whole-rule
  propensities.
- **Converter**: any stoichiometry converts. Reactions of at most two
  molecules keep their 0.6.0 mapping (unit patterns `x`, `y`; converted
  bytes of every earlier model unchanged); a reaction of three or more
  molecules gets one pattern per species with its stoichiometry as `count`
  (omitted when 1), labeled `x`, `y` for two species and `x1`, `x2`, … for
  more. W4 (order 8) converts. The order refusal of `CRN-R016` is gone;
  more than 64 product molecules are still refused.
- **Python** (`python/crn_craft/_rules_loop.py`): the self-check of
  `rules.compile_loop` probes a horizon with the default method (`tree`);
  for a `crn-craft.model/1.0` model it now probes with `static_direct`,
  since `tree` refuses higher-order rules. 0.x models keep the old probe.
  The higher-order case of `check_compiled_loop.py` found this (it failed
  with `CRN-R020` from the probe before the fix).
- Docs: rustdoc of the model, `Event::reactants`, `LogEntry::reactants`
  (through `concat!`, so the result schemas are unchanged), the generators;
  `docs/diagnostics.md` (`CRN-R013`, `CRN-R016`, `CRN-R020`).

## Tests written first

- `crates/crn-craft/tests/higher_order.rs` (8 tests), written before the
  implementation; all 8 failed first:
  `count_is_a_model_1_0_field` (the exact 0.6.0 `CRN-E003` for 0.1–0.4 from
  `check_model` and `prepare`, the 0.x three-pattern refusal, the 1.0 range
  1…2³¹ − 1, `CRN-E003` for −1, 2.5, null and a string);
  `higher_order_patterns_must_be_attribute_free_without_where_or_factor`
  (each refusal with its path, labels unique over three patterns, a count-1
  attributed pattern with `where`, `factor`, `copy` and an attribute
  expression beside a count-2 pattern runs: 2A + C → 2C keeps
  40 − A = 2(C − 5));
  `count_two_is_the_pair_of_two_identical_patterns` (`count: 2` and two
  identical patterns give identical results and event logs on `tree` and
  `static_direct`; three identical patterns equal `count: 3` apart from the
  event log); `every_static_method_runs_higher_order_rules_and_tree_and_mtdna_refuse`
  (`static_direct` with Ziggurat and inverse waits, `next_reaction`, `tau`,
  `hybrid` on 3A → ∅ (A ≡ 1 mod 3 at every report) and 2A + B → C
  (A − 2B = 0); `tree`, the default method, `mtdna_direct` and `binomial`
  refused at `/method`; the rates kernel refused; the event loop emits
  `combination_rate`); `propensities_count_combinations_and_whole_propensities_need_every_count`
  (mean event count over a short horizon against c · C(40, 3) · t within 5
  standard errors; a whole propensity fires exactly while A ≥ 3 and B ≥ 1);
  `overflowing_combination_counts_fail_the_replicate_with_a_coded_diagnostic`
  (the network engine's two overflow cases on replicates: C(2000, 1000) ×
  1e-300 runs; count 2³¹ − 1 of 2⁵³ molecules at rate 5e-324 fails with
  `CRN-R030` at `/rules/0/rate`; an unavailable count never evaluates its
  rate); `the_converter_maps_any_stoichiometry_with_count`;
  `the_rules_model_schema_has_no_count_and_model_1_0_has_it`.
- `rules/combination.rs` unit tests (5), in the file before it was wired
  in, so they first ran with the implementation: binomials up to n = 120
  within m + 1 roundings and exact for k ≤ 2 (declared there); the network
  engine's two overflow tests ported (`high_order_and_small_coefficient_…`,
  `…true_overflow_is_…_bounded_in_work`, now with a count of ratios: under
  4,000 for count 2³¹ − 1, also with the smallest factor); invalid inputs;
  and bit equality with the network engine's propensity on 20 cases of
  order 3 to 1000 (where the network engine refuses an overflow, the result
  is +inf). The first run of the bit test failed because the network engine
  refuses C(2000, 1000) at rate ≥ 1e-292 (test fixed to expect +inf there),
  and the binomial test because its exact u128 reference overflowed above
  n = 120 (range reduced, checked multiplication).
- The exact-law check `validation/merge/check_higher_order.py` failed first
  (both models refused by the M4a converter; kept as
  `laws/higher-order-before-m4b.json`).
- With the implementation: the generated-loop fixture
  `crates/crn-craft-rules/tests/generated/loop_higher_order.rules.json`
  (3A → ∅, 2A + B → C with an expression rate reading an observable, a whole
  propensity over 3A + B, a rule of three patterns, sources, decay, a stop)
  with its generated loop and
  `higher_order_rules_are_generated_and_replay_the_interpreter` (golden
  source; 3 Static Direct policies × 3 observation sets × 1 and 2 threads ×
  5 horizons and budgets, byte-identical to the interpreter, the stop
  fires); `validation/check_compiled_loop.py` compiles the same model
  through `rules.compile_loop` and requires byte-identical results with 1
  and 2 threads and both waiting times. The fixture first had `dc = 0.2`;
  its stop (nC ≥ 60) could not fire at that steady state, so `dc` became
  0.05 before the first replay run. `flat::tests::index_based_lowering_…`
  covers W4 (higher-order) since the converter accepts it.
- Updated: `tests/convert.rs` (`higher_order_reactions_convert_with_counts`
  replaces the refusal test), the CLI convert test (refusal case changed in
  M4a), the `model-1.0` schema snapshot.

## Measured

| Gate (declared in RESULTS.md) | Comparisons | Failures | max \|z\| |
| --- | ---: | ---: | ---: |
| Exact-law check, 3A → ∅ and 2A + B → C, `static_direct` and `next_reaction`, against 0.6.0 FSP | 272 (24 at t = 0) | **0** | 2.83 |
| W4 law gate (order 8) | 192 (32 at t = 0) | **0** | 1.94 |
| Regression run, all twelve models | 4,814 (940 at t = 0) | **0** | 3.88 |

The regression run's power sums and event counts equal those of the earlier
runs for every model. The first W4 gate run was stopped by a host restart
before it wrote a result and was rerun with the declared seeds.

Time each converted workload takes in the law gate (200,000 replicates,
every path kept, conversion and decoding included, Python backend, 2
threads, engine `1f66966`; informative, the speed gate is M4d):

| Model | Wall time | Events |
| --- | ---: | ---: |
| W1 | 13.2 s | |
| W2 | 67.0 s | |
| W3 | 158.4 s | |
| W4 | 532.8 s (558.1 s in the W4-only run) | 7.43e8 |
| W5 | 212.5 s (221.2 s in M4a) | 4.08e8 |
| W6 | 253.1 s | |
| W7 | 112.0 s (121.9 s in M4a) | 6.76e8 |
| enzyme-mechanism | 2.9 s | |
| immigration-death | 2.2 s | |
| michaelis-menten | 2.3 s | |
| repressilator | 58.3 s | |
| vilar-noise | 510.2 s | |
| all twelve | 32 min 6 s | |

Event rates in this run: W4 1.4e6 events per second (2 threads, decoding
included), W5 1.9e6, W7 6.0e6. W4 is the slowest per event: each of its 33
higher-order channels has all eight species as reactants, so every event
recomputes all 33 whole propensities (Michaelis–Menten expressions). CPU
time was not measured; M4d measures the speed gate against the 0.6.0
baseline without path keeping and decoding.

## Checks

On `1f66966` plus the `_rules_loop.py` fix (wheel `_native` `85673a83…`, the
fix copied into the environment), under the CPU lock.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 117 (5 new), `higher_order` 8 (new), `limits_v1` 5, `convert` 10, `contracts_v1` 7, CLI 14, I1 (`stored_requests` 2, `static_direct_defaults` 9), I3 (`generated_loop` 26: the 16 old fixtures byte-identical, 1 new; `generated_kernel` 3), all others |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass |
| Schemas against `crn-craft describe KIND` | `rules-model`, `rules-run`, `rules-result`, `rules-check` and the other 1.0 schemas byte-identical; `model-1.0` regenerated (`count`, three descriptions) |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| `check_dependencies.py` | pass; `-O` identical |
| `check_contracts.py` | 90 of 90; `-O` identical |
| `check_rules.py` | 47 of 47; `-O` identical |
| `check_mtdna.py` | passed; `-O` identical |
| `check_rules_python.py` | 24 checks pass; `-O` identical |
| `check_compiled_loop.py --require-cargo` | before the fix: exit 1 (the new higher-order case, see above); after: 54 rows (35 identical to the interpreter, 4 of them higher-order, 18 refusals, 1 cancellation); `-O` identical |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` | 22 blocks pass, before and after the fix; `-O` identical |
| `check_builder.py` | 18 tests, 0 failed; `-O` identical |
| `check_hybrid.py` | 0 failed; `-O` identical |
| `check_installed.py` (`python -I`, outside the tree) | 17 checks pass, before and after the fix; `-O` identical |
| `check_higher_order.py` | 272 comparisons, 0 failures; `-O` report identical |

## Limitations

- Event logs and observers list a higher-order rule's reactants once per
  pattern with count 1 (`LogEntry`, `Event::reactants`); the multiplicity
  is the model's `count`. `count: 2` on its own lists two entries, as the
  pair of two patterns.
- The combination count of a pattern with count k takes min(k, n − k)
  ratios per evaluation; a model with large counts and few molecules above
  them (for example k = 10⁶, n = 2 · 10⁶ at a rate small enough not to
  overflow) pays that per rate evaluation. The network engine had the same
  cost.
- `tau` treats a higher-order channel as critical within ten firings of
  losing a reactant (count × 10), a direct generalization of its 0.6.0 rule;
  no leap-accuracy study was made for higher-order models (M5c runs the tau
  gates).
- `binomial` does not run higher-order rules (its scope is single-parent
  unit births and deaths); generated rates kernels refuse them; `thinning`
  is not tested with them.
- The Python builder still refuses more than two reactant molecules
  (`["A", "A", "A"]`): accepting them needs `model/1.0` output and the
  `count` mapping in `builder.py`, more than a one-line change, so it stays
  with M6a.

## Deviations from the plan

- `Shape::Multi` holds `Vec<(class, u64)>` (one entry per pattern), not
  `Vec<(usize, u32)>`: patterns over one type add their counts, which can
  exceed `u32`, and the per-pattern list gives the event participants. The
  lowering merges them per species.
- The plan lists `binomial` among the methods that support `Multi`; it
  refuses such models through its existing scope rule, since every `Multi`
  channel takes at least three molecules.
- Generated rates kernels (`generate_rules_kernel`, not named in the plan)
  refuse higher-order rules with `CRN-R020` at `/method`, as they refuse
  whole-rule propensities; generated event loops support them.
- A `count` in a 0.x model is reported with the exact 0.6.0 unknown-field
  diagnostic. One difference remains: in a 0.x document with a `count` and
  another strict-parse error, 0.6.0 reported whichever came first in the
  document; now the other parse error is reported first.
- The converter keeps the 0.6.0 mapping for reactions of up to two molecules
  (`{"A": 2}` stays two patterns), so the converted bytes of every earlier
  model are unchanged; `count` appears only from three molecules on.
- `check_laws.py` prints per-model wall times on stderr (M4a).

## Migration notes

- `crn-craft.model/1.0` patterns accept `count` (1 to 2,147,483,647) and
  rules any number of patterns, for attribute-free types without `where`
  or `factor`; `rate` gives rate × Π C(n, k), `propensity` applies while
  every count is available. `count: 2` equals two identical patterns.
- Methods for higher-order rules: `static_direct`, `next_reaction`, `tau`,
  `hybrid` and generated event loops; `tree` (the default when `method` is
  absent), `mtdna_direct`, `binomial` and generated rates kernels refuse
  them (`CRN-R020` at `/method`).
- `crn-craft convert` and `crn_craft.convert` accept any stoichiometry
  (W4 converts); reactions of three or more molecules use `count` and
  labels `x`, `y` or `x1`, `x2`, ….
- New public function `crn_craft::rules::combination_rate` (also through
  the `crn-craft-rules` facade), which generated event loops call.
- `rules-model/0.1`–`0.4` documents, their schema and their diagnostics are
  unchanged.
