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

## What changed

- `rules/model.rs`: new public limits `MAX_TYPES_1_0` (4,096),
  `MAX_PARAMETERS_1_0` (16,384), `MAX_OBSERVABLES_1_0` (4,096) and
  `MAX_RULES_1_0` (65,536), applied to `crn-craft.model/1.0` documents; the
  0.x constants (64, 1,024, 256, 1,024) and their `CRN-R015` messages are
  unchanged for `rules-model/0.1`–`0.4`. `MAX_OBSERVABLES` stays the tracker
  limit of every request version (256). `Compiled` records the model
  version. The class lookup of a reactant pattern searches only the classes
  of the pattern's type and compares without cloning expressions (it was a
  scan over every class with clones per rule).
- `rules/flat.rs` (static lowering): limits by version, 256 tuples and
  4,096 channels before 1.0, 4,096 tuples and 65,536 channels in 1.0, with
  the 0.x messages (`static backend supports up to N instances`,
  `static channel limit exceeded`). Index-based construction replaces every
  quadratic step:
  - product slots: a hash index by (type, attribute bits) instead of a scan
    over every species per product and channel;
  - class members and observable terms: species lists per type instead of
    scanning every species per class and per observable (and a
    `contains` over the type's observables per species);
  - channel changes: a sorted merge of the reactant and product entries
    instead of a dense delta vector of every species per channel;
  - species read by a rate: once per rule, from the observables the rate
    reads (new `Expr::read_list`), instead of a flag vector over every
    observable and every species per channel;
  - `affected`: the sorted union of the readers of each changed species
    (readers built once) instead of testing every channel against every
    channel;
  - `changed`: the sorted union of the projections holding each changed
    species instead of testing every projection's terms per channel.
  The tables are the same values in the same order; see the equivalence test
  below. Deferred ("cold") projections: the overflow bound per weight is now
  f64::MAX / (max(1024, terms) · 2^53), which is the old bound for every
  projection of up to 1,024 terms (all projections before 1.0).
- `rules/hybrid.rs`: the per-replicate flags "this channel's rate reads a
  continuous variable" come from the rule's read list instead of a flag
  vector over every observable per channel.
- The tree engine's tables (`types[t].classes`, `types[t].observables`,
  class and observable trees) were already built per type; its per-event
  work is linear in rules and observables. The generated event loop has no
  table sized by channels squared; its source is still capped at 1 MiB
  (dense coupling repeats rate code per dependency row). The generated rates
  kernel's `AFFECTED` table has the size of the `affected` lists.
- Converter: `CRN-R016` size refusals now use the 1.0 limits, which equal
  the network contract's own (4,096 species, 16,384 parameters, 65,536
  reactions), so a valid network model no longer reaches them; they stay as
  guards. W5 and W7 convert. Population limits and the other refusals are
  unchanged.
- Discovery and schemas: `describe()` gains `contracts_1_0.limits`; the
  `model-1.0` schema states the 1.0 limits in the `types`, `observables` and
  `rules` descriptions (written by `crn-craft describe model-1.0`); the
  `rules-*` schemas are byte-identical. `docs/diagnostics.md` lists both
  sets of limits under `CRN-R015` and the shorter `CRN-R016` list.
- `validation/merge/check_laws.py` prints each model's wall time on stderr
  (the report stays free of timings).

## Tests written first

`crates/crn-craft/tests/limits_v1.rs` (5 tests) before the implementation:
4 failed, 1 passed (the 0.x limits test, a regression guard that holds
before and after).

- `a_4096_species_8190_reaction_chain_compiles_in_under_two_seconds`:
  `prepare` with `static_direct` and `next_reaction`, then a short run
  (conservation of 1,000 molecules along the chain, 8,190 rule counters,
  4,096 observables). The 2 s assertion applies to optimized builds only.
- `model_1_0_accepts_its_limits_and_refuses_one_more`: 4,096 types and
  observables, 65,536 rules and 16,384 parameters check valid; one more of
  each gives `CRN-R015` at its path with the 1.0 limit in the message.
- `older_model_versions_keep_their_limits_and_diagnostics`: the 0.6.0
  messages for 65 types, 1,025 rules and parameters, 257 observables.
- `static_lowering_limits_follow_the_model_version`: 256/257 tuples and
  4,096/4,352 channels for `rules-model/0.4`, 4,096/4,097 tuples and
  65,536/69,632 channels for 1.0.
- `the_converter_writes_models_up_to_the_1_0_limits`: W5, W7 and a
  4,096-species network model convert to valid 1.0 models.

Added with the implementation: `rules::flat::tests::
index_based_lowering_equals_the_quadratic_construction` keeps the 0.6.0
construction as test code and compares changes, observable terms (bits),
`affected` and `changed` for every closed model of the repository (rules
models, generated-loop fixtures, converted examples and workloads). Updated:
`tests/convert.rs` (the old size refusals become "keeps its network
diagnostic" and "converts"; W5/W7 no longer refused) and the CLI convert
test (the refusal case is now a population limit).

## Measured

Compile-time target (declared above, 2 s): `prepare` of the chain, release
build, under the CPU lock, three runs.

| Method | Run 1 | Run 2 | Run 3 |
| --- | ---: | ---: | ---: |
| `static_direct` | 0.140 s | 0.122 s | 0.120 s |
| `next_reaction` | 0.123 s | 0.162 s | 0.126 s |

**Target met** (at most 0.162 s against 2 s).

Law gate for W5 and W7 (declared in RESULTS.md before the run): 3,520
comparisons, **0 failures**, largest |z| 3.88 (W5) and 2.91 (W7). Wall time
in the law gate (conversion, 200,000 replicates with every path kept,
decoding; Python backend, 2 threads): **W5 221.2 s, W7 121.9 s**; events
4.08e8 and 6.76e8, matching the 0.6.0 reference counts to 0.01%.

## Checks

On `66560f4` (wheel `_native` `32d62d00…`), under the CPU lock.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 112 (1 new), `limits_v1` 5 (new), `convert` 10, CLI 14, `contracts_v1` 7, I1 (`stored_requests` 2, `static_direct_defaults` 9), I3 (`generated_loop` 25, `generated_kernel` 3), all others |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass |
| Schemas against `crn-craft describe KIND` | 8 rules-family schemas: 7 byte-identical, `model-1.0` regenerated (3 descriptions) |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |
| `check_dependencies.py` | pass; `-O` identical |
| `check_contracts.py` | 90 of 90; `-O` identical |
| `check_rules.py` | 47 of 47; `-O` identical |
| `check_mtdna.py` | 151 checks pass; `-O` identical |
| `check_rules_python.py` | 24 checks pass; `-O` identical |
| `check_compiled_loop.py --require-cargo` (`CRN_CRAFT_SOURCE`, `CARGO_TARGET_DIR`) | 50 rows pass; `-O` identical (I3 by the Rust fixtures) |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` | 22 blocks pass; `-O` identical |
| `check_builder.py` | 18 tests, 0 failed; `-O` identical |
| `check_hybrid.py` | 0 failed; `-O` identical |
| `check_installed.py` (`python -I`, outside the tree) | 17 checks pass; `-O` identical |

The validation run was interrupted by a host restart after `check_builder`;
`check_hybrid` and `check_installed` were rerun on the same build afterwards.

## Limitations

- `tau`'s leap condition (`propensity_step`) evaluates every channel's rate
  after one firing of every enabled channel, O(channels²) per leap; it was
  left as it is (I1 covers a stored `tau` request). Large models should use
  the exact methods.
- `affected` lists are as large as the coupling: a rule whose rate reads an
  observable over every species makes every channel of that rule depend on
  every channel. Only the construction stopped being quadratic.
- Admission still reserves 256 + 8 · rules bytes per replicate record, so a
  65,536-rule model admits about 4,000 replicates per request within the
  2 GiB workspace.
- User documentation of the 1.0 limits (rules language reference, manual)
  is task M10b; this task updated rustdoc, the schema and
  `docs/diagnostics.md`.

## Deviations from the plan

- Besides `flat.rs` and the model compiler, `hybrid.rs` had a channels ×
  observables structure; it is index-based now.
- The static lowering limits (tuples and channels) are raised for 1.0 as
  well (4,096 and 65,536, the network engine's species and reaction limits);
  the plan names only the model limits, but a 4,096-species chain cannot be
  lowered within 256 tuples.

## Migration notes

- `crn-craft.model/1.0` accepts 4,096 types, 16,384 parameters, 4,096
  observables and 65,536 rules; its static methods accept 4,096 instance
  tuples and 65,536 channels. `rules-model/0.1`–`0.4` keep 64, 1,024, 256,
  1,024, 256 and 4,096, with the same diagnostics.
- `crn-craft convert` and `crn_craft.convert` accept every valid network
  model size (W5 and W7 now convert).
