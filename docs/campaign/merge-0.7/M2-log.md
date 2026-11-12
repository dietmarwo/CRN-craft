# M2 log: converter from network models to rules models

Task M2 of [the plan](../../plan-0.7.0-engine-merge.md): convert a
`crn-craft.model/0.1` network model into a rules model, expose the
conversion, and pass the law gate of M0c on every convertible model.

## What changed

- New `crates/crn-craft/src/rules/convert.rs`: `from_network(bytes)`, the
  file variant `from_network_file(path)` (used by the CLI) and
  `CONVERTED_VERSION`. The function parses and validates the network model
  (its errors are returned unchanged; warnings do not block), refuses what
  the rules language cannot hold yet, builds a `rules::model::Model`, checks
  it with the rules compiler and returns pretty-printed JSON with a trailing
  newline. Key order follows the rules contract types; parameters and
  observables are sorted by name (`BTreeMap`). Output contract:
  `crn-craft.rules-model/0.3` in this commit (`model/1.0` from M3).
- Mapping (also in the module docs):

  | Network model | Rules model |
  | --- | --- |
  | `name`, `parameters` | the same names and values |
  | species `X`, `initial_count` n | type `X` without attributes; `initial` entry with count n (none for 0); observable `X: {"type": "X"}` |
  | reactants `{"A": 1}`, `{"A": 1, "B": 1}`, `{"A": 2}` | patterns `x: A`; `x: A, y: B` (sorted by name); `x: A, y: A` |
  | products `{"C": 2}` | product `C` listed twice (sorted by name) |
  | `mass_action` c | `"rate": "c"` |
  | `hill_activation` | `"propensity": "basal + vmax / (1 + (k / R)^n)"` |
  | `hill_repression` | `"propensity": "basal + vmax / (1 + (R / k)^n)"` |
  | `michaelis_menten` | `"propensity": "vmax * (S / (km + S))"` |

  A rate argument naming a parameter is written as that (possibly renamed)
  name; a literal as Rust's shortest round-trip decimal (`{:?}`, e.g.
  `0.5`, `1e-7`, `5e-324`, `1.7976931348623157e308`), which the rules
  lexer reads back to the same bits; `-0` is written `0.0`.
- New diagnostic `CRN-R016` (next free in the rules model range), added to
  `docs/diagnostics.md`. Refused, each at the network document's path:
  more than 64 species (`/species`), 1,024 parameters (`/parameters`) or
  1,024 reactions (`/reactions`); no reaction (`/reactions`; a rules model
  needs a rule); more than two reactant molecules
  (`/reactions/i/reactants`); more than 64 product molecules
  (`/reactions/i/products`); an initial count above 2^53
  (`/species/i/initial_count`); a population limit
  (`/species/i/population_limit`).
- `crn-craft convert MODEL` prints the converted bytes (exit 0) or
  `{"success": false, "diagnostics": [...]}` (exit 1). Python:
  `crn_craft.convert(model) -> bytes` through the new native function
  `convert_network` in `bindings/python/src/rules.rs`; refusals raise
  `CrnError`.
- Discovery: the native CLI lists `convert` and reports
  `crn-craft.describe/0.14` (was 0.13); the rules entry
  (`crn_craft::rules::describe`, `crn_craft.rules.describe()`) gains
  `convert` (`from`, `to`, command, Python name). `validation/check_contracts.py`
  and the CLI test expect the new command list and version.
- `validation/merge/check_laws.py`: the law gate, extendable by `--method`
  (any method object, so selectors in M4c), `--models`, `--expect-refused`;
  Python backend (`--python`, a worker per model) or CLI backend
  (`--binary`, standard library only). It reuses `merge_common` (manifest,
  exact power-sum moments, comparisons) and `references.power_sums`, maps
  species to observables by name, chooses `rules-run/0.6` or (after M3)
  `run/1.0` from the converted model's schema, checks every chunk for
  completion and its decoding against the engine's own statistics, and
  writes every row to a JSON report without timings.

## Facts checked in the code before writing it

- `rates.rs`: `CompiledRate::evaluate` returns Hill and Michaelis–Menten
  values as whole propensities and gates every law on reactant availability
  first (`state[s] < coefficient` gives 0); that is a rules `propensity`
  with patterns. Mass action is c × Π C(x, ν), the rules per-combination
  `rate` (n, n_A·n_B, n(n − 1)/2).
- `omega`: only `validate.rs` reads it (finite and positive); no rate law,
  kinetics, FSP or export code uses it. It is not kept; the description
  records a value other than 1.
- Names: `validate.rs` and the rules language use the same identifier
  syntax (1–64 ASCII letters, digits, `_`, not starting with a digit). The
  rules language also reserves its 14 function names and shares one
  namespace between parameters and observables; network species and
  parameters are already disjoint, so observables named after species never
  collide with parameters. `importers.py` checks ids with the same syntax
  (`_identifier`), writes BioNetGen species and reactions as `S_i` and
  `R_i`, `omega` 1.0, mass action only and no population limits; SBML ids
  and BioNetGen parameter names can still be function names. Policy: a
  reserved name gets `_`, or `_2`, `_3`, … if taken (species and parameters
  against both namespaces, reactions against reactions), in the order model
  name, species, parameters (sorted), reactions. `description` is the only
  free-text field of the rules model, so it lists every renaming
  (`species \`exp\` is \`exp_2\``).
- Population limits are refused rather than written as a `stop`: the
  network engine never applies the update and ends the replicate as
  `population_limit` (incomplete, statistics withheld); a rules stop applies
  the event and counts the replicate as `stopped`, which is complete.

## Rate laws: converted against the network engine

Tests in `rules::convert::tests` evaluate both engines on the same network
JSON: the network engine through `CompiledNetwork::propensity`, the
converted model through the rules compiler and the Static Direct channel
rate (`flat::Network::rate`).

| Law | Declared before the first run | Measured |
| --- | --- | --- |
| Hill, activation and repression; 6 (vmax, basal) × 5 k × 7 n from 0.5 to 100 × 14 counts from 0 to 2^53: 5,880 cases | relative ≤ 1e-14 · max(1, n); where (k/R)^n or (R/k)^n overflows (fraction below about 2^-1024), absolute ≤ vmax · 2^-1020; where the network engine refuses an underflowing positive propensity (`CRN-E042`), the new value ≤ basal + vmax · 2^-1020 | 5,764 cases relative, largest 3.5e-15 per unit of max(1, n); 116 at the floor; x = 0 and vmax = 0 bit-identical |
| Michaelis–Menten, 4 parameter sets × 9 counts | bit-identical (same operation order) | bit-identical |
| Mass action, source and first order, 6 coefficients × 12 counts | bit-identical | bit-identical |
| Mass action, A + B and 2A, same coefficients, counts up to 2^53 | relative ≤ 4 ε | within |
| 2A with c = 1, n = 0…50 | both equal n(n − 1)/2 exactly | equal |

The network model sets no upper limit on n (finite and positive), so the
grid stops at 100. The stable forms lose the network engine's
"combine with log(vmax)" branch only where the induced term is below
vmax · 2^-1020; there the converted propensity is basal (or 0) instead.

## Tests written first

`crates/crn-craft/tests/convert.rs` (10 tests: every example and W1, W2,
W3, W6 converts and passes `check_model`, deterministic and canonical
bytes; the documented mapping; literal round trips; renaming; `omega`;
each refusal; network diagnostics kept; a converted model runs with
`static_direct` near the Poisson law) was written before the converter and
run against a stub returning `Err([])`: 10 of 10 failed. The rate-law unit
tests above compare against the network engine and were written with the
implementation. `crates/crn-craft-cli/tests/cli.rs` gained
`convert_prints_the_converted_model_bytes_or_coded_refusals` and the new
discovery assertions.

## Law gate

Declared in `validation/merge/RESULTS.md` (M0c's law gate, plus the M2 run
declaration committed in `d5fb778` before the run). Run with the Python
backend, `{"kind": "static_direct"}` (rules-run/0.6, Ziggurat waits), 2
threads, all twelve models, `--expect-refused W4,W5,W7`.

| Result | Value |
| --- | --- |
| Converted and compared | 9 models (5 examples, W1, W2, W3, W6), 510 cells |
| Refused as declared | W4 (order 8), W5 (256 species), W7 (96 species), `CRN-R016` |
| Comparisons | 1,102 (1,020 against the reference, 82 against FSP laws), 204 at t = 0 |
| Failures | **0** (0 expected) |
| Largest \|z\| | 3.39 (W2); against FSP 2.28 (W1) |
| Cost | 18 min 36 s wall, 34.7 CPU min, 1.56e10 events |

Per-model rows are in RESULTS.md; every comparison is in
`validation/merge/laws/static_direct.json`. Smoke runs before the gate
(2,000 replicates; immigration-death, michaelis-menten, W1, W4,
repressilator) gave identical power sums and rows with the Python and the
CLI backend, and identical reports normally and under `python -O`; they are
not gate runs.

## Checks

On the M2 tree (`d5fb778` plus the result files), under the CPU lock.

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | pass |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | pass |
| `cargo test --release --workspace --locked` | pass: `crn-craft` lib 111 (3 new), `convert` 10 (new), CLI 13 (1 new), facade, generated kernels and loops (I3), `independent`, `static_direct_defaults` and `stored_requests` (I1), doctests 34 + 1 |
| `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked` | pass (the known output-filename collision warning of the CLI binary only) |
| `cargo build --release --workspace --locked`; wheel built from the root and installed | pass |
| I2 `check_manual.py --binary target/release/crn-craft --require-cargo` | pass, 22 blocks; `-O` identical |
| `check_contracts.py` | 86 of 86; `-O` identical |
| `check_rules.py` | 47 of 47; `-O` identical |
| `check_rules_python.py` | 22 checks pass; `-O` identical |
| `check_installed.py` (`python -I`, outside the tree) | 16 checks pass; `-O` identical |
| `check_dependencies.py` | pass; `-O` identical |
| `ruff format --check --exclude '*.md' .`, `ruff check .` | pass |

## Limitations

- Hill laws: below vmax · 2^-1020 the converted induced term is 0 (see
  above). The rules engine also accepts a positive propensity that
  underflows to 0, which the network engine refuses with `CRN-E042`.
- Literal rate arguments are written as decimals in expression strings;
  they round-trip exactly, but a reader sees `1e-7` where the network
  model had a JSON number.
- Initial counts are written as JSON floats (`50.0`), the rules contract's
  number type.
- Static Direct fails a replicate with `CRN-R031` when a positive channel is
  below 2^-53 of the total, so a converted model with a very weak channel
  needs `tree` (as any rules model does).

## Deviations from the plan

- Hill forms are the overflow-free `basal + vmax / (1 + (k / R)^n)` and
  `basal + vmax / (1 + (R / k)^n)` (coordinator's instruction), not the
  plan's `R^n / (k^n + R^n)` forms.
- More refusals than the plan lists (parameters, products, counts above
  2^53, population limits, no reaction): each is a rules limit or a
  semantic difference found while reading the code.
- `from_network_file` and the `convert` discovery entry are additions; the
  native discovery moved to `describe/0.14`, and `check_contracts.py` was
  updated for it (M7 rewrites that script).

## Migration notes

- New: `crn-craft convert MODEL`, `crn_craft.convert(model)`,
  `crn_craft::rules::convert::from_network`. Converted models are rules
  models; run them with `crn-craft rules` or `crn_craft.rules`.
- The native `crn-craft describe` reports `crn-craft.describe/0.14` and
  lists `convert`; the core and Python `describe()` stay at 0.10.
