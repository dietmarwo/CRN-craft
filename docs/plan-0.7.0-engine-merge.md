# CRN-craft 0.7.0: one engine, one API — execution plan

Date: 2026-11-11. Base: `main` at release 0.6.0 (`43105b2`), repository
`dietmarwo/CRN-craft`. Status: **plan only, nothing executed.**

This plan is written so that an agent can execute it task by task without
making design decisions. Every task names its files, its tests, its gate and
what to do when the gate fails. Where the plan says "read X to confirm", do
that before writing code; if the code disagrees with the plan, follow the
code, and record the difference in the task log.

## 1. Goal and decisions already taken

CRN-craft has two engines with separate model formats, Python APIs and docs:

- the **network engine** (crate `crn-craft`, Python `crn_craft.*`): fixed
  species lists, structured rate laws, several exact SSA variants, FSP,
  analysis, archives, import/export, xoshiro256++;
- the **rules engine** (crate `crn-craft-rules`, Python `crn_craft.rules`):
  attributed species, expression rates, builder, pandas/xarray results,
  hybrid ODEs, generated loops, ChaCha8.

Release 0.7.0 merges them: **one model format, one request and result
contract, one Rust crate, one Python API, one set of docs.**

Decisions by the owner (do not reopen them):

| # | Decision |
| --- | --- |
| D1 | Merge at the Rust level, not only in Python. |
| D2 | Byte-identical replay of network-engine results made before the merge is given up. Statistical equivalence is required instead. |
| D3 | Rules-engine requests (`rules-model` 0.1–0.4, `rules-run` 0.1–0.6) keep replaying byte-identically. This is the main guardrail of the whole project. |
| D4 | One random number generator: ChaCha8, as in the rules engine. xoshiro256++ disappears with the old loops. |
| D5 | One crate `crn-craft` and one Python module `crn_craft`. `crn_craft.rules` stays as an alias so existing user code keeps running. The crate `crn-craft-rules` becomes a thin re-export of `crn-craft`. |
| D6 | The merged release is 0.7.0. |

### Feature request included in 0.7.0 (FR1)

A user asked for two things after 0.6.0. Both are part of this plan:

| Part | Request | Task |
| --- | --- | --- |
| FR1a | `result.statistics` should also hold the median and other percentiles, such as 2.5% and 97.5%. | M4e |
| FR1b | `result.save(file_name)` writes the model details and the simulation results to one file: markdown or JSON, chosen by the file extension. | M6b |

They are placed in the merge because quantiles need a new request and result
version, and `run/1.0` and `result/1.0` are that version. Delivering them
earlier in a 0.6.x release would need `rules-run/0.7` and
`rules-result/0.4`; do that only if the owner asks for it.

## 2. Facts measured on main (0.6.0)

Sizes: network crate 12,130 lines, rules crate 14,205 lines (2,140 of them
tests). The rules crate already depends on the network crate for strict JSON,
diagnostics, hashing and the exponential and Poisson samplers. `crn-craft`
already depends on `rand_chacha`, `rand_core` and `libm`.

What the rules engine cannot do today, and the network engine can:

| Gap | Network engine | Rules engine | Closed in |
| --- | --- | --- | --- |
| Species count | 4,096 species, 65,536 reactions | 64 types, 1,024 rules, 256 observables | M4a |
| Reaction order | any stoichiometry (workload W4 has order 8) | at most two reactant patterns (`Shape::Source/Single/Pair/Cross`) | M4b |
| Event selection | DirectCR, SortingDirect, RSSACR, cached direct | one linear scan (`flat::direct_with`) | M4c |
| Exact distributions | `fsp.rs` | none | M5a |
| Network analysis | `analysis.rs`, `conservation.rs`, `cone.rs`, `linkage.rs` | none | M5b |
| Sealed archives | `archive.rs` | none | M5d |
| ODE right-hand side and Jacobians | `kinetics.rs` (four rate laws) | none | M5e |
| SBML/BioNetGen import, GillesPy2 export | `importers.py`, `export.rs` | none | M5f |

What makes the merge feasible:

- The network engine has four rate laws (`contracts::Rate`): `mass_action`,
  `hill_activation`, `hill_repression`, `michaelis_menten`. The rules
  expression language has `+ - * / ^`, comparisons, `if`, `min`, `max`, `exp`,
  `log`, `sqrt`, `abs`, `floor`, so all four can be written as expressions.
- Both use combinatorial mass action (`A + A` has propensity c·n(n−1)/2).
  The network model has one convention only (`stochastic_counts`).
- The rules engine already lowers closed models to a flat network
  (`flat::Network`: species, channels, `changes`, `affected`), and
  `rules.Model` already writes attribute-free models in that form.

Workload models in `validation/ssa-acceleration/models/`: W1 (1 species),
W2 (33 species, 66 reactions), W3 (repressilator), W4 (8 species, 66
reactions, order 8), W5 (256 species, 510 reactions), W6 (32 species, 992
reactions), W7 (96 species, 126 reactions). W4, W5 and W7 exceed today's
rules limits.

Files that use the network Python API and must be migrated: about 70 scripts
under `validation/`, `examples/tutorial/make_assets.py`,
`tutorials/optimization/{calibration,oscillators,sensitivity}.py`. Rust users
of the network crate API: `crates/crn-craft-cli`, `bindings/python`,
`integrations/jobs`, `integrations/diffsol`,
`validation/parameter-optimization/worker`.

## 3. Target architecture

```
crates/crn-craft/src/
  lib.rs                 describe, schema, check, sha256, root re-exports
  strict_json.rs diagnostic.rs sampling.rs ssa/ziggurat*.rs      kept
  rules/                 the engine (moved from crn-craft-rules/src)
    mod.rs model.rs expr.rs engine.rs tree.rs flat.rs flat_path.rs
    observe.rs leap.rs hybrid.rs continuous.rs thinning.rs atp.rs
    mtdna.rs mtdna_path.rs generated.rs generated_loop.rs loop_host.rs
    convert.rs           network model/0.1 -> model/1.0            (M2)
    select.rs            selectors on flat rates                   (M4c)
    flat_view.rs         read-only view of a lowered model         (M5a)
    dual.rs              dual-number expression evaluation         (M5e)
  fsp.rs analysis.rs conservation.rs cone.rs linkage.rs            ported to flat_view
  archive.rs export.rs kinetics.rs                                 ported
crates/crn-craft-rules/  facade: `pub use crn_craft::rules::*;`
crates/crn-craft-cli/    one command set
```

Deleted at the end (M8): `simulation.rs`, `ssa.rs`, `ssa/selectors.rs`,
`ssa/bounds.rs`, `next_reaction.rs`, `tau.rs` (except the Poisson sampler the
rules engine calls, which moves to `sampling.rs`), `rates.rs`, `network.rs`,
`parameters.rs`, `validate.rs`, `contracts.rs`, `codegen.rs`,
`analysis_v1.rs`.

### Unified contracts

| Contract | Meaning |
| --- | --- |
| `crn-craft.model/1.0` | The rules model (shape of `rules-model/0.4`) plus: higher limits (M4a) and reactant multiplicity (M4b). |
| `crn-craft.run/1.0` | The rules request (shape of `rules-run/0.6`) plus: `static_direct.selector` (M4c), `outputs.quantiles` (M4e), method `fsp` (M5a). |
| `crn-craft.result/1.0` | The rules result (shape of `rules-result/0.3`) plus `statistics.quantiles` (M4e). |
| `crn-craft.report/1.0` | The JSON file written by `Result.save` (M6b): model, request and results in one document. Written and read by Python only. |
| `crn-craft.fsp-result/1.0` | Fields of `fsp-result/0.1`; `species` are the flat species labels. |
| `crn-craft.analysis/1.0` | Fields of the current analysis report, computed from the flat view. |
| `crn-craft.archive/1.0` | Model, request, result and a SHA-256 manifest. |
| Still read | `rules-model/0.1–0.4`, `rules-run/0.1–0.6` (byte-identical results); network `model/0.1` through the converter. |
| No longer read | network `run/0.1`, `ssa-run`, `tau-request`, old `result`, `archive`. |

### Unified Python API

```python
import crn_craft as crn

m = crn.Model("sir")                       # the 0.6.0 builder, unchanged
m.add_species({"S": 999, "I": 1, "R": 0}).add_parameter({"beta": 2e-4, "gamma": 0.01})
m.add_reaction("beta * S * I", ["S", "I"], ["I", "I"])
m.add_reaction("gamma * I", ["I"], ["R"])

result = crn.ensemble(m, end=200.0, report_times=[0, 100, 200], replicates=200).simulate()
result.statistics["median"]                # also "quantiles": {0.025: ..., 0.5: ..., 0.975: ...}
result.save("sir.md")                      # or "sir.json"
dist = m.fsp(end=10.0, report_times=[10.0], bounds={"S": 999, "I": 1000, "R": 1000})
report = m.analyze()
m2 = crn.read_sbml("model.xml", convention="binomial")
```

- Top level: `Model`, `ensemble`, `prepare`, `load`, `compile_loop`,
  `loop_source`, `describe`, `check`, `convert`, `load_example`,
  `read_sbml`, `read_bngl`, `CrnError`, `CompiledLoopError`.
- `crn_craft.rules` exports the same objects under the same names as 0.6.0.
- Removed: `run_spec`, `simulate` (module level), `ssa_spec`, `simulate_ssa`,
  `prepare_ssa`, `prepare_parameters`, `ParameterTemplate`,
  `mass_action_kernel_source`, `tau`, module-level `fsp` and `analyze` taking
  network dictionaries.

### Unified command line

`crn-craft describe [KIND]`, `check MODEL [REQUEST]`,
`run MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]`,
`verify DIRECTORY`, `analyze MODEL`, `convert MODEL`,
`export MODEL gillespy2`. `rules-check` and `rules` stay as aliases of
`check` and `run`. Removed: `simulate`, `simulate-cached`, `ssa`, `tau`,
`fsp`, `analyze-v1`.

## 4. Rules for the executing agent

### 4.1 Invariants: check all of them before every PR

| # | Invariant | How to check |
| --- | --- | --- |
| I1 | Stored rules requests replay byte-identically. | `cargo test --release --workspace --locked` (the `stored_requests` and `static_direct_defaults` tests). Never edit files under a `tests/stored` directory. |
| I2 | The user manual runs unchanged. | `python validation/check_manual.py --binary target/release/crn-craft --require-cargo`, normally and under `-O`, identical reports. Until M10 do not edit `docs/mtdna-user-manual.md` or its expected outputs. |
| I2b | Code written for 0.6.0 against `crn_craft.rules` keeps running with identical output, also after M10. | From M10b on: `validation/compat/manual-0.6.0.md`, a frozen copy of the 0.6.0 user manual, is run by `check_manual.py --manual validation/compat/manual-0.6.0.md` in CI. Never edit that copy. |
| I3 | Generated loop sources are unchanged until M8c. | The checked-in fixtures in `tests/generated/` match (`cargo test`). |
| I4 | Every CI job is green before merging. | See 4.3. |
| I5 | Public items are documented. | `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked`; `#![warn(missing_docs)]` stays. |
| I6 | Statistical honesty. | Declare sizes, seeds, tests, alpha and speed targets in a RESULTS file and commit it **before** the first gate run. Keep failed rows. Never loosen a tolerance after seeing a result. |

If I1, I2 or I3 breaks, stop: the change is wrong, not the test.

### 4.2 Working rules

- One task = one branch = one PR. PRs target the integration branch
  `merge-0.7` (created in M0b), not `main`. Branch names: `merge/M<id>-<slug>`.
- Each fix or feature gets a failing-first test: write the test, see it fail,
  then write the code.
- American spelling. No personal names in code or docs.
- Diagnostic codes are stable. For a new diagnostic, take the next free
  number (`grep -rhoE "CRN-[A-Z][0-9]{3}" crates docs | sort -u | tail`) and
  add it to `docs/diagnostics.md`.
- Schemas are generated from Rust types (`crn-craft describe KIND`); update
  the files in `schemas/` from that output, never by hand. For fields whose
  schema must stay a bare `$ref`, and for unit enum variants, write docs as
  `#[doc = concat!("...")]`; a plain `///` there changes the schema shape.
- A change to a delivered contract's shape needs a new version. The 1.0
  contracts are not delivered until 0.7.0 is released, so they may change
  during this project.
- Never push tags, never publish, never force-push `main` or `merge-0.7`.
- Write a short log per task in `docs/campaign/merge-0.7/M<id>-log.md`: what
  changed, test results, limitations, any deviation from this plan.
- Do not edit `CHANGELOG.md` except in M10d; collect migration notes in the
  task logs.

### 4.3 Commands

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --workspace --locked
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --locked
python -m ruff format --check --exclude '*.md' . && python -m ruff check .
cargo build --release --workspace --locked
maturin build --release --locked --out dist   # from the repository root
python -m pip install --force-reinstall --no-deps dist/*.whl
# then every validation script named in .github/workflows/ci.yml,
# each normally and under `python -O`, with identical reports
```

Merge a PR only when
`gh api "repos/dietmarwo/CRN-craft/commits/<sha>/check-runs?per_page=100&filter=latest"`
shows every run `completed` with conclusion `success`, `skipped` or
`neutral`. Create PRs with `gh api repos/dietmarwo/CRN-craft/pulls --input
body.json`; merge with `gh api -X PUT .../pulls/N/merge -f
merge_method=merge -f sha=<sha>`. PR body sections: Commits, Verification,
Notes for review, Merge order. Use the commit and PR attribution lines your
session provides.

Known pitfalls: run ruff format with `--exclude '*.md'`; build the wheel
from the repository root (not with `-m bindings/python/Cargo.toml`); a macOS
runner shortage shows as jobs cancelled with zero steps, and the
`retry-runner-shortage` workflow reruns them.

### 4.4 When to stop and report to the owner

- An invariant (4.1) cannot be kept.
- A gate fails twice after honest fixes.
- A speed target of M4 is missed.
- The plan needs a design decision it does not contain.
- A task would delete something that a later task still needs.

Report what was tried, the measured numbers and the options. Do not work
around a gate.

## 5. Milestones

Order: M0 → M1 → M2 → M3 → M4 → M5 → M6 → M7 → M8 → M9 → M10 → M11. Tasks
inside M5 are independent of each other. M4e needs only M3, and M6b needs
only M4e and M6a, so the feature request FR1 does not wait for M5. Sizes:
S (under a day), M (one to two days), L (several days).

---

### M0. Preparation

**M0a. Stopgap docs (S).** Open a PR from the existing branch
`docs/two-engines` to `main` and merge it on green. It explains the two
engines in the README until 0.7.0 replaces it.

**M0b. Integration branch (S).** Create `merge-0.7` from `main` and push it.
Confirm that `.github/workflows/ci.yml` runs on pull requests against any
base (it has `on: pull_request` without a branch filter). `main` stays
releasable: a fix for 0.6.x is made on `main`, released from there, and
`main` is then merged into `merge-0.7`.

**M0c. Baseline and references (M).** The old engine will be deleted, so its
speed and its laws must be recorded while it exists.

- Create `validation/merge/RESULTS.md` first, with the declarations below,
  and commit it before running anything.
- Write `validation/merge/baseline.py`. For each workload W1–W7 and each
  model in `models/examples/`, and for each algorithm (`direct`,
  `cached_direct`, `direct_cr`, `sorting_direct`, `rssa_cr`,
  `next_reaction`) with `fast` on and off: run the 0.6.0 network engine
  three times with fixed seeds and record the median wall time, the event
  count and the sampler string. Reuse the horizons and replicate counts of
  `validation/acceleration_workloads.py`. Output:
  `validation/merge/baseline-0.6.0.json`.
- Write `validation/merge/references.py`. It stores, under
  `validation/merge/reference/`:
  - for every model where FSP runs within its limits: the full FSP result
    (states, probabilities, moments) at the workload's report times;
  - for every model: mean and variance of each species at each report time
    from 200,000 replicates of `direct` (fixed seed), with standard errors;
  - the analysis report of every model (`crn-craft analyze`);
  - the ODE right-hand side and Jacobians from `kinetics` at the initial
    state for every model that `kinetics` accepts.
- Declarations in `RESULTS.md`:
  - **Law gate:** for each species and report time, the new engine's mean
    from 200,000 replicates lies within 5 standard errors (combined) of the
    reference mean, and within 5 standard errors of the FSP mean where an
    FSP reference exists. Variances: within 5 standard errors, using the
    fourth-moment estimate of the variance's standard error. With about
    N comparisons in total, at most 0 failures are expected; any failure is
    investigated, not rerun with another seed.
  - **Speed gate (M4):** for each workload, the unified engine's best
    method takes at most 1.25 times the time of the 0.6.0 engine's best
    method on the same host, same thread count, medians of three.
  - **Selector rule (M4c):** composition–rejection is ported. Sorting and
    RSSA-CR are ported only if, in `baseline-0.6.0.json`, one of them beats
    `direct_cr` by at least 1.3× on at least one workload. Otherwise they
    are dropped and the log says so.

Gate: the two JSON outputs exist, are deterministic for fixed seeds (run
twice, compare all fields except times), and are committed.

---

### M1. One crate, no behavior change (M)

Pure move. No function body changes.

1. `git mv crates/crn-craft-rules/src/*.rs crates/crn-craft/src/rules/`,
   with `lib.rs` renamed to `mod.rs`. Add `pub mod rules;` to
   `crates/crn-craft/src/lib.rs`.
2. Inside the moved files, replace `crate::` by `crate::rules::` (or
   `super::`) for paths that pointed at the old crate root, and
   `crn_craft::` by `crate::`. The moved files' crate-level attributes and
   crate docs become module-level.
3. `crates/crn-craft-rules/src/lib.rs` becomes a facade:
   `pub use crn_craft::rules::*;` plus `pub use crn_craft::rules::{atp,
   continuous, model, observe, libm};` so every path that compiled against
   0.6.0 still compiles. Keep its `Cargo.toml` name and version.
4. Leave `crates/crn-craft-rules/tests/` where it is; the tests now run
   through the facade and prove it is complete.
5. The generator keeps emitting `crn_craft_rules::` paths (I3). The Python
   loop driver keeps depending on `crn-craft-rules`.
6. Update the dependency policy read by `validation/check_dependencies.py`
   only if the gate reports a difference.

Gate: all commands of 4.3; I1, I2, I3; `python validation/check_compiled_loop.py
--require-cargo` with `CRN_CRAFT_SOURCE` set to the checkout;
`git diff --stat` shows renames and path edits only.

---

### M2. Converter: network model → rules model (M)

New file `crates/crn-craft/src/rules/convert.rs`:

```rust
/// Converts a `crn-craft.model/0.1` document into a rules model document.
pub fn from_network(model_bytes: &[u8]) -> Result<Vec<u8>, Vec<Diagnostic>>
```

Before coding, read `contracts.rs` (`Model`, `Species`, `Reaction`,
`Scalar`, `Rate`) and `rates.rs` to confirm each formula below, in
particular whether Hill and Michaelis–Menten propensities are multiplied by
reactant counts (the field docs say they are whole propensities in events
per second).

| Network model | Rules model |
| --- | --- |
| `name`, `parameters` | the same; a `Scalar::Literal` is written inline as a number |
| species `X` with `initial_count` n | type `X` without attributes; `initial` entry with count n (omitted when 0); observable `X: {"type": "X"}` |
| reaction reactants `{"A": 1, "B": 1}` | patterns `x: A`, `y: B` |
| reactants `{"A": 2}` | two patterns of type `A` (a pair) |
| `mass_action` with coefficient c | `"rate": "<c>"` (per-combination rate; the engine multiplies by the number of reactant combinations) |
| `hill_activation` | `"propensity": "basal + vmax * R^n / (k^n + R^n)"` with the regulator's observable `R` |
| `hill_repression` | `"propensity": "basal + vmax * k^n / (k^n + R^n)"` |
| `michaelis_menten` | `"propensity": "vmax * S / (km + S)"` |
| products `{"C": 2}` | product `C` listed twice |
| `omega` | checked to be 1 or unused in propensities; read `rates.rs` to confirm it does not rescale, as the field doc says |

- Output schema: `crn-craft.rules-model/0.3` in M2 (changed to `model/1.0`
  in M3). The output is canonical: fixed key order, parameters sorted.
- Models outside today's rules limits (more than two reactant molecules,
  more than 64 species, more than 1,024 reactions) are refused with a new
  coded diagnostic until M4.
- Expose it as `crn-craft convert MODEL` (prints the converted model) and
  as `crn_craft.convert(model) -> bytes`.

Tests first (`crates/crn-craft/tests/convert.rs`): each example model
converts; the converted document passes `rules::check_model`; conversion is
deterministic; each refusal has a test.

Gate: 4.3, I1–I3, and the **law gate** of M0c for every convertible model
(all of `models/examples/`, W1, W2, W3, W6), run with
`method: static_direct`. Script: `validation/merge/check_laws.py`, which
later milestones extend. Record the rows in `validation/merge/RESULTS.md`.

---

### M3. Unified contract names (M)

Add the 1.0 contracts as new versions of the rules contracts, with no new
features yet.

1. `ModelVersion`, `RunVersion`, `ResultVersion`, `CheckVersion` in
   `rules/model.rs` and `rules/mod.rs` each gain a variant for
   `crn-craft.model/1.0`, `crn-craft.run/1.0`, `crn-craft.result/1.0`,
   `crn-craft.check/1.0`. A `model/1.0` document means what the same
   document means as `rules-model/0.4`; `run/1.0` means `rules-run/0.6`
   (including its defaults: Ziggurat waits and "auto" sparse trackers); a
   `run/1.0` request produces `result/1.0`, which has the shape of
   `rules-result/0.3`.
2. A `run/1.0` request accepts `model/1.0` and every `rules-model` version;
   a `rules-run` request does not accept `model/1.0` (coded refusal), so old
   requests cannot reach new features.
3. `describe` lists the new versions. Add `schemas/model-1.0.schema.json`,
   `run-1.0.schema.json`, `result-1.0.schema.json` generated from the types.
4. The converter (M2) now writes `crn-craft.model/1.0`.
5. `rules.Model.to_dict()` in Python keeps writing `rules-model/0.3` or
   `0.4` (I2 depends on it). It switches in M6.

Tests first: for each stored request, the same request rewritten as
`run/1.0` gives a result that is identical except for the schema strings and
the hashes that cover them (write a helper that compares after removing
exactly those fields).

Gate: 4.3, I1–I3, the new equivalence test.

---

### M4. Extend the engine to the network engine's scope

All new behavior is reachable only through `model/1.0` and `run/1.0`.

**M4a. Limits (M).** For `model/1.0`: up to 4,096 types, 65,536 rules,
4,096 observables, 16,384 parameters. Old model versions keep their limits
and diagnostics. Before raising them, read `flat::Network::compile` and
`tree.rs` for structures that are quadratic in types or rules (for example
`affected: Vec<Vec<usize>>` built by comparing every channel with every
channel) and replace them with index-based construction where needed; add a
test that compiles a 4,096-species, 8,190-reaction chain in under 2 s.
Converter limits follow. Gate: 4.3, I1–I3, law gate for W5 and W7.

**M4b. Higher-order reactions (L).** A pattern gains an optional integer
field `count` (default 1, 1 ≤ count ≤ 2,147,483,647), and a rule may have
any number of patterns, under these conditions: `model/1.0` only; every
pattern with `count > 1`, and every pattern in a rule with more than two
patterns, is over an attribute-free type without `where` or `factor`.
Otherwise a coded refusal.

- Semantics with `rate`: propensity = rate × ∏ C(nᵢ, countᵢ) over the
  distinct types (two patterns of the same type add their counts). With
  `propensity`: the stated value while every type has at least its count,
  else 0.
- Implementation: extend `Shape` with a variant `Multi(Vec<(usize, u32)>)`
  used only when the rule is not expressible as `Source`, `Single`, `Pair`
  or `Cross`; existing shapes keep their code paths and arithmetic (I1).
  `flat.rs`: `ChannelRate` gets a general path computing the combination
  count with the overflow-safe scaled product of `rates.rs` (read
  `rates.rs` lines around the "binomial(n,k) may overflow" comment and
  port that routine and its tests).
- Methods: `static_direct`, `next_reaction`, `tau`, `binomial`, `hybrid`
  and generated loops support `Multi`. `tree` and `mtdna_direct` refuse
  with a coded diagnostic at `/method`.
- The converter maps any stoichiometry using `count`.

Tests first: exact-law tests against FSP references for `3A → ∅` and
`2A + B → C` (small state spaces; add them to `validation/merge/`); the
overflow tests ported from `rates.rs`; refusals. Gate: 4.3, I1–I3, law gate
for W4.

**M4c. Selectors (L).** `run/1.0` gains
`method.static_direct.selector`: `"linear"`, `"composition_rejection"`,
`"auto"` (and `"sorting"`, `"rssa_cr"` only if the selector rule of M0c
kept them). Absent means `"auto"` in `run/1.0`; old run versions always
use the linear scan (I1).

- New file `rules/select.rs`. Port the composition–rejection structure and
  its tests from `ssa/selectors.rs`: groups by binary exponent, incremental
  update of the channels in `affected`, rejection inside a group. Random
  numbers come from the rules `Words` stream (ChaCha8).
- `"auto"`: linear below a channel-count threshold, composition–rejection
  at or above it. Determine the threshold from a benchmark over channel
  counts 8…4,096 (declare the benchmark first; store it in
  `validation/merge/selector-threshold.json`), then fix it as a constant.
- The sampler string names the selector, so results with different
  selectors are distinguishable.
- Generated loops: support `linear` only; a request with another selector
  and `compiled=` is refused with `CRN-R020`.

Tests first: the ported unit tests; exact-law tests per selector (reuse the
structure of `validation/check_alternative_ssa.py`: immigration–death and
dimerization laws against exact distributions); a test that `linear` under
`run/1.0` equals the `rules-run/0.6` result apart from schema fields.

Gate: 4.3, I1–I3, law gate for W1–W7 with each selector.

**M4d. Speed gate (M).** Run `validation/merge/compare_speed.py` (write
it): unified engine against `baseline-0.6.0.json` on W1–W7 and the example
models, same host, medians of three. Apply the speed gate of M0c. If it
fails for a workload, profile that workload, fix, rerun once; if it still
fails, stop and report (4.4) with the numbers.

**M4e. Quantile statistics (M) — FR1a.** Median and percentiles of every
observable at every report time, computed in the engine.

Request, `run/1.0` only: `outputs.quantiles`, an optional list of at most
16 strictly increasing probabilities with 0 < p < 1. Absent or empty means
none. Old run versions refuse the field (coded diagnostic).

Result, `result/1.0`: `statistics.quantiles`, a list of
`{"p": 0.025, "values": [[...], ...]}` in request order, where `values` is
indexed like `statistics.mean` (report time, then observable, continuous
variables included). It is `null` whenever `statistics.mean` is `null`
(read the existing rule for withholding statistics in `rules/mod.rs` and
apply the same one).

Definition: the linear-interpolation quantile, NumPy's default
(`numpy.quantile(x, p)`, also called type 7): with the n values sorted,
h = (n − 1)·p, result = x[⌊h⌋] + (h − ⌊h⌋)·(x[⌊h⌋+1] − x[⌊h⌋]). State this
definition in the field's documentation.

Implementation:

- Collect the reported values of every replicate in one buffer of
  replicates × report times × observables `f64`, written at the replicate's
  index, so the content does not depend on the number of threads. After the
  ensemble, sort each (time, observable) column with `f64::total_cmp` and
  evaluate the formula.
- Memory: the buffer counts against the existing workspace limit
  (`MAX_WORKSPACE_BYTES`; read how `prepare` accounts for kept paths and
  add the buffer the same way). A request over the limit is refused at
  preparation with a coded diagnostic at `/outputs/quantiles` whose hint
  says to request fewer report times, fewer replicates or no quantiles.
- Read how statistics are accumulated for compiled loops and for
  `bind_parameters` candidates and make quantiles work there too;
  `loop.verify(run)` must still report byte-identical results. If a path
  cannot support them, refuse with `CRN-R020` and a test, and record it in
  the log.

Python (`rules.py`): `ensemble(..., quantiles=...)` and
`prepare(..., quantiles=...)`. `result.statistics` gains
`"quantiles"`, a dict from probability to a DataFrame shaped like
`statistics["mean"]`, and `"median"`, the same DataFrame as
`quantiles[0.5]`, present when 0.5 was requested. Default: through the
top-level API (`crn.ensemble`) the default is `(0.025, 0.5, 0.975)`;
through `crn_craft.rules` the default is none, so 0.6.0 code produces
0.6.0 results (I2). `quantiles=()` switches them off.

Tests first: (a) Rust: for a small ensemble with every path kept, the
quantiles equal the formula applied to the kept paths, bit for bit;
(b) results identical for 1, 2 and 4 threads; (c) the refusals;
(d) Python `check_rules_python.py`: `statistics["quantiles"][p]` equals
`np.quantile(result.paths.values, p, axis=0)` exactly, and `"median"`
equals `quantiles[0.5]`; (e) an exact-law check: for immigration–death at
a late report time, the median and the 2.5% and 97.5% quantiles of 100,000
replicates each lie within one count of the Poisson quantiles (declare in
`validation/merge/RESULTS.md` first).

Gate: 4.3, I1–I3, the new tests, and a measurement of the added time and
memory for 4,000 cells × 37 reports × 6 observables (record it; no target).

---

### M5. Port the features only the network engine has

Common basis, built in M5a and reused by M5b, M5e and M5f:

```rust
// crates/crn-craft/src/rules/flat_view.rs
/// A read-only view of a closed model lowered to a flat network.
pub struct FlatView { /* built from a Prepared model */ }
impl FlatView {
    pub fn species(&self) -> &[String];            // labels, e.g. "A" or "mtDNA{wt=1}"
    pub fn initial(&self) -> Vec<u64>;
    pub fn channels(&self) -> usize;
    pub fn changes(&self, channel: usize) -> &[(usize, i64)];
    pub fn reactants(&self, channel: usize) -> Vec<(usize, u32)>;
    pub fn propensity(&self, channel: usize, counts: &[u64]) -> Result<f64, Diagnostic>;
}
```

A model that is not closed (new attribute values appear at run time) has no
flat view; FSP, analysis, ODE export and GillesPy2 export refuse it with a
coded diagnostic.

**M5a. Finite state projection (L).** `run/1.0` method
`{"kind": "fsp", "bounds": {"<species label>": max, ...}, "max_states": N,
"tolerance": t}`; result `crn-craft.fsp-result/1.0`. Rewrite `fsp.rs` to
take a `FlatView` instead of a `CompiledNetwork`; the numerical method
(uniformization, error accounting) stays as it is. CLI: `crn-craft run`
dispatches on the method. Tests first: port `crates/crn-craft/tests/fsp.rs`
to converted models. Gate: 4.3, I1–I3; probabilities agree with the M0c FSP
references to 1e-10 absolute (rates are now evaluated by expressions, so
bits may differ); `validation/check_fsp.py` migrated to the new command and
passing, including the DSMTS cases.

**M5b. Analysis (M).** `analysis.rs`, `conservation.rs`, `cone.rs`,
`linkage.rs` take a `FlatView`. `crn-craft analyze MODEL` accepts
`model/1.0`, rules models and (through the converter) `model/0.1`. Drop
`analysis_v1.rs` and the `analyze-v1` command. Gate: reports equal the M0c
analysis references field by field, apart from the schema string;
`check_network.py` and `check_structure.py` migrated and passing.

**M5c. Tau leaping (M).** Keep the rules methods `tau` (with `epsilon`) and
`binomial`. Run the gates of `validation/check_tau.py` against the rules
`tau` on converted models (declare first; same tolerances as the script
has today). If a gate fails, stop and report; do not port the old leap code
without asking. Move `tau::poisson_with` (used by the rules engine) to
`sampling.rs` unchanged. Gate: `check_tau.py` migrated and passing.

**M5d. Archives (M).** `crn-craft run MODEL REQUEST --archive NEW_DIRECTORY`
writes the model bytes, the request bytes, the result and a manifest with
the SHA-256 of each file (`crn-craft.archive/1.0`); `crn-craft verify
DIRECTORY` checks the manifest, the hashes inside the result, and that the
result parses. Reuse the manifest and path-safety code of `archive.rs`
(new directory only, portable relative paths, no overwrite). Tests first:
port `crates/crn-craft/tests/archives.rs`, including the tamper cases.
Gate: 4.3, I1–I3.

**M5e. ODE right-hand side and Jacobians (L).** `kinetics.rs` takes a
`FlatView`. New file `rules/dual.rs`: evaluate an `Expr` with dual numbers
(value plus a vector of partial derivatives) for every operator and
function of `expr.rs`; non-differentiable points (`floor`, comparisons,
`if` conditions) use the derivative of the selected branch and 0 for
`floor`. The right-hand side is Σ changes × propensity with counts as
reals; combination counts use the macroscopic form already documented in
`kinetics.rs` (`y^order / order!`). Jacobians with respect to state and to
named parameters come from the dual parts. Tests first: port
`crates/crn-craft/tests/kinetics.rs`; add finite-difference checks of the
Jacobians for a model with `^`, `/`, `exp` and `max`. Gate: values agree
with the M0c kinetics references to 1e-12 relative; `integrations/diffsol`
builds and its tests pass against the new interface; `check_ode.py` passes.

**M5f. Import and export (M).** `python/crn_craft/importers.py` returns
`Model` objects (build the network dictionary as today, then convert).
`export.rs` (GillesPy2) takes a `FlatView`: a `rate` that is a constant
expression exports as mass action, anything else as a custom propensity
string. Gate: `check_importers.py` and `check_export.py` migrated and
passing.

**M5g. Parameters and workers (M).** The network `ParameterTemplate` is
replaced by the rules engine's parameter binding (what Python calls
`bind_parameters`; read `rules/mod.rs` and `bindings/python/src/rules.rs`
for the Rust entry point). Move `integrations/jobs` and
`validation/parameter-optimization/worker` to it, with converted models.
Gate: `check_jobs.py`, `check_parameter_interface.py`,
`validate_parameters.py` passing; each worker's lockfile builds with
`--locked`.

---

### M6. One Python API (L)

1. `python/crn_craft/__init__.py` exports the unified names of section 3.
   `python/crn_craft/rules.py` and `builder.py` keep their code;
   `crn_craft.rules` keeps exporting everything it exported in 0.6.0 (I2).
2. `Model.to_dict()` writes `crn-craft.model/1.0`, and `ensemble`/`prepare`
   write `run/1.0`, **only when called through the top-level names**, or
   when the model uses a 1.0 feature (`count`, a selector, `fsp`). Called
   through `crn_craft.rules` without such a feature, they write the 0.6.0
   schemas, so the manual's printed hashes and outputs stay the same (I2).
   Implement this with one flag on the model and prepared-run objects.
3. New `Model` methods: `fsp(...)` returning an `FspResult` (attributes:
   `states` DataFrame, `probabilities` DataArray with dimensions
   `("time", "state")`, `mean` DataFrame, `truncation` DataFrame of the
   error fields, `metadata`), `analyze()` returning the report dictionary,
   `to_gillespy2()`. `add_reaction` accepts repeated names in the reactant
   list for higher order (`["A", "A", "A"]`) and maps them to `count`.
4. `crn_craft.load_example(name)` returns a `Model` built from the
   converted packaged examples (store the converted files in
   `python/crn_craft/models/`).
5. Remove the functions listed in section 3 and the PyO3 functions behind
   them (`bindings/python/src/lib.rs`).

(The numbered list above is task **M6a**.)

Tests first: `validation/check_python.py` rewritten for the unified API,
keeping every law and error check it has (same tolerances). Gate: 4.3,
I1–I3, `check_python.py`, `check_rules_python.py`, `check_builder.py`,
`check_hybrid.py`, `check_installed.py` (add: FSP of immigration–death
against the Poisson law; `analyze()` of the enzyme model; a third-order
reaction).

**M6b. `Result.save(path)` (M) — FR1b.** One file with the model details
and the simulation results.

```python
result.save("sir.md")     # markdown report
result.save("sir.json")   # the same content as JSON
```

Signature: `Result.save(self, path, *, paths=True)`. The extension decides
the format (`.md`, `.markdown` or `.json`, case-insensitive); any other
extension raises `ValueError` naming the supported ones. `paths=False`
leaves the single trajectories out. An existing file is replaced; write to
a temporary file in the same directory and rename it, so a failed write
leaves no partial file. The method returns the path as a `pathlib.Path`.
The output is deterministic: the same result gives the same bytes (no
timestamp, no host name, no timings).

Prerequisite: a `Result` must know its model and request. Today
`PreparedEnsemble` keeps the exact `(model, request)` bytes in
`_documents`; pass them to the `Result` in `simulate()` and expose them as
`result.model` and `result.request` (parsed dictionaries). Read how
`bind_parameters` candidates and compiled runs create their results and
make sure the documents there are the ones that actually ran (the bound
parameter values, not the template's).

Markdown layout (follow it exactly; `<...>` are placeholders):

```text
# Model

Name: <model name>

## Reactions
S + I -> I + I, beta * S * I
I -> R, gamma * I

## Start values
S = 999
I = 1
R = 0

## Parameter values
beta = 0.0002
gamma = 0.01

## Stop condition
I == 0

# Simulation
seed = 1
end = 250
replicates = 10
report_times = 0, 50, 100, 150, 200, 250
method = static_direct
sampler = <sampler string>
crn-craft = <package version>
model_sha256 = <hash>
run_sha256 = <hash>

# Results

## Terminations
completed = 10

## rule_events (mean per replicate)
infection = 998.3
recovery = 820.6

## Single trajectories for S, I, R
S000: 999, 998, ...
I000: 1, 2, ...
R000: 0, 0, ...
S001: 999, 999, ...

## Summary statistics
samples = 10

### mean
S: 999, 998.5, ...
I: 1, ...
R: 0, ...

### mean_standard_error
### sample_variance
### median
### 2.5% percentile
### 97.5% percentile
```

Rules for the content:

- **Reactions.** One line per rule in declaration order:
  `<reactants> -> <products>, <rate>`. Reactants and products are type
  names joined by ` + `; an empty side is written `0`. The rate is the
  rule's `propensity` string if it has one, else its `rate` string followed
  by ` (mass action)`. A rule that uses attributes, `where`, `factor` or
  `copy` cannot be written this way: print its name, a colon and the rule
  as compact JSON. Models with attributes also get a `## Species types`
  section listing each type with its attributes.
- **Start values.** One line per type in declaration order, including zero
  counts; for attributed initial entries, one line per entry with its
  attributes in braces.
- **Parameter values.** In sorted order. Omit the section when there are
  none. `## Stop condition` and `## Continuous variables` (name, initial,
  derivative, bounds, `independent`) appear only when the model has them.
- **Simulation.** `seed`, `end`, `replicates`, `report_times`, `method`
  with its options, `sampler`, package version and the two hashes.
- **Terminations.** Counts of each termination that occurred.
- **rule_events.** The mean over replicates of each rule's event count.
- **Single trajectories.** Only when paths were kept and `paths=True`. For
  each kept replicate in ascending id, one line per observable (order
  below), labeled `<observable><id>` with the id zero-padded to the width of
  the largest kept id (at least three digits). The values are that
  replicate's series over the report times; values after a budget stop are
  written `nan`. The heading lists the observables.
- **Summary statistics.** Omit the section, with one line
  `statistics withheld (<reason from the result>)`, when
  `result.statistics` is `None`. Otherwise one subsection per statistic the
  result holds, in this order: `mean`, `mean_standard_error`,
  `sample_variance`, then the quantiles in ascending order. The quantile
  0.5 is titled `### median`; any other p is titled
  `### <100·p>% percentile` with the shortest exact decimal (2.5, 97.5,
  25). Each subsection has one line per observable.
- **Order of observables.** In the trajectories and in every statistics
  subsection, observables appear in the order the model declares its
  species types when every observable is named after a type (builder
  models: `S`, `I`, `R` as declared), followed by continuous variables.
  Otherwise they appear in result order (sorted names).
- **Numbers.** Shortest round-trip form (`repr` of a Python float), except
  that a float with an integral value below 2^53 is written without `.0`.
  List separators are `, `.

JSON layout (`crn-craft.report/1.0`), same information, machine-readable:

```json
{
  "schema": "crn-craft.report/1.0",
  "crn_craft_version": "0.7.0",
  "model": { "...the model document that ran..." },
  "request": { "...the request document that ran..." },
  "result": {
    "model_sha256": "...", "run_sha256": "...", "sampler": "...",
    "report_times": [0, 50],
    "observables": ["I", "R", "S"],
    "terminations": {"completed": 10},
    "rule_events_mean": {"infection": 998.3, "recovery": 820.6},
    "statistics": {
      "samples": 10,
      "mean": {"I": [1, 2.5], "R": [0, 0.1], "S": [999, 997.4]},
      "mean_standard_error": {}, "sample_variance": {},
      "quantiles": [{"p": 0.025, "values": {"I": [], "R": [], "S": []}}]
    },
    "paths": {"replicates": [0, 1], "values": {"I": [[1, 2], [1, 1]]}}
  }
}
```

`statistics` is `null` when withheld, `paths` is `null` when no path was
kept or `paths=False`. Non-finite values are written as the strings
`"nan"`, `"inf"`, `"-inf"`, so the file is strict JSON. Keys are in the
order shown; the file ends with a newline. Add
`crn_craft.load_report(path)` returning the parsed dictionary after
checking the schema string; no schema file is generated for this
Python-level format, and the manual documents it.

Relation to archives (M5d): an archive is the sealed, verifiable record of
a run for replay; `save` is a readable report. Say so in the docstring.

Tests first:

- `validation/check_installed.py`: the SIR model of the request (10
  replicates, seed 1, every path kept, quantiles 0.025, 0.5, 0.975) saved
  as markdown equals a checked-in expected file
  `validation/fixtures/sir-report.md` byte for byte, and saved twice gives
  identical bytes.
- `validation/check_rules_python.py`: the JSON report round-trips: every
  number in it equals the corresponding value of `result.records`,
  `rule_events`, `statistics` and `paths`; `paths=False` gives
  `"paths": null`; withheld statistics give `null` and the markdown line;
  an unsupported extension raises `ValueError`; an attributed model (the
  three-class mtDNA model) and a model with a continuous variable save
  without error and contain their extra sections.
- A failed write (target directory missing) raises `OSError` and leaves no
  temporary file.

Docs (in M10b): a manual section "Saving a report" with the SIR example
and its markdown output, the reference table entry, and one line in
`docs/using-from-r.md` (`result$save("sir.md")`).

Gate: 4.3, I1–I3, the new tests normally and under `-O`.

---

### M7. One command line (M)

Implement the command set of section 3 in `crates/crn-craft-cli/src/main.rs`
and `python/crn_craft/__main__.py`. `rules` and `rules-check` are aliases.
`crn-craft describe` lists the 1.0 contracts first and the still-read old
ones under `legacy`. Gate: `check_contracts.py` rewritten for the new
commands (keep every strictness check: unknown fields, duplicate keys,
size limits, exit codes), passing normally and under `-O`.

---

### M8. Delete the old engine and migrate validation

**M8a. Validation scripts (L).** Migrate each script to the unified
commands and API, keeping its statistical rules. Do them one PR per group:

| Group | Scripts |
| --- | --- |
| Exact laws | `check_ssa.py`, `check_cme.py`, `check_alternative_ssa.py`, `check_alternative_contracts.py`, `mutation_check.py`, `reference_rng.py` (replace the xoshiro reference by the ChaCha8 stream reference) |
| Features | `check_fsp.py`, `check_tau.py`, `check_network.py`, `check_structure.py`, `check_export.py`, `check_importers.py`, `check_ode.py` (if not already done in M5) |
| Python and docs | `check_python.py`, `check_tutorial.py` (with M10a), `check_installed.py` |
| Benchmarks | `benchmark_*.py`, `measure_*.py`, `profile_ssa_acceleration.py`, `acceleration_workloads.py`: port those still meaningful; delete those that measure removed algorithms and say so in the log |
| Comparisons | `compare_gillespy2.py`, `compare_references.py`, `compare_tau_samplers.py`: port. `compare_julia*.py` and the gap-10 scripts: port, but re-measure only if Julia is available on the host |

Update `.github/workflows/ci.yml` step by step so that every job keeps
running real checks; never leave a job that passes because it runs nothing.

**M8b. Delete (M).** Remove the files listed in section 3, the PyO3 and CLI
code that only they needed, the schemas `run`, `result`, `result-v1`,
`ssa-run`, `ssa-result`, `tau-request`, `tau-result`, `archive`,
`archive-check`, `analysis-v1`, `state`, `check`, `fsp-result` (replaced by
the 1.0 files), and the old model reader except what the converter needs
(`contracts::{Model, Species, Reaction, Scalar, Rate}` move into
`convert.rs` as private types). Keep the original example models in
`models/legacy/` as converter test inputs; `models/examples/` holds the
converted `model/1.0` files. Gate: 4.3, I1–I3;
`grep -rn "xoshiro" crates` finds nothing; the line count of
`crates/*/src` is recorded in the log.

**M8c. Generated loops on the single crate (M).** The generator emits
`crn_craft::rules::` paths, and the Python loop driver
(`bindings/python/src/rules_loop_main.rs`, `python/crn_craft/_rules_loop.py`)
depends on `crn-craft` at the installed version. Regenerate the fixtures in
`tests/generated/` and state in the log that only paths changed
(`diff` of old and new fixtures after replacing the path prefix must be
empty). Move `crates/crn-craft-rules/tests/` into `crates/crn-craft/tests/`.
The facade crate keeps one smoke test. From here on I3 refers to the new
fixtures. Gate: 4.3, I1, I2, `check_compiled_loop.py --require-cargo`.

---

### M9. Downstream (M)

- `integrations/jobs`, `integrations/diffsol`,
  `validation/parameter-optimization/worker`: final cleanup after M5e/M5g;
  each has its own lockfile, so run `cargo build --locked` in each.
- `tutorials/optimization/`: `calibration.py`, `oscillators.py`,
  `sensitivity.py`, `rules.py`, `run.py`, `test_science.py` and the lesson
  texts C1–C4 on the unified API. Rerun each lesson's smoke preset and its
  science tests; the study presets only if the host allows.
- `examples/circuit`, `examples/tutorial/make_assets.py`.

Gate: each lesson's smoke preset and `test_science.py` pass;
`check_jobs.py` passes.

---

### M10. Documentation

**M10a. Tutorial (L).** Rewrite `docs/tutorial.md` on the unified API and
command line: first model with the builder, checking and diagnostics,
exact simulation, results as data frames, FSP, tau, budgets and
cancellation, reproducibility, analysis, import and export, choosing a
method. `validation/check_tutorial.py` runs every block; regenerate the
screenshots and figures with `examples/tutorial/make_assets.py`.

**M10b. Manual and references (M).** First freeze compatibility: copy
`docs/mtdna-user-manual.md` as it is on `main` at 0.6.0 to
`validation/compat/manual-0.6.0.md`, give `check_manual.py` a `--manual`
option, and add a CI step that runs every block of the frozen copy and
compares its output (I2b). The only allowed differences are the two blocks
whose expected text prints a generated-loop cache path or a timing, which
the checker already treats as variable. This step proves that a 0.6.0 user
script runs on 0.7.0 with the same results. Then
`docs/mtdna-user-manual.md` may
change: use the top-level names (`import crn_craft as crn`), say once that
`crn_craft.rules` is an alias, update printed schema strings and hashes by
running the blocks. `docs/rules-language.md` becomes the model reference
for `model/1.0` (keep the file name, add a redirect note); `contracts.md`,
`result-contract.md`, `archive-contract.md`, `diagnostics.md`,
`direct-ssa.md`, `alternative-ssa.md`, `reusable-parameters.md`,
`network-api.md`, `using-from-r.md`, `python/README.md` updated or merged.
`validation/check_manual.py` passes with the new expected outputs.

**M10c. Architecture and reports (M).** `architecture.md`, `ai-context.md`,
`AGENTS.md`, `technical-report.md`, `extending.md`,
`specialized-kernels.md`. Comparison documents: re-measure GillesPy2 (it is
pip-installable) and update `comparison-gillespy2.md`. For the Julia
comparisons, re-measure only if Julia is available and memory allows;
otherwise put a dated note at the top of each document saying the numbers
were measured with the 0.6.0 network engine. Never mix old and new numbers
in one table.

**M10d. README and changelog (M).** README: replace the "Two engines"
section by a short description of the single engine and API. `CHANGELOG.md`
`[Unreleased]`: Added (including quantile statistics and `Result.save`),
Changed, Removed, and a Migration section with
tables for (a) Python network API → unified API, (b) command line,
(c) contracts, (d) Rust paths (`crn_craft_rules::X` → `crn_craft::rules::X`,
still available through the facade), and the statement that results of the
0.6.0 network engine no longer replay byte-identically (D2) while rules
results do (D3).

Gate for M10: `check_tutorial.py`, `check_manual.py`, a link check over
`docs/` and the READMEs (every relative link resolves), 4.3.

---

### M11. Release 0.7.0 (S)

1. Merge `main` into `merge-0.7` one last time; open the PR `merge-0.7` →
   `main`; merge on green.
2. Release-prep PR: version 0.7.0 in the workspace, the binding,
   `pyproject.toml` and the five lockfiles (only the CRN-craft version
   lines change; `cargo metadata --locked` passes in each directory); crate
   README snippets `"0.7"`; CHANGELOG fold with today's date;
   `.github/workflows/publish-crates.yml` publishes `crn-craft`, then
   `crn-craft-rules` (facade, marked deprecated in its README and crate
   docs), then `crn-craft-cli`.
3. Tell the owner the commit to tag. The owner tags `v0.7.0`.
4. After the tag: both publish workflows green; PyPI and crates.io show
   0.7.0; a clean virtual environment passes `check_installed.py` normally
   and under `-O`; `compile_loop` builds against crates.io without
   `CRN_CRAFT_SOURCE`; `examples/r/sir.R` prints the documented output.

## 6. Risks and how the plan handles them

| Risk | Handling |
| --- | --- |
| A rules result changes bytes | I1 on every PR; new behavior only behind `model/1.0` and `run/1.0`. |
| The unified engine is slower on large networks | Baseline before any change (M0c), declared 1.25× gate (M4d), selectors ported before the old loops are deleted (M8b comes after M4). |
| Laws differ subtly after conversion | References stored while the old engine exists (M0c); the law gate runs in M2, M4a, M4b, M4c and after M8. |
| References are lost with the old engine | They are committed files under `validation/merge/reference/`. |
| Scope creep in validation migration | One PR per group (M8a); benchmarks of removed algorithms are deleted, not ported. |
| Julia comparisons cannot be re-measured | Dated note instead of mixed tables (M10c). |
| A 0.6.0 user script breaks or changes its results | `crn_craft.rules` keeps its names, defaults and schemas (M6a item 2, M4e default); I2 until M10, then the frozen manual I2b. |
| `main` needs a fix during the project | `main` stays releasable; fixes go there and are merged into `merge-0.7` (M0b). |

## 7. Deliverables checklist

- [ ] `validation/merge/` with baseline, references, RESULTS and scripts
- [ ] One crate with `rules/`, facade crate, unified CLI
- [ ] Contracts `model/1.0`, `run/1.0`, `result/1.0`, `fsp-result/1.0`, `analysis/1.0`, `archive/1.0` with schemas
- [ ] Converter and `models/legacy/`
- [ ] Limits, higher-order reactions, selectors
- [ ] Quantile statistics (`median`, percentiles) and `Result.save` for markdown and JSON (FR1)
- [ ] FSP, analysis, archives, ODE export, import/export on the flat view
- [ ] Unified Python API with `crn_craft.rules` alias
- [ ] Old engine deleted; CI jobs all still meaningful
- [ ] Tutorial, manual, references, README, CHANGELOG with migration
- [ ] Task logs in `docs/campaign/merge-0.7/`
- [ ] Release 0.7.0 verified from clean installs
