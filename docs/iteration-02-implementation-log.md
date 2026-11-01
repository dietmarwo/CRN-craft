# CRN-craft iteration 2 implementation log

**Date:** 1 November 2026

**Milestone:** M1a — compiled networks, propensity evaluation and analysis

**Version:** `0.1.0-dev.2`

**Repository:** https://github.com/dietmarwo/CRN-craft

**GitHub baseline:** `8134f510546c4b0711b1e5c21855f18d9b53e80a` (main, merged PR #1)

## Outcome

Implemented and tested the next planned iteration against the actual GitHub source. CRN-craft now compiles validated models into immutable indexed networks, evaluates all four supported rate laws, provides checked atomic reaction updates, and exposes `analyze MODEL [STATE]`. The analyzer includes bounded exact conservation laws and proven initial-class count bounds. Simulation remains unavailable; Direct SSA is M2.

The delivery is a Git mail patch ZIP with application instructions and checksums, plus this standalone Markdown log. It is based on the exact commit above, not the older first-iteration archive. No remote push, tag, package publication or deployment was performed.

## Baseline review

Cloned GitHub main, inspected its history, repository instructions, revised plan, source, CI and the supplied iteration-1 review. Main included the reviewed fixes: forced LF text checkout, serde_json float_roundtrip and its literal regression, corrected depth guidance, stronger population/FSP/output-budget cases, unsupported-contract cases and license/manifest cleanup. All 41 baseline Rust tests passed before changes. A later remote-ref check still reported the same main commit.

Kept these fixes intact. The original model/run/check/result 0.1 schema files remain byte-for-byte unchanged. No existing example model or run hash was edited.

The next milestone in the plan was M1a. Targets were written in `docs/decisions/iteration-02-validation.md` before implementation and numerical checks. Exact conservation was included to meet the known-enzyme-law acceptance gate; linkage classes and exporters remain M1b. No solver or integration infrastructure was added.

## Changes

| Area | Implemented behavior | Main files |
|---|---|---|
| Compilation | Validates typed models, resolves indices/parameters once, owns immutable state, sorts sparse terms by declaration order | `src/network.rs` |
| Stoichiometry | Sparse reactants/products and net changes equal products minus reactants; zero changes omitted | `src/network.rs` |
| Dependencies | Conservative reads include all availability gates and custom regulators/substrates; sorted unique invalidation lists; explicit 1,000,000-edge cap | `src/network.rs` |
| Rate evaluation | Mass action with binomial factors; Hill activation/repression; Michaelis–Menten; availability checked before arithmetic | `src/rates.rs` |
| Numerical handling | Scaled mass-action products, stable logarithmic Hill fractions, MM division before multiplication; non-finite and entirely underflowed positive rates fail explicitly | `src/rates.rs` |
| Count operations | Exact i64 state, complete named states, safety-limit checks and transactional reaction updates; no RNG or clock | `src/network.rs` |
| Conservation | Checked integer elimination and rational back substitution; exact basis verification, initial totals and bounds from nonnegative rows | `src/conservation.rs` |
| Analysis contract | Raw input hashes, initial/evaluation count distinction, structured failure and explicit unavailable-conservation status | `src/analysis.rs`, CLI |
| Documentation/CI | Network API, updated rate card, diagnostics, migration notes, current agent context, acceptance evidence, independent oracle and mutation checks in CI | `docs/`, `validation/`, workflow |

Rust source paths in the table are under `crates/crn-craft/` unless marked CLI. The workspace adds only pinned pure-Rust `libm = 0.2.15` to its external dependencies. Existing locked dependencies are unchanged. The dependency policy explicitly admits libm while retaining its closed allow-list and native-build checks.

### Semantics worth preserving

- A mass-action `2A` channel uses `c*x*(x-1)/2`; a custom law is already the complete propensity and receives no extra binomial factor.
- Reactants gate every law, including catalysts and constant custom propensities. Missing reactants return zero before possibly overflowing arithmetic.
- Dependencies describe propensity invalidation. A future scheduler must separately reschedule a fired channel.
- `fire_reaction` validates every proposed count before modifying any count. Failure cannot partially consume reactants.
- Conservation coefficients and totals are decimal strings to preserve i128 values in JSON. Status is `exact`, `size_limit` or `arithmetic_overflow`; unavailable status never implies absence of laws.
- Count bounds concern the model's initial conservation class. Arbitrary supplied evaluation states are allowed within declared safety limits and are not claimed reachable.
- `check` retains structural warning behavior. `analyze` suppresses a growth warning only when a law proves an initial-class bound.

## CLI and contract migration

```bash
cargo build --release --workspace --locked
./target/release/crn-craft describe
./target/release/crn-craft analyze models/examples/enzyme-mechanism.model.json
./target/release/crn-craft analyze models/examples/immigration-death.model.json models/examples/immigration-death.state.json
```

The named-state example sets A=5 and returns immigration/death propensities `[4.0, 1.0]`. The enzyme example exposes E+ES and S+ES+P with bounds on all four species. Complete outputs are retained in `validation/iteration-02/enzyme-analysis.json` and `named-state-analysis.json`.

New state and analysis contracts are `crn-craft.state/0.1` and `crn-craft.analysis/0.1`. Discovery advances to `crn-craft.describe/0.2` and advertises `analyze`, rate evaluation and analysis limits. Clients consuming discovery must recognize this version. Existing model/run/check/result consumers need no shape migration. Simulation methods remain marked contract-only.

Analysis exits 0 on successful structure/rate evaluation, including an explicitly unavailable conservation sub-result; it exits 1 on invalid input, I/O or numerical errors. No command emits a simulation result. Diagnostics CRN-E040–E044 and CRN-W010–W011 document new compile/state/numerical/update/conservation outcomes.

## Validation executed

| Gate | Result |
|---|---|
| Baseline regression before implementation | 41 Rust tests passed |
| Rust 1.85.0, 1.88.0 and stable 1.99.1 | 73 release tests passed on each; fmt and clippy with warnings denied passed |
| Python schema/CLI checks | 48 passed normally and under `-O`; reports byte-identical |
| Independent exact/Decimal oracle | 218 checks passed normally and under `-O`; reports byte-identical |
| Random-network coverage | 200 seeded networks; exact stoichiometry, reads, dependencies, rank, nullity, residuals and rational span |
| Numerical references | 1,012 evaluated rates; maximum relative error 2.3706671727190152e-14 against the predeclared 2e-12 tolerance plus eight minimum subnormals |
| Source mutation tests | 4/4 compiled mutants caught by real test failures |
| Dependency/formatting gates | Dependency checks normal and `-O` identical; five policy probes passed; Ruff format/lint passed |
| Fresh source patch checkout | Clean apply on baseline, autocrlf=true with all tracked files LF/empty, 73 Rust tests + 48 contract checks + 218 oracle checks passed |

Rust tests include exact homodimer/heterodimer results for counts 0–50, invalid states and indices, subnormal/overflow behavior, custom availability gates, atomic update failures and exact count updates above 2^53. Both sides of the new conservation species/reaction limits and dependency-edge limit are tested. Existing reviewed parser/hash/limit regressions remain.

The independent oracle uses Python integer combinations, 100-digit Decimal arithmetic and Fraction Gauss–Jordan elimination; it does not reuse Rust rate or elimination code. The seeded fixture stream SHA-256 is `15ca5e78084d38599d7e25fb72acad6564d18f4a001cbf95b80d3e9d6bf43f06`.

Four source mutations remove a binomial divisor, reverse a product sign, drop a custom-rate read or skip reactant gating. Each mutant compiled and failed named tests; compile/infrastructure failures are not counted as catches. The oracle also rejects three deliberately corrupted output structures.

Evidence is under `validation/iteration-02/`; M0 evidence remains historical and untouched. Local runs used locked Cargo dependencies and Linux x86_64, Python 3.12.14, Ruff 0.11.13 and jsonschema 4.26.0. Offline validation followed the initial fetch of pinned libm. CI was extended but this delivery does not claim a new remote GitHub Actions run.

## Development failures and fixes

1. The first test run found an incorrect manually entered steep-Hill reference literal. An independent 100-digit Decimal calculation established `0.44471580577189446`; the implementation returned `0.4447158057718944`. Corrected the expected literal, retained the failed run and reference calculation, and added the case to the independent oracle. No tolerance was relaxed.
2. Stable Rust 1.99.1 flagged an indexed elimination loop with needless_range_loop. Replaced it with disjoint mutable slices and zipped iterators. Retained the failed clippy log. All three toolchains and the independent oracle passed afterward.
3. Git patch whitespace checking identified surplus blank lines at the end of captured logs. Removed those trailing blank lines and Markdown hard-break whitespace; retained test/diagnostic content. Final patch whitespace checking passed.

No failed numerical fixtures were removed from the reported acceptance set.

## Limits and next iteration

This milestone does not include SSA/FSP/ODE execution, a random generator, scheduling, linkage classes, expression ASTs, exporters, bindings or ecosystem workers. There is no throughput, stochastic reproducibility or biological validation claim.

Exact conservation is capped at 64 species and 256 reactions and can stop on checked i128 overflow. The reported basis spans the rational conservation space; it is not a full integer-lattice or nonnegative-cone basis. Bounds from its nonnegative rows can be incomplete. Propensities use f64, so conversion of counts above 2^53 can round even though count storage and updates remain exact. Cross-platform bit identity and native Windows/macOS/ARM execution were not tested; LF checkout behavior was tested on Linux with core.autocrlf=true.

Next is the remaining M1b work: linkage classes and pinned reference exporters with explicit coefficient translation, plus expression nodes only for a concrete fixture. M2 then adds a separately tested Direct SSA baseline, live execution budgets, random streams and result/archive semantics.

## Applying the delivery

The ZIP contains `0001-m1a-network-analysis.patch`, this log, `APPLY.md`, checksums and final patch-verification evidence. Start from a clean branch at the exact baseline; `APPLY.md` gives the commands. The final delivery records the resulting commit/tree and validates the exact mail patch independently. Do not apply it blindly to an unrelated or dirty checkout.
