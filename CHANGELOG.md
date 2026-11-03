# Changelog

All notable changes are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) with breaking contract changes
called out. Every contract is versioned on its own (`crn-craft.<kind>/<n>`).

## [Unreleased]

## [0.2.0] - 2026-11-03

The first release: the crates `crn-craft`, `crn-craft-rules` and
`crn-craft-cli` on crates.io and the Python package `crn-craft` on PyPI, with
abi3 wheels for CPython 3.10+ on Linux x86_64, Windows x64 and macOS
x86_64/arm64 and a source distribution. The Python bindings, the diffsol
integration and the job worker stay source-only workspaces. The development
iterations that make up this release follow.

### Highlights

- Exact SSA: legacy checked Direct and dependency-cached Direct archives, and
  `direct`, `cached_direct`, `direct_cr`, `sorting_direct` and `rssa_cr` with
  inverse or Ziggurat waits (`crn-craft.ssa-run/0.3`), where `fast`, emulating
  JumpProcesses.jl, is the default and `fast: false` the checked policy, plus
  Gibson–Bruck next reaction (`crn-craft.ssa-run/0.4`, checked policy).
- Sparse finite state projection, explicitly approximate Poisson leaping,
  network analysis with conservation laws and cone enumeration.
- Experimental rules engine for dynamic, attributed species (native CLI and
  Python `crn_craft.rules`), with exact-law gates for mtDNA and aggregation
  models: parallel ensembles, event trackers and callbacks, Static Direct,
  next reaction, a specialized fixed-population mtDNA sampler, approximate tau
  and binomial steps, ATP-coupled thinning and generated Rust kernels.
- Python bindings with prepared simulations, reusable parameter templates,
  trajectory diagnostics and strict SBML/BioNetGen import subsets; optional
  diffsol BDF with forward sensitivities; fcmaes job applications.
- Optimization tutorials C1–C4 and measured comparisons with GillesPy2 and
  Catalyst/JumpProcesses.
- Packaging: `python-release.yml` builds and tests every wheel and the source
  distribution with the installed-package acceptance and publishes a `v*` tag
  to PyPI through Trusted Publishing; `publish-crates.yml` verifies the three
  crates from their own archives and publishes them in dependency order.
  `scripts/package_version.py` and `scripts/check_release.py` check versions,
  crate metadata and lockfiles. The manual `wheels.yml` workflow is replaced.
- Archives record `cargo_lock_sha256` from a build-script hash of the lockfile
  (the packaged crate's own, or the checkout's root one) instead of
  `include_bytes!` of a path outside the crate, which no packaged crate could
  build. Checkout builds record the same hash as before; see
  `docs/archive-contract.md`.

### mtDNA ensembles and faster generic rules

#### Changed: faster generic rules simulation (S1–S6)

- Static Direct, next reaction, ATP thinning and the leap methods cache rule rates and shared vacancy bases, recompute only the channels and observables an event affects, defer report-only observables, keep the live-instance count incrementally and evaluate the stop and predicates from compiled bytecode. Existing requests return byte-identical results.
- `rules-run/0.4` adds two opt-ins on `static_direct`: `waiting_time: "ziggurat"` and `sparse_trackers: true`. Each has its own sampler identity; the result stays `rules-result/0.2`.
- `generate_rules_kernel` emits a Rust kernel for one closed model with constant or vacancy rates. `execute_with_rules_kernel` runs it under the Static Direct scheduler with byte-identical results. It refuses kernels generated for other model bytes or by another generator version, because the fingerprint covers the emitted source. Rust only; it is not available from Python.
- [docs/specialized-kernels.md](docs/specialized-kernels.md): when a model-specific kernel is worth it, and the review checklist. Review log: [docs/campaign/fast-generic-rules-review-log-2026-11-03.md](docs/campaign/fast-generic-rules-review-log-2026-11-03.md).

#### Added: mtDNA ensembles, event observation and alternative rules methods (A0–A6)

- Three-class mtDNA fixtures (`models/rules/mtdna-three-class*.rules.json`) with independent CME gates, and a usage guide in [docs/mtdna-rules.md](docs/mtdna-rules.md).
- `crn_craft.rules`: prepared parallel rules ensembles from Python with read-only NumPy record columns, cooperative cancellation, budgets, `bind_parameters` candidates, live and batch event callbacks, trackers (crossings, episodes, extrema, integrals, first passage) and a bounded event log. The Rust `EventObserver` interface exposes the same committed events.
- Rules methods: `static_direct`, `mtdna_direct` (exact grouped Direct for the fixed three-class model, with an explained Static Direct fallback, `CRN-R051`), `next_reaction` (Gibson–Bruck), approximate `tau` and `binomial`, and `thinning` for the ATP-coupled model (`rules-model/0.2` `atp` block).
- Contracts: `rules-run/0.2` (method, trackers, event log) and `0.3` (adds `mtdna_direct`) answer with `rules-result/0.2`; `rules-model/0.2` adds `atp`. `rules-run/0.1`, `rules-model/0.1` and the tree method keep their meaning and result bytes, and reject the new fields even when empty or null. Network `ssa-run/0.4` is 0.3 plus the `next_reaction` aggregator (checked policy and inverse waits only) and answers with `ssa-result/0.4`.
- Discovery: core and Python `describe/0.10`, native CLI `describe/0.11`, with the 0.4 request default and the rules methods.
- Migration: stored requests of earlier versions are unchanged. `crn.ssa_spec(aggregator="next_reaction")` writes 0.4 with `fast=False`; other aggregators still get 0.3. A present `method`, `event_log`, `trackers` or `atp` must hold a value; `null` is rejected.
- Review log: [docs/campaign/mtdna-a0-a6-review-log-2026-11-03.md](docs/campaign/mtdna-a0-a6-review-log-2026-11-03.md).

### Alternative exact SSA methods and the Julia-gap campaign (`0.2.0-dev.1`)

#### Changed: `fast` is the default (`crn-craft.ssa-run/0.3`)

- New request version `crn-craft.ssa-run/0.3` with the same fields as 0.2; `method.fast` defaults to `true`. It answers with `crn-craft.ssa-result/0.3`, which always states `method.fast`. Set `"fast": false` for the checked exact policy (ChaCha8 words, lost-increment and 53-bit resolution checks, stop polls before every random word). `fast` must be a boolean; `null` is rejected.
- `crn.ssa_spec` writes 0.3 requests and now defaults to `fast=True`; pass `fast=False` for the checked policy. Non-boolean values raise `TypeError`.
- Migration: stored `ssa-run/0.1` and `0.2` requests keep their meaning and return the same bytes as before (0.1 rejects `fast`, 0.2 keeps it off unless set). Python callers that relied on the old default get different (equally exact) paths unless they pass `fast=False`; the checked replay evidence scripts now pass it explicitly. With the same seed, a 0.3 request with `fast` on reproduces a 0.2 request with `fast: true`, and a 0.3 request with `fast: false` reproduces a 0.1 request; only the schema strings, `method.fast` and `run_sha256` differ.
- Discovery: core and Python `describe/0.8`, native CLI `describe/0.9`; `alternative_exact_ssa.fast` lists the default per request version.

#### Added: simulator optimization tutorials (contributed by Astra)

- Four lessons in [tutorials/optimization](tutorials/optimization/README.md) with native bindings, independent moment, BDF, generator and finite-state checks, disjoint validation and bounded smoke and study presets. SSA runs use the default `fast` policy; C1's `--include-checked` adds a checked arm and every row records policy, sampler and request schema.
- Review fixes: C3's cross-backend control uses its own seeds instead of the training seed that selected the design, a horizon whose spectral grid contains the 20-unit target period, and per-path oscillatory fractions; C1 records the training loss of the true rates and whether a fit stayed at its truth-derived start; C2 reports local condition numbers and a terminal-product twin instead of a dimension-bound rank, and treats failed native solves as undefined points; C4 reports each kernel's exact-law optimum. Options a chapter does not use are rejected, summaries carry a run state, and `results/verified.json` is written by `record.py` from complete runs. The C3 reference (fcmaes-rust `crn_reference`) takes the initial state from the request and rejects unknown fields.

#### Added: dynamic, attributed species (experimental rules engine, K0/K1)

- New crate `crn-craft-rules` and CLI commands `rules-check` and `rules`: species types with numeric attributes, instances created and removed on the fly, rule schemata over zero to two reactant patterns with product-form propensities, observables (`sum`, `distinct`, `max`, `min`), `new_id()` clones, random product attributes and stop conditions. Contracts `crn-craft.rules-model/0.1`, `rules-check/0.1`, `rules-run/0.1`, `rules-result/0.1`; the native CLI's discovery moves to `describe/0.8` (core and Python stay at 0.7). Reference: `docs/rules-language.md`.
- K1 gates (`validation/check_rules.py`, CI job `rules`): 30 exact-law comparisons over 100,000 replicates per fixture pass (largest distance 0.0034 against a bound of 0.0074), including Kowald's fixed-population mtDNA model and constant and multiplicative aggregation; all six wrong-law power checks reject. After an independent review, identical-pattern pairs use exact pair trees (sampler `chacha8-u53/rules-direct/2`), statistics stream in replicate order with admission limits, and overflows and undefined conditions fail instead of passing silently; the gates were rerun unchanged.

#### Added: reusable parameter bindings and generated kernels (imported, by Astra)

- `parameters::ParameterTemplate` (Rust) and `crn.prepare_parameters` (Python) compile a model's topology once and bind parameter vectors or exact name maps per candidate; each binding keeps byte snapshots for replay. Optional `codegen::mass_action_source` / `crn.mass_action_kernel_source` emit a Rust propensity kernel for `ssa::execute_with_kernel`, guarded by `CompiledNetwork::topology_sha256`. Results are byte-identical to fresh preparation. See `docs/reusable-parameters.md` and the optimization comparison with rebop.
- Review fixes: generated kernels fall back to the new `CompiledNetwork::propensity_in_scheduler_state` instead of the state-scanning `propensity` (a 1,500-reaction review model went from about 40 s to 4.2 s with `direct`); ordinary preparation again checks admission before compiling the dependency graph (`CRN-E060` before `CRN-E040`, as before the import); binding mismatches are `CRN-E014` with missing and unknown names, invalid Python seeds `CRN-E020`; Python `bind` rejects sets and turns float overflow into `CrnError`.

#### Measured: fast mode versus Julia at one and eight threads

- `docs/comparison-julia-8-threads.md`: same-host measurement of the ten `fast` configurations against five JumpProcesses 10.0.1 methods on W1–W7, original and 256-path batches, at one and eight workers (2,240 timed observations; 6,514 terminal-marginal checks passed). Geometric mean of best fast / best Julia: 0.69x (original) and 1.04x (256 paths) at eight workers; 1.15x and 1.12x at one worker (derived on import from the same observations), replacing the earlier serial projection of 1.02x. The harness and raw evidence were not part of the import; all summary statistics were recomputed from the reported observations.

#### Added: the `fast` switch (`crn-craft.ssa-run/0.2`)

- `method.fast` (default `false`) emulates JumpProcesses.jl: Xoshiro256++ words per replicate, running totals without lost-increment or 53-bit resolution checks, and stop polls every 256 event-loop iterations. New methods run 1.09–4.3x faster with it; it passes its own 210-comparison exact-law gate, now also in CI.
- `ssa-run/0.2` and `ssa-result/0.2`; 0.1 requests are still accepted, reject `fast`, and return byte-identical results. Discovery moves to `describe/0.7`. Python: `crn.ssa_spec(..., fast=True)`.
- Benchmarks time CRN-craft with `fast` and report the checked default alongside: the GillesPy2 comparison was rerun (fast cached Direct 1.5–8.2x and the best fast method 2.9–9.2x faster than GillesPy2's C++ solver on W1–W4); the Julia comparison gained a projected `fast` table (geometric mean 1.02x).

#### Imported: ten further Julia-gap iterations (C01–C10, by Astra)

- Kept C05 (lower-bound reuse), C06 (shared exact custom-law ranges) and C09 (compiled elementary mass-action kernels); same-host geometric-mean gap to Julia 1.54x to 1.44x. Evidence in `validation/julia-gap10/`.

#### Added (imported acceleration patch)

- `crn-craft.ssa-run/0.1` and `crn-craft.ssa-result/0.1` with `direct`, `cached_direct`, `direct_cr`, `sorting_direct` and `rssa_cr` aggregators, inverse or Ziggurat (Julia 1.13.3 tables, MIT) waits, and the `ssa` command in both CLIs; discovery `describe/0.6` lists them under `alternative_exact_ssa`.
- Reusable prepared handles in Python: `crn.prepare`, `crn.prepare_ssa`, `crn.ssa_spec`, `crn.simulate_ssa`.
- `docs/comparison-catalyst-jumpprocesses.md` on seven workloads, with the evidence in `validation/ssa-acceleration/`.

#### Performance (declared in docs/decisions/julia-gap-2026-11-02.md)

- Kept: ChaCha8 AVX2 dispatch (A02), exact Hill-rate memo by regulator count (A03, W3 about 2x on every method, legacy included), cached composition envelope (A05), sums-only lazily checked rate tree (A06), and one word for composition member and acceptance (A08; DirectCR and RSSACR sampler version 2, 210/210 exact-law checks).
- Rejected and reverted with evidence in `validation/julia-gap/`: A01 (poll per loop), A04 (redundant RSSACR evaluation), A07 (additions-only ordered sums), A09 (batched tree rebuild), A10 (bucket buffer policy).
- DirectCR, SortingDirect and RSSACR 1.07–2.4x faster than the import, cached_direct 1.01–2.05x; projected gap to JumpProcesses 1.11–1.83x (geometric mean 1.46x, from 2.02x measured). Legacy Direct and cached results are byte-identical.

#### Fixed

- SortingDirect no longer fails a replicate when its ordered scan sum rounds one ulp below the tree total; the threshold is redrawn.
- The alternative-SSA output budget is applied once, to the `ssa-result` JSON, instead of after a legacy archive allowance.
- `python -m crn_craft ssa` exits 1 for rejected or incomplete runs and reports diagnostics like the native CLI; `prepare`/`prepare_ssa` convert non-JSON native errors to `CrnError`.
- The native CLI `--version` milestone matches `describe`; the notices no longer point to a file that does not exist here.

### Campaign review and performance iterations P11–P20 (`0.2.0-dev.1`)

#### Performance (declared in docs/campaign/performance-p11-p20.md)

- P11: Direct selects the channel by binary search over the prefix sums recorded during the strict sum.
- P12: `simulate-cached` resumes the strict sums at the lowest changed channel (1.59x on the sparse 66-reaction network).
- P14: one scheduler poll per event; P16: per-worker buffers reused across replicates.
- P17, P18: FSP skips exact-zero rows, then stores the generator by column and fuses the series accumulation (W8b 2.34x against P10, bitwise identical probabilities).
- P19: tau leaping draws leap means above 16 with Hörmann's PTRS; the algorithm identifier becomes `/0.3` (W9 312x).
- Rejected and reverted with evidence: P13 (skip bitwise-unchanged rates) and P20 (single archive serialization). P15 (flat save buffers) gave no speedup and is kept as the memory-accounting fix.
- Final source: `simulate` 1.4–2.0x faster than 5a07334 on the frozen workloads, `simulate-cached` 2.1–3.3x faster than `simulate` on the 12- and sparse 66-reaction networks; all SSA digests unchanged.

#### Fixed

- A negative count from firing is reported as `CRN-E043`, not as a population-limit stop.
- Saved states use one flat buffer per path, so the published memory reservation is an upper bound (it was exceeded for many save points).
- FSP enforces `max_output_bytes` on the exact serialized size (`CRN-E070`) and sums the omitted series weight explicitly.
- Tau leaping caps leap means at 10^7, uses its own random-stream domain, accepts `--cancel-file`, and its CGP factors are tested.
- `sampled_first_passage` reports only grid-supported bounds; `oscillation` needs 16 samples after burn-in and handles constant series.
- Importers reject BioNetGen fixed species, missing SBML stoichiometry or units, `fast` reactions, identifier collisions and pathological literals within one second (`CRN-E090`).

#### Added

- `docs/comparison-gillespy2.md` and `validation/compare_gillespy2.py`: CRN-craft against GillesPy2 1.8.3 on W1–W4 (`simulate-cached` 1.1–3.3x faster than its C++ solver, 43–126x faster than its NumPy solver), with equivalence checks against each other and exact laws.
- README banner and an index linking every Markdown document.
- Python: Ctrl-C cancels native calls; `cancel_file=` for `simulate`, `fsp` and `tau`; `CrnError` with parsed diagnostics; `__version__`; Python CLI `--version`, `--cancel-file` for `fsp`/`tau` and native-CLI output formatting.
- `docs/tutorial.md` with screenshots and figures generated by `examples/tutorial/make_assets.py` and checked by `validation/check_tutorial.py`; `docs/architecture.md`; `docs/technical-report.md`.
- CI job building the wheel and running the Python, importer, tutorial and installed-package checkers normally and under `-O`.
- Two mutation-gate mutants for the lost-increment checks of the prefix-sum and resumed-sum passes (45 in total).
- `validation/benchmark_iterations.py` (paired iteration timings) and `validation/compare_tau_samplers.py`.

### Iteration 4 / M2 (`0.1.0-dev.4`)

#### Added

- Implement Direct SSA for existing count propensities, exact save-time boundaries and absorbing states.
- Pin ChaCha8 streams and owned uniform/wait/channel mappings; reject unresolved dynamic ranges and clock stagnation.
- Add bounded ordered parallel ensembles, reproducible mean/sample variance and mean SE.
- Enforce event/population/deadline/output/workspace limits and cooperative API/file cancellation; retain every replicate's accounting and withhold partial-ensemble statistics.
- Emit result 0.2 and checksummed archive 0.1; add `simulate`, `verify`, generated schemas and migration notes. Preserve result-v1 and existing inputs/analysis/export.
- Add the predeclared 300,000-path analytic campaign, independent scalar RNG and trajectory oracles, archive-corruption tests and nine new source mutants.
- Preserve all three iteration-3 review regressions and the 17-mutant baseline. No FSP or comparative performance claim.
- Review tests for lost positive increments and underflowed waits at the sampler boundary, retained-trajectory IDs, existing and live cancel files, tampered partial-evidence archives and the archive byte budget, and CLI `simulate`/`verify` exit codes. Eight probe mutants had passed every delivered test.
- `validation/check_cme.py`: four 100,000-path fixtures with second-order, Hill and Michaelis–Menten channels against exact master-equation laws (12 DKW tests, family alpha 0.001, declared in `docs/decisions/iteration-04-review-cme.md`). It runs in CI normally and under `-O`.
- Nine review mutants in the mutation gate, now 35.

#### Fixed

- `--cancel-file` no longer slows simulation about tenfold, and a wall-time limit no longer doubles the cost per event. Every poll used to stat the file, which costs about a microsecond against roughly 100 ns per event for small models, and read the clock. A per-worker poller now reads the cancellation flag at every poll and the wall clock and cancel file at every 64th poll, the file at most once per millisecond. A found file cancels all workers, and a worker keeps an observed stop. Numerical results are unchanged.

### Iteration 3 / M1b (`0.1.0-dev.3`)

#### Added

- Exact nonnegative conservation-cone rays with deterministic work/ray caps, checked i128 arithmetic and explicit unavailable statuses.
- Effective initial-class bounds that resolve the reviewed mixed-sign-basis gap when cone enumeration completes.
- Reaction complexes, linkage/strong/terminal classes, weak reversibility and structural deficiency, with iterative graph walks.
- GillesPy2 1.8.3 elementary mass-action Python export, explicit homodimer conversion, collision-free aliases and unsupported-feature diagnostics.
- Independent support-enumeration/reachability oracle for 300 seeded networks, real target propensity/SSA-smoke checks and seven additional source mutants.
- Tests and mutants for gcd reduction of combined cone rays, a caller ray cap below the number of initial orthant rays, and linkage classes joined only through converging edges, which three probe mutants had passed; the structure oracle gains the `primitive_combination` model.

#### Changed

- Analysis advances to 0.2; `analyze-v1` preserves 0.1. Discovery advances to 0.3; export starts at 0.1. Model/run/check/result/state schemas remain unchanged.
- W001 suppression uses exact cone bounds when available and retains sound basis bounds on unavailable outcomes.
- CI covers structure/oracle checks under normal Python and `-O`, 17 mutants after a passing baseline, and optional pinned reference validation.

#### Limits

- Cone enumeration can grow exponentially; size alone does not ensure completion. No partial cone is labelled complete.
- Export supports mass-action reactant order at most two, no population safety limits and exactly representable initial counts up to 2^53. It is for reference SSA, not target ODE semantics or CRN-craft budget/range equivalence.
- Direct SSA and all CRN-craft simulation remain M2 or later.

### Iteration 2 / M1a (`0.1.0-dev.2`)

#### Added

- Immutable compiled networks, sparse stoichiometry, reactant gates and propensity dependency graphs.
- Checked mass-action, Hill activation/repression and Michaelis–Menten evaluation using pinned pure-Rust libm 0.2.15; atomic i64 reaction updates.
- `analyze MODEL [STATE]`, new state/analysis 0.1 schemas, exact checked small-network conservation with explicit cap/overflow status and initial-class bounds.
- Independent 200-network Fraction/Decimal oracle, source mutation tests, numerical edge/atomicity regressions, documentation and retained evidence.
- Tests for the `log1p` Hill path near x = k, for weight division in conservation bounds and for binomial symmetry at high reactant orders, which three probe mutants had passed.

#### Changed

- Discovery is `crn-craft.describe/0.2` and advertises implemented rate evaluation; model/run/check/result 0.1 contracts are unchanged.
- Analyzer growth warnings are suppressed only by proven initial-class bounds. Structural `check` semantics are unchanged.
- CI includes the independent network oracle normally and under `-O`, source mutation checks and network tests in the LF checkout job.
- The mutation check runs the unmutated suite first and requires it to pass, then runs seven mutants.

#### Limitations

- No scheduler, RNG, SSA/FSP/ODE, linkage classes, expression AST, exporters or ecosystem integrations yet.
- Conservation uses checked i128 within 64 species / 256 reactions; its nonnegative basis heuristic is not a complete bound solver.
- f64 propensity evaluation is tested against declared tolerances, not certified correctly rounded or cross-platform identical.

### Iteration 1 / M0 (`0.1.0-dev.1`)

#### Added

- M0 Rust workspace with strict model/run contracts and reserved result envelope.
- JSON duplicate-key, unknown-field, domain, identity, budget and truncation checks.
- Structured `describe` and `check` CLI commands with explicit unavailable simulation status.
- Generated schemas, complete examples, negative fixtures and regression/acceptance suites.
- Local validation evidence, CI configuration, revised plan and agent/release documentation.
- Tests for unsupported run versions, methods and aggregators, and for both sides of the JSON nesting limit.
- Tests for both sides of the population-limit, FSP-bound and minimum-output-budget rules, which three probe mutants had passed.
- A CI job that checks out the repository with `core.autocrlf=true` and runs the contract suite on it.

#### Changed

- CI cancels superseded runs and runs the dependency gate once, normally and under Python `-O`, comparing the two reports.
- `SOURCE_MANIFEST.json` from the delivery archive is not tracked; git history identifies the files.

#### Fixed

- Decimal literals parse to the correctly rounded double (serde_json `float_roundtrip`). The default parser read about a quarter of random literals one ulp away, including shortest round-trip forms such as `3.983236103458967e-14`.
- Text files keep LF line endings in every checkout (`.gitattributes`), so exact-byte model hashes hold under Git for Windows' `core.autocrlf=true`.
- The nesting-limit hint states the rule the parser applies: depth at most 48, with the root at depth 0.

#### Limitations

- No propensity evaluator, network analyzer, SSA, FSP, ODE, bindings, importers or ecosystem worker yet.
- No simulator correctness, throughput or biological validation is claimed.

### 0.2.0-dev.1 campaign

- Sparse reachable FSP with escaping sink, stable split uniformization and separate error accounting.
- Exact dependency-cached Direct with original ordered streams/results; external cancellation checkpoint.
- Isolated Python ABI3 package, trajectory diagnostics and strict SBML/BioNetGen import subsets.
- Native fcmaes job applications, fixed/disjoint seed accounting and actual shared SDK/runner acceptance.
- Optional diffsol BDF, analytic Jacobians/forward sensitivities and global-to-SQP fitting.
- Experimental non-negative Poisson leaping with explicit approximate request/result schemas.
- Ten performance experiments and same-host reference/quality evidence; no package publishing.

[Unreleased]: https://github.com/dietmarwo/CRN-craft/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/dietmarwo/CRN-craft/releases/tag/v0.2.0
