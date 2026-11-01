# Iteration 3 / M1b validation

Date: 1 November 2026. Base: `c167add312a271671bdc31c95162b5b5c1227ee1` (merged PR #2). The untouched baseline passed 76 Rust tests, 48 contract checks and 218 network-oracle checks before implementation.

## Final local gates

| Gate | Result | Evidence |
|---|---|---|
| Rust 1.85.0 | fmt, clippy with warnings denied, 97 release tests | `rust-1.85-*.log` |
| Rust 1.88.0 | fmt, clippy with warnings denied, 97 release tests | `rust-1.88-*.log` |
| Stable 1.99.1 | fmt, clippy with warnings denied, 97 release tests | `rust-stable-*.log` |
| Schema/CLI acceptance | 58 passed normally and under `-O`; byte-identical reports | `contracts*.json` |
| Preserved M1a oracle | 218 passed normally and under `-O`; unchanged fixture stream and numerical results | `network*.json` |
| New cone/graph oracle | 308 passed normally and under `-O`; byte-identical reports | `structure*.json` |
| Actual GillesPy2 1.8.3 | 6 grouped checks, 15,717 propensity comparisons and two-trajectory smoke test; reports identical under `-O` | `export*.json` |
| Source mutation gate | Passing unmutated network/structure suites required; 14/14 mutants caught | `mutations.json` |
| Dependency gate | Passed normally and under `-O`, identical reports; no new external Rust dependencies | `dependencies*.json` |
| Python quality | Ruff format/lint passed | `python-format.log`, `python-lint.log` |

Rust test breakdown: 35 contracts, 33 reviewed network regressions, 19 new structure/export tests and 10 CLI tests. The final mail patch is rebuilt and tested separately on a fresh baseline checkout; its post-commit verification lives in the downloadable ZIP to avoid a self-referential patch artifact.

## Independent exact results

The structural oracle uses support enumeration and Fraction Gauss–Jordan nullspaces, independently of Rust's incremental double description. It checked 300 seeded networks and four targeted models, then rejected four corrupted output probes. Complete primitive ray sets, totals, residuals, bounds, reaction complexes, linkage/SCC/terminal classes and deficiency matched exactly. The cone returned 372 checked extreme rays. Every tested small cone completed exactly within the fixed limits.

The prior basis heuristic bounded 507 species across the model corpus; complete cone analysis bounded 521, a gain of 14. The review counterexample gained all four missing bounds. The corpus differs from the review's 300-network sample, so its counts are reported separately.

Seed: 20261103. Fixture stream SHA-256: `704c72253352754584904bb700702326a992bc6f78fbbe4ff3b803da0dec94ee`. Targets and tolerances were declared in `docs/decisions/iteration-03-validation.md` before the implementation. No gate was loosened after observing results.

The old M1a independent corpus remains seed 20261101 and SHA-256 `15ca5e78084d38599d7e25fb72acad6564d18f4a001cbf95b80d3e9d6bf43f06`: 1,012 rate evaluations, maximum relative error 2.3706671727190152e-14, against 2e-12 plus eight minimum subnormals.

## Export verification

Installed GillesPy2 1.8.3 actually constructs the generated models. Its own generated propensity expressions are compared with integer-combination references and the Rust CLI over all 51×51 count pairs for six elementary channels, plus the shipped immigration–death and enzyme models. The count grid is exact; maximum target relative error on the remaining examples is 1.1328806373726087e-16. Non-unit source Omega, homodimer factorials, aliases, catalysts, zeros and parameter substitution are covered.

The reference SSA smoke uses NumPySSASolver, seed 20261103, two trajectories and save times 0, 1, 2. It checks finite nonnegative integer populations and completion; it is not a CRN-craft solver accuracy or throughput experiment. Runtime package versions and the installed target reaction-source SHA-256 are in `reference-environment.json`. Optional Python reference dependencies are separate from the pure-Rust runtime.

## Failures retained and resolved

- Rust 1.85 clippy flagged a comparison chain; replaced it with `Ordering` matching. See `development-1.85-clippy.log`.
- Stable clippy flagged a cloned one-element test slice; replaced it with `slice::from_ref`. See `development-stable-clippy.log`.
- One later local oracle invocation encountered a build-output executable with mode 0644. Restored its execute permission and reran the checker; no source or tolerance change. See `development-executable-permission.log`. The independent fresh build verifies the source delivery separately.
- The fourteenth mutant initially had an ambiguous text anchor matching both cone and fallback bounds. The harness rejected it; the anchor was made unique and the entire gate rerun successfully. See `development-mutation-anchor.log`.
- All numerical/structural acceptance runs completed without failed fixture rows. Intentional source mutation failures are recorded by named test in `mutations.json`; a build error or pre-existing failing baseline is never scored as a catch.

An extra regression uses a high-weight chain followed by a source equation: cone arithmetic overflows before the last equation, while exact linear elimination still succeeds and proves an isolated species bound. The analysis preserves that bound and labels the cone unavailable. Its dedicated source mutant is the fourteenth mutation.

## Scope and limitations

Linux x86_64 local execution only. GitHub workflow additions are configured, not claimed as a new hosted CI result. The original model/run/check/result/state schemas remain byte-identical, and generated `analysis-v1` is byte-identical to the prior analysis schema. Current analysis is 0.2, export 0.1 and discovery 0.3.

Cone completeness applies only to `exact` outcomes. Size, intermediate-ray, work and checked-i128 limits can prevent completion; no partial cone is emitted. Bounds concern nonnegative linear conservation laws in the initial class, not arbitrary evaluation states or complete dynamical reachability. Graph deficiency is structural. Reference export is restricted to elementary mass action and does not preserve target ODE semantics, live safety stops or all floating-point behavior. No CRN-craft stochastic solver or performance comparison is implemented.

## Repository review, 1 November 2026

The delivered mail patch `0001-m1b-cones-linkage-export.patch` (SHA-256 `936e3d7c…`) applied cleanly to `c167add` and produced tree `2f56d4f`. The delivery's gates reproduced on Linux x86_64 with Rust 1.97.0, Python 3.13.18 and the pinned GillesPy2 1.8.3 stack:

- 97 release tests;
- 58 contract, 218 network, 308 structure and 6 reference-export checks, with every report byte-identical to the delivered one and under `-O`;
- the dependency gate;
- 14 of 14 mutants caught.

Independent checks run in review:

| Check | Result |
|---|---|
| `analyze-v1` against the previous release | Byte-identical output to the `c167add` binary's `analyze` on all 1514 models: every shipped model and 1500 random networks with mass-action, Hill and Michaelis–Menten rates |
| Effective bounds against a linear program | For each species of every random network, `initial_class_upper_bounds` equals ⌊min T·l subject to Nᵀl = 0, l ≥ 0, lᵢ = 1⌋ (HiGHS via SciPy), or null when infeasible: 6783 species, 0 mismatches. The cone bounds are therefore the best bounds that any non-negative linear conservation law gives |
| GillesPy2 homodimer convention | GillesPy2 1.8.3 `Reaction._create_mass_action` builds `k * X * (X - 1) / vol` for 2X, confirming the exporter's halving |
| Cone runtime | 40 random networks with 16–64 species and 8–256 reactions: all exact, at most 0.021 s per `analyze` |
| Mutation probe of 20 cone, linkage, export and analysis rules | 17 caught by the delivered Rust tests. Survivors: combined rays not reduced to primitive form, the initial ray-cap check removed, linkage classes from outgoing edges only. The structure oracle caught only the third. Each now has a Rust test, and the first also an oracle model (`primitive_combination`) |

After the review:

- 100 release tests: 35 contracts, 33 network, 22 structure, 10 CLI;
- the structure oracle passes 309 checks; the fixture stream SHA-256 is now `1373b12d6c9601fe27e3d5424ff0fce260c9708af5face695308366df9cefa98`, and 525 species are bounded by the cone against 511 by the basis;
- the mutation gate catches 17 of 17 mutants after a passing baseline.
