# Iteration 2 / M1a validation results

Date: 1 November 2026. Baseline: GitHub `dietmarwo/CRN-craft` main at `8134f510546c4b0711b1e5c21855f18d9b53e80a` (merged PR #1). The baseline's 41 Rust tests passed before implementation. This evidence supplements the historical M0 report in `validation/RESULTS.md`; it does not replace it.

## Executed local checks

| Check | Result | Evidence |
|---|---|---|
| Rust 1.85.0 | fmt, clippy with warnings denied, 73 release tests passed | `rust-1.85-*.log` |
| Rust 1.88.0 | fmt, clippy with warnings denied, 73 release tests passed | `rust-1.88-*.log` |
| Installed stable 1.99.1 | fmt, clippy with warnings denied, 73 release tests passed | `rust-stable-*.log` |
| Schema/CLI acceptance | 48 passed, 0 failed; normal and `-O` reports byte-identical | `contracts.json`, `contracts-optimized.json` |
| Independent network oracle | 218 checks passed, 0 failed; normal and `-O` reports byte-identical | `network.json`, `network-optimized.json` |
| Source mutation tests | 4/4 compiled mutants rejected by actual failing tests | `mutations.json` |
| Dependency policy | Passed normally and under `-O`, byte-identical; five policy mutation probes retained | `dependencies.json`, `dependencies-optimized.json` |
| Python formatting/lint | Ruff format and lint passed | `python-format.log`, `python-lint.log` |
| Fresh baseline checkout | autocrlf=true; LF preserved; 73 Rust + 48 contract + 218 network checks passed | `fresh-checkout.json`, `fresh-checkout-*.log/json` |
| Original schemas | Model/run/check/result 0.1 files unchanged; all six snapshots match CLI-generated shapes | Contract checks and git diff |

All Cargo validation used `--locked`; local builds also used `--offline` after the pinned libm dependency was fetched. Tests run on Linux x86_64. Environment versions are recorded in `environment.json`. GitHub Actions was extended with these checks, but the new workflow has not been run remotely by this delivery.

## Independent numerical and structural evidence

The oracle generated 200 deterministic small networks with seed 20261101. It also checked four shipped models and numerical boundary cases, then rejected three deliberately corrupted outputs. Across successful numeric fixtures it compared 1,012 rates with Python integer `math.comb` and 100-digit Decimal formulas using the parsed f64 parameter values. Maximum relative error was **2.3706671727190152e-14**, below the predeclared **2e-12** relative tolerance plus eight minimum subnormals absolute tolerance. Exact zeros and unrepresentable values use explicit separate checks.

Fixture stream SHA-256: `15ca5e78084d38599d7e25fb72acad6564d18f4a001cbf95b80d3e9d6bf43f06`.

Sparse stoichiometry, reads and dependencies matched direct Python calculations. Fraction Gauss–Jordan elimination checked exact rank, nullity, residuals, basis independence and rational span against the Rust integer elimination. Enzyme laws E+ES and S+ES+P were found. Bound totals and the absence of bounds from mixed-sign rows were checked independently. This covers a deterministic small-network sample, not exhaustive verification.

Rust tests check homodimer and heterodimer rates exactly for counts 0–50 with coefficients 1 and 0.5; availability before arithmetic; tiny coefficients with overflowing standalone binomials; Hill threshold/extreme arguments; MM overflow avoidance; underflow/overflow diagnostics; immutable compilation; i64 count updates above 2^53; atomic failures; and both sides of conservation species/reaction caps and the dependency-edge cap.

The source mutation harness copies the checkout into a temporary directory. Each mutant must compile and fail named network tests: missing binomial divisor, wrong product sign, missing custom-rate read, and missing availability gate. An infrastructure or compile failure is not counted as a caught mutant.

## Development failures retained

- `development-tests.log`: the first network test run rejected a manually entered steep-Hill expected value. The 100-digit independent Decimal calculation in `reference-values.json` established the corrected reference, `0.44471580577189446`; the implementation returned `0.4447158057718944`. Only the incorrect reference literal changed. The predeclared tolerance was not relaxed. The independent oracle includes this case.
- `development-stable-clippy.log`: stable 1.99.1 rejected an indexed row-elimination loop as needless_range_loop. The loop now uses disjoint row slices and zipped iterators. All three toolchains and the independent oracle passed afterward.
- Intentional failed source mutations are retained in `mutations.json` with actual failed-test names. No failed numerical rows were excluded from the oracle report.

Log normalization removes only surplus blank lines at EOF for patch whitespace hygiene; diagnostic/test content is retained.

## Delivery verification and limits

The source candidate is applied to a fresh baseline checkout with `core.autocrlf=true` and validated independently; see `fresh-checkout.json` and `fresh-checkout-*.log/json`. The final mail patch is checked/applied separately and compared with the final commit tree; its verification is included with the downloadable delivery.

This iteration does not implement time integration or sample stochastic trajectories. It provides no SSA/FSP accuracy, throughput, cross-platform identity or biological validation claim. Conservation uses checked i128 and a 64-species/256-reaction cap. Its nonnegative-basis heuristic does not find every possible bound. Rates use f64; integer counts remain exact, but conversion to f64 can round above 2^53. Native Windows/macOS/ARM execution was not tested; LF behavior was tested using Git's autocrlf setting on Linux.

## Repository review, 1 November 2026

The delivered mail patch `0001-m1a-network-analysis.patch` (SHA-256 `f6d12484…`) applied cleanly to `8134f51` and produced tree `2301d8c`, as its instructions expect. The delivery's checks reproduced on Linux x86_64 with Rust 1.97.0 and Python 3.13.18. There were 73 release tests, 48 contract checks and 218 oracle checks, and the contract and oracle reports were byte-identical to `contracts.json` and `network.json` here. Four of four mutants were caught.

The review added independent probes and tests:

| Check | Result |
|---|---|
| Rates against 80-digit mpmath | 2943 propensities from 20 analyze runs: mass action up to order 60 and counts near 2^62, coefficients down to 1e-300; Hill with n up to 60 and k near x; Michaelis–Menten with km from 1e-100 to 1e20. Maximum relative error 6.1e-14; none above 2e-12 |
| Mutation probe of 20 rate, network, conservation and analysis rules | 17 caught by the delivered tests. Survivors: Hill without `log1p` near x = k, a bound not divided by its weight, mass action without binomial symmetry. Each now has a test, and `mutation_check.py` runs all seven mutants after requiring a passing unmutated suite |
| Bound coverage on 300 random closed conversion networks | Of 875 species whose bound a small integer combination (multipliers −3…3) of the reported basis proves, the analyzer reported 810. The 65 missing bounds (7.4%) come from mixed-sign basis rows that no single non-negative row can repair. This matches the documented limit of the non-negative-row heuristic |

After the review: 76 release tests (35 contract, 33 network, 8 CLI), 48 contract checks, 218 oracle checks and 7 of 7 mutants caught. GitHub Actions runs the full workflow on the pull request.
