# Iteration 3 / M1b scope and acceptance, declared before implementation

Date: 1 November 2026. Base: GitHub main `c167add` (PR #2).

Implement nonnegative conservation-cone enumeration, linkage/strong-linkage
classes and structural deficiency, and one executable pinned reference exporter
for elementary stochastic mass action (GillesPy2 1.8.3). Expressions are still
unnecessary for the fixtures. No CRN-craft simulation is implemented.

Keep the reviewed M1a core and its seven-mutant baseline gate. Add a complete
extreme-ray description of {l >= 0 : N^T l = 0} when checked arithmetic and work
caps permit it. Use orthant/hyperplane double description with combinatorial
adjacency; exact primitive i128 rays, checked residuals/totals and deterministic
ordering. Exhausting size, ray, pair/adjacency-work or arithmetic limits must
report unavailable, never a partial cone labelled complete. Retain existing
linear-basis bounds as a sound fallback; prefer exact cone bounds when present.
The old linear-basis API and schema remain available for compatibility. A
complete cone proves all bounds obtainable from nonnegative linear conservation
laws; it is not a reachability or boundedness decision procedure for the CRN.

The new analysis shape is version 0.2; preserve the old 0.1 schema and provide
legacy analysis explicitly. Discovery moves to 0.3; export report starts at 0.1.
Model/run/check/result/state contracts and exact-byte example hashes stay intact.

Predeclared gates:

- Preserve 76 baseline Rust tests, 48 contract checks, 218 independent M1a checks
  and all seven source mutants; baseline must pass before mutation scoring.
- Compare exact cone rays with a separate support-enumeration/Fraction nullspace
  oracle on 300 seeded networks of at most 7 species. Check complete ray sets,
  primitive nonnegative coefficients, residuals, initial totals and all bounds.
  Include the review's mixed-sign-basis counterexample, weighted totals, open
  networks, disconnected/no-reaction cases, duplicate and reversed reactions.
- Verify both sides of size, ray and deterministic work caps. Arithmetic overflow
  must preserve the linear analysis, expose a diagnostic and emit no partial rays.
- Linkage/SCC/terminal classes match an independent Python reachability oracle on
  the same 300 networks. Pin empty/null/source/sink, catalyst and disconnected
  cases, weak reversibility and n_complexes - n_linkage - rank deficiency. No
  mass-action deficiency theorem is asserted for arbitrary custom kinetics.
- Test exported GillesPy2 models using the installed pinned package, not a stub:
  inspect stoichiometry and compare its generated propensities at counts 0-50
  against integer combinations and the Rust evaluator. Include sources, first
  order, homodimers, heterodimers, non-unit Omega, catalysts, zero coefficients,
  symbol collisions and exact parameter substitution. Declare float limitations.
  Reject custom rates, order >2, population safety limits and initial counts
  above 2^53 before producing source; export no false semantic equivalence.
- Exported code must compile/import/build, and one short reference SSA smoke run
  must complete with valid nonnegative counts. It is an exporter smoke test,
  not validation or benchmarking of an implemented CRN-craft simulator.
- New mutations must exercise cone pairing, adjacency, bound division, linkage
  direction/complex identity and homodimer coefficient conversion.
- fmt/clippy/release checks on 1.85, 1.88 and installed stable; Python checks
  normal/-O with identical fixed-work reports; dependency gate and LF checkout;
  actual fresh application/rebuild of the final patch. Retain failed rows.

Numerical comparisons retain the predeclared M1a tolerance: 2e-12 relative plus
8 minimum subnormals absolute. Exact combinatorial cases must be exact where
specified. No stochastic performance claim is made. Sources, package identity,
limits, failures and final evidence will be recorded in the iteration log.
