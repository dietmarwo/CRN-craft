# Iteration 2 / M1a: scope and acceptance, declared before implementation

Date: 1 November 2026. Base: `8134f510546c4b0711b1e5c21855f18d9b53e80a` (merged PR #1).

Implement immutable indexed networks; sparse reactant/product/net columns;
reactant-availability gates; mass-action, Hill activation/repression and
Michaelis–Menten propensities; dependency graphs; checked transactional reaction
updates; and `analyze MODEL [STATE]`. Add exact conservation analysis for small
networks, with explicit resource/overflow outcomes. Linkage classes, expressions,
exporters and all time integration remain later work.

The existing model/run/check/result shapes stay at 0.1. New named-state and
analysis output contracts start at 0.1. Discovery advances to
`crn-craft.describe/0.2` when it advertises the new command and capabilities.
Simulation remains unavailable. An arbitrary supplied evaluation state is not
asserted to be reachable from the model's initial state.

Predeclared acceptance:

- Homodimer and heterodimer propensities: exact integer/half-integer results
  for counts 0–50 and coefficients 1 and 0.5, including unavailable reactants.
- General mass action: compare with Python's exact integer `math.comb` and
  rational coefficients, including tiny coefficients with large combinations.
- Hill/MM: compare with independently evaluated high-precision Decimal
  formulas. Relative tolerance 2e-12, plus 8 minimum subnormals absolute;
  independently checked zero/unrepresentable results use explicit rules.
- 200 seeded small networks, including conservative conversion networks:
  exact sparse stoichiometry, reads and dependencies against a direct Python
  oracle; exact null-space residuals, dimension and rational span/rank checks.
  Known enzyme laws must include the span of E+ES and S+ES+P.
- Exact conservation uses checked i128 arithmetic and a fixed dense-size cap.
  Exceeding either must return a labelled unavailable result, never a partial
  basis presented as complete. Non-negative laws can prove initial-class
  count bounds; signed laws alone cannot.
- Invalid state length/names/counts, safety limits, unavailable and zero-rate
  channels, invalid reaction indices, overflow and failed updates are tested.
  Failed updates must leave the entire state unchanged.
- Mutations removing a combinatorial divisor, changing a stoichiometric sign,
  dropping custom-rate reads or skipping reactant gates must be caught by tests.
- Existing review regressions, strict JSON and LF checkout behavior remain.
  fmt/clippy/release tests on 1.85, 1.88 and installed stable; independent Python
  checks normally and under -O; generated schemas and dependency gate; fresh
  patch application and validation.

Use pinned pure-Rust libm for transcendental/scaling operations. Finite positive
mathematical rates which overflow or round entirely to zero must fail explicitly;
no infinite rates or silently absorbing underflowed channels are accepted.
Numerical evaluation still uses f64 and does not certify correctly rounded
propensities or cross-platform identity. This is not an SSA accuracy or speed
campaign. Preserve failed checks and fixes in the implementation log.
