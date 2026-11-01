# M1b: conservation cones and reaction graphs

`analyze MODEL [STATE]` now emits analysis 0.2. It adds `nonnegative_conservation`, `linkage` and effective `initial_class_upper_bounds` inside `analysis`. The previous `conservation` block remains the linear-basis result, including its original basis-only bounds. `analyze-v1` preserves the complete 0.1 output and warning policy.

## Nonnegative conservation

The new calculation enumerates the primitive extreme rays of the cone `{l >= 0 : N^T l = 0}`. On `exact`, those rays generate all nonnegative linear conservation laws. Each row has integer-string coefficients and an exact initial total. For species i, the minimum `floor(l^T x0 / l_i)` over rays with positive l_i gives the tightest bound obtainable from such a law. A nonnegative combination cannot improve the minimum: its ratio is a weighted average of the contributing ray ratios.

This closes a gap in using only nonnegative rows of a rational nullspace basis. Consider the shipped fixture:

```bash
crn-craft analyze models/analysis/cone-gap.model.json
```

Its reactions are X2 → X0 + 2X1 and X3 → 2X2 + X0; initial counts are (0, 0, 0, 1). The previous basis rows are mixed-sign and prove no individual count bounds. The cone exposes these two primitive rays:

| Ray | Initial total | Species bounds supplied |
|---|---:|---|
| (0, 1, 2, 4) | 4 | X1 ≤ 4, X2 ≤ 2, X3 ≤ 1 |
| (1, 0, 1, 3) | 3 | X0 ≤ 3, X2 ≤ 3, X3 ≤ 1 |

Effective bounds are (3, 4, 2, 1), so current analysis can suppress the corresponding structural growth warnings. Legacy analysis keeps its previous conservative warnings.

The implementation starts with orthant coordinate rays and intersects reaction hyperplanes, retaining zero rays and combining adjacent opposite-sign rays. Support bitsets implement combinatorial adjacency; normalization and residual checks use checked integer arithmetic. This is an independent implementation of the standard double-description construction and adjacency criterion described in Burton, [Algorithm 3.1 and Lemma 3.4](https://arxiv.org/pdf/0808.4050). No source code from Regina or cddlib is included.

## Limits and result interpretation

| Limit or status | Behavior |
|---|---|
| At most 64 species / 256 reactions | Larger inputs return `size_limit` for the cone |
| At most 4,096 intermediate rays | Exhaustion returns `ray_limit` |
| At most 2,000,000 deterministic work units | Exhaustion returns `work_limit` |
| Checked i128 coordinates, products, sums and totals | Overflow returns `arithmetic_overflow` |
| Exact residual self-check | A nonzero residual returns `verification_failed` |

No unavailable outcome returns partial rays or cone-derived bounds. The independent linear-basis result remains available, and effective bounds fall back to its sound bounds. An unavailable cone does not imply no conservation laws. Ray enumeration can grow exponentially; the species cap alone does not guarantee completion.

Work counts ray visits, sparse dot-product terms, candidate pairs, adjacency probes, combined coordinates, verification terms and initial-total terms. Fixed dimensions bound normalization and map operations; this is a deterministic operation budget, not a wall-time guarantee. `analyze_cone_with_limits` allows smaller ray/work limits; values above the hard caps are clamped. Public reports expose the applied limits and consumed work.

On `exact`, a null bound means no nonnegative linear conservation law bounds that species. It does **not** prove dynamical unboundedness. Rate zeros, reachability restrictions or other nonlinear invariants can still constrain it. All reported bounds concern the initial conservation class; an arbitrary supplied evaluation state need not lie in that class. f64 rates, exact i64 updates and CRN-E042 underflow policy remain unchanged.

## Reaction graph

Complexes are the complete reactant/product multisets, with coefficients and catalysts retained. Their indices follow first occurrence while visiting reaction sides in declaration order. The empty complex is included when it appears in a source, sink or empty null reaction. Isolated declared species do not create a complex. Duplicate channels share complex nodes but keep one edge per reaction. Zero-rate and gated channels remain in this structural graph.

The graph reports weak linkage classes, strong linkage classes, terminal strong classes and weak reversibility. Components and member indices are sorted. Iterative traversal avoids recursion depth failures. Deficiency is `number_of_complexes - number_of_linkage_classes - exact_rank(N)` and is null when the exact rank is unavailable. The public `analyze_linkage(network, exact_rank)` API requires a previously established rank, or None; current CLI analysis obtains it from checked linear conservation analysis.

These terms follow the [Catalyst network-analysis definitions](https://docs.sciml.ai/Catalyst/stable/api/network_analysis_api/). No deficiency-zero theorem, equilibrium, stability or mass-action assumption is inferred from a structural report. Custom rates are allowed in graph analysis.

## Independent validation

The Python oracle enumerates every species support for small networks, computes its Fraction nullspace, and accepts exactly positive one-dimensional supports. This differs from the Rust incremental cone algorithm. It compares complete primitive ray sets, totals, bounds and exact residuals. A separate transitive-closure oracle checks all graph components and deficiency on 300 seeded networks and targeted cases.

Rust tests additionally cover both sides of size/ray/work caps (including a caller cap below the number of initial orthant rays), gcd reduction of combined rays to primitive laws, linkage classes joined only through converging edges, checked arithmetic and total overflow, redundant rays, weighted bounds, duplicate/reversed equations, empty/open/null networks and a 4,096-complex chain. Deliberately corrupted outputs and compiled source mutants must fail after a passing unmutated baseline. See `validation/iteration-03/RESULTS.md` for the executed evidence.
