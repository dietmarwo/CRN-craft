# Compiled networks and analysis

The core provides `CompiledNetwork::compile(&Model)`, `analyze_bytes` and `analyze_files`. Compilation validates even directly constructed models, resolves parameter values and species references once, and owns an immutable snapshot. Editing the source model afterward does not alter the compiled network. Compile again for another parameter set. Compilation errors are diagnostics; successful compilation retains structural warnings.

## Indices, sparse columns and dependencies

Species and reaction indices are zero-based declaration order. Every reaction contains sorted sparse `reactants`, `products` and `net_changes` terms, each `{species, coefficient}`. Net coefficients are products minus reactants with zeros omitted. These are columns of the species-by-reaction stoichiometric matrix N. They are not rows or propensity multipliers.

`reads` includes every reactant gate and the Hill regulator or MM substrate. It is conservative even for zero rate coefficients. `dependencies[j]` is the sorted, unique list of reactions whose propensities can change when reaction j changes counts. It is built from nonzero net changes intersecting reads. A catalyst can gate a reaction without being changed by it. A null reaction has no propensity invalidations. A future event scheduler must separately reschedule its fired channel; these lists only describe rate invalidation.

Compilation caps the graph at 1,000,000 edges and returns CRN-E040 on overflow. Existing parser/model caps still apply. No edges are silently dropped.

## State evaluation and updates

`initial_counts`, `species`, `reactions`, `dependencies` and `warnings` return immutable views. `named_state` maps a complete name-to-count map to declaration order. `validate_state`, `propensity` and `propensities` reject invalid lengths, negative counts and exceeded declared limits. `fire_reaction(j, &mut counts)` requires an enabled channel, checks **every** proposed i64 update and population limit, then commits all changes. Any error leaves the state untouched. Counts above 2^53 remain exact in storage and updates.

These operations have no random generator, event selection, elapsed time, thread scheduling or simulation result. Propensity formulas use f64 and the limits in the [rate-law card](model-cards/rate-laws.md). A validated model may still have an unrepresentable rate at a particular state.

`analyze MODEL [STATE]` uses initial counts unless given a strict `crn-craft.state/0.1` document. Its `model_sha256` and `state_sha256` identify raw supplied bytes, not normalized JSON. The optional state need not lie in the initial conservation class. The output separates `initial_counts` and `evaluation_counts`; it makes no reachability claim.

## Linear-basis conservation (legacy-compatible)

The analyzer computes the rational left nullspace of N by checked integer row elimination with gcd reduction, followed by checked rational back substitution. Primitive integer basis rows satisfy lᵀN=0. A limited sequence of elementary basis combinations exposes nonnegative rows where available; it is **not** an enumeration of the nonnegative conservation cone. Each returned row is verified against the original sparse columns. Basis rows span the rational conservation space, not necessarily the full integer lattice of laws.

Coefficients and `initial_total = lᵀx₀` are decimal strings so i128 integers survive JSON consumers that use doubles. Species order matches the network. Rank means rank(N), and an exact basis has `species_count - rank` independent rows. Zero-reaction networks return the identity basis. Full-rank networks return an exact empty basis.

| Status | Meaning |
|---|---|
| `exact` | Complete rational basis, exact integer residuals and initial totals within checked i128 arithmetic |
| `size_limit` | More than 64 species or 256 reactions; no basis/rank computed, CRN-W010 |
| `arithmetic_overflow` | An intermediate checked i128 operation exceeds its range; no partial result retained, CRN-W011 |

Conservation size/overflow outcomes leave successfully compiled structure and rate evaluation available. `success` remains true, with an explicit status and diagnostics. They do not mean that no conservation law exists.

For each nonnegative basis row, every positive coefficient supplies `floor(initial_total / coefficient)` as an upper bound on that species **in the initial conservation class**. The minimum available bound is reported as an integer string. A signed row alone supplies no count bound. A `null` bound means “not proved by this basis,” not “unbounded.” Even a complete rational basis can miss bounds requiring other nonnegative combinations. Initial-class bounds do not constrain an arbitrary supplied evaluation state.

For the explicit enzyme example, the analyzer exposes E+ES and S+ES+P, yielding bounds on all four species. `analyze` suppresses the structural growth warning only for a species with a proven initial-class bound. `check` keeps its existing conservative warnings and performs no conservation analysis.

## Validation and next use

The independent oracle uses Fraction Gauss–Jordan elimination, a different algorithm from the Rust integer elimination. It checks rank, nullity, residuals, basis independence, rational span and bound totals on 200 seeded small networks. Tests also cover known enzyme laws, a weighted law (A + 2B), signed laws, no reactions, full rank, caps and arithmetic overflow. Seven source mutation tests, run only after the unmutated suite passes, exercise stoichiometric signs, combinatorial divisors, custom reads, availability gates, the `log1p` Hill path near x = k, weight division in bounds and binomial symmetry.

M1b adds a separate exact nonnegative-cone calculation, effective initial-class bounds, linkage analysis and a reference exporter; see [M1b analysis](m1b-analysis.md). The linear-basis fields above are retained unchanged for legacy consumers. M2 must add a separately validated Direct SSA scheduler, RNG stream design, live budgets and result/archive contracts; an atomic reaction update is only one building block.
