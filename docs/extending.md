# Extending the project

1. State the supported mathematical/contract scope and the milestone. Read the revised plan and create a model card before implementing a numerical algorithm.
2. Add typed input fields/variants with strict unknown-field handling, validation, codes and field paths. Update a delivered schema version when its shape changes. Generate snapshots from `describe`; do not hand-edit JSON Schema.
3. Add an independently justified acceptance fixture, including failure cases. A parser change must retain duplicate-key, numeric and depth checks. A rate change must include unavailable-reactant behavior and dependency reads.
4. Keep numeric compilation/algorithms in the library; the CLI handles files and JSON output. Python and job adapters will call the same core.
5. Run the documented checks and retain evidence with toolchain identity. Never introduce a reference simulator as a core dependency.

For M1 rate implementations, add exact small-count mass-action checks and finite-domain/overflow tests. For M2 Direct, add sampler vectors, analytic marginal tests, save-time boundaries, absorbing states and budget/cancellation tests. For later aggregators, reuse the full independent correctness suite before interpreting performance.

The dependency policy is intentionally closed: a new transitive package fails `check_dependencies.py` until reviewed and listed. Check package sources and `links` metadata before updating the allowlist. This gate is a review aid, not a comprehensive supply-chain audit.
