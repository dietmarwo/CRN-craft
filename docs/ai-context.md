# Agent context: CRN-craft M1a

CRN-craft models well-mixed chemical reaction networks. This version supports strict model/run validation, capability discovery, compiled network operations, propensity evaluation and bounded exact conservation analysis. No command runs SSA, FSP, ODEs or optimization.

Start with `crn-craft describe`; discovery is now `crn-craft.describe/0.2`. Use `describe model`, `describe run`, `describe state` and `describe analysis` for canonical schemas. The original model/run/check/result 0.1 shapes are unchanged. Select an example from `models/examples`. After editing a model, compute SHA-256 on its exact bytes and update its run reference; keep LF line endings. Keep the model in or below the run directory.

Run `check MODEL RUN` for structural and domain checks; it does not evaluate rates. Run `analyze MODEL [STATE]` for indices, sparse stoichiometry, propensity invalidation lists, rates and conservation status. State files name every species and may specify an arbitrary valid count vector. Initial-class conservation bounds are not bounds on arbitrary supplied states or a reachability certificate. An analysis success with `size_limit` or `arithmetic_overflow` still has no conservation result; inspect the status and diagnostics.

Identifiers use ASCII letters/digits/underscore, beginning with a letter or underscore, up to 64 bytes. Species/parameter names must be disjoint. Stoichiometry entries are positive integers; omit zeros. Parameter references are exact strings. Mass-action coefficients are stochastic count coefficients; Ω does not rescale them. Use counts and seconds. Hill/MM rates are complete propensities, gated by reactant availability. Read `docs/model-cards/rate-laws.md` and `docs/network-api.md` before translating input or building a scheduler. Dependencies invalidate propensities only; a future scheduler must separately reschedule a fired channel.

Read `docs/diagnostics.md` for stable codes and RFC 6901 pointers. Syntax/type errors stop parsing; semantic checks collect errors in stable order. Compiling directly constructed Rust models still validates them. Count updates are checked and atomic. Never invent a result for an unavailable solver.

Useful current prompt:

> Adapt an example and evaluate its rates at initial counts and one named state. Explain the stochastic count convention, availability gates, conservation status and the difference between initial-class bounds and arbitrary evaluation states. Preserve the model/run byte hash. State that no simulation has been performed.

Next implementation prompt:

> Implement the remaining M1b scope from the revised plan: linkage classes and explicit reference exporters with pinned reference versions and coefficient conventions. Add expression nodes only if a concrete fixture needs them. Preserve the M1a independent oracle, mutation catches, strict schemas and existing review regressions. Keep simulation unavailable until M2's Direct SSA acceptance gate is implemented.
