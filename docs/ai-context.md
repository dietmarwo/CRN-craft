# Agent context: CRN-craft M1b

Current version: 0.1.0-dev.3. Capabilities: strict model/run validation; immutable compiled networks; all four propensity laws; exact bounded linear/cone analysis; reaction linkage graphs; elementary GillesPy2 reference export. No CRN-craft SSA, FSP, ODE or optimization execution exists.

Start with `crn-craft describe` (discovery 0.3). `check MODEL [RUN]` keeps its structural/domain semantics. `analyze MODEL [STATE]` emits analysis 0.2; `analyze-v1` preserves 0.1. State, model, run, check and reserved-result contracts remain 0.1. Use `describe analysis`, `describe analysis-v1` and `describe export` for generated schemas.

Current analysis has a linear-basis `conservation` block and an independent `nonnegative_conservation` cone. Prefer the effective `analysis.initial_class_upper_bounds`: complete cone bounds when exact, otherwise sound basis bounds. Always inspect cone status. Size/work/ray/overflow outcomes never imply no laws or unboundedness. Bounds refer to the initial class, not arbitrary evaluation states. Graphs retain catalysts, empty complexes and zero-rate channels; deficiency is structural and does not imply stability.

Read `docs/network-api.md`, `docs/m1b-analysis.md`, `docs/model-cards/rate-laws.md` and `docs/reference-export.md`. Count updates remain exact i64 and atomic. Propensities use f64, gated before arithmetic. Mass-action coefficients include combinatorial factorials; custom laws are complete propensities. Omega does not automatically rescale coefficients. Dependencies invalidate propensities only; future schedulers must separately reschedule a fired channel. Positive-to-zero rate underflow remains CRN-E042.

Model/run hashes identify exact bytes. Keep LF line endings, recompute run hashes only after intentional model edits, and keep model paths below the run directory. Never silently discard unknown fields, duplicate keys, safety limits or unsupported target features.

Useful current prompt:

> Analyze models/analysis/cone-gap.model.json and explain why exact cone rays find bounds missed by a mixed-sign nullspace basis. Compare analyze with analyze-v1. Export the elementary enzyme model to pinned GillesPy2 and verify its stochastic rate convention, keeping CRN-craft simulation marked unavailable.

Next implementation prompt:

> Implement M2's first Direct SSA baseline with a predeclared numerical acceptance plan, pinned samplers, event/time/population/output budgets, deterministic ensembles and verifiable archives. Preserve current rate underflow errors in the initial strict solver. Test analytic marginals, save-time and absorption behavior, thread-count numerical identity and failed/interrupted archives. Retain all existing oracle/mutation/review regressions. Do not infer solver correctness or speed from the reference-export smoke test.
