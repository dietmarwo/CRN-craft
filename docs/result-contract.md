# Results and claims

M0 emits a `crn-craft.check/0.1` preflight report, not a simulation result. `valid: true` means no input diagnostics of severity error; `execution_available` is always false. `scope` is structural/domain validation. Raw model/run byte hashes identify the checked input snapshots.

The `crn-craft.result/0.1` Rust type/schema is a **reserved envelope only**, to make termination, assessment and evidence distinct from the start. It has no emitter, numerical acceptance checker or archive reader in this iteration. Do not manufacture instances as experiment evidence. M2 defines and versions the complete result/archive contract before use.

Planned termination states distinguish completion, event limit, wall-time limit, population limit, output limit, cancellation and failure. Planned assessments distinguish a fully completed exact-method ensemble from incomplete evidence, failed checks and an FSP approximation. Exact-method simulation still has Monte Carlo error and floating-point limitations.

Future summary logic must check relationships between termination, completed/requested work, missing replicates and evidence. Successful execution alone cannot justify a numerical assessment. State-dependent truncation can bias the completed subset. Report censored/incomplete paths and avoid presenting their conditional mean as an unbiased full-ensemble result.

The planned archive preserves exact consumed input bytes, numerical payloads, work, build/environment identity and hashes. Normalized inputs have separate hashes. A completion marker is written last. Timings/thread count/provenance are not part of numerical byte-identity comparisons. No archive durability or reproducible random stream is claimed in M0.
