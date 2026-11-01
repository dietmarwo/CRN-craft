# Results and claims

`check` emits a `crn-craft.check/0.1` preflight report, not a simulation result. `valid: true` means no input diagnostics of severity error; `execution_available` is always false. `scope` is structural/domain validation. Raw model/run byte hashes identify the checked input snapshots.

The `crn-craft.result/0.1` Rust type/schema is a **reserved envelope only**, to make termination, assessment and evidence distinct from the start. It has no emitter, numerical acceptance checker or archive reader in this iteration. Do not manufacture instances as experiment evidence. M2 defines and versions the complete result/archive contract before use.

Planned termination states distinguish completion, event limit, wall-time limit, population limit, output limit, cancellation and failure. Planned assessments distinguish a fully completed exact-method ensemble from incomplete evidence, failed checks and an FSP approximation. Exact-method simulation still has Monte Carlo error and floating-point limitations.

Future summary logic must check relationships between termination, completed/requested work, missing replicates and evidence. Successful execution alone cannot justify a numerical assessment. State-dependent truncation can bias the completed subset. Report censored/incomplete paths and avoid presenting their conditional mean as an unbiased full-ensemble result.

The planned archive preserves exact consumed input bytes, numerical payloads, work, build/environment identity and hashes. Normalized inputs have separate hashes. A completion marker is written last. Timings/thread count/provenance are not part of numerical byte-identity comparisons. No archive durability or reproducible random stream is claimed in M0.

## M1a analysis reports

`analyze` emits a separate `crn-craft.analysis/0.1` contract with structure, evaluated propensities and bounded exact conservation. `success` describes this analysis; `simulation_available` is false. A conservation size/overflow warning leaves successful structure/rate results available and explicitly labels conservation unavailable. Neither `check` nor `analyze` emits the reserved simulation result envelope. See `docs/contracts.md` and `docs/network-api.md`.

## M1b structural and export reports

Current analysis is 0.2 and includes cone, graph and effective bounds. `analyze-v1` keeps 0.1 output. Export has a separate 0.1 report containing generated Python source and declared target limitations. Neither command executes CRN-craft trajectories or emits a simulator result envelope.

## M2 execution supersedes the reserved envelope

The sections above describe the historical M0/M1 reports. M2 now emits
`crn-craft.result/0.2`; the never-emitted 0.1 schema remains `result-v1`.
Direct SSA and archive verification are implemented. Read
[archive-contract.md](archive-contract.md) for the current payload, accounting,
completion/assessment relationships, verification and migration details.
