# Agent context: CRN-craft

Start with README.md, [architecture](architecture.md), `crn-craft describe`
(`crn-craft.describe/1.0`) and the [contracts](../reference/contracts.md). Release 0.7.0
merged the former network engine and rules engine into one; its
[implementation log](../campaign/merge-0.7/README.md) records every step and
deviation. Older campaign logs describe the engines as they were then.

## The one engine

- One crate, `crn-craft`; the engine is the module `crn_craft::rules`. The
  crate `crn-craft-rules` is a deprecated facade (`pub use
  crn_craft::rules::*`) that keeps 0.6.0 paths compiling. The command line is
  `crn-craft-cli`; Python is `bindings/python` plus `python/crn_craft`.
  Default engine MSRV 1.85. `integrations/jobs` (1.88), `integrations/diffsol`
  (1.89) and `validation/parameter-optimization/worker` are separate
  workspaces with their own lockfiles; keep them out of the default
  dependency gate.
- One model format: `crn-craft.model/1.0`. Attribute-free types are
  ordinary species; attributed types, patterns with `where` and `factor`,
  observables, whole-rule `propensity` expressions, continuous variables and
  a pattern `count` for higher-order reactions are in the
  [model reference](../reference/model-language.md). A network `crn-craft.model/0.1`
  is read only through the converter (`crn-craft convert`,
  `crn_craft.convert`).
- One request and result: `crn-craft.run/1.0` and `crn-craft.result/1.0`;
  `fsp` requests give `crn-craft.fsp-result/1.0`. The `rules-model` 0.1–0.4
  and `rules-run` 0.1–0.6 documents are still read and replay byte for byte.
  A `model/1.0` model needs a `run/1.0` request.
- Methods: `tree` (the default), `static_direct`, `next_reaction`,
  `mtdna_direct`, `thinning`, `tau`, `binomial`, `hybrid`, and `fsp`.
  Static Direct has the selectors `linear`, `composition_rejection`,
  `sorting`, `rssa_cr` and `auto` (the `run/1.0` default)
  ([exact methods](../reference/methods.md)). One random number generator: ChaCha8,
  one stream per replicate.
- Features on the flat view of a closed model (`crn_craft::rules::FlatView`):
  finite state projection, analysis (`analysis/1.0`), ODE kinetics for
  `crn-diffsol`, GillesPy2 export. Models without a closed set of species are
  refused with `CRN-R017`.
- Python: `import crn_craft as crn` writes the 1.0 contracts, asks for the
  quantiles 0.025, 0.5 and 0.975 and leaves the selector to `auto`.
  `crn_craft.rules` is the 0.6.0 API with its defaults (the `rules-*`
  versions, the linear scan, no quantiles). `Result.save` writes markdown or
  `crn-craft.report/1.0`. Parameter sweeps use `bind_parameters`
  ([parameter sweeps](../reusable-parameters.md)).

## Where things are

- Engine sources: `crates/crn-craft/src/rules/` (model, expressions,
  lowering, methods, selectors, flat view, converter, generators) and
  `crates/crn-craft/src/` (FSP, analysis, kinetics, export, archives, strict
  JSON, diagnostics, samplers, discovery and check).
- Tests: `crates/crn-craft/tests/`. Stored requests under `tests/stored/`
  (invariant I1) and generated-code fixtures under `tests/generated/` (I3);
  never edit either by hand.
- Schemas: `schemas/`, written from `crn-craft describe KIND`.
- Validation scripts: `validation/`. The merge gates, the 0.6.0 references
  and their declarations: `validation/merge/` (`RESULTS.md`).
- The frozen 0.6.0 user manual, run in CI: `validation/compat/manual-0.6.0.md`.

## Rules for changes

- Strict input everywhere: unknown fields and duplicate keys are errors;
  diagnostic codes are stable. A new code takes the next free number and goes
  into [diagnostics](../reference/diagnostics.md).
- Schemas come from Rust types. A change to a delivered contract's shape
  needs a new version, and the older versions keep their meaning.
- A new feature goes behind `model/1.0` or `run/1.0`; the stored `rules-run`
  requests must keep replaying byte for byte
  (`stored_requests.rs`, `static_direct_defaults.rs`).
- A new numerical policy or selector gets its own sampler identity and an
  exact-law check before any timing. `validation/merge/check_laws.py` runs
  the law gate against the committed references.
- Declare targets, sizes, seeds and tolerances before running; keep failed
  rows; never loosen a rule after seeing a result.
- Profile with Static Direct, then the generated event loop, before adding a
  specialized kernel ([checklist](specialized-kernels.md)).
- Run `cargo fmt`, clippy with `-D warnings`, `cargo test --release
  --workspace --locked`, rustdoc with `-D warnings`, ruff, and the
  validation scripts that CI names, each normally and under `python -O` with
  identical reports.

Approximate methods (`tau`, `binomial`, `hybrid`) carry their own sampler
identities; FSP reports leakage and error terms and never renormalizes.
Importers accept explicit subsets and conventions only. Incomplete ensembles
withhold statistics; wall time and cancellation are outside fixed-work
reproducibility. A successful run is not biological validation.

The Python package is on PyPI (`pip install crn-craft`) and the crates
`crn-craft`, `crn-craft-rules` and `crn-craft-cli` (`cargo install
crn-craft-cli --locked`) are on crates.io; the diffsol integration and the
job worker build from the repository.

## Optimization lessons

Use [C1–C4](../../tutorials/optimization/README.md) for runnable calibration,
optional BDF and SQP fitting, MODE and MAP-Elites oscillator searches, and
bounded rules inference. They call the library and the command line. Their
SSA runs use Static Direct with its default selector; C1's
`--include-next-reaction` adds a second exact method with its own stream.
[Validation](../../tutorials/optimization/VALIDATION.md) separates completed
checks from fitted-model quality.
