//! Exact stochastic simulation of dynamic, attributed species.
//!
//! **Status: experimental.** The rules engine of
//! [CRN-craft](https://github.com/dietmarwo/CRN-craft) is validated against
//! exact laws (master equations and Poisson laws, see the
//! [K1 gates](https://github.com/dietmarwo/CRN-craft/blob/main/validation/rules/RESULTS.md)),
//! but its language is still small and its contracts may gain versions. The
//! `tau`, `binomial` and `hybrid` methods are approximate and say so in their
//! sampler identity; `hybrid` integrates continuous variables between events
//! at the request's tolerances. Generated kernels and event loops are Rust
//! source; from Python, `crn_craft.rules.compile_loop` builds a generated
//! loop against the `crn-craft` crate and runs it as a separate program.
//!
//! A rules model declares species *types* with numeric attributes. An
//! *instance* is a type plus one value per attribute and a copy number.
//! Instances appear when a rule produces a new attribute vector and disappear
//! when their count reaches zero, so the species need not be listed in
//! advance: mitochondrial DNA clones, polymers with a length or clusters with
//! a size. Rules are schemata over zero, one or two reactant patterns. A
//! rule's propensity is a global rate over parameters and observables times
//! local factors over each matched instance's attributes. Observables reduce
//! the instances of one type to a number; an optional stop predicate over
//! them ends a replicate early. Models whose species are fixed can use the
//! network engine in [`crn_craft`] instead.
//!
//! # Quickstart
//!
//! Prepare a model and a run request once, then execute the ensemble. Both
//! are JSON documents; the request binds the model by the SHA-256 of its exact
//! bytes.
//!
//! ```
//! use crn_craft_rules::{Control, Termination, execute, prepare};
//! use serde_json::json;
//!
//! // Immigration at rate 10 and death at rate 1 per copy, from zero copies.
//! let model = serde_json::to_vec(&json!({
//!     "schema": "crn-craft.rules-model/0.1",
//!     "name": "immigration_death",
//!     "parameters": {"k_in": 10.0, "k_out": 1.0},
//!     "types": [{"name": "X"}],
//!     "initial": [],
//!     "observables": {"X": {"type": "X"}},
//!     "rules": [
//!         {"name": "immigration", "reactants": [], "rate": "k_in",
//!          "products": [{"type": "X"}]},
//!         {"name": "death", "reactants": [{"label": "x", "type": "X"}],
//!          "rate": "k_out", "products": []}
//!     ]
//! }))?;
//! let request = serde_json::to_vec(&json!({
//!     "schema": "crn-craft.rules-run/0.1",
//!     "model": {"path": "immigration-death.rules.json",
//!               "sha256": crn_craft::sha256(&model)},
//!     "time": {"end": 5.0, "report_times": [0.0, 1.0, 5.0]},
//!     "ensemble": {"replicates": 200, "root_seed": 7},
//!     "outputs": {"statistics": ["mean", "variance"],
//!                 "keep_paths": 1, "keep_instances": 0},
//!     "budget": {"max_events_per_replicate": 100000,
//!                "max_sampling_steps_per_replicate": 1000000,
//!                "max_instances": 1, "threads": 2,
//!                "max_output_bytes": 1048576}
//! }))?;
//!
//! let prepared = prepare(&model, &request).map_err(|d| format!("{d:?}"))?;
//! let result = execute(&prepared, &Control::default()).map_err(|d| format!("{d:?}"))?;
//!
//! assert_eq!(result.termination, Termination::Completed);
//! assert_eq!(result.observables, ["X"]);
//! // Rows are report times, columns observables. The exact law at t = 5 is
//! // Poisson with mean 10 (1 - e^-5), about 9.93.
//! let mean = result.statistics.as_ref().and_then(|s| s.mean.as_ref()).unwrap();
//! assert!((mean[2][0] - 9.93).abs() < 1.0);
//! assert_eq!(result.paths[0].values.len(), 3);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Entry points
//!
//! - [`check_model`] and [`check_model_file`] validate a model and return a
//!   [`RulesCheck`] report.
//! - [`prepare`] validates a model and a [`RunRequest`] and returns a
//!   reusable [`Prepared`] run. [`execute`] runs it as often as needed;
//!   [`Prepared::bind_parameters`] derives candidates for parameter sweeps.
//! - [`solve`] and [`solve_files`] prepare and execute in one call.
//! - [`execute_with_observer`] calls an [`observe::EventObserver`] for every
//!   committed event on the worker threads. Trackers and event logs
//!   ([`observe::TrackerSpec`], [`observe::LogSpec`]) summarize events inside
//!   the [`RulesResult`] instead.
//! - [`generate_rules_kernel`] with [`execute_with_rules_kernel`], and
//!   [`generate_rules_loop`] with [`execute_with_rules_loop`], emit and run
//!   ahead-of-time Rust code for a Static Direct model.
//! - [`schema`] returns the JSON Schemas of the contracts; [`describe`]
//!   returns the discovery entry of the native CLI.
//!
//! The model types live in [`model`], observation in [`observe`], the ATP
//! coupling of rules-model/0.2 in [`atp`] and the continuous variables of
//! rules-model/0.3 and 0.4 with the `hybrid` method in [`continuous`].
//!
//! # Contracts
//!
//! Every document names its contract in its `schema` field. Unknown fields
//! and duplicate keys are rejected, and a present optional field must hold a
//! value (`null` is rejected).
//!
//! | Contract | Rust type | Content |
//! | --- | --- | --- |
//! | `crn-craft.rules-model/0.1` | [`Model`] | Types, parameters, initial instances, observables, rules and stop |
//! | `crn-craft.rules-model/0.2` | [`Model`] | 0.1 plus the [`atp`] block |
//! | `crn-craft.rules-model/0.3` | [`Model`] | 0.2 plus whole-rule `propensity` expressions ([`model::Rule::propensity`]) and [`continuous`] variables |
//! | `crn-craft.rules-model/0.4` | [`Model`] | 0.3 plus `independent` continuous variables ([`continuous::ContinuousSpec::independent`]), integrated once per prepared run by `hybrid` |
//! | `crn-craft.rules-check/0.1` | [`RulesCheck`] | Validity, model hash and diagnostics |
//! | `crn-craft.rules-run/0.1` | [`RunRequest`] | Time grid, ensemble, outputs and budget; always the `tree` method |
//! | `crn-craft.rules-run/0.2` | [`RunRequest`] | 0.1 plus `method`, `trackers` and `event_log` |
//! | `crn-craft.rules-run/0.3` | [`RunRequest`] | 0.2 plus the `mtdna_direct` method |
//! | `crn-craft.rules-run/0.4` | [`RunRequest`] | 0.3 plus the Static Direct options `waiting_time` and `sparse_trackers` |
//! | `crn-craft.rules-run/0.5` | [`RunRequest`] | 0.4 plus the `hybrid` method and its tolerances |
//! | `crn-craft.rules-run/0.6` | [`RunRequest`] | 0.5 with new Static Direct defaults: Ziggurat waits, and sparse integration of every eligible tracker |
//! | `crn-craft.rules-result/0.1` | [`RulesResult`] | Result of a rules-run/0.1 request |
//! | `crn-craft.rules-result/0.2` | [`RulesResult`] | Result of a rules-run/0.2–0.4 request, with `observations` |
//! | `crn-craft.rules-result/0.3` | [`RulesResult`] | Result of a rules-run/0.5 or 0.6 request: 0.2 plus `method` and `continuous` |
//!
//! # Methods
//!
//! [`Method`] selects the algorithm. Every method except `tree` first lowers
//! the model to a fixed network: the initial tuples plus the products with
//! constant attributes (at most 256 tuples), and at most 4096 channels. A
//! model whose products can leave that set, for example through random
//! attributes, `new_id()` or attribute arithmetic, runs only on `tree`;
//! `prepare` reports `CRN-R020` at `/method` otherwise. Every method samples
//! whole-rule propensities (rules-model/0.3); the static methods need each
//! of their reactant types to lower to one tuple, `mtdna_direct` falls back
//! to Static Direct with `CRN-R051`, and generated rates kernels refuse them
//! with `CRN-R020`. Generated event loops emit every rate, propensity and
//! observable that Static Direct runs.
//!
//! | `kind` | Exact | Scope |
//! | --- | --- | --- |
//! | `tree` (default) | yes | Every rules model without ATP or continuous variables; dynamic instances over sum trees |
//! | `static_direct` | yes | Closed models; Ziggurat waits and sparse trackers, opt-in in rules-run/0.4 and 0.5 and the default from 0.6 |
//! | `mtdna_direct` | yes | The fixed three-class vacancy model; falls back to Static Direct with `CRN-R051` |
//! | `next_reaction` | yes | Closed models, with indexed clocks ([DOI 10.1021/jp993732q](https://doi.org/10.1021/jp993732q)) |
//! | `thinning` | yes | Required by, and only for, rules-model/0.2 models with ATP |
//! | `tau` | no | Closed models; adaptive Poisson leaps ([DOI 10.1063/1.2159468](https://doi.org/10.1063/1.2159468)) |
//! | `hybrid` | no | Required by, and only for, closed models with [`continuous`] variables (also runs closed models without); Dormand–Prince 5(4) between events at the request's tolerances |
//! | `binomial` | no | Closed models with single-parent unit births and deaths |
//!
//! # Reproducibility
//!
//! Replicate `i` draws from its own ChaCha8 stream, [`stream`]`(root_seed,
//! i)`, and the ensemble is reduced in replicate order. The numbers of a
//! completed result therefore depend only on the exact model bytes, the
//! request's settings and seed, and the sampler identity in
//! [`RulesResult::sampler`]. They do not depend on `budget.threads`, which
//! changes only `run_sha256`. Timing, cancellation and wall-time limits can
//! change which replicates complete.
//!
//! # Diagnostics
//!
//! Errors are [`crn_craft::diagnostic::Diagnostic`] values with stable codes
//! and an RFC 6901 pointer into the offending document:
//!
//! | Code | Meaning |
//! | --- | --- |
//! | `CRN-E001`–`CRN-E004` | Malformed JSON, duplicate keys, unknown fields or wrong types, size limits |
//! | `CRN-E030` | A model or request file cannot be read |
//! | `CRN-R010`–`CRN-R015` | Invalid names, unknown types, expressions, rule shapes, initial values and model size or ATP limits |
//! | `CRN-R020` | Invalid run request, method, tracker or log, budget reservation or generated-code fingerprint |
//! | `CRN-R030` | A negative or non-finite rate, factor, weight or attribute, an undefined stop, or an overflow |
//! | `CRN-R031` | A positive weight below 2⁻⁵³ of its total, which a 53-bit uniform cannot select |
//! | `CRN-R040`–`CRN-R045` | Event, sampling, instance, wall-time, cancellation and output limits |
//! | `CRN-R050` | Warning: the ensemble is incomplete and statistics are withheld |
//! | `CRN-R051` | Warning: a requested specialization did not apply |
//!
//! # Related crates and documentation
//!
//! - [`crn-craft`](https://docs.rs/crn-craft): the fixed-network engine and
//!   the shared diagnostics, strict JSON parsing and contract types.
//! - [`crn-craft-cli`](https://crates.io/crates/crn-craft-cli): the native
//!   CLI, with `crn-craft rules-check`, `crn-craft rules` and
//!   `crn-craft describe rules-run`.
//! - [`crn-craft` on PyPI](https://pypi.org/project/crn-craft/): Python
//!   ensembles through `crn_craft.rules`.
//! - [Rules language](https://github.com/dietmarwo/CRN-craft/blob/main/docs/rules-language.md):
//!   the model language, propensities, expressions and numerical policy.
//! - [mtDNA ensembles](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-rules.md)
//!   and the [mtDNA user manual](https://github.com/dietmarwo/CRN-craft/blob/main/docs/mtdna-user-manual.md):
//!   methods, trackers, event logs, callbacks and later contract versions.
//! - [Rules kernel checklist](https://github.com/dietmarwo/CRN-craft/blob/main/docs/specialized-kernels.md):
//!   generated rules kernels and event loops.
//! - [Architecture](https://github.com/dietmarwo/CRN-craft/blob/main/docs/architecture.md)
//!   and [diagnostic codes](https://github.com/dietmarwo/CRN-craft/blob/main/docs/diagnostics.md).

#![warn(missing_docs)]

// This crate is a facade. Since 0.7.0 the rules engine is the module
// `crn_craft::rules` of the `crn-craft` crate, and this glob re-export keeps
// every path that compiled against crn-craft-rules 0.6.0: the items, the
// public modules `atp`, `continuous`, `model` and `observe`, and the hidden
// `libm` re-export that kernels and loops generated by 0.6.0 call. Since
// 0.7.0 generated code names `crn_craft::rules::` and needs only `crn-craft`;
// code generated by 0.6.0 keeps compiling through this facade.
// `tests/facade.rs` lists the 0.6.0 paths.
pub use crn_craft::rules::*;
