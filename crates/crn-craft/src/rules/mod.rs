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
//! loop against this crate and runs it as a separate program.
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
//! network engine in [`crn_craft`](crate) instead.
//!
//! # Quickstart
//!
//! Prepare a model and a run request once, then execute the ensemble. Both
//! are JSON documents; the request binds the model by the SHA-256 of its exact
//! bytes.
//!
//! ```
//! use crn_craft::rules::{Control, Termination, execute, prepare};
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
//! | `crn-craft.model/1.0` | [`Model`] | The unified model contract: means what rules-model/0.4 means; runs only with run/1.0 |
//! | `crn-craft.run/1.0` | [`RunRequest`] | The unified request contract: the fields and defaults of rules-run/0.6 plus quantile statistics ([`Outputs::quantiles`]); accepts model/1.0 and every rules-model version |
//! | `crn-craft.result/1.0` | [`RulesResult`] | Result of a run/1.0 request, with the fields of rules-result/0.3 plus [`Statistics::quantiles`] |
//! | `crn-craft.check/1.0` | [`RulesCheck`] | Check report of a model/1.0 document, with the fields of rules-check/0.1 |
//!
//! [`convert::from_network`] writes a network model (`crn-craft.model/0.1`)
//! as a `crn-craft.model/1.0` model.
//!
//! # Methods
//!
//! [`Method`] selects the algorithm. Every method except `tree` first lowers
//! the model to a fixed network: the initial tuples plus the products with
//! constant attributes (at most 256 tuples), and at most 4096 channels
//! (4,096 tuples and 65,536 channels for a `crn-craft.model/1.0` model). A
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
//! Errors are [`crn_craft::diagnostic::Diagnostic`](crate::diagnostic::Diagnostic) values with
//! stable codes and an RFC 6901 pointer into the offending document:
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

// Some contract fields and enum variants are documented with
// `#[doc = concat!(..)]` instead of `///`. schemars copies literal doc
// strings into the generated JSON Schemas, where a description would wrap a
// bare `$ref` in `allOf` or split a string enum into `oneOf` entries, which
// changes the committed schema shapes under `schemas/`. schemars ignores the
// `concat!` form; rustdoc and `missing_docs` treat it as a doc comment.

pub mod atp;
mod combination;
pub use combination::combination_rate;
pub mod continuous;
pub mod convert;
mod engine;
mod expr;
mod flat;
pub mod flat_view;
pub use flat_view::FlatView;
mod generated;
mod generated_loop;
pub use engine::{LoopHalt, LoopHost};
pub use generated::{RulesKernel, execute_with_rules_kernel, generate_rules_kernel};
pub use generated_loop::{RulesLoop, execute_with_rules_loop, generate_rules_loop};
/// The math library emitted stop predicates call, so that a generated loop
/// needs no dependency of its own.
#[doc(hidden)]
pub use libm;
pub mod model;
mod mtdna;
pub mod observe;
mod quantile;
mod tree;
use observe::{EventObserver, LogSpec, NoObserver, Observation, TrackerSpec};

use crate::contracts::{Ensemble, ModelRef, Statistic};
use crate::diagnostic::Diagnostic;
use rand_chacha::ChaCha8Rng;
use rand_core::SeedableRng;
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Condvar, Mutex};
use std::time::{Duration, Instant};

pub use model::Model;

/// Sampler identity of the `tree` method, reported in
/// [`RulesResult::sampler`] for every rules-run/0.1 request and for
/// `{"kind": "tree"}`.
///
/// Version 2 selects identical-pattern pairs exactly from pair trees.
/// Version 1, used only for the first exact-law gate run, drew two molecules
/// and rejected a draw of the same one; it is no longer produced.
pub const SAMPLER: &str = "chacha8-u53/rules-direct/2";
/// Largest `ensemble.replicates` a request may ask for; more is `CRN-R020`.
pub const MAX_REPLICATES: u64 = 1_000_000;
/// Largest number of entries in `time.report_times`; more is `CRN-R020`.
pub const MAX_REPORT_TIMES: usize = 100_000;
/// Largest `budget.threads`. The engine starts `min(threads, replicates)`
/// workers.
pub const MAX_WORKERS: u32 = 64;
/// Hard cap of 1 GiB on the compact JSON of a result, applied on top of
/// `budget.max_output_bytes`.
pub const MAX_RESULT_BYTES: u64 = 1 << 30;
/// Upper bound of 2 GiB on the reservation [`prepare`] computes for records,
/// in-flight and kept series, kept final instances and observations. It is a
/// conservative admission estimate (`CRN-R020` at `/budget` when exceeded),
/// not a cap on allocator or OS memory.
pub const MAX_WORKSPACE_BYTES: u64 = 2 << 30;
/// Largest number of entries in `outputs.quantiles` (`crn-craft.run/1.0`).
pub const MAX_QUANTILES: usize = 16;

/// The random stream of one replicate.
///
/// ChaCha8 keyed by the SHA-256 of the domain `crn-craft.rules/0.1` followed
/// by `root_seed` as eight little-endian bytes, with the replicate index as
/// stream id. Every method draws replicate `replicate` from this stream, so a
/// replicate can be replayed alone; the sampler identity in
/// [`RulesResult::sampler`] says how the words are used.
///
/// ```
/// use rand_core::RngCore;
///
/// let mut a = crn_craft::rules::stream(17, 3);
/// let mut b = crn_craft::rules::stream(17, 3);
/// assert_eq!(a.next_u64(), b.next_u64());
/// assert_ne!(
///     crn_craft::rules::stream(17, 4).next_u64(),
///     crn_craft::rules::stream(18, 3).next_u64(),
/// );
/// ```
pub fn stream(root_seed: u64, replicate: u64) -> ChaCha8Rng {
    let mut digest = Sha256::new();
    digest.update(b"crn-craft.rules/0.1");
    digest.update(root_seed.to_le_bytes());
    let mut rng = ChaCha8Rng::from_seed(digest.finalize().into());
    rng.set_stream(replicate);
    rng
}

/// Cooperative cancellation, by flag or by the existence of a file.
///
/// Workers poll the control when a replicate starts and every 1024 events.
/// A cancelled run ends its running replicates as
/// [`Termination::Cancelled`] (`CRN-R044`), records the replicates it never
/// started with `started: false`, and withholds statistics. The flag is never
/// cleared, so use a fresh `Control` for the next run. [`Control::default`]
/// is not cancelled and watches no file.
///
/// ```
/// # use crn_craft::rules::{Control, Termination, execute, prepare};
/// # use serde_json::json;
/// # let model = serde_json::to_vec(&json!({
/// #     "schema": "crn-craft.rules-model/0.1", "name": "immigration_death",
/// #     "parameters": {"k_in": 10.0, "k_out": 1.0},
/// #     "types": [{"name": "X"}], "initial": [],
/// #     "observables": {"X": {"type": "X"}},
/// #     "rules": [
/// #         {"name": "immigration", "reactants": [], "rate": "k_in", "products": [{"type": "X"}]},
/// #         {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "k_out", "products": []}
/// #     ]
/// # })).unwrap();
/// # let request = serde_json::to_vec(&json!({
/// #     "schema": "crn-craft.rules-run/0.1",
/// #     "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(&model)},
/// #     "time": {"end": 5.0, "report_times": [0.0, 5.0]},
/// #     "ensemble": {"replicates": 8, "root_seed": 1},
/// #     "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
/// #     "budget": {"max_events_per_replicate": 100000, "max_sampling_steps_per_replicate": 1000000,
/// #                "max_instances": 1, "threads": 1, "max_output_bytes": 1048576}
/// # })).unwrap();
/// let prepared = prepare(&model, &request).unwrap();
/// let control = Control::default();
/// control.cancel(); // for example from another thread or a signal handler
/// let result = execute(&prepared, &control).unwrap();
/// assert_eq!(result.termination, Termination::Cancelled);
/// assert!(result.statistics.is_none());
/// assert_eq!(result.diagnostics[0].code, "CRN-R050");
/// ```
#[derive(Debug, Default)]
pub struct Control {
    cancelled: AtomicBool,
    /// Cancel when this file exists. The check runs at the same polls as the
    /// flag; it reads file system metadata, so per-event callers should use
    /// [`Control::cancel_flagged`].
    pub cancel_file: Option<PathBuf>,
}
impl Control {
    /// Sets the cancellation flag. Safe to call from any thread while a run
    /// is in progress.
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }
    /// The in-memory flag only. Per-event callers use this to avoid a file
    /// system check on every event; the engine polls the file separately.
    pub fn cancel_flagged(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
    /// Whether the flag is set or the cancel file exists. Finding the file
    /// also sets the flag.
    pub fn is_cancelled(&self) -> bool {
        if self.cancelled.load(Ordering::Relaxed) {
            return true;
        }
        if self.cancel_file.as_ref().is_some_and(|p| p.exists()) {
            self.cancel();
            return true;
        }
        false
    }
}

// ---------------------------------------------------------------------------
// Run request

/// The version of a run request, its `schema` field.
///
/// Versions 0.2 to 0.5 each add fields to the previous one; using a field
/// before its version is `CRN-R020` at `/schema`. Version 0.6 adds no field:
/// it changes what an absent `static_direct` option means (Ziggurat waits,
/// sparse trackers where eligible), and a stored request of an earlier
/// version keeps its meaning. Version 0.1 returns rules-result/0.1, versions
/// 0.2 to 0.4 rules-result/0.2, and 0.5 and 0.6 rules-result/0.3.
#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq, PartialOrd, Ord,
)]
pub enum RunVersion {
    // Variant docs use `concat!` to stay out of the schema's string enum.
    #[doc = concat!(
        "`crn-craft.rules-run/0.1`: the `tree` method only, without ",
        "`method`, `trackers` or `event_log`."
    )]
    #[serde(rename = "crn-craft.rules-run/0.1")]
    V01,
    #[doc = concat!(
        "`crn-craft.rules-run/0.2`: adds `method`, `trackers` and `event_log`."
    )]
    #[serde(rename = "crn-craft.rules-run/0.2")]
    V02,
    #[doc = concat!(
        "`crn-craft.rules-run/0.3`: adds the `mtdna_direct` method."
    )]
    #[serde(rename = "crn-craft.rules-run/0.3")]
    V03,
    #[doc = concat!(
        "`crn-craft.rules-run/0.4`: adds the Static Direct options ",
        "`waiting_time` and `sparse_trackers`."
    )]
    #[serde(rename = "crn-craft.rules-run/0.4")]
    V04,
    #[doc = concat!(
        "`crn-craft.rules-run/0.5`: adds the `hybrid` method; returns ",
        "rules-result/0.3."
    )]
    #[serde(rename = "crn-craft.rules-run/0.5")]
    V05,
    #[doc = concat!(
        "`crn-craft.rules-run/0.6`: the fields of 0.5; Static Direct ",
        "defaults to Ziggurat waits and to sparse integration of every ",
        "eligible tracker. Returns rules-result/0.3."
    )]
    #[serde(rename = "crn-craft.rules-run/0.6")]
    V06,
    #[doc = concat!(
        "`crn-craft.run/1.0`: the unified request contract, with the fields ",
        "and defaults of 0.6 plus `outputs.quantiles`. It accepts ",
        "`crn-craft.model/1.0` and every rules-model version and returns ",
        "`crn-craft.result/1.0`."
    )]
    #[serde(rename = "crn-craft.run/1.0")]
    V10,
}

/// The simulated horizon and the report grid of a run request.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Time {
    /// Simulated time horizon, positive and finite, in the time unit of the
    /// model's rates.
    pub end: f64,
    /// Strictly increasing times in [0, end] at which observables are
    /// reported (right-continuous: events at a report time are included).
    /// Between 1 and 100,000 entries.
    pub report_times: Vec<f64>,
}

/// What a result keeps besides the per-replicate records.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Outputs {
    /// Ensemble statistics per report time and observable, each listed at
    /// most once: `mean`, and `variance` for the sample variance and the
    /// standard error of the mean. They are withheld unless every replicate
    /// completed or stopped.
    pub statistics: Vec<Statistic>,
    /// Observable series of replicates 0..keep_paths, at most `replicates`.
    pub keep_paths: u64,
    /// Final instances of replicates 0..keep_instances, at most
    /// `replicates`.
    pub keep_instances: u64,
    /// `crn-craft.run/1.0` only (earlier versions refuse the field with
    /// `CRN-R020` at `/schema`): probabilities of quantile statistics, at
    /// most 16, strictly increasing, each in (0, 1); absent or empty means
    /// none. The result's `statistics.quantiles` then holds, for every
    /// report time and observable, the linear-interpolation quantile of the
    /// values of all replicates, NumPy's default (Hyndman and Fan type 7):
    /// with the n values sorted ascending, `h = (n - 1) * p`, `j = floor(h)`
    /// and `g = h - j`, the quantile is `x[j] + g * (x[j+1] - x[j])`,
    /// evaluated as `numpy.quantile(x, p)` evaluates it
    /// (`x[j+1] - (1 - g) * (x[j+1] - x[j])` for `g >= 0.5`, `x[n-1]` for
    /// `h >= n - 1`, NaN if a value is NaN), so the values equal NumPy's bit
    /// for bit. Quantiles are withheld with the other statistics. They keep
    /// replicates x report times x observables values in memory during the
    /// run, counted against the 2 GiB workspace limit (`CRN-R020` at
    /// `/outputs/quantiles` above it).
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "Vec<f64>", length(max = 16))]
    pub quantiles: Option<Vec<f64>>,
}

impl Outputs {
    /// The requested quantile probabilities; an absent list and an empty
    /// list are equal.
    #[must_use]
    pub fn quantiles(&self) -> &[f64] {
        self.quantiles.as_deref().unwrap_or(&[])
    }
}

/// Explicit work, memory and output limits of a run request.
///
/// A replicate that reaches a per-replicate limit ends with the matching
/// termination and diagnostic; the other replicates continue.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    /// Events one replicate may execute; reaching it ends the replicate with
    /// `event_limit` (`CRN-R040`).
    pub max_events_per_replicate: u64,
    /// Random 64-bit words one replicate may draw, counting every word,
    /// including product draws; exhausting it ends the replicate with
    /// `sampling_limit` (`CRN-R041`). The leap methods also count their
    /// iterations against it, and `hybrid` its integration steps.
    pub max_sampling_steps_per_replicate: u64,
    /// Distinct instances alive at once, over all types; positive. Exceeding
    /// it ends the replicate with `instance_limit` (`CRN-R042`). On `tree`,
    /// so does a count above 2^53; the static methods fail the replicate
    /// with `CRN-R030` instead.
    pub max_instances: u64,
    /// Wall-clock limit for the whole run in seconds, positive and finite
    /// when present. Reaching it ends running replicates with
    /// `wall_time_limit` (`CRN-R043`) and starts no others.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_s: Option<f64>,
    /// Worker threads, 1 to 64. Completed results do not depend on it.
    #[schemars(range(min = 1, max = 64))]
    pub threads: u32,
    /// Size limit of the result as compact JSON plus a trailing newline,
    /// capped at 1 GiB. Must hold at least 160 bytes per replicate
    /// (`CRN-R020`). A larger result drops statistics, paths, observations
    /// and final instances and ends as `output_limit` (`CRN-R045`); if the
    /// records alone still do not fit, the run fails with `CRN-R045`.
    #[schemars(range(min = 1))]
    pub max_output_bytes: u64,
}

/// How Static Direct draws waiting times (rules-run/0.4 and later).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RulesWaitingTime {
    // Variant docs use `concat!` to stay out of the schema's string enum.
    #[doc = concat!(
        "Inversion of one open uniform, `-ln(1 - u) / total`; the default ",
        "before rules-run/0.6."
    )]
    Inverse,
    #[doc = concat!(
        "The Ziggurat exponential sampler. It consumes the stream ",
        "differently, so results differ from `inverse` and the sampler ",
        "identity says `ziggurat`. The default from rules-run/0.6."
    )]
    Ziggurat,
}

/// The simulation method of a run request, tagged by `kind`.
///
/// Absent means `tree`. Each method reports its own sampler identity in the
/// result's `sampler`; the same seed gives different trajectories under
/// different samplers. Every method except `tree` needs a model that lowers
/// to fixed instance tuples (`CRN-R020` at `/method` otherwise).
///
#[doc = concat!(
    "See [`RulesResult::sampler`] and the [crate documentation](crate#methods) ",
    "for which models each method accepts."
)]
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    // Empty struct variants, so `deny_unknown_fields` applies to every kind:
    // serde ignores extra keys next to the tag of a unit variant.
    /// Exact direct method over dynamic instances: rules are selected by
    /// their propensities, then instances from sum trees. Supports random
    /// product attributes and `new_id()`. Sampler
    /// `chacha8-u53/rules-direct/2`.
    Tree {},
    /// Exact direct method over the lowered network of fixed instance tuples
    /// and channels. Every positive channel must be at least 2^-53 of the
    /// total (`CRN-R031`). Sampler `chacha8-u53/rules-static-direct/1`.
    /// Ziggurat waits change the prefix to `chacha8-u53-ziggurat`, and
    /// effective sparse trackers change the suffix to
    /// `rules-static-direct-sparse-trackers/1`. From rules-run/0.6 both are
    /// the default: `waiting_time: "inverse"` with `sparse_trackers: false`
    /// runs the default of rules-run/0.5.
    StaticDirect {
        /// rules-run/0.4 and later. Absent selects inverse waits before
        /// rules-run/0.6 and Ziggurat waits from 0.6. `null` is rejected.
        #[serde(
            default,
            deserialize_with = "present",
            skip_serializing_if = "Option::is_none"
        )]
        #[schemars(with = "RulesWaitingTime")]
        waiting_time: Option<RulesWaitingTime>,
        /// rules-run/0.4 and later. `true` integrates unweighted
        /// single-class count trackers without first-passage predicates over
        /// their constant spans, which can change the last bits of an
        /// integral. Without such a tracker it changes nothing and warns
        /// with `CRN-R051`. `false` integrates every tracker at every event.
        /// Absent means `false` before rules-run/0.6; from 0.6 it integrates
        /// the eligible trackers sparsely and, without one, runs the
        /// per-event path without a warning. `null` is rejected.
        #[serde(
            default,
            deserialize_with = "present",
            skip_serializing_if = "Option::is_none"
        )]
        #[schemars(with = "bool")]
        sparse_trackers: Option<bool>,
    },
    /// rules-run/0.3 and later. Exact grouped direct method for the fixed
    /// three-class mtDNA vacancy model, sampler
    /// `chacha8-u53-ziggurat/rules-mtdna-direct/1`. A model or parameter set
    /// it cannot certify runs Static Direct instead, with that sampler
    /// identity and a `CRN-R051` warning naming the reason.
    MtdnaDirect {},
    /// Exact next reaction method with indexed clocks over the lowered
    /// channels. Sampler `chacha8-u53/rules-next-reaction/1`.
    NextReaction {},
    /// Exact thinning with the analytic ATP flow of a rules-model/0.2 model;
    /// required for ATP models and refused for others. Sampler
    /// `chacha8-u53/rules-atp-thinning/1`.
    Thinning {},
    /// Approximate adaptive Poisson leaps with exact fallback steps. Sampler
    /// `chacha8-u53/rules-approximate-cgp-poisson/2`.
    Tau {
        /// Leap tolerance, in (0, 0.1] (`CRN-R020` otherwise): the step
        /// selection bounds the mean and spread of each species and rate
        /// change in one leap by this fraction of its current value.
        epsilon: f64,
    },
    /// Approximate death-then-birth splitting on the fixed lattice
    /// `k * step`, for single-parent unit births and deaths only. Sampler
    /// `chacha8-u53/rules-approximate-binomial-split/2`.
    Binomial {
        /// Lattice step in simulated time, finite and positive (`CRN-R020`
        /// otherwise).
        step: f64,
    },
    /// rules-run/0.5 and later. Piecewise-deterministic simulation of a
    /// closed model with continuous variables (rules-model/0.3): between
    /// events the variables and the integrated total propensity follow an
    /// adaptive Dormand–Prince 5(4) integration, and the next event happens
    /// where the integrated propensity reaches a unit exponential threshold,
    /// located on the dense output. Approximate at the given tolerances;
    /// sampler `chacha8-u53/rules-approximate-hybrid-dopri54/1`. Also runs
    /// models without continuous variables. Refuses ATP models (use
    /// `thinning`) and models that do not lower to fixed tuples.
    Hybrid {
        /// Relative tolerance of the error control, in [1e-12, 1e-2]
        /// (`CRN-R020` otherwise).
        rtol: f64,
        /// Absolute tolerance of the error control, in (0, 1e-2]
        /// (`CRN-R020` otherwise). It applies to every variable and to the
        /// integrated propensity, which runs from 0 to a unit exponential.
        atol: f64,
    },
}

/// A run request: `crn-craft.rules-run/0.1` to `0.6`.
///
/// It binds one model by hash and states the time grid, the ensemble, the
/// outputs and an explicit budget. Unknown fields and duplicate keys are
/// rejected.
///
#[doc = concat!(
    "Parse and validate it together with its model bytes with [`prepare`] ",
    "or [`solve`]; [`RunRequest::trackers`] reads the tracker list.\n\n",
    "The same type reads `crn-craft.run/1.0`, the unified request contract, ",
    "which means what `rules-run/0.6` means, adds quantile statistics ",
    "([`Outputs::quantiles`]), also accepts ",
    "`crn-craft.model/1.0` models and returns `crn-craft.result/1.0`. A ",
    "`rules-run` request refuses a `model/1.0` model (`CRN-R020` at `/schema`)."
)]
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunRequest {
    // Docs of `$ref` fields use `concat!` so the schema keeps a bare `$ref`.
    #[doc = concat!("The request version, one of [`RunVersion`].")]
    pub schema: RunVersion,
    #[doc = concat!(
        "The model by portable relative path below the request's directory ",
        "and by the SHA-256 of its exact bytes. [`prepare`] checks the hash ",
        "and the path syntax; [`solve_files`] also resolves the path."
    )]
    pub model: ModelRef,
    #[doc = concat!("Horizon and report grid.")]
    pub time: Time,
    #[doc = concat!(
        "Replicates (1 to 1,000,000) and the root seed of their streams; ",
        "see [`stream`]."
    )]
    pub ensemble: Ensemble,
    #[doc = concat!("Statistics, kept paths and kept final instances.")]
    pub outputs: Outputs,
    #[doc = concat!("Work, instance, thread and output limits.")]
    pub budget: Budget,
    /// Event trackers, at most 256, each summarizing one observable per
    /// replicate in the result's `observations`. rules-run/0.2 and later.
    /// Not allowed in 0.1, even when empty.
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "Vec<TrackerSpec>")]
    pub trackers: Option<Vec<TrackerSpec>>,
    /// The simulation method. rules-run/0.2 and later; absent selects
    /// `tree`. `null` is rejected.
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "Method")]
    pub method: Option<Method>,
    /// Retains selected events of each replicate in the result's
    /// `observations`. rules-run/0.2 and later. `null` is rejected.
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "LogSpec")]
    pub event_log: Option<LogSpec>,
}

impl RunRequest {
    /// The requested trackers; an absent list and an empty list are equal.
    #[must_use]
    pub fn trackers(&self) -> &[TrackerSpec] {
        self.trackers.as_deref().unwrap_or(&[])
    }
}

/// A present optional field must hold a value. `null` is rejected rather than
/// read as absent, so a field's presence can be checked against the version.
pub(crate) fn present<'de, D, T>(d: D) -> Result<Option<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    T::deserialize(d).map(Some)
}

// ---------------------------------------------------------------------------
// Result

/// The version of a result, its `schema` field: rules-result/0.1 answers a
/// rules-run/0.1 request, rules-result/0.2 answers rules-run/0.2 to 0.4 and
/// rules-result/0.3 answers rules-run/0.5 and 0.6.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ResultVersion {
    // Variant docs use `concat!` to stay out of the schema's string enum.
    #[doc = concat!("`crn-craft.rules-result/0.1`, for rules-run/0.1 requests.")]
    #[serde(rename = "crn-craft.rules-result/0.1")]
    V01,
    #[doc = concat!(
        "`crn-craft.rules-result/0.2`, for rules-run/0.2 to 0.4 requests; ",
        "it can carry `observations`."
    )]
    #[serde(rename = "crn-craft.rules-result/0.2")]
    V02,
    #[doc = concat!(
        "`crn-craft.rules-result/0.3`, for rules-run/0.5 and 0.6 requests: ",
        "0.2 plus `method` and `continuous`."
    )]
    #[serde(rename = "crn-craft.rules-result/0.3")]
    V03,
    #[doc = concat!(
        "`crn-craft.result/1.0`, for `crn-craft.run/1.0` requests, with the ",
        "fields of rules-result/0.3 plus `statistics.quantiles`."
    )]
    #[serde(rename = "crn-craft.result/1.0")]
    V10,
}

/// How a replicate, or the whole ensemble, ended.
///
/// Values are ordered by severity, from `completed` to `output_limit`. The
/// ensemble reports `completed` when every replicate completed or stopped,
/// and otherwise the most severe replicate termination.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    // Variant docs other than `Stopped` use `concat!` to stay out of the
    // schema's string enum. `Stopped` keeps its `///` comment, which already
    // gives it a `oneOf` entry of its own.
    #[doc = concat!(
        "The replicate ran to the horizon `time.end`, possibly holding its ",
        "state after the last channel switched off."
    )]
    Completed,
    /// The model's stop condition held, ATP was depleted under thinning, a
    /// continuous variable reached its lower bound with `on_bound: "stop"`
    /// under hybrid, or an observer stopped the replicate; the state is held
    /// to the end of the grid.
    Stopped,
    #[doc = concat!(
        "`max_events_per_replicate` was reached (`CRN-R040`)."
    )]
    EventLimit,
    #[doc = concat!(
        "`max_sampling_steps_per_replicate` was reached (`CRN-R041`)."
    )]
    SamplingLimit,
    #[doc = concat!(
        "`max_instances` was exceeded, or on `tree` a count passed 2^53 ",
        "(`CRN-R042`)."
    )]
    InstanceLimit,
    #[doc = concat!(
        "`wall_time_s` elapsed (`CRN-R043`). Replicates not started by ",
        "then are recorded with `started: false`."
    )]
    WallTimeLimit,
    #[doc = concat!(
        "The run was cancelled through its [`Control`] (`CRN-R044`). ",
        "Replicates not started by then are recorded with `started: false`."
    )]
    Cancelled,
    #[doc = concat!(
        "A numerical error ended the replicate, for example a negative or ",
        "non-finite rate (`CRN-R030`) or an unresolvable weight ",
        "(`CRN-R031`); the record's `diagnostic` says which."
    )]
    Failed,
    #[doc = concat!(
        "Ensemble only: the result exceeded `max_output_bytes`, so ",
        "statistics, paths, observations and final instances were dropped ",
        "(`CRN-R045`)."
    )]
    OutputLimit,
}
impl Termination {
    /// Whether the replicate counts toward statistics: `completed` or
    /// `stopped`.
    pub fn complete(self) -> bool {
        matches!(self, Self::Completed | Self::Stopped)
    }
}

/// A species type of the model, as listed in the result's `types`.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TypeHeader {
    /// The type name from the model.
    pub name: String,
    /// The type's attribute names, in the order of an instance's
    /// `attributes` values.
    pub attributes: Vec<String>,
}

/// One instance: a type, its attribute values and its copy number.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct InstanceRecord {
    /// Index into `types`.
    pub type_index: usize,
    /// One value per attribute of the type, in declared order.
    pub attributes: Vec<f64>,
    /// Copy number, at least 1 for a final instance. Event log participants
    /// always have count 1.
    pub count: u64,
}

/// The accounting of one replicate. Every replicate has a record, including
/// replicates that a cancellation or the wall-time limit kept from starting.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Record {
    /// Replicate index, 0 to `requested_replicates - 1`; it selects the
    /// random stream.
    pub id: u64,
    /// False for a replicate that never ran because the run was cancelled or
    /// hit its wall-time limit first; its counters are then zero.
    pub started: bool,
    // Docs of `$ref` fields use `concat!` so the schema keeps a bare `$ref`.
    #[doc = concat!(
        "How this replicate ended; `stopped` for the model's stop condition ",
        "or an observer's [`observe::ObserverAction::Stop`]."
    )]
    pub termination: Termination,
    /// Events executed. Approximate methods count every firing in a leap.
    pub events: u64,
    /// End of the horizon, the time the stop condition fired, or the time a
    /// limit or failure ended the path.
    pub stop_time: f64,
    /// Report times filled with observable values: all of them for a
    /// completed or stopped replicate, the ones before the halt otherwise.
    pub saved_points: usize,
    /// Events per rule, in rule order.
    pub rule_events: Vec<u64>,
    /// Random words, including product draws; the final event's product
    /// draws can take it past the sampling budget that stopped the path.
    pub rng_words: u64,
    /// Largest number of distinct instances alive at once. On `tree` it
    /// includes the instance whose creation exceeded `max_instances`; the
    /// static methods reject that event before counting it.
    pub peak_instances: u64,
    /// Instances created by events (excluding the initial state).
    pub created_instances: u64,
    /// The limit or numerical error that ended the replicate; `null` for a
    /// completed or stopped replicate.
    pub diagnostic: Option<Diagnostic>,
}

/// The observable series of one kept replicate.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PathSeries {
    /// Replicate index, below `keep_paths`.
    pub id: u64,
    /// report time × observable: one row per saved report time, with the
    /// columns in the order of `observables`.
    pub values: Vec<Vec<f64>>,
}

/// The instances present at the end of one kept replicate.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FinalInstances {
    /// Replicate index, below `keep_instances`.
    pub id: u64,
    /// Instances with a positive count, sorted by type and attributes.
    pub instances: Vec<InstanceRecord>,
}

/// Ensemble statistics per report time and observable, accumulated in
/// replicate order so that they do not depend on the thread count.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Statistics {
    /// Replicates included: all of them, since statistics are withheld
    /// unless every replicate completed or stopped.
    pub samples: u64,
    /// report time × observable.
    pub mean: Option<Vec<Vec<f64>>>,
    /// Sample variance (divisor `samples - 1`), report time × observable;
    /// present when `variance` was requested and `samples` is at least 2.
    pub sample_variance: Option<Vec<Vec<f64>>>,
    /// Standard error of the mean, `sqrt(sample_variance / samples)`, under
    /// the same conditions as `sample_variance`.
    pub mean_standard_error: Option<Vec<Vec<f64>>>,
    /// `crn-craft.result/1.0`: one entry per probability of the request's
    /// `outputs.quantiles`, in request order, each with the
    /// linear-interpolation quantile (NumPy's default, type 7; see the
    /// request field) of every report time and observable, indexed like
    /// `mean`. Present only when the request asked for quantiles; withheld
    /// with the other statistics.
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "Vec<Quantile>")]
    pub quantiles: Option<Vec<Quantile>>,
}

/// One quantile statistic of a `crn-craft.result/1.0` result.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Quantile {
    /// The probability, as requested.
    pub p: f64,
    /// report time × observable: the linear-interpolation quantile of the
    /// replicates' values (NumPy's default, type 7).
    pub values: Vec<Vec<f64>>,
}

/// The result of a rules run: `crn-craft.rules-result/0.1`, `0.2` or `0.3`.
///
/// A result is returned even when replicates hit limits or fail: check
/// `termination`, the per-replicate records and `diagnostics`. Matrices are
/// indexed by report time, then observable.
///
#[doc = concat!(
    "Returned by [`execute`], [`solve`] and the generated-code entry points. ",
    "Serialize it with `serde_json` to obtain the JSON document. A ",
    "`crn-craft.run/1.0` request returns `crn-craft.result/1.0`, with the ",
    "fields of rules-result/0.3 plus [`Statistics::quantiles`]."
)]
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RulesResult {
    // Docs of `$ref` fields use `concat!` so the schema keeps a bare `$ref`.
    #[doc = concat!("The result version, see [`ResultVersion`].")]
    pub schema: ResultVersion,
    /// SHA-256 of the exact model bytes, lowercase hexadecimal.
    pub model_sha256: String,
    /// SHA-256 of the exact request bytes, lowercase hexadecimal. A
    /// parameter-bound candidate hashes its re-serialized request.
    pub run_sha256: String,
    /// The sampler identity of the method that ran, for example
    /// `chacha8-u53/rules-direct/2`. Replays compare equal only under the
    /// same identity.
    pub sampler: String,
    /// The request's root seed.
    pub root_seed: u64,
    /// The model's types, in declaration order.
    pub types: Vec<TypeHeader>,
    /// Rule names, in declaration order; `rule_events` and rule indices
    /// follow it.
    pub rules: Vec<String>,
    /// Observable names in column order: the model's observables sorted by
    /// name, then `ATP` for an ATP model or the continuous variables in
    /// declaration order.
    pub observables: Vec<String>,
    /// The request's report times, one per matrix row.
    pub report_times: Vec<f64>,
    /// `completed` when every replicate completed or stopped; otherwise the
    /// most severe replicate termination.
    pub termination: Termination,
    /// The request's `ensemble.replicates`.
    pub requested_replicates: u64,
    /// Replicates that completed or stopped.
    pub completed_replicates: u64,
    /// Replicates that stopped, a subset of `completed_replicates`.
    pub stopped_replicates: u64,
    /// Events over all replicates.
    pub events: u64,
    /// One record per replicate, in replicate order.
    pub replicates: Vec<Record>,
    /// Observable series of the started replicates below `keep_paths`.
    pub paths: Vec<PathSeries>,
    /// Final instances of the started replicates below `keep_instances`.
    pub final_instances: Vec<FinalInstances>,
    /// Withheld unless every replicate completed or stopped.
    pub statistics: Option<Statistics>,
    /// Ensemble-level diagnostics: `CRN-R050` for an incomplete ensemble,
    /// `CRN-R051` for a specialization that did not apply and `CRN-R045`
    /// for a trimmed result. Replicate errors are in their records.
    pub diagnostics: Vec<Diagnostic>,
    /// rules-result/0.2: tracker results and event logs, one entry per
    /// observed replicate in replicate order. Omitted when empty.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub observations: Vec<Observation>,
    #[doc = concat!(
        "rules-result/0.3: the request's method, including the tolerances ",
        "of `hybrid`, so the result records the approximation it ran at. ",
        "Omitted before 0.3."
    )]
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "Method")]
    pub method: Option<Method>,
    /// rules-result/0.3: the names of the continuous observable columns, in
    /// column order (continuous variables, or `ATP` of an ATP model). Their
    /// values are reported from the integrated or analytic flow, not only
    /// at events. Omitted when empty and before 0.3.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub continuous: Vec<String>,
}

impl Default for Method {
    fn default() -> Self {
        Self::Tree {}
    }
}

// ---------------------------------------------------------------------------

/// A validated run, ready to execute.
///
/// [`prepare`] parses both documents, validates and compiles the model,
/// checks the request against it, reserves memory and output, lowers the
/// model for the static methods and certifies `mtdna_direct` once. A
/// `Prepared` value is immutable: run it with [`execute`],
/// [`execute_with_observer`], [`execute_with_rules_kernel`] or
/// [`execute_with_rules_loop`] as often as needed. Each run starts from
/// fresh state and gives the same completed results. Use
/// [`Prepared::bind_parameters`] for other parameter values.
pub struct Prepared {
    source: Model,
    mtdna: Option<mtdna::Kernel>,
    /// Why `mtdna_direct` runs Static Direct instead, when it does.
    mtdna_fallback: Option<&'static str>,
    /// `sparse_trackers: true` was requested but no tracker is eligible.
    sparse_unused: bool,
    /// The request's method with the Static Direct defaults of its version.
    method: Method,
    compiled: model::Compiled,
    observations: observe::Spec,
    flat: Option<flat::Network>,
    /// `hybrid`: the shared solution of the independent continuous
    /// variables (rules-model/0.4), solved once here.
    independent: Option<engine::Shared>,
    request: RunRequest,
    model_sha256: String,
    run_sha256: String,
}

impl Prepared {
    /// Build an independently reusable candidate with the same run settings.
    /// Rules parameters can affect attributes, projections and even topology,
    /// so each candidate is revalidated and lowered once before simulation.
    ///
    /// The candidate is a full [`prepare`] of the re-serialized model, with
    /// `overrides` applied, and of the request bound to it. Its
    /// `model_sha256` and `run_sha256` hash those bytes, not the original
    /// files, so they can differ from the template's hashes even for an
    /// empty map. Specializations (`mtdna_direct`, `sparse_trackers`) are
    /// decided again for the new values, and a generated kernel or loop for
    /// the template is refused for the candidate.
    ///
    /// # Errors
    ///
    /// `CRN-R020` at `/parameters` for a name the model does not declare or a
    /// non-finite value; otherwise any diagnostic [`prepare`] reports for the
    /// new values.
    ///
    /// ```
    /// # use crn_craft::rules::{Control, execute, prepare};
    /// # use serde_json::json;
    /// # use std::collections::BTreeMap;
    /// # let model = serde_json::to_vec(&json!({
    /// #     "schema": "crn-craft.rules-model/0.1", "name": "immigration_death",
    /// #     "parameters": {"k_in": 10.0, "k_out": 1.0},
    /// #     "types": [{"name": "X"}], "initial": [],
    /// #     "observables": {"X": {"type": "X"}},
    /// #     "rules": [
    /// #         {"name": "immigration", "reactants": [], "rate": "k_in", "products": [{"type": "X"}]},
    /// #         {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "k_out", "products": []}
    /// #     ]
    /// # })).unwrap();
    /// # let request = serde_json::to_vec(&json!({
    /// #     "schema": "crn-craft.rules-run/0.1",
    /// #     "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(&model)},
    /// #     "time": {"end": 5.0, "report_times": [0.0, 5.0]},
    /// #     "ensemble": {"replicates": 200, "root_seed": 1},
    /// #     "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
    /// #     "budget": {"max_events_per_replicate": 100000, "max_sampling_steps_per_replicate": 1000000,
    /// #                "max_instances": 1, "threads": 1, "max_output_bytes": 1048576}
    /// # })).unwrap();
    /// let template = prepare(&model, &request).unwrap();
    /// let control = Control::default();
    /// let mean_at_end = |k_in: f64| {
    ///     let candidate = template
    ///         .bind_parameters(&BTreeMap::from([("k_in".to_string(), k_in)]))
    ///         .unwrap();
    ///     let result = execute(&candidate, &control).unwrap();
    ///     result.statistics.unwrap().mean.unwrap()[1][0]
    /// };
    /// // The stationary mean is k_in / k_out.
    /// assert!((mean_at_end(20.0) - 19.9).abs() < 2.0);
    /// assert!((mean_at_end(5.0) - 5.0).abs() < 1.0);
    ///
    /// let typo = BTreeMap::from([("k_inn".to_string(), 1.0)]);
    /// let diagnostics = template.bind_parameters(&typo).err().unwrap();
    /// assert_eq!(diagnostics[0].code, "CRN-R020");
    /// ```
    pub fn bind_parameters(
        &self,
        overrides: &BTreeMap<String, f64>,
    ) -> Result<Self, Vec<Diagnostic>> {
        let mut model = self.source.clone();
        for (name, value) in overrides {
            let Some(target) = model.parameters.get_mut(name).filter(|_| value.is_finite()) else {
                return Err(vec![request_error(
                    "/parameters",
                    format!("unknown parameter or nonfinite value: {name}"),
                    "Use declared parameter names and finite values.",
                )]);
            };
            *target = *value;
        }
        let model_bytes = serde_json::to_vec(&model).expect("validated serializable model");
        let mut request = self.request.clone();
        request.model.sha256 = crate::sha256(&model_bytes);
        let request_bytes = serde_json::to_vec(&request).expect("validated serializable request");
        prepare(&model_bytes, &request_bytes)
    }
}

fn request_error(path: &str, message: impl Into<String>, hint: &str) -> Diagnostic {
    Diagnostic::error("CRN-R020", path, message, hint)
}

fn validate_request(request: &RunRequest, model_sha256: &str) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    if !crate::validate::valid_relative_path(&request.model.path) {
        out.push(request_error("/model/path", "non-portable or non-relative model path", "Use a slash-separated path below the request directory without '.', '..', empty components, drive prefixes or backslashes."));
    }
    if request.model.sha256 != model_sha256 {
        out.push(request_error(
            "/model/sha256",
            "model SHA-256 does not match the supplied model",
            "Recompute the hash of the exact model bytes.",
        ));
    }
    if request.schema == RunVersion::V01
        && (request.trackers.is_some() || request.event_log.is_some() || request.method.is_some())
    {
        out.push(request_error(
            "/schema",
            "trackers, method and event_log require rules-run/0.2",
            "Select the version 0.2 request contract or remove these fields.",
        ));
    }
    if request.method == Some(Method::MtdnaDirect {}) && request.schema < RunVersion::V03 {
        out.push(request_error(
            "/schema",
            "mtdna_direct requires rules-run/0.3 or later",
            "Select the version 0.3 or 0.4 request contract.",
        ));
    }
    if matches!(request.method, Some(Method::Hybrid { .. })) && request.schema < RunVersion::V05 {
        out.push(request_error(
            "/schema",
            "hybrid requires rules-run/0.5",
            "Select the version 0.5 request contract.",
        ));
    }
    if matches!(
        request.method,
        Some(
            Method::StaticDirect {
                waiting_time: Some(_),
                ..
            } | Method::StaticDirect {
                sparse_trackers: Some(_),
                ..
            }
        )
    ) && request.schema < RunVersion::V04
    {
        out.push(request_error(
            "/schema",
            "Static Direct options require rules-run/0.4",
            "Select the version 0.4 request contract.",
        ));
    }
    if request.trackers().len() > model::MAX_OBSERVABLES {
        out.push(request_error(
            "/trackers",
            "too many trackers",
            "Use at most 256 trackers.",
        ));
    }
    let time = &request.time;
    if !(time.end.is_finite() && time.end > 0.0) {
        out.push(request_error(
            "/time/end",
            "end must be positive and finite",
            "Use a positive finite horizon.",
        ));
    }
    if time.report_times.is_empty() || time.report_times.len() > MAX_REPORT_TIMES {
        out.push(request_error(
            "/time/report_times",
            "between 1 and 100000 report times are required",
            "Give a sorted report grid.",
        ));
    }
    for (i, &t) in time.report_times.iter().enumerate() {
        let ordered = i == 0 || t > time.report_times[i - 1];
        if !(t.is_finite() && t >= 0.0 && t <= time.end && ordered) {
            out.push(request_error(
                &format!("/time/report_times/{i}"),
                "report times must be strictly increasing within [0, end]",
                "Sort the grid and keep it inside the horizon.",
            ));
            break;
        }
    }
    let replicates = request.ensemble.replicates;
    if replicates == 0 || replicates > MAX_REPLICATES {
        out.push(request_error(
            "/ensemble/replicates",
            "replicates must be between 1 and 1000000",
            "Use a supported ensemble size.",
        ));
    }
    if request.outputs.keep_paths > replicates || request.outputs.keep_instances > replicates {
        out.push(request_error(
            "/outputs",
            "cannot keep more replicates than requested",
            "Keep at most `replicates` paths and instance sets.",
        ));
    }
    if let Some(quantiles) = &request.outputs.quantiles {
        if request.schema < RunVersion::V10 {
            out.push(request_error(
                "/schema",
                "outputs.quantiles requires crn-craft.run/1.0",
                "Select the crn-craft.run/1.0 request contract or remove outputs.quantiles.",
            ));
        } else if quantiles.len() > MAX_QUANTILES {
            out.push(request_error(
                "/outputs/quantiles",
                "at most 16 quantiles may be requested",
                "Request at most 16 probabilities.",
            ));
        } else if let Some(i) = quantiles
            .iter()
            .enumerate()
            .position(|(i, &p)| !(p > 0.0 && p < 1.0 && (i == 0 || p > quantiles[i - 1])))
        {
            out.push(request_error(
                &format!("/outputs/quantiles/{i}"),
                "quantile probabilities must be strictly increasing within (0, 1)",
                "Give increasing probabilities such as [0.025, 0.5, 0.975].",
            ));
        }
    }
    let mut seen = std::collections::BTreeSet::new();
    if !request
        .outputs
        .statistics
        .iter()
        .all(|s| seen.insert(s.clone()))
    {
        out.push(request_error(
            "/outputs/statistics",
            "duplicate statistic",
            "List each statistic once.",
        ));
    }
    let budget = &request.budget;
    if budget.threads == 0 || budget.threads > MAX_WORKERS {
        out.push(request_error(
            "/budget/threads",
            "threads must be between 1 and 64",
            "Use a supported worker count.",
        ));
    }
    if budget.max_output_bytes == 0 {
        out.push(request_error(
            "/budget/max_output_bytes",
            "output budget must be positive",
            "Give an explicit output budget.",
        ));
    }
    if budget.max_instances == 0 {
        out.push(request_error(
            "/budget/max_instances",
            "instance budget must be positive",
            "Give an explicit instance budget.",
        ));
    }
    if budget
        .wall_time_s
        .is_some_and(|w| !(w.is_finite() && w > 0.0))
    {
        out.push(request_error(
            "/budget/wall_time_s",
            "wall time must be positive and finite",
            "Omit it or give a positive number of seconds.",
        ));
    }
    out
}

/// The version of a model check report, its `schema` field.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum CheckVersion {
    // Variant docs use `concat!` to stay out of the schema's string enum.
    #[doc = concat!("`crn-craft.rules-check/0.1`.")]
    #[serde(rename = "crn-craft.rules-check/0.1")]
    V01,
    #[doc = concat!(
        "`crn-craft.check/1.0`, the report of a `crn-craft.model/1.0` document, ",
        "with the fields of rules-check/0.1."
    )]
    #[serde(rename = "crn-craft.check/1.0")]
    V10,
}

/// The report of a model check: `crn-craft.rules-check/0.1`.
///
/// It covers rules-model/0.1 to 0.4 documents. The check stops at the first
/// stage that fails (parsing, size limits, names and types, then the rest of
/// the model), so `diagnostics` lists the errors of that stage.
///
#[doc = concat!(
    "Returned by [`check_model`] and [`check_model_file`], and printed by ",
    "`crn-craft rules-check`. A `crn-craft.model/1.0` document gets a ",
    "`crn-craft.check/1.0` report with the same fields."
)]
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RulesCheck {
    // Docs of `$ref` fields use `concat!` so the schema keeps a bare `$ref`.
    #[doc = concat!("The report version, see [`CheckVersion`].")]
    pub schema: CheckVersion,
    /// True when the model parses and compiles without diagnostics.
    pub valid: bool,
    /// SHA-256 of the model bytes, when they could be read.
    pub model_sha256: Option<String>,
    /// Errors with stable codes and JSON Pointer paths into the model:
    /// `CRN-E001` to `CRN-E004` from parsing, `CRN-R010` to `CRN-R015` from
    /// validation, `CRN-E030` for an unreadable file.
    pub diagnostics: Vec<Diagnostic>,
}

/// Parses and checks a model; every expression is compiled and name-checked.
///
/// This validates the model alone. A run request adds its own checks in
/// [`prepare`], for example whether a static method can lower the model.
///
/// ```
/// use crn_craft::rules::check_model;
///
/// let model = br#"{
///     "schema": "crn-craft.rules-model/0.1",
///     "name": "decay",
///     "parameters": {"k": 0.5},
///     "types": [{"name": "A"}],
///     "initial": [{"type": "A", "count": 100}],
///     "observables": {"A": {"type": "A"}},
///     "rules": [{"name": "decay", "reactants": [{"label": "a", "type": "A"}],
///                "rate": "k", "products": []}]
/// }"#;
/// let report = check_model(model);
/// assert!(report.valid);
/// assert_eq!(report.model_sha256.as_deref(), Some(crn_craft::sha256(model).as_str()));
///
/// // A product of an undeclared type is reported with a JSON Pointer.
/// let broken = String::from_utf8_lossy(model).replace(r#""products": []"#,
///                                                      r#""products": [{"type": "B"}]"#);
/// let report = check_model(broken.as_bytes());
/// assert!(!report.valid);
/// assert_eq!(report.diagnostics[0].code, "CRN-R011");
/// assert_eq!(report.diagnostics[0].path, "/rules/0/products/0/type");
/// ```
pub fn check_model(model_bytes: &[u8]) -> RulesCheck {
    let parsed = parse_model(model_bytes);
    // The report version follows the model's: also for a model/1.0
    // document that does not parse.
    let unified = match &parsed {
        Ok(model) => model.schema == model::ModelVersion::V10,
        Err(_) => crate::schema_field(model_bytes).as_deref() == Some("crn-craft.model/1.0"),
    };
    let diagnostics = parsed
        .map_err(|e| vec![e])
        .and_then(|model| model::compile(&model).map(|_| ()))
        .err()
        .unwrap_or_default();
    RulesCheck {
        schema: if unified {
            CheckVersion::V10
        } else {
            CheckVersion::V01
        },
        valid: diagnostics.is_empty(),
        model_sha256: Some(crate::sha256(model_bytes)),
        diagnostics,
    }
}

/// Validates a model and a run request and prepares the run.
///
/// Both arguments are the exact JSON bytes; the request's `model.sha256` must
/// be the SHA-256 of `model_bytes`. Its `model.path` must be a portable
/// relative path, but is not resolved here (see [`solve_files`]).
/// Preparation checks the request against its version, compiles the model,
/// compiles trackers and the event log, admits the memory and output
/// reservations, and lowers the model for every method except `tree`.
///
/// # Errors
///
/// All diagnostics of the first failing stage. A document that does not
/// parse gives one `CRN-E001` to `CRN-E004` diagnostic; the request checks
/// and the model validation are reported together (`CRN-R020` for the
/// request, `CRN-R010` to `CRN-R015` for the model); trackers, the event
/// log, the budget reservation, the method's parameters and the lowering
/// each fail with `CRN-R020`.
///
/// ```
/// use crn_craft::rules::{Control, execute, prepare};
/// use serde_json::json;
///
/// # let model = serde_json::to_vec(&json!({
/// #     "schema": "crn-craft.rules-model/0.1", "name": "immigration_death",
/// #     "parameters": {"k_in": 10.0, "k_out": 1.0},
/// #     "types": [{"name": "X"}], "initial": [],
/// #     "observables": {"X": {"type": "X"}},
/// #     "rules": [
/// #         {"name": "immigration", "reactants": [], "rate": "k_in", "products": [{"type": "X"}]},
/// #         {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "k_out", "products": []}
/// #     ]
/// # })).unwrap();
/// // Static Direct with a tracker on X and a bounded log of death events.
/// let request = serde_json::to_vec(&json!({
///     "schema": "crn-craft.rules-run/0.2",
///     "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(&model)},
///     "time": {"end": 2.0, "report_times": [0.0, 1.0, 2.0]},
///     "ensemble": {"replicates": 4, "root_seed": 42},
///     "outputs": {"statistics": [], "keep_paths": 0, "keep_instances": 4},
///     "budget": {"max_events_per_replicate": 10000, "max_sampling_steps_per_replicate": 100000,
///                "max_instances": 1, "threads": 1, "max_output_bytes": 1048576},
///     "method": {"kind": "static_direct"},
///     "trackers": [{"observable": "X", "keep_episodes": 2}],
///     "event_log": {"max_entries": 3, "rules": ["death"]}
/// }))
/// .unwrap();
/// let prepared = prepare(&model, &request).unwrap();
/// let result = execute(&prepared, &Control::default()).unwrap();
/// assert_eq!(result.sampler, "chacha8-u53/rules-static-direct/1");
///
/// let first = &result.observations[0];
/// let tracker = &first.trackers[0];
/// assert_eq!(tracker.observable, "X");
/// assert!(tracker.maximum >= tracker.minimum);
/// assert!(first.log.len() <= 3);
/// assert!(first.log.iter().all(|entry| entry.rule == Some(1)));
///
/// // Other model bytes, even a trailing newline, no longer match the hash.
/// let mut edited = model.clone();
/// edited.push(b'\n');
/// let diagnostics = prepare(&edited, &request).err().unwrap();
/// assert_eq!((diagnostics[0].code.as_str(), diagnostics[0].path.as_str()),
///            ("CRN-R020", "/model/sha256"));
/// ```
pub fn prepare(model_bytes: &[u8], request_bytes: &[u8]) -> Result<Prepared, Vec<Diagnostic>> {
    let model: Model = parse_model(model_bytes).map_err(|e| vec![e])?;
    let request: RunRequest = crate::strict_json::parse(request_bytes).map_err(|e| {
        if crate::fsp::v1::is_fsp_request(request_bytes) {
            vec![crate::fsp::v1::not_a_simulation()]
        } else {
            vec![e]
        }
    })?;
    let model_sha256 = crate::sha256(model_bytes);
    let mut diagnostics = validate_request(&request, &model_sha256);
    if model.schema == model::ModelVersion::V10 && request.schema != RunVersion::V10 {
        // Old requests cannot reach the features of the 1.0 contracts.
        diagnostics.push(request_error(
            "/schema",
            "a crn-craft.model/1.0 model runs only with a crn-craft.run/1.0 request",
            "Select the crn-craft.run/1.0 request contract.",
        ));
    }
    let compiled = match model::compile(&model) {
        Ok(c) => Some(c),
        Err(mut d) => {
            diagnostics.append(&mut d);
            None
        }
    };
    if !diagnostics.is_empty() {
        return Err(diagnostics);
    }
    let compiled = compiled.expect("compiled when valid");
    let mut observations = observe::Spec::compile(
        &model,
        &compiled,
        request.trackers(),
        request.event_log.as_ref(),
    )
    .map_err(|(path, e)| {
        vec![request_error(
            path,
            e,
            "Use declared observable and rule names.",
        )]
    })?;
    match request.method {
        Some(Method::Tau { epsilon })
            if !epsilon.is_finite() || epsilon <= 0.0 || epsilon > 0.1 =>
        {
            return Err(vec![request_error(
                "/method/epsilon",
                "epsilon must be in (0, 0.1]",
                "Set an explicit tolerance.",
            )]);
        }
        Some(Method::Binomial { step }) if !step.is_finite() || step <= 0.0 => {
            return Err(vec![request_error(
                "/method/step",
                "step must be finite and positive",
                "Specify a fixed step.",
            )]);
        }
        Some(Method::Hybrid { rtol, .. }) if !(1e-12..=1e-2).contains(&rtol) => {
            return Err(vec![request_error(
                "/method/rtol",
                "rtol must be in [1e-12, 1e-2]",
                "Set an explicit relative tolerance, for example 1e-6.",
            )]);
        }
        Some(Method::Hybrid { atol, .. }) if !(atol > 0.0 && atol <= 1e-2) => {
            return Err(vec![request_error(
                "/method/atol",
                "atol must be in (0, 1e-2]",
                "Set an explicit absolute tolerance, for example 1e-9.",
            )]);
        }
        _ => {}
    }
    if let Some(diagnostic) = higher_order_method(&compiled, request.method.unwrap_or_default()) {
        return Err(vec![diagnostic]);
    }
    let hybrid = matches!(request.method, Some(Method::Hybrid { .. }));
    if !compiled.continuous.is_empty() && !hybrid {
        return Err(vec![request_error(
            "/method",
            "continuous variables require the hybrid method",
            "Select {\"kind\": \"hybrid\"} with rules-run/0.5.",
        )]);
    }
    if compiled.atp.is_some() != matches!(request.method, Some(Method::Thinning {})) {
        let message = if compiled.atp.is_some() {
            "ATP models require the thinning method"
        } else {
            "thinning requires an ATP rules model"
        };
        return Err(vec![request_error(
            "/method",
            message,
            "Use thinning with an ATP rules model.",
        )]);
    }
    let mut flat = match request.method.unwrap_or_default() {
        Method::Tree {} => None,
        Method::StaticDirect { .. }
        | Method::MtdnaDirect {}
        | Method::NextReaction {}
        | Method::Thinning {}
        | Method::Tau { .. }
        | Method::Binomial { .. }
        | Method::Hybrid { .. } => Some(flat::Network::compile(&compiled).map_err(|e| {
            vec![request_error(
                "/method",
                e,
                "Use a statically closed model or select tree Direct.",
            )]
        })?),
    };
    if matches!(request.method, Some(Method::Binomial { .. }))
        && flat.as_ref().is_some_and(|n| !n.supports_binomial())
    {
        return Err(vec![request_error(
            "/method",
            "binomial splitting supports single-parent unit births and deaths",
            "Select tau or exact SSA for other channels.",
        )]);
    }
    let sizes = observation_sizes(&compiled, request.method, flat.as_ref());
    let reserved = observations.reservation(&sizes);
    let reservation = reserved.per_replicate() * u128::from(request.ensemble.replicates);
    let budget = MAX_WORKSPACE_BYTES.min(request.budget.max_output_bytes);
    if reservation > u128::from(budget) {
        return Err(vec![reservation_error(
            &reserved,
            request.ensemble.replicates,
            reservation,
            request.budget.max_output_bytes,
        )]);
    }
    let diagnostics = admission(&compiled, &request, reservation);
    if !diagnostics.is_empty() {
        return Err(diagnostics);
    }
    let mut sparse_unused = false;
    if let Some(net) = &mut flat {
        let sparse = sparse_request(&request);
        if sparse != Sparse::Off {
            observations.enable_sparse(net);
            // No eligible tracker: the run is the per-event one, and an
            // explicit request says so.
            sparse_unused = sparse == Sparse::Required && !observations.sparse_effective();
        }
        net.plan_observables(&compiled, &observations);
    }
    let (mtdna, mtdna_fallback) = match (request.method, flat.as_ref()) {
        (Some(Method::MtdnaDirect {}), Some(net)) => {
            match mtdna::Kernel::compile(net, &compiled, &observations) {
                Ok(kernel) => (Some(kernel), None),
                Err(reason) => (None, Some(reason)),
            }
        }
        _ => (None, None),
    };
    // The independent subsystem does not depend on the replicate, so one
    // solution at the request's tolerances serves every replicate and
    // thread.
    let independent = match request.method {
        Some(Method::Hybrid { rtol, atol }) => engine::Shared::solve(
            &compiled,
            rtol,
            atol,
            request.time.end,
            request.budget.max_sampling_steps_per_replicate,
        ),
        _ => None,
    };
    Ok(Prepared {
        source: model,
        mtdna,
        mtdna_fallback,
        sparse_unused,
        method: effective_method(&request),
        compiled,
        flat,
        independent,
        observations,
        request,
        model_sha256,
        run_sha256: crate::sha256(request_bytes),
    })
}

/// Parses a rules model strictly. A pattern `count` before
/// `crn-craft.model/1.0` gets the strict parser's unknown-field error of
/// 0.6.0.
fn parse_model(model_bytes: &[u8]) -> Result<Model, Diagnostic> {
    let model: Model = crate::strict_json::parse(model_bytes)?;
    match model::legacy_count(&model) {
        Some(diagnostic) => Err(diagnostic),
        None => Ok(model),
    }
}

/// `CRN-R020` at `/method` for a higher-order rule on a method that cannot
/// sample it: `tree` (also the default) and `mtdna_direct`.
fn higher_order_method(compiled: &model::Compiled, method: Method) -> Option<Diagnostic> {
    let name = match method {
        Method::Tree {} => "tree",
        Method::MtdnaDirect {} => "mtdna_direct",
        _ => return None,
    };
    let rule = compiled
        .rules
        .iter()
        .find(|r| matches!(r.shape, model::Shape::Multi(_)))?;
    Some(request_error(
        "/method",
        format!(
            "rule `{}` takes more than two reactant molecules, which {name} does not sample",
            rule.name
        ),
        "Select static_direct, next_reaction, tau or hybrid for higher-order rules.",
    ))
}

/// How a request asks for sparse tracker integration on Static Direct.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Sparse {
    /// Per-event integration of every tracker.
    Off,
    /// The rules-run/0.6 default: sparse where eligible, silently.
    Auto,
    /// `sparse_trackers: true`: sparse where eligible, `CRN-R051` otherwise.
    Required,
}

fn sparse_request(request: &RunRequest) -> Sparse {
    match request.method {
        Some(Method::StaticDirect {
            sparse_trackers: Some(true),
            ..
        }) => Sparse::Required,
        Some(Method::StaticDirect {
            sparse_trackers: None,
            ..
        }) if request.schema >= RunVersion::V06 => Sparse::Auto,
        _ => Sparse::Off,
    }
}

/// The method that runs: the request's method, or `tree`, with an absent
/// Static Direct `waiting_time` filled in by the default of the request's
/// version. Every other method is the request's own.
fn effective_method(request: &RunRequest) -> Method {
    match request.method.unwrap_or_default() {
        Method::StaticDirect {
            waiting_time: None,
            sparse_trackers,
        } if request.schema >= RunVersion::V06 => Method::StaticDirect {
            waiting_time: Some(RulesWaitingTime::Ziggurat),
            sparse_trackers,
        },
        method => method,
    }
}

fn workers(request: &RunRequest) -> u64 {
    u64::from(request.budget.threads.clamp(1, MAX_WORKERS)).min(request.ensemble.replicates)
}

/// Replicates in flight or waiting for ordered reduction, per worker.
const WINDOW_PER_WORKER: u64 = 4;

/// The sizes that bound one log entry of this model and method.
fn observation_sizes(
    compiled: &model::Compiled,
    method: Option<Method>,
    flat: Option<&flat::Network>,
) -> observe::Sizes {
    observe::Sizes {
        attributes: compiled
            .types
            .iter()
            .map(|t| t.attributes.len())
            .max()
            .unwrap_or(0),
        observables: compiled.observables.len(),
        rules: compiled.rules.len(),
        participants: compiled
            .rules
            .iter()
            .map(|r| r.reactant_types.len() + r.products.len())
            .max()
            .unwrap_or(0),
        leap_channels: match method {
            Some(Method::Tau { .. } | Method::Binomial { .. }) => flat.map(|n| n.channels.len()),
            _ => None,
        },
    }
}

/// `CRN-R020` at `/budget` for an observation reservation above the output
/// budget or the workspace limit. The message itemizes the reservation; the
/// hint names the smallest `max_output_bytes` that admits it.
fn reservation_error(
    reserved: &observe::Reservation,
    replicates: u64,
    total: u128,
    max_output_bytes: u64,
) -> Diagnostic {
    let mut parts = Vec::new();
    if reserved.entries > 0 {
        parts.push(format!(
            "{} log entries x {} bytes",
            reserved.entries, reserved.entry
        ));
    }
    if reserved.trackers > 0 {
        parts.push(format!("{} tracker bytes", reserved.trackers));
    }
    let message = format!(
        "observation reservation exceeds budget: {replicates} replicates x {} bytes ({}) = {total} bytes, above max_output_bytes = {max_output_bytes}",
        reserved.per_replicate(),
        parts.join(" + "),
    );
    let hint = if total <= u128::from(MAX_WORKSPACE_BYTES) {
        format!(
            "Set max_output_bytes to at least {total}, or reduce replicates, event_log.max_entries or keep_episodes."
        )
    } else {
        format!(
            "The reservation exceeds the workspace limit of {MAX_WORKSPACE_BYTES} bytes; reduce replicates, event_log.max_entries or keep_episodes."
        )
    };
    request_error("/budget", message, &hint)
}

/// Workspace and output reservations from sizes alone, checked before any run.
fn admission(
    model: &model::Compiled,
    request: &RunRequest,
    observation_bytes: u128,
) -> Vec<Diagnostic> {
    let n = u128::from(request.ensemble.replicates);
    let series = request.time.report_times.len() as u128 * model.observables.len() as u128 * 8;
    let record = 256 + 8 * model.rules.len() as u128;
    let attributes = model
        .types
        .iter()
        .map(|t| t.attributes.len())
        .max()
        .unwrap_or(0) as u128;
    let instances = u128::from(request.outputs.keep_instances)
        * u128::from(request.budget.max_instances)
        * (8 * attributes + 48);
    let in_flight = u128::from(workers(request)) * u128::from(1 + WINDOW_PER_WORKER);
    let workspace = n * record
        + (in_flight + u128::from(request.outputs.keep_paths)) * series
        + instances
        + observation_bytes;
    // Lower bound of the compact result: a record is at least 160 bytes.
    let evidence = n * 160;
    let mut out = Vec::new();
    if workspace > u128::from(MAX_WORKSPACE_BYTES) {
        out.push(request_error(
            "/budget",
            format!("execution reservation of {workspace} bytes exceeds {MAX_WORKSPACE_BYTES}"),
            "Use fewer replicates, report times, kept paths or kept instance sets, or a smaller max_instances.",
        ));
    } else if !request.outputs.quantiles().is_empty() {
        // Quantiles keep every replicate's series, plus one column to sort.
        let quantiles = n * series + n * 8;
        if workspace + quantiles > u128::from(MAX_WORKSPACE_BYTES) {
            out.push(request_error(
                "/outputs/quantiles",
                format!(
                    "the quantile buffer of {n} replicates x {} report times x {} observables x 8 bytes, with its sort column {quantiles} bytes, brings the execution reservation to {} bytes, above {MAX_WORKSPACE_BYTES}",
                    request.time.report_times.len(),
                    model.observables.len(),
                    workspace + quantiles
                ),
                "Request fewer report times, fewer replicates or no quantiles.",
            ));
        }
    }
    if evidence > u128::from(request.budget.max_output_bytes.min(MAX_RESULT_BYTES)) {
        out.push(request_error(
            "/budget/max_output_bytes",
            format!("{n} replicate records need at least {evidence} output bytes"),
            "Increase max_output_bytes or use fewer replicates.",
        ));
    }
    out
}

fn record(id: u64, path: &engine::PathOutput) -> Record {
    Record {
        id,
        started: true,
        termination: path.termination,
        events: path.events,
        stop_time: path.stop_time,
        saved_points: path.saved_points,
        rule_events: path.rule_events.clone(),
        rng_words: path.rng_words,
        peak_instances: path.peak_instances,
        created_instances: path.created_instances,
        diagnostic: path.diagnostic.clone(),
    }
}

fn fits_json(value: &impl Serialize, limit: u64) -> bool {
    struct Counter {
        used: u64,
        limit: u64,
    }
    impl std::io::Write for Counter {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            let n = bytes.len() as u64;
            if n > self.limit.saturating_sub(self.used) {
                return Err(std::io::Error::other("output limit"));
            }
            self.used += n;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    // The budget covers the compact JSON and a trailing newline.
    let mut counter = Counter {
        used: 0,
        limit: limit.saturating_sub(1),
    };
    serde_json::to_writer(&mut counter, value).is_ok()
}

/// Ordered reduction: replicates are folded in id order (Welford), so the
/// statistics do not depend on the thread count; only kept series stay, and
/// with quantiles every complete replicate's series at its index.
struct Reducer {
    next: u64,
    pending: BTreeMap<u64, engine::PathOutput>,
    records: Vec<Record>,
    observations: Vec<Observation>,
    paths: Vec<PathSeries>,
    finals: Vec<FinalInstances>,
    samples: u64,
    mean: Vec<f64>,
    m2: Vec<f64>,
    /// Replicate × report time × observable, for quantiles; empty when the
    /// request asks for none.
    quantiles: Vec<f64>,
    halted: bool,
}

impl Reducer {
    #[allow(clippy::cast_precision_loss)]
    fn fold(&mut self, id: u64, mut path: engine::PathOutput, keep_paths: u64, columns: usize) {
        self.records.push(record(id, &path));
        if let Some(obs) = path.observation.take() {
            self.observations.push(obs);
        }
        if path.termination.complete() {
            self.samples += 1;
            let k = self.samples as f64;
            for ((m, s), &x) in self.mean.iter_mut().zip(&mut self.m2).zip(&path.values) {
                let delta = x - *m;
                *m += delta / k;
                *s += delta * (x - *m);
            }
            if !self.quantiles.is_empty() {
                let cells = self.mean.len();
                let start = usize::try_from(id).expect("replicate index fits usize") * cells;
                self.quantiles[start..start + cells].copy_from_slice(&path.values[..cells]);
            }
        }
        if id < keep_paths {
            self.paths.push(PathSeries {
                id,
                values: if columns == 0 {
                    vec![Vec::new(); path.saved_points]
                } else {
                    path.values[..path.saved_points * columns]
                        .chunks(columns)
                        .map(<[f64]>::to_vec)
                        .collect()
                },
            });
        }
        if let Some(instances) = path.instances {
            self.finals.push(FinalInstances { id, instances });
        }
    }
}

#[allow(clippy::cast_precision_loss)]
fn statistics(
    wanted: &[Statistic],
    probabilities: &[f64],
    reducer: &Reducer,
    rows: usize,
    columns: usize,
) -> Statistics {
    let n = reducer.samples;
    let chunk = |flat: Vec<f64>| -> Vec<Vec<f64>> {
        if columns == 0 {
            vec![Vec::new(); reducer.mean.len()]
        } else {
            flat.chunks(columns).map(<[f64]>::to_vec).collect()
        }
    };
    let variance: Vec<f64> = reducer
        .m2
        .iter()
        .map(|s| s / n.saturating_sub(1) as f64)
        .collect();
    let has_variance = wanted.contains(&Statistic::Variance) && n > 1;
    Statistics {
        samples: n,
        mean: wanted
            .contains(&Statistic::Mean)
            .then(|| chunk(reducer.mean.clone())),
        mean_standard_error: has_variance
            .then(|| chunk(variance.iter().map(|v| libm::sqrt(v / n as f64)).collect())),
        sample_variance: has_variance.then(|| chunk(variance)),
        // Statistics are reported only when every replicate completed, so
        // the buffer holds all of them.
        quantiles: (!probabilities.is_empty()).then(|| {
            quantile::statistics(
                &reducer.quantiles,
                usize::try_from(n).expect("replicate count fits usize"),
                rows,
                columns,
                probabilities,
            )
        }),
    }
}

/// Runs the ensemble. Replicates are independent keyed streams, so results do
/// not depend on the thread count.
///
/// The engine starts `min(budget.threads, replicates)` workers, each taking
/// the next replicate, and reduces finished replicates in replicate order.
/// Replicates that hit a limit or fail do not make this call fail: they are
/// reported in [`RulesResult::replicates`], the ensemble
/// [`RulesResult::termination`] names the most severe of them, and
/// statistics are withheld with a `CRN-R050` warning. A cancellation or the
/// wall-time limit stops the remaining replicates from starting.
///
/// # Errors
///
/// `CRN-R045` at `/budget/max_output_bytes` when even the replicate records,
/// without statistics, paths, observations and final instances, exceed the
/// output budget. A result that fits only after that trimming is returned
/// with termination `output_limit`.
#[allow(clippy::missing_panics_doc, clippy::too_many_lines)]
pub fn execute(prepared: &Prepared, control: &Control) -> Result<RulesResult, Diagnostic> {
    execute_with_observer(prepared, control, &NoObserver)
}

/// Observe committed events on each native worker. Ordering is per replicate.
///
/// Like [`execute`], but `observer` sees every committed event: one per
/// exact step, or one per leap of the approximate methods. It runs on the
/// worker that simulates the replicate, after the state update and the
/// trackers and before the stop condition. Events of one replicate arrive in
/// order; events of different replicates can interleave. Returning
/// [`observe::ObserverAction::Stop`] ends that replicate as `stopped` at the
/// event's time, which counts as complete for statistics. An observer that
/// always continues leaves the serialized result byte-identical to
/// [`execute`]. Observation consumes no random words.
///
/// # Errors
///
/// As for [`execute`].
///
/// # Panics
///
/// `on_event` must not panic. A panic is not caught: it ends that worker
/// without a result for its replicate, and the remaining workers can then
/// wait for that result indefinitely, so the call may never return. Catch
/// panics inside the observer and return
/// [`observe::ObserverAction::Stop`] instead.
///
/// ```
/// use crn_craft::rules::observe::{Event, EventObserver, ObserverAction};
/// use crn_craft::rules::{Control, Termination, execute_with_observer, prepare};
/// use serde_json::json;
/// use std::sync::atomic::{AtomicU64, Ordering};
///
/// /// Counts events and stops a replicate once X reaches 5 copies.
/// struct StopAtFive {
///     events: AtomicU64,
/// }
///
/// impl EventObserver for StopAtFive {
///     fn on_event(&self, event: &Event<'_>) -> ObserverAction {
///         self.events.fetch_add(1, Ordering::Relaxed);
///         if event.observables[0] >= 5.0 {
///             ObserverAction::Stop
///         } else {
///             ObserverAction::Continue
///         }
///     }
/// }
///
/// # let model = serde_json::to_vec(&json!({
/// #     "schema": "crn-craft.rules-model/0.1", "name": "immigration_death",
/// #     "parameters": {"k_in": 10.0, "k_out": 1.0},
/// #     "types": [{"name": "X"}], "initial": [],
/// #     "observables": {"X": {"type": "X"}},
/// #     "rules": [
/// #         {"name": "immigration", "reactants": [], "rate": "k_in", "products": [{"type": "X"}]},
/// #         {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "k_out", "products": []}
/// #     ]
/// # })).unwrap();
/// # let request = serde_json::to_vec(&json!({
/// #     "schema": "crn-craft.rules-run/0.1",
/// #     "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(&model)},
/// #     "time": {"end": 5.0, "report_times": [0.0, 5.0]},
/// #     "ensemble": {"replicates": 16, "root_seed": 1},
/// #     "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
/// #     "budget": {"max_events_per_replicate": 100000, "max_sampling_steps_per_replicate": 1000000,
/// #                "max_instances": 1, "threads": 2, "max_output_bytes": 1048576}
/// # })).unwrap();
/// let prepared = prepare(&model, &request).unwrap();
/// let observer = StopAtFive { events: AtomicU64::new(0) };
/// let result = execute_with_observer(&prepared, &Control::default(), &observer).unwrap();
///
/// assert_eq!(observer.events.load(Ordering::Relaxed), result.events);
/// assert_eq!(result.stopped_replicates, result.requested_replicates);
/// assert!(result.replicates.iter().all(|r| r.termination == Termination::Stopped));
/// // A stopped replicate holds its state: X is 5 at the end of the grid.
/// assert_eq!(result.statistics.unwrap().mean.unwrap()[1][0], 5.0);
/// ```
pub fn execute_with_observer<O: EventObserver>(
    prepared: &Prepared,
    control: &Control,
    observer: &O,
) -> Result<RulesResult, Diagnostic> {
    execute_inner(prepared, control, observer, None, None)
}

fn execute_inner<O: EventObserver>(
    prepared: &Prepared,
    control: &Control,
    observer: &O,
    generated: Option<&dyn RulesKernel>,
    generated_loop: Option<&dyn RulesLoop>,
) -> Result<RulesResult, Diagnostic> {
    let model = &prepared.compiled;
    let request = &prepared.request;
    let start = Instant::now();
    let deadline = request
        .budget
        .wall_time_s
        .map(|s| start + Duration::from_secs_f64(s.min(1e9)));
    let limits = engine::Limits {
        max_events: request.budget.max_events_per_replicate,
        max_words: request.budget.max_sampling_steps_per_replicate,
        max_instances: request.budget.max_instances,
        deadline,
    };
    let n = request.ensemble.replicates;
    let columns = model.observables.len();
    let rows = request.time.report_times.len();
    let keep_paths = request.outputs.keep_paths;
    let keep_instances = request.outputs.keep_instances;
    let next = AtomicU64::new(0);
    let halt = AtomicBool::new(false);
    let workers = workers(request);
    let window = workers * WINDOW_PER_WORKER;
    let reducer = Mutex::new(Reducer {
        next: 0,
        pending: BTreeMap::new(),
        records: Vec::with_capacity(usize::try_from(n).unwrap_or(0)),
        observations: Vec::new(),
        paths: Vec::new(),
        finals: Vec::new(),
        samples: 0,
        mean: vec![0.0; rows * columns],
        m2: vec![0.0; rows * columns],
        quantiles: if request.outputs.quantiles().is_empty() {
            Vec::new()
        } else {
            // Admitted by `prepare` against the workspace limit.
            vec![0.0; usize::try_from(n).expect("replicate count fits usize") * rows * columns]
        },
        halted: false,
    });
    let reduced = Condvar::new();
    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| {
                loop {
                    if halt.load(Ordering::Relaxed) {
                        break;
                    }
                    let id = next.fetch_add(1, Ordering::Relaxed);
                    if id >= n {
                        break;
                    }
                    {
                        // Do not run ahead of the ordered reduction.
                        let mut r = reducer.lock().expect("no poisoned workers");
                        while id >= r.next + window && !r.halted {
                            r = reduced.wait(r).expect("no poisoned workers");
                        }
                        if r.halted {
                            break;
                        }
                    }
                    let path = engine::Path {
                        model,
                        report_times: &request.time.report_times,
                        end: request.time.end,
                        limits: &limits,
                        control,
                        keep_instances: id < keep_instances,
                        observations: &prepared.observations,
                        flat: prepared.flat.as_ref(),
                        generated,
                        generated_loop,
                        mtdna: prepared.mtdna.as_ref(),
                        method: prepared.method,
                        independent: prepared.independent.as_ref(),
                        id,
                    };
                    let output = path.run(&mut stream(request.ensemble.root_seed, id), observer);
                    let stop = matches!(
                        output.termination,
                        Termination::Cancelled | Termination::WallTimeLimit
                    );
                    let mut r = reducer.lock().expect("no poisoned workers");
                    r.pending.insert(id, output);
                    while let Some(path) = {
                        let key = r.next;
                        r.pending.remove(&key)
                    } {
                        let id = r.next;
                        r.fold(id, path, keep_paths, columns);
                        r.next += 1;
                    }
                    if stop {
                        halt.store(true, Ordering::Relaxed);
                        r.halted = true;
                    }
                    reduced.notify_all();
                }
            });
        }
    });
    let mut reducer = reducer.into_inner().expect("no poisoned workers");
    // After a halt, finished replicates beyond a gap are reduced in id order;
    // the gap's never-started replicates get unstarted records below.
    let control_stop = if control.is_cancelled() {
        Termination::Cancelled
    } else {
        Termination::WallTimeLimit
    };
    let rest: Vec<_> = std::mem::take(&mut reducer.pending).into_iter().collect();
    let mut rest = rest.into_iter().peekable();
    while reducer.next < n {
        let id = reducer.next;
        if let Some((_, path)) = rest.next_if(|(i, _)| *i == id) {
            reducer.fold(id, path, keep_paths, columns);
        } else {
            reducer.records.push(Record {
                id,
                started: false,
                termination: control_stop,
                events: 0,
                stop_time: 0.0,
                saved_points: 0,
                rule_events: vec![0; model.rules.len()],
                rng_words: 0,
                peak_instances: 0,
                created_instances: 0,
                diagnostic: None,
            });
        }
        reducer.next += 1;
    }
    let replicates = std::mem::take(&mut reducer.records);
    let complete = replicates.iter().all(|r| r.termination.complete());
    let termination = replicates
        .iter()
        .map(|r| {
            if r.termination == Termination::Stopped {
                Termination::Completed
            } else {
                r.termination
            }
        })
        .max()
        .unwrap_or(Termination::Failed);
    // Quantiles are statistics too: withheld by the same rule.
    let statistics = (complete
        && !(request.outputs.statistics.is_empty() && request.outputs.quantiles().is_empty()))
    .then(|| {
        statistics(
            &request.outputs.statistics,
            request.outputs.quantiles(),
            &reducer,
            rows,
            columns,
        )
    });
    let paths = std::mem::take(&mut reducer.paths);
    let final_instances = std::mem::take(&mut reducer.finals);
    let mut diagnostics = Vec::new();
    if !complete {
        diagnostics.push(Diagnostic::warning(
            "CRN-R050",
            "/replicates",
            "ensemble incomplete; statistics withheld",
            "See the per-replicate terminations and diagnostics.",
        ));
    }
    if let Some(reason) = prepared.mtdna_fallback {
        diagnostics.push(Diagnostic::warning(
            "CRN-R051",
            "/method",
            format!("mtdna_direct ran Static Direct: {reason}"),
            "The result is exact; check `sampler` before comparing replays.",
        ));
    }
    if prepared.sparse_unused {
        diagnostics.push(Diagnostic::warning(
            "CRN-R051",
            "/method/sparse_trackers",
            "sparse_trackers changed nothing: no tracker is an unweighted single-class count without a first-passage predicate",
            "Per-event integration ran; the sampler identity is the default one.",
        ));
    }
    // The identity names what ran: the effective waits, and sparse
    // integration only where a tracker was eligible.
    let ziggurat = matches!(
        prepared.method,
        Method::StaticDirect {
            waiting_time: Some(RulesWaitingTime::Ziggurat),
            ..
        }
    );
    let sparse = prepared.observations.sparse_effective();
    let mut result = RulesResult {
        schema: match request.schema {
            RunVersion::V01 => ResultVersion::V01,
            RunVersion::V02 | RunVersion::V03 | RunVersion::V04 => ResultVersion::V02,
            RunVersion::V05 | RunVersion::V06 => ResultVersion::V03,
            RunVersion::V10 => ResultVersion::V10,
        },
        model_sha256: prepared.model_sha256.clone(),
        run_sha256: prepared.run_sha256.clone(),
        sampler: match prepared.method {
            Method::Tree {} => SAMPLER,
            Method::MtdnaDirect {} if prepared.mtdna.is_some() => {
                "chacha8-u53-ziggurat/rules-mtdna-direct/1"
            }
            Method::StaticDirect { .. } => match (ziggurat, sparse) {
                (true, true) => "chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1",
                (true, false) => "chacha8-u53-ziggurat/rules-static-direct/1",
                (false, true) => "chacha8-u53/rules-static-direct-sparse-trackers/1",
                (false, false) => "chacha8-u53/rules-static-direct/1",
            },
            Method::MtdnaDirect {} => "chacha8-u53/rules-static-direct/1",
            Method::NextReaction {} => "chacha8-u53/rules-next-reaction/1",
            Method::Thinning {} => "chacha8-u53/rules-atp-thinning/1",
            Method::Tau { .. } => "chacha8-u53/rules-approximate-cgp-poisson/2",
            Method::Binomial { .. } => "chacha8-u53/rules-approximate-binomial-split/2",
            Method::Hybrid { .. } => "chacha8-u53/rules-approximate-hybrid-dopri54/1",
        }
        .into(),
        root_seed: request.ensemble.root_seed,
        types: model
            .types
            .iter()
            .map(|t| TypeHeader {
                name: t.name.clone(),
                attributes: t.attributes.clone(),
            })
            .collect(),
        rules: model.rules.iter().map(|r| r.name.clone()).collect(),
        observables: model.observables.iter().map(|o| o.name.clone()).collect(),
        report_times: request.time.report_times.clone(),
        termination,
        requested_replicates: n,
        completed_replicates: replicates
            .iter()
            .filter(|r| r.termination.complete())
            .count() as u64,
        stopped_replicates: replicates
            .iter()
            .filter(|r| r.termination == Termination::Stopped)
            .count() as u64,
        events: replicates.iter().map(|r| r.events).sum(),
        replicates,
        paths,
        final_instances,
        statistics,
        diagnostics,
        observations: reducer.observations,
        method: (request.schema >= RunVersion::V05).then(|| request.method.unwrap_or_default()),
        continuous: if request.schema >= RunVersion::V05 {
            model
                .atp
                .iter()
                .map(|_| "ATP".to_string())
                .chain(model.continuous.iter().map(|v| v.name.clone()))
                .collect()
        } else {
            Vec::new()
        },
    };
    let limit = request.budget.max_output_bytes.min(MAX_RESULT_BYTES);
    if !fits_json(&result, limit) {
        result.termination = Termination::OutputLimit;
        result.statistics = None;
        result.paths.clear();
        result.observations.clear();
        result.final_instances.clear();
        result.diagnostics.push(Diagnostic::error(
            "CRN-R045",
            "/outputs",
            "rules output exceeds its byte budget",
            "Increase max_output_bytes or keep fewer paths; statistics, paths and instances have been removed.",
        ));
        if !fits_json(&result, limit) {
            return Err(Diagnostic::error(
                "CRN-R045",
                "/budget/max_output_bytes",
                "replicate records exceed the output budget",
                "Increase max_output_bytes.",
            ));
        }
    }
    Ok(result)
}

/// Prepares and executes in one call: [`prepare`] followed by [`execute`].
///
/// # Errors
///
/// The diagnostics of [`prepare`], or the single `CRN-R045` error of
/// [`execute`].
pub fn solve(
    model_bytes: &[u8],
    request_bytes: &[u8],
    control: &Control,
) -> Result<RulesResult, Vec<Diagnostic>> {
    execute(&prepare(model_bytes, request_bytes)?, control).map_err(|e| vec![e])
}

pub(crate) fn read(path: &Path) -> Result<Vec<u8>, Diagnostic> {
    let error = |message: String| {
        Diagnostic::error(
            "CRN-E030",
            "",
            message,
            "Provide a readable regular JSON file.",
        )
    };
    let meta = std::fs::metadata(path).map_err(|e| error(format!("{}: {e}", path.display())))?;
    if !meta.is_file() {
        return Err(error(format!("{} is not a regular file", path.display())));
    }
    if meta.len() > crate::strict_json::MAX_DOCUMENT_BYTES as u64 {
        return Err(Diagnostic::error(
            "CRN-E004",
            "",
            format!("{} exceeds the document size limit", path.display()),
            "Keep input documents below 16 MiB.",
        ));
    }
    std::fs::read(path).map_err(|e| error(format!("{}: {e}", path.display())))
}

/// Runs a request file against the model it references, which must be the
/// supplied model file below the request's directory.
///
/// The request's `model.path` is resolved against the request file's
/// directory and must name the same file as `model_path`, inside that
/// directory. This is the `crn-craft rules MODEL REQUEST` command.
///
/// # Errors
///
/// `CRN-E030` for a file that cannot be read or is not a regular file,
/// `CRN-E004` for a file above 16 MiB, `CRN-E001` to `CRN-E004` for a
/// request that does not parse, `CRN-R020` at `/model/path` when the request
/// does not reference `model_path`, and otherwise the errors of [`solve`].
pub fn solve_files(
    model_path: &Path,
    request_path: &Path,
    control: &Control,
) -> Result<RulesResult, Vec<Diagnostic>> {
    let model = read(model_path).map_err(|e| vec![e])?;
    let bytes = read(request_path).map_err(|e| vec![e])?;
    let request: RunRequest = crate::strict_json::parse(&bytes).map_err(|e| vec![e])?;
    if !references_model(model_path, request_path, &request.model.path) {
        return Err(vec![request_error(
            "/model/path",
            "rules request does not identify the supplied model below its directory",
            "Use the referenced model within the request directory.",
        )]);
    }
    solve(&model, &bytes, control)
}

/// Whether a request file's `model.path` names `model_path`, inside the
/// request file's directory.
pub(crate) fn references_model(model_path: &Path, request_path: &Path, relative: &str) -> bool {
    let parent = request_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let resolved = (
        parent.canonicalize(),
        parent.join(relative).canonicalize(),
        model_path.canonicalize(),
    );
    crate::validate::valid_relative_path(relative)
        && matches!(resolved, (Ok(ref base), Ok(ref reference), Ok(ref supplied)) if reference.starts_with(base) && reference == supplied)
}

/// Reads a model file and checks it like [`check_model`]; this is the
/// `crn-craft rules-check MODEL` command.
///
/// A file that cannot be read gives an invalid report with one `CRN-E030`
/// diagnostic (or `CRN-E004` above 16 MiB) and no `model_sha256`.
///
/// ```
/// use std::path::Path;
///
/// let report = crn_craft::rules::check_model_file(Path::new("no/such/model.rules.json"));
/// assert!(!report.valid);
/// assert_eq!(report.model_sha256, None);
/// assert_eq!(report.diagnostics[0].code, "CRN-E030");
/// ```
pub fn check_model_file(model_path: &Path) -> RulesCheck {
    match read(model_path) {
        Ok(bytes) => check_model(&bytes),
        Err(e) => RulesCheck {
            schema: CheckVersion::V01,
            valid: false,
            model_sha256: None,
            diagnostics: vec![e],
        },
    }
}

/// Discovery entry for the native CLI's `describe` (`describe/0.14`).
///
/// The native CLI embeds it as `dynamic_species`; Python's
/// `crn_craft.rules.describe()` returns it. It names the status, the
/// commands, the accepted contract versions (`accepted_models`,
/// `accepted_requests`), the request version that introduced each option
/// or default, the methods (and which of them are approximate), the
/// converter from network models (`convert`), the 1.0 contracts
/// (`contracts_1_0`), the tree sampler identity and the model and run
/// limits.
///
/// ```
/// let entry = crn_craft::rules::describe();
/// assert_eq!(entry["status"], "experimental");
/// assert_eq!(entry["static_direct_defaults_request"], "crn-craft.rules-run/0.6");
/// assert!(entry["methods"].as_array().unwrap().contains(&"static_direct".into()));
/// assert_eq!(entry["limits"]["threads"], crn_craft::rules::MAX_WORKERS);
/// ```
pub fn describe() -> Value {
    serde_json::json!({
        "status": "experimental",
        "crate": "crn-craft-rules",
        "commands": ["rules-check", "rules"],
        "model": "crn-craft.rules-model/0.1",
        "accepted_models": ["crn-craft.rules-model/0.1","crn-craft.rules-model/0.2","crn-craft.rules-model/0.3","crn-craft.rules-model/0.4"],
        "check": "crn-craft.rules-check/0.1",
        "request": "crn-craft.rules-run/0.1",
        "result": "crn-craft.rules-result/0.1",
        "extended_request": "crn-craft.rules-run/0.2",
        "extended_result": "crn-craft.rules-result/0.2",
        "static_direct_options_request": "crn-craft.rules-run/0.4",
        "static_direct_defaults_request": "crn-craft.rules-run/0.6",
        "accepted_requests": [
            "crn-craft.rules-run/0.1", "crn-craft.rules-run/0.2", "crn-craft.rules-run/0.3",
            "crn-craft.rules-run/0.4", "crn-craft.rules-run/0.5", "crn-craft.rules-run/0.6"
        ],
        "mtdna_direct_request": "crn-craft.rules-run/0.3",
        "methods": ["tree", "static_direct", "mtdna_direct", "next_reaction", "tau", "binomial", "thinning", "hybrid"],
        "hybrid_request": "crn-craft.rules-run/0.5",
        "hybrid_result": "crn-craft.rules-result/0.3",
        "python": "crn_craft.rules",
        "approximate_methods": ["tau", "binomial", "hybrid"],
        "convert": {
            "command": "convert", "python": "crn_craft.convert",
            "from": "crn-craft.model/0.1", "to": convert::converted_schema()
        },
        "contracts_1_0": {
            "model": "crn-craft.model/1.0", "request": "crn-craft.run/1.0",
            "result": "crn-craft.result/1.0", "check": "crn-craft.check/1.0",
            "fsp": {"method": "fsp", "result": "crn-craft.fsp-result/1.0", "schema": "fsp-result-1.0", "python": "crn_craft.rules.fsp"},
            "means": {
                "model": "crn-craft.rules-model/0.4", "request": "crn-craft.rules-run/0.6",
                "result": "crn-craft.rules-result/0.3", "check": "crn-craft.rules-check/0.1"
            },
            "request_accepts_models": [
                "crn-craft.rules-model/0.1", "crn-craft.rules-model/0.2", "crn-craft.rules-model/0.3",
                "crn-craft.rules-model/0.4", "crn-craft.model/1.0"
            ],
            "schemas": ["model-1.0", "run-1.0", "result-1.0", "check-1.0"],
            "limits": {
                "types": model::MAX_TYPES_1_0, "parameters": model::MAX_PARAMETERS_1_0,
                "observables": model::MAX_OBSERVABLES_1_0, "rules": model::MAX_RULES_1_0
            },
            "quantiles": {
                "request": "outputs.quantiles", "result": "statistics.quantiles",
                "max": MAX_QUANTILES, "definition": "numpy-linear-type-7"
            }
        },
        "sampler": SAMPLER,
        "reference": "docs/rules-language.md",
        "limits": {
            "types": model::MAX_TYPES, "attributes": model::MAX_ATTRIBUTES,
            "parameters": model::MAX_PARAMETERS, "observables": model::MAX_OBSERVABLES,
            "rules": model::MAX_RULES, "products": model::MAX_PRODUCTS,
            "initial_entries": model::MAX_INITIAL, "expression_bytes": model::MAX_EXPRESSION_BYTES,
            "replicates": MAX_REPLICATES, "report_times": MAX_REPORT_TIMES,
            "threads": MAX_WORKERS, "result_bytes": MAX_RESULT_BYTES
        }
    })
}

/// JSON Schemas generated from the Rust contract types.
///
/// `kind` is `rules-model` ([`Model`]), `rules-check` ([`RulesCheck`]),
/// `rules-run` ([`RunRequest`]) or `rules-result` ([`RulesResult`]), or one
/// of the 1.0 contracts generated from the same types: `model-1.0`,
/// `check-1.0`, `run-1.0` and `result-1.0`. Any other kind gives `None`;
/// `fsp-result-1.0` comes from [`crate::schema`]. The `run-1.0` schema
/// lists two request shapes under `oneOf`: the simulation request
/// ([`RunRequest`], `definitions/RunRequest`) and the finite state
/// projection request ([`crate::fsp::v1::FspRequest`]).
/// Each schema lists only its own versions in its `schema` field: the
/// `rules-*` schemas the 0.x versions, the 1.0 schemas the 1.0 version,
/// with a description of that contract. The schemas are JSON Schema draft
/// 07 and are committed under `schemas/` in the repository. They express
/// the shape and some local bounds; [`prepare`] enforces the remaining
/// rules.
///
/// ```
/// let run = crn_craft::rules::schema("rules-run").unwrap();
/// assert_eq!(run["$schema"], "http://json-schema.org/draft-07/schema#");
/// assert!(run["properties"]["budget"].is_object());
/// assert!(crn_craft::rules::schema("model").is_none());
/// let unified = crn_craft::rules::schema("run-1.0").unwrap();
/// assert_eq!(unified["definitions"]["RunVersion"]["enum"][0], "crn-craft.run/1.0");
/// ```
pub fn schema(kind: &str) -> Option<Value> {
    let (schema, definition, unified) = match kind {
        "rules-model" => (schema_for!(Model), "ModelVersion", None),
        "rules-check" => (schema_for!(RulesCheck), "CheckVersion", None),
        "rules-run" => (schema_for!(RunRequest), "RunVersion", None),
        "rules-result" => (schema_for!(RulesResult), "ResultVersion", None),
        "model-1.0" => (schema_for!(Model), "ModelVersion", Some(MODEL_1_0)),
        "check-1.0" => (schema_for!(RulesCheck), "CheckVersion", Some(CHECK_1_0)),
        "run-1.0" => (schema_for!(RunRequest), "RunVersion", Some(RUN_1_0)),
        "result-1.0" => (schema_for!(RulesResult), "ResultVersion", Some(RESULT_1_0)),
        _ => return None,
    };
    let mut value = serde_json::to_value(schema).ok()?;
    let definitions = value.get_mut("definitions")?.as_object_mut()?;
    if unified.is_none() {
        // A rules-* document cannot hold the fields of the 1.0 contracts.
        for (name, property) in ONLY_1_0 {
            if let Some(properties) = definitions
                .get_mut(*name)
                .and_then(|d| d.get_mut("properties"))
                .and_then(Value::as_object_mut)
            {
                properties.remove(*property);
            }
        }
        for name in ONLY_1_0_DEFINITIONS {
            definitions.remove(*name);
        }
    }
    let version = definitions.get_mut(definition)?;
    version
        .get_mut("enum")?
        .as_array_mut()?
        .retain(|v| v.as_str().is_some_and(is_unified) == unified.is_some());
    if let Some((description, version_description)) = unified {
        version["description"] = version_description.into();
        value["description"] = description.into();
        // The rules-model text keeps a rustdoc path of 0.6.0 (see the field).
        if let Some(text) = value.pointer_mut("/properties/continuous/description") {
            *text = text
                .as_str()?
                .replace(
                    "See [`crate::continuous`].",
                    "See `crn_craft::rules::continuous`.",
                )
                .into();
        }
        // Field texts that state the limits of the 0.x model versions.
        let fields: &[(&str, &str)] = if kind == "model-1.0" {
            &MODEL_1_0_FIELDS
        } else {
            &[]
        };
        for (pointer, text) in fields {
            if let Some(description) = value.pointer_mut(pointer) {
                *description = (*text).into();
            }
        }
    }
    if kind == "run-1.0" {
        // The fsp method has a request shape of its own.
        return crate::fsp::v1::unified_run_schema(value);
    }
    Some(value)
}

/// Descriptions of `crn-craft.model/1.0` fields whose `rules-model` text
/// states a 0.x limit.
const MODEL_1_0_FIELDS: [(&str, &str); 6] = [
    (
        "/definitions/Rule/properties/reactants/description",
        "Reactant patterns, any number (see `count`). Reactants are consumed; list a `copy` product to keep one.",
    ),
    (
        "/definitions/Pattern/properties/label/description",
        "Name of the matched instance in this pattern's expressions and in the rule's products, for example `x` in `x.size`. The patterns of a rule need different labels.",
    ),
    (
        "/definitions/Pattern/properties/type/description",
        "Name of a declared type. Two different patterns over the same type are not supported (`CRN-R013`); identical ones form pairs, or with more molecules a higher-order rule.",
    ),
    (
        "/properties/types/description",
        "Species types, 1 to 4096, with unique names.",
    ),
    (
        "/properties/observables/description",
        "Observables by name, at most 4096. Results list them sorted by name.",
    ),
    (
        "/properties/rules/description",
        "Rules, 1 to 65536, with unique names. Results count events per rule in this order.",
    ),
];

/// Properties (definition, property) of the shared contract types that only
/// the 1.0 contracts have: the `rules-*` schemas leave them out.
const ONLY_1_0: &[(&str, &str)] = &[
    ("Outputs", "quantiles"),
    ("Statistics", "quantiles"),
    ("Pattern", "count"),
];
/// Definitions that only the 1.0 schemas use.
const ONLY_1_0_DEFINITIONS: &[&str] = &["Quantile"];

/// Whether a `schema` string names a 1.0 contract.
fn is_unified(version: &str) -> bool {
    version.ends_with("/1.0")
}

/// Descriptions of a 1.0 schema and of its version definition.
type Descriptions = (&'static str, &'static str);

const MODEL_1_0: Descriptions = (
    "A model: `crn-craft.model/1.0`, the unified model contract. A document means what the same document means as `crn-craft.rules-model/0.4`.\n\nUnknown fields and duplicate keys are rejected. Validation also requires at least one type and one rule. It runs only with `crn-craft.run/1.0` requests; `crn-craft convert` writes network models in this contract.",
    "The version of a model, its `schema` field.",
);
const CHECK_1_0: Descriptions = (
    "The report of a model check of a `crn-craft.model/1.0` document: `crn-craft.check/1.0`, with the fields of `crn-craft.rules-check/0.1`.\n\nThe check stops at the first stage that fails (parsing, size limits, names and types, then the rest of the model), so `diagnostics` lists the errors of that stage.",
    "The version of a model check report, its `schema` field.",
);
const RUN_1_0: Descriptions = (
    "A run request: `crn-craft.run/1.0`, the unified request contract, with the fields and defaults of `crn-craft.rules-run/0.6` plus `outputs.quantiles`.\n\nIt binds one model, `crn-craft.model/1.0` or any `crn-craft.rules-model` version, by hash and states the time grid, the ensemble, the outputs and an explicit budget. Unknown fields and duplicate keys are rejected. It returns `crn-craft.result/1.0`.",
    "The version of a run request, its `schema` field.",
);
const RESULT_1_0: Descriptions = (
    "The result of a `crn-craft.run/1.0` request: `crn-craft.result/1.0`, with the fields of `crn-craft.rules-result/0.3` plus `statistics.quantiles`.\n\nA result is returned even when replicates hit limits or fail: check `termination`, the per-replicate records and `diagnostics`. Matrices are indexed by report time, then observable.",
    "The version of a result, its `schema` field.",
);

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_feedback;
#[cfg(test)]
mod tests_review;
