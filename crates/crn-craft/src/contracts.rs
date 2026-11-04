//! Versioned data contracts. Validation is a separate, mandatory step.
//!
//! These Rust types define the fixed-network input contracts and generate
//! their JSON Schemas (`crn-craft describe model`, `describe run`,
//! `describe result-v1`), so the doc comments below double as schema
//! descriptions:
//!
//! | Schema | Type | Role |
//! |---|---|---|
//! | `crn-craft.model/0.1` | [`Model`] | Species, parameters and reactions with rate laws |
//! | `crn-craft.run/0.1` | [`Run`] | Direct SSA or FSP request bound to one model file |
//! | `crn-craft.result/0.1` | [`RunResult`] | Reserved envelope; no solver emits it |
//!
//! Unknown fields are rejected everywhere. Read documents with
//! [`crate::parse_model`] and [`crate::parse_run`], which also reject
//! duplicate keys, and validate them with [`crate::check_bytes`] or the
//! functions in [`crate::validate`]; deserializing with plain `serde_json`
//! skips both safeguards. Other requests reuse parts of these types:
//! [`crate::ssa::Request`] embeds [`ModelRef`], [`TimeGrid`], [`Ensemble`] and
//! [`Outputs`], and [`crate::tau::Request`] embeds a whole [`Run`]. Field
//! semantics are specified in
//! [docs/contracts.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/contracts.md).
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of a model document.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ModelVersion {
    #[doc = concat!("`crn-craft.model/0.1`, the only model version.")]
    #[serde(rename = "crn-craft.model/0.1")]
    V01,
}

/// Schema tag of a fixed-network run request.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum RunVersion {
    #[doc = concat!("`crn-craft.run/0.1`, the only run version.")]
    #[serde(rename = "crn-craft.run/0.1")]
    V01,
}

/// Schema tag of the reserved result envelope, which no solver emits.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ResultVersion {
    #[doc = concat!("`crn-craft.result/0.1`; executed runs report `crn-craft.result/0.2`.")]
    #[serde(rename = "crn-craft.result/0.1")]
    V01,
}

/// How a model's numbers are read: integer molecule counts and complete stochastic propensities.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Convention {
    #[doc = concat!("`stochastic_counts`: populations are molecule counts, time is in seconds and")]
    #[doc = concat!("each rate law gives a channel's full propensity in events per second.")]
    StochasticCounts,
}

/// A fixed reaction network: species with initial counts, named parameters and reactions
/// with rate laws.
/// Declaration order fixes every species and reaction index in results and analyses.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Model {
    #[doc = concat!("Always `crn-craft.model/0.1`.")]
    pub schema: ModelVersion,
    /// Model identifier: 1–64 ASCII letters, digits or underscores, not starting with a digit.
    pub name: String,
    #[doc = concat!("Always `stochastic_counts`.")]
    pub convention: Convention,
    /// Counts per chosen concentration unit. Does not rescale SSA coefficients.
    #[schemars(schema_with = "positive_number")]
    pub omega: f64,
    /// Finite named values that rate laws reference by name; at most 16384. Names follow the
    /// identifier rules and must not collide with species names.
    pub parameters: BTreeMap<String, f64>,
    /// Species in declaration order, 1–4096 entries with unique names.
    pub species: Vec<Species>,
    /// Reactions in declaration order, at most 65536; an empty list gives an absorbing model.
    pub reactions: Vec<Reaction>,
}

/// One species: a unique name, its initial molecule count and an optional safety limit.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Species {
    /// Identifier, unique among species and distinct from every parameter name.
    pub name: String,
    /// Initial molecule count, a non-negative 64-bit integer; counts above 2^53 stay exact.
    #[schemars(range(min = 0))]
    pub initial_count: i64,
    /// A safety stop, never a reflecting boundary or an FSP truncation bound.
    /// When present it is at least `initial_count`. An update that would exceed it is never
    /// applied; a stochastic replicate stops with `population_limit` instead.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub population_limit: Option<i64>,
}

/// One reaction channel: stoichiometry by species name and a rate law.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Reaction {
    /// Identifier, unique among reactions.
    pub name: String,
    /// Consumed species and coefficients from 1 to 2147483647; omit zeros. An empty map is a
    /// source. A channel can fire only while every reactant count covers its coefficient.
    pub reactants: BTreeMap<String, i64>,
    /// Produced species and coefficients from 1 to 2147483647; omit zeros. An empty map is a sink.
    pub products: BTreeMap<String, i64>,
    #[doc = concat!("Propensity law of the channel.")]
    pub rate: Rate,
}

/// A literal or the exact name of a model parameter; no string evaluation.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum Scalar {
    /// A finite number used as given.
    Literal(f64),
    /// The name of an entry in the model's `parameters`; an unknown name is `CRN-E014`.
    Parameter(String),
}

/// Custom laws are complete propensities, not mass-action multipliers.
/// All laws require sufficient reactants before firing (implemented in M1/M2).
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Rate {
    /// `mass_action`: c × Π binomial(x, ν) over the reactants, with coefficient c ≥ 0, so `2A`
    /// gives c·x(x−1)/2 and an empty reactant side gives c.
    MassAction {
        #[doc = concat!("Stochastic rate constant c ≥ 0 in count units; `omega` does not")]
        #[doc = concat!("rescale it.")]
        coefficient: Scalar,
    },
    /// `hill_activation`: basal + vmax·xⁿ/(kⁿ + xⁿ), with x the regulator count;
    /// vmax, basal ≥ 0 and k, n > 0.
    HillActivation {
        /// Species whose count x drives the law.
        regulator: String,
        #[doc = concat!("Largest induced propensity, ≥ 0, in events per second.")]
        vmax: Scalar,
        #[doc = concat!("Propensity added at every state, ≥ 0, in events per second.")]
        basal: Scalar,
        #[doc = concat!("Half-saturation count, > 0.")]
        k: Scalar,
        #[doc = concat!("Hill exponent, > 0; it need not be an integer.")]
        n: Scalar,
    },
    /// `hill_repression`: basal + vmax·kⁿ/(kⁿ + xⁿ), with x the regulator count;
    /// vmax, basal ≥ 0 and k, n > 0.
    HillRepression {
        /// Species whose count x drives the law.
        regulator: String,
        #[doc = concat!("Propensity at x = 0 above `basal`, ≥ 0, in events per second.")]
        vmax: Scalar,
        #[doc = concat!("Propensity added at every state, ≥ 0, in events per second.")]
        basal: Scalar,
        #[doc = concat!("Half-repression count, > 0.")]
        k: Scalar,
        #[doc = concat!("Hill exponent, > 0; it need not be an integer.")]
        n: Scalar,
    },
    /// `michaelis_menten`: vmax·x/(km + x), with x the substrate count; vmax ≥ 0 and km > 0.
    MichaelisMenten {
        /// Species whose count enters the law.
        substrate: String,
        #[doc = concat!("Saturated propensity, ≥ 0, in events per second.")]
        vmax: Scalar,
        #[doc = concat!("Michaelis constant in molecule counts, > 0.")]
        km: Scalar,
    },
}

/// Reference from a request to its model file: a relative path and the SHA-256 of the exact
/// model bytes.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelRef {
    /// Portable, relative to the run file's directory; no parent traversal.
    pub path: String,
    /// SHA-256 of the exact model input bytes, lowercase hexadecimal.
    #[schemars(regex(pattern = "^[0-9a-f]{64}$"))]
    pub sha256: String,
}

/// Solver of a run request, selected by `kind`.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    /// `ssa`: exact Direct stochastic simulation of an ensemble. Needs `ensemble` and
    /// `budget.max_events_per_replicate`.
    Ssa {
        #[doc = concat!("Event selection; `direct` is the only value.")]
        aggregator: Aggregator,
    },
    /// `fsp`: finite state projection of the probability distribution on a bounded box.
    /// Takes no ensemble, event cap or retained trajectories.
    Fsp {
        /// Inclusive upper count bound for every species, by name, each at least the initial
        /// count. Probability that leaves the box moves to an absorbing sink.
        bounds: BTreeMap<String, i64>,
        /// Most reachable states to enumerate, ≥ 1; execution also caps it at 100000.
        max_states: u64,
        /// Total Poisson-series truncation budget over the whole run, 0 < tolerance < 1.
        /// It does not bound floating-point error.
        tolerance: f64,
    },
}

/// Event selection of a `run/0.1` SSA request.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Aggregator {
    #[doc = concat!("`direct`: the Direct method. `simulate-cached` runs its dependency-cached")]
    #[doc = concat!("variant, which draws the same random words and returns the same numbers.")]
    Direct,
}

/// A fixed-network execution request bound to one model file by relative path and SHA-256.
/// `simulate`, `simulate-cached` and `fsp` read it, and tau requests embed it.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Run {
    #[doc = concat!("Always `crn-craft.run/0.1`.")]
    pub schema: RunVersion,
    #[doc = concat!("The model file this request belongs to.")]
    pub model: ModelRef,
    #[doc = concat!("`ssa` or `fsp`.")]
    pub method: Method,
    #[doc = concat!("Horizon and save grid, in seconds.")]
    pub time: TimeGrid,
    /// Replicate count and root seed; required for `ssa` and absent for `fsp`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensemble: Option<Ensemble>,
    #[doc = concat!("Statistics and retained trajectories to report.")]
    pub outputs: Outputs,
    #[doc = concat!("Event, time, thread and output limits.")]
    pub budget: Budget,
}

/// Time horizon and the times at which states are reported, in seconds.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TimeGrid {
    /// Horizon in seconds, finite and positive.
    #[schemars(schema_with = "positive_number")]
    pub end: f64,
    /// Explicit sorted grid including 0 and end; no accumulated step rounding.
    /// 2–1000000 strictly increasing values. The state saved at a time includes every event
    /// at or before that time.
    pub save_times: Vec<f64>,
}

/// Independent replicates of a stochastic run.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Ensemble {
    /// Number of independent replicates, ≥ 1, with ids 0 to replicates − 1.
    #[schemars(range(min = 1))]
    pub replicates: u64,
    /// Root seed. Each replicate draws from its own stream, derived from this seed and the
    /// replicate id, so results do not depend on the thread count.
    pub root_seed: u64,
}

/// Ensemble statistic to report per species and save time.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Statistic {
    #[doc = concat!("`mean`: the sample mean, with its standard error when there are two or more")]
    #[doc = concat!("replicates.")]
    Mean,
    #[doc = concat!("`variance`: the unbiased sample variance (denominator n − 1); needs at least")]
    #[doc = concat!("two replicates.")]
    Variance,
}

/// What a run reports besides the per-replicate accounting; at least one statistic or trajectory.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Outputs {
    /// Statistics to compute, each at most once. Stochastic methods report them only when
    /// every replicate completed.
    pub statistics: Vec<Statistic>,
    /// How many replicates, starting from id 0, return their saved states; at most
    /// `replicates`, and 0 for FSP.
    pub keep_trajectories: u64,
}

/// Resource limits of a run request.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    /// Required for SSA; omitted for FSP, which has max_states.
    /// A replicate that reaches the cap stops with `event_limit`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_events_per_replicate: Option<u64>,
    /// Optional wall-clock limit in seconds, finite and positive. It is polled cooperatively,
    /// so it is not a hard deadline; unfinished replicates stop with `wall_time_limit`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_s: Option<f64>,
    /// Worker allowance, 1–1024. Execution uses at most 64 workers and no more than the
    /// replicate count; completed results do not depend on it.
    #[schemars(range(min = 1, max = 1024))]
    pub threads: u32,
    /// Byte budget for the result, at least 8 × save_times × (1 + species × (statistics +
    /// kept trajectories)). For an archive it covers every archive file; execution never
    /// uses more than 256 MiB.
    #[schemars(range(min = 1))]
    pub max_output_bytes: u64,
}

/// Termination states planned for the reserved `result/0.1` envelope.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    #[doc = concat!("`completed`: every replicate reached the horizon.")]
    Completed,
    #[doc = concat!("`event_limit`: a replicate reached its event cap.")]
    EventLimit,
    #[doc = concat!("`wall_time_limit`: the wall-clock limit expired.")]
    WallTimeLimit,
    #[doc = concat!("`population_limit`: an update would have exceeded a population limit.")]
    PopulationLimit,
    #[doc = concat!("`output_limit`: the result did not fit its byte budget.")]
    OutputLimit,
    #[doc = concat!("`cancelled`: cancellation was requested.")]
    Cancelled,
    #[doc = concat!("`failed`: a numerical or state-update check failed.")]
    Failed,
}

/// Assessments planned for the reserved `result/0.1` envelope.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Assessment {
    #[doc = concat!("`exact_ensemble`: a fully completed exact-method ensemble.")]
    ExactEnsemble,
    #[doc = concat!("`fsp_approximation`: a finite state projection.")]
    FspApproximation,
    #[doc = concat!("`incomplete`: some requested work did not finish.")]
    Incomplete,
    #[doc = concat!("`failed_checks`: a numerical check failed.")]
    FailedChecks,
}

/// Reserved result envelope; M0 never fabricates a simulation result.
/// No solver emits it: executed runs report `crn-craft.result/0.2` instead.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunResult {
    #[doc = concat!("Always `crn-craft.result/0.1`.")]
    pub schema: ResultVersion,
    /// SHA-256 of the exact model bytes, lowercase hexadecimal.
    pub model_sha256: String,
    /// SHA-256 of the exact run bytes, lowercase hexadecimal.
    pub run_sha256: String,
    #[doc = concat!("How the run ended.")]
    pub termination: Termination,
    #[doc = concat!("What the result may be used for.")]
    pub assessment: Assessment,
    /// Replicates the run asked for.
    pub requested_replicates: u64,
    /// Replicates that reached the horizon.
    pub completed_replicates: u64,
    /// Events executed over all replicates.
    pub events: u64,
    /// Checks that support the assessment.
    pub evidence: Vec<Evidence>,
}

/// One check supporting an assessment of the reserved `result/0.1` envelope.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Evidence {
    /// Identifier of the check.
    pub check: String,
    /// Outcome of the check.
    pub passed: bool,
    /// Explanation of the outcome for readers.
    pub detail: String,
}

fn positive_number(generator: &mut schemars::r#gen::SchemaGenerator) -> schemars::schema::Schema {
    let mut schema = f64::json_schema(generator).into_object();
    schema.number().exclusive_minimum = Some(0.0);
    schema.into()
}
