//! The network model contract `crn-craft.model/0.1`, which only the
//! converter reads: its types (they generate the schema that `crn-craft
//! describe model` prints), strict parsing and validation.
//!
//! Until release 0.7.0 these were `crn_craft::contracts` and
//! `crn_craft::validate`, shared with the network engine.
use crate::diagnostic::{Diagnostic, pointer};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Schema tag of a model document.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub(super) enum ModelVersion {
    #[doc = concat!("`crn-craft.model/0.1`, the only model version.")]
    #[serde(rename = "crn-craft.model/0.1")]
    V01,
}

/// How a model's numbers are read: integer molecule counts and complete stochastic propensities.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub(super) enum Convention {
    #[doc = concat!("`stochastic_counts`: populations are molecule counts, time is in seconds and")]
    #[doc = concat!("each rate law gives a channel's full propensity in events per second.")]
    StochasticCounts,
}

/// A fixed reaction network: species with initial counts, named parameters and reactions
/// with rate laws.
/// Declaration order fixes every species and reaction index in results and analyses.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub(super) struct Model {
    #[doc = concat!("Always `crn-craft.model/0.1`.")]
    pub(super) schema: ModelVersion,
    /// Model identifier: 1–64 ASCII letters, digits or underscores, not starting with a digit.
    pub(super) name: String,
    #[doc = concat!("Always `stochastic_counts`.")]
    pub(super) convention: Convention,
    /// Counts per chosen concentration unit. Does not rescale SSA coefficients.
    #[schemars(schema_with = "positive_number")]
    pub(super) omega: f64,
    /// Finite named values that rate laws reference by name; at most 16384. Names follow the
    /// identifier rules and must not collide with species names.
    pub(super) parameters: BTreeMap<String, f64>,
    /// Species in declaration order, 1–4096 entries with unique names.
    pub(super) species: Vec<Species>,
    /// Reactions in declaration order, at most 65536; an empty list gives an absorbing model.
    pub(super) reactions: Vec<Reaction>,
}

/// One species: a unique name, its initial molecule count and an optional safety limit.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub(super) struct Species {
    /// Identifier, unique among species and distinct from every parameter name.
    pub(super) name: String,
    /// Initial molecule count, a non-negative 64-bit integer; counts above 2^53 stay exact.
    #[schemars(range(min = 0))]
    pub(super) initial_count: i64,
    /// A safety stop, never a reflecting boundary or an FSP truncation bound.
    /// When present it is at least `initial_count`. An update that would exceed it is never
    /// applied; a stochastic replicate stops with `population_limit` instead.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(super) population_limit: Option<i64>,
}

/// One reaction channel: stoichiometry by species name and a rate law.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub(super) struct Reaction {
    /// Identifier, unique among reactions.
    pub(super) name: String,
    /// Consumed species and coefficients from 1 to 2147483647; omit zeros. An empty map is a
    /// source. A channel can fire only while every reactant count covers its coefficient.
    pub(super) reactants: BTreeMap<String, i64>,
    /// Produced species and coefficients from 1 to 2147483647; omit zeros. An empty map is a sink.
    pub(super) products: BTreeMap<String, i64>,
    #[doc = concat!("Propensity law of the channel.")]
    pub(super) rate: Rate,
}

/// A literal or the exact name of a model parameter; no string evaluation.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub(super) enum Scalar {
    /// A finite number used as given.
    Literal(f64),
    /// The name of an entry in the model's `parameters`; an unknown name is `CRN-E014`.
    Parameter(String),
}

/// Custom laws are complete propensities, not mass-action multipliers.
/// All laws require sufficient reactants before firing (implemented in M1/M2).
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(super) enum Rate {
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

/// Most species a model may declare (`CRN-E004`); at least one is required.
const MAX_SPECIES: usize = 4096;
/// Most parameters a model may declare (`CRN-E004`).
const MAX_PARAMETERS: usize = 16384;
/// Most reactions a model may declare (`CRN-E004`); zero is allowed.
const MAX_REACTIONS: usize = 65536;

/// Strictly parses a network model: syntax and shape, not semantics.
/// A shape error (`CRN-E003`) in a document whose `schema` names a contract
/// of the unified engine, a `crn-craft.rules-*` or a 1.0 contract, becomes
/// `CRN-E032` at `/schema`: such a document never parses as a network model,
/// and the strict parser would report the first field it does not know.
pub(super) fn parse(bytes: &[u8]) -> Result<Model, Diagnostic> {
    crate::strict_json::parse(bytes).map_err(|error| {
        if error.code != "CRN-E003" {
            return error;
        }
        match crate::schema_field(bytes) {
            Some(schema) if schema.starts_with("crn-craft.rules-") || schema.ends_with("/1.0") => {
                Diagnostic::error(
                    "CRN-E032",
                    "/schema",
                    format!(
                        "`{schema}` is not a network model; the converter reads `crn-craft.model/0.1` network models"
                    ),
                    "This model needs no conversion: check it with `crn-craft check` and run it with `crn-craft run`. `crn-craft convert` writes a network model as a crn-craft.model/1.0 model.",
                )
            }
            _ => error,
        }
    })
}

/// The generated JSON Schema of the network model (`crn-craft describe model`).
pub(super) fn schema() -> serde_json::Value {
    serde_json::to_value(schemars::schema_for!(Model)).expect("a schema serializes")
}

fn identifier(name: &str) -> bool {
    let mut bytes = name.bytes();
    name.len() <= 64
        && bytes
            .next()
            .is_some_and(|x| x.is_ascii_alphabetic() || x == b'_')
        && bytes.all(|x| x.is_ascii_alphanumeric() || x == b'_')
}

fn name_check(name: &str, path: &str, out: &mut Vec<Diagnostic>) {
    if !identifier(name) {
        out.push(Diagnostic::error(
            "CRN-E010",
            path,
            "invalid identifier",
            "Use 1–64 ASCII letters, digits or underscores, starting with a letter or underscore.",
        ));
    }
}

fn species_ref(name: &str, path: &str, names: &BTreeSet<&str>, out: &mut Vec<Diagnostic>) {
    if !names.contains(name) {
        out.push(Diagnostic::error(
            "CRN-E012",
            path,
            format!("unknown species {name:?}"),
            "Declare the species or correct its name.",
        ));
    }
}

fn scalar(
    value: &Scalar,
    path: &str,
    positive: bool,
    parameters: &BTreeMap<String, f64>,
    out: &mut Vec<Diagnostic>,
) {
    let number = match value {
        Scalar::Literal(x) => *x,
        Scalar::Parameter(name) => match parameters.get(name) {
            Some(x) => *x,
            None => {
                out.push(Diagnostic::error(
                    "CRN-E014",
                    path,
                    format!("unknown parameter {name:?}"),
                    "Declare the parameter or provide a numeric literal.",
                ));
                return;
            }
        },
    };
    if !number.is_finite() || number < 0.0 || (positive && number == 0.0) {
        out.push(Diagnostic::error(
            "CRN-E015",
            path,
            "rate argument outside its finite numeric domain",
            if positive {
                "Use a finite, strictly positive value."
            } else {
                "Use a finite, non-negative value."
            },
        ));
    }
}

/// The errors of a parsed network model, in a deterministic order: the
/// structural and domain checks of release 0.6.0 without its warnings
/// (`crn-craft analyze` reports those for the converted model).
///
/// Checks identifiers (1–64 ASCII letters, digits or underscores, not starting
/// with a digit), unique names, a species/parameter namespace without
/// collisions, known species and parameters, stoichiometry in 1..=2147483647,
/// finite parameters, rate arguments in their domains (coefficients, `vmax`
/// and `basal` non-negative; `k`, `km` and `n` positive), a finite positive
/// `omega`, non-negative initial counts and population limits at least the
/// initial count.
pub(super) fn validate(model: &Model) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    name_check(&model.name, "/name", &mut out);
    if !model.omega.is_finite() || model.omega <= 0.0 {
        out.push(Diagnostic::error(
            "CRN-E015",
            "/omega",
            "invalid system size",
            "Use a finite positive counts-per-concentration scale.",
        ));
    }
    if model.species.is_empty()
        || model.species.len() > MAX_SPECIES
        || model.parameters.len() > MAX_PARAMETERS
        || model.reactions.len() > MAX_REACTIONS
    {
        out.push(Diagnostic::error(
            "CRN-E004",
            "",
            "model collection limit exceeded or no species",
            "Use 1–4096 species, at most 16384 parameters and 65536 reactions.",
        ));
    }
    let mut names = BTreeSet::new();
    for (i, item) in model.species.iter().enumerate() {
        let base = format!("/species/{i}");
        name_check(&item.name, &pointer(&base, "name"), &mut out);
        if !names.insert(item.name.as_str()) {
            out.push(Diagnostic::error(
                "CRN-E011",
                &pointer(&base, "name"),
                "duplicate species name",
                "Give each species a unique name.",
            ));
        }
        if item.initial_count < 0 {
            out.push(Diagnostic::error(
                "CRN-E016",
                &pointer(&base, "initial_count"),
                "negative initial population",
                "Use a non-negative i64 count.",
            ));
        }
        if item
            .population_limit
            .is_some_and(|limit| limit < 0 || limit < item.initial_count)
        {
            out.push(Diagnostic::error(
                "CRN-E016",
                &pointer(&base, "population_limit"),
                "population limit below the initial population",
                "Use a non-negative limit at least as large as initial_count.",
            ));
        }
    }
    for (name, value) in &model.parameters {
        let path = pointer("/parameters", name);
        name_check(name, &path, &mut out);
        if names.contains(name.as_str()) {
            out.push(Diagnostic::error(
                "CRN-E011",
                &path,
                "parameter and species share a name",
                "Keep parameter and species namespaces disjoint.",
            ));
        }
        if !value.is_finite() {
            out.push(Diagnostic::error(
                "CRN-E015",
                &path,
                "non-finite parameter",
                "Use a finite number.",
            ));
        }
    }
    let mut reaction_names = BTreeSet::new();
    for (i, reaction) in model.reactions.iter().enumerate() {
        let base = format!("/reactions/{i}");
        name_check(&reaction.name, &pointer(&base, "name"), &mut out);
        if !reaction_names.insert(&reaction.name) {
            out.push(Diagnostic::error(
                "CRN-E011",
                &pointer(&base, "name"),
                "duplicate reaction name",
                "Give each reaction a unique name.",
            ));
        }
        for (label, side) in [
            ("reactants", &reaction.reactants),
            ("products", &reaction.products),
        ] {
            for (name, value) in side {
                let path = pointer(&pointer(&base, label), name);
                species_ref(name, &path, &names, &mut out);
                if *value <= 0 || *value > i64::from(i32::MAX) {
                    out.push(Diagnostic::error(
                        "CRN-E013",
                        &path,
                        "invalid stoichiometry",
                        "Use positive integers up to 2147483647; omit zero entries.",
                    ));
                }
            }
        }
        let base = pointer(&base, "rate");
        match &reaction.rate {
            Rate::MassAction { coefficient } => {
                scalar(
                    coefficient,
                    &pointer(&base, "coefficient"),
                    false,
                    &model.parameters,
                    &mut out,
                );
            }
            Rate::HillActivation {
                regulator,
                vmax,
                basal,
                k,
                n,
            }
            | Rate::HillRepression {
                regulator,
                vmax,
                basal,
                k,
                n,
            } => {
                species_ref(regulator, &pointer(&base, "regulator"), &names, &mut out);
                for (key, value, positive) in [
                    ("vmax", vmax, false),
                    ("basal", basal, false),
                    ("k", k, true),
                    ("n", n, true),
                ] {
                    scalar(
                        value,
                        &pointer(&base, key),
                        positive,
                        &model.parameters,
                        &mut out,
                    );
                }
            }
            Rate::MichaelisMenten {
                substrate,
                vmax,
                km,
            } => {
                species_ref(substrate, &pointer(&base, "substrate"), &names, &mut out);
                scalar(
                    vmax,
                    &pointer(&base, "vmax"),
                    false,
                    &model.parameters,
                    &mut out,
                );
                scalar(km, &pointer(&base, "km"), true, &model.parameters, &mut out);
            }
        }
    }
    out
}

fn positive_number(generator: &mut schemars::r#gen::SchemaGenerator) -> schemars::schema::Schema {
    let mut schema = f64::json_schema(generator).into_object();
    schema.number().exclusive_minimum = Some(0.0);
    schema.into()
}
