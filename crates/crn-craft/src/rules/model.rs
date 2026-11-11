//! `crn-craft.rules-model/0.1`: species types with numeric attributes, rule
//! schemata, observables and an optional stop condition.
//!
//! [`Model`] is the JSON document; `crn-craft.rules-model/0.2` adds the
//! optional [`crate::rules::atp`] block, `0.3` whole-rule propensities
//! ([`Rule::propensity`]) and continuous variables, and `0.4` independent
//! continuous variables ([`crate::rules::continuous::ContinuousSpec::independent`]). A model is parsed strictly and validated by
//! [`crate::rules::check_model`] or [`crate::rules::prepare`], which compile every
//! expression. The [rules language
//! reference](https://github.com/dietmarwo/CRN-craft/blob/main/docs/rules-language.md)
//! defines the propensities, the expression syntax and which names each
//! field may use. In brief:
//!
//! - Names, including the model's `name`, are ASCII identifiers
//!   (`[A-Za-z_][A-Za-z0-9_]*`, at most 64 bytes) that are not function
//!   names. Parameters and observables share one namespace.
//! - Expressions use `|| && == != < <= > >= + - * / ^`, unary `-` and `!`,
//!   and `min max exp log sqrt abs floor ceil if`. Product attributes may
//!   also draw `bernoulli(p)`, `choice(p, a, b)`, `uniform(a, b)`,
//!   `normal(m, s)` and `new_id()`.
//! - A rule with no reactant fires at its rate g; with one pattern at
//!   g × Σ h(x) n(x); with two identical patterns at g times the sum over
//!   unordered pairs of distinct molecules; with patterns over two types at
//!   g × (Σ h n)(Σ h n). For one instance of `A` this is mass action with
//!   the combinatorial count of distinct reactant molecules: n for `A`,
//!   n(n − 1)/2 for `A + A` and n_A n_B for `A + B`.
//! - A rule with a `propensity` P (rules-model/0.3) fires at P while every
//!   reactant pattern can take a distinct molecule, and at 0 otherwise.
//!
//! The limits below are checked during validation and reported as
//! `CRN-R015` (`CRN-R010` for names, `CRN-R014` for counts).
//!
//! ```
//! use crn_craft::rules::model::{Model, Reduce};
//!
//! // Clusters merge pairwise; the observables count clusters and mass.
//! // `check_model` and `prepare` parse strictly; plain `serde_json` does not
//! // reject duplicate keys.
//! let model: Model = serde_json::from_str(r#"{
//!     "schema": "crn-craft.rules-model/0.1",
//!     "name": "aggregation",
//!     "parameters": {"k": 0.01},
//!     "types": [{"name": "C", "attributes": ["size"]}],
//!     "initial": [{"type": "C", "attributes": {"size": 1}, "count": 12}],
//!     "observables": {
//!         "clusters": {"type": "C"},
//!         "mass": {"type": "C", "weight": "x.size"},
//!         "largest": {"type": "C", "weight": "x.size", "reduce": "max"}
//!     },
//!     "rules": [{
//!         "name": "merge",
//!         "reactants": [{"label": "a", "type": "C"}, {"label": "b", "type": "C"}],
//!         "rate": "k",
//!         "products": [{"type": "C", "attributes": {"size": "a.size + b.size"}}]
//!     }],
//!     "stop": "clusters == 1"
//! }"#).unwrap();
//! assert_eq!(model.observables["largest"].reduce, Reduce::Max);
//! assert!(crn_craft::rules::check_model(&serde_json::to_vec(&model).unwrap()).valid);
//! ```

use crate::diagnostic::Diagnostic;
use crate::rules::expr::{Expr, Predicate, Scope};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Most species types in a model.
pub const MAX_TYPES: usize = 64;
/// Most attributes of one type.
pub const MAX_ATTRIBUTES: usize = 64;
/// Most parameters in a model.
pub const MAX_PARAMETERS: usize = 1024;
/// Most observables in a model; also the most trackers in a run request
/// (`CRN-R020` there).
pub const MAX_OBSERVABLES: usize = 256;
/// Most rules in a model.
pub const MAX_RULES: usize = 1024;
/// Most products of one rule.
pub const MAX_PRODUCTS: usize = 64;
/// Most entries in `initial` (before equal attribute vectors merge).
pub const MAX_INITIAL: usize = 100_000;
/// Longest expression in bytes, including tracker `first_passage`
/// predicates.
pub const MAX_EXPRESSION_BYTES: usize = 4096;
/// Longest name in bytes (`CRN-R010`).
pub const MAX_NAME_BYTES: usize = 64;
/// Counts and identities stay exact integers in f64.
///
/// Initial counts may be at most 2^53 (`CRN-R014`). Identities from
/// `new_id()` must stay below it: a model that calls `new_id()` with an
/// initial attribute of 2^53 - 1 or more is rejected (`CRN-R014`), and
/// reaching 2^53 during a run fails the replicate (`CRN-R030`).
pub const MAX_COUNT: u64 = 1 << 53;

/// The version of a rules model, its `schema` field.
///
/// Versions are ordered; each adds fields to the previous one.
#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq, PartialOrd, Ord,
)]
pub enum ModelVersion {
    // Variant docs use `concat!` to stay out of the schema's string enum;
    // see the note at the top of the crate root.
    #[doc = concat!("`crn-craft.rules-model/0.1`: the rules language without ATP.")]
    #[serde(rename = "crn-craft.rules-model/0.1")]
    V01,
    #[doc = concat!(
        "`crn-craft.rules-model/0.2`: 0.1 plus the optional `atp` block, ",
        "which requires the `thinning` method."
    )]
    #[serde(rename = "crn-craft.rules-model/0.2")]
    V02,
    #[doc = concat!(
        "`crn-craft.rules-model/0.3`: 0.2 plus whole-rule `propensity` ",
        "expressions in rules and `continuous` variables, which require the ",
        "`hybrid` method."
    )]
    #[serde(rename = "crn-craft.rules-model/0.3")]
    V03,
    #[doc = concat!(
        "`crn-craft.rules-model/0.4`: 0.3 plus `independent` continuous ",
        "variables, which `hybrid` integrates once per prepared run."
    )]
    #[serde(rename = "crn-craft.rules-model/0.4")]
    V04,
}

/// A rules model: `crn-craft.rules-model/0.1`, `0.2`, `0.3` or `0.4`.
///
/// Unknown fields and duplicate keys are rejected. Validation also requires
/// at least one type and one rule.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Model {
    // Docs of `$ref` fields use `concat!` so the schema keeps a bare `$ref`.
    #[doc = concat!("The model version, see [`ModelVersion`].")]
    pub schema: ModelVersion,
    /// The model's name, an ASCII identifier of at most 64 bytes.
    pub name: String,
    /// Free text for readers; not interpreted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// Finite numbers by name. Names share one namespace with observables.
    #[serde(default)]
    pub parameters: BTreeMap<String, f64>,
    /// Species types, 1 to 64, with unique names.
    pub types: Vec<SpeciesType>,
    /// The initial instances, at most 100,000 entries. Entries with equal
    /// attribute vectors merge; an empty list starts from no instances.
    pub initial: Vec<Initial>,
    /// Observables by name, at most 256. Results list them sorted by name.
    #[serde(default)]
    pub observables: BTreeMap<String, Observable>,
    /// Rules, 1 to 1024, with unique names. Results count events per rule in
    /// this order.
    pub rules: Vec<Rule>,
    /// Parameters and observables. Once true, the state is held and the
    /// replicate is marked `stopped`. It is checked at time 0 and after each
    /// event; a NaN value fails the replicate (`CRN-R030`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stop: Option<String>,
    /// ATP coupling, rules-model/0.2 only (`CRN-R015` in 0.1); a model with
    /// it runs only with the `thinning` method. `null` is rejected.
    #[serde(
        default,
        deserialize_with = "crate::rules::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "crate::rules::atp::AtpSpec")]
    pub atp: Option<crate::rules::atp::AtpSpec>,
    /// Continuous variables, rules-model/0.3 and later (`CRN-R015` before),
    /// at most 64: ODEs between events whose values rule rates and
    /// propensities may read. A model with them runs only with the `hybrid`
    /// method, and it cannot also have `atp` (`CRN-R015`). `null` is
    /// rejected. See [`crate::continuous`].
    // schemars copies this doc text into the rules-model schema, so it keeps
    // the 0.6.0 path, which rustdoc cannot resolve inside `crn_craft`.
    #[allow(rustdoc::broken_intra_doc_links)]
    #[serde(
        default,
        deserialize_with = "crate::rules::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "Vec<crate::rules::continuous::ContinuousSpec>")]
    pub continuous: Option<Vec<crate::rules::continuous::ContinuousSpec>>,
}

/// A species type: a name and its numeric attributes.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SpeciesType {
    /// Type name, unique in the model.
    pub name: String,
    /// Numeric attributes; an instance is a type plus one value per attribute.
    /// At most 64 unique names; values are finite and `-0` is stored as `0`.
    #[serde(default)]
    pub attributes: Vec<String>,
}

/// A number, or an expression string.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum Value {
    /// A literal number; must be finite.
    Number(f64),
    /// An expression, evaluated with the names the field allows.
    Expression(String),
}

/// One entry of the initial state: a type, every attribute value and a copy
/// number.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Initial {
    /// Name of a declared type.
    #[serde(rename = "type")]
    pub species_type: String,
    /// Every attribute of the type; parameters only.
    #[serde(default)]
    pub attributes: BTreeMap<String, Value>,
    /// Nonnegative integer; parameters only.
    pub count: Value,
}

/// How an observable reduces the matching instances to one number.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Reduce {
    /// Count-weighted sum over matching instances.
    #[default]
    Sum,
    /// Sum of weights over present matching instances, each counted once
    /// whatever its copy number (weight 1 counts distinct instances).
    Distinct,
    /// Largest weight over present matching instances; 0 when none.
    Max,
    /// Smallest weight over present matching instances; 0 when none.
    Min,
}

/// A named reduction over the present instances of one type that satisfy
/// `where`.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Observable {
    /// Name of a declared type.
    #[serde(rename = "type")]
    pub species_type: String,
    /// Name for the instance in `where` and `weight`; default `x`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    /// Predicate over parameters and the instance's attributes, for example
    /// `x.wt == 0`; default true. A NaN value is an error, not true.
    #[serde(rename = "where", default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    /// Parameters and the instance's attributes; default 1. Weights must be
    /// finite and may be negative.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub weight: Option<String>,
    /// `sum` (default), `distinct`, `max` or `min`.
    #[serde(default)]
    pub reduce: Reduce,
}

/// A rule schema: one rule stands for the reactions of every matching
/// instance.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    /// Rule name, unique in the model; event logs and `rule_events` refer to
    /// rules by their index in the model.
    pub name: String,
    /// Zero, one or two patterns. Reactants are consumed; list a `copy`
    /// product to keep one.
    pub reactants: Vec<Pattern>,
    /// Global factor: parameters and observables only. It is evaluated only
    /// while the reactant weight is positive; a negative or non-finite value
    /// fails the replicate (`CRN-R030`). Required in rules-model/0.1 and 0.2
    /// (`CRN-E003` when missing); from 0.3 a rule gives exactly one of
    /// `rate` and `propensity` (`CRN-R013`). `null` is rejected.
    #[serde(
        default,
        deserialize_with = "crate::rules::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "String")]
    pub rate: Option<String>,
    /// rules-model/0.3 and later (`CRN-R015` before): the rule's whole
    /// propensity, an expression over parameters and observables, instead
    /// of a `rate` per reactant tuple. It applies while every reactant type has enough
    /// molecules (one per pattern, two for two patterns over one type) and
    /// is 0 otherwise; a firing picks distinct reactant molecules uniformly.
    /// The reactant patterns take no `where` or `factor` (`CRN-R013`), and
    /// the expression cannot read their labels (`CRN-R012`). A negative or
    /// non-finite value fails the replicate (`CRN-R030`). `null` is
    /// rejected.
    #[serde(
        default,
        deserialize_with = "crate::rules::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "String")]
    pub propensity: Option<String>,
    /// Up to 64 products. A product equal to a reactant cancels its
    /// consumption.
    pub products: Vec<Product>,
}

/// A reactant pattern: the instances of one type that may react, and their
/// local factor.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Pattern {
    /// Name of the matched instance in this pattern's expressions and in the
    /// rule's products, for example `x` in `x.size`. The two patterns of a
    /// rule need different labels.
    pub label: String,
    /// Name of a declared type. Two different patterns over the same type
    /// are not supported (`CRN-R013`); identical ones form pairs.
    #[serde(rename = "type")]
    pub species_type: String,
    /// Parameters and this reactant's attributes.
    #[serde(rename = "where", default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    /// Local factor h(x) >= 0: parameters and this reactant's attributes.
    /// Default 1; an instance that fails `where` has h = 0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub factor: Option<String>,
}

/// A product of a rule: a new instance built from constants, parameters,
/// reactant attributes and random draws.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Product {
    /// Name of a declared type.
    #[serde(rename = "type")]
    pub species_type: String,
    /// Label of a reactant of the same type whose attributes are copied.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub copy: Option<String>,
    /// Overrides (with `copy`) or all attributes (without): parameters,
    /// reactant attributes and random draws.
    #[serde(default)]
    pub attributes: BTreeMap<String, Value>,
}

// ---------------------------------------------------------------------------
// Compiled form

#[derive(Debug)]
pub(crate) struct TypeInfo {
    pub name: String,
    pub attributes: Vec<String>,
    pub classes: Vec<usize>,
    pub observables: Vec<usize>,
}

#[derive(Debug)]
pub(crate) struct ClassInfo {
    pub type_id: usize,
    pub condition: Option<Expr>,
    pub factor: Option<Expr>,
    pub squares: bool,
    pub pointer: String,
}

#[derive(Debug)]
pub(crate) struct ObservableInfo {
    pub name: String,
    pub condition: Option<Expr>,
    pub weight: Option<Expr>,
    pub reduce: Reduce,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Shape {
    Source,
    Single(usize),
    /// Two identical patterns: unordered pairs of distinct molecules.
    Pair(usize),
    /// Two patterns over different types.
    Cross(usize, usize),
}

#[derive(Debug)]
pub(crate) enum Source {
    Copy { reactant: usize, index: usize },
    Expression(Expr),
}

#[derive(Debug)]
pub(crate) struct ProductInfo {
    pub type_id: usize,
    pub attributes: Vec<Source>,
}

#[derive(Debug)]
pub(crate) struct RuleInfo {
    pub name: String,
    /// The rule's `rate`, or its `propensity` when `whole`.
    pub rate: Expr,
    /// A whole-rule propensity: the rate applies once while the reactants are
    /// available instead of per reactant tuple.
    pub whole: bool,
    pub shape: Shape,
    /// Reactant type ids in label order.
    pub reactant_types: Vec<usize>,
    pub products: Vec<ProductInfo>,
}

impl RuleInfo {
    /// The model field that holds the rule's rate expression.
    pub fn field(&self) -> &'static str {
        if self.whole { "propensity" } else { "rate" }
    }

    /// JSON Pointer to the rate expression of rule `index`.
    pub fn pointer(&self, index: usize) -> String {
        format!("/rules/{index}/{}", self.field())
    }
}

#[derive(Debug)]
pub(crate) struct Compiled {
    pub params: Vec<f64>,
    pub types: Vec<TypeInfo>,
    /// (type, attributes, count), merged and without zero counts.
    pub initial: Vec<(usize, Vec<f64>, u64)>,
    pub first_id: f64,
    pub classes: Vec<ClassInfo>,
    pub observables: Vec<ObservableInfo>,
    pub rules: Vec<RuleInfo>,
    pub stop: Option<Expr>,
    pub stop_program: Option<Predicate>,
    pub atp: Option<crate::rules::atp::Atp>,
    /// Continuous variables; their values follow the model's observables in
    /// `observables`, in declaration order.
    pub continuous: Vec<crate::rules::continuous::Variable>,
}

const FUNCTIONS: [&str; 14] = [
    "min",
    "max",
    "exp",
    "log",
    "sqrt",
    "abs",
    "floor",
    "ceil",
    "if",
    "bernoulli",
    "choice",
    "uniform",
    "normal",
    "new_id",
];

pub(crate) fn invalid(
    code: &str,
    path: &str,
    message: impl Into<String>,
    hint: &str,
) -> Diagnostic {
    Diagnostic::error(code, path, message, hint)
}

fn check_name(name: &str, path: &str, out: &mut Vec<Diagnostic>) {
    let mut chars = name.chars();
    let valid = name.len() <= MAX_NAME_BYTES
        && chars
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
        && !FUNCTIONS.contains(&name);
    if !valid {
        out.push(invalid(
            "CRN-R010",
            path,
            format!("invalid name `{name}`"),
            "Use an ASCII identifier of at most 64 bytes ([A-Za-z_][A-Za-z0-9_]*) that is not a function name.",
        ));
    }
}

pub(crate) fn pointer(parent: &str, key: impl std::fmt::Display) -> String {
    crate::diagnostic::pointer(parent, key)
}

struct Names {
    params: BTreeMap<String, usize>,
    observables: BTreeMap<String, usize>,
}

pub(crate) fn parse_expr(
    source: &str,
    scope: &Scope,
    path: &str,
    out: &mut Vec<Diagnostic>,
) -> Option<Expr> {
    if source.len() > MAX_EXPRESSION_BYTES {
        out.push(invalid(
            "CRN-R015",
            path,
            "expression is longer than 4096 bytes",
            "Split the model or shorten the expression.",
        ));
        return None;
    }
    match Expr::parse(source, scope) {
        Ok(e) => Some(e),
        Err(message) => {
            out.push(invalid(
                "CRN-R012",
                path,
                format!("{message} in `{source}`"),
                "See docs/rules-language.md for the names and functions each field may use.",
            ));
            None
        }
    }
}

/// Deterministic expression over parameters only.
fn constant(
    value: &Value,
    names: &Names,
    params: &[f64],
    path: &str,
    out: &mut Vec<Diagnostic>,
) -> Option<f64> {
    match value {
        Value::Number(x) => Some(*x),
        Value::Expression(source) => {
            let scope = Scope {
                params: &names.params,
                observables: None,
                reactants: &[],
                random: false,
            };
            let e = parse_expr(source, &scope, path, out)?;
            Some(e.eval(&crate::rules::expr::Env {
                params,
                observables: &[],
                attributes: [&[], &[]],
            }))
        }
    }
}

fn normalize(x: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x }
}

pub(crate) fn compile(model: &Model) -> Result<Compiled, Vec<Diagnostic>> {
    let mut out = Vec::new();
    let limit = |out: &mut Vec<Diagnostic>, n: usize, max: usize, path: &str, what: &str| {
        if n > max {
            out.push(invalid(
                "CRN-R015",
                path,
                format!("{n} {what} exceed the limit of {max}"),
                "Reduce the model size.",
            ));
        }
    };
    limit(&mut out, model.types.len(), MAX_TYPES, "/types", "types");
    limit(
        &mut out,
        model.parameters.len(),
        MAX_PARAMETERS,
        "/parameters",
        "parameters",
    );
    limit(
        &mut out,
        model.observables.len(),
        MAX_OBSERVABLES,
        "/observables",
        "observables",
    );
    limit(&mut out, model.rules.len(), MAX_RULES, "/rules", "rules");
    limit(
        &mut out,
        model.initial.len(),
        MAX_INITIAL,
        "/initial",
        "initial entries",
    );
    if !out.is_empty() {
        return Err(out);
    }
    check_name(&model.name, "/name", &mut out);
    if model.types.is_empty() {
        out.push(invalid(
            "CRN-R011",
            "/types",
            "no species types",
            "Declare at least one type.",
        ));
    }

    // Names: parameters and observables share a namespace.
    let mut names = Names {
        params: BTreeMap::new(),
        observables: BTreeMap::new(),
    };
    let mut params = Vec::new();
    for (i, (name, value)) in model.parameters.iter().enumerate() {
        let path = pointer("/parameters", name);
        check_name(name, &path, &mut out);
        if !value.is_finite() {
            out.push(invalid(
                "CRN-R014",
                &path,
                "parameter is not finite",
                "Use a finite number.",
            ));
        }
        names.params.insert(name.clone(), i);
        params.push(normalize(*value));
    }
    for (i, name) in model.observables.keys().enumerate() {
        let path = pointer("/observables", name);
        check_name(name, &path, &mut out);
        if names.params.contains_key(name) {
            out.push(invalid(
                "CRN-R010",
                &path,
                format!("observable `{name}` shadows a parameter"),
                "Rename the observable or the parameter.",
            ));
        }
        names.observables.insert(name.clone(), i);
    }
    // Continuous variables share the namespace and follow the observables.
    let specs = model.continuous.as_deref().unwrap_or(&[]);
    if model.continuous.is_some() && model.schema < ModelVersion::V03 {
        out.push(invalid(
            "CRN-R015",
            "/continuous",
            "continuous variables require rules-model/0.3",
            "Select the 0.3 model contract.",
        ));
    }
    limit(
        &mut out,
        specs.len(),
        crate::rules::continuous::MAX_VARIABLES,
        "/continuous",
        "continuous variables",
    );
    for (k, spec) in specs.iter().enumerate() {
        let path = pointer(&pointer("/continuous", k), "name");
        check_name(&spec.name, &path, &mut out);
        if names.params.contains_key(&spec.name) || names.observables.contains_key(&spec.name) {
            out.push(invalid(
                "CRN-R010",
                &path,
                format!(
                    "continuous variable `{}` repeats a parameter, observable or variable name",
                    spec.name
                ),
                "Give every parameter, observable and continuous variable its own name.",
            ));
        }
        names
            .observables
            .insert(spec.name.clone(), model.observables.len() + k);
    }

    let mut types: Vec<TypeInfo> = Vec::new();
    let mut type_ids = BTreeMap::new();
    for (t, ty) in model.types.iter().enumerate() {
        let path = pointer("/types", t);
        check_name(&ty.name, &pointer(&path, "name"), &mut out);
        if type_ids.insert(ty.name.clone(), t).is_some() {
            out.push(invalid(
                "CRN-R010",
                &pointer(&path, "name"),
                format!("duplicate type `{}`", ty.name),
                "Type names must be unique.",
            ));
        }
        limit(
            &mut out,
            ty.attributes.len(),
            MAX_ATTRIBUTES,
            &pointer(&path, "attributes"),
            "attributes",
        );
        for (a, attribute) in ty.attributes.iter().enumerate() {
            let p = pointer(&pointer(&path, "attributes"), a);
            check_name(attribute, &p, &mut out);
            if ty.attributes[..a].contains(attribute) {
                out.push(invalid(
                    "CRN-R010",
                    &p,
                    format!("duplicate attribute `{attribute}`"),
                    "Attribute names must be unique within a type.",
                ));
            }
        }
        types.push(TypeInfo {
            name: ty.name.clone(),
            attributes: ty.attributes.clone(),
            classes: Vec::new(),
            observables: Vec::new(),
        });
    }
    let type_of = |name: &str, path: &str, out: &mut Vec<Diagnostic>| -> Option<usize> {
        let id = type_ids.get(name).copied();
        if id.is_none() {
            out.push(invalid(
                "CRN-R011",
                path,
                format!("unknown type `{name}`"),
                "Declare the type under /types.",
            ));
        }
        id
    };
    if !out.is_empty() {
        return Err(out);
    }

    // Initial state.
    let mut merged: BTreeMap<(usize, Vec<u64>), (Vec<f64>, u64)> = BTreeMap::new();
    let mut largest = 0.0f64;
    for (i, entry) in model.initial.iter().enumerate() {
        let path = pointer("/initial", i);
        let Some(t) = type_of(&entry.species_type, &pointer(&path, "type"), &mut out) else {
            continue;
        };
        let attributes = &types[t].attributes;
        for key in entry.attributes.keys() {
            if !attributes.contains(key) {
                out.push(invalid(
                    "CRN-R013",
                    &pointer(&pointer(&path, "attributes"), key),
                    format!("type `{}` has no attribute `{key}`", types[t].name),
                    "Remove the attribute or declare it on the type.",
                ));
            }
        }
        let mut values = Vec::with_capacity(attributes.len());
        for name in attributes {
            let p = pointer(&pointer(&path, "attributes"), name);
            let Some(value) = entry.attributes.get(name) else {
                out.push(invalid(
                    "CRN-R013",
                    &p,
                    format!("missing attribute `{name}`"),
                    "Give every attribute of the type.",
                ));
                continue;
            };
            if let Some(x) = constant(value, &names, &params, &p, &mut out) {
                if x.is_finite() {
                    largest = largest.max(x);
                    values.push(normalize(x));
                } else {
                    out.push(invalid(
                        "CRN-R014",
                        &p,
                        "attribute is not finite",
                        "Use a finite value.",
                    ));
                }
            }
        }
        let p = pointer(&path, "count");
        let count = constant(&entry.count, &names, &params, &p, &mut out);
        let Some(count) = count else { continue };
        #[allow(clippy::cast_precision_loss)]
        if !(count >= 0.0 && count.fract() == 0.0 && count <= MAX_COUNT as f64) {
            out.push(invalid(
                "CRN-R014",
                &p,
                format!("count {count} is not an integer in [0, 2^53]"),
                "Use a nonnegative integer count.",
            ));
            continue;
        }
        if values.len() != attributes.len() {
            continue;
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = count as u64;
        if count == 0 {
            continue;
        }
        let key = (t, values.iter().map(|x| x.to_bits()).collect());
        let entry = merged.entry(key).or_insert((values, 0));
        entry.1 += count;
        if entry.1 > MAX_COUNT {
            out.push(invalid(
                "CRN-R014",
                &p,
                "merged initial count exceeds 2^53",
                "Reduce the initial counts.",
            ));
        }
    }
    let initial: Vec<_> = merged
        .into_iter()
        .map(|((t, _), (v, n))| (t, v, n))
        .collect();
    // new_id() values exceed every initial attribute value (and 0).
    let first_id = libm::floor(largest) + 1.0;

    // Observables.
    let mut observables = Vec::new();
    for (name, obs) in &model.observables {
        let path = pointer("/observables", name);
        let Some(t) = type_of(&obs.species_type, &pointer(&path, "type"), &mut out) else {
            continue;
        };
        let label = obs.label.clone().unwrap_or_else(|| "x".into());
        check_name(&label, &pointer(&path, "label"), &mut out);
        let reactants = [(label.as_str(), types[t].attributes.as_slice())];
        let scope = Scope {
            params: &names.params,
            observables: None,
            reactants: &reactants,
            random: false,
        };
        let condition = obs
            .condition
            .as_ref()
            .and_then(|s| parse_expr(s, &scope, &pointer(&path, "where"), &mut out));
        let weight = obs
            .weight
            .as_ref()
            .and_then(|s| parse_expr(s, &scope, &pointer(&path, "weight"), &mut out));
        types[t].observables.push(observables.len());
        observables.push(ObservableInfo {
            name: name.clone(),
            condition,
            weight,
            reduce: obs.reduce,
        });
    }

    let global_scope = Scope {
        params: &names.params,
        observables: Some(&names.observables),
        reactants: &[],
        random: false,
    };
    let stop = model
        .stop
        .as_ref()
        .and_then(|s| parse_expr(s, &global_scope, "/stop", &mut out));
    if let Some(stop) = &stop {
        // A stop is checked at events; continuous values change in between.
        let mut flags = vec![false; model.observables.len() + specs.len()];
        stop.reads(&mut flags);
        if let Some(k) = flags[model.observables.len()..].iter().position(|&f| f) {
            out.push(invalid(
                "CRN-R012",
                "/stop",
                format!("stop reads continuous variable `{}`", specs[k].name),
                "Stop on a continuous variable with its lower_bound and on_bound \"stop\".",
            ));
        }
    }

    // Rules.
    let mut classes: Vec<ClassInfo> = Vec::new();
    let mut rules = Vec::new();
    let mut rule_names = BTreeMap::new();
    for (r, rule) in model.rules.iter().enumerate() {
        let path = pointer("/rules", r);
        check_name(&rule.name, &pointer(&path, "name"), &mut out);
        if rule_names.insert(rule.name.clone(), r).is_some() {
            out.push(invalid(
                "CRN-R010",
                &pointer(&path, "name"),
                format!("duplicate rule `{}`", rule.name),
                "Rule names must be unique.",
            ));
        }
        let whole = rule.propensity.is_some();
        let rate = match (&rule.rate, &rule.propensity) {
            (Some(rate), None) => {
                parse_expr(rate, &global_scope, &pointer(&path, "rate"), &mut out)
            }
            (None, Some(_)) | (Some(_), Some(_)) if model.schema < ModelVersion::V03 => {
                out.push(invalid(
                    "CRN-R015",
                    &pointer(&path, "propensity"),
                    "propensity requires rules-model/0.3",
                    "Select the 0.3 model contract.",
                ));
                None
            }
            (None, Some(propensity)) => parse_expr(
                propensity,
                &global_scope,
                &pointer(&path, "propensity"),
                &mut out,
            ),
            (Some(_), Some(_)) => {
                out.push(invalid(
                    "CRN-R013",
                    &path,
                    "a rule gives both rate and propensity",
                    "Give `rate` per reactant tuple or the whole `propensity`, not both.",
                ));
                None
            }
            (None, None) if model.schema < ModelVersion::V03 => {
                // A required field before 0.3: report it as the parser did.
                out.push(invalid(
                    "CRN-E003",
                    &path,
                    "missing field `rate`",
                    "Compare this field with the schema emitted by describe.",
                ));
                None
            }
            (None, None) => {
                out.push(invalid(
                    "CRN-R013",
                    &path,
                    "a rule needs a rate or a propensity",
                    "Give `rate` per reactant tuple or the whole `propensity`.",
                ));
                None
            }
        };
        if whole {
            for (i, p) in rule.reactants.iter().enumerate() {
                for (key, present) in [
                    ("where", p.condition.is_some()),
                    ("factor", p.factor.is_some()),
                ] {
                    if present {
                        out.push(invalid(
                            "CRN-R013",
                            &pointer(&pointer(&pointer(&path, "reactants"), i), key),
                            format!("a propensity rule's reactants take no `{key}`"),
                            "Write the dependence into the propensity, or use `rate` with reactant patterns.",
                        ));
                    }
                }
            }
        }
        if rule.reactants.len() > 2 {
            out.push(invalid(
                "CRN-R013",
                &pointer(&path, "reactants"),
                "more than two reactant patterns",
                "Split the rule; rules-model/0.1 supports zero to two reactants.",
            ));
            continue;
        }
        if rule.products.len() > MAX_PRODUCTS {
            out.push(invalid(
                "CRN-R015",
                &pointer(&path, "products"),
                "more than 64 products",
                "Reduce the number of products.",
            ));
            continue;
        }
        // Patterns.
        let mut pattern_types = Vec::new();
        let mut compiled = Vec::new();
        for (i, p) in rule.reactants.iter().enumerate() {
            let pp = pointer(&pointer(&path, "reactants"), i);
            check_name(&p.label, &pointer(&pp, "label"), &mut out);
            let Some(t) = type_of(&p.species_type, &pointer(&pp, "type"), &mut out) else {
                continue;
            };
            let reactants = [(p.label.as_str(), types[t].attributes.as_slice())];
            let scope = Scope {
                params: &names.params,
                observables: None,
                reactants: &reactants,
                random: false,
            };
            let condition = p
                .condition
                .as_ref()
                .and_then(|s| parse_expr(s, &scope, &pointer(&pp, "where"), &mut out));
            let factor = p
                .factor
                .as_ref()
                .and_then(|s| parse_expr(s, &scope, &pointer(&pp, "factor"), &mut out));
            pattern_types.push(t);
            compiled.push((t, condition, factor, pp));
        }
        if compiled.len() != rule.reactants.len() {
            continue;
        }
        if rule.reactants.len() == 2 && rule.reactants[0].label == rule.reactants[1].label {
            out.push(invalid(
                "CRN-R010",
                &format!("{path}/reactants/1/label"),
                "duplicate reactant label",
                "Give the two reactants different labels.",
            ));
            continue;
        }
        let same = |a: &(usize, Option<Expr>, Option<Expr>, String),
                    b: &(usize, Option<Expr>, Option<Expr>, String)| {
            a.0 == b.0
                && Expr::same_option(a.1.as_ref(), b.1.as_ref())
                && Expr::same_option(a.2.as_ref(), b.2.as_ref())
        };
        let mut class_of = |c: (usize, Option<Expr>, Option<Expr>, String), squares: bool| {
            let key = (c.0, c.1.clone(), c.2.clone(), String::new());
            if let Some(i) = classes.iter().position(|k| {
                same(
                    &(
                        k.type_id,
                        k.condition.clone(),
                        k.factor.clone(),
                        String::new(),
                    ),
                    &key,
                )
            }) {
                classes[i].squares |= squares;
                return i;
            }
            types[c.0].classes.push(classes.len());
            classes.push(ClassInfo {
                type_id: c.0,
                condition: c.1,
                factor: c.2,
                squares,
                pointer: c.3,
            });
            classes.len() - 1
        };
        let shape = match compiled.len() {
            0 => Shape::Source,
            1 => Shape::Single(class_of(compiled.pop().expect("one pattern"), false)),
            _ => {
                let second = compiled.pop().expect("two patterns");
                let first = compiled.pop().expect("two patterns");
                if same(&first, &second) {
                    Shape::Pair(class_of(first, true))
                } else if first.0 == second.0 {
                    out.push(invalid(
                        "CRN-R013",
                        &pointer(&path, "reactants"),
                        "two different patterns over the same type",
                        "Use identical patterns (same where and factor) or patterns over different types; rules-model/0.1 has no other same-type pairs.",
                    ));
                    continue;
                } else {
                    Shape::Cross(class_of(first, false), class_of(second, false))
                }
            }
        };
        // Products.
        let labels: Vec<(&str, &[String])> = rule
            .reactants
            .iter()
            .zip(&pattern_types)
            .map(|(p, &t)| (p.label.as_str(), types[t].attributes.as_slice()))
            .collect();
        let scope = Scope {
            params: &names.params,
            observables: None,
            reactants: &labels,
            random: true,
        };
        let mut products = Vec::new();
        for (j, product) in rule.products.iter().enumerate() {
            let pp = pointer(&pointer(&path, "products"), j);
            let Some(t) = type_of(&product.species_type, &pointer(&pp, "type"), &mut out) else {
                continue;
            };
            let copy = match &product.copy {
                None => None,
                Some(label) => {
                    let found = rule.reactants.iter().position(|p| &p.label == label);
                    match found {
                        Some(k) if pattern_types[k] == t => Some(k),
                        Some(_) => {
                            out.push(invalid(
                                "CRN-R013",
                                &pointer(&pp, "copy"),
                                format!("reactant `{label}` has a different type"),
                                "Copy only a reactant of the product's type.",
                            ));
                            continue;
                        }
                        None => {
                            out.push(invalid(
                                "CRN-R013",
                                &pointer(&pp, "copy"),
                                format!("unknown reactant label `{label}`"),
                                "Copy one of the rule's reactant labels.",
                            ));
                            continue;
                        }
                    }
                }
            };
            for key in product.attributes.keys() {
                if !types[t].attributes.contains(key) {
                    out.push(invalid(
                        "CRN-R013",
                        &pointer(&pointer(&pp, "attributes"), key),
                        format!("type `{}` has no attribute `{key}`", types[t].name),
                        "Remove the attribute or declare it on the type.",
                    ));
                }
            }
            let mut sources = Vec::new();
            for (index, name) in types[t].attributes.iter().enumerate() {
                let ap = pointer(&pointer(&pp, "attributes"), name);
                match (product.attributes.get(name), copy) {
                    (Some(Value::Number(x)), _) => {
                        if x.is_finite() {
                            sources.push(Source::Expression(Expr::constant(normalize(*x))));
                        } else {
                            out.push(invalid(
                                "CRN-R014",
                                &ap,
                                "attribute is not finite",
                                "Use a finite value.",
                            ));
                        }
                    }
                    (Some(Value::Expression(s)), _) => {
                        if let Some(e) = parse_expr(s, &scope, &ap, &mut out) {
                            sources.push(Source::Expression(e));
                        }
                    }
                    (None, Some(reactant)) => sources.push(Source::Copy { reactant, index }),
                    (None, None) => out.push(invalid(
                        "CRN-R013",
                        &ap,
                        format!("missing attribute `{name}`"),
                        "Give every attribute, or copy a reactant and override some.",
                    )),
                }
            }
            products.push(ProductInfo {
                type_id: t,
                attributes: sources,
            });
        }
        if let Some(rate) = rate {
            rules.push(RuleInfo {
                name: rule.name.clone(),
                rate,
                whole,
                shape,
                reactant_types: pattern_types,
                products,
            });
        }
    }
    if model.rules.is_empty() {
        out.push(invalid(
            "CRN-R013",
            "/rules",
            "no rules",
            "Declare at least one rule.",
        ));
    }
    let uses_new_id = rules.iter().any(|r| {
        r.products.iter().any(|p| {
            p.attributes
                .iter()
                .any(|a| matches!(a, Source::Expression(e) if e.uses_new_id()))
        })
    });
    if uses_new_id && first_id >= crate::rules::expr::ID_LIMIT {
        out.push(invalid(
            "CRN-R014",
            "/initial",
            format!(
                "initial attribute values up to {largest} leave no exact identities for new_id()"
            ),
            "Keep initial attribute values below 2^53 - 1 when rules call new_id().",
        ));
    }
    let atp = model.atp.as_ref().and_then(|spec| {
        if model.schema == ModelVersion::V01 {
            out.push(invalid("CRN-R015", "/atp", "ATP requires rules-model/0.2", "Select the 0.2 model contract."));
            return None;
        }
        match crate::rules::atp::Atp::compile(spec, &observables, &rules) {
            Ok(atp) => {
                observables.push(ObservableInfo { name: "ATP".into(), condition: None, weight: None, reduce: Reduce::Sum });
                Some(atp)
            }
            Err(message) => {
                out.push(invalid("CRN-R015", "/atp", message, "Use declared observable and replication-rule names with finite ATP parameters."));
                None
            }
        }
    });
    if model.atp.is_some() && model.continuous.is_some() {
        out.push(invalid(
            "CRN-R015",
            "/continuous",
            "a model has either an atp block or continuous variables",
            "Write ATP as a continuous variable, or keep the legacy atp block alone.",
        ));
    }
    let continuous = crate::rules::continuous::Variable::compile_all(
        specs,
        model.schema,
        &names.params,
        &names.observables,
        &params,
        model.observables.len(),
        &mut out,
    );
    for variable in &continuous {
        observables.push(ObservableInfo {
            name: variable.name.clone(),
            condition: None,
            weight: None,
            reduce: Reduce::Sum,
        });
    }
    if out.is_empty() {
        let stop_program = stop
            .as_ref()
            .map(|e| Predicate::compile(e, observables.len()));
        Ok(Compiled {
            params,
            types,
            initial,
            first_id,
            classes,
            observables,
            rules,
            stop_program,
            stop,
            atp,
            continuous,
        })
    } else {
        Err(out)
    }
}
