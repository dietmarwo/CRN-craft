//! `crn-craft.rules-model/0.1`: species types with numeric attributes, rule
//! schemata, observables and an optional stop condition.

use crate::expr::{Expr, Scope};
use crn_craft::diagnostic::Diagnostic;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MAX_TYPES: usize = 64;
pub const MAX_ATTRIBUTES: usize = 64;
pub const MAX_PARAMETERS: usize = 1024;
pub const MAX_OBSERVABLES: usize = 256;
pub const MAX_RULES: usize = 1024;
pub const MAX_PRODUCTS: usize = 64;
pub const MAX_INITIAL: usize = 100_000;
pub const MAX_EXPRESSION_BYTES: usize = 4096;
pub const MAX_NAME_BYTES: usize = 64;
/// Counts and identities stay exact integers in f64.
pub const MAX_COUNT: u64 = 1 << 53;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ModelVersion {
    #[serde(rename = "crn-craft.rules-model/0.1")]
    V01,
    #[serde(rename = "crn-craft.rules-model/0.2")]
    V02,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Model {
    pub schema: ModelVersion,
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// Finite numbers by name. Names share one namespace with observables.
    #[serde(default)]
    pub parameters: BTreeMap<String, f64>,
    pub types: Vec<SpeciesType>,
    pub initial: Vec<Initial>,
    #[serde(default)]
    pub observables: BTreeMap<String, Observable>,
    pub rules: Vec<Rule>,
    /// Parameters and observables. Once true, the state is held and the
    /// replicate is marked `stopped`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stop: Option<String>,
    /// rules-model/0.2 only. `null` is rejected.
    #[serde(
        default,
        deserialize_with = "crate::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "crate::atp::AtpSpec")]
    pub atp: Option<crate::atp::AtpSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SpeciesType {
    pub name: String,
    /// Numeric attributes; an instance is a type plus one value per attribute.
    #[serde(default)]
    pub attributes: Vec<String>,
}

/// A number, or an expression string.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum Value {
    Number(f64),
    Expression(String),
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Initial {
    #[serde(rename = "type")]
    pub species_type: String,
    /// Every attribute of the type; parameters only.
    #[serde(default)]
    pub attributes: BTreeMap<String, Value>,
    /// Nonnegative integer; parameters only.
    pub count: Value,
}

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

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Observable {
    #[serde(rename = "type")]
    pub species_type: String,
    /// Name for the instance in `where` and `weight`; default `x`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(rename = "where", default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    /// Parameters and the instance's attributes; default 1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub weight: Option<String>,
    #[serde(default)]
    pub reduce: Reduce,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    pub name: String,
    /// Zero, one or two patterns. Reactants are consumed; list a `copy`
    /// product to keep one.
    pub reactants: Vec<Pattern>,
    /// Global factor: parameters and observables only.
    pub rate: String,
    pub products: Vec<Product>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Pattern {
    pub label: String,
    #[serde(rename = "type")]
    pub species_type: String,
    /// Parameters and this reactant's attributes.
    #[serde(rename = "where", default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    /// Local factor h(x) >= 0: parameters and this reactant's attributes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub factor: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Product {
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
    pub rate: Expr,
    pub shape: Shape,
    /// Reactant type ids in label order.
    pub reactant_types: Vec<usize>,
    pub products: Vec<ProductInfo>,
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
    pub atp: Option<crate::atp::Atp>,
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

fn invalid(code: &str, path: &str, message: impl Into<String>, hint: &str) -> Diagnostic {
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

fn pointer(parent: &str, key: impl std::fmt::Display) -> String {
    crn_craft::diagnostic::pointer(parent, key)
}

struct Names {
    params: BTreeMap<String, usize>,
    observables: BTreeMap<String, usize>,
}

fn parse_expr(source: &str, scope: &Scope, path: &str, out: &mut Vec<Diagnostic>) -> Option<Expr> {
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
            Some(e.eval(&crate::expr::Env {
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
        let rate = parse_expr(&rule.rate, &global_scope, &pointer(&path, "rate"), &mut out);
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
    if uses_new_id && first_id >= crate::expr::ID_LIMIT {
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
        match crate::atp::Atp::compile(spec, &observables, &rules) {
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
    if out.is_empty() {
        Ok(Compiled {
            params,
            types,
            initial,
            first_id,
            classes,
            observables,
            rules,
            stop,
            atp,
        })
    } else {
        Err(out)
    }
}
