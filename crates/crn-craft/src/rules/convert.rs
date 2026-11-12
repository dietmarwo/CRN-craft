//! Conversion of network models (`crn-craft.model/0.1`) into rules models.
//!
//! [`from_network`] writes a `crn-craft.model/1.0` document with one
//! attribute-free species type per network species and one rule per
//! reaction, so the converted model has the same species, reactions and
//! propensities; run it with `crn-craft.run/1.0` requests:
//!
//! | Network model | Rules model |
//! | --- | --- |
//! | `name`, `parameters` | the same names and values |
//! | species `X` with `initial_count` n | type `X` without attributes, an `initial` entry with count n (none when n is 0) and the observable `X: {"type": "X"}` |
//! | reactants `{"A": 1}`, `{"A": 1, "B": 1}`, `{"A": 2}` | the patterns `x: A`; `x: A, y: B`; `x: A, y: A` (a pair) |
//! | three or more reactant molecules, for example `{"A": 2, "B": 1}` or `{"A": 1, "B": 1, "C": 1}` | one pattern per species with its stoichiometry as `count` (omitted when 1): `x: A count 2, y: B`; three or more species are labeled `x1`, `x2`, … |
//! | products `{"C": 2}` | the product `C` listed twice |
//! | `mass_action` with coefficient c | `"rate": "c"`, per reactant combination: c, c·n, c·n_A·n_B, c·n(n − 1)/2 or c·Π C(nᵢ, kᵢ) |
//! | `hill_activation` with regulator R | `"propensity": "basal + vmax / (1 + (k / R)^n)"` |
//! | `hill_repression` with regulator R | `"propensity": "basal + vmax / (1 + (R / k)^n)"` |
//! | `michaelis_menten` with substrate S | `"propensity": "vmax * (S / (km + S))"` |
//!
//! A rate argument that names a parameter is written as that name, a
//! literal as the shortest decimal that parses back to the same `f64`.
//! Hill and Michaelis–Menten laws are whole propensities that apply while
//! every reactant has its stoichiometric count, which is what a rules
//! `propensity` means. The Hill forms are the algebraically equal fractions
//! that cannot overflow: R = 0 gives `k / 0 = inf` and an induced term of 0,
//! as in the network engine. `omega` does not scale any propensity, so it is
//! not kept; the description records a value other than 1.
//!
//! Names follow the same identifier rules in both languages, except that the
//! rules language reserves its function names (`min`, `exp`, `if`, ...). Such
//! a name gets the suffix `_`, or `_2`, `_3`, ... if that is taken; the
//! converted model's `description` lists every renaming.
//!
//! The output is canonical: pretty-printed JSON with a trailing newline, keys
//! in the order of the rules contract types, parameters and observables
//! sorted by name. Models the rules language cannot hold yet are refused
//! with `CRN-R016`: no reaction, more than 64 product molecules in a
//! reaction, an initial count above 2^53 or a population limit. The size limits of `crn-craft.model/1.0` (4,096 types,
//! 16,384 parameters, 65,536 rules) equal those of the network contract.
//!
//! ```
//! use crn_craft::rules::{check_model, convert::from_network};
//!
//! let network = br#"{"schema": "crn-craft.model/0.1", "name": "decay",
//!   "convention": "stochastic_counts", "omega": 1.0, "parameters": {"k": 0.5},
//!   "species": [{"name": "A", "initial_count": 100}],
//!   "reactions": [{"name": "decay", "reactants": {"A": 1}, "products": {},
//!                  "rate": {"kind": "mass_action", "coefficient": "k"}}]}"#;
//! let rules = from_network(network).unwrap();
//! assert!(check_model(&rules).valid);
//! let model: serde_json::Value = serde_json::from_slice(&rules).unwrap();
//! assert_eq!(model["rules"][0]["rate"], "k");
//! assert_eq!(model["observables"]["A"]["type"], "A");
//! ```
use crate::contracts::{self, Rate, Scalar};
use crate::diagnostic::{Diagnostic, Severity, pointer};
use crate::rules::model::{
    FUNCTIONS, Initial, MAX_COUNT, MAX_PARAMETERS_1_0, MAX_PRODUCTS, MAX_RULES_1_0, MAX_TYPES_1_0,
    Model, ModelVersion, Observable, Pattern, Product, Reduce, Rule, SpeciesType, Value,
};
use std::collections::{BTreeMap, BTreeSet};

/// The contract version of a converted model.
pub const CONVERTED_VERSION: ModelVersion = ModelVersion::V10;

/// The `schema` string of the network models the converter reads.
pub(crate) const NETWORK_SCHEMA: &str = "crn-craft.model/0.1";

/// The `schema` string of a converted model.
pub(crate) fn converted_schema() -> serde_json::Value {
    serde_json::to_value(CONVERTED_VERSION).expect("a version serializes")
}

fn refusal(path: &str, message: impl Into<String>, hint: &str) -> Diagnostic {
    Diagnostic::error("CRN-R016", path, message, hint)
}

/// Converts a `crn-craft.model/0.1` document into a rules model document.
///
/// The input is parsed strictly and validated like `crn-craft check`; the
/// output passes [`crate::rules::check_model`]. See the [module
/// documentation](self) for the mapping and the canonical form.
///
/// # Errors
///
/// The parse and validation errors of the network model (`CRN-E001` to
/// `CRN-E016`; warnings do not block conversion), or `CRN-R016` at the
/// network document's path for each feature the rules language cannot hold.
pub fn from_network(model_bytes: &[u8]) -> Result<Vec<u8>, Vec<Diagnostic>> {
    let network = crate::parse_model(model_bytes).map_err(|d| vec![d])?;
    let errors: Vec<_> = crate::validate::validate_model(&network)
        .into_iter()
        .filter(|d| d.severity == Severity::Error)
        .collect();
    if !errors.is_empty() {
        return Err(errors);
    }
    let refusals = limits(&network);
    if !refusals.is_empty() {
        return Err(refusals);
    }
    let model = convert(&network);
    // The mapping is total on the checked domain; this guards it.
    crate::rules::model::compile(&model)?;
    let mut bytes = serde_json::to_vec_pretty(&model).expect("finite rules model");
    bytes.push(b'\n');
    Ok(bytes)
}

/// Reads a network model file and converts it like [`from_network`]; this is
/// the `crn-craft convert MODEL` command.
///
/// # Errors
///
/// `CRN-E030` for a file that cannot be read or is not a regular file,
/// `CRN-E004` for a file above 16 MiB, and otherwise the errors of
/// [`from_network`].
pub fn from_network_file(path: &std::path::Path) -> Result<Vec<u8>, Vec<Diagnostic>> {
    from_network(&crate::rules::read(path).map_err(|d| vec![d])?)
}

fn limits(network: &contracts::Model) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    // The model/1.0 limits equal the network contract's own (4,096 species,
    // 16,384 parameters, 65,536 reactions), so a valid network model passes
    // these checks; they guard the mapping if either limit changes.
    if network.species.len() > MAX_TYPES_1_0 {
        out.push(refusal(
            "/species",
            format!(
                "{} species exceed the rules limit of {MAX_TYPES_1_0} types",
                network.species.len()
            ),
            "Run this model with the network commands, or reduce it to at most 4096 species.",
        ));
    }
    if network.parameters.len() > MAX_PARAMETERS_1_0 {
        out.push(refusal(
            "/parameters",
            format!(
                "{} parameters exceed the rules limit of {MAX_PARAMETERS_1_0}",
                network.parameters.len()
            ),
            "Run this model with the network commands, or reduce it to at most 16384 parameters.",
        ));
    }
    if network.reactions.len() > MAX_RULES_1_0 {
        out.push(refusal(
            "/reactions",
            format!(
                "{} reactions exceed the rules limit of {MAX_RULES_1_0} rules",
                network.reactions.len()
            ),
            "Run this model with the network commands, or reduce it to at most 65536 reactions.",
        ));
    }
    if network.reactions.is_empty() {
        out.push(refusal(
            "/reactions",
            "a rules model needs at least one rule",
            "Add a reaction; a model without reactions never changes its state.",
        ));
    }
    for (i, species) in network.species.iter().enumerate() {
        let base = pointer("/species", i);
        if u64::try_from(species.initial_count).is_ok_and(|n| n > MAX_COUNT) {
            out.push(refusal(
                &pointer(&base, "initial_count"),
                "initial count above 2^53, the largest exact rules count",
                "Run this model with the network commands, or rescale the counts.",
            ));
        }
        if species.population_limit.is_some() {
            out.push(refusal(
                &pointer(&base, "population_limit"),
                "population limits have no rules equivalent",
                "Remove the limit, or run this model with the network commands; a rules stop condition ends a replicate after the event instead of refusing it.",
            ));
        }
    }
    for (i, reaction) in network.reactions.iter().enumerate() {
        let base = pointer("/reactions", i);
        let products: i64 = reaction.products.values().sum();
        if products > MAX_PRODUCTS as i64 {
            out.push(refusal(
                &pointer(&base, "products"),
                format!("{products} product molecules exceed the rules limit of {MAX_PRODUCTS}"),
                "Run this model with the network commands, or split the reaction.",
            ));
        }
    }
    out
}

/// Deterministic renaming of names the rules language reserves.
struct Names {
    species: BTreeMap<String, String>,
    parameters: BTreeMap<String, String>,
    reactions: BTreeMap<String, String>,
    model: String,
    /// `(kind, old, new)` in declaration order: model, species, parameters
    /// (sorted), reactions.
    renamed: Vec<(&'static str, String, String)>,
}

impl Names {
    fn new(network: &contracts::Model) -> Self {
        // Species and parameters share the rules namespace of observables
        // and parameters; reactions and the model name have their own.
        let mut shared: BTreeSet<String> = network
            .species
            .iter()
            .map(|s| s.name.clone())
            .chain(network.parameters.keys().cloned())
            .collect();
        let mut reactions: BTreeSet<String> =
            network.reactions.iter().map(|r| r.name.clone()).collect();
        let mut renamed = Vec::new();
        let mut rename = |kind: &'static str, name: &str, taken: &mut BTreeSet<String>| {
            if !FUNCTIONS.contains(&name) {
                return name.to_string();
            }
            let mut candidate = format!("{name}_");
            let mut suffix = 1;
            while taken.contains(&candidate) {
                suffix += 1;
                candidate = format!("{name}_{suffix}");
            }
            taken.insert(candidate.clone());
            renamed.push((kind, name.to_string(), candidate.clone()));
            candidate
        };
        let model = rename("model", &network.name, &mut BTreeSet::new());
        let species = network
            .species
            .iter()
            .map(|s| (s.name.clone(), rename("species", &s.name, &mut shared)))
            .collect();
        let parameters = network
            .parameters
            .keys()
            .map(|p| (p.clone(), rename("parameter", p, &mut shared)))
            .collect();
        let reactions = network
            .reactions
            .iter()
            .map(|r| (r.name.clone(), rename("reaction", &r.name, &mut reactions)))
            .collect();
        Self {
            species,
            parameters,
            reactions,
            model,
            renamed,
        }
    }

    fn scalar(&self, value: &Scalar) -> String {
        match value {
            Scalar::Parameter(name) => self.parameters[name].clone(),
            Scalar::Literal(x) => number(*x),
        }
    }
}

/// The shortest decimal that parses back to `x`, in a form the rules
/// expression language reads (validated literals are finite and
/// non-negative; `-0` is written as `0.0`).
fn number(x: f64) -> String {
    if x == 0.0 {
        "0.0".into()
    } else {
        format!("{x:?}")
    }
}

fn description(network: &contracts::Model, names: &Names) -> String {
    let mut text = String::from("Converted from a crn-craft.model/0.1 document.");
    if network.omega.to_bits() != 1.0f64.to_bits() {
        text.push_str(&format!(
            " Its omega, {:?}, does not scale propensities and is not kept.",
            network.omega
        ));
    }
    if !names.renamed.is_empty() {
        let list: Vec<String> = names
            .renamed
            .iter()
            .map(|(kind, old, new)| format!("{kind} `{old}` is `{new}`"))
            .collect();
        text.push_str(&format!(
            " Renamed because the rules language reserves function names: {}.",
            list.join("; ")
        ));
    }
    text
}

fn convert(network: &contracts::Model) -> Model {
    let names = Names::new(network);
    let species = |name: &String| names.species[name].clone();
    let types = network
        .species
        .iter()
        .map(|s| SpeciesType {
            name: species(&s.name),
            attributes: Vec::new(),
        })
        .collect();
    #[allow(clippy::cast_precision_loss, clippy::cast_sign_loss)]
    let initial = network
        .species
        .iter()
        .filter(|s| s.initial_count > 0)
        .map(|s| Initial {
            species_type: species(&s.name),
            attributes: BTreeMap::new(),
            // Counts up to 2^53 (checked) are exact in f64.
            count: Value::Number(s.initial_count as u64 as f64),
        })
        .collect();
    let observables = network
        .species
        .iter()
        .map(|s| {
            (
                species(&s.name),
                Observable {
                    species_type: species(&s.name),
                    label: None,
                    condition: None,
                    weight: None,
                    reduce: Reduce::Sum,
                },
            )
        })
        .collect();
    let rules = network
        .reactions
        .iter()
        .map(|reaction| {
            let order: i64 = reaction.reactants.values().sum();
            let mut reactants = Vec::new();
            if order <= 2 {
                // One count-1 pattern per molecule, as before model/1.0.
                for (name, &count) in &reaction.reactants {
                    for _ in 0..count {
                        let label = ["x", "y"][reactants.len()];
                        reactants.push(Pattern {
                            label: label.into(),
                            species_type: species(name),
                            condition: None,
                            factor: None,
                            count: None,
                        });
                    }
                }
            } else {
                // A higher-order rule: one pattern per species with its count
                // (validated stoichiometry is 1 to 2^31 - 1).
                let many = reaction.reactants.len() > 2;
                for (i, (name, &count)) in reaction.reactants.iter().enumerate() {
                    let label = if many {
                        format!("x{}", i + 1)
                    } else {
                        ["x", "y"][i].to_string()
                    };
                    reactants.push(Pattern {
                        label,
                        species_type: species(name),
                        condition: None,
                        factor: None,
                        count: (count > 1).then(|| u64::try_from(count).expect("validated")),
                    });
                }
            }
            let mut products = Vec::new();
            for (name, &count) in &reaction.products {
                for _ in 0..count {
                    products.push(Product {
                        species_type: species(name),
                        copy: None,
                        attributes: BTreeMap::new(),
                    });
                }
            }
            let s = |value: &Scalar| names.scalar(value);
            let (rate, propensity) = match &reaction.rate {
                Rate::MassAction { coefficient } => (Some(s(coefficient)), None),
                Rate::HillActivation {
                    regulator,
                    vmax,
                    basal,
                    k,
                    n,
                } => (
                    None,
                    Some(format!(
                        "{} + {} / (1 + ({} / {})^{})",
                        s(basal),
                        s(vmax),
                        s(k),
                        species(regulator),
                        s(n)
                    )),
                ),
                Rate::HillRepression {
                    regulator,
                    vmax,
                    basal,
                    k,
                    n,
                } => (
                    None,
                    Some(format!(
                        "{} + {} / (1 + ({} / {})^{})",
                        s(basal),
                        s(vmax),
                        species(regulator),
                        s(k),
                        s(n)
                    )),
                ),
                Rate::MichaelisMenten {
                    substrate,
                    vmax,
                    km,
                } => (
                    None,
                    Some(format!(
                        "{} * ({} / ({} + {}))",
                        s(vmax),
                        species(substrate),
                        s(km),
                        species(substrate)
                    )),
                ),
            };
            Rule {
                name: names.reactions[&reaction.name].clone(),
                reactants,
                rate,
                propensity,
                products,
            }
        })
        .collect();
    Model {
        schema: CONVERTED_VERSION,
        name: names.model.clone(),
        description: Some(description(network, &names)),
        parameters: network
            .parameters
            .iter()
            .map(|(name, &value)| (names.parameters[name].clone(), value))
            .collect(),
        types,
        initial,
        observables,
        rules,
        stop: None,
        atp: None,
        continuous: None,
    }
}

#[cfg(test)]
mod tests {
    //! The converted propensities against the network engine's rate laws.
    use super::from_network;
    use crate::rules::expr::Env;
    use crate::rules::{flat, model};
    use serde_json::{Value, json};

    /// Declared in docs/campaign/merge-0.7/M2-log.md: relative tolerance of
    /// a converted Hill propensity against `rates::hill`, per unit of the
    /// Hill exponent n (with n below 1 counted as 1).
    const HILL_RELATIVE_PER_EXPONENT: f64 = 1e-14;
    /// Declared there too: where the fraction kⁿ/(kⁿ + xⁿ) or its
    /// complement falls below about 2^-1024, the converted form gives 0 and
    /// the network engine a value below vmax · 2^-1020.
    const HILL_ABSOLUTE_PER_VMAX: f64 = 8.900_295_434_028_806e-308; // 2^-1020
    /// Declared there too: mass action over pairs, relative.
    const PAIR_RELATIVE: f64 = 4.0 * f64::EPSILON;

    fn network(species: &[(&str, i64)], reaction: &Value) -> Value {
        json!({
            "schema": "crn-craft.model/0.1", "name": "law",
            "convention": "stochastic_counts", "omega": 1.0, "parameters": {},
            "species": species.iter().map(|(name, count)| json!({"name": name, "initial_count": count})).collect::<Vec<_>>(),
            "reactions": [reaction]
        })
    }

    /// The network engine's propensity of reaction 0 and the converted
    /// rules model's channel 0 at the same counts.
    fn both(model: &Value, counts: &[i64]) -> (Result<f64, String>, f64) {
        let bytes = serde_json::to_vec(model).unwrap();
        let old = crate::CompiledNetwork::compile(&crate::parse_model(&bytes).unwrap())
            .unwrap()
            .propensity(0, counts)
            .map_err(|d| d.message);
        let converted = from_network(&bytes).unwrap();
        let rules: model::Model = serde_json::from_slice(&converted).unwrap();
        let compiled = model::compile(&rules).unwrap();
        let net = flat::Network::compile(&compiled).unwrap();
        assert_eq!(net.species.len(), counts.len(), "one tuple per species");
        let counts: Vec<u64> = counts.iter().map(|&n| n as u64).collect();
        let mut values = vec![0.0; compiled.observables.len()];
        net.refresh(&counts, &mut values).ok().unwrap();
        let env = Env {
            params: &compiled.params,
            observables: &values,
            attributes: [&[], &[]],
        };
        let new = net.rate(0, &counts, &env).ok().unwrap();
        (old, new)
    }

    #[test]
    #[allow(clippy::cast_precision_loss)]
    fn hill_propensities_match_the_network_law() {
        let counts = [
            0_i64,
            1,
            2,
            3,
            7,
            39,
            40,
            41,
            1000,
            8191,
            8192,
            1 << 20,
            1 << 40,
            1 << 53,
        ];
        let mut compared = 0;
        let mut floor = 0;
        // Largest relative difference per unit of max(1, n).
        let mut worst = 0.0_f64;
        for repression in [false, true] {
            for (vmax, basal) in [
                (216.0, 0.216),
                (3.0, 0.0),
                (1e300, 0.0),
                (5.0, 1e-300),
                (1e-3, 0.0),
                (0.0, 2.0),
            ] {
                for k in [1.0, 7.5, 40.0, 1000.0, 1e6] {
                    for n in [0.5_f64, 1.0, 2.0, 3.7, 10.0, 60.0, 100.0] {
                        let kind = if repression {
                            "hill_repression"
                        } else {
                            "hill_activation"
                        };
                        let model = network(
                            &[("R", 1), ("G", 1)],
                            &json!({"name": "law", "reactants": {}, "products": {"G": 1},
                                    "rate": {"kind": kind, "regulator": "R", "vmax": vmax,
                                             "basal": basal, "k": k, "n": n}}),
                        );
                        for &x in &counts {
                            let (old, new) = both(&model, &[x, 1]);
                            assert!(new.is_finite() && new >= 0.0);
                            let tolerance = HILL_RELATIVE_PER_EXPONENT * n.max(1.0);
                            match old {
                                Ok(old) => {
                                    let difference = (new - old).abs();
                                    if difference <= tolerance * old {
                                        compared += 1;
                                        if old > 0.0 {
                                            worst = worst.max(difference / old / n.max(1.0));
                                        }
                                    } else {
                                        assert!(
                                            difference <= vmax * HILL_ABSOLUTE_PER_VMAX,
                                            "{kind} vmax {vmax} basal {basal} k {k} n {n} x {x}: {old} vs {new}"
                                        );
                                        floor += 1;
                                    }
                                }
                                // The network engine refuses a positive
                                // propensity that underflows to 0.
                                Err(message) => {
                                    assert!(message.contains("underflows"), "{message}");
                                    assert!(new <= basal + vmax * HILL_ABSOLUTE_PER_VMAX);
                                    floor += 1;
                                }
                            }
                            // x = 0 and vmax = 0 are exact in both forms.
                            if x == 0 || vmax == 0.0 {
                                let exact = basal + if repression { vmax } else { 0.0 };
                                assert_eq!(new.to_bits(), exact.to_bits());
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "Hill: {compared} cases within the relative tolerance (largest {worst:e} per unit of n), {floor} at the floor"
        );
        // Most cases are within the relative tolerance.
        assert!(
            compared > 20 * floor,
            "{compared} relative, {floor} at the floor"
        );
    }

    #[test]
    fn michaelis_menten_is_bit_identical() {
        for (vmax, km) in [(5.0, 20.0), (1e300, 1e-3), (1e-3, 1e6), (0.0, 1.0)] {
            let model = network(
                &[("S", 1), ("P", 1)],
                &json!({"name": "law", "reactants": {"S": 1}, "products": {"P": 1},
                        "rate": {"kind": "michaelis_menten", "substrate": "S",
                                 "vmax": vmax, "km": km}}),
            );
            for x in [0, 1, 2, 19, 20, 21, 1000, 1 << 40, 1 << 53] {
                let (old, new) = both(&model, &[x, 1]);
                assert_eq!(old.unwrap().to_bits(), new.to_bits(), "{vmax} {km} {x}");
            }
        }
    }

    #[test]
    #[allow(clippy::cast_precision_loss)]
    fn mass_action_matches_the_network_kernels() {
        let counts = [
            0_i64,
            1,
            2,
            3,
            4,
            5,
            10,
            50,
            1000,
            1 << 26,
            1 << 40,
            1 << 53,
        ];
        for c in [1e-7, 0.1, 0.5, 1.0, 3.0, 2.5e10] {
            for (species, reactants, exact) in [
                (vec![("A", 1)], json!({}), true),
                (vec![("A", 1)], json!({"A": 1}), true),
                (vec![("A", 1), ("B", 1)], json!({"A": 1, "B": 1}), false),
                (vec![("A", 1)], json!({"A": 2}), false),
            ] {
                let model = network(
                    &species,
                    &json!({"name": "law", "reactants": reactants, "products": {},
                            "rate": {"kind": "mass_action", "coefficient": c}}),
                );
                for &x in &counts {
                    for &y in &counts[..if species.len() == 2 { counts.len() } else { 1 }] {
                        let state: Vec<i64> = if species.len() == 2 {
                            vec![x, y]
                        } else {
                            vec![x]
                        };
                        let (old, new) = both(&model, &state);
                        let old = old.unwrap();
                        if exact {
                            assert_eq!(old.to_bits(), new.to_bits(), "{c} {reactants} {state:?}");
                        } else {
                            assert!(
                                (new - old).abs() <= PAIR_RELATIVE * old,
                                "{c} {reactants} {state:?}: {old} vs {new}"
                            );
                        }
                    }
                }
            }
        }
        // The pair convention n(n - 1)/2 equals the network's dimer kernel
        // exactly at small counts.
        let model = network(
            &[("A", 1)],
            &json!({"name": "law", "reactants": {"A": 2}, "products": {},
                    "rate": {"kind": "mass_action", "coefficient": 1.0}}),
        );
        for n in 0..=50_i64 {
            let (old, new) = both(&model, &[n]);
            assert_eq!(old.unwrap(), (n * (n - 1).max(0) / 2) as f64);
            assert_eq!(new, (n * (n - 1).max(0) / 2) as f64);
        }
    }
}
