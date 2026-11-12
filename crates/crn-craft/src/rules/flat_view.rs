//! A read-only view of a closed model lowered to a flat network.
//!
//! The static methods (`static_direct`, `next_reaction`, `tau`, ...) run a
//! closed rules model as a fixed network: one *species* per instance tuple
//! (a type and its attribute values) that can ever hold molecules, and one
//! *channel* per rule and choice of reactant tuples. [`FlatView`] exposes
//! that network to the features that need a fixed species list: finite
//! state projection ([`crate::fsp`]), network analysis, ODE right-hand sides
//! and exporters.
//!
//! The view is built from the same lowering the Static Direct engine runs,
//! and [`FlatView::propensity`] evaluates a channel with the engine's own
//! channel-rate function, so it gives exactly the value Static Direct uses
//! in that state, including higher-order (`count`) channels and whole-rule
//! `propensity` expressions.
//!
//! A model that is not closed (its products can create attribute tuples
//! that are not known before the run), a model with an `atp` block or with
//! continuous variables, and a model above the static lowering limits has
//! no flat view: [`FlatView::new`] refuses it with `CRN-R017`.
//!
//! ```
//! use crn_craft::rules::FlatView;
//! use serde_json::json;
//!
//! let model = serde_json::to_vec(&json!({
//!     "schema": "crn-craft.model/1.0", "name": "dimerization",
//!     "parameters": {"k": 0.5},
//!     "types": [{"name": "A"}, {"name": "B"}],
//!     "initial": [{"type": "A", "count": 10}],
//!     "observables": {"A": {"type": "A"}, "B": {"type": "B"}},
//!     "rules": [{"name": "bind", "reactants": [{"label": "x", "type": "A", "count": 2}],
//!                "rate": "k", "products": [{"type": "B"}]}]
//! }))
//! .unwrap();
//! let view = FlatView::from_bytes(&model).unwrap();
//! assert_eq!(view.species(), ["A", "B"]);
//! assert_eq!(view.initial(), [10, 0]);
//! assert_eq!(view.changes(0), [(0, -2), (1, 1)]);
//! assert_eq!(view.reactants(0), [(0, 2)]);
//! // k * C(10, 2) = 0.5 * 45
//! assert_eq!(view.propensity(0, &[10, 0]).unwrap(), 22.5);
//! ```
use super::engine::Halt;
use super::expr::{Env, Expr, Slot, python_literal};
use super::flat::{Network, RateCache};
use super::model::{self, Compiled, MAX_COUNT, Model, Reduce};
use crate::diagnostic::Diagnostic;
use std::collections::HashMap;

/// The diagnostic code of a model that has no flat view.
pub const NO_FLAT_VIEW: &str = "CRN-R017";

/// A read-only view of a closed model lowered to a flat network.
///
/// Species are numbered in *view order*: by type in the model's declaration
/// order, and within a type in the order the lowering finds the tuples (the
/// initial entries first, then the tuples products create, by rule). Only
/// tuples that can hold molecules are species: a type that has no initial
/// instance and that no rule produces has none, except in a view built with
/// [`FlatView::with_every_type`], where every type without attributes has
/// its species. Channels keep the order of
/// the lowering: by rule, then by choice of reactant tuples. All species
/// indices this view takes or returns are view indices; all channel indices
/// are channel numbers of the lowered network.
pub struct FlatView {
    source: Model,
    compiled: Compiled,
    net: Network,
    labels: Vec<String>,
    index: HashMap<String, usize>,
    initial: Vec<u64>,
    /// View index to network species.
    order: Vec<usize>,
    /// Network species to view index.
    position: Vec<usize>,
    identity: bool,
    reactants: Vec<Vec<(usize, u64)>>,
    products: Vec<Vec<(usize, u64)>>,
    changes: Vec<Vec<(usize, i64)>>,
    reads: Vec<Vec<usize>>,
    /// The projection of each model observable in the lowered network.
    projections: Vec<Option<usize>>,
}

impl std::fmt::Debug for FlatView {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FlatView")
            .field("model", &self.source.name)
            .field("species", &self.labels)
            .field("channels", &self.changes.len())
            .finish_non_exhaustive()
    }
}

fn refusal(path: &str, message: impl Into<String>, hint: &str) -> Diagnostic {
    Diagnostic::error(NO_FLAT_VIEW, path, message, hint)
}

fn halt(error: Halt) -> Diagnostic {
    match error {
        Halt::Limit(_, diagnostic) | Halt::Numerical(diagnostic) => diagnostic,
    }
}

/// Parses a rules model strictly, as `prepare` does: a pattern `count`
/// before `crn-craft.model/1.0` gets the 0.6.0 unknown-field error.
pub(crate) fn parse(model_bytes: &[u8]) -> Result<Model, Diagnostic> {
    let model: Model = crate::strict_json::parse(model_bytes)
        .map_err(|error| super::network_model_error(error, model_bytes))?;
    match model::legacy_count(&model) {
        Some(diagnostic) => Err(diagnostic),
        None => Ok(model),
    }
}

/// An attribute value in a species label: integers below 2^53 without a
/// fraction, other values in Rust's shortest round-trip form (`0.5`,
/// `1e-7`).
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn number(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < MAX_COUNT as f64 {
        format!("{}", x as i64)
    } else {
        format!("{x:?}")
    }
}

/// Merged (species, count) entries, sorted by species.
fn merged(entries: impl IntoIterator<Item = (usize, u64)>) -> Vec<(usize, u64)> {
    let mut list: Vec<(usize, u64)> = entries.into_iter().collect();
    list.sort_by_key(|&(s, _)| s);
    let mut out: Vec<(usize, u64)> = Vec::with_capacity(list.len());
    for (s, k) in list {
        match out.last_mut() {
            Some(last) if last.0 == s => last.1 += k,
            _ => out.push((s, k)),
        }
    }
    out
}

impl FlatView {
    /// Parses model bytes strictly and builds their view; see [`Self::new`].
    ///
    /// # Errors
    ///
    /// `CRN-E001` to `CRN-E004` for a document that does not parse, and the
    /// errors of [`Self::new`].
    pub fn from_bytes(model_bytes: &[u8]) -> Result<Self, Vec<Diagnostic>> {
        Self::new(&parse(model_bytes).map_err(|d| vec![d])?)
    }

    /// Validates and compiles a rules model (any `crn-craft.rules-model`
    /// version or `crn-craft.model/1.0`) and lowers it as the static methods
    /// do.
    ///
    /// # Errors
    ///
    /// The model's validation diagnostics, or `CRN-R017` when it has no flat
    /// view: an `atp` block (at `/atp`), continuous variables (at
    /// `/continuous`), or a lowering failure (at `/rules`: products that can
    /// create new attribute tuples, random product attributes, a whole-rule
    /// propensity over a type with several tuples, the static tuple and
    /// channel limits).
    pub fn new(model: &Model) -> Result<Self, Vec<Diagnostic>> {
        Self::lower(model, false)
    }

    /// Like [`Self::new`], but every type without attributes is a species,
    /// also one that never holds molecules (its count stays 0), and each
    /// rule over such types has its channels. A converted network model
    /// keeps every species and reaction of the network this way: network
    /// analysis and the GillesPy2 export use this view, finite state
    /// projection uses [`Self::new`].
    ///
    /// # Errors
    ///
    /// As [`Self::new`].
    pub fn with_every_type(model: &Model) -> Result<Self, Vec<Diagnostic>> {
        Self::lower(model, true)
    }

    /// Builds the view of [`Self::with_every_type`] from the bytes of any
    /// model the engine reads: a network `crn-craft.model/0.1` (converted
    /// first with [`super::convert::from_network`]), a rules model or a
    /// `crn-craft.model/1.0`.
    ///
    /// # Errors
    ///
    /// The converter's diagnostics for a network model (its validation
    /// errors, or `CRN-R016`), `CRN-E001` to `CRN-E004` for a document that
    /// does not parse, and the errors of [`Self::new`].
    pub fn of_any_model(model_bytes: &[u8]) -> Result<Self, Vec<Diagnostic>> {
        if crate::schema_field(model_bytes).as_deref() == Some("crn-craft.model/0.1") {
            let converted = super::convert::from_network(model_bytes)?;
            return Self::with_every_type(&parse(&converted).map_err(|d| vec![d])?);
        }
        Self::with_every_type(&parse(model_bytes).map_err(|d| vec![d])?)
    }

    fn lower(model: &Model, every_type: bool) -> Result<Self, Vec<Diagnostic>> {
        let mut compiled = model::compile(model)?;
        let hint = "FSP, analysis, ODE and export need a closed model of counts; simulate this model with its own method instead.";
        if compiled.atp.is_some() {
            return Err(vec![refusal(
                "/atp",
                "a model with an ATP block has no flat view",
                hint,
            )]);
        }
        if !compiled.continuous.is_empty() {
            return Err(vec![refusal(
                "/continuous",
                "a model with continuous variables has no flat view",
                hint,
            )]);
        }
        if every_type {
            // A zero-count initial entry makes the lowering list the type's
            // one tuple even when no rule produces it.
            for (ty, info) in compiled.types.iter().enumerate() {
                if info.attributes.is_empty() && !compiled.initial.iter().any(|e| e.0 == ty) {
                    compiled.initial.push((ty, Vec::new(), 0));
                }
            }
        }
        let net = Network::compile(&compiled).map_err(|e| {
            vec![refusal(
                "/rules",
                format!("the model has no flat view: {e}"),
                "Use a closed model: every product attribute must be fixed by the parameters or copied from a reactant tuple of the lowered network.",
            )]
        })?;
        Ok(Self::build(model.clone(), compiled, net))
    }

    fn build(source: Model, compiled: Compiled, net: Network) -> Self {
        let mut order: Vec<usize> = (0..net.species.len()).collect();
        order.sort_by_key(|&s| net.species[s].ty);
        let mut position = vec![0; order.len()];
        for (v, &s) in order.iter().enumerate() {
            position[s] = v;
        }
        let identity = order.iter().enumerate().all(|(v, &s)| v == s);
        let labels: Vec<String> = order
            .iter()
            .map(|&s| {
                let species = &net.species[s];
                let info = &compiled.types[species.ty];
                if info.attributes.is_empty() {
                    info.name.clone()
                } else {
                    let values: Vec<String> = info
                        .attributes
                        .iter()
                        .zip(&species.attrs)
                        .map(|(name, &x)| format!("{name}={}", number(x)))
                        .collect();
                    format!("{}{{{}}}", info.name, values.join(","))
                }
            })
            .collect();
        let index = labels
            .iter()
            .enumerate()
            .map(|(i, l)| (l.clone(), i))
            .collect();
        let initial = order.iter().map(|&s| net.species[s].initial).collect();
        let at = |s: usize| position[s];
        let reactants = net
            .channels
            .iter()
            .map(|c| match &c.multi {
                Some(list) => merged(list.iter().map(|&(s, k)| (at(s), k))),
                None => merged(c.reactants.iter().map(|&s| (at(s), 1))),
            })
            .collect();
        let products = net
            .channels
            .iter()
            .map(|c| merged(c.products.iter().map(|&s| (at(s), 1))))
            .collect();
        let changes = net
            .channels
            .iter()
            .map(|c| {
                let mut list: Vec<(usize, i64)> =
                    c.changes.iter().map(|&(s, d)| (at(s), d)).collect();
                list.sort_by_key(|&(s, _)| s);
                list
            })
            .collect();
        let mut projections = vec![None; compiled.observables.len()];
        for (p, (o, _, _)) in net.observables.iter().enumerate() {
            projections[*o] = Some(p);
        }
        // The species each rule's rate reads through its observables, as the
        // lowering computes them for `affected`.
        let rule_reads: Vec<Vec<usize>> = compiled
            .rules
            .iter()
            .map(|rule| {
                rule.rate
                    .read_list()
                    .into_iter()
                    .filter_map(|o| projections[o])
                    .flat_map(|p| net.observables[p].2.iter().map(|&(s, _)| s))
                    .collect()
            })
            .collect();
        let reads = net
            .channels
            .iter()
            .map(|c| {
                let mut list: Vec<usize> = rule_reads[c.rule]
                    .iter()
                    .chain(&c.reactants)
                    .map(|&s| at(s))
                    .collect();
                list.sort_unstable();
                list.dedup();
                list
            })
            .collect();
        Self {
            source,
            compiled,
            net,
            labels,
            index,
            initial,
            order,
            position,
            identity,
            reactants,
            products,
            changes,
            reads,
            projections,
        }
    }

    /// The source model.
    #[must_use]
    pub fn model(&self) -> &Model {
        &self.source
    }

    /// Species labels in view order: the type name for a type without
    /// attributes (`A`), otherwise the type name and every attribute in
    /// declaration order (`mtDNA{wt=1,mut=0}`). Values are written as
    /// integers when they are integers below 2^53, else in Rust's shortest
    /// round-trip form (`0.5`, `1e-7`).
    #[must_use]
    pub fn species(&self) -> &[String] {
        &self.labels
    }

    /// The view index of a species label.
    #[must_use]
    pub fn index(&self, label: &str) -> Option<usize> {
        self.index.get(label).copied()
    }

    /// The type of each species, as an index into the model's `types`.
    #[must_use]
    pub fn species_type(&self, species: usize) -> usize {
        self.net.species[self.order[species]].ty
    }

    /// The attribute values of a species, in the type's declaration order.
    #[must_use]
    pub fn attributes(&self, species: usize) -> &[f64] {
        &self.net.species[self.order[species]].attrs
    }

    /// Initial counts in view order.
    #[must_use]
    pub fn initial(&self) -> &[u64] {
        &self.initial
    }

    /// The number of channels.
    #[must_use]
    pub fn channels(&self) -> usize {
        self.changes.len()
    }

    /// The rule a channel belongs to, as an index into the model's `rules`.
    #[must_use]
    pub fn rule(&self, channel: usize) -> usize {
        self.net.channels[channel].rule
    }

    /// The name of a rule.
    #[must_use]
    pub fn rule_name(&self, rule: usize) -> &str {
        &self.compiled.rules[rule].name
    }

    /// The number of rules (channels may be fewer or more).
    #[must_use]
    pub fn rules(&self) -> usize {
        self.compiled.rules.len()
    }

    /// The net change of one firing: (species, change) sorted by species,
    /// zero changes left out.
    #[must_use]
    pub fn changes(&self, channel: usize) -> &[(usize, i64)] {
        &self.changes[channel]
    }

    /// The reactants of a channel: (species, molecules taken) sorted by
    /// species. A pattern with `count` k, or two patterns over one tuple,
    /// take several molecules of one species.
    #[must_use]
    pub fn reactants(&self, channel: usize) -> &[(usize, u64)] {
        &self.reactants[channel]
    }

    /// The products of a channel: (species, molecules made) sorted by
    /// species.
    #[must_use]
    pub fn products(&self, channel: usize) -> &[(usize, u64)] {
        &self.products[channel]
    }

    /// Whether the channel's rule has a whole-rule `propensity`, which
    /// applies once while every reactant count is available. Otherwise the
    /// rule's `rate` is multiplied by the number of reactant combinations,
    /// Π C(n, k), and by [`Self::factor`].
    #[must_use]
    pub fn whole(&self, channel: usize) -> bool {
        self.net.channels[channel].whole
    }

    /// The product of the reactant patterns' local factors (`factor`
    /// expressions over the tuples' attributes); 1 for patterns without one.
    #[must_use]
    pub fn factor(&self, channel: usize) -> f64 {
        self.net.channels[channel].factor
    }

    /// The species whose counts the channel's propensity reads, sorted: its
    /// reactants and the species of every observable its rule's rate or
    /// propensity reads.
    #[must_use]
    pub fn reads(&self, channel: usize) -> &[usize] {
        &self.reads[channel]
    }

    /// The channels whose propensity can change when `channel` fires, sorted:
    /// those that read a species the firing changes. Empty for a channel
    /// without net change. This is the list Static Direct updates after an
    /// event.
    #[must_use]
    pub fn dependencies(&self, channel: usize) -> &[usize] {
        &self.net.affected[channel]
    }

    /// The source text of a rule's `rate`, or of its `propensity` when the
    /// rule has one ([`Self::whole`]).
    #[must_use]
    pub fn rate_source(&self, rule: usize) -> &str {
        let rule = &self.source.rules[rule];
        rule.propensity
            .as_deref()
            .or(rule.rate.as_deref())
            .unwrap_or_default()
    }

    /// The value of a rule's rate when it depends on parameters only (a
    /// constant for the run); `None` when it reads observables.
    #[must_use]
    pub fn rate_constant(&self, rule: usize) -> Option<f64> {
        match &self.net.rates[rule] {
            super::expr::RateKernel::Constant(k) => Some(*k),
            _ => None,
        }
    }

    /// The model observables a rule's rate or propensity reads, sorted, as
    /// indices into [`Self::observable_names`].
    #[must_use]
    pub fn rate_observables(&self, rule: usize) -> Vec<usize> {
        self.compiled.rules[rule].rate.read_list()
    }

    /// The compiled rate or propensity expression of a rule. Its parameter
    /// slots follow [`Self::parameter_names`] and its observable slots
    /// [`Self::observable_names`]. The kinetics evaluate it with dual
    /// numbers for their Jacobians.
    #[must_use]
    pub(crate) fn rate_expression(&self, rule: usize) -> &Expr {
        &self.compiled.rules[rule].rate
    }

    /// A rule's rate or propensity as a Python expression over the species
    /// `S0`, `S1`, … (view indices), for the GillesPy2 export: parameter
    /// values inlined, each observable as the weighted sum of its species,
    /// `^` as `**`. The error names the first construct GillesPy2 cannot
    /// evaluate in both its solvers (comparisons, logic, functions, an
    /// observable that is not a sum).
    pub(crate) fn python_rate(&self, rule: usize) -> Result<String, String> {
        let names = self.observable_names();
        let values = self.parameter_values();
        self.compiled.rules[rule]
            .rate
            .to_python(&mut |slot| match slot {
                Slot::Param(i) => python_literal(values[i]),
                Slot::Observable(o) => {
                    let (reduce, terms) = self.observable(o);
                    if reduce != Reduce::Sum {
                        return Err(format!("the observable `{}`, which is not a sum", names[o]));
                    }
                    if terms.is_empty() {
                        return Ok("0.0".into());
                    }
                    let parts: Vec<String> = terms
                        .iter()
                        .map(|&(s, w)| {
                            if w == 1.0 {
                                Ok(format!("S{s}"))
                            } else {
                                Ok(format!("{} * S{s}", python_literal(w)?))
                            }
                        })
                        .collect::<Result<_, String>>()?;
                    Ok(match parts.as_slice() {
                        [one] if !one.contains(' ') => one.clone(),
                        _ => format!("({})", parts.join(" + ")),
                    })
                }
                Slot::Attribute { .. } => Err("a reactant attribute".into()),
            })
    }

    /// Parameter names in the order expressions index them (sorted).
    #[must_use]
    pub fn parameter_names(&self) -> Vec<&str> {
        self.source.parameters.keys().map(String::as_str).collect()
    }

    /// Parameter values in [`Self::parameter_names`] order.
    #[must_use]
    pub fn parameter_values(&self) -> &[f64] {
        &self.compiled.params
    }

    /// Observable names in the order expressions index them (sorted).
    #[must_use]
    pub fn observable_names(&self) -> Vec<&str> {
        self.compiled
            .observables
            .iter()
            .map(|o| o.name.as_str())
            .collect()
    }

    /// How an observable reduces the species: its reduction and its
    /// (species, weight) terms in increasing species order. An observable
    /// over a type without species has no terms.
    #[must_use]
    pub fn observable(&self, observable: usize) -> (Reduce, Vec<(usize, f64)>) {
        let reduce = self.compiled.observables[observable].reduce;
        let Some(p) = self.projections[observable] else {
            return (reduce, Vec::new());
        };
        let mut terms: Vec<(usize, f64)> = self.net.observables[p]
            .2
            .iter()
            .map(|&(s, w)| (self.position[s], w))
            .collect();
        terms.sort_by_key(|&(s, _)| s);
        (reduce, terms)
    }

    /// Whether the model has a `stop` condition.
    #[must_use]
    pub fn has_stop(&self) -> bool {
        self.compiled.stop_program.is_some()
    }

    /// The propensity of one channel at `counts` (view order), in events per
    /// unit time: exactly the value Static Direct uses in that state.
    ///
    /// # Errors
    ///
    /// `CRN-R030` at `/counts` for a state of the wrong length or with a
    /// count above 2^53, at `/observables` when an observable overflows, and
    /// at the rule's `rate` or `propensity` when the value is negative or
    /// not finite (the replicate failure of Static Direct).
    pub fn propensity(&self, channel: usize, counts: &[u64]) -> Result<f64, Diagnostic> {
        let mut evaluator = self.evaluator();
        evaluator.load(counts)?;
        evaluator.channel(channel)
    }

    /// Every channel's propensity at `counts`; see [`Self::propensity`].
    ///
    /// # Errors
    ///
    /// As [`Self::propensity`], for the first channel that fails.
    pub fn propensities(&self, counts: &[u64]) -> Result<Vec<f64>, Diagnostic> {
        let mut evaluator = self.evaluator();
        evaluator.load(counts)?;
        (0..self.channels()).map(|c| evaluator.channel(c)).collect()
    }

    /// Whether the model's `stop` condition holds at `counts`; `false` for a
    /// model without one. A stopped replicate holds its state, so a stop
    /// state is absorbing.
    ///
    /// # Errors
    ///
    /// As [`Self::propensity`] for the state, and `CRN-R030` at `/stop` when
    /// the condition is NaN.
    pub fn stops(&self, counts: &[u64]) -> Result<bool, Diagnostic> {
        let mut evaluator = self.evaluator();
        evaluator.load(counts)?;
        evaluator.stops()
    }

    /// Scratch space for evaluating many states.
    pub(crate) fn evaluator(&self) -> Evaluator<'_> {
        Evaluator {
            view: self,
            counts: vec![0; self.order.len()],
            values: vec![0.0; self.compiled.observables.len()],
            cache: self.net.rate_cache(),
        }
    }
}

/// Evaluates channels of a [`FlatView`] at one loaded state.
pub(crate) struct Evaluator<'a> {
    view: &'a FlatView,
    /// The loaded state in the lowered network's species order.
    counts: Vec<u64>,
    /// Model observables at the loaded state.
    values: Vec<f64>,
    cache: RateCache,
}

impl Evaluator<'_> {
    /// Loads a state given in view order, with the observables the engine
    /// computes from it.
    pub(crate) fn load(&mut self, counts: &[u64]) -> Result<(), Diagnostic> {
        let view = self.view;
        if counts.len() != view.order.len() {
            return Err(Diagnostic::error(
                "CRN-R030",
                "/counts",
                format!(
                    "a state of {} counts for {} species",
                    counts.len(),
                    view.order.len()
                ),
                "Give one count per species in view order.",
            ));
        }
        if let Some(v) = counts.iter().position(|&n| n > MAX_COUNT) {
            return Err(Diagnostic::error(
                "CRN-R030",
                &crate::diagnostic::pointer("/counts", &view.labels[v]),
                "count above 2^53",
                "Counts are exact integers in f64; keep them at most 2^53.",
            ));
        }
        if view.identity {
            self.counts.copy_from_slice(counts);
        } else {
            for (&s, &n) in view.order.iter().zip(counts) {
                self.counts[s] = n;
            }
        }
        view.net
            .refresh(&self.counts, &mut self.values)
            .map_err(halt)?;
        self.cache.reset();
        Ok(())
    }

    /// The propensity of a channel at the loaded state.
    pub(crate) fn channel(&mut self, channel: usize) -> Result<f64, Diagnostic> {
        let env = Env {
            params: &self.view.compiled.params,
            observables: &self.values,
            attributes: [&[], &[]],
        };
        self.cache
            .channel(&self.view.net, channel, &self.counts, &env)
            .map_err(halt)
    }

    /// The stop condition at the loaded state.
    pub(crate) fn stops(&self) -> Result<bool, Diagnostic> {
        let Some(e) = self.view.compiled.stop_program.as_ref() else {
            return Ok(false);
        };
        let value = e.eval(&Env {
            params: &self.view.compiled.params,
            observables: &self.values,
            attributes: [&[], &[]],
        });
        if value.is_nan() {
            return Err(Diagnostic::error(
                "CRN-R030",
                "/stop",
                "stop condition is undefined (NaN)",
                "Guard divisions so the stop condition is defined in every state.",
            ));
        }
        Ok(value != 0.0)
    }
}
