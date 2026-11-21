//! Scheduled count interventions and fixed reaction completions for closed networks.
use super::{Method, RunRequest, RunVersion, Selector, flat::Network, model, request_error};
use crate::diagnostic::Diagnostic;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// External actions and delayed completions in a `crn-craft.run/1.1` request.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Protocol {
    /// Unit of action times, the horizon and model rates. No implicit conversion.
    #[serde(
        default,
        deserialize_with = "super::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "String")]
    pub time_unit: Option<String>,
    /// Atomic groups in nondecreasing time order.
    #[serde(default)]
    pub actions: Vec<Action>,
    /// Completion attached to the initiation of a named reaction.
    #[serde(default)]
    pub delays: Vec<Delay>,
    /// Maximum number of outstanding completions, including queue bookkeeping.
    #[serde(default = "default_pending")]
    pub max_pending: usize,
}

impl Default for Protocol {
    fn default() -> Self {
        Self {
            time_unit: None,
            actions: Vec::new(),
            delays: Vec::new(),
            max_pending: default_pending(),
        }
    }
}

const fn default_pending() -> usize {
    4096
}

/// An atomic action group. Count operations refer to the pre-group counts.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Action {
    /// Simulated time in the same unit as the run.
    pub time: f64,
    /// Replace named parameters for all later events, until changed again.
    #[serde(default)]
    pub parameters: BTreeMap<String, f64>,
    /// Signed count increments.
    #[serde(default)]
    pub change: BTreeMap<String, i64>,
    /// Assign absolute counts.
    #[serde(default)]
    pub set: BTreeMap<String, u64>,
    /// Keep each copy independently with the given probability.
    #[serde(default)]
    pub partition: BTreeMap<String, f64>,
}

/// Fixed-delay completion. With an owner, the firing creates one busy molecule
/// and completion consumes that molecule. Loss of the busy molecule cancels it.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Delay {
    /// Initiation reaction, which still applies its immediate stoichiometry.
    pub reaction: String,
    /// Positive delay in the same time unit as the run.
    pub duration: f64,
    /// Counts produced on completion.
    #[serde(default)]
    pub products: BTreeMap<String, u64>,
    /// An optional busy species, private to the delay rule.
    #[serde(
        default,
        deserialize_with = "super::present",
        skip_serializing_if = "Option::is_none"
    )]
    #[schemars(with = "String")]
    pub owner: Option<String>,
}

pub(super) struct Phase {
    pub model: model::Compiled,
    pub network: Network,
}

pub(super) struct Boundary {
    pub time: f64,
    pub operations: Vec<(usize, Operation)>,
    pub phase: usize,
}

pub(super) enum Operation {
    Change(i64),
    Set(u64),
    Partition(f64),
}

pub(super) struct Completion {
    pub duration: f64,
    pub owner: Option<usize>,
    pub products: Vec<(usize, u64)>,
}

pub(super) struct Schedule {
    pub boundaries: Vec<Boundary>,
    pub phases: Vec<Phase>,
    pub delays: Vec<Option<Completion>>,
    pub owned: BTreeSet<usize>,
    pub max_pending: usize,
}

fn error(path: &str, message: impl Into<String>) -> Diagnostic {
    request_error(
        path,
        message,
        "Use a closed count model and supported protocol actions.",
    )
}

pub(super) fn validate(request: &RunRequest) -> Result<(), Diagnostic> {
    if request.schema == RunVersion::V11 && request.protocol.is_none() {
        return Err(error("/protocol", "run/1.1 requires a protocol"));
    }
    let Some(protocol) = &request.protocol else {
        return Ok(());
    };
    if request.schema != RunVersion::V11 {
        return Err(error("/schema", "protocol requires crn-craft.run/1.1"));
    }
    if !matches!(
        request.method,
        Some(
            Method::StaticDirect {
                selector: None | Some(Selector::Linear),
                ..
            } | Method::NextReaction {}
        )
    ) {
        return Err(error(
            "/method",
            "protocol supports linear Static Direct or next reaction",
        ));
    }
    if protocol.actions.len() > 4096
        || protocol.delays.len() > 4096
        || !(1..=4096).contains(&protocol.max_pending)
    {
        return Err(error(
            "/protocol",
            "protocol action or completion capacity exceeded",
        ));
    }
    if protocol
        .time_unit
        .as_ref()
        .is_some_and(|s| s.trim().is_empty() || s.len() > 64)
    {
        return Err(error("/protocol/time_unit", "invalid time unit label"));
    }
    let mut previous = 0.0;
    for (i, action) in protocol.actions.iter().enumerate() {
        if !action.time.is_finite() || action.time < previous || action.time > request.time.end {
            return Err(error(
                &format!("/protocol/actions/{i}/time"),
                "action time must be ordered and within the horizon",
            ));
        }
        previous = action.time;
    }
    Ok(())
}

fn index(
    model: &model::Model,
    net: &Network,
    label: &str,
    path: &str,
) -> Result<usize, Diagnostic> {
    net.species
        .iter()
        .position(|s| model.types[s.ty].name == label)
        .ok_or_else(|| error(path, format!("unknown species {label:?}")))
}

fn operations(
    action: &Action,
    source: &model::Model,
    net: &Network,
    path: &str,
) -> Result<Vec<(usize, Operation)>, Diagnostic> {
    let mut ops = BTreeMap::new();
    let mut add = |label: &str, operation| {
        let i = index(source, net, label, path)?;
        if ops.insert(i, operation).is_some() {
            return Err(error(path, format!("conflicting operations on {label:?}")));
        }
        Ok(())
    };
    for (label, &n) in &action.change {
        add(label, Operation::Change(n))?;
    }
    for (label, &n) in &action.set {
        if n > model::MAX_COUNT {
            return Err(error(path, "absolute count exceeds the model limit"));
        }
        add(label, Operation::Set(n))?;
    }
    for (label, &p) in &action.partition {
        if !(0.0..=1.0).contains(&p) {
            return Err(error(path, "partition probability must be in [0, 1]"));
        }
        add(label, Operation::Partition(p))?;
    }
    Ok(ops.into_iter().collect())
}

pub(super) fn include_empty_types(compiled: &mut model::Compiled) {
    for ty in 0..compiled.types.len() {
        if !compiled.initial.iter().any(|entry| entry.0 == ty) {
            compiled.initial.push((ty, Vec::new(), 0));
        }
    }
}

impl Schedule {
    pub fn build(
        protocol: &Protocol,
        source: &model::Model,
        base: &Network,
        observations: &super::observe::Spec,
    ) -> Result<Self, Diagnostic> {
        let mut source = source.clone();
        let mut phases = Vec::new();
        let mut boundaries = Vec::with_capacity(protocol.actions.len());
        // Distinct parameter changes require their own compiled rate constants.
        let work = protocol
            .actions
            .iter()
            .filter(|a| !a.parameters.is_empty())
            .count()
            .saturating_mul(base.species.len() + base.channels.len());
        if work > 1_000_000 {
            return Err(error(
                "/protocol",
                "phase lowering exceeds its work capacity",
            ));
        }
        for (i, action) in protocol.actions.iter().enumerate() {
            let path = format!("/protocol/actions/{i}");
            let operations = operations(action, &source, base, &path)?;
            if !action.parameters.is_empty() {
                for (name, &value) in &action.parameters {
                    let target = source
                        .parameters
                        .get_mut(name)
                        .ok_or_else(|| error(&path, format!("unknown parameter {name:?}")))?;
                    if !value.is_finite() {
                        return Err(error(&path, "non-finite parameter"));
                    }
                    *target = value;
                }
                let mut model = model::compile(&source)
                    .map_err(|d| error(&path, format!("invalid phase: {}", d[0].message)))?;
                include_empty_types(&mut model);
                let mut network = Network::compile(&model).map_err(|e| error(&path, e))?;
                if !same_layout(base, &network) {
                    return Err(error(
                        &path,
                        "a kinetic edit changes initial counts or channel topology",
                    ));
                }
                network.plan_observables(&model, observations);
                phases.push(Phase { model, network });
            }
            boundaries.push(Boundary {
                time: action.time,
                operations,
                phase: phases.len(),
            });
        }

        let (delays, owned) = completions(protocol, &source, base)?;
        Ok(Self {
            boundaries,
            phases,
            delays,
            owned,
            max_pending: protocol.max_pending,
        })
    }
}

fn same_layout(a: &Network, b: &Network) -> bool {
    a.species.len() == b.species.len()
        && a.channels.len() == b.channels.len()
        && a.species
            .iter()
            .zip(&b.species)
            .all(|(a, b)| (a.ty, a.initial) == (b.ty, b.initial))
        && a.channels.iter().zip(&b.channels).all(|(a, b)| {
            a.rule == b.rule
                && a.reactants == b.reactants
                && a.multi == b.multi
                && a.products == b.products
        })
}

type CompletionSpec = (Vec<Option<Completion>>, BTreeSet<usize>);

fn completions(
    protocol: &Protocol,
    source: &model::Model,
    base: &Network,
) -> Result<CompletionSpec, Diagnostic> {
    let mut delays: Vec<_> = (0..base.channels.len()).map(|_| None).collect();
    let mut owned = BTreeSet::new();
    for (i, spec) in protocol.delays.iter().enumerate() {
        let path = format!("/protocol/delays/{i}");
        let mut matching = base
            .channels
            .iter()
            .enumerate()
            .filter(|(_, c)| source.rules[c.rule].name == spec.reaction);
        let Some((channel, _)) = matching.next() else {
            return Err(error(&path, "unknown initiation reaction"));
        };
        if matching.next().is_some()
            || delays[channel].is_some()
            || !spec.duration.is_finite()
            || spec.duration <= 0.0
        {
            return Err(error(
                &path,
                "initiation needs one channel, one positive delay and one declaration",
            ));
        }
        let owner = spec
            .owner
            .as_deref()
            .map(|name| index(source, base, name, &path))
            .transpose()?;
        if let Some(s) = owner {
            if base.species[s].initial != 0 || !owned.insert(s) {
                return Err(error(
                    &path,
                    "busy species must start empty and belong to one delay",
                ));
            }
            for (k, c) in base.channels.iter().enumerate() {
                let change = c
                    .changes
                    .iter()
                    .find(|&(id, _)| *id == s)
                    .map_or(0, |&(_, n)| n);
                if (k == channel && change != 1) || (k != channel && change > 0) {
                    return Err(error(
                        &path,
                        "only the initiation can make the busy species, one at a time",
                    ));
                }
            }
        }
        let products = spec
            .products
            .iter()
            .map(|(name, &n)| {
                if n == 0 || n > model::MAX_COUNT {
                    return Err(error(&path, "invalid completion count"));
                }
                Ok((index(source, base, name, &path)?, n))
            })
            .collect::<Result<Vec<_>, Diagnostic>>()?;
        delays[channel] = Some(Completion {
            duration: spec.duration,
            owner,
            products,
        });
    }
    if delays
        .iter()
        .flatten()
        .any(|d| d.products.iter().any(|(s, _)| owned.contains(s)))
    {
        return Err(error(
            "/protocol/delays",
            "a completion cannot produce a busy species",
        ));
    }
    Ok((delays, owned))
}
