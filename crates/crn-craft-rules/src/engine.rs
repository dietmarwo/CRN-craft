//! Exact direct-method simulation of one replicate.

#[path = "flat_path.rs"]
mod flat_path;
#[path = "leap.rs"]
mod leap;
#[path = "loop_host.rs"]
mod loop_host;
pub use loop_host::{LoopHalt, LoopHost};
#[path = "mtdna_path.rs"]
mod mtdna_path;
#[path = "thinning.rs"]
mod thinning;

use crate::expr::{Draw, Env};
use crate::model::{Compiled, MAX_COUNT, Reduce, Shape, Source};
use crate::observe::{Event, EventObserver, Monitor, Observation, ObserverAction, Spec};
use crate::tree::{MaxTree, PairChoice, PairTree, SumTree};
use crate::{Control, InstanceRecord, Termination};
use crn_craft::diagnostic::Diagnostic;
use crn_craft::sampling::uniform_closed_open;
use rand_core::RngCore;
use std::collections::HashMap;
use std::time::Instant;

/// A positive weight below this fraction of its total cannot be resolved by a
/// 53-bit uniform.
pub(crate) const RESOLUTION: f64 = 1.0 / 9_007_199_254_740_992.0;
const POLL_EVENTS: u64 = 1024;

pub(crate) struct Limits {
    pub max_events: u64,
    pub max_words: u64,
    pub max_instances: u64,
    pub deadline: Option<Instant>,
}

pub(crate) struct PathOutput {
    pub termination: Termination,
    pub events: u64,
    pub stop_time: f64,
    pub rule_events: Vec<u64>,
    pub rng_words: u64,
    pub peak_instances: u64,
    pub created_instances: u64,
    pub saved_points: usize,
    /// report time × observable, row-major; rows past `saved_points` are 0.
    pub values: Vec<f64>,
    pub diagnostic: Option<Diagnostic>,
    pub observation: Option<Observation>,
    pub instances: Option<Vec<InstanceRecord>>,
}

pub(crate) enum Halt {
    Limit(Termination, Diagnostic),
    Numerical(Diagnostic),
}

pub(crate) fn numerical(path: &str, message: impl Into<String>, hint: &str) -> Halt {
    Halt::Numerical(Diagnostic::error("CRN-R030", path, message, hint))
}

struct TypeState {
    stride: usize,
    attrs: Vec<f64>,
    counts: Vec<u64>,
    free: Vec<usize>,
    index: HashMap<Box<[u64]>, usize>,
}

struct ClassState {
    h: Vec<f64>,
    /// Classes of one-pattern and cross-type rules only.
    single: Option<SumTree>,
    /// Classes of identical-pair rules; also serves single selection.
    pair: Option<PairTree>,
}

impl ClassState {
    fn weight(&self) -> f64 {
        match (&self.single, &self.pair) {
            (Some(t), _) => t.total(),
            (None, Some(t)) => t.total_weight(),
            (None, None) => unreachable!("class without tree"),
        }
    }

    /// Instance holding the weighted molecule at `target` in [0, weight).
    fn select(&self, target: f64) -> usize {
        match (&self.single, &self.pair) {
            (Some(t), _) => t.select(target),
            (None, Some(t)) => t.descend(1, target),
            (None, None) => unreachable!("class without tree"),
        }
    }

    /// Every positive instance weight is resolvable by a 53-bit uniform. For
    /// pair trees the smallest factor bounds every h·n from below.
    fn resolved(&self) -> bool {
        match (&self.single, &self.pair) {
            (Some(t), _) => t.min_positive() >= t.total() * RESOLUTION,
            (None, Some(t)) => t.min_factor() >= t.total_weight() * RESOLUTION,
            (None, None) => unreachable!("class without tree"),
        }
    }
}

struct ObservableState {
    /// Weight per slot; NaN when the instance does not match.
    w: Vec<f64>,
    sum: Option<SumTree>,
    extreme: Option<MaxTree>,
}

struct State<'m> {
    model: &'m Compiled,
    types: Vec<TypeState>,
    classes: Vec<ClassState>,
    observables: Vec<ObservableState>,
    values: Vec<f64>,
    live: u64,
    peak: u64,
    created: u64,
    max_instances: u64,
}

fn grow(v: &mut Vec<f64>, slot: usize, fill: f64) {
    if v.len() <= slot {
        v.resize(slot + 1, fill);
    }
}

/// A `where` condition; an undefined (NaN) value is an error, not true.
fn condition(
    expr: Option<&crate::expr::Expr>,
    env: &Env,
    path: &str,
    attrs: &[f64],
) -> Result<bool, Halt> {
    let Some(e) = expr else { return Ok(true) };
    let value = e.eval(env);
    if value.is_nan() {
        return Err(numerical(
            path,
            format!("condition is undefined (NaN) for attributes {attrs:?}"),
            "Guard divisions and logarithms so conditions are defined for every instance.",
        ));
    }
    Ok(value != 0.0)
}

impl<'m> State<'m> {
    fn new(model: &'m Compiled, max_instances: u64) -> Self {
        Self {
            model,
            types: model
                .types
                .iter()
                .map(|t| TypeState {
                    stride: t.attributes.len(),
                    attrs: Vec::new(),
                    counts: Vec::new(),
                    free: Vec::new(),
                    index: HashMap::new(),
                })
                .collect(),
            classes: model
                .classes
                .iter()
                .map(|c| ClassState {
                    h: Vec::new(),
                    single: (!c.squares).then(|| SumTree::new(true)),
                    pair: c.squares.then(PairTree::new),
                })
                .collect(),
            observables: model
                .observables
                .iter()
                .map(|o| ObservableState {
                    w: Vec::new(),
                    sum: matches!(o.reduce, Reduce::Sum | Reduce::Distinct)
                        .then(|| SumTree::new(false)),
                    extreme: matches!(o.reduce, Reduce::Max | Reduce::Min).then(MaxTree::new),
                })
                .collect(),
            values: vec![0.0; model.observables.len()],
            live: 0,
            peak: 0,
            created: 0,
            max_instances,
        }
    }

    fn attributes(&self, t: usize, slot: usize) -> &[f64] {
        let ts = &self.types[t];
        &ts.attrs[slot * ts.stride..(slot + 1) * ts.stride]
    }

    /// Local factors and observable weights of a new instance.
    fn evaluate(&mut self, t: usize, slot: usize) -> Result<(), Halt> {
        let model = self.model;
        let ts = &self.types[t];
        let attrs = &ts.attrs[slot * ts.stride..(slot + 1) * ts.stride];
        let env = Env {
            params: &model.params,
            observables: &[],
            attributes: [attrs, &[]],
        };
        for &c in &model.types[t].classes {
            let info = &model.classes[c];
            let matched = condition(
                info.condition.as_ref(),
                &env,
                &format!("{}/where", info.pointer),
                attrs,
            )?;
            let h = if matched {
                info.factor.as_ref().map_or(1.0, |e| e.eval(&env))
            } else {
                0.0
            };
            if !(h.is_finite() && h >= 0.0) {
                return Err(numerical(
                    &format!("{}/factor", info.pointer),
                    format!("factor is {h} for attributes {attrs:?}"),
                    "Local factors must be finite and nonnegative for every instance the pattern matches.",
                ));
            }
            let state = &mut self.classes[c];
            grow(&mut state.h, slot, 0.0);
            state.h[slot] = h;
        }
        for &o in &model.types[t].observables {
            let info = &model.observables[o];
            let at = crn_craft::diagnostic::pointer("/observables", &info.name);
            let matched = condition(info.condition.as_ref(), &env, &format!("{at}/where"), attrs)?;
            let w = if matched {
                info.weight.as_ref().map_or(1.0, |e| e.eval(&env))
            } else {
                f64::NAN
            };
            if matched && !w.is_finite() {
                return Err(numerical(
                    &crn_craft::diagnostic::pointer("/observables", &info.name),
                    format!("weight is {w} for attributes {attrs:?}"),
                    "Observable weights must be finite.",
                ));
            }
            let state = &mut self.observables[o];
            grow(&mut state.w, slot, f64::NAN);
            state.w[slot] = w;
        }
        Ok(())
    }

    #[allow(clippy::cast_precision_loss)]
    fn set_leaves(&mut self, t: usize, slot: usize) -> Result<(), Halt> {
        let model = self.model;
        let count = self.types[t].counts[slot];
        let n = count as f64;
        for &c in &model.types[t].classes {
            let state = &mut self.classes[c];
            let h = state.h[slot];
            let weight = h * n;
            if !weight.is_finite() {
                return Err(numerical(
                    &model.classes[c].pointer,
                    "pattern weight overflows",
                    "Reduce the factor or the population.",
                ));
            }
            if let Some(single) = state.single.as_mut() {
                single.set(slot, weight);
            }
            if let Some(pair) = state.pair.as_mut() {
                pair.set(slot, h, count);
            }
        }
        for &o in &model.types[t].observables {
            let state = &mut self.observables[o];
            let w = state.w[slot];
            let matched = !w.is_nan() && n > 0.0;
            if let Some(sum) = state.sum.as_mut() {
                let leaf = match (matched, model.observables[o].reduce) {
                    (false, _) => 0.0,
                    (true, Reduce::Distinct) => w,
                    (true, _) => w * n,
                };
                if !leaf.is_finite() {
                    return Err(numerical(
                        &crn_craft::diagnostic::pointer("/observables", &model.observables[o].name),
                        "observable weight overflows",
                        "Reduce the weight or the population.",
                    ));
                }
                sum.set(slot, leaf);
            }
            if let Some(extreme) = state.extreme.as_mut() {
                let signed = if model.observables[o].reduce == Reduce::Min {
                    -w
                } else {
                    w
                };
                extreme.set(slot, matched.then_some(signed));
            }
        }
        Ok(())
    }

    fn add(&mut self, t: usize, attrs: &[f64], k: u64) -> Result<(), Halt> {
        let key: Box<[u64]> = attrs.iter().map(|x| x.to_bits()).collect();
        if let Some(&slot) = self.types[t].index.get(&key) {
            let count = &mut self.types[t].counts[slot];
            *count += k;
            if *count > MAX_COUNT {
                return Err(Halt::Limit(
                    Termination::InstanceLimit,
                    Diagnostic::error(
                        "CRN-R042",
                        "/budget/max_instances",
                        "an instance count exceeds 2^53",
                        "Counts must stay exact integers; reduce the population.",
                    ),
                ));
            }
            return self.set_leaves(t, slot);
        }
        self.live += 1;
        self.peak = self.peak.max(self.live);
        if self.live > self.max_instances {
            return Err(Halt::Limit(
                Termination::InstanceLimit,
                Diagnostic::error(
                    "CRN-R042",
                    "/budget/max_instances",
                    "instance limit reached",
                    "Increase max_instances; incomplete ensembles withhold statistics.",
                ),
            ));
        }
        let ts = &mut self.types[t];
        let slot = if let Some(slot) = ts.free.pop() {
            ts.attrs[slot * ts.stride..(slot + 1) * ts.stride].copy_from_slice(attrs);
            ts.counts[slot] = k;
            slot
        } else {
            ts.attrs.extend_from_slice(attrs);
            ts.counts.push(k);
            ts.counts.len() - 1
        };
        ts.index.insert(key, slot);
        self.evaluate(t, slot)?;
        self.set_leaves(t, slot)
    }

    fn remove(&mut self, t: usize, slot: usize, k: u64) -> Result<(), Halt> {
        let ts = &mut self.types[t];
        ts.counts[slot] -= k;
        if ts.counts[slot] == 0 {
            let key: Box<[u64]> = ts.attrs[slot * ts.stride..(slot + 1) * ts.stride]
                .iter()
                .map(|x| x.to_bits())
                .collect();
            ts.index.remove(&key);
            ts.free.push(slot);
            self.live -= 1;
        }
        self.set_leaves(t, slot)
    }

    fn refresh(&mut self) -> Result<(), Halt> {
        for (o, state) in self.observables.iter().enumerate() {
            let value = if let Some(sum) = &state.sum {
                sum.total()
            } else {
                let top = state.extreme.as_ref().expect("extreme tree").top();
                if top == f64::NEG_INFINITY {
                    0.0
                } else if self.model.observables[o].reduce == Reduce::Min {
                    -top
                } else {
                    top
                }
            };
            if !value.is_finite() {
                return Err(numerical(
                    &crn_craft::diagnostic::pointer(
                        "/observables",
                        &self.model.observables[o].name,
                    ),
                    "observable total overflows",
                    "Rescale the weights; observables must stay finite.",
                ));
            }
            self.values[o] = value;
        }
        Ok(())
    }

    fn global_env(&self) -> Env<'_> {
        Env {
            params: &self.model.params,
            observables: &self.values,
            attributes: [&[], &[]],
        }
    }

    /// Rule propensities in rule order and their sequential sum.
    fn propensities(&self, a: &mut [f64]) -> Result<f64, Halt> {
        let env = self.global_env();
        let mut total = 0.0;
        for (r, rule) in self.model.rules.iter().enumerate() {
            let weight = match rule.shape {
                Shape::Source => 1.0,
                Shape::Single(c) => self.classes[c].weight(),
                Shape::Pair(c) => self.classes[c]
                    .pair
                    .as_ref()
                    .expect("pair class keeps a pair tree")
                    .total_pairs(),
                Shape::Cross(c, d) => self.classes[c].weight() * self.classes[d].weight(),
            };
            if !weight.is_finite() {
                return Err(numerical(
                    &format!("/rules/{r}"),
                    format!("reactant weight of `{}` overflows", rule.name),
                    "Reduce factors or populations.",
                ));
            }
            a[r] = if weight > 0.0 {
                let rate = rule.rate.eval(&env);
                if !(rate.is_finite() && rate >= 0.0) {
                    return Err(numerical(
                        &format!("/rules/{r}/rate"),
                        format!(
                            "rate of `{}` is {rate} at observables {:?}",
                            rule.name, self.values
                        ),
                        "Rates must be finite and nonnegative wherever the rule has reactants.",
                    ));
                }
                rate * weight
            } else {
                0.0
            };
            if !a[r].is_finite() {
                return Err(numerical(
                    &format!("/rules/{r}"),
                    format!("propensity of `{}` overflows", rule.name),
                    "Reduce rates or populations.",
                ));
            }
            total += a[r];
        }
        if !total.is_finite() {
            return Err(numerical(
                "/rules",
                "total propensity overflows",
                "Reduce rates or populations.",
            ));
        }
        // Every positive weight must be resolvable by a 53-bit uniform.
        for (r, rule) in self.model.rules.iter().enumerate() {
            if a[r] <= 0.0 {
                continue;
            }
            let unresolved = a[r] < total * RESOLUTION
                || match rule.shape {
                    Shape::Source => false,
                    Shape::Single(c) => !self.classes[c].resolved(),
                    Shape::Cross(c, d) => {
                        !self.classes[c].resolved() || !self.classes[d].resolved()
                    }
                    Shape::Pair(c) => {
                        // Every pair weight h(x)h(y) is at least the smallest
                        // factor squared.
                        let tree = self.classes[c].pair.as_ref().expect("pair tree");
                        let m = tree.min_factor();
                        !self.classes[c].resolved() || m * m < tree.total_pairs() * RESOLUTION
                    }
                };
            if unresolved {
                return Err(Halt::Numerical(Diagnostic::error(
                    "CRN-R031",
                    &format!("/rules/{r}"),
                    format!(
                        "a positive weight of `{}` is below 2^-53 of its total",
                        rule.name
                    ),
                    "Rescale rates or factors; weights this small cannot be selected exactly with 53-bit uniforms.",
                )));
            }
        }
        Ok(total)
    }

    fn instances(&self) -> Vec<InstanceRecord> {
        let mut out = Vec::new();
        for (t, ts) in self.types.iter().enumerate() {
            for (slot, &count) in ts.counts.iter().enumerate() {
                if count > 0 {
                    out.push(InstanceRecord {
                        type_index: t,
                        attributes: self.attributes(t, slot).to_vec(),
                        count,
                    });
                }
            }
        }
        out.sort_by(|a, b| {
            a.type_index.cmp(&b.type_index).then_with(|| {
                a.attributes
                    .iter()
                    .zip(&b.attributes)
                    .map(|(x, y)| x.total_cmp(y))
                    .find(|o| o.is_ne())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        });
        out
    }
}

pub(crate) struct Words<'r, R: RngCore> {
    rng: &'r mut R,
    used: u64,
    max: u64,
}
#[cfg(test)]
impl<'r, R: RngCore> Words<'r, R> {
    pub(crate) fn for_test(rng: &'r mut R, max: u64) -> Self {
        Self { rng, used: 0, max }
    }
}

impl<R: RngCore> Words<'_, R> {
    pub(crate) fn open(&mut self) -> Result<f64, Halt> {
        loop {
            let u = self.uniform()?;
            if u > 0.0 {
                return Ok(u);
            }
        }
    }

    pub(crate) fn uniform(&mut self) -> Result<f64, Halt> {
        Ok(uniform_closed_open(self.word()?))
    }

    pub(crate) fn word(&mut self) -> Result<u64, Halt> {
        if self.used >= self.max {
            return Err(Halt::Limit(
                Termination::SamplingLimit,
                Diagnostic::error(
                    "CRN-R041",
                    "/budget/max_sampling_steps_per_replicate",
                    "sampling work limit reached",
                    "Increase the sampling-step budget; incomplete ensembles withhold statistics.",
                ),
            ));
        }
        self.used += 1;
        Ok(self.rng.next_u64())
    }
}

pub(crate) struct Path<'a> {
    pub model: &'a Compiled,
    pub report_times: &'a [f64],
    pub end: f64,
    pub limits: &'a Limits,
    pub control: &'a Control,
    pub keep_instances: bool,
    pub observations: &'a Spec,
    pub id: u64,
    pub flat: Option<&'a crate::flat::Network>,
    pub generated: Option<&'a dyn crate::RulesKernel>,
    pub generated_loop: Option<&'a dyn crate::RulesLoop>,
    pub mtdna: Option<&'a crate::mtdna::Kernel>,
    pub method: crate::Method,
}

impl Path<'_> {
    pub(crate) fn run<R: RngCore + 'static, O: EventObserver>(
        &self,
        rng: &mut R,
        observer: &O,
    ) -> PathOutput {
        let model = self.model;
        let n_obs = model.observables.len();
        let mut out = PathOutput {
            termination: Termination::Completed,
            events: 0,
            stop_time: 0.0,
            rule_events: vec![0; model.rules.len()],
            rng_words: 0,
            peak_instances: 0,
            created_instances: 0,
            saved_points: 0,
            values: vec![0.0; self.report_times.len() * n_obs],
            diagnostic: None,
            observation: None,
            instances: None,
        };
        let mut words = Words {
            rng,
            used: 0,
            max: self.limits.max_words,
        };
        let mut time = 0.0;
        let mut monitor = self
            .observations
            .enabled()
            .then(|| Monitor::new(self.observations, self.id));
        let mut values = vec![0.0; n_obs];
        let result = if let Some(network) = self.flat {
            let mut counts = network.initial();
            let result = self.simulate_flat(
                network,
                &mut counts,
                &mut values,
                &mut words,
                &mut out,
                &mut time,
                &mut monitor,
                observer,
            );
            if self.keep_instances {
                out.instances = Some(network.instances(&counts));
            }
            result
        } else {
            let mut state = State::new(model, self.limits.max_instances);
            let result = self.simulate(
                &mut state,
                &mut words,
                &mut out,
                &mut time,
                &mut monitor,
                observer,
            );
            out.peak_instances = state.peak;
            out.created_instances = state.created;
            if self.keep_instances {
                out.instances = Some(state.instances());
            }
            values = state.values;
            result
        };
        out.rng_words = words.used;
        match result {
            Ok(stopped) => {
                // The state is held to the end of the grid.
                while out.saved_points < self.report_times.len() {
                    let row = out.saved_points * n_obs;
                    out.values[row..row + n_obs].copy_from_slice(&values);
                    out.saved_points += 1;
                }
                if let Some(at) = stopped {
                    out.termination = Termination::Stopped;
                    out.stop_time = at;
                } else {
                    out.stop_time = self.end;
                }
            }
            Err(Halt::Limit(termination, diagnostic)) => {
                out.termination = termination;
                out.diagnostic = Some(diagnostic);
                out.stop_time = time;
            }
            Err(Halt::Numerical(diagnostic)) => {
                out.termination = Termination::Failed;
                out.diagnostic = Some(diagnostic);
                out.stop_time = time;
            }
        }
        if let Some(observation) = monitor.and_then(|m| m.finish(out.stop_time)) {
            if observation
                .trackers
                .iter()
                .any(|t| !t.integral.is_finite() || !t.omitted_duration.is_finite())
            {
                out.termination = Termination::Failed;
                out.diagnostic = Some(Diagnostic::error(
                    "CRN-R030",
                    "/trackers",
                    "tracker integral or duration exceeds numerical range",
                    "Rescale the observable or the time horizon.",
                ));
            } else {
                out.observation = Some(observation);
            }
        }
        out
    }

    fn stop_now(&self, state: &State) -> Result<bool, Halt> {
        self.stop_values(&state.values)
    }

    fn stop_values(&self, values: &[f64]) -> Result<bool, Halt> {
        let Some(e) = self.model.stop_program.as_ref() else {
            return Ok(false);
        };
        let value = e.eval(&Env {
            params: &self.model.params,
            observables: values,
            attributes: [&[], &[]],
        });
        if value.is_nan() {
            return Err(numerical(
                "/stop",
                format!("stop condition is undefined (NaN) at observables {values:?}"),
                "Guard divisions so the stop condition is defined in every state.",
            ));
        }
        Ok(value != 0.0)
    }

    fn poll(&self) -> Result<(), Halt> {
        if self.control.is_cancelled() {
            return Err(Halt::Limit(
                Termination::Cancelled,
                Diagnostic::error(
                    "CRN-R044",
                    "",
                    "cancelled",
                    "Rerun to obtain a complete ensemble.",
                ),
            ));
        }
        if self.limits.deadline.is_some_and(|d| Instant::now() >= d) {
            return Err(Halt::Limit(
                Termination::WallTimeLimit,
                Diagnostic::error(
                    "CRN-R043",
                    "/budget/wall_time_s",
                    "wall-time limit reached",
                    "Increase wall_time_s; incomplete ensembles withhold statistics.",
                ),
            ));
        }
        Ok(())
    }

    /// Returns the stop time if the stop condition fired.
    #[allow(clippy::too_many_lines, clippy::cast_precision_loss)]
    fn simulate<R: RngCore, O: EventObserver>(
        &self,
        state: &mut State,
        words: &mut Words<R>,
        out: &mut PathOutput,
        time: &mut f64,
        monitor: &mut Option<Monitor<'_>>,
        observer: &O,
    ) -> Result<Option<f64>, Halt> {
        let model = self.model;
        let n_obs = model.observables.len();
        self.poll()?;
        for (t, attrs, count) in &model.initial {
            state.add(*t, attrs, *count)?;
        }
        state.refresh()?;
        if let Some(m) = monitor {
            m.start(&model.params, &state.values);
        }
        if self.stop_now(state)? {
            return Ok(Some(0.0));
        }
        let mut a = vec![0.0; model.rules.len()];
        let mut reactants: [Vec<f64>; 2] = [Vec::new(), Vec::new()];
        let mut products: Vec<(usize, Vec<f64>)> = Vec::new();
        let mut next_id = model.first_id;
        loop {
            if out.events % POLL_EVENTS == 0 && out.events > 0 {
                self.poll()?;
            }
            let total = state.propensities(&mut a)?;
            if total <= 0.0 {
                return Ok(None);
            }
            let u = 1.0 - words.uniform()?;
            let next = *time + (-libm::log(u) / total);
            while out.saved_points < self.report_times.len()
                && self.report_times[out.saved_points] < next
            {
                let row = out.saved_points * n_obs;
                out.values[row..row + n_obs].copy_from_slice(&state.values);
                out.saved_points += 1;
            }
            if next > self.end {
                return Ok(None);
            }
            if out.events >= self.limits.max_events {
                return Err(Halt::Limit(
                    Termination::EventLimit,
                    Diagnostic::error(
                        "CRN-R040",
                        "/budget/max_events_per_replicate",
                        "event limit reached",
                        "Increase max_events_per_replicate; incomplete ensembles withhold statistics.",
                    ),
                ));
            }
            // Rule.
            let target = words.uniform()? * total;
            let mut cumulative = 0.0;
            let mut chosen = None;
            for (r, &ar) in a.iter().enumerate() {
                if ar > 0.0 {
                    chosen = Some(r);
                    cumulative += ar;
                    if target < cumulative {
                        break;
                    }
                }
            }
            let r = chosen.expect("positive total has a positive rule");
            let rule = &model.rules[r];
            // Reactant instances.
            let mut slots = [usize::MAX; 2];
            match rule.shape {
                Shape::Source => {}
                Shape::Single(c) => {
                    let class = &state.classes[c];
                    slots[0] = class.select(words.uniform()? * class.weight());
                }
                Shape::Pair(c) => {
                    // Exact: a node's pair weight splits into pairs inside its
                    // children and pairs across them; no rejection.
                    let tree = state.classes[c].pair.as_ref().expect("pair tree");
                    match tree.select(words.uniform()? * tree.total_pairs()) {
                        PairChoice::Same(i) => slots = [i, i],
                        PairChoice::Across(l, r) => {
                            slots[0] = tree.descend(l, words.uniform()? * tree.node_weight(l));
                            slots[1] = tree.descend(r, words.uniform()? * tree.node_weight(r));
                        }
                    }
                }
                Shape::Cross(c, d) => {
                    let first = &state.classes[c];
                    slots[0] = first.select(words.uniform()? * first.weight());
                    let second = &state.classes[d];
                    slots[1] = second.select(words.uniform()? * second.weight());
                }
            }
            for (k, &t) in rule.reactant_types.iter().enumerate() {
                reactants[k].clear();
                reactants[k].extend_from_slice(state.attributes(t, slots[k]));
            }
            // Products are evaluated before the reactants are consumed.
            products.truncate(rule.products.len());
            while products.len() < rule.products.len() {
                products.push((0, Vec::new()));
            }
            {
                let env = Env {
                    params: &model.params,
                    observables: &[],
                    attributes: [&reactants[0], &reactants[1]],
                };
                let mut id_overflow = false;
                let mut draw = Draw {
                    rng: &mut *words.rng,
                    next_id: &mut next_id,
                    words: &mut words.used,
                    id_overflow: &mut id_overflow,
                };
                for (j, product) in rule.products.iter().enumerate() {
                    let (t, values) = &mut products[j];
                    *t = product.type_id;
                    values.clear();
                    for (index, source) in product.attributes.iter().enumerate() {
                        let x = match source {
                            Source::Copy { reactant, index } => reactants[*reactant][*index],
                            Source::Expression(e) => e.eval_random(&env, &mut draw),
                        };
                        if *draw.id_overflow {
                            return Err(numerical(
                                &format!("/rules/{r}/products/{j}"),
                                "new_id() reached 2^53",
                                "Identities must stay exact integers; start initial identities lower.",
                            ));
                        }
                        if !x.is_finite() {
                            return Err(numerical(
                                &format!(
                                    "/rules/{r}/products/{j}/attributes/{}",
                                    model.types[product.type_id].attributes[index]
                                ),
                                format!("product attribute is {x}"),
                                "Product attributes must be finite.",
                            ));
                        }
                        values.push(if x == 0.0 { 0.0 } else { x });
                    }
                }
                if *draw.words > words.max {
                    return Err(Halt::Limit(
                        Termination::SamplingLimit,
                        Diagnostic::error(
                            "CRN-R041",
                            "/budget/max_sampling_steps_per_replicate",
                            "sampling work limit reached",
                            "Increase the sampling-step budget; incomplete ensembles withhold statistics.",
                        ),
                    ));
                }
            }
            // A product equal to a reactant cancels its consumption, so a kept
            // reactant (`copy` without overrides) is never removed and re-created.
            let mut kept = [false; 2];
            let mut skip = [false; crate::model::MAX_PRODUCTS];
            for (j, (t, values)) in products.iter().enumerate() {
                for (k, &rt) in rule.reactant_types.iter().enumerate() {
                    if !kept[k]
                        && rt == *t
                        && reactants[k]
                            .iter()
                            .zip(values)
                            .all(|(x, y)| x.to_bits() == y.to_bits())
                    {
                        kept[k] = true;
                        skip[j] = true;
                        break;
                    }
                }
            }
            for (k, &rt) in rule.reactant_types.iter().enumerate() {
                if !kept[k] {
                    state.remove(rt, slots[k], 1)?;
                }
            }
            for (j, (t, values)) in products.iter().enumerate() {
                if skip[j] {
                    continue;
                }
                let before = state.live;
                state.add(*t, values, 1)?;
                if state.live > before {
                    state.created += 1;
                }
            }
            *time = next;
            out.events += 1;
            out.rule_events[r] += 1;
            state.refresh()?;
            let stopped = if O::ENABLED || monitor.is_some() {
                let participants =
                    O::ENABLED || monitor.as_ref().is_some_and(|m| m.keeps_exact_event(r));
                let (rs, ps) = if participants {
                    let rs = rule
                        .reactant_types
                        .iter()
                        .copied()
                        .zip(&reactants)
                        .map(|(t, a)| (t, a.as_slice()))
                        .collect();
                    let ps = products.iter().map(|(t, a)| (*t, a.as_slice())).collect();
                    (rs, ps)
                } else {
                    (Vec::new(), Vec::new())
                };
                let event = Event {
                    replicate: self.id,
                    time: next,
                    rule: Some(r),
                    reactants: &rs,
                    products: &ps,
                    observables: &state.values,
                    firings: &[],
                };
                if let Some(m) = monitor {
                    m.observe(&event, &model.params);
                }
                O::ENABLED && observer.on_event(&event) == ObserverAction::Stop
            } else {
                false
            };
            if stopped || self.stop_now(state)? {
                return Ok(Some(next));
            }
        }
    }
}
