//! Propensity bounds of RSSA-CR over population intervals
//! (Thanh, Priami & Zunino, J. Chem. Phys. 141, 134116 (2014)), for the
//! lowered channels of the rules engine.
//!
//! Ported from the 0.6.0 network engine (`ssa/bounds.rs`). Every species has
//! an interval around its count; a channel's bounds enclose the value the
//! engine computes for it ([`Network::value_with`] with the global factor of
//! [`Network::global_fresh`]) at every integer state of the box. They are
//! recomputed only for the readers of a species whose count leaves its
//! interval.
//!
//! - The count factors are monotone in each count. For channels of at most
//!   two molecules their arithmetic is a product of correctly rounded
//!   operations on nonnegative values, which is monotone itself, so the
//!   corners are exact bounds. A scaled combination count
//!   ([`crate::rules::combination_rate`]) is within (2m + 1) roundings of
//!   the exact count for m ratios; its corners are widened by twice that.
//! - A rule whose rate reads one species is evaluated at every count of that
//!   species' interval (at most 129), the engine's own arithmetic, so its
//!   range over the box is exact, whatever the expression.
//! - A rule whose rate reads two or more species pins them: their intervals
//!   are their counts, and the bounds are the value at the current state.
//! - A box where a value is not finite, negative, or too small for a 53-bit
//!   acceptance draw pins the channel's reads for the rest of the path
//!   (the network engine's point fallback).
use crate::rules::engine::{Halt, RESOLUTION};
use crate::rules::expr::{Env, RateKernel};
use crate::rules::flat::{Network, RateCache};
use crate::rules::select::{Composition, RateTree};

/// Per-network tables, built once per prepared run.
pub(crate) struct Plan {
    /// Reactant molecules each channel needs: `(species, count)`.
    needs: Vec<Vec<(usize, u64)>>,
    /// The single species a rule's rate reads, if it reads exactly one.
    single: Vec<Option<usize>>,
    /// Species read by rules that read two or more species.
    pinned: Vec<bool>,
    /// For each channel, the species its value reads: reactants and the
    /// rule's reads, sorted.
    reads: Vec<Vec<usize>>,
    /// The range cache of each rule's global factor. Rules whose rates
    /// compile to the same program (and so read the same projections and
    /// compute the same value in every state) share one; every other rule
    /// has its own.
    range_slot: Vec<usize>,
    range_slots: usize,
}

impl Plan {
    pub(crate) fn new(net: &Network) -> Self {
        let mut pinned = vec![false; net.species.len()];
        for reads in &net.rule_reads {
            if reads.len() > 1 {
                for &s in reads {
                    pinned[s] = true;
                }
            }
        }
        let needs: Vec<_> = (0..net.channels.len()).map(|i| net.needs(i)).collect();
        let reads = net
            .channels
            .iter()
            .zip(&needs)
            .map(|(c, needs)| {
                let mut reads: Vec<usize> = needs
                    .iter()
                    .map(|&(s, _)| s)
                    .chain(net.rule_reads[c.rule].iter().copied())
                    .collect();
                reads.sort_unstable();
                reads.dedup();
                reads
            })
            .collect();
        let mut shared = std::collections::HashMap::new();
        let mut range_slots = 0;
        let range_slot = net
            .rates
            .iter()
            .map(|kernel| {
                let mut fresh = || {
                    range_slots += 1;
                    range_slots - 1
                };
                match kernel {
                    RateKernel::Expression { program, .. } => match program.key() {
                        Some(key) => *shared.entry(key).or_insert_with(fresh),
                        None => fresh(),
                    },
                    _ => fresh(),
                }
            })
            .collect();
        Self {
            needs,
            single: net
                .rule_reads
                .iter()
                .map(|r| if r.len() == 1 { Some(r[0]) } else { None })
                .collect(),
            pinned,
            reads,
            range_slot,
            range_slots,
        }
    }
}

/// The interval of a species with `count` molecules: a point when pinned or
/// empty, otherwise ± clamp(count / 10, 4, 64).
fn interval(count: u64, pinned: bool) -> (u64, u64) {
    if pinned || count == 0 {
        return (count, count);
    }
    let width = (count / 10).clamp(4, 64);
    (count.saturating_sub(width), count.saturating_add(width))
}

/// A rule's global factor over one species' counts `start..=end`.
#[derive(Clone, Copy)]
struct Range {
    start: u64,
    end: u64,
    minimum: f64,
    maximum: f64,
    positive: f64,
}

/// Lower, upper and smallest positive value of a channel over the box.
type Bound = (f64, f64, f64);

const UNIT_ROUNDOFF: f64 = 1.0 / 9_007_199_254_740_992.0;

/// The next float above a nonnegative `x` (infinity and NaN unchanged).
fn next_up(x: f64) -> f64 {
    if !x.is_finite() {
        x
    } else if x <= 0.0 {
        f64::from_bits(1)
    } else {
        f64::from_bits(x.to_bits() + 1)
    }
}

/// The next float below a positive finite `x`.
fn next_down(x: f64) -> f64 {
    f64::from_bits(x.to_bits() - 1)
}

/// Per-path state.
#[derive(Default)]
pub(crate) struct Bounds {
    low: Vec<u64>,
    high: Vec<u64>,
    pinned: Vec<bool>,
    /// Lower bound of each channel.
    pub(crate) lower: Vec<f64>,
    /// Channels with a positive lower bound.
    pub(crate) positive_lower: usize,
    marks: Vec<u32>,
    generation: u32,
    dirty: Vec<usize>,
    ranges: Vec<Option<Range>>,
    saved: Vec<f64>,
    /// The current counts of a channel's reactants while its corners are
    /// evaluated.
    corner_saved: Vec<u64>,
}

impl Bounds {
    pub(crate) fn reset(&mut self, net: &Network, plan: &Plan, counts: &[u64]) {
        self.pinned.clone_from(&plan.pinned);
        self.low.clear();
        self.high.clear();
        for (s, &count) in counts.iter().enumerate() {
            let (low, high) = interval(count, self.pinned[s]);
            self.low.push(low);
            self.high.push(high);
        }
        let n = net.channels.len();
        self.lower.clear();
        self.lower.resize(n, 0.0);
        self.positive_lower = 0;
        self.marks.clear();
        self.marks.resize(n, 0);
        self.generation = 0;
        self.dirty.clear();
        self.ranges.clear();
        self.ranges.resize(plan.range_slots, None);
    }

    /// The range of channel `i`'s global factor over the box, or `None` if
    /// some value there is not finite or negative.
    fn global_range(
        &mut self,
        net: &Network,
        plan: &Plan,
        i: usize,
        counts: &mut [u64],
        values: &mut [f64],
        params: &[f64],
    ) -> Option<(f64, f64, f64)> {
        if let Some(k) = net.constant_global(i) {
            let positive = if k > 0.0 { k } else { f64::INFINITY };
            return (k.is_finite() && k >= 0.0).then_some((k, k, positive));
        }
        let rule = net.channels[i].rule;
        let Some(s) = plan.single[rule].filter(|&s| self.low[s] < self.high[s]) else {
            // Every species the rate reads is a point: its current value.
            let env = Env {
                params,
                observables: values,
                attributes: [&[], &[]],
            };
            let g = net.global_fresh(i, &env);
            let positive = if g > 0.0 { g } else { f64::INFINITY };
            return (g.is_finite() && g >= 0.0).then_some((g, g, positive));
        };
        // Counts below the channel's own need of `s` give it no value.
        let need = plan.needs[i]
            .iter()
            .find(|&&(t, _)| t == s)
            .map_or(0, |&(_, k)| k);
        let start = self.low[s].max(need);
        let end = self.high[s];
        let slot = plan.range_slot[rule];
        if let Some(range) = self.ranges[slot] {
            if range.start == start && range.end == end {
                return Some((range.minimum, range.maximum, range.positive));
            }
        }
        let projections = &net.rule_projections[rule];
        self.saved.clear();
        self.saved
            .extend(projections.iter().map(|&p| values[net.observables[p].0]));
        let current = counts[s];
        let mut range = Range {
            start,
            end,
            minimum: f64::INFINITY,
            maximum: 0.0,
            positive: f64::INFINITY,
        };
        let mut valid = true;
        for v in start..=end {
            counts[s] = v;
            for &p in projections {
                let x = net.project(p, counts);
                valid &= x.is_finite();
                values[net.observables[p].0] = x;
            }
            if !valid {
                break;
            }
            let env = Env {
                params,
                observables: values,
                attributes: [&[], &[]],
            };
            let g = net.global_fresh(i, &env);
            if !(g.is_finite() && g >= 0.0) {
                valid = false;
                break;
            }
            range.minimum = range.minimum.min(g);
            range.maximum = range.maximum.max(g);
            if g > 0.0 {
                range.positive = range.positive.min(g);
            }
        }
        counts[s] = current;
        for (&p, &x) in projections.iter().zip(&self.saved) {
            values[net.observables[p].0] = x;
        }
        if !valid {
            return None;
        }
        self.ranges[slot] = Some(range);
        Some((range.minimum, range.maximum, range.positive))
    }

    /// Bounds of channel `i` over the box, or `None` when it must be pinned.
    #[allow(clippy::too_many_arguments, clippy::cast_precision_loss)]
    fn bound(
        &mut self,
        net: &Network,
        plan: &Plan,
        i: usize,
        counts: &mut [u64],
        values: &mut [f64],
        params: &[f64],
    ) -> Option<Bound> {
        let needs = &plan.needs[i];
        if needs.iter().any(|&(s, k)| self.high[s] < k) {
            return Some((0.0, 0.0, f64::INFINITY));
        }
        let gated = needs.iter().any(|&(s, k)| self.low[s] < k);
        let (g_low, g_high, g_positive) =
            self.global_range(net, plan, i, counts, values, params)?;
        if g_high == 0.0 {
            return Some((0.0, 0.0, f64::INFINITY));
        }
        let mut saved = std::mem::take(&mut self.corner_saved);
        saved.clear();
        saved.extend(needs.iter().map(|&(s, _)| counts[s]));
        let (low, high) = (&self.low, &self.high);
        let at = |global: f64, corner: &dyn Fn(usize, u64) -> u64, counts: &mut [u64]| {
            for &(s, k) in needs {
                counts[s] = corner(s, k);
            }
            net.value_with(i, global, counts)
        };
        let upper = at(g_high, &|s, _| high[s], counts);
        let positive = if g_positive.is_finite() {
            at(g_positive, &|s, k| low[s].max(k), counts)
        } else {
            f64::INFINITY
        };
        let lower = if gated {
            0.0
        } else {
            at(g_low, &|s, _| low[s], counts)
        };
        let ratios: u64 = needs
            .iter()
            .map(|&(s, k)| k.min(high[s].saturating_sub(k)))
            .sum();
        for (&(s, _), &count) in needs.iter().zip(&saved) {
            counts[s] = count;
        }
        self.corner_saved = saved;
        let (lower, upper, positive) = if net.scaled_combination(i) {
            // Within (2m + 1) roundings of the exact, monotone count times
            // the global factor; widen by twice that, rounded outward.
            let slack = 4.0 * (2.0 * ratios as f64 + 2.0) * UNIT_ROUNDOFF;
            let down = |x: f64| {
                if x < f64::MIN_POSITIVE * 1e20 {
                    0.0
                } else {
                    next_down(x * (1.0 - slack))
                }
            };
            (
                down(lower),
                next_up(upper * (1.0 + slack)) + f64::from_bits(2),
                if positive.is_finite() {
                    down(positive)
                } else {
                    positive
                },
            )
        } else {
            (lower, upper, positive)
        };
        let usable = upper.is_finite()
            && upper <= f64::MAX / 4.0
            && lower >= 0.0
            && lower <= upper
            && (upper == 0.0 || positive >= upper * RESOLUTION);
        usable.then_some((lower, upper, positive))
    }

    /// Pins channel `i`'s reads to their counts and returns its exact value
    /// at the current state as a point bound.
    #[allow(clippy::too_many_arguments)]
    fn point(
        &mut self,
        net: &Network,
        plan: &Plan,
        i: usize,
        counts: &[u64],
        values: &[f64],
        params: &[f64],
        cache: &mut RateCache,
    ) -> Result<Bound, Halt> {
        for &s in &plan.reads[i] {
            self.pinned[s] = true;
            self.low[s] = counts[s];
            self.high[s] = counts[s];
        }
        let env = Env {
            params,
            observables: values,
            attributes: [&[], &[]],
        };
        let a = cache.channel(net, i, counts, &env)?;
        Ok((a, a, if a > 0.0 { a } else { f64::INFINITY }))
    }

    /// After the first call (`fired` is `None`, every channel) or an event of
    /// channel `fired`: new intervals for the species that left theirs, and
    /// new bounds for their readers, in channel order, staged in `tree` and
    /// `composition`. `cache` must be reset for the current state.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn refresh(
        &mut self,
        net: &Network,
        plan: &Plan,
        fired: Option<usize>,
        counts: &mut [u64],
        values: &mut [f64],
        params: &[f64],
        cache: &mut RateCache,
        tree: &mut RateTree,
        composition: &mut Composition,
    ) -> Result<(), Halt> {
        self.dirty.clear();
        if let Some(fired) = fired {
            self.generation = self.generation.wrapping_add(1);
            if self.generation == 0 {
                self.marks.fill(0);
                self.generation = 1;
            }
            for &(s, _) in &net.channels[fired].changes {
                if counts[s] >= self.low[s] && counts[s] <= self.high[s] {
                    continue;
                }
                (self.low[s], self.high[s]) = interval(counts[s], self.pinned[s]);
                for &channel in &net.readers[s] {
                    if self.marks[channel] != self.generation {
                        self.marks[channel] = self.generation;
                        self.dirty.push(channel);
                    }
                }
            }
            // Fixed channel order makes the pinning fallback deterministic.
            self.dirty.sort_unstable();
        } else {
            self.dirty.extend(0..net.channels.len());
        }
        for k in 0..self.dirty.len() {
            let i = self.dirty[k];
            let (lower, upper, positive) = match self.bound(net, plan, i, counts, values, params) {
                Some(bound) => bound,
                None => self.point(net, plan, i, counts, values, params, cache)?,
            };
            self.positive_lower -= usize::from(self.lower[i] > 0.0);
            self.positive_lower += usize::from(lower > 0.0);
            self.lower[i] = lower;
            tree.set(i, upper, positive);
            composition.update(i, upper)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    //! The bounds against brute force: every integer state of the box,
    //! evaluated with the engine's own channel arithmetic.
    use super::*;
    use crate::rules::model;
    use rand_core::RngCore;
    use serde_json::{Value, json};

    fn network(value: &Value) -> Network {
        let source: model::Model = serde_json::from_value(value.clone()).unwrap();
        let compiled = model::compile(&source).unwrap_or_else(|d| panic!("{d:?}"));
        Network::compile(&compiled).unwrap()
    }

    fn params(value: &Value) -> Vec<f64> {
        let source: model::Model = serde_json::from_value(value.clone()).unwrap();
        model::compile(&source).unwrap().params
    }

    /// Two attribute-free types X and G (and P with an attribute), with the
    /// given rules.
    fn model(rules: Value, parameters: Value) -> Value {
        json!({
            "schema": "crn-craft.model/1.0",
            "name": "bounds",
            "parameters": parameters,
            "types": [{"name": "X"}, {"name": "G"}, {"name": "P", "attributes": ["size"]}],
            "initial": [{"type": "X", "count": 3}, {"type": "G", "count": 3},
                        {"type": "P", "attributes": {"size": 1}, "count": 2},
                        {"type": "P", "attributes": {"size": 3}, "count": 1}],
            "observables": {"X": {"type": "X"}, "G": {"type": "G"},
                            "mass": {"type": "P", "weight": "x.size"},
                            "big": {"type": "P", "where": "x.size > 2"}},
            "rules": rules
        })
    }

    /// The engine's value of every channel at `counts`: projections
    /// refreshed, then the rate cache.
    fn exact(net: &Network, params: &[f64], counts: &[u64], n_obs: usize) -> Vec<Option<f64>> {
        let mut values = vec![0.0; n_obs];
        if net.refresh(counts, &mut values).is_err() {
            return vec![None; net.channels.len()];
        }
        let env = Env {
            params,
            observables: &values,
            attributes: [&[], &[]],
        };
        let mut cache = net.rate_cache();
        (0..net.channels.len())
            .map(|i| cache.channel(net, i, counts, &env).ok())
            .collect()
    }

    /// Every bound of every channel of `value` encloses the engine's value at
    /// every state of random boxes, and its positive minimum is at most every
    /// positive value. Returns (interval bounds, pinned bounds) checked.
    fn enclose(value: &Value, seed: u64, boxes: usize) -> (usize, usize) {
        let net = network(value);
        let params = params(value);
        let plan = Plan::new(&net);
        let n_obs = value["observables"].as_object().unwrap().len();
        let species = net.species.len();
        let mut rng = crate::sampling::stream(seed, 1);
        let (mut checked, mut pinned) = (0, 0);
        for _ in 0..boxes {
            // A box per species, within 0..=12, and a current state inside.
            let mut low = vec![0u64; species];
            let mut high = vec![0u64; species];
            let mut counts = vec![0u64; species];
            for s in 0..species {
                let a = rng.next_u64() % 10;
                let b = a + rng.next_u64() % 5;
                low[s] = a;
                high[s] = b;
                counts[s] = a + rng.next_u64() % (b - a + 1);
            }
            let mut values = vec![0.0; n_obs];
            assert!(net.refresh(&counts, &mut values).is_ok());
            let mut bounds = Bounds::default();
            bounds.reset(&net, &plan, &counts);
            for s in 0..species {
                if !bounds.pinned[s] {
                    bounds.low[s] = low[s];
                    bounds.high[s] = high[s];
                }
            }
            let (low, high) = (bounds.low.clone(), bounds.high.clone());
            let mut found = Vec::new();
            for i in 0..net.channels.len() {
                let before = (counts.clone(), values.clone());
                let bound = bounds.bound(&net, &plan, i, &mut counts, &mut values, &params);
                assert_eq!(before.0, counts, "counts restored");
                assert!(
                    before
                        .1
                        .iter()
                        .zip(&values)
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    "values restored"
                );
                match bound {
                    Some(_) => checked += 1,
                    None => pinned += 1,
                }
                found.push(bound);
            }
            // Every state of the box, the other species at every value of
            // theirs too.
            let mut state = low.clone();
            loop {
                let values = exact(&net, &params, &state, n_obs);
                for (i, bound) in found.iter().enumerate() {
                    let (Some((lower, upper, positive)), Some(a)) = (*bound, values[i]) else {
                        continue;
                    };
                    assert!(
                        lower <= a && a <= upper,
                        "channel {i} at {state:?}: {lower} <= {a} <= {upper}, box {low:?}..{high:?}"
                    );
                    if a > 0.0 {
                        assert!(
                            positive <= a,
                            "channel {i} at {state:?}: minimum {positive} > {a}"
                        );
                    }
                }
                let mut s = 0;
                while s < species {
                    if state[s] < high[s] {
                        state[s] += 1;
                        break;
                    }
                    state[s] = low[s];
                    s += 1;
                }
                if s == species {
                    break;
                }
            }
        }
        (checked, pinned)
    }

    #[test]
    fn every_integer_state_is_enclosed_in_mass_action_boxes() {
        let rules = json!([
            {"name": "source", "reactants": [], "rate": "k", "products": [{"type": "X"}]},
            {"name": "single", "reactants": [{"label": "x", "type": "X"}], "rate": "0.3", "products": []},
            {"name": "pair", "reactants": [{"label": "x", "type": "X"}, {"label": "y", "type": "X"}],
             "rate": "0.1", "products": [{"type": "G"}]},
            {"name": "cross", "reactants": [{"label": "x", "type": "X"}, {"label": "y", "type": "G"}],
             "rate": "0.03", "products": [{"type": "G"}]},
            {"name": "attributed", "reactants": [{"label": "x", "type": "G"},
                                                 {"label": "p", "type": "P", "factor": "p.size * 0.7"}],
             "rate": "1.3", "products": [{"type": "P", "copy": "p"}]},
            {"name": "third", "reactants": [{"label": "x", "type": "X", "count": 3}],
             "rate": "0.002", "products": []},
            {"name": "two_one", "reactants": [{"label": "x", "type": "X", "count": 2},
                                              {"label": "y", "type": "G"}],
             "rate": "0.07", "products": [{"type": "X"}]},
            {"name": "three", "reactants": [{"label": "x", "type": "X"}, {"label": "y", "type": "G"},
                                            {"label": "z", "type": "X"}],
             "rate": "0.011", "products": [{"type": "P", "attributes": {"size": 1}}]}
        ]);
        let value = model(rules, json!({"k": 2.5}));
        let (checked, pinned) = enclose(&value, 7, 300);
        assert_eq!(pinned, 0);
        assert!(checked > 2000, "{checked}");
    }

    #[test]
    fn every_integer_state_is_enclosed_in_expression_boxes() {
        // Michaelis–Menten and Hill laws of one species (enumerated), whole
        // propensities, a weighted filtered observable, higher-order rules
        // with expression rates, and steep laws whose neighbors underflow,
        // overflow or turn negative (pinned).
        let rules = json!([
            {"name": "mm", "reactants": [{"label": "x", "type": "X"}],
             "propensity": "vmax * X / (km + X)", "products": [{"type": "G"}]},
            {"name": "mm_rate", "reactants": [{"label": "x", "type": "G"}],
             "rate": "vmax / (km + X)", "products": []},
            {"name": "hill_act", "reactants": [], "propensity": "basal + vmax * G^n / (k^n + G^n)",
             "products": [{"type": "X"}]},
            {"name": "hill_rep", "reactants": [{"label": "x", "type": "G"}, {"label": "y", "type": "G"}],
             "propensity": "basal + vmax * k^n / (k^n + X^n)", "products": [{"type": "X"}]},
            {"name": "vacancy", "reactants": [{"label": "x", "type": "X"}],
             "rate": "max(0, 9 - G) / (1 + G)", "products": [{"type": "G"}]},
            {"name": "filtered", "reactants": [{"label": "x", "type": "X"}],
             "rate": "0.5 * big + 1", "products": []},
            {"name": "steep", "reactants": [], "propensity": "exp(-200 * X) + 1e-300",
             "products": []},
            {"name": "huge", "reactants": [{"label": "x", "type": "G"}],
             "rate": "exp(150 * G)", "products": []},
            {"name": "negative", "reactants": [{"label": "x", "type": "X"}],
             "rate": "2 - X", "products": []},
            {"name": "order", "reactants": [{"label": "x", "type": "X", "count": 2},
                                            {"label": "y", "type": "G"}],
             "rate": "1 / (1 + X)", "products": []},
            {"name": "whole_order", "reactants": [{"label": "x", "type": "X", "count": 3}],
             "propensity": "km * X", "products": []}
        ]);
        let value = model(
            rules,
            json!({"vmax": 4.0, "km": 2.5, "basal": 0.1, "k": 3.0, "n": 2.5}),
        );
        let net = network(&value);
        assert!(net.rule_reads.iter().all(|r| r.len() <= 1));
        let (checked, pinned) = enclose(&value, 9, 300);
        assert!(checked > 2000, "{checked}");
        assert!(pinned > 100, "{pinned}");
    }

    #[test]
    fn every_integer_state_is_enclosed_when_rates_read_several_species() {
        let rules = json!([
            {"name": "two_reads", "reactants": [{"label": "x", "type": "X"}],
             "propensity": "0.2 * X * G / (X + G + 1)", "products": []},
            {"name": "mass_rate", "reactants": [{"label": "x", "type": "X"}],
             "rate": "0.1 * mass + big", "products": []},
            {"name": "pair", "reactants": [{"label": "x", "type": "G"}, {"label": "y", "type": "G"}],
             "rate": "0.1", "products": []},
            {"name": "third", "reactants": [{"label": "x", "type": "X", "count": 3}],
             "rate": "0.002", "products": []}
        ]);
        let value = model(rules, json!({}));
        let net = network(&value);
        let plan = Plan::new(&net);
        assert!(plan.pinned.iter().all(|&p| p));
        let (checked, pinned) = enclose(&value, 10, 200);
        assert_eq!(pinned, 0);
        assert!(checked >= 800, "{checked}");
    }

    /// Near the count limit consecutive combination counts differ by about
    /// their rounding error, so the scaled product need not be monotone; the
    /// widened corners enclose every state of boxes there.
    #[test]
    fn scaled_combination_bounds_enclose_states_near_the_count_limit() {
        for (k, rate) in [
            (3u64, "1e-40"),
            (4, "1e-55"),
            (5, "1e-70"),
            (8, "1e-120"),
            (2, "1e-20"),
        ] {
            let value = json!({
                "schema": "crn-craft.model/1.0",
                "name": "large",
                "parameters": {},
                "types": [{"name": "X"}],
                "initial": [{"type": "X", "count": 1}],
                "observables": {"X": {"type": "X"}},
                "rules": [{"name": "order", "reactants": [{"label": "x", "type": "X", "count": k}],
                           "rate": rate, "products": []}]
            });
            let net = network(&value);
            let params = params(&value);
            let plan = Plan::new(&net);
            for top in [
                (1u64 << 52) + 977,
                (1u64 << 53) - 3,
                (1u64 << 51) - 1,
                1u64 << 40,
            ] {
                let mut counts = vec![top - 64];
                let mut values = vec![0.0];
                assert!(net.refresh(&counts, &mut values).is_ok());
                let mut bounds = Bounds::default();
                bounds.reset(&net, &plan, &counts);
                assert_eq!(bounds.high[0], top);
                let (lower, upper, positive) = bounds
                    .bound(&net, &plan, 0, &mut counts, &mut values, &params)
                    .unwrap();
                for x in bounds.low[0]..=bounds.high[0] {
                    let a = exact(&net, &params, &[x], 1)[0].unwrap();
                    assert!(lower <= a && a <= upper && positive <= a, "k {k} x {x}");
                }
            }
        }
    }

    #[test]
    fn enumerated_ranges_are_reused_per_rule_and_interval() {
        let value = model(
            json!([
                {"name": "mm", "reactants": [{"label": "x", "type": "G"}],
                 "propensity": "4 * X / (2 + X)", "products": []}
            ]),
            json!({}),
        );
        let net = network(&value);
        let params = params(&value);
        let plan = Plan::new(&net);
        let mut counts = vec![30, 50, 2, 1];
        let mut values = vec![0.0; 4];
        assert!(net.refresh(&counts, &mut values).is_ok());
        let mut bounds = Bounds::default();
        bounds.reset(&net, &plan, &counts);
        assert_eq!((bounds.low[0], bounds.high[0]), (26, 34));
        let first = bounds
            .bound(&net, &plan, 0, &mut counts, &mut values, &params)
            .unwrap();
        assert!(bounds.ranges[0].is_some());
        // The law is increasing in X: the bounds are its values at 26 and 34
        // (G = 50 keeps the whole propensity available over its interval).
        assert_eq!(first.1, 4.0 * 34.0 / (2.0 + 34.0));
        assert_eq!(first.0, 4.0 * 26.0 / (2.0 + 26.0));
        assert_eq!(first.2, first.0);
        let again = bounds
            .bound(&net, &plan, 0, &mut counts, &mut values, &params)
            .unwrap();
        assert_eq!(first, again);
    }

    /// Rules whose rates compile to the same program share one range cache
    /// (M4f), so a model with many copies of one law (W4) enumerates it once
    /// per interval; the bounds do not change.
    #[test]
    fn rules_with_one_rate_program_share_their_ranges() {
        let value = model(
            json!([
                {"name": "a", "reactants": [{"label": "x", "type": "G"}],
                 "propensity": "4 * X / (2 + X)", "products": []},
                {"name": "b", "reactants": [{"label": "x", "type": "X"}],
                 "propensity": "k * X / (2 + X)", "products": []},
                {"name": "c", "reactants": [{"label": "x", "type": "G"}],
                 "propensity": "4 * X / (3 + X)", "products": []},
                {"name": "d", "reactants": [{"label": "x", "type": "G"}],
                 "rate": "0.5", "products": []},
                {"name": "e", "reactants": [{"label": "x", "type": "G"}],
                 "rate": "0.5", "products": []}
            ]),
            json!({"k": 4}),
        );
        let net = network(&value);
        let params = params(&value);
        let plan = Plan::new(&net);
        // a and b compute 4 X / (2 + X); c another law; constants unshared.
        assert_eq!(plan.range_slot[0], plan.range_slot[1]);
        assert_ne!(plan.range_slot[0], plan.range_slot[2]);
        assert_ne!(plan.range_slot[3], plan.range_slot[4]);
        assert_eq!(plan.range_slots, 4);
        let mut counts = vec![30, 50, 2, 1];
        let mut values = vec![0.0; 4];
        assert!(net.refresh(&counts, &mut values).is_ok());
        let mut shared = Bounds::default();
        shared.reset(&net, &plan, &counts);
        let mut alone = Bounds::default();
        for i in 0..net.channels.len() {
            // A fresh cache for each channel: nothing reused.
            alone.reset(&net, &plan, &counts);
            let expected = alone.bound(&net, &plan, i, &mut counts, &mut values, &params);
            let got = shared.bound(&net, &plan, i, &mut counts, &mut values, &params);
            assert_eq!(got, expected, "channel {i}");
        }
    }

    #[test]
    fn rules_that_read_two_species_pin_them_and_refresh_on_every_change() {
        let value = model(
            json!([
                {"name": "ratio", "reactants": [], "propensity": "X * G / (1 + X + G)", "products": []},
                {"name": "decay", "reactants": [{"label": "x", "type": "X"}], "rate": "1", "products": []}
            ]),
            json!({}),
        );
        let net = network(&value);
        let params = params(&value);
        let plan = Plan::new(&net);
        assert!(plan.pinned[0] && plan.pinned[1] && !plan.pinned[2]);
        let mut counts = vec![30, 20, 2, 1];
        let mut values = vec![0.0; 4];
        assert!(net.refresh(&counts, &mut values).is_ok());
        let mut bounds = Bounds::default();
        bounds.reset(&net, &plan, &counts);
        let mut tree = RateTree::default();
        tree.reset(2);
        let mut composition = Composition::default();
        composition.reset(2);
        let mut cache = net.rate_cache();
        bounds
            .refresh(
                &net,
                &plan,
                None,
                &mut counts,
                &mut values,
                &params,
                &mut cache,
                &mut tree,
                &mut composition,
            )
            .ok()
            .unwrap();
        let ratio = 30.0 * 20.0 / (1.0 + 30.0 + 20.0);
        assert_eq!(tree.total().ok(), Some(ratio + 30.0));
        assert_eq!(bounds.lower, vec![ratio, 30.0]);
        // One X decays: X is pinned, so both channels are refreshed.
        counts[0] = 29;
        assert!(net.refresh(&counts, &mut values).is_ok());
        cache.reset();
        bounds
            .refresh(
                &net,
                &plan,
                Some(1),
                &mut counts,
                &mut values,
                &params,
                &mut cache,
                &mut tree,
                &mut composition,
            )
            .ok()
            .unwrap();
        let ratio = 29.0 * 20.0 / (1.0 + 29.0 + 20.0);
        assert_eq!(tree.total().ok(), Some(ratio + 29.0));
    }

    #[test]
    fn neighboring_failures_pin_reads_without_masking_current_state_errors() {
        // 2 - X is negative above X = 2: an interval around X = 1 is pinned,
        // and the point bound is the value at X = 1.
        let value = model(
            json!([
                {"name": "negative", "reactants": [{"label": "x", "type": "G"}],
                 "rate": "2 - X", "products": []}
            ]),
            json!({}),
        );
        let net = network(&value);
        let params = params(&value);
        let plan = Plan::new(&net);
        let mut counts = vec![1, 3, 2, 1];
        let mut values = vec![0.0; 4];
        assert!(net.refresh(&counts, &mut values).is_ok());
        let mut bounds = Bounds::default();
        bounds.reset(&net, &plan, &counts);
        let mut tree = RateTree::default();
        tree.reset(1);
        let mut composition = Composition::default();
        composition.reset(1);
        let mut cache = net.rate_cache();
        bounds
            .refresh(
                &net,
                &plan,
                None,
                &mut counts,
                &mut values,
                &params,
                &mut cache,
                &mut tree,
                &mut composition,
            )
            .ok()
            .unwrap();
        assert_eq!(tree.total().ok(), Some(3.0));
        assert!(bounds.pinned[0] && bounds.pinned[1]);
        // At X = 5 the current state itself is invalid: the replicate fails.
        let mut counts = vec![5, 3, 2, 1];
        assert!(net.refresh(&counts, &mut values).is_ok());
        bounds.reset(&net, &plan, &counts);
        let mut cache = net.rate_cache();
        let Err(Halt::Numerical(error)) = bounds.refresh(
            &net,
            &plan,
            None,
            &mut counts,
            &mut values,
            &params,
            &mut cache,
            &mut tree,
            &mut composition,
        ) else {
            panic!("invalid rate accepted");
        };
        assert_eq!(error.code, "CRN-R030");
    }

    #[test]
    fn the_fresh_global_factor_is_the_cached_one() {
        // Vacancy rates (shared and single-channel), expressions, constants.
        let value = json!({
            "schema": "crn-craft.model/1.0",
            "name": "globals",
            "parameters": {"K": 40, "c": 0.3},
            "types": [{"name": "A", "attributes": ["s"]}, {"name": "B"}],
            "initial": [{"type": "A", "attributes": {"s": 1}, "count": 10},
                        {"type": "A", "attributes": {"s": 2}, "count": 7},
                        {"type": "B", "count": 5}],
            "observables": {"N": {"type": "A"}, "D": {"type": "A", "weight": "x.s"},
                            "B": {"type": "B"}},
            "rules": [
                {"name": "rep", "reactants": [{"label": "x", "type": "A", "factor": "x.s"}],
                 "rate": "max(0, K - N) / D * c", "products": [{"type": "A", "copy": "x"},
                                                               {"type": "A", "copy": "x"}]},
                {"name": "one", "reactants": [{"label": "x", "type": "B"}],
                 "rate": "max(0, K - N) / D * 2", "products": [{"type": "B"}, {"type": "B"}]},
                {"name": "expr", "reactants": [{"label": "x", "type": "B"}],
                 "rate": "sqrt(B) + log(1 + N)", "products": []},
                {"name": "deg", "reactants": [{"label": "x", "type": "A"}], "rate": "c", "products": []}
            ]
        });
        let net = network(&value);
        let params = params(&value);
        let mut rng = crate::sampling::stream(3, 3);
        for _ in 0..200 {
            let counts: Vec<u64> = (0..net.species.len())
                .map(|_| rng.next_u64() % 30)
                .collect();
            let mut values = vec![0.0; 3];
            assert!(net.refresh(&counts, &mut values).is_ok());
            let env = Env {
                params: &params,
                observables: &values,
                attributes: [&[], &[]],
            };
            let mut cache = net.rate_cache();
            for i in 0..net.channels.len() {
                let Ok(cached) = cache.channel(&net, i, &counts, &env) else {
                    continue;
                };
                let fresh = net.value_with(i, net.global_fresh(i, &env), &counts);
                assert_eq!(
                    cached.to_bits(),
                    fresh.to_bits(),
                    "channel {i} at {counts:?}"
                );
            }
        }
    }
}
