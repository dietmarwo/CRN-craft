//! Fixed-instance lowering and direct SSA over the resulting channels.
use crate::rules::InstanceRecord;
use crate::rules::engine::{Clock, Halt, Limits, Words, numerical};
use crate::rules::expr::{Env, RateKernel};
use crate::rules::model::{Compiled, MAX_COUNT, ModelVersion, Reduce, Shape, Source};
use rand_core::RngCore;
use std::collections::HashMap;

pub(crate) struct Species {
    pub ty: usize,
    pub attrs: Vec<f64>,
    pub initial: u64,
}

pub(crate) struct Channel {
    pub rule: usize,
    /// One species per reactant molecule; for a higher-order rule one per
    /// pattern (the event participants).
    pub reactants: Vec<usize>,
    /// A higher-order rule's (species, count) list, sorted by species with
    /// equal species merged; `None` for the 0.6.0 shapes.
    pub multi: Option<Vec<(usize, u64)>>,
    pub products: Vec<usize>,
    pub changes: Vec<(usize, i64)>,
    pub(crate) factor: f64,
    /// A whole-rule propensity: the rate applies once while the reactants
    /// are present instead of per reactant tuple.
    pub whole: bool,
}

pub(crate) type Projection = (usize, Reduce, Vec<(usize, f64)>);

pub(crate) struct Network {
    pub species: Vec<Species>,
    pub channels: Vec<Channel>,
    pub(crate) observables: Vec<Projection>,
    pub(crate) rates: Vec<RateKernel>,
    pub affected: Vec<Vec<usize>>,
    plans: Vec<ChannelRate>,
    shared_rates: usize,
    vacancy_slots: Vec<Option<usize>>,
    shared_vacancies: usize,
    pub stop_affected: Vec<bool>,
    pub observe_affected: Vec<bool>,
    pub cold: Vec<usize>,
    pub(crate) changed: Vec<Vec<usize>>,
    hot_changed: Vec<Vec<usize>>,
    /// The channels whose propensity reads each species (its reactants and
    /// the species its rule's rate reads), in increasing order.
    pub(crate) readers: Vec<Vec<usize>>,
    /// The species each rule's rate reads through its observables.
    pub(crate) rule_reads: Vec<Vec<usize>>,
    /// The projections (indices into `observables`) each rule's rate reads.
    pub(crate) rule_projections: Vec<Vec<usize>>,
}

enum Factor {
    Constant(f64),
    Vacancy {
        slot: usize,
        limit: f64,
        count: usize,
        denominator: usize,
    },
    General(Option<usize>),
}

/// The reactant count factor of a channel: the number of reactant tuples
/// (for a whole-rule propensity 1 while every reactant is present), decided
/// once per channel.
#[derive(Clone, Copy)]
enum Count {
    /// No reactants: 1.
    Source,
    /// One molecule of a species: its count.
    Single(usize),
    /// Two molecules of one species: n (n - 1) / 2.
    Pair(usize),
    /// One molecule of each of two species: the product of the counts.
    Cross(usize, usize),
    /// Whole-rule propensities: 1 while the molecules are present.
    WholeSingle(usize),
    WholePair(usize),
    WholeCross(usize, usize),
    /// A higher-order channel: its reactants are in `Channel::multi`.
    HigherOrder,
}

struct ChannelRate {
    factor: Factor,
    count: Count,
    weight: f64,
    rule: usize,
    whole: bool,
}

impl ChannelRate {
    /// The reactant count factor of a channel of at most two molecules:
    /// the number of reactant tuples, or for a whole-rule propensity 1 while
    /// every reactant is present.
    #[allow(clippy::cast_precision_loss)]
    #[inline]
    fn count_weight(&self, counts: &[u64]) -> f64 {
        match self.count {
            Count::Source => 1.0,
            Count::Single(s) => counts[s] as f64,
            Count::Pair(s) => counts[s] as f64 * counts[s].saturating_sub(1) as f64 * 0.5,
            Count::Cross(s, t) => counts[s] as f64 * counts[t] as f64,
            // A whole-rule propensity needs one molecule per pattern.
            Count::WholeSingle(s) => f64::from(u8::from(counts[s] > 0)),
            Count::WholePair(s) => f64::from(u8::from(counts[s] > 1)),
            Count::WholeCross(s, t) => f64::from(u8::from(counts[s] > 0 && counts[t] > 0)),
            Count::HigherOrder => unreachable!("higher-order channels use their own path"),
        }
    }

    fn higher_order(&self) -> bool {
        matches!(self.count, Count::HigherOrder)
    }
}

/// Static lowering limits: instance tuples and channels of a
/// rules-model/0.1 to 0.4 model, and of a `crn-craft.model/1.0` model.
const MAX_TUPLES: usize = 256;
const MAX_CHANNELS: usize = 4096;
const MAX_TUPLES_1_0: usize = 4096;
const MAX_CHANNELS_1_0: usize = 65_536;

/// The lowered species by type and attribute bits, so a product finds its
/// slot without a scan over every species.
#[derive(Default)]
struct SpeciesIndex(HashMap<(usize, Vec<u64>), usize>);

impl SpeciesIndex {
    fn key(ty: usize, attrs: &[f64]) -> (usize, Vec<u64>) {
        (ty, attrs.iter().map(|x| x.to_bits()).collect())
    }

    /// The first species of type `ty` whose attributes have the same bits.
    fn slot(&self, ty: usize, attrs: &[f64]) -> Option<usize> {
        self.0.get(&Self::key(ty, attrs)).copied()
    }

    /// Records `slot` unless an equal species is already recorded.
    fn insert(&mut self, ty: usize, attrs: &[f64], slot: usize) {
        self.0.entry(Self::key(ty, attrs)).or_insert(slot);
    }
}

/// Sorted, merged species changes of one firing: -1 per reactant molecule,
/// +1 per product, zero sums left out.
fn changes_of(
    reactants: impl IntoIterator<Item = (usize, i64)>,
    products: &[usize],
) -> Vec<(usize, i64)> {
    let mut delta: Vec<(usize, i64)> = reactants
        .into_iter()
        .chain(products.iter().map(|&s| (s, 1)))
        .collect();
    delta.sort_by_key(|&(s, _)| s);
    let mut out: Vec<(usize, i64)> = Vec::with_capacity(delta.len());
    for (s, d) in delta {
        match out.last_mut() {
            Some(last) if last.0 == s => last.1 += d,
            _ => out.push((s, d)),
        }
    }
    out.retain(|&(_, d)| d != 0);
    out
}

/// The sorted union of `lists[k]` over `keys`.
fn union(lists: &[Vec<usize>], keys: impl IntoIterator<Item = usize>) -> Vec<usize> {
    let mut out: Vec<usize> = keys
        .into_iter()
        .flat_map(|k| lists[k].iter().copied())
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}

fn local(
    env: &Env<'_>,
    condition: Option<&crate::rules::expr::Expr>,
    weight: Option<&crate::rules::expr::Expr>,
) -> Result<Option<f64>, String> {
    let matched = condition.map_or(1.0, |e| e.eval(env));
    if matched.is_nan() {
        return Err("undefined fixed-instance predicate".into());
    }
    if matched == 0.0 {
        return Ok(None);
    }
    let w = weight.map_or(1.0, |e| e.eval(env));
    if w.is_finite() {
        Ok(Some(w))
    } else {
        Err("nonfinite fixed-instance weight".into())
    }
}

impl Network {
    /// Lowers a closed model. Every table is built from indices (species by
    /// type, readers by species, projections by species), so the work grows
    /// with the size of the tables rather than with the square of the
    /// number of species or channels.
    #[allow(clippy::too_many_lines)]
    pub fn compile(model: &Compiled) -> Result<Self, String> {
        let (max_tuples, max_channels) = if model.version == ModelVersion::V10 {
            (MAX_TUPLES_1_0, MAX_CHANNELS_1_0)
        } else {
            (MAX_TUPLES, MAX_CHANNELS)
        };
        let too_many = || format!("static backend supports up to {max_tuples} instances");
        let mut species: Vec<_> = model
            .initial
            .iter()
            .map(|(ty, attrs, count)| Species {
                ty: *ty,
                attrs: attrs.clone(),
                initial: *count,
            })
            .collect();
        if species.len() > max_tuples {
            return Err(too_many());
        }
        let mut index = SpeciesIndex::default();
        for (i, s) in species.iter().enumerate() {
            index.insert(s.ty, &s.attrs, i);
        }
        let env = Env {
            params: &model.params,
            observables: &[],
            attributes: [&[], &[]],
        };
        for rule in &model.rules {
            for p in &rule.products {
                for a in &p.attributes {
                    if matches!(a, Source::Expression(e) if !e.deterministic()) {
                        return Err("random product attributes require the dynamic engine".into());
                    }
                }
                if p.attributes
                    .iter()
                    .all(|a| matches!(a, Source::Expression(e) if e.parameter_only()))
                {
                    let attrs: Vec<_> = p
                        .attributes
                        .iter()
                        .map(|a| {
                            let Source::Expression(e) = a else {
                                unreachable!()
                            };
                            normalize(e.eval(&env))
                        })
                        .collect();
                    if attrs.iter().any(|x| !x.is_finite()) {
                        return Err("nonfinite fixed product attribute".into());
                    }
                    if index.slot(p.type_id, &attrs).is_none() {
                        if species.len() == max_tuples {
                            return Err(too_many());
                        }
                        index.insert(p.type_id, &attrs, species.len());
                        species.push(Species {
                            ty: p.type_id,
                            attrs,
                            initial: 0,
                        });
                    }
                }
            }
        }
        if species.len() > max_tuples {
            return Err(too_many());
        }
        // Species of each type, in increasing order.
        let mut by_type = vec![Vec::new(); model.types.len()];
        for (i, s) in species.iter().enumerate() {
            by_type[s.ty].push(i);
        }
        let classes = model
            .classes
            .iter()
            .map(|c| {
                by_type[c.type_id]
                    .iter()
                    .filter_map(|&i| {
                        let s = &species[i];
                        let env = Env {
                            attributes: [&s.attrs, &[]],
                            ..env
                        };
                        match local(&env, c.condition.as_ref(), c.factor.as_ref()) {
                            Ok(None | Some(0.0)) => None,
                            Ok(Some(w)) if w > 0.0 => Some(Ok((i, w))),
                            Ok(Some(_)) => Some(Err("negative fixed-instance factor".into())),
                            Err(e) => Some(Err(e)),
                        }
                    })
                    .collect::<Result<Vec<_>, String>>()
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut channels = Vec::new();
        for (r, rule) in model.rules.iter().enumerate() {
            // `reactants` holds one species per molecule, or per pattern with
            // `multi` the merged (species, count) list of a higher-order rule.
            let mut add = |reactants: &[usize],
                           multi: Option<Vec<(usize, u64)>>,
                           factor|
             -> Result<(), String> {
                if channels.len() >= max_channels {
                    return Err("static channel limit exceeded".into());
                }
                // Rules of more than two patterns have attribute-free
                // reactants only.
                let mut attrs = [&[][..]; 2];
                for (i, s) in reactants.iter().take(2).enumerate() {
                    attrs[i] = &species[*s].attrs;
                }
                let env = Env {
                    attributes: attrs,
                    ..env
                };
                let mut products = Vec::new();
                for p in &rule.products {
                    let values: Vec<_> = p
                        .attributes
                        .iter()
                        .map(|a| {
                            normalize(match a {
                                Source::Copy { reactant, index } => attrs[*reactant][*index],
                                Source::Expression(e) => e.eval(&env),
                            })
                        })
                        .collect();
                    let i = index
                        .slot(p.type_id, &values)
                        .ok_or("products can create new attribute tuples")?;
                    products.push(i);
                }
                let changes = match &multi {
                    None => changes_of(reactants.iter().map(|&s| (s, -1)), &products),
                    Some(list) => changes_of(
                        list.iter()
                            .map(|&(s, k)| (s, -i64::try_from(k).expect("validated count"))),
                        &products,
                    ),
                };
                channels.push(Channel {
                    rule: r,
                    reactants: reactants.to_vec(),
                    multi,
                    changes,
                    products,
                    factor,
                    whole: rule.whole,
                });
                Ok(())
            };
            if rule.whole {
                // Firing picks molecules uniformly, which one channel per
                // tuple could express only with the type's total count.
                let several = |c: usize| classes[c].len() > 1;
                let single = match &rule.shape {
                    Shape::Source => true,
                    Shape::Single(c) | Shape::Pair(c) => !several(*c),
                    Shape::Cross(c, d) => !several(*c) && !several(*d),
                    Shape::Multi(patterns) => patterns.iter().all(|&(c, _)| !several(c)),
                };
                if !single {
                    return Err(format!(
                        "the propensity of rule `{}` reads a type with several instance tuples; static methods support whole-rule propensities over types with one tuple",
                        rule.name
                    ));
                }
            }
            match &rule.shape {
                Shape::Source => add(&[], None, 1.0)?,
                &Shape::Single(c) => {
                    for &(i, w) in &classes[c] {
                        add(&[i], None, w)?;
                    }
                }
                &Shape::Cross(c, d) => {
                    for &(i, w) in &classes[c] {
                        for &(j, v) in &classes[d] {
                            add(&[i, j], None, w * v)?;
                        }
                    }
                }
                &Shape::Pair(c) => {
                    if classes[c].len() > 1 {
                        return Err("static identical pairs require one matching instance".into());
                    }
                    if let Some(&(i, w)) = classes[c].first() {
                        add(&[i, i], None, w * w)?;
                    }
                }
                Shape::Multi(patterns) => {
                    // One channel per choice of a member for every pattern.
                    // Patterns of one class are over an attribute-free type,
                    // whose class has at most one member, so they agree.
                    if patterns.iter().all(|&(c, _)| !classes[c].is_empty()) {
                        let mut choice = vec![0; patterns.len()];
                        loop {
                            let members: Vec<(usize, f64)> = patterns
                                .iter()
                                .zip(&choice)
                                .map(|(&(c, _), &j)| classes[c][j])
                                .collect();
                            let participants: Vec<usize> =
                                members.iter().map(|&(s, _)| s).collect();
                            let mut merged: Vec<(usize, u64)> = participants
                                .iter()
                                .zip(patterns)
                                .map(|(&s, &(_, k))| (s, k))
                                .collect();
                            merged.sort_by_key(|&(s, _)| s);
                            merged.dedup_by(|next, kept| {
                                let same = next.0 == kept.0;
                                if same {
                                    kept.1 += next.1;
                                }
                                same
                            });
                            let factor = members.iter().map(|&(_, w)| w).product();
                            add(&participants, Some(merged), factor)?;
                            // The next choice, last pattern fastest.
                            let mut k = patterns.len();
                            loop {
                                if k == 0 {
                                    break;
                                }
                                k -= 1;
                                choice[k] += 1;
                                if choice[k] < classes[patterns[k].0].len() {
                                    break;
                                }
                                choice[k] = 0;
                            }
                            if choice.iter().all(|&j| j == 0) {
                                break;
                            }
                        }
                    }
                }
            }
        }
        // The type of each model observable (ATP and continuous variables
        // have none).
        let mut owner = vec![None; model.observables.len()];
        for (t, info) in model.types.iter().enumerate() {
            for &o in &info.observables {
                owner[o] = Some(t);
            }
        }
        let mut observables = Vec::new();
        // The projection of each model observable.
        let mut projection = vec![None; model.observables.len()];
        for (o, info) in model.observables.iter().enumerate() {
            if model.atp.as_ref().is_some_and(|atp| atp.index == o)
                || model.continuous.iter().any(|v| v.index == o)
            {
                continue;
            }
            let mut terms = Vec::new();
            for &i in owner[o].map_or(&[][..], |t| by_type[t].as_slice()) {
                let s = &species[i];
                let env = Env {
                    attributes: [&s.attrs, &[]],
                    ..env
                };
                if let Some(w) = local(&env, info.condition.as_ref(), info.weight.as_ref())? {
                    terms.push((i, w));
                }
            }
            projection[o] = Some(observables.len());
            observables.push((o, info.reduce, terms));
        }
        // Projections read by each rule's rate.
        let rule_projections: Vec<Vec<usize>> = model
            .rules
            .iter()
            .map(|rule| {
                rule.rate
                    .read_list()
                    .into_iter()
                    .filter_map(|o| projection[o])
                    .collect()
            })
            .collect();
        // Species read by each rule's rate through its observables.
        let rule_reads: Vec<Vec<usize>> = rule_projections
            .iter()
            .map(|read| {
                let mut out: Vec<usize> = read
                    .iter()
                    .flat_map(|&p| observables[p].2.iter().map(|&(s, _)| s))
                    .collect();
                out.sort_unstable();
                out.dedup();
                out
            })
            .collect();
        // The channels whose rate reads each species, in increasing order.
        let mut readers = vec![Vec::new(); species.len()];
        for (j, c) in channels.iter().enumerate() {
            let mut reads: Vec<usize> = rule_reads[c.rule]
                .iter()
                .copied()
                .chain(c.reactants.iter().copied())
                .collect();
            reads.sort_unstable();
            reads.dedup();
            for s in reads {
                readers[s].push(j);
            }
        }
        let affected = channels
            .iter()
            .map(|c| union(&readers, c.changes.iter().map(|&(s, _)| s)))
            .collect();
        // The projections whose terms hold each species.
        let mut holders = vec![Vec::new(); species.len()];
        for (i, (_, _, terms)) in observables.iter().enumerate() {
            for &(s, _) in terms {
                if holders[s].last() != Some(&i) {
                    holders[s].push(i);
                }
            }
        }
        let changed: Vec<Vec<usize>> = channels
            .iter()
            .map(|c| union(&holders, c.changes.iter().map(|&(s, _)| s)))
            .collect();
        let rates: Vec<_> = model
            .rules
            .iter()
            .map(|r| r.rate.kernel(&model.params))
            .collect();
        let mut uses = vec![0usize; rates.len()];
        for c in &channels {
            uses[c.rule] += 1;
        }
        let mut shared_rates = 0;
        let rate_slots: Vec<_> = uses
            .iter()
            .zip(&rates)
            .map(|(&uses, rate)| {
                (uses > 1 && !matches!(rate, RateKernel::Constant(_))).then(|| {
                    let slot = shared_rates;
                    shared_rates += 1;
                    slot
                })
            })
            .collect();
        let mut groups = std::collections::BTreeMap::new();
        for (rate, &uses) in rates.iter().zip(&uses) {
            if uses > 0 {
                if let RateKernel::Vacancy {
                    limit,
                    count,
                    denominator,
                    ..
                } = rate
                {
                    *groups
                        .entry((limit.to_bits(), *count, *denominator))
                        .or_insert(0usize) += 1;
                }
            }
        }
        let slots: std::collections::BTreeMap<_, _> = groups
            .into_iter()
            .filter(|(_, n)| *n > 1)
            .enumerate()
            .map(|(i, (key, _))| (key, i))
            .collect();
        let plans = channels
            .iter()
            .map(|c| {
                let factor = match &rates[c.rule] {
                    RateKernel::Constant(k) => Factor::Constant(*k),
                    RateKernel::Vacancy {
                        limit,
                        count,
                        denominator,
                        ..
                    } => match slots.get(&(limit.to_bits(), *count, *denominator)) {
                        Some(&slot) if uses[c.rule] == 1 => Factor::Vacancy {
                            slot,
                            limit: *limit,
                            count: *count,
                            denominator: *denominator,
                        },
                        _ => Factor::General(rate_slots[c.rule]),
                    },
                    _ => Factor::General(rate_slots[c.rule]),
                };
                let count = match (c.multi.is_some(), c.reactants.as_slice(), c.whole) {
                    (true, _, _) => Count::HigherOrder,
                    (false, [], _) => Count::Source,
                    (false, [s], false) => Count::Single(*s),
                    (false, [s], true) => Count::WholeSingle(*s),
                    (false, [s, t], false) if s == t => Count::Pair(*s),
                    (false, [s, t], true) if s == t => Count::WholePair(*s),
                    (false, [s, t], false) => Count::Cross(*s, *t),
                    (false, [s, t], true) => Count::WholeCross(*s, *t),
                    (false, _, _) => unreachable!("at most two reactant molecules"),
                };
                ChannelRate {
                    factor,
                    count,
                    weight: c.factor,
                    rule: c.rule,
                    whole: c.whole,
                }
            })
            .collect();
        let vacancy_slots = rates
            .iter()
            .map(|rate| match rate {
                RateKernel::Vacancy {
                    limit,
                    count,
                    denominator,
                    ..
                } => slots.get(&(limit.to_bits(), *count, *denominator)).copied(),
                _ => None,
            })
            .collect();
        Ok(Self {
            species,
            channels,
            observables,
            rates,
            affected,
            plans,
            shared_rates,
            vacancy_slots,
            shared_vacancies: slots.len(),
            stop_affected: Vec::new(),
            observe_affected: Vec::new(),
            cold: Vec::new(),
            hot_changed: changed.clone(),
            changed,
            readers,
            rule_reads,
            rule_projections,
        })
    }

    pub fn initial(&self) -> Vec<u64> {
        self.species.iter().map(|s| s.initial).collect()
    }

    pub fn plan_observables(&mut self, model: &Compiled, spec: &crate::rules::observe::Spec) {
        let mut stop = vec![false; model.observables.len()];
        if let Some(expr) = &model.stop_program {
            expr.reads(&mut stop);
        }
        self.stop_affected = self
            .changed
            .iter()
            .map(|row| row.iter().any(|&i| stop[self.observables[i].0]))
            .collect();
        self.observe_affected = if let Some(indices) = spec.sparse_observables() {
            let tracked: Vec<_> = indices.collect();
            self.changed
                .iter()
                .map(|row| {
                    row.iter()
                        .any(|&i| tracked.contains(&self.observables[i].0))
                })
                .collect()
        } else {
            vec![spec.enabled(); self.channels.len()]
        };
        let mut hot = stop;
        spec.reads(&mut hot);
        for rule in &model.rules {
            rule.rate.reads(&mut hot);
        }
        let mut cold = vec![false; self.observables.len()];
        if model.atp.is_none() && model.continuous.is_empty() {
            for (i, (index, _, terms)) in self.observables.iter().enumerate() {
                // Retain eager validation for any projection that might
                // overflow: a deferred sum of at most `bound` products, each
                // at most f64::MAX / bound, stays finite. The bound is 1024
                // for every projection of up to 1024 terms (all of them
                // before model/1.0, whose lowering has at most 256 species).
                #[allow(clippy::cast_precision_loss)]
                let bound = terms.len().max(1024) as f64;
                if !hot[*index]
                    && terms
                        .iter()
                        .all(|&(_, w)| w.abs() <= f64::MAX / (bound * MAX_COUNT as f64))
                {
                    cold[i] = true;
                    self.cold.push(i);
                }
            }
        }
        for row in &mut self.hot_changed {
            row.retain(|&i| !cold[i]);
        }
    }

    #[inline]
    pub fn refresh_channel(
        &self,
        channel: usize,
        counts: &[u64],
        values: &mut [f64],
        full: bool,
    ) -> Result<(), Halt> {
        let changed = if full {
            &self.changed
        } else {
            &self.hot_changed
        };
        let indices = &changed[channel];
        if indices.is_empty() {
            return Ok(());
        }
        self.refresh_indices(indices.iter().copied(), counts, values)
    }

    pub fn refresh_generated(
        &self,
        channel: usize,
        counts: &[u64],
        values: &mut [f64],
        full: bool,
        kernel: &dyn crate::rules::RulesKernel,
    ) -> Result<(), Halt> {
        let changed = if full {
            &self.changed
        } else {
            &self.hot_changed
        };
        for &i in &changed[channel] {
            let value = kernel.projection(i, counts);
            if !value.is_finite() {
                return Err(numerical(
                    "/observables",
                    "observable overflows",
                    "Rescale weights.",
                ));
            }
            values[self.observables[i].0] = value;
        }
        Ok(())
    }

    #[cold]
    #[inline(never)]
    pub fn rate_error(&self, channel: usize) -> Halt {
        let c = &self.channels[channel];
        numerical(
            &format!(
                "/rules/{}/{}",
                c.rule,
                if c.whole { "propensity" } else { "rate" }
            ),
            "invalid fixed-channel propensity",
            "Rates must be finite and nonnegative.",
        )
    }

    pub fn refresh_cold(&self, counts: &[u64], values: &mut [f64]) -> Result<(), Halt> {
        self.refresh_indices(self.cold.iter().copied(), counts, values)
    }

    pub fn refresh_changed(
        &self,
        changed: &[bool],
        counts: &[u64],
        values: &mut [f64],
    ) -> Result<(), Halt> {
        let indices = self
            .observables
            .iter()
            .enumerate()
            .filter_map(|(i, (_, _, terms))| terms.iter().any(|&(s, _)| changed[s]).then_some(i));
        self.refresh_indices(indices, counts, values)
    }

    pub fn refresh(&self, counts: &[u64], values: &mut [f64]) -> Result<(), Halt> {
        self.refresh_indices(0..self.observables.len(), counts, values)
    }

    fn refresh_indices(
        &self,
        indices: impl IntoIterator<Item = usize>,
        counts: &[u64],
        values: &mut [f64],
    ) -> Result<(), Halt> {
        for index in indices {
            let i = &self.observables[index].0;
            let value = self.project(index, counts);
            if !value.is_finite() {
                return Err(numerical(
                    "/observables",
                    "observable overflows",
                    "Rescale weights.",
                ));
            }
            values[*i] = value;
        }
        Ok(())
    }

    /// The value of projection `index` at `counts`, unchecked.
    #[allow(clippy::cast_precision_loss)]
    #[inline]
    pub(crate) fn project(&self, index: usize, counts: &[u64]) -> f64 {
        let (_, reduce, terms) = &self.observables[index];
        match reduce {
            Reduce::Sum => terms.iter().map(|&(s, w)| counts[s] as f64 * w).sum(),
            Reduce::Distinct => terms
                .iter()
                .filter(|&&(s, _)| counts[s] > 0)
                .map(|&(_, w)| w)
                .sum(),
            Reduce::Max => terms
                .iter()
                .filter(|&&(s, _)| counts[s] > 0)
                .map(|&(_, w)| w)
                .reduce(f64::max)
                .unwrap_or(0.0),
            Reduce::Min => terms
                .iter()
                .filter(|&&(s, _)| counts[s] > 0)
                .map(|&(_, w)| w)
                .reduce(f64::min)
                .unwrap_or(0.0),
        }
    }

    /// The reactant molecules channel `i` needs: `(species, count)` pairs,
    /// one per species.
    pub(crate) fn needs(&self, i: usize) -> Vec<(usize, u64)> {
        let c = &self.channels[i];
        match (&c.multi, c.reactants.as_slice()) {
            (Some(list), _) => list.clone(),
            (None, [s, t]) if s == t => vec![(*s, 2)],
            (None, reactants) => reactants.iter().map(|&s| (s, 1)).collect(),
        }
    }

    /// The constant global factor of channel `i`, if its rate is one.
    pub(crate) fn constant_global(&self, i: usize) -> Option<f64> {
        match self.plans[i].factor {
            Factor::Constant(k) => Some(k),
            _ => None,
        }
    }

    /// Whether channel `i`'s value is a scaled combination count
    /// ([`crate::rules::combination_rate`]) rather than a product of at
    /// most three correctly rounded operations.
    pub(crate) fn scaled_combination(&self, i: usize) -> bool {
        let c = &self.plans[i];
        c.higher_order() && !c.whole
    }

    /// The global factor of channel `i` at `env`, computed without the
    /// per-state cache: the value [`RateCache::channel`] multiplies by.
    pub(crate) fn global_fresh(&self, i: usize, env: &Env<'_>) -> f64 {
        let c = &self.plans[i];
        match c.factor {
            Factor::Constant(k) => k,
            Factor::Vacancy {
                limit,
                count,
                denominator,
                ..
            } => {
                let base = (limit - env.observables[count]).max(0.0) / env.observables[denominator];
                let RateKernel::Vacancy { factors, .. } = &self.rates[c.rule] else {
                    unreachable!()
                };
                factors.iter().fold(base, |x, f| x * f)
            }
            Factor::General(_) => {
                let kernel = &self.rates[c.rule];
                if self.vacancy_slots[c.rule].is_some() {
                    let RateKernel::Vacancy {
                        limit,
                        count,
                        denominator,
                        factors,
                    } = kernel
                    else {
                        unreachable!()
                    };
                    let base =
                        (limit - env.observables[*count]).max(0.0) / env.observables[*denominator];
                    factors.iter().fold(base, |x, f| x * f)
                } else {
                    kernel.eval(env)
                }
            }
        }
    }

    /// Channel `i`'s value with global factor `global` at `counts`, with the
    /// arithmetic of [`RateCache::channel`]: 0 while a reactant is missing.
    /// The value is not checked.
    pub(crate) fn value_with(&self, i: usize, global: f64, counts: &[u64]) -> f64 {
        let c = &self.plans[i];
        if c.higher_order() {
            let reactants = self.channels[i]
                .multi
                .as_deref()
                .expect("higher-order channel");
            if reactants.iter().any(|&(s, k)| counts[s] < k) {
                return 0.0;
            }
            return if c.whole {
                global * (1.0 * c.weight)
            } else {
                crate::rules::combination_rate(global, c.weight, counts, reactants)
            };
        }
        let weight = c.count_weight(counts);
        if weight == 0.0 {
            return 0.0;
        }
        global * (weight * c.weight)
    }

    pub fn rate_cache(&self) -> RateCache {
        RateCache {
            rules: vec![f64::NAN; self.shared_rates],
            vacancies: vec![f64::NAN; self.shared_vacancies],
        }
    }

    #[cfg(test)]
    pub fn rate(&self, i: usize, counts: &[u64], env: &Env<'_>) -> Result<f64, Halt> {
        self.rate_cache().channel(self, i, counts, env)
    }

    #[inline]
    pub fn apply(
        &self,
        i: usize,
        counts: &mut [u64],
        live: &mut u64,
        limits: &Limits,
    ) -> Result<u64, Halt> {
        let changes = &self.channels[i].changes;
        let mut next_live = *live;
        let mut created = 0;
        for &(s, d) in changes {
            let Some(next) = counts[s].checked_add_signed(d).filter(|&n| n <= MAX_COUNT) else {
                return Err(count_range());
            };
            next_live = next_live + u64::from(next > 0) - u64::from(counts[s] > 0);
            created += u64::from(next > 0 && counts[s] == 0);
        }
        if next_live > limits.max_instances {
            return Err(instance_limit());
        }
        for &(s, d) in changes {
            // Validated above.
            counts[s] = counts[s].wrapping_add_signed(d);
        }
        *live = next_live;
        Ok(created)
    }

    pub fn instances(&self, counts: &[u64]) -> Vec<InstanceRecord> {
        self.species
            .iter()
            .zip(counts)
            .filter(|&(_, &n)| n > 0)
            .map(|(s, &n)| InstanceRecord {
                type_index: s.ty,
                attributes: s.attrs.clone(),
                count: n,
            })
            .collect()
    }
}

/// Lazy per-state global factors. A channel with no reactants available
/// must still skip evaluating its global factor, even when another channel
/// of the same rule is present in the network.
pub(crate) struct RateCache {
    rules: Vec<f64>,
    vacancies: Vec<f64>,
}

impl RateCache {
    pub fn reset(&mut self) {
        self.rules.fill(f64::NAN);
        self.vacancies.fill(f64::NAN);
    }

    fn rule(&mut self, net: &Network, rule: usize, slot: Option<usize>, env: &Env<'_>) -> f64 {
        if let Some(i) = slot {
            if !self.rules[i].is_nan() {
                return self.rules[i];
            }
        }
        let kernel = &net.rates[rule];
        let value = if let Some(i) = net.vacancy_slots[rule] {
            let RateKernel::Vacancy {
                limit,
                count,
                denominator,
                factors,
            } = kernel
            else {
                unreachable!()
            };
            let base = &mut self.vacancies[i];
            if base.is_nan() {
                *base = (limit - env.observables[*count]).max(0.0) / env.observables[*denominator];
            }
            factors.iter().fold(*base, |x, f| x * f)
        } else {
            kernel.eval(env)
        };
        if let Some(i) = slot {
            self.rules[i] = value;
        }
        value
    }

    /// Channel `i`'s propensity at `counts` and `env`, or the rate error of
    /// its rule if the value is not finite and nonnegative.
    #[inline]
    pub fn channel(
        &mut self,
        net: &Network,
        i: usize,
        counts: &[u64],
        env: &Env<'_>,
    ) -> Result<f64, Halt> {
        match self.value(net, i, counts, env) {
            Some(value) => Ok(value),
            None => Err(net.rate_error(i)),
        }
    }

    /// [`Self::channel`]'s value, `None` for an invalid one.
    #[inline]
    fn value(&mut self, net: &Network, i: usize, counts: &[u64], env: &Env<'_>) -> Option<f64> {
        let c = &net.plans[i];
        if c.higher_order() {
            return self.higher_order(net, i, counts, env);
        }
        let weight = c.count_weight(counts);
        if weight == 0.0 {
            return Some(0.0);
        }
        let global = match c.factor {
            Factor::Constant(k) => k,
            _ => self.global(net, c, env),
        };
        let value = global * (weight * c.weight);
        (value.is_finite() && value >= 0.0).then_some(value)
    }

    /// The global factor of a channel whose reactants are available.
    #[inline(never)]
    fn global(&mut self, net: &Network, c: &ChannelRate, env: &Env<'_>) -> f64 {
        match c.factor {
            Factor::Constant(k) => k,
            Factor::Vacancy {
                slot,
                limit,
                count,
                denominator,
            } => {
                let base = &mut self.vacancies[slot];
                if base.is_nan() {
                    *base =
                        (limit - env.observables[count]).max(0.0) / env.observables[denominator];
                }
                let RateKernel::Vacancy { factors, .. } = &net.rates[c.rule] else {
                    unreachable!()
                };
                factors.iter().fold(*base, |x, f| x * f)
            }
            Factor::General(slot) => self.rule(net, c.rule, slot, env),
        }
    }

    /// A higher-order channel: 0 unless every reactant type has its count;
    /// then the whole propensity, or rate × factor × Π C(n, k) as the
    /// overflow-safe scaled product of [`crate::rules::combination_rate`].
    #[inline(never)]
    fn higher_order(
        &mut self,
        net: &Network,
        i: usize,
        counts: &[u64],
        env: &Env<'_>,
    ) -> Option<f64> {
        let c = &net.plans[i];
        let reactants = net.channels[i]
            .multi
            .as_deref()
            .expect("higher-order channel");
        if reactants.iter().any(|&(s, k)| counts[s] < k) {
            return Some(0.0);
        }
        let global = self.global(net, c, env);
        let value = if c.whole {
            global * (1.0 * c.weight)
        } else {
            crate::rules::combination_rate(global, c.weight, counts, reactants)
        };
        (value.is_finite() && value >= 0.0).then_some(value)
    }
}

pub(crate) fn check_instances(live: u64, limits: &Limits) -> Result<(), Halt> {
    if live <= limits.max_instances {
        return Ok(());
    }
    Err(instance_limit())
}

#[cold]
#[inline(never)]
fn instance_limit() -> Halt {
    Halt::Limit(
        crate::rules::Termination::InstanceLimit,
        crate::diagnostic::Diagnostic::error(
            "CRN-R042",
            "/budget/max_instances",
            "instance limit reached",
            "Increase max_instances.",
        ),
    )
}

#[cold]
#[inline(never)]
fn count_range() -> Halt {
    numerical(
        "/rules",
        "fixed count outside 0..2^53",
        "Reduce population.",
    )
}

fn normalize(x: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x }
}

pub(crate) fn direct<R: RngCore>(
    rates: &[f64],
    time: f64,
    words: &mut Words<'_, R>,
    clock: &mut Clock,
) -> Result<Option<(f64, usize)>, Halt> {
    direct_with(rates, time, words, clock, |words, total| {
        // Open uniform: u = 0 would give a zero wait. Other words are
        // unchanged. time + (-log(1 - u) / total) is time - log(1 - u) / total.
        Ok(-libm::log(1.0 - words.open()?) / total)
    })
}

/// From this many channels on, selection finds the smallest positive rate
/// while summing instead of testing every rate in a second pass. Some
/// positive rate is below the resolution floor exactly when the smallest
/// one is, so the outcome is the same.
pub(crate) const FUSED_SELECTION: usize = 32;

/// The total in `Iterator::sum`'s order and the smallest positive rate
/// (infinite if none is positive).
fn fused_total(rates: &[f64]) -> (f64, f64) {
    let mut total = -0.0;
    let mut smallest = f64::INFINITY;
    for &a in rates {
        total += a;
        if a > 0.0 && a < smallest {
            smallest = a;
        }
    }
    (total, smallest)
}

/// The error of a linear selection whose total is not finite (`CRN-R030`) or
/// whose smallest positive rate is below 2⁻⁵³ of the finite total
/// (`CRN-R031`).
#[cold]
pub(crate) fn unresolved(total: f64) -> Halt {
    // Static channels are selected in one flat level; the tree engine
    // resolves rules and instances separately and can still run.
    if total.is_finite() {
        return Halt::Numerical(crate::diagnostic::Diagnostic::error(
            "CRN-R031",
            "/rules",
            "a positive channel is below 53-bit resolution of the total",
            "Use the tree method or rescale rates.",
        ));
    }
    numerical(
        "/rules",
        "unresolvable total propensity",
        "Rescale rates or populations.",
    )
}

/// Linear Static Direct selection. `wait` draws the waiting time for the
/// total; `clock` decides whether a wait that does not advance the rounded
/// time is accepted.
pub(crate) fn direct_with<R: RngCore>(
    rates: &[f64],
    time: f64,
    words: &mut Words<'_, R>,
    clock: &mut Clock,
    wait: impl FnOnce(&mut Words<'_, R>, f64) -> Result<f64, Halt>,
) -> Result<Option<(f64, usize)>, Halt> {
    let (total, smallest) = if rates.len() >= FUSED_SELECTION {
        fused_total(rates)
    } else {
        (rates.iter().sum::<f64>(), f64::NAN)
    };
    if total == 0.0 {
        return Ok(None);
    }
    let floor = total * crate::rules::engine::RESOLUTION;
    let below = if smallest.is_nan() {
        rates.iter().any(|&a| a > 0.0 && a < floor)
    } else {
        smallest < floor
    };
    if !total.is_finite() || below {
        return Err(unresolved(total));
    }
    let wait = wait(words, total)?;
    let next = clock.advance(time, wait)?;
    let mut target = words.uniform()? * total;
    let mut last = None;
    for (i, &a) in rates.iter().enumerate().filter(|&(_, a)| *a > 0.0) {
        last = Some(i);
        if target < a {
            return Ok(Some((next, i)));
        }
        target -= a;
    }
    Ok(last.map(|i| (next, i)))
}

impl Network {
    pub fn supports_binomial(&self) -> bool {
        self.channels.iter().all(|c| {
            // A higher-order channel takes at least three molecules.
            c.multi.is_none()
                && c.reactants.len() == 1
                && c.changes.len() == 1
                && (c.changes[0] == (c.reactants[0], -1) || c.changes[0].1 == 1)
        })
    }
}

#[cfg(test)]
mod tests {
    //! The index-based lowering (M4a) against the quadratic construction it
    //! replaced, on every closed model of the repository.
    use super::*;
    use crate::rules::model;

    /// The tables as the 0.6.0 lowering built them: dense deltas, a scan of
    /// every species per observable, every channel against every channel.
    #[allow(clippy::type_complexity)]
    fn quadratic(
        net: &Network,
        model: &Compiled,
    ) -> (
        Vec<Vec<(usize, i64)>>,
        Vec<Vec<(usize, f64)>>,
        Vec<Vec<usize>>,
        Vec<Vec<usize>>,
    ) {
        let species = &net.species;
        let changes = net
            .channels
            .iter()
            .map(|c| {
                let mut delta = vec![0i64; species.len()];
                match &c.multi {
                    None => {
                        for &s in &c.reactants {
                            delta[s] -= 1;
                        }
                    }
                    Some(list) => {
                        for &(s, k) in list {
                            delta[s] -= i64::try_from(k).unwrap();
                        }
                    }
                }
                for &s in &c.products {
                    delta[s] += 1;
                }
                delta
                    .into_iter()
                    .enumerate()
                    .filter(|(_, d)| *d != 0)
                    .collect()
            })
            .collect();
        let env = Env {
            params: &model.params,
            observables: &[],
            attributes: [&[], &[]],
        };
        let terms = net
            .observables
            .iter()
            .map(|&(o, _, _)| {
                let info = &model.observables[o];
                let mut terms = Vec::new();
                for (i, s) in species.iter().enumerate() {
                    if model.types[s.ty].observables.contains(&o) {
                        let env = Env {
                            attributes: [&s.attrs, &[]],
                            ..env
                        };
                        if let Some(w) =
                            local(&env, info.condition.as_ref(), info.weight.as_ref()).unwrap()
                        {
                            terms.push((i, w));
                        }
                    }
                }
                terms
            })
            .collect();
        let reads: Vec<_> = net
            .channels
            .iter()
            .map(|c| {
                let mut obs = vec![false; model.observables.len()];
                model.rules[c.rule].rate.reads(&mut obs);
                let mut flags = vec![false; species.len()];
                for (i, _, terms) in &net.observables {
                    if obs[*i] {
                        for &(s, _) in terms {
                            flags[s] = true;
                        }
                    }
                }
                for &s in &c.reactants {
                    flags[s] = true;
                }
                flags
            })
            .collect();
        let affected = net
            .channels
            .iter()
            .map(|c| {
                reads
                    .iter()
                    .enumerate()
                    .filter_map(|(j, flags)| {
                        (c.changes.iter().any(|&(s, _)| flags[s])).then_some(j)
                    })
                    .collect()
            })
            .collect();
        let changed = net
            .channels
            .iter()
            .map(|c| {
                net.observables
                    .iter()
                    .enumerate()
                    .filter_map(|(i, (_, _, terms))| {
                        terms
                            .iter()
                            .any(|&(s, _)| c.changes.iter().any(|&(x, _)| s == x))
                            .then_some(i)
                    })
                    .collect()
            })
            .collect();
        (changes, terms, affected, changed)
    }

    fn corpus() -> Vec<(String, Vec<u8>)> {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let mut out = Vec::new();
        for folder in ["models/rules", "crates/crn-craft-rules/tests/generated"] {
            for entry in std::fs::read_dir(root.join(folder)).unwrap() {
                let path = entry.unwrap().path();
                if path.to_string_lossy().ends_with(".rules.json") {
                    out.push((path.display().to_string(), std::fs::read(&path).unwrap()));
                }
            }
        }
        let networks = ["models/examples", "validation/ssa-acceleration/models"];
        for folder in networks {
            for entry in std::fs::read_dir(root.join(folder)).unwrap() {
                let path = entry.unwrap().path();
                if path.to_string_lossy().ends_with(".model.json") {
                    let bytes = std::fs::read(&path).unwrap();
                    if let Ok(rules) = crate::rules::convert::from_network(&bytes) {
                        out.push((path.display().to_string(), rules));
                    }
                }
            }
        }
        out.sort();
        out
    }

    #[test]
    fn index_based_lowering_equals_the_quadratic_construction() {
        let mut lowered = 0;
        for (name, bytes) in corpus() {
            let source: model::Model = serde_json::from_slice(&bytes).unwrap();
            let Ok(compiled) = model::compile(&source) else {
                continue;
            };
            let Ok(net) = Network::compile(&compiled) else {
                continue;
            };
            let (changes, terms, affected, changed) = quadratic(&net, &compiled);
            for (i, c) in net.channels.iter().enumerate() {
                assert_eq!(c.changes, changes[i], "{name}: changes of channel {i}");
            }
            for (i, (_, _, t)) in net.observables.iter().enumerate() {
                let bits = |v: &[(usize, f64)]| -> Vec<(usize, u64)> {
                    v.iter().map(|&(s, w)| (s, w.to_bits())).collect()
                };
                assert_eq!(bits(t), bits(&terms[i]), "{name}: terms of projection {i}");
            }
            assert_eq!(net.affected, affected, "{name}: affected");
            assert_eq!(net.changed, changed, "{name}: changed");
            lowered += 1;
        }
        // Every workload, example and closed rules model of the repository.
        assert!(lowered >= 30, "{lowered} models lowered");
    }
}
