//! Fixed-instance lowering and direct SSA over the resulting channels.
use crate::InstanceRecord;
use crate::engine::{Halt, Limits, Words, numerical};
use crate::expr::{Env, RateKernel};
use crate::model::{Compiled, MAX_COUNT, Reduce, Shape, Source};
use rand_core::RngCore;

pub(crate) struct Species {
    pub ty: usize,
    pub attrs: Vec<f64>,
    pub initial: u64,
}

pub(crate) struct Channel {
    pub rule: usize,
    pub reactants: Vec<usize>,
    pub products: Vec<usize>,
    pub changes: Vec<(usize, i64)>,
    pub(crate) factor: f64,
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

struct ChannelRate {
    factor: Factor,
    reactants: [usize; 2],
    arity: usize,
    weight: f64,
    rule: usize,
}

fn same(a: &[f64], b: &[f64]) -> bool {
    a.iter()
        .map(|x| x.to_bits())
        .eq(b.iter().map(|x| x.to_bits()))
}

fn slot(species: &[Species], ty: usize, attrs: &[f64]) -> Option<usize> {
    species
        .iter()
        .position(|s| s.ty == ty && same(&s.attrs, attrs))
}

fn local(
    env: &Env<'_>,
    condition: Option<&crate::expr::Expr>,
    weight: Option<&crate::expr::Expr>,
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
    pub fn compile(model: &Compiled) -> Result<Self, String> {
        let mut species: Vec<_> = model
            .initial
            .iter()
            .map(|(ty, attrs, count)| Species {
                ty: *ty,
                attrs: attrs.clone(),
                initial: *count,
            })
            .collect();
        if species.len() > 256 {
            return Err("static backend supports up to 256 instances".into());
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
                    if slot(&species, p.type_id, &attrs).is_none() {
                        if species.len() == 256 {
                            return Err("static backend supports up to 256 instances".into());
                        }
                        species.push(Species {
                            ty: p.type_id,
                            attrs,
                            initial: 0,
                        });
                    }
                }
            }
        }
        if species.len() > 256 {
            return Err("static backend supports up to 256 instances".into());
        }
        let classes = model
            .classes
            .iter()
            .map(|c| {
                species
                    .iter()
                    .enumerate()
                    .filter(|(_, s)| s.ty == c.type_id)
                    .filter_map(|(i, s)| {
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
            let mut add = |reactants: &[usize], factor| -> Result<(), String> {
                if channels.len() >= 4096 {
                    return Err("static channel limit exceeded".into());
                }
                let mut attrs = [&[][..]; 2];
                for (i, s) in reactants.iter().enumerate() {
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
                    let i = slot(&species, p.type_id, &values)
                        .ok_or("products can create new attribute tuples")?;
                    products.push(i);
                }
                let mut delta = vec![0i64; species.len()];
                for &s in reactants {
                    delta[s] -= 1;
                }
                for &s in &products {
                    delta[s] += 1;
                }
                channels.push(Channel {
                    rule: r,
                    reactants: reactants.to_vec(),
                    products,
                    factor,
                    changes: delta
                        .into_iter()
                        .enumerate()
                        .filter(|(_, d)| *d != 0)
                        .collect(),
                });
                Ok(())
            };
            match rule.shape {
                Shape::Source => add(&[], 1.0)?,
                Shape::Single(c) => {
                    for &(i, w) in &classes[c] {
                        add(&[i], w)?;
                    }
                }
                Shape::Cross(c, d) => {
                    for &(i, w) in &classes[c] {
                        for &(j, v) in &classes[d] {
                            add(&[i, j], w * v)?;
                        }
                    }
                }
                Shape::Pair(c) => {
                    if classes[c].len() > 1 {
                        return Err("static identical pairs require one matching instance".into());
                    }
                    if let Some(&(i, w)) = classes[c].first() {
                        add(&[i, i], w * w)?;
                    }
                }
            }
        }
        let mut observables = Vec::new();
        for (o, info) in model.observables.iter().enumerate() {
            if model.atp.as_ref().is_some_and(|atp| atp.index == o) {
                continue;
            }
            let mut terms = Vec::new();
            for (i, s) in species.iter().enumerate() {
                if model.types[s.ty].observables.contains(&o) {
                    let env = Env {
                        attributes: [&s.attrs, &[]],
                        ..env
                    };
                    if let Some(w) = local(&env, info.condition.as_ref(), info.weight.as_ref())? {
                        terms.push((i, w));
                    }
                }
            }
            observables.push((o, info.reduce, terms));
        }
        let reads: Vec<_> = channels
            .iter()
            .map(|c| {
                let mut obs = vec![false; model.observables.len()];
                model.rules[c.rule].rate.reads(&mut obs);
                let mut flags = vec![false; species.len()];
                for (i, _, terms) in &observables {
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
        let affected = channels
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
                let mut reactants = [0; 2];
                reactants[..c.reactants.len()].copy_from_slice(&c.reactants);
                ChannelRate {
                    factor,
                    reactants,
                    arity: c.reactants.len(),
                    weight: c.factor,
                    rule: c.rule,
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
        let changed: Vec<Vec<usize>> = channels
            .iter()
            .map(|c| {
                observables
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
        })
    }

    pub fn initial(&self) -> Vec<u64> {
        self.species.iter().map(|s| s.initial).collect()
    }

    pub fn plan_observables(&mut self, model: &Compiled, spec: &crate::observe::Spec) {
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
        if model.atp.is_none() {
            for (i, (index, _, terms)) in self.observables.iter().enumerate() {
                // Retain eager validation for any projection that might
                // overflow. Summing up to 256 products leaves ample headroom.
                if !hot[*index]
                    && terms
                        .iter()
                        .all(|&(_, w)| w.abs() <= f64::MAX / (1024.0 * MAX_COUNT as f64))
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
        self.refresh_indices(changed[channel].iter().copied(), counts, values)
    }

    pub fn refresh_generated(
        &self,
        channel: usize,
        counts: &[u64],
        values: &mut [f64],
        full: bool,
        kernel: &dyn crate::RulesKernel,
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

    pub fn rate_error(&self, channel: usize) -> Halt {
        numerical(
            &format!("/rules/{}/rate", self.channels[channel].rule),
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

    #[allow(clippy::cast_precision_loss)]
    fn refresh_indices(
        &self,
        indices: impl IntoIterator<Item = usize>,
        counts: &[u64],
        values: &mut [f64],
    ) -> Result<(), Halt> {
        for index in indices {
            let (i, reduce, terms) = &self.observables[index];
            let value = match reduce {
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
            };
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
            let next = counts[s]
                .checked_add_signed(d)
                .filter(|&n| n <= MAX_COUNT)
                .ok_or_else(|| {
                    numerical(
                        "/rules",
                        "fixed count outside 0..2^53",
                        "Reduce population.",
                    )
                })?;
            next_live = next_live + u64::from(next > 0) - u64::from(counts[s] > 0);
            created += u64::from(next > 0 && counts[s] == 0);
        }
        check_instances(next_live, limits)?;
        for &(s, d) in changes {
            counts[s] = counts[s].checked_add_signed(d).expect("validated update");
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

    #[allow(clippy::cast_precision_loss)]
    pub fn channel(
        &mut self,
        net: &Network,
        i: usize,
        counts: &[u64],
        env: &Env<'_>,
    ) -> Result<f64, Halt> {
        let c = &net.plans[i];
        let [s, t] = c.reactants;
        let weight = match c.arity {
            0 => 1.0,
            1 => counts[s] as f64,
            _ if s == t => counts[s] as f64 * counts[s].saturating_sub(1) as f64 * 0.5,
            _ => counts[s] as f64 * counts[t] as f64,
        };
        if weight == 0.0 {
            return Ok(0.0);
        }
        let global = match c.factor {
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
        };
        let value = global * (weight * c.weight);
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(net.rate_error(i))
        }
    }
}

pub(crate) fn check_instances(live: u64, limits: &Limits) -> Result<(), Halt> {
    if live <= limits.max_instances {
        return Ok(());
    }
    Err(Halt::Limit(
        crate::Termination::InstanceLimit,
        crn_craft::diagnostic::Diagnostic::error(
            "CRN-R042",
            "/budget/max_instances",
            "instance limit reached",
            "Increase max_instances.",
        ),
    ))
}

fn normalize(x: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x }
}

pub(crate) fn direct<R: RngCore>(
    rates: &[f64],
    time: f64,
    words: &mut Words<'_, R>,
) -> Result<Option<(f64, usize)>, Halt> {
    direct_with(rates, time, words, |words, total| {
        // Open uniform: u = 0 would give a zero wait. Other words are unchanged.
        Ok(time - libm::log(1.0 - words.open()?) / total)
    })
}

pub(crate) fn direct_with<R: RngCore>(
    rates: &[f64],
    time: f64,
    words: &mut Words<'_, R>,
    next_time: impl FnOnce(&mut Words<'_, R>, f64) -> Result<f64, Halt>,
) -> Result<Option<(f64, usize)>, Halt> {
    let total: f64 = rates.iter().sum();
    if total == 0.0 {
        return Ok(None);
    }
    if !total.is_finite()
        || rates
            .iter()
            .any(|&a| a > 0.0 && a < total * crate::engine::RESOLUTION)
    {
        // Static channels are selected in one flat level; the tree engine
        // resolves rules and instances separately and can still run.
        if total.is_finite() {
            return Err(Halt::Numerical(crn_craft::diagnostic::Diagnostic::error(
                "CRN-R031",
                "/rules",
                "a positive channel is below 53-bit resolution of the total",
                "Use the tree method or rescale rates.",
            )));
        }
        return Err(numerical(
            "/rules",
            "unresolvable total propensity",
            "Rescale rates or populations.",
        ));
    }
    let next = next_time(words, total)?;
    if next <= time {
        return Err(numerical(
            "/time",
            "clock failed to advance",
            "Rescale rates and time.",
        ));
    }
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
            c.reactants.len() == 1
                && c.changes.len() == 1
                && (c.changes[0] == (c.reactants[0], -1) || c.changes[0].1 == 1)
        })
    }
}
