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
        Ok(Self {
            species,
            channels,
            observables,
            rates: model
                .rules
                .iter()
                .map(|r| r.rate.kernel(&model.params))
                .collect(),
            affected,
        })
    }

    pub fn initial(&self) -> Vec<u64> {
        self.species.iter().map(|s| s.initial).collect()
    }

    #[allow(clippy::cast_precision_loss)]
    pub fn refresh(&self, counts: &[u64], values: &mut [f64]) -> Result<(), Halt> {
        for (i, reduce, terms) in &self.observables {
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

    #[allow(clippy::cast_precision_loss)]
    pub fn rate(&self, i: usize, counts: &[u64], env: &Env<'_>) -> Result<f64, Halt> {
        let c = &self.channels[i];
        let weight = match c.reactants.as_slice() {
            [] => 1.0,
            &[s] => counts[s] as f64,
            &[s, t] if s == t => counts[s] as f64 * counts[s].saturating_sub(1) as f64 * 0.5,
            &[s, t] => counts[s] as f64 * counts[t] as f64,
            _ => unreachable!(),
        };
        if weight == 0.0 {
            return Ok(0.0);
        }
        let value = self.rates[c.rule].eval(env) * (weight * c.factor);
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(numerical(
                &format!("/rules/{}/rate", c.rule),
                "invalid fixed-channel propensity",
                "Rates must be finite and nonnegative.",
            ))
        }
    }

    pub fn apply(&self, i: usize, counts: &mut [u64], limits: &Limits) -> Result<(), Halt> {
        let changes = &self.channels[i].changes;
        let mut live = counts.iter().filter(|&&n| n > 0).count() as u64;
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
            live = live + u64::from(next > 0) - u64::from(counts[s] > 0);
        }
        check_instances(live, limits)?;
        for &(s, d) in changes {
            counts[s] = counts[s].checked_add_signed(d).expect("validated update");
        }
        Ok(())
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
    let total: f64 = rates.iter().sum();
    if total == 0.0 {
        return Ok(None);
    }
    if !total.is_finite()
        || rates
            .iter()
            .any(|&a| a > 0.0 && a < total * crate::engine::RESOLUTION)
    {
        return Err(numerical(
            "/rules",
            "unresolvable total propensity",
            "Rescale rates or populations.",
        ));
    }
    let next = time - libm::log(1.0 - words.uniform()?) / total;
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
