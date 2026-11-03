//! Preparation-time recognition for the three-class fixed-population sampler.
use crate::engine::RESOLUTION;
use crate::expr::RateKernel;
use crate::flat::Network;
use crate::model::{Compiled, Reduce};

const DEATH: [usize; 3] = [0, 4, 6];
const BIRTH: [usize; 3] = [1, 5, 7];
const MUTATION: [usize; 2] = [2, 3];

pub(crate) struct Kernel {
    pub capacity: u64,
    pub delta: f64,
    pub weights: [f64; 3],
    pub birth: [f64; 3],
    pub wt_unchanged: f64,
    pub wt_through_neutral: f64,
    pub threshold: Option<f64>,
    pub channels: [usize; 8],
    pub species: [usize; 3],
    pub count_trackers: Option<[Vec<usize>; 3]>,
}

// Slots ordered by the event representation in the specialized loop.
const SHAPES: [(usize, usize, i64); 8] = [
    (0, 0, -1),
    (0, 0, 1),
    (0, 1, 1),
    (0, 2, 1),
    (1, 1, -1),
    (1, 1, 1),
    (2, 2, -1),
    (2, 2, 1),
];

impl Kernel {
    pub(crate) fn compile(
        net: &Network,
        model: &Compiled,
        spec: &crate::observe::Spec,
    ) -> Option<Self> {
        if net.species.len() != 3 || net.channels.len() != SHAPES.len() {
            return None;
        }
        // Identify the parent of both mutations, then the two destinations in
        // static species order. This is independent of model and rule names.
        let wt = net
            .channels
            .iter()
            .find_map(|c| match (&c.reactants[..], &c.changes[..]) {
                ([parent], [(target, 1)]) if parent != target => Some(*parent),
                _ => None,
            })?;
        let mut order = (0..3).filter(|&i| i != wt);
        let species = [wt, order.next()?, order.next()?];
        let mut channels = [0; 8];
        for (slot, (parent, target, change)) in SHAPES.into_iter().enumerate() {
            channels[slot] = net.channels.iter().position(|c| {
                c.reactants == [species[parent]]
                    && c.changes == [(species[target], change)]
                    && c.factor == 1.0
            })?;
        }
        let rate = |slot: usize| &net.rates[net.channels[channels[slot]].rule];
        let RateKernel::Constant(delta) = *rate(DEATH[0]) else {
            return None;
        };
        if !delta.is_finite()
            || delta < 0.0
            || DEATH
                .into_iter()
                .any(|i| !matches!(rate(i), RateKernel::Constant(d) if *d == delta))
        {
            return None;
        }
        let RateKernel::Vacancy {
            limit,
            count,
            denominator,
            ..
        } = *rate(BIRTH[0])
        else {
            return None;
        };
        if !(1.0..=1_000_000.0).contains(&limit)
            || limit.fract() != 0.0
            || net.species.iter().map(|s| s.initial).sum::<u64>() > limit as u64
        {
            return None;
        }
        let projection = |index| {
            let (_, reduce, terms) = net.observables.iter().find(|(i, _, _)| *i == index)?;
            if !matches!(reduce, Reduce::Sum) {
                return None;
            }
            Some(species.map(|s| {
                terms
                    .iter()
                    .find(|&&(i, _)| i == s)
                    .map_or(0.0, |&(_, w)| w)
            }))
        };
        if projection(count)? != [1.0; 3] {
            return None;
        }
        let weights = projection(denominator)?;
        if weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return None;
        }
        // Certify finite report projections even when they are not needed at
        // each event. Leave more general reductions to Static Direct.
        for &(i, _, _) in &net.observables {
            if !(projection(i)?.iter().map(|w| w.abs()).sum::<f64>() * limit).is_finite() {
                return None;
            }
        }
        let threshold = if let Some(stop) = &model.stop {
            let (mutants, _, _) = net.observables.iter().find(|(i, _, _)| {
                projection(*i) == Some([0.0, 1.0, 1.0]) && stop.fraction_stop(count, *i).is_some()
            })?;
            Some(stop.fraction_stop(count, *mutants)?)
        } else {
            None
        };
        let den_min = weights.into_iter().fold(f64::INFINITY, f64::min);
        let den_max = limit * weights.into_iter().fold(0.0, f64::max);
        let low = 1.0 / den_max;
        let high = limit / den_min;
        if low < f64::MIN_POSITIVE || high > f64::MAX / 16.0 || den_max > f64::MAX / 16.0 {
            return None;
        }
        let mut lower = if delta > 0.0 { delta } else { f64::INFINITY };
        let mut coefficient = |slot| {
            let RateKernel::Vacancy {
                limit: k,
                count: n,
                denominator: d,
                factors,
            } = rate(slot)
            else {
                return None;
            };
            if (*k, *n, *d) != (limit, count, denominator) {
                return None;
            }
            let mut bounds = [low, high, 1.0];
            for &factor in factors {
                if !factor.is_finite() || factor < 0.0 {
                    return None;
                }
                let positive = factor > 0.0 && bounds[0] > 0.0;
                bounds = bounds.map(|x| x * factor);
                if !bounds[1].is_finite()
                    || !bounds[2].is_finite()
                    || (positive && bounds[0].min(bounds[2]) < f64::MIN_POSITIVE)
                {
                    return None;
                }
            }
            if bounds[0] > 0.0 {
                lower = lower.min(bounds[0]);
            }
            if !(bounds[1] * limit).is_finite() {
                return None;
            }
            Some(bounds[2])
        };
        let plain = coefficient(BIRTH[0])?;
        let neutral = coefficient(MUTATION[0])?;
        let advantaged = coefficient(MUTATION[1])?;
        let birth = [
            plain + neutral + advantaged,
            coefficient(BIRTH[1])?,
            coefficient(BIRTH[2])?,
        ];
        let max_birth = birth
            .into_iter()
            .zip(weights)
            .map(|(rate, weight)| rate / weight)
            .fold(0.0, f64::max);
        let upper = limit * (delta + max_birth);
        // Reserve rounding margin for grouped cumulative thresholds and for
        // reordering the folded products. Very stiff candidates fall back.
        if !upper.is_finite()
            || lower < (16.0 * RESOLUTION * upper).max(f64::MIN_POSITIVE)
            || !birth[0].is_finite()
        {
            return None;
        }
        Some(Self {
            capacity: limit as u64,
            delta,
            weights,
            birth,
            wt_unchanged: plain,
            wt_through_neutral: plain + neutral,
            threshold,
            channels,
            species,
            count_trackers: spec.count_trackers(net),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Env;
    use crate::observe::Spec;

    #[test]
    fn grouped_hazards_match_every_original_channel() {
        let source: crate::Model = serde_json::from_str(include_str!(
            "../../../models/rules/mtdna-three-class.rules.json"
        ))
        .unwrap();
        let model = crate::model::compile(&source).unwrap();
        let net = Network::compile(&model).unwrap();
        let spec = Spec::compile(&source, &model, &[], None).unwrap();
        let kernel = Kernel::compile(&net, &model, &spec).unwrap();
        for classes in [
            [1000, 0, 0],
            [999, 0, 0],
            [4, 0, 5],
            [0, 1, 0],
            [13, 71, 109],
        ] {
            let mut counts = [0; 3];
            for (slot, count) in kernel.species.into_iter().zip(classes) {
                counts[slot] = count;
            }
            let mut values = vec![0.0; model.observables.len()];
            net.refresh(&counts, &mut values).ok().unwrap();
            let env = Env {
                params: &model.params,
                observables: &values,
                attributes: [&[], &[]],
            };
            let [w, n, a] = classes.map(|v| v as f64);
            let [ww, nw, aw] = kernel.weights;
            let scale = (kernel.capacity as f64 - w - n - a) / (w * ww + n * nw + a * aw);
            let grouped = [
                w * kernel.delta,
                w * scale * kernel.wt_unchanged,
                w * scale * (kernel.wt_through_neutral - kernel.wt_unchanged),
                w * scale * (kernel.birth[0] - kernel.wt_through_neutral),
                n * kernel.delta,
                n * scale * kernel.birth[1],
                a * kernel.delta,
                a * scale * kernel.birth[2],
            ];
            for (rate, channel) in grouped.into_iter().zip(kernel.channels) {
                let original = net.rate(channel, &counts, &env).ok().unwrap();
                assert!(
                    (rate - original).abs() <= 1e-10 * original,
                    "{rate} != {original}"
                );
            }
        }
    }
}
