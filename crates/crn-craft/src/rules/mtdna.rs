//! Preparation-time recognition for the three-class fixed-population sampler.
use crate::rules::engine::RESOLUTION;
use crate::rules::expr::{Env, Expr, Func, Node, Op, RateKernel, Slot};
use crate::rules::flat::Network;
use crate::rules::model::{Compiled, Reduce};

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

const TOPOLOGY: &str = "the model is not three classes with the eight fixed-population channels";
const DEGRADATION: &str = "degradation is not one constant rate shared by all classes";
const VACANCY: &str = "replication does not use one shared vacancy law max(0, K - N) / D";
const CAPACITY: &str = "K is not an integer from 1 to 1,000,000 at or above the initial population";
const PROJECTIONS: &str = "observables are not finite weighted sums of the three classes";
const STOP: &str = "stop is neither absent nor the guarded mutant fraction";
const RANGE: &str = "hazards span more than the certified floating-point range";
const WHOLE: &str = "the model has whole-rule propensities";

/// The law of one channel, as the certificate reads it from the model.
enum Law {
    /// `delta` per molecule of the parent class.
    Death(f64),
    /// `max(0, limit - N) / D` times the factors, multiplied in order, per
    /// molecule of the parent class; `N` is the total count and `D` the sum
    /// of the counts with `weights`.
    Birth {
        limit: f64,
        weights: [f64; 3],
        factors: Vec<f64>,
    },
}

/// An exact linear form over the three classes (in the kernel's class
/// order): a constant and one coefficient per class count.
#[derive(Clone, Copy, PartialEq)]
struct Linear {
    constant: f64,
    counts: [f64; 3],
}

/// `a + b`, if it is exact (the TwoSum error is zero and the sum finite).
fn exact_sum(a: f64, b: f64) -> Option<f64> {
    let s = a + b;
    let bb = s - a;
    let error = (a - (s - bb)) + (b - bb);
    (s.is_finite() && error == 0.0).then_some(s)
}

/// `a * b`, if it is exact (the fused residual is zero and the product
/// finite).
fn exact_product(a: f64, b: f64) -> Option<f64> {
    let p = a * b;
    (p.is_finite() && a.mul_add(b, -p) == 0.0).then_some(p)
}

/// Reads whole-rule propensities and stops of the builder form: a model
/// whose classes are attribute-free species, whose propensities write the
/// vacancy law and the parent count, and whose stop writes the mutant
/// fraction over sums of counts.
struct Reader<'a> {
    params: &'a [f64],
    /// The weights of each Sum projection over the classes, by observable.
    projections: Vec<Option<[f64; 3]>>,
}

impl Reader<'_> {
    fn constant(&self, node: &Node) -> Option<f64> {
        let expr = Expr { root: node.clone() };
        expr.parameter_only().then(|| {
            expr.eval(&Env {
                params: self.params,
                observables: &[],
                attributes: [&[], &[]],
            })
        })
    }

    /// The exact linear form of `node`: sums, differences, negations and
    /// products with constants of observables that are weighted sums of the
    /// classes, every coefficient computed exactly.
    fn linear(&self, node: &Node) -> Option<Linear> {
        if let Some(c) = self.constant(node) {
            return c.is_finite().then_some(Linear {
                constant: c,
                counts: [0.0; 3],
            });
        }
        let combine = |a: Linear, b: Linear, sign: f64| -> Option<Linear> {
            let mut counts = [0.0; 3];
            for (k, count) in counts.iter_mut().enumerate() {
                *count = exact_sum(a.counts[k], sign * b.counts[k])?;
            }
            Some(Linear {
                constant: exact_sum(a.constant, sign * b.constant)?,
                counts,
            })
        };
        let scale = |a: Linear, c: f64| -> Option<Linear> {
            let mut counts = [0.0; 3];
            for (count, &a) in counts.iter_mut().zip(&a.counts) {
                *count = exact_product(a, c)?;
            }
            Some(Linear {
                constant: exact_product(a.constant, c)?,
                counts,
            })
        };
        match node {
            Node::Var(Slot::Observable(i)) => Some(Linear {
                constant: 0.0,
                counts: (*self.projections.get(*i)?)?,
            }),
            Node::Neg(a) => scale(self.linear(a)?, -1.0),
            Node::Bin(Op::Add, a, b) => combine(self.linear(a)?, self.linear(b)?, 1.0),
            Node::Bin(Op::Sub, a, b) => combine(self.linear(a)?, self.linear(b)?, -1.0),
            Node::Bin(Op::Mul, a, b) => match (self.constant(a), self.constant(b)) {
                (Some(c), None) => scale(self.linear(b)?, c),
                (None, Some(c)) => scale(self.linear(a)?, c),
                _ => None,
            },
            _ => None,
        }
    }

    /// The factors of a product: `a * b` and `a / b` in any association,
    /// numerator factors first (in order), then denominator factors.
    fn factors<'n>(node: &'n Node, top: &mut Vec<&'n Node>, bottom: &mut Vec<&'n Node>) {
        match node {
            Node::Bin(Op::Mul, a, b) => {
                Self::factors(a, top, bottom);
                Self::factors(b, top, bottom);
            }
            Node::Bin(Op::Div, a, b) => {
                Self::factors(a, top, bottom);
                Self::factors(b, bottom, top);
            }
            _ => top.push(node),
        }
    }

    /// The law of a whole-rule propensity of a channel whose parent class
    /// is `parent`: constants times the parent's count (a death), or the
    /// vacancy law `max(0, K - N) / D` times constants times the parent's
    /// count (a birth). Equal to the propensity as a real-valued function of
    /// the counts; `None` for any other propensity.
    fn whole(&self, expr: &Expr, parent: usize) -> Option<Law> {
        let (mut top, mut bottom) = (Vec::new(), Vec::new());
        Self::factors(&expr.root, &mut top, &mut bottom);
        let mut parent_count = Linear {
            constant: 0.0,
            counts: [0.0; 3],
        };
        parent_count.counts[parent] = 1.0;
        let mut factors = Vec::new();
        let mut counted = false;
        let mut vacancy = None;
        for node in top {
            if let Some(c) = self.constant(node) {
                factors.push(c);
            } else if !counted && self.linear(node) == Some(parent_count) {
                counted = true;
            } else if vacancy.is_none() {
                // max(0, K - N) or max(K - N, 0).
                let Node::Call(Func::Max, args) = node else {
                    return None;
                };
                let free = match (self.constant(&args[0]), self.constant(&args[1])) {
                    (Some(z), None) if z.to_bits() == 0f64.to_bits() => &args[1],
                    (None, Some(z)) if z.to_bits() == 0f64.to_bits() => &args[0],
                    _ => return None,
                };
                let free = self.linear(free)?;
                if free.counts != [-1.0; 3] {
                    return None;
                }
                vacancy = Some(free.constant);
            } else {
                return None;
            }
        }
        if !counted {
            return None;
        }
        match (vacancy, bottom.as_slice()) {
            (None, []) => {
                // A death: one constant, the rate per molecule.
                let [delta] = factors[..] else {
                    return None;
                };
                Some(Law::Death(delta))
            }
            (Some(limit), [denominator]) => {
                let d = self.linear(denominator)?;
                (d.constant.to_bits() == 0f64.to_bits()).then_some(Law::Birth {
                    limit,
                    weights: d.counts,
                    factors,
                })
            }
            _ => None,
        }
    }

    /// The fraction `f` of a stop `T > 0 && P >= f * T` (or `T * f`), with
    /// `T` a sum of observables that are integer-weighted sums of the
    /// classes and equal to the total count, and `P` such a sum equal to the
    /// mutant count. Every partial sum is then an integer of at most
    /// `bound` in magnitude, so the stop computes the total and the mutant
    /// count exactly, as the fast path does.
    fn stop(&self, expr: &Expr, bound: f64) -> Option<f64> {
        let Node::Bin(Op::And, guard, test) = &expr.root else {
            return None;
        };
        let Node::Bin(Op::Gt, total, zero) = guard.as_ref() else {
            return None;
        };
        if self.constant(zero).map(f64::to_bits) != Some(0f64.to_bits()) {
            return None;
        }
        let Node::Bin(Op::Ge, part, product) = test.as_ref() else {
            return None;
        };
        let Node::Bin(Op::Mul, a, b) = product.as_ref() else {
            return None;
        };
        let (fraction, again) = match (self.constant(a), self.constant(b)) {
            (Some(f), None) => (f, b.as_ref()),
            (None, Some(f)) => (f, a.as_ref()),
            _ => return None,
        };
        if !(0.0..=1.0).contains(&fraction) {
            return None;
        }
        let integral = |node: &Node, expected: [f64; 3]| -> Option<()> {
            let mut weight = 0.0;
            if !self.integer_sum(node, &mut weight) || weight * bound > 9_007_199_254_740_992.0 {
                return None;
            }
            (self.linear(node)?
                == Linear {
                    constant: 0.0,
                    counts: expected,
                })
            .then_some(())
        };
        integral(total, [1.0; 3])?;
        integral(again, [1.0; 3])?;
        integral(part, [0.0, 1.0, 1.0])?;
        Some(fraction)
    }

    /// Whether `node` is a sum of observables with integer weights; adds the
    /// sum of their absolute weights to `weight`.
    fn integer_sum(&self, node: &Node, weight: &mut f64) -> bool {
        match node {
            Node::Var(Slot::Observable(i)) => {
                let Some(Some(w)) = self.projections.get(*i) else {
                    return false;
                };
                *weight += w.iter().map(|x| x.abs()).sum::<f64>();
                w.iter().all(|x| x.fract() == 0.0)
            }
            Node::Bin(Op::Add, a, b) => self.integer_sum(a, weight) && self.integer_sum(b, weight),
            _ => false,
        }
    }
}

impl Kernel {
    /// Recognize the fixed three-class model, or say why Static Direct runs.
    ///
    /// The model file's form (per-instance rates over the classes of one
    /// attributed type, the vacancy law through the observables `N` and
    /// `D`) is recognized for every request version. With `builder` (a
    /// request of `crn-craft.run/1.0`) a model the file form refuses is read
    /// again as the builder form ([`Self::builder_form`]); a model without
    /// whole-rule propensities that both refuse keeps the file form's
    /// reason.
    pub(crate) fn compile(
        net: &Network,
        model: &Compiled,
        spec: &crate::rules::observe::Spec,
        builder: bool,
    ) -> Result<Self, &'static str> {
        match Self::file_form(net, model, spec) {
            Err(reason) if builder => {
                let whole = net.channels.iter().any(|c| c.whole);
                Self::builder_form(net, model, spec)
                    .map_err(|other| if whole { other } else { reason })
            }
            result => result,
        }
    }

    /// The parent and the two destinations of the mutations in static
    /// species order, and the channel of each slot of [`SHAPES`].
    fn topology(net: &Network) -> Result<([usize; 8], [usize; 3]), &'static str> {
        if net.species.len() != 3 || net.channels.len() != SHAPES.len() {
            return Err(TOPOLOGY);
        }
        // Identify the parent of both mutations, then the two destinations in
        // static species order. This is independent of model and rule names.
        let wt = net
            .channels
            .iter()
            .find_map(|c| match (&c.reactants[..], &c.changes[..]) {
                ([parent], [(target, 1)]) if parent != target => Some(*parent),
                _ => None,
            })
            .ok_or(TOPOLOGY)?;
        let mut order = (0..3).filter(|&i| i != wt);
        let species = [
            wt,
            order.next().ok_or(TOPOLOGY)?,
            order.next().ok_or(TOPOLOGY)?,
        ];
        let mut channels = [0; 8];
        for (slot, (parent, target, change)) in SHAPES.into_iter().enumerate() {
            channels[slot] = net
                .channels
                .iter()
                .position(|c| {
                    c.reactants == [species[parent]]
                        && c.changes == [(species[target], change)]
                        && c.factor == 1.0
                })
                .ok_or(TOPOLOGY)?;
        }
        Ok((channels, species))
    }

    /// The model file's form, as 0.7.0 recognized it.
    fn file_form(
        net: &Network,
        model: &Compiled,
        spec: &crate::rules::observe::Spec,
    ) -> Result<Self, &'static str> {
        if net.channels.iter().any(|c| c.whole) {
            return Err(WHOLE);
        }
        let (channels, species) = Self::topology(net)?;
        let rate = |slot: usize| &net.rates[net.channels[channels[slot]].rule];
        let RateKernel::Constant(delta) = *rate(DEATH[0]) else {
            return Err(DEGRADATION);
        };
        if !delta.is_finite()
            || delta < 0.0
            || DEATH
                .into_iter()
                .any(|i| !matches!(rate(i), RateKernel::Constant(d) if *d == delta))
        {
            return Err(DEGRADATION);
        }
        let RateKernel::Vacancy {
            limit,
            count,
            denominator,
            ..
        } = *rate(BIRTH[0])
        else {
            return Err(VACANCY);
        };
        if !(1.0..=1_000_000.0).contains(&limit)
            || limit.fract() != 0.0
            || net.species.iter().map(|s| s.initial).sum::<u64>() > limit as u64
        {
            return Err(CAPACITY);
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
        if projection(count).ok_or(VACANCY)? != [1.0; 3] {
            return Err(VACANCY);
        }
        let weights = projection(denominator).ok_or(VACANCY)?;
        if weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return Err(VACANCY);
        }
        // Certify finite report projections even when they are not needed at
        // each event. Leave more general reductions to Static Direct.
        for &(i, _, _) in &net.observables {
            let weights = projection(i).ok_or(PROJECTIONS)?;
            if !(weights.iter().map(|w| w.abs()).sum::<f64>() * limit).is_finite() {
                return Err(PROJECTIONS);
            }
        }
        let threshold = if let Some(stop) = &model.stop {
            let (mutants, _, _) = net
                .observables
                .iter()
                .find(|(i, _, _)| {
                    projection(*i) == Some([0.0, 1.0, 1.0])
                        && stop.fraction_stop(count, *i).is_some()
                })
                .ok_or(STOP)?;
            Some(stop.fraction_stop(count, *mutants).ok_or(STOP)?)
        } else {
            None
        };
        let factors = |slot: usize| match rate(slot) {
            RateKernel::Vacancy {
                limit: k,
                count: n,
                denominator: d,
                factors,
            } if (*k, *n, *d) == (limit, count, denominator) => Ok(factors.as_slice()),
            _ => Err(VACANCY),
        };
        Self::certify(
            net, spec, channels, species, delta, limit, weights, threshold, factors,
        )
    }

    /// The builder form: the same chemistry written with attribute-free
    /// classes and whole-rule propensities, with or without the parameters
    /// inlined, and a takeover stop over sums of counts. Each channel's law
    /// is read from its rule ([`Reader::whole`] for a whole-rule propensity,
    /// the rate kernel otherwise); the laws must be the file form's: one
    /// death rate, one vacancy law with one capacity and one weight vector.
    fn builder_form(
        net: &Network,
        model: &Compiled,
        spec: &crate::rules::observe::Spec,
    ) -> Result<Self, &'static str> {
        let (channels, species) = Self::topology(net)?;
        // A whole-rule channel of one molecule: its propensity applies while
        // the parent class is present.
        if channels.iter().any(|&c| net.channels[c].multi.is_some()) {
            return Err(TOPOLOGY);
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
        let reader = Reader {
            params: &model.params,
            projections: (0..model.observables.len()).map(projection).collect(),
        };
        let law = |slot: usize| -> Result<Law, &'static str> {
            let channel = &net.channels[channels[slot]];
            let birth = SHAPES[slot].2 == 1;
            let wrong = if birth { VACANCY } else { DEGRADATION };
            if channel.whole {
                return reader
                    .whole(&model.rules[channel.rule].rate, SHAPES[slot].0)
                    .ok_or(wrong);
            }
            match &net.rates[channel.rule] {
                &RateKernel::Constant(d) if !birth => Ok(Law::Death(d)),
                RateKernel::Vacancy {
                    limit,
                    count,
                    denominator,
                    factors,
                } if birth && projection(*count) == Some([1.0; 3]) => Ok(Law::Birth {
                    limit: *limit,
                    weights: projection(*denominator).ok_or(VACANCY)?,
                    factors: factors.clone(),
                }),
                _ => Err(wrong),
            }
        };
        let laws = (0..SHAPES.len()).map(law).collect::<Result<Vec<_>, _>>()?;
        let Law::Death(delta) = laws[DEATH[0]] else {
            return Err(DEGRADATION);
        };
        if !delta.is_finite()
            || delta < 0.0
            || DEATH
                .into_iter()
                .any(|i| !matches!(laws[i], Law::Death(d) if d == delta))
        {
            return Err(DEGRADATION);
        }
        let Law::Birth { limit, weights, .. } = laws[BIRTH[0]] else {
            return Err(VACANCY);
        };
        if !(1.0..=1_000_000.0).contains(&limit)
            || limit.fract() != 0.0
            || net.species.iter().map(|s| s.initial).sum::<u64>() > limit as u64
        {
            return Err(CAPACITY);
        }
        if weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return Err(VACANCY);
        }
        // Certify finite report projections even when they are not needed at
        // each event. Leave more general reductions to Static Direct.
        for &(i, _, _) in &net.observables {
            let weights = projection(i).ok_or(PROJECTIONS)?;
            if !(weights.iter().map(|w| w.abs()).sum::<f64>() * limit).is_finite() {
                return Err(PROJECTIONS);
            }
        }
        let threshold = match &model.stop {
            Some(stop) => {
                // The file form's stop over the observables N and mutants,
                // or the stop over sums of counts.
                let file = net.observables.iter().find_map(|&(n, _, _)| {
                    let (m, _, _) = net.observables.iter().find(|(m, _, _)| {
                        projection(n) == Some([1.0; 3])
                            && projection(*m) == Some([0.0, 1.0, 1.0])
                            && stop.fraction_stop(n, *m).is_some()
                    })?;
                    stop.fraction_stop(n, *m)
                });
                Some(file.or_else(|| reader.stop(stop, limit)).ok_or(STOP)?)
            }
            None => None,
        };
        let factors = |slot: usize| match &laws[slot] {
            Law::Birth {
                limit: k,
                weights: w,
                factors,
            } if k.to_bits() == limit.to_bits()
                && w.map(f64::to_bits) == weights.map(f64::to_bits) =>
            {
                Ok(factors.as_slice())
            }
            _ => Err(VACANCY),
        };
        Self::certify(
            net, spec, channels, species, delta, limit, weights, threshold, factors,
        )
    }

    /// The checks and constants that both forms share, from the death rate,
    /// the vacancy law and the stop's threshold; `factors` gives the factors
    /// of a birth slot's law, multiplied in order.
    #[allow(clippy::too_many_arguments)]
    fn certify<'f>(
        net: &Network,
        spec: &crate::rules::observe::Spec,
        channels: [usize; 8],
        species: [usize; 3],
        delta: f64,
        limit: f64,
        weights: [f64; 3],
        threshold: Option<f64>,
        factors: impl Fn(usize) -> Result<&'f [f64], &'static str>,
    ) -> Result<Self, &'static str> {
        let den_min = weights.into_iter().fold(f64::INFINITY, f64::min);
        let den_max = limit * weights.into_iter().fold(0.0, f64::max);
        let low = 1.0 / den_max;
        let high = limit / den_min;
        if low < f64::MIN_POSITIVE || high > f64::MAX / 16.0 || den_max > f64::MAX / 16.0 {
            return Err(RANGE);
        }
        // Smallest positive coefficient of each replication outcome.
        let mut smallest = f64::INFINITY;
        let mut coefficient = |slot| {
            let factors = factors(slot)?;
            let mut bounds = [low, high, 1.0];
            for &factor in factors {
                if !factor.is_finite() || factor < 0.0 {
                    return Err(RANGE);
                }
                let positive = factor > 0.0 && bounds[0] > 0.0;
                bounds = bounds.map(|x| x * factor);
                if !bounds[1].is_finite()
                    || !bounds[2].is_finite()
                    || (positive && bounds[0].min(bounds[2]) < f64::MIN_POSITIVE)
                {
                    return Err(RANGE);
                }
            }
            if bounds[2] > 0.0 {
                smallest = smallest.min(bounds[2]);
            }
            if !(bounds[1] * limit).is_finite() {
                return Err(RANGE);
            }
            Ok(bounds[2])
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
        let max_weight = weights.into_iter().fold(0.0, f64::max);
        // Resolution is a property of one state: every positive hazard must be
        // at least 16 * 2^-53 of the total in the same state. With N in the
        // population, V = K - N >= 1 vacancies and D <= max_weight * N, a
        // class with count n has replication hazard n * c * V / D and the
        // total is at most N * delta + V * max_birth. The ratio for n = 1 is
        // smallest at V = 1, so c / (K * max_weight * (K * delta + max_birth))
        // bounds every replication outcome. A degradation hazard n * delta
        // against a total of at most K * max(delta, max_birth) gives the
        // other bound. The margin covers grouped thresholds and reordered
        // products.
        let total_rate = limit * delta + max_birth;
        let floor = 16.0 * RESOLUTION;
        let births_resolved =
            !smallest.is_finite() || smallest / (limit * max_weight * total_rate) >= floor;
        let deaths_resolved = delta == 0.0 || delta / (limit * delta.max(max_birth)) >= floor;
        if !total_rate.is_finite()
            || !(limit * delta.max(max_birth)).is_finite()
            || !births_resolved
            || !deaths_resolved
            || !birth[0].is_finite()
        {
            return Err(RANGE);
        }
        Ok(Self {
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
    use crate::rules::expr::Env;
    use crate::rules::observe::Spec;

    /// The kernel of `source`, or the reason it is refused.
    fn kernel(source: &serde_json::Value, builder: bool) -> Result<Kernel, &'static str> {
        let source: crate::rules::Model = serde_json::from_value(source.clone()).unwrap();
        let model = crate::rules::model::compile(&source).unwrap();
        let net = Network::compile(&model).unwrap();
        let spec = Spec::compile(&source, &model, &[], None).unwrap();
        Kernel::compile(&net, &model, &spec, builder)
    }

    /// Every channel's grouped hazard equals the channel's propensity as
    /// Static Direct evaluates it, at states across the population range.
    fn hazards_match(source: &serde_json::Value, builder: bool) {
        let parsed: crate::rules::Model = serde_json::from_value(source.clone()).unwrap();
        let model = crate::rules::model::compile(&parsed).unwrap();
        let net = Network::compile(&model).unwrap();
        let kernel = kernel(source, builder).unwrap();
        for classes in [
            [1000, 0, 0],
            [999, 0, 0],
            [4, 0, 5],
            [0, 1, 0],
            [13, 71, 109],
            [0, 0, 1000],
            [500, 250, 250],
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
                    "{rate} != {original} at {classes:?}"
                );
            }
        }
    }

    fn file_model() -> serde_json::Value {
        serde_json::from_str(include_str!(
            "../../../../models/rules/mtdna-three-class.rules.json"
        ))
        .unwrap()
    }

    /// The three-class model as `three_class()` of guide chapter 12 writes
    /// it (`python/crn_craft` builder, `model/1.0`), with species `names`,
    /// `free` the vacancy term and `stop` the stop, if any.
    fn builder_model(names: [&str; 3], free: &str, stop: Option<&str>) -> serde_json::Value {
        let [w, m, a] = names;
        let rule = |name: &str, parent: &str, propensity: String, products: [&str; 2]| {
            serde_json::json!({
                "name": name,
                "reactants": [{"label": "x", "type": parent}],
                "propensity": propensity,
                "products": products.map(|p| serde_json::json!({"type": p})),
            })
        };
        let death = |name: &str, parent: &str| {
            serde_json::json!({
                "name": name,
                "reactants": [{"label": "x", "type": parent}],
                "propensity": format!("delta * {parent}"),
                "products": [],
            })
        };
        let mut model = serde_json::json!({
            "schema": "crn-craft.model/1.0",
            "name": "mtdna_three_class",
            "parameters": {"K": 1000.0, "delta": 0.0028881132523331052, "sw": 1.0,
                           "sm": 1.5, "mu": 3e-05, "fracAdv": 0.2},
            "types": [{"name": w}, {"name": m}, {"name": a}],
            "initial": [{"type": w, "count": 1000}],
            "observables": {w: {"type": w}, m: {"type": m}, a: {"type": a}},
            "rules": [
                rule("replication_wt", w, format!("{free} * sw * (1 - mu) * {w}"), [w, w]),
                rule("mutation_mt", w, format!("{free} * sw * mu * (1 - fracAdv) * {w}"), [w, m]),
                rule("mutation_mtAdv", w, format!("{free} * sw * mu * fracAdv * {w}"), [w, a]),
                rule("replication_mt", m, format!("{free} * sw * {m}"), [m, m]),
                rule("replication_mtAdv", a, format!("{free} * sm * {a}"), [a, a]),
                death("degradation_wt", w),
                death("degradation_mt", m),
                death("degradation_mtAdv", a),
            ],
        });
        if let Some(stop) = stop {
            model["stop"] = stop.into();
        }
        model
    }

    const FREE: &str = "max(0, K - W - M - A) / (W + sm * (M + A))";
    const TAKEOVER: &str = "W + M + A > 0 && M + A >= 0.9 * (W + M + A)";

    #[test]
    fn grouped_hazards_match_every_original_channel() {
        hazards_match(&file_model(), false);
    }

    /// E1: the builder form of guide chapter 12 is certified under
    /// `run/1.0`, with the file form's kernel to the bit (so its draws are
    /// the file form's), and refused for the `rules-run` versions.
    #[test]
    fn the_builder_form_is_certified_with_the_file_forms_kernel() {
        let builder = builder_model(["W", "M", "A"], FREE, Some(TAKEOVER));
        assert_eq!(kernel(&builder, false).err(), Some(WHOLE));
        let mine = kernel(&builder, true).unwrap();
        let mut file = file_model();
        file["stop"] = "N > 0 && mutants >= 0.9 * N".into();
        let theirs = kernel(&file, false).unwrap();
        assert_eq!(mine.capacity, theirs.capacity);
        let bits = |k: &Kernel| {
            [
                k.delta,
                k.weights[0],
                k.weights[1],
                k.weights[2],
                k.birth[0],
                k.birth[1],
                k.birth[2],
                k.wt_unchanged,
                k.wt_through_neutral,
                k.threshold.unwrap(),
            ]
            .map(f64::to_bits)
        };
        assert_eq!(bits(&mine), bits(&theirs));
        hazards_match(&builder, true);
        // The file form is unchanged under run/1.0.
        assert!(kernel(&file, true).is_ok());
    }

    /// E1: variants of the builder form that are the same law.
    #[test]
    fn builder_variants_that_are_the_same_law_are_certified() {
        let free_xyz = "max(0, K - X - Y - Z) / (X + sm * (Y + Z))";
        let inlined = "max(0, 1000 - W - M - A) / (W + 1.5 * (M + A))";
        let reordered = "max(K - (A + M + W), 0) / (sm * A + sm * M + W)";
        for (names, free, stop) in [
            (["W", "M", "A"], FREE, None),
            (["W", "M", "A"], FREE, Some(TAKEOVER)),
            (
                ["X", "Y", "Z"],
                free_xyz,
                Some("X + Y + Z > 0 && Y + Z >= 0.9 * (X + Y + Z)"),
            ),
            (
                ["W", "M", "A"],
                inlined,
                Some("W + M + A > 0 && M + A >= (W + M + A) * 0.5"),
            ),
            (
                ["W", "M", "A"],
                reordered,
                Some("A + W + M > 0 && A + M >= 1 * (W + M + A)"),
            ),
        ] {
            let model = builder_model(names, free, stop);
            assert!(kernel(&model, true).is_ok(), "{free} {stop:?}");
            hazards_match(&model, true);
        }
        // Species declared in another order, factors in another order and
        // association, the death rate written after the count.
        let mut model = builder_model(["W", "M", "A"], FREE, Some(TAKEOVER));
        model["types"] = serde_json::json!([{"name": "A"}, {"name": "W"}, {"name": "M"}]);
        model["rules"][0]["propensity"] =
            "W * (sw * (1 - mu)) * max(0, K - W - M - A) / (W + sm * (M + A))".into();
        model["rules"][5]["propensity"] = "W * delta".into();
        assert!(kernel(&model, true).is_ok());
        hazards_match(&model, true);
    }

    /// E1: near-misses of the builder form are refused, each with its
    /// reason.
    #[test]
    fn builder_near_misses_are_refused() {
        let refused = |model: &serde_json::Value| kernel(model, true).err();
        let base = || builder_model(["W", "M", "A"], FREE, Some(TAKEOVER));
        // A changed coefficient in one rule's vacancy weights.
        let mut model = base();
        model["rules"][3]["propensity"] =
            "max(0, K - W - M - A) / (W + 1.4 * (M + A)) * sw * M".into();
        assert_eq!(refused(&model), Some(VACANCY));
        // The vacancy without one class.
        let model = builder_model(
            ["W", "M", "A"],
            "max(0, K - W - M) / (W + sm * (M + A))",
            Some(TAKEOVER),
        );
        assert_eq!(refused(&model), Some(VACANCY));
        // Another capacity in one rule.
        let mut model = base();
        model["rules"][4]["propensity"] =
            "max(0, 999 - W - M - A) / (W + sm * (M + A)) * sm * A".into();
        assert_eq!(refused(&model), Some(VACANCY));
        // max(1, ...) is not the vacancy.
        let model = builder_model(
            ["W", "M", "A"],
            "max(1, K - W - M - A) / (W + sm * (M + A))",
            None,
        );
        assert_eq!(refused(&model), Some(VACANCY));
        // A replication without the parent's count, or with another count.
        let mut model = base();
        model["rules"][3]["propensity"] = format!("{FREE} * sw").into();
        assert_eq!(refused(&model), Some(VACANCY));
        let mut model = base();
        model["rules"][3]["propensity"] = format!("{FREE} * sw * W").into();
        assert_eq!(refused(&model), Some(VACANCY));
        // A replication squared in its count.
        let mut model = base();
        model["rules"][3]["propensity"] = format!("{FREE} * sw * M * M").into();
        assert_eq!(refused(&model), Some(VACANCY));
        // Another death rate for one class; a death without its count.
        let mut model = base();
        model["rules"][6]["propensity"] = "2 * delta * M".into();
        assert_eq!(refused(&model), Some(DEGRADATION));
        let mut model = base();
        model["rules"][7]["propensity"] = "delta".into();
        assert_eq!(refused(&model), Some(DEGRADATION));
        let mut model = base();
        model["rules"][6]["propensity"] = "0.0028881132523331057 * M".into();
        assert_eq!(refused(&model), Some(DEGRADATION));
        // A fourth reaction.
        let mut model = base();
        model["rules"]
            .as_array_mut()
            .unwrap()
            .push(serde_json::json!({
                "name": "back_mutation",
                "reactants": [{"label": "x", "type": "M"}],
                "propensity": "1e-6 * M",
                "products": [{"type": "W"}],
            }));
        assert_eq!(refused(&model), Some(TOPOLOGY));
        // A capacity below the initial population.
        let mut model = base();
        model["parameters"]["K"] = 999.0.into();
        assert_eq!(refused(&model), Some(CAPACITY));
        // Stops that are not the guarded mutant fraction over counts.
        for stop in [
            "W + M + A > 0 && M + A > 0.9 * (W + M + A)",
            "W + M + A > 0 && M >= 0.9 * (W + M + A)",
            "W + M > 0 && M + A >= 0.9 * (W + M)",
            "W + M + A > 0 && M + A >= 1.5 * (W + M + A)",
            "M + A >= 0.9 * (W + M + A)",
            "W + M + A > 0 && (M + A) / (W + M + A) >= 0.9",
        ] {
            let model = builder_model(["W", "M", "A"], FREE, Some(stop));
            assert_eq!(refused(&model), Some(STOP), "{stop}");
        }
    }
}
