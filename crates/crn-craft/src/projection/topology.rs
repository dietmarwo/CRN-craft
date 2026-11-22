//! Enumeration of the reachable states and edges of a prepared projection.
use super::{Interrupt, Operation, budget, input};
use crate::diagnostic::Diagnostic;
use crate::rules::FlatView;
use std::collections::{HashMap, HashSet};

/// Most states of a prepared projection, whatever `max_states` says.
pub const MAX_STATES: usize = 250_000;
/// Most generator edges of a prepared projection, whatever `max_edges` says.
pub const MAX_EDGES: usize = 2_000_000;
/// The state-table rule: (states + 1) * (128 + 40 * species) bytes at most.
pub const STATE_TABLE_BYTES: u64 = 256 << 20;
/// The largest exact count, 2^53.
pub const MAX_COUNT: u64 = 1 << 53;
/// The largest `min(k, n - k)` of a reactant binomial coefficient.
const MAX_CHOOSE: u64 = 1024;

/// One channel of the lowered network as the enumeration reads it.
#[derive(Debug, Clone)]
pub struct ChannelShape {
    /// Reactants: (species, molecules taken), sorted by species.
    pub needs: Vec<(usize, u64)>,
    /// Net change of one firing: (species, change), sorted by species.
    pub changes: Vec<(usize, i64)>,
    /// Whether the channel's rule has a whole-rule propensity (multiplicity 1).
    pub whole: bool,
}

impl ChannelShape {
    /// The channels of a flat view, in channel order.
    #[must_use]
    pub fn of_view(view: &FlatView) -> Vec<Self> {
        (0..view.channels())
            .map(|c| Self {
                needs: view.reactants(c).to_vec(),
                changes: view.changes(c).to_vec(),
                whole: view.whole(c),
            })
            .collect()
    }
}

/// What a projection enumerates: bounds, the initial support, channels and
/// the count operations of scheduled interventions.
#[derive(Debug, Clone)]
pub struct Spec {
    /// Species labels in state order (used in messages).
    pub species: Vec<String>,
    /// Inclusive upper bound per species, at most 2^53.
    pub bounds: Vec<u64>,
    /// The states of the initial law, in the order they are admitted (the
    /// caller sorts them). Counts are signed so that invalid ones are
    /// refused here.
    pub initial: Vec<Vec<i64>>,
    /// The lowered channels.
    pub channels: Vec<ChannelShape>,
    /// The count operations of each scheduled action that has any, in
    /// schedule order; within an action, in the order they apply.
    pub interventions: Vec<Vec<(usize, Operation)>>,
    /// Most states, `1..=MAX_STATES`.
    pub max_states: usize,
    /// Most edges, `1..=MAX_EDGES`.
    pub max_edges: usize,
}

/// The reachable states inside the bounds and the generator edges between
/// them. The state, edge and multiplicity order is the discovery order of a
/// breadth-first search from the sorted initial support, with each
/// intervention's outcomes admitted after the states reached before it.
#[derive(Debug, Clone)]
pub struct Topology {
    species: Vec<String>,
    bounds: Vec<u64>,
    width: usize,
    count: usize,
    states: Vec<u64>,
    index: HashMap<Box<[u64]>, usize>,
    src: Vec<usize>,
    dst: Vec<usize>,
    channel: Vec<usize>,
    multiplicity: Vec<f64>,
    max_states: usize,
}

/// A state is inside the bounds (an index) or beyond them (the sink).
enum Admitted {
    State(usize),
    Sink,
}

struct Builder<'a> {
    spec: &'a Spec,
    retained: Option<&'a HashSet<Vec<u64>>>,
    width: usize,
    count: usize,
    states: Vec<u64>,
    index: HashMap<Box<[u64]>, usize>,
}

impl Builder<'_> {
    fn len(&self) -> usize {
        self.count
    }

    fn admit(&mut self, state: &[i128]) -> Result<Admitted, Diagnostic> {
        if state.len() != self.width || state.iter().any(|&n| n < 0) {
            return Err(input("invalid initial or reachable counts"));
        }
        if state
            .iter()
            .zip(&self.spec.bounds)
            .any(|(&n, &bound)| n > i128::from(bound))
        {
            return Ok(Admitted::Sink);
        }
        // Inside the bounds every count fits in u64.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let key: Box<[u64]> = state.iter().map(|&n| n as u64).collect();
        if self.retained.is_some_and(|states| !states.contains(&*key)) {
            return Ok(Admitted::Sink);
        }
        if let Some(&i) = self.index.get(&key) {
            return Ok(Admitted::State(i));
        }
        let len = self.len();
        if len >= self.spec.max_states {
            return Err(budget("reachable-state budget exceeded"));
        }
        // Reserve the state table and the indices used for lookups.
        if (len as u64 + 1) * (128 + 40 * self.width as u64) > STATE_TABLE_BYTES {
            return Err(budget(
                "projection state table exceeds the workspace budget",
            ));
        }
        self.states.extend_from_slice(&key);
        self.index.insert(key, len);
        self.count += 1;
        Ok(Admitted::State(len))
    }

    fn state(&self, i: usize) -> &[u64] {
        &self.states[i * self.width..(i + 1) * self.width]
    }
}

/// The outcomes of a count operation on `count` molecules, as the range
/// `first..=last`; `first` is negative when a decrement exceeds the count.
fn outcomes(count: u64, operation: Operation) -> (i128, i128) {
    let count = i128::from(count);
    match operation {
        Operation::Partition(p) if p > 0.0 && p < 1.0 => (0, count),
        Operation::Partition(p) => {
            let n = if p == 0.0 { 0 } else { count };
            (n, n)
        }
        Operation::Set(n) => (i128::from(n), i128::from(n)),
        Operation::Change(d) => (count + i128::from(d), count + i128::from(d)),
    }
}

impl Topology {
    /// Enumerates the projection of `spec`.
    ///
    /// # Errors
    ///
    /// `CRN-E071` for an initial state outside the bounds, an invalid count
    /// or reactant combinatorics beyond `min(k, n - k) = 1024`, at
    /// `/projection/overflow` for a multiplicity above the largest double;
    /// `CRN-E070` for the state, edge, state-table and enumeration budgets
    /// and for an interruption.
    pub fn enumerate(spec: &Spec, interrupt: &Interrupt<'_>) -> Result<Self, Diagnostic> {
        Self::enumerate_retained(spec, None, interrupt)
    }

    /// Enumerate within a finite retained set. Omitted transitions flow
    /// straight to the existing absorbing sink.
    ///
    /// # Errors
    /// Invalid retained states, or any error of ordinary enumeration.
    #[allow(clippy::too_many_lines)]
    pub fn enumerate_retained(
        spec: &Spec,
        retained: Option<&HashSet<Vec<u64>>>,
        interrupt: &Interrupt<'_>,
    ) -> Result<Self, Diagnostic> {
        let width = spec.bounds.len();
        if retained.is_some_and(|states| {
            states.is_empty()
                || states.len() > spec.max_states
                || states.iter().any(|s| {
                    s.len() != width || s.iter().zip(&spec.bounds).any(|(n, bound)| n > bound)
                })
        }) {
            return Err(input("invalid retained projection states"));
        }
        let mut builder = Builder {
            spec,
            retained,
            width,
            count: 0,
            states: Vec::new(),
            index: HashMap::new(),
        };
        let max_edges = spec.max_edges as i128;
        let mut src = Vec::new();
        let mut dst = Vec::new();
        let mut channel = Vec::new();
        let mut multiplicity = Vec::new();
        let mut touched: i128 = 0;
        let mut scratch: Vec<i128> = Vec::with_capacity(width);
        for state in &spec.initial {
            scratch.clear();
            scratch.extend(state.iter().map(|&n| i128::from(n)));
            if let Admitted::Sink = builder.admit(&scratch)? {
                return Err(input("initial law lies outside the projection"));
            }
        }
        let mut cursor = 0;
        for changes in std::iter::once(None).chain(spec.interventions.iter().map(Some)) {
            if let Some(changes) = changes.filter(|c| !c.is_empty()) {
                // Every state reached before this action may be its input.
                let reached = builder.len();
                for s in 0..reached {
                    let state = builder.state(s).to_vec();
                    let ranges: Vec<(i128, i128)> = changes
                        .iter()
                        .map(|&(i, op)| outcomes(state[i], op))
                        .collect();
                    if ranges.iter().any(|&(first, _)| first < 0) {
                        continue;
                    }
                    let total = ranges.iter().fold(1u128, |acc, &(first, last)| {
                        acc.saturating_mul((last - first + 1).unsigned_abs())
                    });
                    let room = max_edges * 4 - touched;
                    if room < 0 || total > room.unsigned_abs() {
                        return Err(budget("intervention expansion exceeds work budget"));
                    }
                    // Odometer over the outcome ranges, the last fastest.
                    let mut values: Vec<i128> = ranges.iter().map(|r| r.0).collect();
                    'product: loop {
                        touched += 1;
                        if touched > max_edges * 4 {
                            return Err(budget("intervention expansion exceeds work budget"));
                        }
                        if touched % 4096 == 0 {
                            interrupt.check()?;
                        }
                        scratch.clear();
                        scratch.extend(state.iter().map(|&n| i128::from(n)));
                        for (&(i, _), &v) in changes.iter().zip(&values) {
                            scratch[i] = v;
                        }
                        builder.admit(&scratch)?;
                        let mut d = values.len();
                        loop {
                            if d == 0 {
                                break 'product;
                            }
                            d -= 1;
                            if values[d] < ranges[d].1 {
                                values[d] += 1;
                                break;
                            }
                            values[d] = ranges[d].0;
                        }
                    }
                }
            }
            while cursor < builder.len() {
                let state = builder.state(cursor).to_vec();
                for (c, shape) in spec.channels.iter().enumerate() {
                    touched += 1;
                    if touched > max_edges * 25 {
                        return Err(budget("projection enumeration exceeds work budget"));
                    }
                    if touched % 4096 == 0 {
                        interrupt.check()?;
                    }
                    if shape.needs.iter().any(|&(j, n)| state[j] < n) {
                        continue;
                    }
                    scratch.clear();
                    scratch.extend(state.iter().map(|&n| i128::from(n)));
                    let mut escaped = false;
                    for &(j, delta) in &shape.changes {
                        scratch[j] += i128::from(delta);
                        if scratch[j] > i128::from(spec.bounds[j]) {
                            escaped = true;
                            break;
                        }
                    }
                    if !escaped
                        && scratch
                            .iter()
                            .zip(&state)
                            .all(|(&a, &b)| a == i128::from(b))
                    {
                        continue;
                    }
                    if src.len() >= spec.max_edges {
                        return Err(budget("projection edge budget exceeded"));
                    }
                    let target = if escaped {
                        usize::MAX
                    } else {
                        match builder.admit(&scratch)? {
                            Admitted::State(j) => j,
                            Admitted::Sink => usize::MAX,
                        }
                    };
                    src.push(cursor);
                    dst.push(target);
                    channel.push(c);
                    multiplicity.push(if shape.whole {
                        1.0
                    } else {
                        combinations(&state, &shape.needs)?
                    });
                }
                cursor += 1;
            }
        }
        let n = builder.len();
        for d in &mut dst {
            if *d == usize::MAX {
                *d = n;
            }
        }
        Ok(Self {
            species: spec.species.clone(),
            bounds: spec.bounds.clone(),
            width,
            count: n,
            states: builder.states,
            index: builder.index,
            src,
            dst,
            channel,
            multiplicity,
            max_states: spec.max_states,
        })
    }

    /// The number of projected states (the sink, index `len()`, excluded).
    #[must_use]
    pub fn len(&self) -> usize {
        self.count
    }

    /// Whether no state is projected (never after a successful enumeration).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// The number of species of a state.
    #[must_use]
    pub fn width(&self) -> usize {
        self.width
    }

    /// Species labels in state order.
    #[must_use]
    pub fn species(&self) -> &[String] {
        &self.species
    }

    /// Inclusive bound per species.
    #[must_use]
    pub fn bounds(&self) -> &[u64] {
        &self.bounds
    }

    /// The configured state budget (the intervention workspace limit).
    #[must_use]
    pub fn max_states(&self) -> usize {
        self.max_states
    }

    /// All states, row-major (`len() * width()` counts).
    #[must_use]
    pub fn states(&self) -> &[u64] {
        &self.states
    }

    /// The counts of state `i`.
    #[must_use]
    pub fn state(&self, i: usize) -> &[u64] {
        &self.states[i * self.width..(i + 1) * self.width]
    }

    /// The index of a state, `None` outside the projection.
    #[must_use]
    pub fn index_of(&self, counts: &[u64]) -> Option<usize> {
        self.index.get(counts).copied()
    }

    /// Edge sources, in edge order (nondecreasing).
    #[must_use]
    pub fn sources(&self) -> &[usize] {
        &self.src
    }

    /// Edge targets; `len()` is the sink.
    #[must_use]
    pub fn targets(&self) -> &[usize] {
        &self.dst
    }

    /// The lowered channel of each edge.
    #[must_use]
    pub fn channels(&self) -> &[usize] {
        &self.channel
    }

    /// The reactant multiplicity of each edge: `prod C(n_i, k_i)` rounded
    /// once to the nearest double, or 1 for a whole-rule propensity.
    #[must_use]
    pub fn multiplicities(&self) -> &[f64] {
        &self.multiplicity
    }
}

/// `prod_i C(state[i], k_i)` over the reactants, exact, rounded once to the
/// nearest double (ties to even), as Python's `float(int)` rounds.
fn combinations(state: &[u64], needs: &[(usize, u64)]) -> Result<f64, Diagnostic> {
    let mut product = Big::one();
    let mut overflow = false;
    for &(i, k) in needs {
        let n = state[i];
        let choose = k.min(n - k);
        if choose > MAX_CHOOSE {
            return Err(input(
                "reactant combinatorics exceed probability solver range",
            ));
        }
        if overflow {
            continue;
        }
        // C(n, c) = prod_{j=1..c} (n - c + j) / j; each partial product times
        // the earlier factors is an integer.
        for j in 1..=choose {
            product.mul_small(n - choose + j);
            product.div_exact_small(j);
        }
        // Beyond 2^1100 no later factor (each >= 1) brings it back below
        // the largest double.
        if product.bits() > 1100 {
            overflow = true;
        }
    }
    match product.to_f64() {
        Some(value) if !overflow => Ok(value),
        _ => Err(Diagnostic::error(
            "CRN-E071",
            "/projection/overflow",
            "int too large to convert to float",
            "Lower the bounds of species that react with themselves.",
        )),
    }
}

/// A nonnegative integer in little-endian 64-bit limbs.
struct Big(Vec<u64>);

impl Big {
    fn one() -> Self {
        Self(vec![1])
    }

    #[allow(clippy::cast_possible_truncation)]
    fn mul_small(&mut self, factor: u64) {
        let mut carry: u128 = 0;
        for limb in &mut self.0 {
            let value = u128::from(*limb) * u128::from(factor) + carry;
            *limb = value as u64;
            carry = value >> 64;
        }
        if carry > 0 {
            self.0.push(carry as u64);
        }
        self.trim();
    }

    #[allow(clippy::cast_possible_truncation)]
    fn div_exact_small(&mut self, divisor: u64) {
        let mut remainder: u128 = 0;
        for limb in self.0.iter_mut().rev() {
            let value = (remainder << 64) | u128::from(*limb);
            *limb = (value / u128::from(divisor)) as u64;
            remainder = value % u128::from(divisor);
        }
        debug_assert_eq!(remainder, 0);
        self.trim();
    }

    fn trim(&mut self) {
        while self.0.len() > 1 && self.0.last() == Some(&0) {
            self.0.pop();
        }
    }

    fn bits(&self) -> u64 {
        let top = *self.0.last().unwrap_or(&0);
        (self.0.len() as u64 - 1) * 64 + u64::from(64 - top.leading_zeros())
    }

    /// Correctly rounded conversion; `None` when it overflows.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        clippy::cast_precision_loss
    )]
    fn to_f64(&self) -> Option<f64> {
        let bits = self.bits();
        if bits <= 64 {
            return Some(self.0[0] as f64);
        }
        // The top 64 bits with every lower bit folded into the lowest one:
        // u64 -> f64 rounds to nearest even, and the folded bit lies below
        // the rounding position, so it acts as the sticky bit.
        let shift = bits - 64;
        let (word, offset) = ((shift / 64) as usize, (shift % 64) as u32);
        let mut top = self.0[word] >> offset;
        if offset > 0 {
            top |= self.0[word + 1] << (64 - offset);
        }
        let lower_nonzero = self.0[..word].iter().any(|&l| l != 0)
            || (offset > 0 && self.0[word] & ((1u64 << offset) - 1) != 0);
        if lower_nonzero {
            top |= 1;
        }
        let value = libm::scalbn(top as f64, shift as i32);
        value.is_finite().then_some(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combinations_round_like_python() {
        // C(10, 2) * C(7, 3) = 45 * 35
        assert_eq!(combinations(&[10, 7], &[(0, 2), (1, 3)]).unwrap(), 1575.0);
        // C(2^53, 2) = 2^52 (2^53 - 1): 105 bits, rounds to 2^105 - 2^52 exactly
        // representable? 2^52 * (2^53 - 1) has 53 significant bits.
        let n = 1u64 << 53;
        let value = combinations(&[n], &[(0, 2)]).unwrap();
        assert_eq!(value, 4_503_599_627_370_496.0 * 9_007_199_254_740_991.0);
        // C(2^53, 3) = 2^53 (2^53 - 1)(2^53 - 2) / 6, compared with an
        // independent u128-free computation: the f64 product of the exact
        // factors, rounded at most three times, within 2 ulp.
        let value = combinations(&[n], &[(0, 3)]).unwrap();
        let approx = (n as f64) * ((n - 1) as f64) * ((n - 2) as f64) / 6.0;
        assert!((value - approx).abs() <= 2.0 * approx * f64::EPSILON);
        // Too large for a double.
        let error = combinations(&[n], &[(0, 1000)]).unwrap_err();
        assert_eq!(error.path, "/projection/overflow");
        let error = combinations(&[5000], &[(0, 2000)]).unwrap_err();
        assert_eq!(error.code, "CRN-E071");
    }

    #[test]
    fn big_rounding_ties_to_even() {
        // 2^64 + 2^11: exactly halfway between 2^64 and 2^64 + 2^12 -> even (2^64).
        let big = Big(vec![1 << 11, 1]);
        assert_eq!(big.to_f64().unwrap(), 18_446_744_073_709_551_616.0);
        // One more unit is above the half: rounds up.
        let big = Big(vec![(1 << 11) + 1, 1]);
        assert_eq!(big.to_f64().unwrap(), 18_446_744_073_709_551_616.0 + 4096.0);
        // A sticky bit below the 64-bit window: 2^128 + 2^75 + 1.
        let big = Big(vec![1, 1 << 11, 1]);
        let exact_half = libm::ldexp(1.0, 128);
        assert!(big.to_f64().unwrap() > exact_half);
    }
}
