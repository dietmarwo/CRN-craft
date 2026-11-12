//! Event selection of Static Direct over the lowered channels
//! (`crn-craft.run/1.0`, `method.selector`).
//!
//! Ported from the 0.6.0 network engine (`ssa/selectors.rs`), with the
//! random words of the rules stream ([`Words`]):
//!
//! - [`RateTree`]: pairwise sums of the channel rates, so the total depends
//!   only on the current rates, with the strict checks of the linear scan
//!   (a non-finite total and a positive rate below 2⁻⁵³ of the total are
//!   refused) evaluated lazily;
//! - [`Composition`]: composition–rejection over groups by binary exponent
//!   (Slepoy, Thompson & Plimpton, J. Chem. Phys. 128, 205101 (2008)): a
//!   group is chosen by its envelope mass, then a member uniformly and
//!   accepted with probability rate / group ceiling;
//! - [`sorting_sample`]: sorting direct (McCollum et al., Computational
//!   Biology and Chemistry 30 (2006), 39–49): a linear scan in an order that
//!   moves the fired channel one place forward.
//!
//! RSSA-CR's propensity bounds are in `rssa.rs`.
use crate::rules::engine::{Halt, RESOLUTION, Words, numerical};
use crate::rules::flat::Network;
use crate::rules::{RunVersion, Selector};
use rand_core::RngCore;

/// The selector that runs, resolved from the request at preparation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    Linear,
    CompositionRejection,
    Sorting,
    RssaCr,
}

impl Kind {
    /// The selector's part of the sampler identity; empty for the linear
    /// scan, whose identities are those of 0.6.0.
    pub(crate) fn sampler_suffix(self) -> &'static str {
        match self {
            Self::Linear => "",
            Self::CompositionRejection => "-cr",
            Self::Sorting => "-sorting",
            Self::RssaCr => "-rssa-cr",
        }
    }
}

/// `auto` picks composition–rejection or RSSA-CR from this many channels on
/// (declared benchmark in `validation/merge/selector-threshold.json`).
pub(crate) const AUTO_CHANNELS: usize = 64;
/// `auto` picks RSSA-CR over composition–rejection when the channels'
/// dependency lists hold at least this many channels on average.
pub(crate) const AUTO_RSSA_AFFECTED: f64 = 16.0;
/// `auto` picks RSSA-CR from this many channels on when at least half of
/// the channels evaluate an expression or a combination count.
pub(crate) const AUTO_RSSA_EXPRESSIONS: usize = 64;

/// Structural features of a lowered network that `auto` reads.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Features {
    pub channels: usize,
    /// The mean length of the channels' dependency lists (`affected`).
    pub mean_affected: f64,
    /// The fraction of channels whose rate is an expression of observables
    /// or a scaled combination count, rather than a constant.
    pub expressions: f64,
    /// Every rule's rate reads at most one species, so every RSSA-CR bound
    /// is an interval bound rather than a point.
    pub interval_bounds: bool,
}

impl Features {
    pub(crate) fn of(net: &Network) -> Self {
        let channels = net.channels.len();
        let affected: usize = net.affected.iter().map(Vec::len).sum();
        let expressions = (0..channels)
            .filter(|&i| net.constant_global(i).is_none() || net.scaled_combination(i))
            .count();
        Self {
            channels,
            mean_affected: affected as f64 / channels.max(1) as f64,
            expressions: expressions as f64 / channels.max(1) as f64,
            interval_bounds: net.rule_reads.iter().all(|r| r.len() <= 1),
        }
    }

    /// The `auto` rule.
    pub(crate) fn auto(self) -> Kind {
        let rssa = self.interval_bounds
            && ((self.expressions >= 0.5 && self.channels >= AUTO_RSSA_EXPRESSIONS)
                || (self.channels >= AUTO_CHANNELS && self.mean_affected >= AUTO_RSSA_AFFECTED));
        if rssa {
            Kind::RssaCr
        } else if self.channels >= AUTO_CHANNELS {
            Kind::CompositionRejection
        } else {
            Kind::Linear
        }
    }
}

/// The selector a Static Direct request runs: always the linear scan before
/// `crn-craft.run/1.0`; there the requested one, `auto` when absent.
pub(crate) fn resolve(
    version: RunVersion,
    requested: Option<Selector>,
    net: Option<&Network>,
) -> Kind {
    if version < RunVersion::V10 {
        return Kind::Linear;
    }
    match requested.unwrap_or(Selector::Auto) {
        Selector::Linear => Kind::Linear,
        Selector::CompositionRejection => Kind::CompositionRejection,
        Selector::Sorting => Kind::Sorting,
        Selector::RssaCr => Kind::RssaCr,
        Selector::Auto => net.map_or(Kind::Linear, |net| Features::of(net).auto()),
    }
}

fn unresolvable() -> Halt {
    numerical(
        "/rules",
        "unresolvable total propensity",
        "Rescale rates or populations.",
    )
}

fn below_resolution() -> Halt {
    Halt::Numerical(crate::diagnostic::Diagnostic::error(
        "CRN-R031",
        "/rules",
        "a positive channel is below 53-bit resolution of the total",
        "Use the tree method or rescale rates.",
    ))
}

#[derive(Clone, Copy, Debug)]
struct Node {
    sum: f64,
    minimum: f64,
    bad: bool,
}

impl Default for Node {
    fn default() -> Self {
        Self {
            sum: 0.0,
            minimum: f64::INFINITY,
            bad: false,
        }
    }
}

fn leaf(rate: f64, minimum_positive: f64) -> Node {
    Node {
        sum: rate,
        minimum: minimum_positive,
        bad: !rate.is_finite()
            || rate < 0.0
            || minimum_positive <= 0.0
            || minimum_positive.is_nan(),
    }
}

fn parent(a: Node, b: Node) -> Node {
    let sum = a.sum + b.sum;
    Node {
        sum,
        minimum: a.minimum.min(b.minimum),
        bad: a.bad
            || b.bad
            || !sum.is_finite()
            || (a.sum > 0.0 && sum <= b.sum)
            || (b.sum > 0.0 && sum <= a.sum),
    }
}

/// Pairwise sums of the channel rates, rebuilt along the ancestor chains of
/// the changed leaves; no accumulated add/subtract drift.
///
/// Only the sums are maintained per node. [`RateTree::total`] decides
/// exactly as a tree that also carries per-node minima and loss flags: it
/// rejects an invalid leaf, a non-finite sum, a pairwise addition that loses
/// a positive operand (`CRN-R030`), and a smallest positive leaf minimum below
/// 2⁻⁵³ of the total (`CRN-R031`). `floor` never exceeds a current positive
/// leaf or leaf minimum. When it is at least 2⁻⁵² of a normal total, every
/// positive pairwise operand exceeds half an ulp of any sum (all sums are at
/// most the total), so no addition can lose it, and the resolution test
/// passes; otherwise the checked tree is evaluated in full.
///
/// Leaves are staged with [`RateTree::set`]; [`RateTree::commit`] (also
/// called by [`RateTree::total`]) recomputes the ancestors of the staged
/// leaves one level at a time, each ancestor once, which gives the sums of a
/// per-leaf update in fewer additions when the changed channels are close
/// together.
#[derive(Default)]
pub(crate) struct RateTree {
    sums: Vec<f64>,
    minima: Vec<f64>,
    base: usize,
    invalid: usize,
    floor: f64,
    checked: Vec<Node>,
    staged: Vec<usize>,
    level: Vec<usize>,
}

impl RateTree {
    pub(crate) fn reset(&mut self, count: usize) {
        self.base = count.max(1).next_power_of_two();
        self.sums.clear();
        self.sums.resize(self.base * 2, 0.0);
        self.minima.clear();
        self.minima.resize(self.base, f64::INFINITY);
        self.invalid = 0;
        self.floor = f64::INFINITY;
        self.staged.clear();
    }

    #[inline]
    pub(crate) fn get(&self, index: usize) -> f64 {
        self.sums[self.base + index]
    }

    /// Stages leaf `index` with `rate` and the smallest positive value the
    /// channel can take (`rate` itself for an exact rate, infinite for none).
    #[inline]
    pub(crate) fn set(&mut self, index: usize, rate: f64, minimum_positive: f64) {
        let i = self.base + index;
        self.invalid -= usize::from(leaf(self.sums[i], self.minima[index]).bad);
        self.invalid += usize::from(leaf(rate, minimum_positive).bad);
        self.sums[i] = rate;
        self.minima[index] = minimum_positive;
        let positive = if rate > 0.0 { rate } else { f64::INFINITY };
        self.floor = self.floor.min(positive.min(minimum_positive));
        self.staged.push(index);
    }

    /// Recomputes the ancestors of the staged leaves.
    pub(crate) fn commit(&mut self) {
        if self.staged.is_empty() {
            return;
        }
        if self.base == 1 {
            // The only leaf is the root.
            self.staged.clear();
            return;
        }
        let mut level = std::mem::take(&mut self.level);
        level.clear();
        self.staged.sort_unstable();
        for &index in &self.staged {
            let node = (self.base + index) / 2;
            if level.last() != Some(&node) {
                level.push(node);
            }
        }
        self.staged.clear();
        loop {
            for &node in &level {
                self.sums[node] = self.sums[2 * node] + self.sums[2 * node + 1];
            }
            if level[0] == 1 {
                break;
            }
            let mut kept = 0;
            for k in 0..level.len() {
                let node = level[k] / 2;
                if kept == 0 || level[kept - 1] != node {
                    level[kept] = node;
                    kept += 1;
                }
            }
            level.truncate(kept);
        }
        self.level = level;
    }

    /// The total of the committed leaves, with the checks of the linear scan.
    pub(crate) fn total(&mut self) -> Result<f64, Halt> {
        self.commit();
        let root = self.sums[1];
        if self.invalid == 0 && root.is_normal() && self.floor >= root * (1.0 / (1u64 << 52) as f64)
        {
            return Ok(root);
        }
        self.checked_total()
    }

    #[cold]
    fn checked_total(&mut self) -> Result<f64, Halt> {
        let base = self.base;
        self.checked.resize(2 * base, Node::default());
        for i in 0..base {
            self.checked[base + i] = leaf(self.sums[base + i], self.minima[i]);
        }
        for i in (1..base).rev() {
            self.checked[i] = parent(self.checked[2 * i], self.checked[2 * i + 1]);
        }
        let root = self.checked[1];
        debug_assert_eq!(root.sum.to_bits(), self.sums[1].to_bits());
        let positive = (0..base)
            .map(|i| self.sums[base + i])
            .filter(|&rate| rate > 0.0)
            .fold(root.minimum, f64::min);
        self.floor = positive;
        // As the linear scan: a non-finite total or rate is `CRN-R030`, a
        // finite total that loses or cannot resolve a positive rate
        // `CRN-R031`.
        if self.invalid > 0 || !root.sum.is_finite() {
            return Err(unresolvable());
        }
        if root.bad || (root.minimum.is_finite() && root.minimum / root.sum < RESOLUTION) {
            return Err(below_resolution());
        }
        Ok(root.sum)
    }
}

const GROUPS: usize = 2098;

/// The group of a positive finite rate: its binary exponent, with
/// subnormals by their leading bit.
fn group(rate: f64) -> usize {
    let bits = rate.to_bits();
    let exponent = (bits >> 52) as usize;
    if exponent == 0 {
        (63 - bits.leading_zeros()) as usize
    } else {
        exponent + 51
    }
}

/// The power of two above every rate of a group.
fn ceiling(group: usize) -> f64 {
    let exponent = group as i32 + 1 - 1074;
    if exponent > 1023 {
        f64::MAX
    } else if exponent >= -1022 {
        f64::from_bits(((exponent + 1023) as u64) << 52)
    } else {
        f64::from_bits(1u64 << (exponent + 1074))
    }
}

/// Groups up to this size take the member and the acceptance uniform from
/// one random word. Lemire's product `w * n` gives the member in its high 64
/// bits; after his rejection step, the low 64 bits are, given the member,
/// uniform on a lattice of spacing n / 2⁶⁴ ≤ 2⁻⁵⁴. Their top 53 bits are the
/// acceptance uniform, so the acceptance probability errs by at most
/// 2⁻⁵³ + 2⁻⁵⁴.
const JOINT_MEMBERS: usize = 1024;

fn joint<R: RngCore>(
    members: &[usize],
    bound: f64,
    words: &mut Words<'_, R>,
) -> Result<(usize, f64), Halt> {
    let n = members.len() as u64;
    let cutoff = n.wrapping_neg() % n;
    loop {
        let product = u128::from(words.word()?) * u128::from(n);
        let low = product as u64;
        if low >= cutoff {
            let uniform = crate::sampling::uniform_closed_open(low);
            return Ok((members[(product >> 64) as usize], uniform * bound));
        }
    }
}

/// A uniform index below `n` (Lemire's method).
fn index<R: RngCore>(n: usize, words: &mut Words<'_, R>) -> Result<usize, Halt> {
    let n = n as u64;
    let cutoff = n.wrapping_neg() % n;
    loop {
        let product = u128::from(words.word()?) * u128::from(n);
        if product as u64 >= cutoff {
            return Ok((product >> 64) as usize);
        }
    }
}

/// Composition–rejection: channels grouped by the binary exponent of their
/// rate (or RSSA-CR upper bound).
#[derive(Default)]
pub(crate) struct Composition {
    groups: Vec<Vec<usize>>,
    /// The non-empty groups, increasing.
    active: Vec<usize>,
    /// (group, position) of each channel; `usize::MAX` for a zero rate.
    locations: Vec<(usize, usize)>,
    /// Running envelope sums over `active`, rebuilt only after a membership
    /// change, with the same additions in the same order as a fresh scan.
    cumulative: Vec<f64>,
    stale: bool,
}

impl Composition {
    pub(crate) fn reset(&mut self, count: usize) {
        if self.groups.is_empty() {
            self.groups.resize_with(GROUPS, Vec::new);
        }
        for id in self.active.drain(..) {
            self.groups[id].clear();
        }
        self.locations.clear();
        self.locations.resize(count, (usize::MAX, 0));
        self.stale = true;
    }

    /// Moves channel `index` to the group of `value`, a finite nonnegative
    /// rate.
    #[inline]
    pub(crate) fn update(&mut self, index: usize, value: f64) -> Result<(), Halt> {
        if !value.is_finite() || value < 0.0 {
            return Err(unresolvable());
        }
        let next = if value == 0.0 {
            usize::MAX
        } else {
            group(value)
        };
        let (old, pos) = self.locations[index];
        if next == old {
            return Ok(());
        }
        self.stale = true;
        if old != usize::MAX {
            self.groups[old].swap_remove(pos);
            if let Some(&moved) = self.groups[old].get(pos) {
                self.locations[moved].1 = pos;
            }
            // Bound retained capacity by live memberships, rather than by
            // every exponent bucket a channel has visited along a path.
            if self.groups[old].capacity() > 4 * self.groups[old].len().max(4) {
                self.groups[old].shrink_to_fit();
            }
            if self.groups[old].is_empty() {
                let at = self.active.binary_search(&old).expect("active group");
                self.active.remove(at);
            }
        }
        if next == usize::MAX {
            self.locations[index] = (next, 0);
        } else {
            if self.groups[next].is_empty() {
                let at = self
                    .active
                    .binary_search(&next)
                    .expect_err("inactive group");
                self.active.insert(at, next);
            }
            self.locations[index] = (next, self.groups[next].len());
            self.groups[next].push(index);
        }
        Ok(())
    }

    fn rebuild(&mut self) -> Result<(), Halt> {
        let largest = *self.active.last().ok_or_else(unresolvable)?;
        let scale = ceiling(largest);
        self.cumulative.clear();
        let mut sum = 0.0;
        for &id in &self.active {
            let weight = (ceiling(id) / scale) * self.groups[id].len() as f64;
            let next = sum + weight;
            if weight <= 0.0 || next <= sum {
                return Err(below_resolution());
            }
            sum = next;
            self.cumulative.push(sum);
        }
        self.stale = false;
        Ok(())
    }

    /// A channel with probability `tree.get(i)` / total over the members.
    pub(crate) fn sample<R: RngCore>(
        &mut self,
        tree: &RateTree,
        words: &mut Words<'_, R>,
    ) -> Result<usize, Halt> {
        if self.stale {
            self.rebuild()?;
        }
        let sum = *self.cumulative.last().ok_or_else(unresolvable)?;
        loop {
            let threshold = words.uniform()? * sum;
            if threshold >= sum {
                return Err(numerical(
                    "/rules",
                    "composition threshold rounded to its total",
                    "Rescale rates or populations.",
                ));
            }
            // The first running sum above the threshold, as in a scan.
            let k = self.cumulative.partition_point(|&c| c <= threshold);
            let id = *self.active.get(k).ok_or_else(unresolvable)?;
            let members = &self.groups[id];
            let bound = ceiling(id);
            let (channel, threshold) = if members.len() <= JOINT_MEMBERS {
                joint(members, bound, words)?
            } else {
                let channel = members[index(members.len(), words)?];
                (channel, words.uniform()? * bound)
            };
            if threshold >= bound {
                return Err(numerical(
                    "/rules",
                    "composition acceptance rounded to its bound",
                    "Rescale rates or populations.",
                ));
            }
            if threshold < tree.get(channel) {
                return Ok(channel);
            }
        }
    }
}

/// Sorting direct: a scan of `rates` in `order`; the fired channel moves one
/// position forward. The threshold scales the tree's pairwise total, but the
/// scan adds rates in the learned order, and the two sums may differ in the
/// last bit. A threshold at or above the sequential sum is redrawn, so
/// accepted thresholds are uniform below the sum the scan partitions.
pub(crate) fn sorting_sample<R: RngCore>(
    order: &mut [usize],
    rates: &[f64],
    total: f64,
    words: &mut Words<'_, R>,
) -> Result<usize, Halt> {
    loop {
        let threshold = words.uniform()? * total;
        if threshold >= total {
            return Err(numerical(
                "/rules",
                "sorting threshold rounded to its total",
                "Rescale rates or populations.",
            ));
        }
        let mut sum = 0.0;
        for (position, &channel) in order.iter().enumerate() {
            let rate = rates[channel];
            let next = sum + rate;
            if rate > 0.0 && next <= sum {
                return Err(below_resolution());
            }
            sum = next;
            if threshold < sum {
                if position > 0 {
                    order.swap(position, position - 1);
                }
                return Ok(channel);
            }
        }
        if sum == 0.0 {
            return Err(unresolvable());
        }
    }
}

#[cfg(test)]
mod tests {
    //! Ported from the network engine's `ssa/selectors.rs`, on the rules
    //! engine's random words.
    use super::*;
    use crate::rules::engine::Words;
    use rand_core::RngCore;

    #[test]
    fn composition_matches_declared_categorical_dkw_bound_after_membership_changes() {
        let rates = [0.0, 0.125, 0.5, 1.0, 2.0, 7.0, 8.0, 50.0];
        let mut tree = RateTree::default();
        let mut selector = Composition::default();
        tree.reset(rates.len());
        selector.reset(rates.len());
        // Exercise removal, group migration, zero activation and swap-removal.
        for (i, &rate) in rates.iter().enumerate() {
            selector.update(i, 8.0).ok().unwrap();
            selector.update(i, rate).ok().unwrap();
            tree.set(i, rate, if rate > 0.0 { rate } else { f64::INFINITY });
        }
        let total = tree.total().ok().unwrap();
        let mut rng = crate::sampling::stream(4242, 777);
        let mut words = Words::for_test(&mut rng, 20_000_000);
        let mut counts = [0u64; 8];
        const N: u64 = 1_000_000;
        for _ in 0..N {
            let index = selector.sample(&tree, &mut words).ok().unwrap();
            counts[index] += 1;
        }
        assert_eq!(counts[0], 0);
        let bound = (libm::log(2.0 * 100.0 / 0.001) / (2.0 * N as f64)).sqrt();
        let (mut observed, mut expected) = (0.0, 0.0);
        for (&count, &rate) in counts.iter().zip(&rates) {
            observed += count as f64 / N as f64;
            expected += rate / total;
            assert!((observed - expected).abs() <= bound);
        }
    }

    #[test]
    fn large_groups_draw_the_member_and_the_acceptance_separately() {
        // 1,500 channels in one group (above JOINT_MEMBERS) and one alone.
        let mut rates = vec![1.5; 1500];
        rates.push(1000.0);
        let mut tree = RateTree::default();
        let mut selector = Composition::default();
        tree.reset(rates.len());
        selector.reset(rates.len());
        for (i, &rate) in rates.iter().enumerate() {
            selector.update(i, rate).ok().unwrap();
            tree.set(i, rate, rate);
        }
        let total = tree.total().ok().unwrap();
        let mut rng = crate::sampling::stream(4242, 778);
        let mut words = Words::for_test(&mut rng, 20_000_000);
        const N: u64 = 400_000;
        let mut big = 0u64;
        for _ in 0..N {
            if selector.sample(&tree, &mut words).ok().unwrap() == 1500 {
                big += 1;
            }
        }
        let (p, n) = (1000.0 / total, N as f64);
        let observed = big as f64 / n;
        assert!(
            (observed - p).abs() <= 5.0 * (p * (1.0 - p) / n).sqrt(),
            "{observed} {p}"
        );
    }

    struct Fixed(Vec<u64>);
    impl RngCore for Fixed {
        fn next_u64(&mut self) -> u64 {
            self.0.remove(0)
        }
        fn next_u32(&mut self) -> u32 {
            unreachable!()
        }
        fn fill_bytes(&mut self, _: &mut [u8]) {
            unreachable!()
        }
        fn try_fill_bytes(&mut self, _: &mut [u8]) -> Result<(), rand_core::Error> {
            unreachable!()
        }
    }

    #[test]
    fn sorting_selection_redraws_a_threshold_beyond_the_sequential_sum() {
        // The pairwise total exceeds the ordered sum by one ulp.
        let rates = [
            1.5667571478300881,
            0.2652157528648588,
            1.4018558426681513,
            0.6054530008643468,
            0.40195711670041256,
        ];
        let mut tree = RateTree::default();
        tree.reset(rates.len());
        for (i, &rate) in rates.iter().enumerate() {
            tree.set(i, rate, rate);
        }
        let total = tree.total().ok().unwrap();
        let sequential = rates.iter().fold(0.0, |sum, &rate| sum + rate);
        assert!(sequential < total, "{sequential} {total}");
        // The largest 53-bit uniform overruns; the next draw selects channel 0.
        let mut rng = Fixed(vec![u64::MAX, 0]);
        let mut words = Words::for_test(&mut rng, 10);
        let mut order: Vec<usize> = (0..rates.len()).collect();
        let chosen = sorting_sample(&mut order, &rates, total, &mut words);
        assert_eq!(chosen.ok(), Some(0));
        assert_eq!(words.used(), 2);
    }

    #[test]
    fn sorting_moves_the_fired_channel_forward() {
        let rates = [1.0, 1.0, 1.0];
        let mut order = vec![0, 1, 2];
        // A threshold in the third channel's interval: 2.5 of 3.
        let word = ((2.5 / 3.0) * (1u64 << 53) as f64) as u64;
        let mut rng = Fixed(vec![word << 11]);
        let mut words = Words::for_test(&mut rng, 10);
        assert_eq!(
            sorting_sample(&mut order, &rates, 3.0, &mut words).ok(),
            Some(2)
        );
        assert_eq!(order, [0, 2, 1]);
    }

    #[test]
    fn exponent_groups_cover_normal_subnormal_and_extreme_values() {
        for x in [
            f64::from_bits(1),
            f64::from_bits(15),
            f64::MIN_POSITIVE,
            0.5,
            1.0,
            1.999,
            2.0,
            f64::MAX,
        ] {
            let top = ceiling(group(x));
            assert!(top >= x && top / 2.0 <= x);
        }
    }

    /// The previous implementation: every node carries its checked state.
    struct Reference {
        nodes: Vec<Node>,
        base: usize,
    }

    impl Reference {
        fn new(count: usize) -> Self {
            let base = count.max(1).next_power_of_two();
            Self {
                nodes: vec![Node::default(); 2 * base],
                base,
            }
        }
        fn set(&mut self, index: usize, rate: f64, minimum_positive: f64) {
            let mut i = self.base + index;
            self.nodes[i] = leaf(rate, minimum_positive);
            while i > 1 {
                i /= 2;
                self.nodes[i] = parent(self.nodes[2 * i], self.nodes[2 * i + 1]);
            }
        }
        fn total(&self) -> Option<f64> {
            let root = self.nodes[1];
            (!(root.bad
                || (root.minimum.is_finite()
                    && root.minimum / root.sum < 1.0 / (1u64 << 53) as f64)))
                .then_some(root.sum)
        }
    }

    #[test]
    fn lazily_checked_tree_decides_exactly_like_the_fully_checked_tree() {
        let mut rng = crate::sampling::stream(4242, 31);
        let specials = [
            0.0,
            f64::MIN_POSITIVE,
            f64::from_bits(1),
            1.0,
            1.0 + f64::EPSILON,
            0.5 * f64::EPSILON,
            f64::EPSILON,
            2.0f64.powi(-53),
            f64::MAX,
            f64::INFINITY,
            f64::NAN,
            -1.0,
            1e-300,
            1e300,
        ];
        let rate_of_index = |index: usize| (index + 1) as f64;
        let (mut accepted, mut rejected) = (0, 0);
        for case in 0..400 {
            let count = 1 + (rng.next_u64() % 40) as usize;
            let mut tree = RateTree::default();
            tree.reset(count);
            let mut reference = Reference::new(count);
            for _ in 0..200 {
                // One to three leaves per commit, as one event's updates.
                let updates = 1 + (rng.next_u64() % 3) as usize;
                for _ in 0..updates {
                    let index = (rng.next_u64() % count as u64) as usize;
                    let word = rng.next_u64();
                    let wild = case % 4 == 0;
                    let span: u64 = if case % 5 == 0 { 120 } else { 40 };
                    let scale = 2.0f64.powi(-(case % 3) * 330);
                    let rate = match word % 8 {
                        0 if wild => specials[(word >> 8) as usize % specials.len()],
                        0 => rate_of_index(index) * scale,
                        1 => 0.0,
                        // Wide dynamic ranges, sometimes beyond 53-bit resolution.
                        2 | 3 => {
                            f64::from_bits(0x3ff0_0000_0000_0000 + (word >> 12))
                                * 2.0f64.powi(((word >> 3) % span) as i32 - span as i32 / 2)
                                * scale
                        }
                        _ => (word >> 11) as f64 / (1u64 << 53) as f64 * 10.0 * scale,
                    };
                    // Mostly the rate itself; sometimes a smaller or odd bound.
                    let minimum = match (word >> 60) % 6 {
                        0 => rate * 0.25,
                        1 if wild => specials[(word >> 20) as usize % specials.len()],
                        _ if rate > 0.0 => rate,
                        _ => f64::INFINITY,
                    };
                    tree.set(index, rate, minimum);
                    reference.set(index, rate, minimum);
                    assert_eq!(tree.get(index).to_bits(), rate.to_bits());
                }
                let expected = reference.total();
                match (tree.total(), expected) {
                    (Ok(a), Some(b)) => {
                        assert_eq!(a.to_bits(), b.to_bits());
                        accepted += 1;
                    }
                    (Err(_), None) => rejected += 1,
                    (a, b) => panic!("case {case}: {:?} versus {b:?}", a.ok()),
                }
                // The committed sums are the reference's.
                for (node, reference) in reference.nodes.iter().enumerate().skip(1) {
                    assert_eq!(tree.sums[node].to_bits(), reference.sum.to_bits());
                }
            }
        }
        assert!(
            accepted > 10_000 && rejected > 10_000,
            "{accepted} {rejected}"
        );
    }

    #[test]
    fn tree_rebuilds_without_incremental_cancellation_drift() {
        let mut tree = RateTree::default();
        tree.reset(3);
        for _ in 0..1000 {
            tree.set(0, 1e12, 1e12);
            tree.set(1, 0.1, 0.1);
            tree.set(2, 0.2, 0.2);
            tree.set(0, 0.0, f64::INFINITY);
            assert_eq!(tree.total().ok().unwrap(), 0.1 + 0.2);
        }
        tree.set(0, 1e100, 1e100);
        let Err(Halt::Numerical(error)) = tree.total() else {
            panic!("resolved");
        };
        assert_eq!(error.code, "CRN-R031");
    }

    #[test]
    fn a_single_channel_tree_is_its_leaf() {
        let mut tree = RateTree::default();
        tree.reset(1);
        tree.set(0, 3.5, 3.5);
        assert_eq!(tree.total().ok(), Some(3.5));
        tree.set(0, 0.0, f64::INFINITY);
        assert_eq!(tree.total().ok(), Some(0.0));
    }
}
