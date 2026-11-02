//! Composition/rejection with power-of-two envelopes and deterministic sums.
//! Algorithm family: Slepoy, Thompson & Plimpton, JCP 128, 205101 (2008).
//! This implementation chooses groups by envelope mass and rejects globally.
//! The accepted index has probability a_i / sum(a); waits use the actual sum.
use super::{Bail, Draws, numeric};
use crate::diagnostic::Diagnostic;
use rand_core::RngCore;

/// McCollum et al., Computational Biology and Chemistry 30 (2006), 39–49.
/// The fired channel moves one position forward; the order resets per path.
pub(super) fn sorting_sample(
    order: &mut [usize],
    rates: &[f64],
    total: f64,
    rng: &mut impl RngCore,
    draws: &mut Draws,
) -> Result<usize, Bail> {
    let threshold = draws.closed(rng)? * total;
    if threshold >= total {
        return Err(numeric("sorting threshold rounded to total").into());
    }
    let mut sum = 0.0;
    for (position, &channel) in order.iter().enumerate() {
        if position % 64 == 0 {
            draws.poll()?;
        }
        let next = sum + rates[channel];
        if rates[channel] > 0.0 && next <= sum {
            return Err(numeric("sorting sum loses a positive channel").into());
        }
        sum = next;
        if threshold < sum {
            if position > 0 {
                order.swap(position, position - 1);
            }
            return Ok(channel);
        }
    }
    Err(numeric("sorting selection cannot resolve the checked total").into())
}

#[derive(Clone, Copy)]
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

/// Pairwise sums rebuild the ancestor chain; no accumulated add/subtract drift.
///
/// Only the sums are maintained per node. `total` must decide exactly as a
/// tree that also carries per-node minima and loss flags (`Node`): reject an
/// invalid leaf, a non-finite sum, a pairwise addition that loses a positive
/// operand, or a smallest positive leaf below 2^-53 of the total. `floor` never
/// exceeds a current positive leaf or leaf minimum. When it is at least 2^-52
/// of a normal total, every positive pairwise operand exceeds half an ulp of
/// any sum (all sums are at most the total), so no addition can lose it, and
/// the resolution test passes; otherwise the checked tree is evaluated in full.
#[derive(Default)]
pub(super) struct RateTree {
    sums: Vec<f64>,
    minima: Vec<f64>,
    base: usize,
    len: usize,
    invalid: usize,
    floor: f64,
    checked: Vec<Node>,
    level: Vec<usize>,
}
impl RateTree {
    pub(super) fn reset(&mut self, count: usize) {
        self.len = count;
        self.base = count.max(1).next_power_of_two();
        self.sums.resize(self.base * 2, 0.0);
        self.sums.fill(0.0);
        self.minima.resize(self.base, f64::INFINITY);
        self.minima.fill(f64::INFINITY);
        self.invalid = 0;
        self.floor = f64::INFINITY;
    }
    pub(super) fn get(&self, index: usize) -> f64 {
        self.sums[self.base + index]
    }
    /// Writes one leaf; its ancestors are stale until `rebuild` covers it.
    pub(super) fn set_leaf(&mut self, index: usize, rate: f64, minimum_positive: f64) {
        let i = self.base + index;
        self.invalid -= usize::from(leaf(self.sums[i], self.minima[index]).bad);
        self.invalid += usize::from(leaf(rate, minimum_positive).bad);
        self.sums[i] = rate;
        self.minima[index] = minimum_positive;
        let positive = if rate > 0.0 { rate } else { f64::INFINITY };
        self.floor = self.floor.min(positive.min(minimum_positive));
    }
    fn chain(&mut self, index: usize) {
        let mut i = self.base + index;
        while i > 1 {
            i /= 2;
            self.sums[i] = self.sums[2 * i] + self.sums[2 * i + 1];
        }
    }
    #[cfg(test)]
    pub(super) fn set(&mut self, index: usize, rate: f64, minimum_positive: f64) {
        self.set_leaf(index, rate, minimum_positive);
        self.chain(index);
    }
    pub(super) fn rebuild_all(&mut self) {
        for i in (1..self.base).rev() {
            self.sums[i] = self.sums[2 * i] + self.sums[2 * i + 1];
        }
    }
    /// Recomputes every ancestor of the given leaves once, a level at a time.
    /// Each node is the sum of its two children, so the result equals one
    /// `set` per leaf. Ascending leaves keep each level free of duplicates.
    pub(super) fn rebuild(&mut self, leaves: &[usize]) {
        if leaves.len() <= 4 {
            for &k in leaves {
                self.chain(k);
            }
            return;
        }
        if 2 * leaves.len() >= self.len {
            self.rebuild_all();
            return;
        }
        let mut level = std::mem::take(&mut self.level);
        level.clear();
        for &k in leaves {
            let p = (self.base + k) / 2;
            if level.last() != Some(&p) {
                level.push(p);
            }
        }
        loop {
            for &i in &level {
                self.sums[i] = self.sums[2 * i] + self.sums[2 * i + 1];
            }
            if level.first() == Some(&1) {
                break;
            }
            let mut kept = 0;
            for k in 0..level.len() {
                let p = level[k] / 2;
                if kept == 0 || level[kept - 1] != p {
                    level[kept] = p;
                    kept += 1;
                }
            }
            level.truncate(kept);
        }
        self.level = level;
    }
    pub(super) fn total(&mut self) -> Result<f64, Diagnostic> {
        let root = self.sums[1];
        if self.invalid == 0 && root.is_normal() && self.floor >= root * (1.0 / (1u64 << 52) as f64)
        {
            return Ok(root);
        }
        self.checked_total()
    }
    #[cold]
    fn checked_total(&mut self) -> Result<f64, Diagnostic> {
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
        if root.bad
            || (root.minimum.is_finite() && root.minimum / root.sum < 1.0 / (1u64 << 53) as f64)
        {
            return Err(numeric(
                "pairwise total is non-finite or cannot resolve a positive channel",
            ));
        }
        Ok(root.sum)
    }
}

const GROUPS: usize = 2098;
fn group(rate: f64) -> usize {
    let bits = rate.to_bits();
    let exponent = (bits >> 52) as usize;
    if exponent == 0 {
        (63 - bits.leading_zeros()) as usize
    } else {
        exponent + 51
    }
}
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

#[derive(Default)]
pub(super) struct Composition {
    groups: Vec<Vec<usize>>,
    active: Vec<usize>,
    locations: Vec<(usize, usize)>,
    /// Running envelope sums over `active`, rebuilt only after a membership
    /// change with the same additions in the same order as a fresh scan.
    cumulative: Vec<f64>,
    stale: bool,
}
impl Composition {
    pub(super) fn reset(&mut self, count: usize) {
        if self.groups.is_empty() {
            self.groups.resize_with(GROUPS, Vec::new);
        }
        for id in self.active.drain(..) {
            self.groups[id].clear();
            self.groups[id].shrink_to_fit();
        }
        self.locations.resize(count, (usize::MAX, 0));
        self.locations.fill((usize::MAX, 0));
        self.stale = true;
    }
    pub(super) fn update(&mut self, index: usize, value: f64) -> Result<(), Diagnostic> {
        if !value.is_finite() || value < 0.0 {
            return Err(numeric("invalid composition rate"));
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
            // Bound retained capacity by live memberships, rather than by all
            // exponent buckets a reaction has visited over a long trajectory.
            if self.groups[old].capacity() > 4 * self.groups[old].len() {
                self.groups[old].shrink_to_fit();
            }
            if self.groups[old].is_empty() {
                let at = self.active.binary_search(&old).unwrap();
                self.active.remove(at);
            }
        }
        if next == usize::MAX {
            self.locations[index] = (next, 0);
        } else {
            if self.groups[next].is_empty() {
                let at = self.active.binary_search(&next).unwrap_err();
                self.active.insert(at, next);
            }
            self.locations[index] = (next, self.groups[next].len());
            self.groups[next].push(index);
        }
        Ok(())
    }
    fn rebuild(&mut self, draws: &mut Draws) -> Result<(), Bail> {
        let largest = *self
            .active
            .last()
            .ok_or_else(|| numeric("empty composition selection"))?;
        let scale = ceiling(largest);
        self.cumulative.clear();
        let mut sum = 0.0;
        for (k, &id) in self.active.iter().enumerate() {
            if k % 64 == 0 {
                draws.poll()?;
            }
            let weight = (ceiling(id) / scale) * self.groups[id].len() as f64;
            let next = sum + weight;
            if weight <= 0.0 || next <= sum {
                return Err(numeric("composition envelope loses a positive group").into());
            }
            sum = next;
            self.cumulative.push(sum);
        }
        self.stale = false;
        Ok(())
    }
    pub(super) fn sample(
        &mut self,
        tree: &RateTree,
        rng: &mut impl RngCore,
        draws: &mut Draws,
    ) -> Result<usize, Bail> {
        if self.stale {
            self.rebuild(draws)?;
        }
        let sum = *self
            .cumulative
            .last()
            .ok_or_else(|| numeric("empty composition selection"))?;
        loop {
            draws.work.proposals += 1;
            let threshold = draws.closed(rng)? * sum;
            if threshold >= sum {
                return Err(numeric("composition threshold rounded to its total").into());
            }
            // First running sum above the threshold, as in a sequential scan.
            let k = self.cumulative.partition_point(|&c| c <= threshold);
            let id = *self
                .active
                .get(k)
                .ok_or_else(|| numeric("composition group selection failed"))?;
            let members = &self.groups[id];
            let index = members[draws.index(rng, members.len())?];
            let bound = ceiling(id);
            let threshold = draws.closed(rng)? * bound;
            if threshold >= bound {
                return Err(numeric("composition acceptance rounded to its bound").into());
            }
            if threshold < tree.get(index) {
                return Ok(index);
            }
            draws.poll()?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn composition_matches_declared_categorical_dkw_bound_after_membership_changes() {
        let rates = [0.0, 0.125, 0.5, 1.0, 2.0, 7.0, 8.0, 50.0];
        let mut tree = RateTree::default();
        let mut selector = Composition::default();
        tree.reset(rates.len());
        selector.reset(rates.len());
        // Exercise removal, group migration, zero activation and swap-removal.
        for (i, &rate) in rates.iter().enumerate() {
            selector.update(i, 8.0).unwrap();
            selector.update(i, rate).unwrap();
            tree.set(i, rate, if rate > 0.0 { rate } else { f64::INFINITY });
        }
        let total = tree.total().unwrap();
        let control = crate::simulation::Control::default();
        let mut poller = crate::simulation::Poller::new(&control, std::time::Instant::now(), None);
        let mut work = super::super::Work::default();
        let mut draws = Draws {
            limit: 20_000_000,
            work: &mut work,
            poller: &mut poller,
        };
        let mut rng = crate::sampling::stream(4242, 777);
        let mut counts = [0u64; 8];
        const N: u64 = 1_000_000;
        for _ in 0..N {
            let index = selector.sample(&tree, &mut rng, &mut draws).ok().unwrap();
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
                let expected = reference.total();
                match (tree.total(), expected) {
                    (Ok(a), Some(b)) => {
                        assert_eq!(a.to_bits(), b.to_bits());
                        accepted += 1;
                    }
                    (Err(_), None) => rejected += 1,
                    (a, b) => panic!("case {case}: {a:?} versus {b:?}"),
                }
                assert_eq!(tree.get(index).to_bits(), rate.to_bits());
            }
        }
        assert!(
            accepted > 20_000 && rejected > 20_000,
            "{accepted} {rejected}"
        );
    }

    #[test]
    fn batched_ancestor_rebuilds_equal_one_update_per_leaf() {
        let mut rng = crate::sampling::stream(4242, 32);
        for count in [1, 2, 3, 5, 8, 9, 31, 64, 100, 257] {
            let (mut single, mut batched) = (RateTree::default(), RateTree::default());
            single.reset(count);
            batched.reset(count);
            for _ in 0..200 {
                let mut leaves: Vec<usize> = (0..count)
                    .filter(|_| rng.next_u64() % (1 + count as u64 / 4) == 0)
                    .collect();
                if rng.next_u64() % 5 == 0 {
                    leaves = (0..count).collect();
                }
                for &k in &leaves {
                    let rate = (rng.next_u64() >> 11) as f64 / (1u64 << 40) as f64;
                    single.set(k, rate, if rate > 0.0 { rate } else { f64::INFINITY });
                    batched.set_leaf(k, rate, if rate > 0.0 { rate } else { f64::INFINITY });
                }
                batched.rebuild(&leaves);
                let bits = |t: &RateTree| t.sums.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                assert_eq!(bits(&single), bits(&batched), "{count} {leaves:?}");
            }
        }
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
            assert_eq!(tree.total().unwrap(), 0.1 + 0.2);
        }
        tree.set(0, 1e100, 1e100);
        assert!(tree.total().is_err());
    }
}
