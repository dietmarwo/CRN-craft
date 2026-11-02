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

/// Pairwise sums rebuild the ancestor chain; no accumulated add/subtract drift.
#[derive(Default)]
pub(super) struct RateTree {
    nodes: Vec<Node>,
    base: usize,
    len: usize,
}
impl RateTree {
    pub(super) fn reset(&mut self, count: usize) {
        self.len = count;
        self.base = count.max(1).next_power_of_two();
        self.nodes.resize(self.base * 2, Node::default());
        self.nodes.fill(Node::default());
    }
    pub(super) fn get(&self, index: usize) -> f64 {
        self.nodes[self.base + index].sum
    }
    pub(super) fn set(&mut self, index: usize, rate: f64, minimum_positive: f64) {
        let mut i = self.base + index;
        self.nodes[i] = Node {
            sum: rate,
            minimum: minimum_positive,
            bad: !rate.is_finite()
                || rate < 0.0
                || minimum_positive <= 0.0
                || minimum_positive.is_nan(),
        };
        while i > 1 {
            i /= 2;
            let a = self.nodes[2 * i];
            let b = self.nodes[2 * i + 1];
            let sum = a.sum + b.sum;
            self.nodes[i] = Node {
                sum,
                minimum: a.minimum.min(b.minimum),
                bad: a.bad
                    || b.bad
                    || !sum.is_finite()
                    || (a.sum > 0.0 && sum <= b.sum)
                    || (b.sum > 0.0 && sum <= a.sum),
            };
        }
    }
    pub(super) fn total(&self) -> Result<f64, Diagnostic> {
        let root = self.nodes[1];
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
    pub(super) fn sample(
        &self,
        tree: &RateTree,
        rng: &mut impl RngCore,
        draws: &mut Draws,
    ) -> Result<usize, Bail> {
        let largest = *self
            .active
            .last()
            .ok_or_else(|| numeric("empty composition selection"))?;
        let scale = ceiling(largest);
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
        }
        loop {
            draws.work.proposals += 1;
            let threshold = draws.closed(rng)? * sum;
            if threshold >= sum {
                return Err(numeric("composition threshold rounded to its total").into());
            }
            let mut cumulative = 0.0;
            let mut selected = None;
            for (k, &id) in self.active.iter().enumerate() {
                if k % 64 == 0 {
                    draws.poll()?;
                }
                cumulative += (ceiling(id) / scale) * self.groups[id].len() as f64;
                if threshold < cumulative {
                    selected = Some(id);
                    break;
                }
            }
            let id = selected.ok_or_else(|| numeric("composition group selection failed"))?;
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
