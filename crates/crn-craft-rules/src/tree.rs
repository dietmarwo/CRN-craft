//! Complete binary trees over instance slots.
//!
//! Internal nodes are recomputed from their children on every leaf update, so
//! totals carry no accumulated drift. Capacity is a power of two and doubles
//! when a type gains more slots; the tree is then rebuilt from its leaves.

/// Pairwise sums with optional minimum positive leaf (for the resolution
/// check of selection trees).
#[derive(Debug, Clone)]
pub(crate) struct SumTree {
    capacity: usize,
    sums: Vec<f64>,
    min_positive: Option<Vec<f64>>,
}

impl SumTree {
    pub(crate) fn new(track_min: bool) -> Self {
        Self {
            capacity: 1,
            sums: vec![0.0; 2],
            min_positive: track_min.then(|| vec![f64::INFINITY; 2]),
        }
    }

    pub(crate) fn total(&self) -> f64 {
        self.sums[1]
    }

    /// Smallest positive leaf, or infinity when none.
    pub(crate) fn min_positive(&self) -> f64 {
        self.min_positive.as_ref().map_or(f64::INFINITY, |m| m[1])
    }

    #[cfg(test)]
    pub(crate) fn leaf(&self, slot: usize) -> f64 {
        if slot < self.capacity {
            self.sums[self.capacity + slot]
        } else {
            0.0
        }
    }

    fn grow(&mut self, slots: usize) {
        if slots <= self.capacity {
            return;
        }
        let mut capacity = self.capacity;
        while capacity < slots {
            capacity *= 2;
        }
        let mut sums = vec![0.0; 2 * capacity];
        sums[capacity..capacity + self.capacity]
            .copy_from_slice(&self.sums[self.capacity..2 * self.capacity]);
        let mut min = self.min_positive.as_ref().map(|old| {
            let mut m = vec![f64::INFINITY; 2 * capacity];
            m[capacity..capacity + self.capacity]
                .copy_from_slice(&old[self.capacity..2 * self.capacity]);
            m
        });
        for i in (1..capacity).rev() {
            sums[i] = sums[2 * i] + sums[2 * i + 1];
            if let Some(m) = min.as_mut() {
                m[i] = m[2 * i].min(m[2 * i + 1]);
            }
        }
        self.capacity = capacity;
        self.sums = sums;
        self.min_positive = min;
    }

    pub(crate) fn set(&mut self, slot: usize, value: f64) {
        self.grow(slot + 1);
        let mut i = self.capacity + slot;
        self.sums[i] = value;
        if let Some(m) = self.min_positive.as_mut() {
            m[i] = if value > 0.0 { value } else { f64::INFINITY };
        }
        while i > 1 {
            i /= 2;
            self.sums[i] = self.sums[2 * i] + self.sums[2 * i + 1];
            if let Some(m) = self.min_positive.as_mut() {
                m[i] = m[2 * i].min(m[2 * i + 1]);
            }
        }
    }

    /// Slot whose cumulative interval contains `target` in `[0, total)`.
    /// Never returns a zero leaf while the total is positive.
    pub(crate) fn select(&self, mut target: f64) -> usize {
        let mut i = 1;
        while i < self.capacity {
            let left = self.sums[2 * i];
            let right = self.sums[2 * i + 1];
            if target < left || right <= 0.0 {
                i *= 2;
            } else {
                target -= left;
                i = 2 * i + 1;
            }
        }
        i - self.capacity
    }
}

/// Pair weights of one pattern over its instances. Each node stores the
/// pattern weight S = Σ h·n of its leaves, the weight of unordered pairs of
/// distinct molecules inside it, P = P_left + P_right + S_left·S_right (a
/// leaf holds h²·n(n−1)/2), and the smallest positive local factor h. Pair
/// totals are sums of nonnegative terms, without the cancellation of
/// ½[(Σ h n)² − Σ h² n].
#[derive(Debug, Clone)]
pub(crate) struct PairTree {
    capacity: usize,
    sums: Vec<f64>,
    pairs: Vec<f64>,
    min_factor: Vec<f64>,
}

/// Where a selected pair lies.
pub(crate) enum PairChoice {
    /// Two molecules of one instance.
    Same(usize),
    /// Molecules of two subtrees: (node of the first, node of the second).
    Across(usize, usize),
}

impl PairTree {
    pub(crate) fn new() -> Self {
        Self {
            capacity: 1,
            sums: vec![0.0; 2],
            pairs: vec![0.0; 2],
            min_factor: vec![f64::INFINITY; 2],
        }
    }

    pub(crate) fn total_weight(&self) -> f64 {
        self.sums[1]
    }

    pub(crate) fn total_pairs(&self) -> f64 {
        self.pairs[1]
    }

    /// Smallest positive h over instances that are present.
    pub(crate) fn min_factor(&self) -> f64 {
        self.min_factor[1]
    }

    fn combine(&mut self, i: usize) {
        let (l, r) = (2 * i, 2 * i + 1);
        self.sums[i] = self.sums[l] + self.sums[r];
        self.pairs[i] = self.pairs[l] + self.pairs[r] + self.sums[l] * self.sums[r];
        self.min_factor[i] = self.min_factor[l].min(self.min_factor[r]);
    }

    fn grow(&mut self, slots: usize) {
        if slots <= self.capacity {
            return;
        }
        let mut capacity = self.capacity;
        while capacity < slots {
            capacity *= 2;
        }
        let old = self.capacity;
        let mut next = Self {
            capacity,
            sums: vec![0.0; 2 * capacity],
            pairs: vec![0.0; 2 * capacity],
            min_factor: vec![f64::INFINITY; 2 * capacity],
        };
        next.sums[capacity..capacity + old].copy_from_slice(&self.sums[old..2 * old]);
        next.pairs[capacity..capacity + old].copy_from_slice(&self.pairs[old..2 * old]);
        next.min_factor[capacity..capacity + old].copy_from_slice(&self.min_factor[old..2 * old]);
        for i in (1..capacity).rev() {
            next.combine(i);
        }
        *self = next;
    }

    /// Leaf for an instance with factor `h` and count `n`.
    #[allow(clippy::cast_precision_loss)]
    pub(crate) fn set(&mut self, slot: usize, h: f64, n: u64) {
        self.grow(slot + 1);
        let mut i = self.capacity + slot;
        let count = n as f64;
        self.sums[i] = h * count;
        // n(n−1)/2 is exact for n ≤ 2^26 and correctly rounded beyond.
        self.pairs[i] = if n >= 2 {
            h * h * (count * (count - 1.0) * 0.5)
        } else {
            0.0
        };
        self.min_factor[i] = if n > 0 && h > 0.0 { h } else { f64::INFINITY };
        while i > 1 {
            i /= 2;
            self.combine(i);
        }
    }

    /// Locates the pair for `target` in [0, total_pairs).
    pub(crate) fn select(&self, mut target: f64) -> PairChoice {
        let mut i = 1;
        while i < self.capacity {
            let (l, r) = (2 * i, 2 * i + 1);
            let across = self.sums[l] * self.sums[r];
            if target < self.pairs[l] {
                i = l;
                continue;
            }
            target -= self.pairs[l];
            if target < across || self.pairs[r] <= 0.0 {
                if across > 0.0 {
                    return PairChoice::Across(l, r);
                }
                i = l;
                continue;
            }
            target -= across;
            i = r;
        }
        PairChoice::Same(i - self.capacity)
    }

    /// Weight of a subtree node, for the second stage of an `Across` choice.
    pub(crate) fn node_weight(&self, node: usize) -> f64 {
        self.sums[node]
    }

    /// Leaf below `node` holding the weighted molecule at `target` in
    /// [0, node_weight(node)).
    pub(crate) fn descend(&self, mut node: usize, mut target: f64) -> usize {
        while node < self.capacity {
            let (l, r) = (2 * node, 2 * node + 1);
            if target < self.sums[l] || self.sums[r] <= 0.0 {
                node = l;
            } else {
                target -= self.sums[l];
                node = r;
            }
        }
        node - self.capacity
    }
}

/// Maximum (or, negated, minimum) over present leaves.
#[derive(Debug, Clone)]
pub(crate) struct MaxTree {
    capacity: usize,
    nodes: Vec<f64>,
}

impl MaxTree {
    pub(crate) fn new() -> Self {
        Self {
            capacity: 1,
            nodes: vec![f64::NEG_INFINITY; 2],
        }
    }

    pub(crate) fn top(&self) -> f64 {
        self.nodes[1]
    }

    fn grow(&mut self, slots: usize) {
        if slots <= self.capacity {
            return;
        }
        let mut capacity = self.capacity;
        while capacity < slots {
            capacity *= 2;
        }
        let mut nodes = vec![f64::NEG_INFINITY; 2 * capacity];
        nodes[capacity..capacity + self.capacity]
            .copy_from_slice(&self.nodes[self.capacity..2 * self.capacity]);
        for i in (1..capacity).rev() {
            nodes[i] = nodes[2 * i].max(nodes[2 * i + 1]);
        }
        self.capacity = capacity;
        self.nodes = nodes;
    }

    /// `None` removes the leaf.
    pub(crate) fn set(&mut self, slot: usize, value: Option<f64>) {
        self.grow(slot + 1);
        let mut i = self.capacity + slot;
        self.nodes[i] = value.unwrap_or(f64::NEG_INFINITY);
        while i > 1 {
            i /= 2;
            self.nodes[i] = self.nodes[2 * i].max(self.nodes[2 * i + 1]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sums_grow_and_select() {
        let mut t = SumTree::new(true);
        for (slot, w) in [(0, 1.0), (5, 2.0), (2, 0.5), (9, 4.0)] {
            t.set(slot, w);
        }
        assert_eq!(t.total(), 7.5);
        assert_eq!(t.min_positive(), 0.5);
        assert_eq!(t.select(0.0), 0);
        assert_eq!(t.select(0.99), 0);
        assert_eq!(t.select(1.0), 2);
        assert_eq!(t.select(1.5), 5);
        assert_eq!(t.select(3.49), 5);
        assert_eq!(t.select(3.5), 9);
        assert_eq!(t.select(7.4999), 9);
        // A target at the total (rounding) never lands on an empty leaf.
        assert_eq!(t.select(7.5), 9);
        t.set(9, 0.0);
        assert_eq!(t.select(3.5), 5);
        assert_eq!(t.leaf(9), 0.0);
        assert_eq!(t.leaf(100), 0.0);
        t.set(2, 0.0);
        assert_eq!(t.min_positive(), 1.0);
    }

    #[test]
    fn pair_weights_and_selection() {
        let mut t = PairTree::new();
        // Factors 1 (count 3), 2 (count 1), 5 (count 2).
        t.set(0, 1.0, 3);
        t.set(1, 2.0, 1);
        t.set(4, 5.0, 2);
        let (s, q) = (3.0 + 2.0 + 10.0, 3.0 + 4.0 + 50.0);
        assert_eq!(t.total_weight(), s);
        assert_eq!(t.total_pairs(), 0.5 * (s * s - q));
        assert_eq!(t.min_factor(), 1.0);
        // Two-stage selection as in the engine, against exact pair weights
        // h(x)h(y) over unordered pairs of distinct molecules.
        use rand_core::RngCore;
        let mut rng = crate::stream(5, 0);
        let mut u = || (rng.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0;
        let draws = 400_000;
        let mut counts = std::collections::BTreeMap::new();
        for _ in 0..draws {
            let pair = match t.select(u() * t.total_pairs()) {
                PairChoice::Same(i) => (i, i),
                PairChoice::Across(a, b) => {
                    let i = t.descend(a, u() * t.node_weight(a));
                    let j = t.descend(b, u() * t.node_weight(b));
                    (i.min(j), i.max(j))
                }
            };
            *counts.entry(pair).or_insert(0u32) += 1;
        }
        let total = t.total_pairs();
        // Slots 0 (h 1, n 3), 1 (h 2, n 1), 4 (h 5, n 2).
        for (pair, weight) in [
            ((0, 0), 3.0),
            ((0, 1), 6.0),
            ((0, 4), 30.0),
            ((1, 4), 20.0),
            ((4, 4), 25.0),
        ] {
            let p = weight / total;
            let f = f64::from(counts.get(&pair).copied().unwrap_or(0)) / f64::from(draws);
            let se = (p * (1.0 - p) / f64::from(draws)).sqrt();
            assert!((f - p).abs() < 5.0 * se, "{pair:?}: {f} vs {p}");
        }
        assert_eq!(counts.len(), 5);
        t.set(0, 1.0, 0);
        assert_eq!(t.min_factor(), 2.0);
        t.set(1, 2.0, 0);
        t.set(4, 5.0, 1);
        assert_eq!(t.total_pairs(), 0.0);
    }

    #[test]
    fn maxima() {
        let mut t = MaxTree::new();
        assert_eq!(t.top(), f64::NEG_INFINITY);
        t.set(3, Some(2.0));
        t.set(0, Some(5.0));
        t.set(7, Some(-1.0));
        assert_eq!(t.top(), 5.0);
        t.set(0, None);
        assert_eq!(t.top(), 2.0);
    }
}
