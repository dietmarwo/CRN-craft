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
