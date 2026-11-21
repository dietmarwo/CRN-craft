//! Exact local-event sampling of a two-genotype turnover chain.
use crate::diagnostic::Diagnostic;
use crate::rules::{self, Control};
use crate::sampling::uniform_closed_open;
use crate::ssa::exponential_from_words;
use rand_core::RngCore;
use serde::{Deserialize, Serialize};

const MAX_COUNT: u64 = 1 << 53;
const MAX_LABELS: u64 = 8_000_000;
// A positive weight must stay above 16 × 2⁻⁵³ of the total for the 53-bit
// selection draw (docs/developer/specialized-kernels.md).
const RESOLUTION: f64 = 16.0 / MAX_COUNT as f64;

/// A chain with count feedback, nearest-neighbor hopping and neutral clone IDs.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Turnover {
    /// Target wild-type copy count.
    pub target: f64,
    /// Coefficient of linear birth feedback with clipping at zero.
    pub feedback: f64,
    /// Mutant weight in the feedback expression.
    pub mutant_weight: f64,
    /// Per-copy loss rate.
    pub death: f64,
    /// Per-copy hop rate in each available direction.
    pub hop: f64,
    /// Probability that a birth makes a new mutant clone.
    pub mutation: f64,
    /// Keep the leftmost region's initial counts as a reservoir.
    pub fixed_left: bool,
}

/// Endpoint accounting of an exact chain simulation.
#[derive(Debug, Serialize)]
pub struct TurnoverResult {
    /// Effective sampler identity.
    pub sampler: &'static str,
    /// Endpoint wild-type population by region.
    pub wild: Vec<u64>,
    /// Endpoint mutant population by region.
    pub mutant: Vec<u64>,
    /// Mutant clone IDs by region when requested.
    pub labels: Option<Vec<Vec<u64>>>,
    /// Executed stochastic events.
    pub events: u64,
    /// Births creating a new clone.
    pub mutations: u64,
}

fn error_with(code: &'static str, message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        code,
        "/turnover",
        message,
        "Check chain rates, counts, time and event budget.",
    )
}

fn error(message: impl Into<String>) -> Diagnostic {
    error_with("CRN-R030", message)
}

struct Tree {
    base: usize,
    sums: Vec<f64>,
}

impl Tree {
    fn new(n: usize) -> Self {
        let base = n.next_power_of_two();
        Self {
            base,
            sums: vec![0.0; 2 * base],
        }
    }

    fn update(&mut self, site: usize, rate: f64) {
        let mut i = self.base + site;
        self.sums[i] = rate;
        while i > 1 {
            i /= 2;
            self.sums[i] = self.sums[2 * i] + self.sums[2 * i + 1];
        }
    }

    fn pick(&self, mut draw: f64) -> usize {
        let mut i = 1;
        while i < self.base {
            let left = i * 2;
            if draw < self.sums[left] || self.sums[left + 1] == 0.0 {
                i = left;
            } else {
                draw -= self.sums[left];
                i = left + 1;
            }
        }
        i - self.base
    }
}

struct Chain<'a> {
    spec: &'a Turnover,
    wild: Vec<u64>,
    mutant: Vec<u64>,
    labels: Option<Vec<Vec<u64>>>,
    tree: Tree,
    birth: Vec<f64>,
    rates: Vec<f64>,
    left: Vec<f64>,
    right: Vec<f64>,
    minimum_weight: f64,
    labels_live: u64,
}

impl<'a> Chain<'a> {
    fn new(spec: &'a Turnover, counts: &[[u64; 2]], labels: bool) -> Self {
        let n = counts.len();
        Self {
            spec,
            wild: counts.iter().map(|v| v[0]).collect(),
            mutant: counts.iter().map(|v| v[1]).collect(),
            labels: labels.then(|| counts.iter().map(|v| vec![0; v[1] as usize]).collect()),
            tree: Tree::new(n),
            birth: vec![0.0; n],
            rates: vec![0.0; n],
            left: (0..n).map(|i| if i > 0 { spec.hop } else { 0.0 }).collect(),
            minimum_weight: f64::INFINITY,
            labels_live: if labels {
                counts.iter().map(|v| v[1]).sum()
            } else {
                0
            },
            right: (0..n)
                .map(|i| if i + 1 < n { spec.hop } else { 0.0 })
                .collect(),
        }
    }

    fn reservoir(&self, i: usize) -> bool {
        self.spec.fixed_left && i == 0
    }

    #[allow(clippy::cast_precision_loss)]
    fn update(&mut self, i: usize) -> Result<(), Diagnostic> {
        let pop = self.wild[i] + self.mutant[i];
        let b = if self.reservoir(i) {
            0.0
        } else {
            (self.spec.death
                + self.spec.feedback
                    * (self.spec.target
                        - self.wild[i] as f64
                        - self.spec.mutant_weight * self.mutant[i] as f64))
                .max(0.0)
        };
        self.birth[i] = b;
        let total = (b + self.death_at(i)) + (self.left[i] + self.right[i]);
        self.rates[i] = total;
        let weight = pop as f64 * total;
        if pop > MAX_COUNT || !weight.is_finite() {
            return Err(error(
                "non-finite hazard or population exceeds the exact integer range",
            ));
        }
        if pop > 0 {
            let inv = 1.0 / total;
            for rate in [b, self.death_at(i), self.left[i], self.right[i]] {
                if rate > 0.0 && rate * inv < RESOLUTION {
                    return Err(error("a positive local event is below RNG resolution"));
                }
            }
        }
        if weight > 0.0 {
            self.minimum_weight = self.minimum_weight.min(weight);
        }
        self.tree.update(i, weight);
        Ok(())
    }

    fn death_at(&self, i: usize) -> f64 {
        if self.reservoir(i) {
            0.0
        } else {
            self.spec.death
        }
    }

    fn remove(&mut self, site: usize, copy: u64) -> (bool, u64) {
        if copy < self.mutant[site] {
            if self.reservoir(site) {
                return (true, 0);
            }
            self.mutant[site] -= 1;
            self.labels_live -= u64::from(self.labels.is_some());
            (
                true,
                self.labels
                    .as_mut()
                    .map_or(0, |labels| labels[site].swap_remove(copy as usize)),
            )
        } else {
            if !self.reservoir(site) {
                self.wild[site] -= 1;
            }
            (false, 0)
        }
    }

    fn add(&mut self, site: usize, mutant: bool, label: u64) -> Result<(), Diagnostic> {
        if self.reservoir(site) {
            return Ok(());
        }
        if let Some(labels) = &mut self.labels {
            if mutant {
                if self.labels_live >= MAX_LABELS {
                    return Err(error("clone storage budget exceeded"));
                }
                labels[site].push(label);
                self.labels_live += 1;
            }
        }
        if mutant {
            self.mutant[site] += 1;
        } else {
            self.wild[site] += 1;
        }
        Ok(())
    }
}

#[inline]
fn uniform(rng: &mut impl RngCore) -> f64 {
    uniform_closed_open(rng.next_u64())
}

/// Simulate one exact endpoint. Only the selected region and the destination
/// of a hop need their hazards updated. The virtual left reservoir exchanges
/// molecules at its fixed starting counts.
///
/// # Errors
/// Rejects unsupported rates, work limits, cancellations and populations.
#[allow(clippy::too_many_lines, clippy::cast_precision_loss)]
pub fn simulate(
    spec: &Turnover,
    initial: &[[u64; 2]],
    end: f64,
    seed: u64,
    max_events: u64,
    labels: bool,
    control: &Control,
) -> Result<TurnoverResult, Diagnostic> {
    if initial.is_empty()
        || initial.len() > 4096
        || !end.is_finite()
        || end < 0.0
        || ![
            spec.target,
            spec.feedback,
            spec.mutant_weight,
            spec.death,
            spec.hop,
            spec.mutation,
        ]
        .iter()
        .all(|&v| v.is_finite() && v >= 0.0)
        || spec.mutation > 1.0
        || initial
            .iter()
            .any(|&[w, m]| w > MAX_COUNT || m > MAX_COUNT - w || (labels && m > 2_000_000))
    {
        return Err(error("invalid chain setup"));
    }
    if labels && initial.iter().map(|v| u128::from(v[1])).sum::<u128>() > u128::from(MAX_LABELS) {
        return Err(error("initial clone storage exceeds the label budget"));
    }
    let mut chain = Chain::new(spec, initial, labels);
    let mut rng = rules::stream(seed, 0);
    for i in 0..initial.len() {
        chain.update(i)?;
    }
    let mut time = 0.0;
    let mut carry = 0.0;
    let mut events = 0;
    let mut mutations = 0;
    while time < end {
        if events % 4096 == 0 && control.is_cancelled() {
            return Err(error_with("CRN-R044", "chain simulation cancelled"));
        }
        let total = chain.tree.sums[1];
        if total == 0.0 {
            break;
        }
        if !total.is_normal() {
            return Err(error("invalid chain total hazard"));
        }
        if chain.minimum_weight < total * RESOLUTION {
            chain.minimum_weight = chain.tree.sums[chain.tree.base..]
                .iter()
                .filter(|&&rate| rate > 0.0)
                .copied()
                .fold(f64::INFINITY, f64::min);
            if chain.minimum_weight < total * RESOLUTION {
                return Err(error("a region is below RNG resolution"));
            }
        }
        let unit = exponential_from_words(|_| Ok::<_, std::convert::Infallible>(rng.next_u64()))
            .expect("random source is infallible");
        // Compensated waiting-time accumulation for very long paths. A
        // positive wait can be smaller than the ulp of the absolute clock.
        let increment = unit / total - carry;
        let next = time + increment;
        carry = (next - time) - increment;
        if next > end {
            break;
        }
        if events == max_events {
            return Err(error_with(
                "CRN-R040",
                format!("event budget {max_events} reached at time {time}"),
            ));
        }
        let i = chain.tree.pick(uniform(&mut rng) * total);
        let pop = chain.wild[i] + chain.mutant[i];
        let copy = ((uniform(&mut rng) * pop as f64) as u64).min(pop - 1);
        let draw = uniform(&mut rng) * chain.rates[i];
        let birth = chain.birth[i];
        let loss = birth + chain.death_at(i);
        if draw < birth {
            let new = spec.mutation > 0.0 && uniform(&mut rng) < spec.mutation;
            let (mutant, label) = if new {
                mutations += 1;
                (true, mutations)
            } else if copy < chain.mutant[i] {
                (
                    true,
                    chain
                        .labels
                        .as_ref()
                        .map_or(0, |labels| labels[i][copy as usize]),
                )
            } else {
                (false, 0)
            };
            chain.add(i, mutant, label)?;
        } else {
            let (mutant, label) = chain.remove(i, copy);
            if draw >= loss {
                let destination = if chain.right[i] == 0.0 || draw < loss + chain.left[i] {
                    i - 1
                } else {
                    i + 1
                };
                chain.add(destination, mutant, label)?;
                chain.update(destination)?;
            }
        }
        chain.update(i)?;
        time = next;
        events += 1;
    }
    Ok(TurnoverResult {
        sampler: "chacha8-u53-ziggurat/spatial-turnover/2",
        wild: chain.wild,
        mutant: chain.mutant,
        labels: chain.labels,
        events,
        mutations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(hop: f64) -> Turnover {
        Turnover {
            target: 20.0,
            feedback: 0.01,
            mutant_weight: 0.5,
            death: 0.07,
            hop,
            mutation: 0.0,
            fixed_left: false,
        }
    }

    fn hops_only() -> Turnover {
        Turnover {
            feedback: 0.0,
            death: 0.0,
            target: 0.0,
            ..spec(0.3)
        }
    }

    #[test]
    fn tree_pick_skips_zero_weights() {
        let mut tree = Tree::new(5);
        for (site, rate) in [0.0, 3.0, 0.0, 0.0, 5.0].into_iter().enumerate() {
            tree.update(site, rate);
        }
        assert_eq!(tree.sums[1], 8.0);
        assert_eq!(tree.pick(0.0), 1);
        assert_eq!(tree.pick(2.999), 1);
        assert_eq!(tree.pick(3.0), 4);
        assert_eq!(tree.pick(f64::from_bits(8.0_f64.to_bits() - 1)), 4);
    }

    #[test]
    fn hops_conserve_each_genotype_and_keep_labels() {
        let initial = [[7, 3], [0, 5], [4, 0], [2, 2]];
        let result = simulate(
            &hops_only(),
            &initial,
            50.0,
            11,
            1_000_000,
            true,
            &Control::default(),
        )
        .expect("hops only");
        assert!(result.events > 0);
        assert_eq!(result.wild.iter().sum::<u64>(), 13);
        assert_eq!(result.mutant.iter().sum::<u64>(), 10);
        let labels = result.labels.expect("labels were requested");
        for (ids, &count) in labels.iter().zip(&result.mutant) {
            assert_eq!(ids.len() as u64, count);
        }
        assert_eq!(result.mutations, 0);
    }

    #[test]
    fn fixed_left_keeps_the_reservoir_counts() {
        let reservoir = Turnover {
            fixed_left: true,
            ..spec(0.3)
        };
        let initial = [[10, 10], [20, 0], [20, 0]];
        let result = simulate(
            &reservoir,
            &initial,
            30.0,
            5,
            1_000_000,
            false,
            &Control::default(),
        )
        .expect("reservoir chain");
        assert!(result.events > 0);
        assert_eq!((result.wild[0], result.mutant[0]), (10, 10));
    }

    #[test]
    fn zero_horizon_and_zero_rates_do_nothing() {
        let initial = [[3, 1], [2, 2]];
        let at_zero = simulate(&spec(0.3), &initial, 0.0, 1, 10, false, &Control::default())
            .expect("zero horizon");
        assert_eq!(at_zero.events, 0);
        assert_eq!((at_zero.wild, at_zero.mutant), (vec![3, 2], vec![1, 2]));
        let frozen = Turnover {
            hop: 0.0,
            ..hops_only()
        };
        let still = simulate(&frozen, &initial, 100.0, 1, 10, false, &Control::default())
            .expect("no hazards");
        assert_eq!(still.events, 0);
    }

    #[test]
    fn the_same_seed_gives_the_same_path() {
        let initial = [[20, 0], [15, 5], [20, 0]];
        let run = |seed| {
            let r = simulate(
                &spec(0.2),
                &initial,
                40.0,
                seed,
                1_000_000,
                false,
                &Control::default(),
            )
            .expect("chain");
            (r.wild, r.mutant, r.events)
        };
        assert_eq!(run(3), run(3));
        assert_ne!(run(3), run(4));
    }

    #[test]
    fn invalid_setups_are_numerical_errors() {
        let initial = [[5, 5]];
        let control = Control::default();
        let code = |spec: &Turnover, counts: &[[u64; 2]], end: f64| {
            simulate(spec, counts, end, 1, 10, false, &control)
                .expect_err("invalid setup")
                .code
        };
        let too_likely = Turnover {
            mutation: 1.5,
            ..spec(0.1)
        };
        let not_a_number = Turnover {
            death: f64::NAN,
            ..spec(0.1)
        };
        assert_eq!(code(&too_likely, &initial, 1.0), "CRN-R030");
        assert_eq!(code(&not_a_number, &initial, 1.0), "CRN-R030");
        assert_eq!(code(&spec(0.1), &[], 1.0), "CRN-R030");
        assert_eq!(code(&spec(0.1), &vec![[1, 0]; 4097], 1.0), "CRN-R030");
        assert_eq!(code(&spec(0.1), &initial, f64::INFINITY), "CRN-R030");
    }

    #[test]
    fn budget_and_cancellation_have_their_own_codes() {
        let initial = [[20, 0], [20, 0]];
        let budget = simulate(
            &spec(0.2),
            &initial,
            1.0e6,
            1,
            10,
            false,
            &Control::default(),
        )
        .expect_err("event budget");
        assert_eq!(budget.code, "CRN-R040");
        let control = Control::default();
        control.cancel();
        let cancelled =
            simulate(&spec(0.2), &initial, 1.0, 1, 1_000, false, &control).expect_err("cancelled");
        assert_eq!(cancelled.code, "CRN-R044");
    }
}
