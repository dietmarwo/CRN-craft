//! CGP adaptive Poisson leaps with one competing critical event.
//! Cao, Gillespie & Petzold (2006), DOI 10.1063/1.2159468: the species-based
//! bound and, for expression rates, the propensity condition of their
//! section III with first differences (Gillespie & Petzold 2003,
//! DOI 10.1063/1.1613254).
use super::*;
use crate::rules::flat::{Network, RateCache, check_instances, direct};
use crate::rules::model::MAX_COUNT;

struct LeapWorkspace {
    cache: RateCache,
    rates: Vec<f64>,
    critical: Vec<f64>,
    proposed: Vec<i128>,
    changed: Vec<bool>,
    delta: Vec<(f64, f64)>,
    firings: Vec<u64>,
    remaining: Vec<u64>,
    remaining_probability: Vec<f64>,
    /// Per channel: mean and variance of the rate change per unit time, and
    /// the smallest nonzero change from one firing.
    change: Vec<(f64, f64, f64)>,
    shifted_counts: Vec<u64>,
    shifted_values: Vec<f64>,
    /// Rate cache for the shifted state of one firing; reset per firing.
    shifted_cache: RateCache,
    row: Vec<f64>,
}

impl LeapWorkspace {
    fn new(net: &Network) -> Self {
        let m = net.channels.len();
        let n = net.species.len();
        Self {
            cache: net.rate_cache(),
            rates: vec![0.0; m],
            critical: vec![0.0; m],
            proposed: vec![0; n],
            changed: vec![false; n],
            delta: vec![(0.0, 0.0); n],
            firings: vec![0; m],
            remaining: vec![0; n],
            remaining_probability: vec![1.0; n],
            change: vec![(0.0, 0.0, f64::INFINITY); m],
            shifted_counts: vec![0; n],
            shifted_values: Vec::new(),
            shifted_cache: net.rate_cache(),
            row: Vec::new(),
        }
    }

    fn rates(&mut self, net: &Network, counts: &[u64], env: &Env<'_>) -> Result<(), Halt> {
        self.cache.reset();
        for (i, a) in self.rates.iter_mut().enumerate() {
            *a = self.cache.channel(net, i, counts, env)?;
        }
        Ok(())
    }

    fn begin(&mut self, counts: &[u64]) {
        for (v, &n) in self.proposed.iter_mut().zip(counts) {
            *v = i128::from(n);
        }
        self.firings.fill(0);
    }

    fn fire(&mut self, net: &Network, channel: usize, count: u64) {
        self.firings[channel] = count;
        for &(s, d) in &net.channels[channel].changes {
            self.proposed[s] += i128::from(d) * i128::from(count);
        }
    }

    fn valid(&self) -> bool {
        self.proposed
            .iter()
            .all(|&x| (0..=i128::from(MAX_COUNT)).contains(&x))
    }

    /// Rate changes a_j(x + nu_k) - a_j(x) of every channel j after one
    /// firing of channel k, written to `row`. Rates may read global
    /// observables, so these are first differences of the rates.
    fn differences(
        &mut self,
        net: &Network,
        counts: &[u64],
        params: &[f64],
        k: usize,
        row: &mut [f64],
    ) -> Result<(), Halt> {
        self.shifted_counts.copy_from_slice(counts);
        for &(s, d) in &net.channels[k].changes {
            // An enabled channel has its reactants; products only grow.
            self.shifted_counts[s] = self.shifted_counts[s].saturating_add_signed(d);
        }
        net.refresh(&self.shifted_counts, &mut self.shifted_values)?;
        let env = Env {
            params,
            observables: &self.shifted_values,
            attributes: [&[], &[]],
        };
        self.shifted_cache.reset();
        for (j, f) in row.iter_mut().enumerate() {
            *f = self
                .shifted_cache
                .channel(net, j, &self.shifted_counts, &env)?
                - self.rates[j];
        }
        Ok(())
    }

    /// CGP criticality for rates: a channel is critical if ten firings of an
    /// enabled channel can exhaust its rate, as a reactant count within ten
    /// firings of zero makes a mass-action channel critical. Then the
    /// largest step for which every channel's expected rate change and its
    /// standard deviation stay within epsilon * a_j. Zero means that one
    /// firing already exceeds that allowance.
    fn propensity_step(
        &mut self,
        net: &Network,
        counts: &[u64],
        values: &[f64],
        params: &[f64],
        epsilon: f64,
    ) -> Result<f64, Halt> {
        let m = net.channels.len();
        self.shifted_values.resize(values.len(), 0.0);
        let mut row = std::mem::take(&mut self.row);
        row.resize(m, 0.0);
        for k in 0..m {
            if self.rates[k] == 0.0 {
                continue;
            }
            self.differences(net, counts, params, k, &mut row)?;
            for (j, &f) in row.iter().enumerate() {
                if f < 0.0 && self.rates[j] <= -10.0 * f {
                    self.critical[j] = self.rates[j];
                }
            }
        }
        self.change.fill((0.0, 0.0, f64::INFINITY));
        for k in 0..m {
            let a_k = self.rates[k];
            if a_k == 0.0 || self.critical[k] != 0.0 {
                continue;
            }
            self.differences(net, counts, params, k, &mut row)?;
            for (j, &f) in row.iter().enumerate() {
                if f != 0.0 {
                    let (mean, variance, floor) = &mut self.change[j];
                    *mean += f * a_k;
                    *variance += f * f * a_k;
                    *floor = floor.min(f.abs());
                }
            }
        }
        self.row = row;
        let mut dt = f64::INFINITY;
        for (j, &(mean, variance, floor)) in self.change.iter().enumerate() {
            if mean == 0.0 && variance == 0.0 {
                continue;
            }
            // As for a submolecule species allowance: if one firing moves a
            // rate by more than epsilon of itself, a frozen-rate leap is
            // unsuitable, and the caller takes an exact step. The CGP floor
            // c_j would accept such a change; for a vacancy rate with a few
            // free places it biases the population (review attempt 1).
            let allowance = epsilon * self.rates[j];
            if allowance < floor {
                return Ok(0.0);
            }
            if mean != 0.0 {
                dt = dt.min(allowance / mean.abs());
            }
            if variance > 0.0 {
                dt = dt.min(allowance * allowance / variance);
            }
        }
        Ok(dt)
    }

    #[allow(clippy::cast_precision_loss, clippy::too_many_arguments)]
    fn tau<R: RngCore>(
        &mut self,
        net: &Network,
        counts: &[u64],
        values: &[f64],
        params: &[f64],
        time: f64,
        target: f64,
        epsilon: f64,
        words: &mut Words<'_, R>,
    ) -> Result<f64, Halt> {
        self.delta.fill((0.0, 0.0));
        for (i, c) in net.channels.iter().enumerate() {
            let a = self.rates[i];
            // Low-count births are critical too, so a new mutant is visible
            // even if its parent population is large.
            let critical = c
                .changes
                .iter()
                .any(|&(s, d)| counts[s] <= 10 * d.unsigned_abs())
                || c.reactants.iter().any(|&s| counts[s] <= 10);
            self.critical[i] = if critical { a } else { 0.0 };
            if !critical {
                for &(s, d) in &c.changes {
                    let d = d as f64;
                    self.delta[s].0 += d * a;
                    self.delta[s].1 += d * d * a;
                }
            }
        }
        let total: f64 = self.rates.iter().sum();
        let mut dt = target - time;
        let noncritical = self
            .rates
            .iter()
            .zip(&self.critical)
            .filter_map(|(&a, &c)| (c == 0.0).then_some(a))
            .fold(0.0, f64::max);
        if noncritical == 0.0 {
            return self.exact(net, counts, time, target, words);
        }
        dt = dt.min(1.0e8 / noncritical);
        for (&n, &(drift, variance)) in counts.iter().zip(&self.delta) {
            // Expression rates use g=2 as the documented heuristic.
            let scale = epsilon * n as f64 / 2.0;
            // The CGP unit floor is a mass-action heuristic. For arbitrary
            // expression rates, a submolecule allowance means a leap is
            // unsuitable even if the reactant is above the critical cutoff.
            if (drift != 0.0 || variance > 0.0) && scale < 1.0 {
                return self.exact(net, counts, time, target, words);
            }
            if drift != 0.0 {
                dt = dt.min(scale / drift.abs());
            }
            if variance > 0.0 {
                dt = dt.min(scale * scale / variance);
            }
        }
        if dt < 1.0 / total {
            return self.exact(net, counts, time, target, words);
        }
        dt = dt.min(self.propensity_step(net, counts, values, params, epsilon)?);
        if dt < 1.0 / total {
            return self.exact(net, counts, time, target, words);
        }
        let critical = direct(&self.critical, time, words)?;
        for _ in 0..20 {
            let proposed = (time + dt).min(target);
            let end = critical.map_or(proposed, |(at, _)| at.min(proposed));
            if end <= time {
                break;
            }
            self.begin(counts);
            for i in 0..self.rates.len() {
                if self.critical[i] == 0.0 {
                    let k =
                        crate::tau::poisson_with(self.rates[i] * (end - time), || words.open())?;
                    self.fire(net, i, k);
                }
            }
            if let Some((at, i)) = critical {
                if at == end {
                    self.fire(net, i, 1);
                }
            }
            if self.valid() {
                return Ok(end);
            }
            dt *= 0.5;
        }
        self.exact(net, counts, time, target, words)
    }

    fn exact<R: RngCore>(
        &mut self,
        net: &Network,
        counts: &[u64],
        time: f64,
        target: f64,
        words: &mut Words<'_, R>,
    ) -> Result<f64, Halt> {
        self.begin(counts);
        if let Some((next, i)) = direct(&self.rates, time, words)? {
            if next <= target {
                self.fire(net, i, 1);
                return Ok(next);
            }
        }
        Ok(target)
    }

    #[allow(clippy::cast_precision_loss)]
    fn binomial<R: RngCore>(
        &mut self,
        net: &Network,
        counts: &[u64],
        values: &[f64],
        dt: f64,
        params: &[f64],
        words: &mut Words<'_, R>,
    ) -> Result<(), Halt> {
        self.begin(counts);
        self.remaining.copy_from_slice(counts);
        self.remaining_probability.fill(1.0);
        self.split_phase(net, counts, dt, -1, words)?;
        let after: Vec<_> = self.proposed.iter().map(|&x| x as u64).collect();
        let mut updated = values.to_vec();
        for ((changed, old), new) in self.changed.iter_mut().zip(counts).zip(&after) {
            *changed = old != new;
        }
        net.refresh_changed(&self.changed, &after, &mut updated)?;
        self.rates(
            net,
            &after,
            &Env {
                params,
                observables: &updated,
                attributes: [&[], &[]],
            },
        )?;
        self.remaining.copy_from_slice(&after);
        self.remaining_probability.fill(1.0);
        self.split_phase(net, &after, dt, 1, words)?;
        if !self.valid() {
            return Err(numerical(
                "/method",
                "binomial count exceeds range",
                "Reduce step or population.",
            ));
        }
        Ok(())
    }

    #[allow(clippy::cast_precision_loss)]
    fn split_phase<R: RngCore>(
        &mut self,
        net: &Network,
        counts: &[u64],
        dt: f64,
        direction: i64,
        words: &mut Words<'_, R>,
    ) -> Result<(), Halt> {
        // Group a parent's competing channels into one multinomial. Cap the
        // total probability at one, preserving the relative channel weights.
        self.delta.fill((0.0, 0.0));
        for (i, c) in net
            .channels
            .iter()
            .enumerate()
            .filter(|(_, c)| c.changes[0].1 == direction)
        {
            self.delta[c.reactants[0]].0 += self.rates[i];
        }
        for (i, c) in net
            .channels
            .iter()
            .enumerate()
            .filter(|(_, c)| c.changes[0].1 == direction)
        {
            let s = c.reactants[0];
            let a = self.rates[i];
            if counts[s] == 0 || a == 0.0 || self.remaining[s] == 0 {
                continue;
            }
            let total = self.delta[s].0;
            if !total.is_finite() {
                return Err(numerical(
                    "/method",
                    "binomial rate sum overflows",
                    "Rescale rates.",
                ));
            }
            // Competing risks with rates frozen over the step: the parent's
            // first event falls in the step with probability 1 - e^(-L dt),
            // L its total per-parent rate, and is channel i with a / total.
            let probability = -libm::expm1(-total * dt / counts[s] as f64) * (a / total);
            let available = self.remaining_probability[s];
            let conditional = if available > 0.0 {
                (probability / available).min(1.0)
            } else {
                0.0
            };
            let k = binomial(self.remaining[s], conditional, words)?;
            self.remaining[s] -= k;
            self.remaining_probability[s] = (available - probability).max(0.0);
            self.fire(net, i, k);
        }
        Ok(())
    }
}

#[allow(clippy::cast_precision_loss)]
fn binomial<R: RngCore>(n: u64, p: f64, words: &mut Words<'_, R>) -> Result<u64, Halt> {
    if p == 0.0 {
        return Ok(0);
    }
    if p == 1.0 {
        return Ok(n);
    }
    let invert = p > 0.5;
    let log_failure = libm::log1p(-if invert { 1.0 - p } else { p });
    let mut left = n;
    let mut successes = 0;
    while left > 0 {
        let gap = (libm::log(words.open()?) / log_failure).floor();
        if gap >= left as f64 {
            break;
        }
        left -= gap as u64 + 1;
        successes += 1;
    }
    Ok(if invert { n - successes } else { successes })
}

impl Path<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn simulate_leaps<R: RngCore, O: EventObserver>(
        &self,
        net: &Network,
        counts: &mut [u64],
        values: &mut [f64],
        words: &mut Words<'_, R>,
        out: &mut PathOutput,
        time: &mut f64,
        monitor: &mut Option<Monitor<'_>>,
        observer: &O,
    ) -> Result<Option<f64>, Halt> {
        let mut work = LeapWorkspace::new(net);
        let mut steps = 0;
        let mut lattice: u64 = 0;
        loop {
            self.poll()?;
            while out.saved_points < self.report_times.len()
                && self.report_times[out.saved_points] <= *time
            {
                let n = values.len();
                let offset = out.saved_points * n;
                out.values[offset..offset + n].copy_from_slice(values);
                out.saved_points += 1;
            }
            if *time >= self.end {
                return Ok(None);
            }
            // A zero-firing leap can advance without drawing a random
            // word. Bound these iterations as well as the random draws.
            if steps >= self.limits.max_words {
                return Err(Halt::Limit(
                    Termination::SamplingLimit,
                    Diagnostic::error(
                        "CRN-R041",
                        "/budget/max_sampling_steps_per_replicate",
                        "leap step limit reached",
                        "Increase the sampling-step budget or step size.",
                    ),
                ));
            }
            steps += 1;
            let target = self
                .report_times
                .get(out.saved_points)
                .copied()
                .unwrap_or(self.end);
            work.rates(
                net,
                counts,
                &Env {
                    params: &self.model.params,
                    observables: values,
                    attributes: [&[], &[]],
                },
            )?;
            let next = match self.method {
                crate::rules::Method::Tau { epsilon } => work.tau(
                    net,
                    counts,
                    values,
                    &self.model.params,
                    *time,
                    target,
                    epsilon,
                    words,
                )?,
                crate::rules::Method::Binomial { step } => {
                    // Steps lie on the lattice k * step, independent of the
                    // report grid; a report between lattice points holds the
                    // state of the last step.
                    lattice += 1;
                    #[allow(clippy::cast_precision_loss)]
                    let next = (lattice as f64 * step).min(self.end);
                    let n = values.len();
                    while out.saved_points < self.report_times.len()
                        && self.report_times[out.saved_points] < next
                    {
                        let offset = out.saved_points * n;
                        out.values[offset..offset + n].copy_from_slice(values);
                        out.saved_points += 1;
                    }
                    work.binomial(net, counts, values, next - *time, &self.model.params, words)?;
                    next
                }
                _ => unreachable!(),
            };
            if next <= *time || !work.valid() {
                return Err(numerical(
                    "/method",
                    "invalid leap",
                    "Reduce the step size.",
                ));
            }
            let events = work
                .firings
                .iter()
                .try_fold(0u64, |acc, &k| acc.checked_add(k))
                .ok_or_else(|| {
                    numerical("/method", "leap event count overflows", "Reduce step size.")
                })?;
            self.check_events(out, events)?;
            let live = work.proposed.iter().filter(|&&n| n > 0).count() as u64;
            check_instances(live, self.limits)?;
            out.created_instances += counts
                .iter()
                .zip(&work.proposed)
                .filter(|&(&old, &n)| old == 0 && n > 0)
                .count() as u64;
            for ((n, &v), changed) in counts.iter_mut().zip(&work.proposed).zip(&mut work.changed) {
                *changed = *n != v as u64;
                *n = v as u64;
            }
            *time = next;
            out.peak_instances = out.peak_instances.max(live);
            out.events += events;
            for (c, &k) in net.channels.iter().zip(&work.firings) {
                out.rule_events[c.rule] += k;
            }
            net.refresh_changed(&work.changed, counts, values)?;
            if events > 0 && (O::ENABLED || monitor.is_some()) {
                let firings: Vec<_> = net
                    .channels
                    .iter()
                    .zip(&work.firings)
                    .filter_map(|(c, &k)| (k > 0).then_some((c.rule, k)))
                    .collect();
                let event = Event {
                    replicate: self.id,
                    time: next,
                    rule: None,
                    reactants: &[],
                    products: &[],
                    observables: values,
                    firings: &firings,
                };
                if let Some(m) = monitor {
                    m.observe(&event, &self.model.params);
                }
                if O::ENABLED && observer.on_event(&event) == ObserverAction::Stop {
                    return Ok(Some(next));
                }
            }
            if self.stop_values(values)? {
                return Ok(Some(next));
            }
        }
    }
}
