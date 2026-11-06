//! The `hybrid` method: jump times by integrated propensities, continuous
//! variables by Dormand–Prince 5(4) with dense output. See
//! [`crate::continuous`] for the method and its references.
use super::*;
use crate::continuous::{Attempt, Dense, OnBound, Stepper, Variable};
use crate::flat::{Network, RateCache};

/// Where a step ends.
#[derive(Clone, Copy, PartialEq)]
enum Stop {
    /// The integrated propensity reaches the threshold: an event.
    Hazard,
    /// Variable `k` reaches its lower bound.
    Bound(usize),
}

/// The jump rates and the right-hand side between two events.
struct Flow<'a> {
    net: &'a Network,
    vars: &'a [Variable],
    params: &'a [f64],
    /// Channels whose rates read a continuous variable.
    dynamic: Vec<bool>,
    /// Rates of the other channels, constant between events.
    fixed: Vec<f64>,
    fixed_total: f64,
    cache: RateCache,
    /// Observable values with the continuous columns of the current stage.
    scratch: Vec<f64>,
    rates: Vec<f64>,
}

fn derivative_error(k: usize, value: f64) -> Halt {
    numerical(
        &format!("/continuous/{k}/derivative"),
        format!("derivative is {value}"),
        "Derivatives must be finite; rescale the variables or tighten the tolerances.",
    )
}

/// The value at which rates, derivatives and reports read variable `k`: a
/// clamped variable never reads below its bound.
fn read(vars: &[Variable], k: usize, x: f64) -> f64 {
    match vars[k].bound {
        Some((lower, OnBound::Clamp)) => x.max(lower),
        _ => x,
    }
}

/// Writes the variables `y` into their observable columns.
fn store(vars: &[Variable], y: &[f64], values: &mut [f64]) {
    for (k, v) in vars.iter().enumerate() {
        values[v.index] = read(vars, k, y[k]);
    }
}

impl<'a> Flow<'a> {
    fn new(net: &'a Network, model: &'a crate::model::Compiled) -> Self {
        let vars = model.continuous.as_slice();
        let dynamic = net
            .channels
            .iter()
            .map(|c| {
                let mut flags = vec![false; model.observables.len()];
                model.rules[c.rule].rate.reads(&mut flags);
                vars.iter().any(|v| flags[v.index])
            })
            .collect();
        Self {
            net,
            vars,
            params: &model.params,
            dynamic,
            fixed: vec![0.0; net.channels.len()],
            fixed_total: 0.0,
            cache: net.rate_cache(),
            scratch: Vec::new(),
            rates: vec![0.0; net.channels.len()],
        }
    }

    fn load(&mut self, values: &[f64], y: &[f64]) {
        self.scratch.clear();
        self.scratch.extend_from_slice(values);
        store(self.vars, y, &mut self.scratch);
    }

    /// Rates of the channels that do not read continuous variables.
    fn fix(&mut self, counts: &[u64], values: &[f64]) -> Result<(), Halt> {
        let env = Env {
            params: self.params,
            observables: values,
            attributes: [&[], &[]],
        };
        self.cache.reset();
        self.fixed_total = 0.0;
        for (i, rate) in self.fixed.iter_mut().enumerate() {
            *rate = if self.dynamic[i] {
                0.0
            } else {
                self.cache.channel(self.net, i, counts, &env)?
            };
            self.fixed_total += *rate;
        }
        Ok(())
    }

    /// dy/dt: the variables' derivatives, then the total propensity.
    fn rhs(
        &mut self,
        counts: &[u64],
        values: &[f64],
        y: &[f64],
        dy: &mut [f64],
    ) -> Result<(), Halt> {
        self.load(values, y);
        let env = Env {
            params: self.params,
            observables: &self.scratch,
            attributes: [&[], &[]],
        };
        self.cache.reset();
        let mut total = self.fixed_total;
        for (i, &dynamic) in self.dynamic.iter().enumerate() {
            if dynamic {
                total += self.cache.channel(self.net, i, counts, &env)?;
            }
        }
        if !total.is_finite() {
            return Err(numerical(
                "/rules",
                "total propensity overflows",
                "Reduce rates or populations.",
            ));
        }
        let m = self.vars.len();
        dy[m] = total;
        for (k, v) in self.vars.iter().enumerate() {
            let mut g = v.derivative.eval(&env);
            if !g.is_finite() {
                return Err(derivative_error(k, g));
            }
            if let Some((lower, OnBound::Clamp)) = v.bound {
                if y[k] <= lower && g < 0.0 {
                    g = 0.0;
                }
            }
            dy[k] = g;
        }
        Ok(())
    }

    /// Every channel's rate at the variables `y`, and their sum.
    fn all_rates(&mut self, counts: &[u64], values: &[f64], y: &[f64]) -> Result<f64, Halt> {
        self.load(values, y);
        let env = Env {
            params: self.params,
            observables: &self.scratch,
            attributes: [&[], &[]],
        };
        self.cache.reset();
        let mut total = 0.0;
        for i in 0..self.rates.len() {
            self.rates[i] = if self.dynamic[i] {
                self.cache.channel(self.net, i, counts, &env)?
            } else {
                self.fixed[i]
            };
            total += self.rates[i];
        }
        Ok(total)
    }

    /// Whether a variable with `on_bound: "stop"` sits at its bound with a
    /// derivative that is not positive; such a variable is set to the bound.
    fn depleted(&self, values: &mut [f64]) -> Result<bool, Halt> {
        let env = Env {
            params: self.params,
            observables: values,
            attributes: [&[], &[]],
        };
        let mut stopped = None;
        for (k, v) in self.vars.iter().enumerate() {
            if let Some((lower, OnBound::Stop)) = v.bound {
                if values[v.index] <= lower {
                    let g = v.derivative.eval(&env);
                    if !g.is_finite() {
                        return Err(derivative_error(k, g));
                    }
                    if g <= 0.0 {
                        stopped = Some((v.index, lower));
                        break;
                    }
                }
            }
        }
        if let Some((index, lower)) = stopped {
            values[index] = lower;
        }
        Ok(stopped.is_some())
    }
}

fn select<R: RngCore>(rates: &[f64], total: f64, words: &mut Words<'_, R>) -> Result<usize, Halt> {
    if !(total.is_finite() && total > 0.0) {
        return Err(numerical(
            "/method",
            format!(
                "the integrated propensity reached its threshold where the total propensity is {total}"
            ),
            "Tighten rtol and atol.",
        ));
    }
    let floor = total * RESOLUTION;
    if rates.iter().any(|&a| a > 0.0 && a < floor) {
        return Err(Halt::Numerical(Diagnostic::error(
            "CRN-R031",
            "/rules",
            "a positive channel is below 53-bit resolution of the total",
            "Rescale rates.",
        )));
    }
    let mut target = words.uniform()? * total;
    let mut last = 0;
    for (i, &a) in rates.iter().enumerate().filter(|&(_, a)| *a > 0.0) {
        last = i;
        if target < a {
            return Ok(i);
        }
        target -= a;
    }
    Ok(last)
}

impl Path<'_> {
    /// Reports before `until` from the dense output of the current step.
    fn save_dense(
        &self,
        out: &mut PathOutput,
        values: &[f64],
        vars: &[Variable],
        dense: &Dense,
        until: f64,
    ) {
        let n = values.len();
        while out.saved_points < self.report_times.len()
            && self.report_times[out.saved_points] < until
        {
            let row = out.saved_points * n;
            out.values[row..row + n].copy_from_slice(values);
            let theta = ((self.report_times[out.saved_points] - dense.t) / dense.h).clamp(0.0, 1.0);
            for (k, v) in vars.iter().enumerate() {
                out.values[row + v.index] = read(vars, k, dense.component(k, theta));
            }
            out.saved_points += 1;
        }
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn simulate_hybrid<R: RngCore, O: EventObserver>(
        &self,
        rtol: f64,
        atol: f64,
        net: &Network,
        counts: &mut [u64],
        values: &mut [f64],
        words: &mut Words<'_, R>,
        out: &mut PathOutput,
        time: &mut f64,
        monitor: &mut Option<Monitor<'_>>,
        observer: &O,
    ) -> Result<Option<f64>, Halt> {
        let mut flow = Flow::new(net, self.model);
        let vars = flow.vars;
        let m = vars.len();
        let mut y = vec![0.0; m + 1];
        let mut stepper = Stepper::new(m + 1, rtol, atol);
        let mut dense = Dense::default();
        let mut snapshot = Vec::with_capacity(values.len());
        let mut live = out.peak_instances;
        let mut steps: u64 = 0;
        let mut predicted: Option<f64> = None;
        loop {
            // A new interval between events: fixed counts, a fresh threshold.
            self.poll()?;
            if flow.depleted(values)? {
                return Ok(Some(*time));
            }
            flow.fix(counts, values)?;
            for (k, v) in vars.iter().enumerate() {
                y[k] = values[v.index];
            }
            y[m] = 0.0;
            let threshold = -libm::log(words.open()?);
            stepper.fresh = false;
            snapshot.clear();
            snapshot.extend_from_slice(values);
            let (stop, until) = {
                let fixed_counts: &[u64] = counts;
                let snapshot = snapshot.as_slice();
                let mut rhs = |_t: f64, state: &[f64], dy: &mut [f64]| {
                    flow.rhs(fixed_counts, snapshot, state, dy)
                };
                loop {
                    if steps >= self.limits.max_words {
                        return Err(Halt::Limit(
                            Termination::SamplingLimit,
                            Diagnostic::error(
                                "CRN-R041",
                                "/budget/max_sampling_steps_per_replicate",
                                "integration step limit reached",
                                "Increase the sampling-step budget or loosen the tolerances.",
                            ),
                        ));
                    }
                    steps += 1;
                    if steps % POLL_EVENTS == 0 {
                        self.poll()?;
                    }
                    let remaining = self.end - *time;
                    if remaining <= 0.0 {
                        // Reports at the horizon read `values`.
                        store(vars, &y, values);
                        return Ok(None);
                    }
                    let mut h = match predicted {
                        Some(h) => h,
                        None => stepper.initial_step(*time, &y, remaining, &mut rhs)?,
                    };
                    let rate = stepper.derivative(*time, &y, &mut rhs)?[m];
                    if rate > 0.0 {
                        // The threshold is expected within a quarter of the step.
                        h = h.min(4.0 * (threshold - y[m]) / rate);
                    }
                    let last = h >= remaining;
                    if last {
                        h = remaining;
                    }
                    if !(h.is_finite() && h > 4.0 * f64::EPSILON * time.abs()) {
                        return Err(numerical(
                            "/method",
                            format!("integration step {h} underflows at time {time}"),
                            "Loosen rtol and atol, or rescale the variables and rates.",
                        ));
                    }
                    match stepper.attempt(*time, &y, h, &mut rhs, &mut dense)? {
                        Attempt::Rejected(smaller) => {
                            predicted = Some(smaller);
                            continue;
                        }
                        Attempt::Accepted { next_h } => predicted = Some(next_h),
                    }
                    // The earliest event inside the step; a bound wins a tie.
                    let mut theta = 1.0;
                    let mut stop = None;
                    if stepper.next[m] >= threshold {
                        theta = dense.crossing_up(m, threshold);
                        stop = Some(Stop::Hazard);
                    }
                    for (k, v) in vars.iter().enumerate() {
                        let Some((lower, action)) = v.bound else {
                            continue;
                        };
                        let crossed = match action {
                            OnBound::Stop => stepper.next[k] <= lower,
                            OnBound::Clamp => stepper.next[k] < lower,
                        };
                        if y[k] > lower && crossed {
                            let at = dense.crossing_down(k, lower);
                            if at <= theta {
                                theta = at;
                                stop = Some(Stop::Bound(k));
                            }
                        }
                    }
                    let until = if stop.is_none() && last {
                        self.end
                    } else {
                        (*time + theta * h).min(self.end)
                    };
                    self.save_dense(out, values, vars, &dense, until);
                    match stop {
                        None => {
                            y.copy_from_slice(&stepper.next);
                            *time = until;
                        }
                        Some(Stop::Bound(k))
                            if matches!(vars[k].bound, Some((_, OnBound::Clamp))) =>
                        {
                            // Continue from the bound with the same threshold.
                            dense.state(theta, &mut y);
                            y[k] = vars[k].bound.expect("bounded variable").0;
                            *time = until;
                            stepper.fresh = false;
                        }
                        Some(stop) => {
                            dense.state(theta, &mut y);
                            break (stop, until);
                        }
                    }
                }
            };
            match stop {
                Stop::Bound(k) => {
                    y[k] = vars[k].bound.expect("bounded variable").0;
                    *time = until;
                    store(vars, &y, values);
                    return Ok(Some(until));
                }
                Stop::Hazard => {
                    if until <= *time {
                        return Err(numerical(
                            "/time",
                            "clock failed to advance",
                            "Rescale rates and time.",
                        ));
                    }
                    *time = until;
                    store(vars, &y, values);
                    let total = flow.all_rates(counts, values, &y)?;
                    let channel = select(&flow.rates, total, words)?;
                    self.apply_flat(net, channel, counts, &mut live, out)?;
                    net.refresh_channel(channel, counts, values, true)?;
                    if self.observe_flat(net, channel, values, until, monitor, observer)
                        || self.stop_values(values)?
                    {
                        return Ok(Some(until));
                    }
                }
            }
        }
    }
}
