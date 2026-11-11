//! The `hybrid` method: jump times by integrated propensities, continuous
//! variables by Dormand–Prince 5(4) with dense output. See
//! [`crate::rules::continuous`] for the method and its references.
use super::*;
use crate::rules::continuous::{Attempt, Dense, OnBound, Stepper, Variable, interpolate};
use crate::rules::flat::{Network, RateCache};

/// Where a step ends.
#[derive(Clone, Copy, PartialEq)]
enum Stop {
    /// The integrated propensity reaches the threshold: an event.
    Hazard,
    /// The dependent variable at state position `j` reaches its lower bound.
    Bound(usize),
    /// The independent variable `k` reaches its stop bound at the end of
    /// the shared solution.
    Independent(usize),
}

/// How the shared solution of the independent variables ends.
#[derive(Clone)]
pub(crate) enum Ending {
    /// At the request's horizon.
    Horizon,
    /// Variable `k` (a model index) reaches its stop bound: every replicate
    /// that gets there stops at the same time.
    Bound(usize),
    /// The solve failed (a non-finite derivative, a step that underflows or
    /// the step budget); every replicate that gets there fails with it.
    Failed(Halt),
}

/// The independent continuous variables (rules-model/0.4), integrated once
/// per prepared run over [0, end] by Dormand–Prince 5(4) at the request's
/// tolerances, with the same step control, bound handling and continuous
/// extension as the per-replicate integration. Replicates read the
/// variables from this dense output, so results do not depend on threads.
pub(crate) struct Shared {
    /// Model indices of the independent variables, in declaration order.
    vars: Vec<usize>,
    /// Start time and size of each accepted step.
    starts: Vec<f64>,
    sizes: Vec<f64>,
    /// Five continuous-extension coefficients per variable per step.
    coefficients: Vec<f64>,
    initial: Vec<f64>,
    /// The solution covers [0, until].
    until: f64,
    ending: Ending,
}

impl Shared {
    /// Integrates the independent variables of `model` over [0, end]; at
    /// most `max_steps` attempted steps, as one replicate's budget.
    #[allow(clippy::too_many_lines)]
    pub(crate) fn solve(
        model: &crate::rules::model::Compiled,
        rtol: f64,
        atol: f64,
        end: f64,
        max_steps: u64,
    ) -> Option<Self> {
        let all = model.continuous.as_slice();
        let vars: Vec<usize> = (0..all.len()).filter(|&k| all[k].independent).collect();
        if vars.is_empty() {
            return None;
        }
        let m = vars.len();
        let initial: Vec<f64> = vars.iter().map(|&k| all[k].initial).collect();
        let mut shared = Self {
            vars,
            starts: Vec::new(),
            sizes: Vec::new(),
            coefficients: Vec::new(),
            initial,
            until: 0.0,
            ending: Ending::Horizon,
        };
        let mut scratch = vec![0.0; model.observables.len()];
        let vars = shared.vars.clone();
        let params = model.params.as_slice();
        let mut rhs = |_t: f64, y: &[f64], dy: &mut [f64]| -> Result<(), Halt> {
            for (j, &k) in vars.iter().enumerate() {
                scratch[all[k].index] = read(all, k, y[j]);
            }
            let env = Env {
                params,
                observables: &scratch,
                attributes: [&[], &[]],
            };
            for (j, &k) in vars.iter().enumerate() {
                let mut g = all[k].derivative.eval(&env);
                if !g.is_finite() {
                    return Err(derivative_error(k, g));
                }
                if let Some((lower, OnBound::Clamp)) = all[k].bound {
                    if y[j] <= lower && g < 0.0 {
                        g = 0.0;
                    }
                }
                dy[j] = g;
            }
            Ok(())
        };
        let mut y = shared.initial.clone();
        let mut time = 0.0f64;
        let mut stepper = Stepper::new(m, rtol, atol);
        let mut dense = Dense::default();
        let mut predicted: Option<f64> = None;
        let mut steps: u64 = 0;
        let mut solve = |shared: &mut Self| -> Result<Ending, Halt> {
            // A variable that starts at its stop bound and cannot rise.
            let mut dy = vec![0.0; m];
            rhs(0.0, &y, &mut dy)?;
            for (j, &k) in vars.iter().enumerate() {
                if let Some((lower, OnBound::Stop)) = all[k].bound {
                    if y[j] <= lower && dy[j] <= 0.0 {
                        return Ok(Ending::Bound(k));
                    }
                }
            }
            loop {
                if steps >= max_steps {
                    return Err(Halt::Limit(
                        Termination::SamplingLimit,
                        Diagnostic::error(
                            "CRN-R041",
                            "/budget/max_sampling_steps_per_replicate",
                            "integration step limit reached for the independent variables",
                            "Increase the sampling-step budget or loosen the tolerances.",
                        ),
                    ));
                }
                steps += 1;
                let remaining = end - time;
                if remaining <= 0.0 {
                    return Ok(Ending::Horizon);
                }
                let mut h = match predicted {
                    Some(h) => h,
                    None => stepper.initial_step(time, &y, remaining, &mut rhs)?,
                };
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
                match stepper.attempt(time, &y, h, &mut rhs, &mut dense)? {
                    Attempt::Rejected(smaller) => {
                        predicted = Some(smaller);
                        continue;
                    }
                    Attempt::Accepted { next_h } => predicted = Some(next_h),
                }
                let mut theta = 1.0;
                let mut stop = None;
                for (j, &k) in vars.iter().enumerate() {
                    let Some((lower, action)) = all[k].bound else {
                        continue;
                    };
                    let crossed = match action {
                        OnBound::Stop => stepper.next[j] <= lower,
                        OnBound::Clamp => stepper.next[j] < lower,
                    };
                    if y[j] > lower && crossed {
                        let at = dense.crossing_down(j, lower);
                        if at < theta || (at == theta && action == OnBound::Stop) {
                            theta = at;
                            stop = Some(j);
                        }
                    }
                }
                let until = if stop.is_none() && last {
                    end
                } else {
                    (time + theta * h).min(end)
                };
                shared.starts.push(dense.t);
                shared.sizes.push(dense.h);
                for j in 0..m {
                    shared.coefficients.extend(dense.coefficients(j));
                }
                time = until;
                match stop {
                    None => y.copy_from_slice(&stepper.next),
                    Some(j) => {
                        let k = vars[j];
                        let (lower, action) = all[k].bound.expect("bounded variable");
                        if action == OnBound::Stop {
                            return Ok(Ending::Bound(k));
                        }
                        // Continue from the bound.
                        dense.state(theta, &mut y);
                        y[j] = lower;
                        stepper.fresh = false;
                    }
                }
            }
        };
        let ending = solve(&mut shared).unwrap_or_else(Ending::Failed);
        shared.until = time;
        shared.ending = ending;
        Some(shared)
    }

    /// Calls `f(j, x)` for every independent variable `j` (in declaration
    /// order) with its value `x` (no clamp applied) at time `t` in
    /// [0, until]: from the last step that starts at or before `t`.
    fn each(&self, t: f64, mut f: impl FnMut(usize, f64)) {
        if self.starts.is_empty() {
            for (j, &x) in self.initial.iter().enumerate() {
                f(j, x);
            }
            return;
        }
        let i = self.starts.partition_point(|&s| s <= t).saturating_sub(1);
        let theta = ((t - self.starts[i]) / self.sizes[i]).clamp(0.0, 1.0);
        let m = self.vars.len();
        for (j, r) in self.coefficients[i * m * 5..(i + 1) * m * 5]
            .chunks_exact(5)
            .enumerate()
        {
            f(j, interpolate([r[0], r[1], r[2], r[3], r[4]], theta));
        }
    }
}

/// The jump rates and the right-hand side between two events.
struct Flow<'a> {
    net: &'a Network,
    vars: &'a [Variable],
    params: &'a [f64],
    /// Model indices of the variables integrated per replicate, in
    /// declaration order: all of them, or the dependent ones when the
    /// independent ones come from `shared`. The state vector holds them,
    /// then the integrated propensity.
    dependent: Vec<usize>,
    shared: Option<&'a Shared>,
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

/// Writes the variables at time `t` into their observable columns: the
/// dependent ones from the state `y`, the independent ones from the shared
/// solution.
fn store(
    vars: &[Variable],
    dependent: &[usize],
    shared: Option<&Shared>,
    y: &[f64],
    t: f64,
    values: &mut [f64],
) {
    for (j, &k) in dependent.iter().enumerate() {
        values[vars[k].index] = read(vars, k, y[j]);
    }
    if let Some(shared) = shared {
        shared.each(t, |j, x| {
            let k = shared.vars[j];
            values[vars[k].index] = read(vars, k, x);
        });
    }
}

impl<'a> Flow<'a> {
    fn new(
        net: &'a Network,
        model: &'a crate::rules::model::Compiled,
        shared: Option<&'a Shared>,
    ) -> Self {
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
            dependent: (0..vars.len())
                .filter(|&k| !(shared.is_some() && vars[k].independent))
                .collect(),
            shared,
            dynamic,
            fixed: vec![0.0; net.channels.len()],
            fixed_total: 0.0,
            cache: net.rate_cache(),
            scratch: Vec::new(),
            rates: vec![0.0; net.channels.len()],
        }
    }

    fn load(&mut self, values: &[f64], y: &[f64], t: f64) {
        self.scratch.clear();
        self.scratch.extend_from_slice(values);
        store(
            self.vars,
            &self.dependent,
            self.shared,
            y,
            t,
            &mut self.scratch,
        );
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

    /// dy/dt at time `t`: the dependent variables' derivatives, then the
    /// total propensity.
    fn rhs(
        &mut self,
        counts: &[u64],
        values: &[f64],
        t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) -> Result<(), Halt> {
        self.load(values, y, t);
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
        let m = self.dependent.len();
        dy[m] = total;
        for (j, &k) in self.dependent.iter().enumerate() {
            let v = &self.vars[k];
            let mut g = v.derivative.eval(&env);
            if !g.is_finite() {
                return Err(derivative_error(k, g));
            }
            if let Some((lower, OnBound::Clamp)) = v.bound {
                if y[j] <= lower && g < 0.0 {
                    g = 0.0;
                }
            }
            dy[j] = g;
        }
        Ok(())
    }

    /// Every channel's rate at time `t` and the variables `y`, and their
    /// sum.
    fn all_rates(
        &mut self,
        counts: &[u64],
        values: &[f64],
        y: &[f64],
        t: f64,
    ) -> Result<f64, Halt> {
        self.load(values, y, t);
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

    /// Whether a dependent variable with `on_bound: "stop"` sits at its
    /// bound with a derivative that is not positive; such a variable is set
    /// to the bound. Independent variables end the shared solution instead.
    fn depleted(&self, values: &mut [f64]) -> Result<bool, Halt> {
        let env = Env {
            params: self.params,
            observables: values,
            attributes: [&[], &[]],
        };
        let mut stopped = None;
        for &k in &self.dependent {
            let v = &self.vars[k];
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
    /// Reports before `until` from the dense output of the current step,
    /// and the independent variables from the shared solution.
    #[allow(clippy::too_many_arguments)]
    fn save_dense(
        &self,
        out: &mut PathOutput,
        values: &[f64],
        vars: &[Variable],
        dependent: &[usize],
        shared: Option<&Shared>,
        dense: &Dense,
        until: f64,
    ) {
        let n = values.len();
        while out.saved_points < self.report_times.len()
            && self.report_times[out.saved_points] < until
        {
            let row = out.saved_points * n;
            out.values[row..row + n].copy_from_slice(values);
            let at = self.report_times[out.saved_points];
            let theta = ((at - dense.t) / dense.h).clamp(0.0, 1.0);
            for (j, &k) in dependent.iter().enumerate() {
                out.values[row + vars[k].index] = read(vars, k, dense.component(j, theta));
            }
            if let Some(shared) = shared {
                shared.each(at, |j, x| {
                    let k = shared.vars[j];
                    out.values[row + vars[k].index] = read(vars, k, x);
                });
            }
            out.saved_points += 1;
        }
    }

    /// The end of the shared solution, reached at `time`: nothing at the
    /// horizon, a stop at an independent bound, or the solve's failure.
    fn shared_end(
        shared: Option<&Shared>,
        vars: &[Variable],
        values: &mut [f64],
        time: f64,
    ) -> Result<Option<f64>, Halt> {
        match shared.map(|s| &s.ending) {
            None | Some(Ending::Horizon) => Ok(None),
            Some(&Ending::Bound(k)) => {
                values[vars[k].index] = vars[k].bound.expect("bounded variable").0;
                Ok(Some(time))
            }
            Some(Ending::Failed(halt)) => Err(halt.clone()),
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
        let shared = self.independent;
        let mut flow = Flow::new(net, self.model, shared);
        let vars = flow.vars;
        let dependent = flow.dependent.clone();
        let m = dependent.len();
        // Replicates integrate to the end of the shared solution.
        let horizon = shared.map_or(self.end, |s| s.until.min(self.end));
        let independent_stop = match shared.map(|s| &s.ending) {
            Some(&Ending::Bound(k)) => Some(k),
            _ => None,
        };
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
            for (j, &k) in dependent.iter().enumerate() {
                y[j] = values[vars[k].index];
            }
            y[m] = 0.0;
            let threshold = -libm::log(words.open()?);
            stepper.fresh = false;
            snapshot.clear();
            snapshot.extend_from_slice(values);
            let (stop, until) = {
                let fixed_counts: &[u64] = counts;
                let snapshot = snapshot.as_slice();
                let mut rhs = |t: f64, state: &[f64], dy: &mut [f64]| {
                    flow.rhs(fixed_counts, snapshot, t, state, dy)
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
                    let remaining = horizon - *time;
                    if remaining <= 0.0 {
                        // Reports at the horizon read `values`.
                        store(vars, &dependent, shared, &y, *time, values);
                        return Self::shared_end(shared, vars, values, *time);
                    }
                    let mut h = match predicted {
                        Some(h) => h,
                        None => stepper.initial_step(*time, &y, remaining, &mut rhs)?,
                    };
                    let rate = stepper.derivative(*time, &y, &mut rhs)?[m];
                    if rate > 0.0 {
                        // The threshold is expected within a quarter of the
                        // step; the cap never goes below a resolvable step.
                        let floor = 16.0 * f64::EPSILON * time.abs().max(1.0);
                        h = h.min((4.0 * (threshold - y[m]) / rate).max(floor));
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
                    // The earliest event inside the step. A stop wins a tie
                    // with the threshold; a clamp loses it, and the variable
                    // then reads the bound from the event on.
                    let mut theta = 1.0;
                    let mut stop = None;
                    if stepper.next[m] >= threshold {
                        theta = dense.crossing_up(m, threshold);
                        stop = Some(Stop::Hazard);
                    }
                    for (j, &k) in dependent.iter().enumerate() {
                        let Some((lower, action)) = vars[k].bound else {
                            continue;
                        };
                        let crossed = match action {
                            OnBound::Stop => stepper.next[j] <= lower,
                            OnBound::Clamp => stepper.next[j] < lower,
                        };
                        if y[j] > lower && crossed {
                            let at = dense.crossing_down(j, lower);
                            if at < theta || (at == theta && action == OnBound::Stop) {
                                theta = at;
                                stop = Some(Stop::Bound(j));
                            }
                        }
                    }
                    // An independent stop bound ends the shared solution: it
                    // stops the replicate at the horizon unless an event or
                    // a dependent stop comes first.
                    if let (true, Some(k)) = (last, independent_stop) {
                        let dependent_stop = matches!(
                            stop,
                            Some(Stop::Bound(j))
                                if matches!(vars[dependent[j]].bound, Some((_, OnBound::Stop)))
                        );
                        if theta == 1.0 && !dependent_stop {
                            stop = Some(Stop::Independent(k));
                        }
                    }
                    let until = if last && matches!(stop, None | Some(Stop::Independent(_))) {
                        horizon
                    } else {
                        (*time + theta * h).min(horizon)
                    };
                    self.save_dense(out, values, vars, &dependent, shared, &dense, until);
                    match stop {
                        None => {
                            y.copy_from_slice(&stepper.next);
                            *time = until;
                        }
                        Some(Stop::Bound(j))
                            if matches!(vars[dependent[j]].bound, Some((_, OnBound::Clamp))) =>
                        {
                            // Continue from the bound with the same threshold.
                            dense.state(theta, &mut y);
                            y[j] = vars[dependent[j]].bound.expect("bounded variable").0;
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
                Stop::Bound(j) => {
                    y[j] = vars[dependent[j]].bound.expect("bounded variable").0;
                    *time = until;
                    store(vars, &dependent, shared, &y, until, values);
                    return Ok(Some(until));
                }
                Stop::Independent(_) => {
                    *time = until;
                    store(vars, &dependent, shared, &y, until, values);
                    return Self::shared_end(shared, vars, values, until);
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
                    store(vars, &dependent, shared, &y, until, values);
                    let total = flow.all_rates(counts, values, &y, until)?;
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
