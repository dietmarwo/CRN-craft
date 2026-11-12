//! Fixed-channel trajectory loop, with the same result and control accounting.
use super::*;
use crate::next_reaction::Clocks;
use crate::rules::flat::{Network, check_instances, direct, direct_with, unresolved};
use crate::rules::rssa::Bounds;
use crate::rules::select::{Composition, Kind, RateTree, sorting_sample};

impl Path<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn simulate_flat<R: RngCore + 'static, O: EventObserver>(
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
        self.poll()?;

        out.peak_instances = counts.iter().filter(|&&n| n > 0).count() as u64;
        check_instances(out.peak_instances, self.limits)?;
        net.refresh(counts, values)?;
        if let Some(atp) = &self.model.atp {
            values[atp.index] = atp.initial;
        }
        for variable in &self.model.continuous {
            values[variable.index] = variable.initial;
        }
        if let Some(m) = monitor {
            m.start(&self.model.params, values);
        }
        if self.stop_values(values)? {
            return Ok(Some(0.0));
        }
        if let crate::rules::Method::Hybrid { rtol, atol } = self.method {
            return self.simulate_hybrid(
                rtol, atol, net, counts, values, words, out, time, monitor, observer,
            );
        }
        if matches!(
            self.method,
            crate::rules::Method::Tau { .. } | crate::rules::Method::Binomial { .. }
        ) {
            return self.simulate_leaps(net, counts, values, words, out, time, monitor, observer);
        }
        if let Some(atp) = &self.model.atp {
            return self.simulate_thinning(
                atp, net, counts, values, words, out, time, monitor, observer,
            );
        }
        if let Some(kernel) = self.mtdna {
            let stopped = self.simulate_mtdna(
                kernel, net, counts, values, words, out, time, monitor, observer,
            )?;
            net.refresh(counts, values)?;
            return Ok(stopped);
        }
        if let Some(kernel) = self.generated_loop {
            return self
                .simulate_generated_loop(kernel, net, counts, values, words, out, time, monitor);
        }
        if self.selector != Kind::Linear {
            return self
                .simulate_selected(net, counts, values, words, out, time, monitor, observer);
        }
        if self.generated.is_none() {
            match self.method {
                crate::rules::Method::StaticDirect {
                    waiting_time: Some(crate::rules::RulesWaitingTime::Ziggurat),
                    ..
                } => {
                    return self.simulate_linear::<R, O, true>(
                        net, counts, values, words, out, time, monitor, observer,
                    );
                }
                crate::rules::Method::StaticDirect { .. } => {
                    return self.simulate_linear::<R, O, false>(
                        net, counts, values, words, out, time, monitor, observer,
                    );
                }
                _ => {}
            }
        }
        let mut live = out.peak_instances;
        let mut cache = net.rate_cache();
        let mut rates = vec![0.0; net.channels.len()];
        let mut clocks = Clocks::default();
        let mut previous = None;
        if matches!(self.method, crate::rules::Method::NextReaction {}) {
            clocks.reset(rates.len());
        }
        let ziggurat = matches!(
            self.method,
            crate::rules::Method::StaticDirect {
                waiting_time: Some(crate::rules::RulesWaitingTime::Ziggurat),
                ..
            }
        );
        let mut clock = Clock::new(self.lenient_clock);
        loop {
            if out.events % POLL_EVENTS == 0 && out.events > 0 {
                self.poll()?;
            }
            let env = Env {
                params: &self.model.params,
                observables: values,
                attributes: [&[], &[]],
            };
            cache.reset();
            let candidate = if matches!(self.method, crate::rules::Method::NextReaction {}) {
                let mut update = |i, fired| {
                    let rate = cache.channel(net, i, counts, &env)?;
                    let fresh = if clocks.needs_draw(i, rate, fired) {
                        Some(-libm::log(1.0 - words.open()?))
                    } else {
                        None
                    };
                    if !clocks.update(i, rate, *time, fresh) {
                        return Err(numerical(
                            "/time",
                            "reaction clock failed to advance",
                            "Rescale rates and time.",
                        ));
                    }
                    Ok(())
                };
                if let Some(i) = previous {
                    update(i, true)?;
                    for &j in &net.affected[i] {
                        if i != j {
                            update(j, false)?;
                        }
                    }
                } else {
                    for i in 0..rates.len() {
                        update(i, false)?;
                    }
                }
                clocks.next()
            } else {
                if let Some(kernel) = self.generated {
                    kernel
                        .rates(previous, counts, values, &mut rates)
                        .map_err(|i| net.rate_error(i))?;
                } else if let Some(channel) = previous {
                    for &i in &net.affected[channel] {
                        rates[i] = cache.channel(net, i, counts, &env)?;
                    }
                } else {
                    for (i, rate) in rates.iter_mut().enumerate() {
                        *rate = cache.channel(net, i, counts, &env)?;
                    }
                }
                if ziggurat {
                    direct_with(&rates, *time, words, &mut clock, |words, total| {
                        let unit = crate::ssa::exponential_from_words(|retry| {
                            if retry {
                                self.poll()?;
                            }
                            words.word()
                        })?;
                        Ok(unit / total)
                    })?
                } else {
                    direct(&rates, *time, words, &mut clock)?
                }
            };
            let Some((next, channel)) = candidate else {
                net.refresh_cold(counts, values)?;
                return Ok(None);
            };
            previous = Some(channel);
            if out.saved_points < self.report_times.len()
                && self.report_times[out.saved_points] < next
            {
                net.refresh_cold(counts, values)?;
                self.save_before(out, values, next);
            }
            if next > self.end {
                net.refresh_cold(counts, values)?;
                return Ok(None);
            }
            self.apply_flat(net, channel, counts, &mut live, out)?;
            *time = next;
            if let Some(kernel) = self.generated {
                net.refresh_generated(channel, counts, values, O::ENABLED, kernel)?;
            } else {
                net.refresh_channel(channel, counts, values, O::ENABLED)?;
            }
            if (O::ENABLED || net.observe_affected[channel])
                && self.observe_flat(net, channel, values, next, monitor, observer)
            {
                return Ok(Some(next));
            }
            if net.stop_affected[channel] {
                match self.stop_values(values) {
                    Ok(false) => {}
                    Ok(true) => {
                        net.refresh_cold(counts, values)?;
                        return Ok(Some(next));
                    }
                    Err(_) => {
                        // The diagnostic includes all projection values.
                        net.refresh_cold(counts, values)?;
                        self.stop_values(values)?;
                    }
                }
            }
        }
    }

    /// Static Direct with the linear scan, without a generated kernel: the
    /// arithmetic, checks and random words of [`direct_with`] in the loop of
    /// [`Self::simulate_flat`], with the running totals kept between events.
    ///
    /// The total is the sequential sum of the channel rates from `-0.0`
    /// (`Iterator::sum`). After an event only the channels in the fired
    /// channel's `affected` list change; they are sorted, so the running
    /// sums before the first of them are the ones a full pass would compute,
    /// and the pass resumes there. The smallest positive rate decides the
    /// resolution check exactly as testing every rate does.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn simulate_linear<R: RngCore + 'static, O: EventObserver, const ZIGGURAT: bool>(
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
        let n = net.channels.len();
        let params = &self.model.params;
        let mut live = out.peak_instances;
        let mut cache = net.rate_cache();
        let mut rates = vec![0.0; n];
        // Running total and smallest positive rate through each channel.
        let mut sums = vec![-0.0; n];
        let mut minima = vec![f64::INFINITY; n];
        let mut clock = Clock::new(self.lenient_clock);
        let mut previous: Option<usize> = None;
        loop {
            if out.events % POLL_EVENTS == 0 && out.events > 0 {
                self.poll()?;
            }
            let env = Env {
                params,
                observables: values,
                attributes: [&[], &[]],
            };
            cache.reset();
            let first = if let Some(channel) = previous {
                let affected = &net.affected[channel];
                for &i in affected {
                    rates[i] = cache.channel(net, i, counts, &env)?;
                }
                affected.first().copied()
            } else {
                for (i, rate) in rates.iter_mut().enumerate() {
                    *rate = cache.channel(net, i, counts, &env)?;
                }
                Some(0)
            };
            if let Some(first) = first {
                let (mut sum, mut smallest) = if first == 0 {
                    (-0.0, f64::INFINITY)
                } else {
                    (sums[first - 1], minima[first - 1])
                };
                for i in first..n {
                    let a = rates[i];
                    sum += a;
                    if a > 0.0 && a < smallest {
                        smallest = a;
                    }
                    sums[i] = sum;
                    minima[i] = smallest;
                }
            }
            let (total, smallest) = match n {
                0 => (-0.0, f64::INFINITY),
                _ => (sums[n - 1], minima[n - 1]),
            };
            if total == 0.0 {
                net.refresh_cold(counts, values)?;
                return Ok(None);
            }
            if !total.is_finite() || smallest < total * crate::rules::engine::RESOLUTION {
                return Err(unresolved(total));
            }
            let wait = if ZIGGURAT {
                crate::ssa::exponential_from_words(|retry| {
                    if retry {
                        self.poll()?;
                    }
                    words.word()
                })? / total
            } else {
                -libm::log(1.0 - words.open()?) / total
            };
            let next = clock.advance(*time, wait)?;
            let mut target = words.uniform()? * total;
            let mut channel = usize::MAX;
            for (i, &a) in rates.iter().enumerate() {
                if a > 0.0 {
                    channel = i;
                    if target < a {
                        break;
                    }
                    target -= a;
                }
            }
            previous = Some(channel);
            if out.saved_points < self.report_times.len()
                && self.report_times[out.saved_points] < next
            {
                net.refresh_cold(counts, values)?;
                self.save_before(out, values, next);
            }
            if next > self.end {
                net.refresh_cold(counts, values)?;
                return Ok(None);
            }
            self.apply_flat(net, channel, counts, &mut live, out)?;
            *time = next;
            net.refresh_channel(channel, counts, values, O::ENABLED)?;
            if (O::ENABLED || net.observe_affected[channel])
                && self.observe_flat(net, channel, values, next, monitor, observer)
            {
                return Ok(Some(next));
            }
            if net.stop_affected[channel] {
                match self.stop_values(values) {
                    Ok(false) => {}
                    Ok(true) => {
                        net.refresh_cold(counts, values)?;
                        return Ok(Some(next));
                    }
                    Err(_) => {
                        // The diagnostic includes all projection values.
                        net.refresh_cold(counts, values)?;
                        self.stop_values(values)?;
                    }
                }
            }
        }
    }

    /// Static Direct with composition–rejection, sorting direct or RSSA-CR
    /// selection (`crn-craft.run/1.0`). Waiting times, saves, the horizon,
    /// events, observation and stops are those of the linear scan; RSSA-CR
    /// also advances the clock on rejected candidates.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn simulate_selected<R: RngCore + 'static, O: EventObserver>(
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
        let kind = self.selector;
        let n = net.channels.len();
        let params = &self.model.params;
        let mut live = out.peak_instances;
        let mut cache = net.rate_cache();
        let mut tree = RateTree::default();
        tree.reset(n);
        let mut composition = Composition::default();
        let mut rates = Vec::new();
        let mut order = Vec::new();
        let mut bounds = Bounds::default();
        match kind {
            Kind::Sorting => {
                rates.resize(n, 0.0);
                order.extend(0..n);
            }
            Kind::RssaCr => {
                composition.reset(n);
                bounds.reset(net, self.rssa.expect("RSSA-CR plan"), counts);
            }
            _ => composition.reset(n),
        }
        let ziggurat = matches!(
            self.method,
            crate::rules::Method::StaticDirect {
                waiting_time: Some(crate::rules::RulesWaitingTime::Ziggurat),
                ..
            }
        );
        let mut clock = Clock::new(self.lenient_clock);
        // `Some` after an event: the state changed since the last pass.
        let mut fired: Option<usize> = None;
        let mut started = false;
        let mut absorbed = false;
        let mut iterations: u64 = 0;
        loop {
            iterations += 1;
            if iterations % POLL_EVENTS == 0 {
                self.poll()?;
            }
            if !started || fired.is_some() {
                cache.reset();
                if kind == Kind::RssaCr {
                    let plan = self.rssa.expect("RSSA-CR plan");
                    bounds.refresh(
                        net,
                        plan,
                        fired,
                        counts,
                        values,
                        params,
                        &mut cache,
                        &mut tree,
                        &mut composition,
                    )?;
                    // Every upper bound may be positive while every
                    // propensity is zero.
                    absorbed = false;
                    if bounds.positive_lower == 0 {
                        let env = Env {
                            params,
                            observables: values,
                            attributes: [&[], &[]],
                        };
                        absorbed = true;
                        for i in 0..n {
                            if cache.channel(net, i, counts, &env)? > 0.0 {
                                absorbed = false;
                                break;
                            }
                        }
                    }
                } else {
                    let env = Env {
                        params,
                        observables: values,
                        attributes: [&[], &[]],
                    };
                    let mut update = |i: usize| -> Result<(), Halt> {
                        let rate = cache.channel(net, i, counts, &env)?;
                        tree.set(i, rate, if rate > 0.0 { rate } else { f64::INFINITY });
                        if kind == Kind::Sorting {
                            rates[i] = rate;
                        } else {
                            composition.update(i, rate)?;
                        }
                        Ok(())
                    };
                    match fired {
                        Some(channel) => {
                            for &i in &net.affected[channel] {
                                update(i)?;
                            }
                        }
                        None => {
                            for i in 0..n {
                                update(i)?;
                            }
                        }
                    }
                }
                started = true;
                fired = None;
            }
            let total = if absorbed { 0.0 } else { tree.total()? };
            if total == 0.0 {
                net.refresh_cold(counts, values)?;
                return Ok(None);
            }
            let wait = if ziggurat {
                let unit = crate::ssa::exponential_from_words(|retry| {
                    if retry {
                        self.poll()?;
                    }
                    words.word()
                })?;
                unit / total
            } else {
                -libm::log(1.0 - words.open()?) / total
            };
            let next = clock.advance(*time, wait)?;
            let channel = if kind == Kind::Sorting {
                sorting_sample(&mut order, &rates, total, words)?
            } else {
                composition.sample(&tree, words)?
            };
            if out.saved_points < self.report_times.len()
                && self.report_times[out.saved_points] < next
            {
                net.refresh_cold(counts, values)?;
                self.save_before(out, values, next);
            }
            if next > self.end {
                net.refresh_cold(counts, values)?;
                return Ok(None);
            }
            if kind == Kind::RssaCr {
                // Thinning against the candidate's bounds: accept below the
                // lower bound, otherwise against the exact propensity.
                let upper = tree.get(channel);
                let threshold = words.uniform()? * upper;
                if threshold >= upper {
                    return Err(numerical(
                        "/rules",
                        "RSSA-CR acceptance rounded to its bound",
                        "Rescale rates or populations.",
                    ));
                }
                let lower = bounds.lower[channel];
                if threshold >= lower {
                    let env = Env {
                        params,
                        observables: values,
                        attributes: [&[], &[]],
                    };
                    let a = cache.channel(net, channel, counts, &env)?;
                    if a < lower || a > upper {
                        return Err(numerical(
                            "/rules",
                            "a propensity left its RSSA-CR bounds",
                            "Select the linear selector and report the model.",
                        ));
                    }
                    if threshold >= a {
                        *time = next;
                        continue;
                    }
                }
            }
            fired = Some(channel);
            self.apply_flat(net, channel, counts, &mut live, out)?;
            *time = next;
            net.refresh_channel(channel, counts, values, O::ENABLED)?;
            if (O::ENABLED || net.observe_affected[channel])
                && self.observe_flat(net, channel, values, next, monitor, observer)
            {
                return Ok(Some(next));
            }
            if net.stop_affected[channel] {
                match self.stop_values(values) {
                    Ok(false) => {}
                    Ok(true) => {
                        net.refresh_cold(counts, values)?;
                        return Ok(Some(next));
                    }
                    Err(_) => {
                        net.refresh_cold(counts, values)?;
                        self.stop_values(values)?;
                    }
                }
            }
        }
    }

    pub(super) fn save_before(&self, out: &mut PathOutput, values: &[f64], time: f64) {
        let n = values.len();
        while out.saved_points < self.report_times.len()
            && self.report_times[out.saved_points] < time
        {
            let row = out.saved_points * n;
            out.values[row..row + n].copy_from_slice(values);
            out.saved_points += 1;
        }
    }

    #[inline]
    pub(crate) fn check_events(&self, out: &PathOutput, events: u64) -> Result<(), Halt> {
        if events <= self.limits.max_events.saturating_sub(out.events) {
            return Ok(());
        }
        Err(event_limit())
    }

    #[inline]
    pub(super) fn apply_flat(
        &self,
        net: &Network,
        channel: usize,
        counts: &mut [u64],
        live: &mut u64,
        out: &mut PathOutput,
    ) -> Result<(), Halt> {
        self.check_events(out, 1)?;
        let c = &net.channels[channel];
        let created = net.apply(channel, counts, live, self.limits)?;
        out.events += 1;
        out.rule_events[c.rule] += 1;
        out.created_instances += created;
        out.peak_instances = out.peak_instances.max(*live);
        Ok(())
    }

    pub(super) fn observe_flat<O: EventObserver>(
        &self,
        net: &Network,
        channel: usize,
        values: &[f64],
        time: f64,
        monitor: &mut Option<Monitor<'_>>,
        observer: &O,
    ) -> bool {
        if !O::ENABLED && monitor.is_none() {
            return false;
        }
        let c = &net.channels[channel];
        let participants = |ids: &[usize]| {
            ids.iter()
                .map(|&s| {
                    let species = &net.species[s];
                    (species.ty, species.attrs.as_slice())
                })
                .collect::<Vec<_>>()
        };
        let (reactants, products) = if O::ENABLED
            || monitor
                .as_ref()
                .is_some_and(|m| m.keeps_exact_event(c.rule))
        {
            (participants(&c.reactants), participants(&c.products))
        } else {
            (Vec::new(), Vec::new())
        };
        let event = Event {
            replicate: self.id,
            time,
            rule: Some(c.rule),
            reactants: &reactants,
            products: &products,
            observables: values,
            firings: &[],
        };
        if let Some(m) = monitor {
            m.observe(&event, &self.model.params);
        }
        O::ENABLED && observer.on_event(&event) == ObserverAction::Stop
    }
}

#[cold]
#[inline(never)]
fn event_limit() -> Halt {
    Halt::Limit(
        Termination::EventLimit,
        Diagnostic::error(
            "CRN-R040",
            "/budget/max_events_per_replicate",
            "event limit reached",
            "Increase max_events_per_replicate; incomplete ensembles withhold statistics.",
        ),
    )
}
