//! Fixed-channel trajectory loop, with the same result and control accounting.
use super::*;
use crate::flat::{Network, check_instances, direct, direct_with};
use crn_craft::next_reaction::Clocks;

impl Path<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn simulate_flat<R: RngCore, O: EventObserver>(
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
        if let Some(m) = monitor {
            m.start(&self.model.params, values);
        }
        if self.stop_values(values)? {
            return Ok(Some(0.0));
        }
        if matches!(
            self.method,
            crate::Method::Tau { .. } | crate::Method::Binomial { .. }
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
        let mut live = out.peak_instances;
        let mut cache = net.rate_cache();
        let mut rates = vec![0.0; net.channels.len()];
        let mut clocks = Clocks::default();
        let mut previous = None;
        if self.method == crate::Method::NextReaction {
            clocks.reset(rates.len());
        }
        let ziggurat = matches!(
            self.method,
            crate::Method::StaticDirect {
                waiting_time: Some(crate::RulesWaitingTime::Ziggurat),
                ..
            }
        );
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
            let candidate = if self.method == crate::Method::NextReaction {
                let mut update = |i, fired| {
                    let rate = cache.channel(net, i, counts, &env)?;
                    let fresh = if clocks.needs_draw(i, rate, fired) {
                        Some(-libm::log(1.0 - words.uniform()?))
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
                    direct_with(&rates, *time, words, |words, total| {
                        let unit = crn_craft::ssa::exponential_from_words(|retry| {
                            if retry {
                                self.poll()?;
                            }
                            words.word()
                        })?;
                        Ok(*time + unit / total)
                    })?
                } else {
                    direct(&rates, *time, words)?
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

    pub(crate) fn check_events(&self, out: &PathOutput, events: u64) -> Result<(), Halt> {
        if events <= self.limits.max_events.saturating_sub(out.events) {
            return Ok(());
        }
        Err(Halt::Limit(
            Termination::EventLimit,
            Diagnostic::error(
                "CRN-R040",
                "/budget/max_events_per_replicate",
                "event limit reached",
                "Increase max_events_per_replicate; incomplete ensembles withhold statistics.",
            ),
        ))
    }

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
        observer.on_event(&event) == ObserverAction::Stop
    }
}
