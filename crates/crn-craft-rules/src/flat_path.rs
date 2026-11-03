//! Fixed-channel trajectory loop, with the same result and control accounting.
use super::*;
use crate::flat::{Network, check_instances, direct};
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
        let mut rates = vec![0.0; net.channels.len()];
        let mut clocks = Clocks::default();
        let mut previous = None;
        if self.method == crate::Method::NextReaction {
            clocks.reset(rates.len());
        }
        loop {
            if out.events % POLL_EVENTS == 0 && out.events > 0 {
                self.poll()?;
            }
            let env = Env {
                params: &self.model.params,
                observables: values,
                attributes: [&[], &[]],
            };
            let candidate = if self.method == crate::Method::NextReaction {
                let mut update = |i, fired| {
                    let rate = net.rate(i, counts, &env)?;
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
                for (i, a) in rates.iter_mut().enumerate() {
                    *a = net.rate(i, counts, &env)?;
                }
                direct(&rates, *time, words)?
            };
            let Some((next, channel)) = candidate else {
                return Ok(None);
            };
            previous = Some(channel);
            self.save_before(out, values, next);
            if next > self.end {
                return Ok(None);
            }
            self.apply_flat(net, channel, counts, out)?;
            *time = next;
            net.refresh(counts, values)?;
            if self.observe_flat(net, channel, values, next, monitor, observer)
                || self.stop_values(values)?
            {
                return Ok(Some(next));
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
        out: &mut PathOutput,
    ) -> Result<(), Halt> {
        self.check_events(out, 1)?;
        let c = &net.channels[channel];
        let created = c
            .changes
            .iter()
            .filter(|&&(s, d)| d > 0 && counts[s] == 0)
            .count() as u64;
        net.apply(channel, counts, self.limits)?;
        out.events += 1;
        out.rule_events[c.rule] += 1;
        out.created_instances += created;
        out.peak_instances = out
            .peak_instances
            .max(counts.iter().filter(|&&n| n > 0).count() as u64);
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
        let reactants = participants(&c.reactants);
        let products = participants(&c.products);
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
