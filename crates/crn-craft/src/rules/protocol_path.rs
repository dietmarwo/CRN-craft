//! A trajectory whose reaction clock competes with a protocol and completions.
use super::flat_path::event_limit;
use super::leap::binomial;
use super::*;
use crate::rules::clocks::Clocks;
use crate::rules::flat::{Network, check_instances, direct, direct_with};
use crate::rules::observe::{ProtocolEvent, ProtocolEventKind};
use crate::rules::protocol::{Completion, Operation, Schedule};
use std::collections::BTreeMap;

// Positive finite times have the same order as their bit patterns.
type Key = (u64, u64);

struct Pending {
    jobs: BTreeMap<Key, (usize, usize)>,
    owners: Vec<Vec<Key>>,
    serial: u64,
}

impl Pending {
    fn new(species: usize) -> Self {
        Self {
            jobs: BTreeMap::new(),
            owners: vec![Vec::new(); species],
            serial: 0,
        }
    }

    fn add(&mut self, time: f64, channel: usize, completion: &Completion) {
        let key = (time.to_bits(), self.serial);
        self.serial += 1;
        let slot = completion.owner.map_or(0, |s| {
            let list = &mut self.owners[s];
            list.push(key);
            list.len() - 1
        });
        self.jobs.insert(key, (channel, slot));
    }

    fn remove(&mut self, key: Key, owner: Option<usize>) {
        let (_, slot) = self.jobs.remove(&key).expect("live completion");
        if let Some(s) = owner {
            let list = &mut self.owners[s];
            list.swap_remove(slot);
            if let Some(moved) = list.get(slot) {
                self.jobs.get_mut(moved).expect("live owner").1 = slot;
            }
        }
    }

    fn cancel<R: RngCore>(
        &mut self,
        s: usize,
        n: u64,
        words: &mut Words<'_, R>,
    ) -> Result<(), Halt> {
        for _ in 0..n {
            let count = self.owners[s].len() as u64;
            if count == 0 {
                return Err(protocol_error(
                    "busy count has no corresponding pending completion",
                ));
            }
            let slot = if count == 1 {
                0
            } else {
                // Rejection of the incomplete residue interval gives a
                // uniform identity even when the owner count is not a power of two.
                let floor = count.wrapping_neg() % count;
                loop {
                    let x = words.word()?;
                    if x >= floor {
                        break (x % count) as usize;
                    }
                }
            };
            self.remove(self.owners[s][slot], Some(s));
        }
        Ok(())
    }
}

fn protocol_error(message: impl Into<String>) -> Halt {
    numerical(
        "/protocol",
        message,
        "Check counts, durations and the pending completion budget.",
    )
}

fn action_error(index: usize, message: String) -> Halt {
    numerical(
        &format!("/protocol/actions/{index}"),
        message,
        "Check the counts that the action group reads and the pending completion budget.",
    )
}

fn count_add(n: u64, change: i64) -> Result<u64, Halt> {
    n.checked_add_signed(change)
        .filter(|&x| x <= MAX_COUNT)
        .ok_or_else(|| protocol_error(format!("count {n} cannot change by {change}")))
}

fn commit(
    counts: &mut [u64],
    next: &[u64],
    limits: &Limits,
    out: &mut PathOutput,
) -> Result<u64, Halt> {
    let live = next.iter().filter(|&&n| n > 0).count() as u64;
    check_instances(live, limits)?;
    out.created_instances += next
        .iter()
        .zip(&*counts)
        .filter(|&(&n, &old)| n > 0 && old == 0)
        .count() as u64;
    out.peak_instances = out.peak_instances.max(live);
    counts.copy_from_slice(next);
    Ok(live)
}

fn operation<R: RngCore>(
    n: u64,
    operation: &Operation,
    words: &mut Words<'_, R>,
) -> Result<u64, Halt> {
    match *operation {
        Operation::Change(change) => count_add(n, change),
        Operation::Set(count) => Ok(count),
        Operation::Partition(p) => binomial(n, p, words),
    }
}

impl Path<'_> {
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn simulate_protocol<R: RngCore, O: EventObserver>(
        &self,
        schedule: &Schedule,
        base: &Network,
        counts: &mut [u64],
        values: &mut [f64],
        words: &mut Words<'_, R>,
        out: &mut PathOutput,
        time: &mut f64,
        monitor: &mut Option<Monitor<'_>>,
        observer: &O,
    ) -> Result<Option<f64>, Halt> {
        use crate::rules::{Method, RulesWaitingTime};
        let mut net = base;
        let mut model = self.model;
        let mut live = counts.iter().filter(|&&n| n > 0).count() as u64;
        out.peak_instances = live;
        check_instances(live, self.limits)?;
        net.refresh(counts, values)?;
        if let Some(monitor) = monitor {
            monitor.start(&model.params, values);
        }
        let mut actions = schedule.boundaries.iter().peekable();
        if actions.peek().is_none_or(|a| a.time != 0.0) && Self::stop_values_with(model, values)? {
            return Ok(Some(0.0));
        }
        let mut pending = Pending::new(counts.len());
        let mut next_counts = counts.to_vec();
        let mut deterministic_events = 0_u64;
        let mut cache = net.rate_cache();
        let mut rates = vec![0.0; net.channels.len()];
        let use_clocks = matches!(self.method, Method::NextReaction {});
        let mut clocks = Clocks::default();
        clocks.reset(rates.len());
        let mut clock = Clock::new(self.lenient_clock);
        let ziggurat = matches!(
            self.method,
            Method::StaticDirect {
                waiting_time: Some(RulesWaitingTime::Ziggurat),
                ..
            }
        );
        let mut previous = None;
        loop {
            if (out.events + deterministic_events) % POLL_EVENTS == 0 {
                self.poll()?;
            }
            let next_action = actions.peek().map_or(f64::INFINITY, |a| a.time);
            let job = pending
                .jobs
                .first_key_value()
                .map(|(&key, &(channel, _))| (key, channel));
            let next_completion = job.map_or(f64::INFINITY, |(key, _)| f64::from_bits(key.0));
            let boundary = next_action.min(next_completion);
            let candidate = if boundary == *time {
                None
            } else {
                cache.reset();
                let env = Env {
                    params: &model.params,
                    observables: values,
                    attributes: [&[], &[]],
                };
                let mut update = |i, fired| -> Result<(), Halt> {
                    let a = cache.channel(net, i, counts, &env)?;
                    rates[i] = a;
                    if use_clocks {
                        let fresh = if clocks.needs_draw(i, a, fired) {
                            Some(-libm::log1p(-words.open()?))
                        } else {
                            None
                        };
                        if !clocks.update(i, a, *time, fresh, self.lenient_clock) {
                            return Err(protocol_error("reaction clock failed to advance"));
                        }
                    }
                    Ok(())
                };
                if let Some(channel) = previous {
                    if use_clocks {
                        update(channel, true)?;
                    }
                    for &i in &net.affected[channel] {
                        if !use_clocks || i != channel {
                            update(i, false)?;
                        }
                    }
                } else {
                    for i in 0..net.channels.len() {
                        update(i, false)?;
                    }
                }
                if use_clocks {
                    clocks.next()
                } else if ziggurat {
                    direct_with(&rates, *time, words, &mut clock, |words, total| {
                        Ok(crate::ssa::exponential_from_words(|retry| {
                            if retry {
                                self.poll()?;
                            }
                            words.word()
                        })? / total)
                    })?
                } else {
                    direct(&rates, *time, words, &mut clock)?
                }
            };
            let tentative = candidate.map_or(f64::INFINITY, |(t, _)| t);
            let next = boundary.min(tentative);
            if self
                .report_times
                .get(out.saved_points)
                .is_some_and(|&t| t < next)
            {
                net.refresh_cold(counts, values)?;
                self.save_before(out, values, next);
            }
            if next > self.end {
                net.refresh_cold(counts, values)?;
                return Ok(None);
            }
            if deterministic_events >= self.limits.max_events.saturating_sub(out.events) {
                return Err(event_limit());
            }
            let before_pending = pending.jobs.len();
            if boundary <= tentative {
                next_counts.copy_from_slice(counts);
                let kind = if next_completion <= next_action {
                    ProtocolEventKind::Completion
                } else {
                    ProtocolEventKind::Intervention
                };
                if next_completion <= next_action {
                    let (key, channel) = job.expect("there is a due completion");
                    let completion = schedule.delays[channel].as_ref().expect("known completion");
                    if let Some(s) = completion.owner {
                        next_counts[s] = count_add(next_counts[s], -1)?;
                    }
                    for &(s, n) in &completion.products {
                        next_counts[s] = count_add(next_counts[s], n as i64)?;
                    }
                    pending.remove(key, completion.owner);
                } else {
                    let action = actions.next().expect("there is a scheduled action");
                    let name = |s: usize| &model.types[net.species[s].ty].name;
                    for (s, op) in &action.operations {
                        // Only a count change can fail on the counts; other errors,
                        // such as the sampling budget of a partition, pass through.
                        let fails = matches!(op, Operation::Change(_));
                        next_counts[*s] = operation(counts[*s], op, words).map_err(|e| {
                            if !fails {
                                return e;
                            }
                            action_error(
                                action.index,
                                format!(
                                    "action group {} at time {}: {:?} has {} copies; the operation \
                                     would leave it negative or above the count limit",
                                    action.index,
                                    action.time,
                                    name(*s),
                                    counts[*s]
                                ),
                            )
                        })?;
                    }
                    for &s in &schedule.owned {
                        let removed = counts[s].checked_sub(next_counts[s]).ok_or_else(|| {
                            action_error(
                                action.index,
                                format!(
                                    "action group {} at time {}: an intervention cannot create \
                                     busy {:?} molecules",
                                    action.index,
                                    action.time,
                                    name(s)
                                ),
                            )
                        })?;
                        pending.cancel(s, removed, words)?;
                    }
                    if action.phase > 0 {
                        let phase = &schedule.phases[action.phase - 1];
                        net = &phase.network;
                        model = &phase.model;
                    }
                    cache = net.rate_cache();
                }
                let log_counts = monitor.as_ref().is_some_and(|m| m.keeps_protocol_event());
                let before = log_counts.then(|| counts.to_vec());
                live = commit(counts, &next_counts, self.limits, out)?;
                *time = next;
                deterministic_events += 1;
                net.refresh(counts, values)?;
                let event = Event {
                    replicate: self.id,
                    time: next,
                    rule: None,
                    reactants: &[],
                    products: &[],
                    observables: values,
                    firings: &[],
                    protocol: Some(ProtocolEvent {
                        kind,
                        pending: pending.jobs.len(),
                        cancelled: if next_completion <= next_action {
                            0
                        } else {
                            before_pending - pending.jobs.len()
                        },
                    }),
                };
                if let Some(m) = monitor {
                    m.observe_protocol(&event, &model.params, before.as_deref(), counts, net);
                }
                if O::ENABLED && observer.on_event(&event) == ObserverAction::Stop {
                    return Ok(Some(next));
                }
                previous = None;
                if use_clocks {
                    clocks.reset(rates.len());
                }
            } else {
                let (_, channel) = candidate.expect("there is a reaction");
                // A firing time of the next-reaction method that equals the
                // current time is a wait below the clock's resolution.
                if use_clocks {
                    clock.arrive(*time, next)?;
                }
                let completion = schedule.delays[channel]
                    .as_ref()
                    .map(|c| (c, next + c.duration));
                if let Some((_, due)) = completion {
                    if !due.is_finite() || due <= next || pending.jobs.len() >= schedule.max_pending
                    {
                        return Err(protocol_error(
                            "completion clock or pending work limit exceeded",
                        ));
                    }
                }
                for &(s, n) in &net.channels[channel].changes {
                    if n < 0 && schedule.owned.contains(&s) {
                        pending.cancel(s, n.unsigned_abs(), words)?;
                    }
                }
                self.apply_flat(net, channel, counts, &mut live, out)?;
                *time = next;
                let cancelled = before_pending - pending.jobs.len();
                let mut protocol = None;
                if let Some((spec, due)) = completion {
                    pending.add(due, channel, spec);
                    protocol = Some(ProtocolEvent {
                        kind: ProtocolEventKind::Initiation,
                        pending: pending.jobs.len(),
                        cancelled,
                    });
                } else if cancelled > 0 {
                    protocol = Some(ProtocolEvent {
                        kind: ProtocolEventKind::Cancellation,
                        pending: pending.jobs.len(),
                        cancelled,
                    });
                }
                net.refresh_channel(channel, counts, values, O::ENABLED)?;
                if (O::ENABLED || net.observe_affected[channel])
                    && self.observe_flat_with_params(
                        net,
                        channel,
                        values,
                        next,
                        monitor,
                        observer,
                        &model.params,
                        protocol,
                    )
                {
                    net.refresh_cold(counts, values)?;
                    return Ok(Some(next));
                }
                previous = Some(channel);
            }
            if Self::stop_values_with(model, values)? {
                net.refresh_cold(counts, values)?;
                return Ok(Some(next));
            }
        }
    }
}
