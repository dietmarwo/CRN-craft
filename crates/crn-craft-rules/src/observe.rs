//! Bounded, read-only event observation. Observation consumes no RNG words.
use crate::InstanceRecord;
use crate::expr::{Env, Expr, Scope};
use crate::model::{Compiled, Model};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObserverAction {
    Continue,
    Stop,
}

/// A committed state change. For a leap, rule is None and multiplicities
/// identify the aggregate firings. Reactants and products describe exact steps.
pub struct Event<'a> {
    pub replicate: u64,
    pub time: f64,
    pub rule: Option<usize>,
    pub reactants: &'a [(usize, &'a [f64])],
    pub products: &'a [(usize, &'a [f64])],
    pub observables: &'a [f64],
    pub firings: &'a [(usize, u64)],
}

pub trait EventObserver: Sync {
    const ENABLED: bool = true;
    fn on_event(&self, event: &Event<'_>) -> ObserverAction;
}

pub(crate) struct NoObserver;
impl EventObserver for NoObserver {
    const ENABLED: bool = false;
    #[inline]
    fn on_event(&self, _: &Event<'_>) -> ObserverAction {
        ObserverAction::Continue
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TrackerSpec {
    pub observable: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub first_passage: Option<String>,
    #[serde(default)]
    pub keep_episodes: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LogSpec {
    pub max_entries: usize,
    /// An empty list logs all rules. Leaps match any listed firing.
    #[serde(default)]
    pub rules: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Episode {
    pub start: f64,
    /// Open at the end of observation when absent.
    pub end: Option<f64>,
    pub peak: f64,
    /// Observed duration, including the observed part of an open episode.
    pub duration: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TrackerResult {
    pub observable: String,
    pub appearances: u64,
    pub extinctions: u64,
    pub last_appearance: Option<f64>,
    pub last_extinction: Option<f64>,
    pub first_passage: Option<f64>,
    pub minimum: f64,
    pub maximum: f64,
    pub integral: f64,
    pub episodes: Vec<Episode>,
    pub omitted_episodes: u64,
    pub omitted_duration: f64,
    pub omitted_peak: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LogEntry {
    pub time: f64,
    pub rule: Option<usize>,
    pub reactants: Vec<InstanceRecord>,
    pub products: Vec<InstanceRecord>,
    pub observables: Vec<f64>,
    pub firings: Vec<(usize, u64)>,
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Observation {
    pub id: u64,
    pub trackers: Vec<TrackerResult>,
    pub log: Vec<LogEntry>,
    pub omitted_log_entries: u64,
}

pub(crate) struct Spec {
    trackers: Vec<(TrackerSpec, usize, Option<Expr>)>,
    log: Option<(usize, Vec<usize>)>,
}

impl Spec {
    pub(crate) fn compile(
        model: &Model,
        compiled: &Compiled,
        trackers: &[TrackerSpec],
        log: Option<&LogSpec>,
    ) -> Result<Self, (&'static str, String)> {
        let names: BTreeMap<_, _> = compiled
            .observables
            .iter()
            .enumerate()
            .map(|(i, o)| (o.name.clone(), i))
            .collect();
        let params = model
            .parameters
            .keys()
            .cloned()
            .enumerate()
            .map(|(i, n)| (n, i))
            .collect();
        let scope = Scope {
            params: &params,
            observables: Some(&names),
            reactants: &[],
            random: false,
        };
        let trackers: Vec<_> = trackers
            .iter()
            .map(|t| {
                let index = *names
                    .get(&t.observable)
                    .ok_or_else(|| format!("unknown observable {}", t.observable))?;
                let predicate = t
                    .first_passage
                    .as_ref()
                    .map(|p| {
                        if p.len() > crate::model::MAX_EXPRESSION_BYTES {
                            return Err("first passage expression exceeds length limit".into());
                        }
                        Expr::parse(p, &scope)
                    })
                    .transpose()?;
                Ok::<_, String>((t.clone(), index, predicate))
            })
            .collect::<Result<_, String>>()
            .map_err(|e| ("/trackers", e))?;
        if let Some(atp) = &compiled.atp {
            for (_, index, predicate) in &trackers {
                let mut flags = vec![false; compiled.observables.len()];
                if let Some(expr) = predicate {
                    expr.reads(&mut flags);
                }
                if *index == atp.index || flags[atp.index] {
                    return Err((
                        "/trackers",
                        "trackers apply to event-driven observables; ATP is continuous".into(),
                    ));
                }
            }
        }
        let log = log
            .map(|l| {
                let indices = l
                    .rules
                    .iter()
                    .map(|name| {
                        compiled
                            .rules
                            .iter()
                            .position(|r| r.name == *name)
                            .ok_or_else(|| format!("unknown event log rule {name}"))
                    })
                    .collect::<Result<_, _>>()?;
                Ok::<_, String>((l.max_entries, indices))
            })
            .transpose()
            .map_err(|e| ("/event_log", e))?;
        Ok(Self { trackers, log })
    }

    pub(crate) fn enabled(&self) -> bool {
        !self.trackers.is_empty() || self.log.is_some()
    }

    /// Direct count trackers can integrate each constant span in one operation.
    pub(crate) fn count_trackers(&self, net: &crate::flat::Network) -> Option<[Vec<usize>; 3]> {
        if self.log.is_some() {
            return None;
        }
        let mut by_species: [Vec<usize>; 3] = Default::default();
        for (tracker, (_, index, predicate)) in self.trackers.iter().enumerate() {
            let (_, reduce, terms) = net.observables.iter().find(|(i, _, _)| i == index)?;
            if predicate.is_some() || !matches!(reduce, crate::model::Reduce::Sum) {
                return None;
            }
            let &[(species, 1.0)] = terms.as_slice() else {
                return None;
            };
            by_species[species].push(tracker);
        }
        Some(by_species)
    }

    pub(crate) fn reservation(&self, attributes: usize, observables: usize, rules: usize) -> u128 {
        let trackers = self
            .trackers
            .iter()
            .map(|(t, _, _)| 256 + t.keep_episodes as u128 * 64)
            .sum::<u128>();
        let entry = 160
            + observables as u128 * 8
            + rules as u128 * 16
            + (2 + crate::model::MAX_PRODUCTS as u128) * (48 + 8 * attributes as u128);
        trackers + self.log.as_ref().map_or(0, |(n, _)| *n as u128 * entry)
    }
}

struct Tracker {
    result: TrackerResult,
    value: f64,
    last_time: f64,
    current: Option<Episode>,
}

impl Tracker {
    fn update(&mut self, value: f64, time: f64, keep_episodes: usize) {
        self.result.integral += self.value * (time - self.last_time);
        if self.value <= 0.0 && value > 0.0 {
            self.result.appearances += 1;
            self.result.last_appearance = Some(time);
            self.current = Some(Episode {
                start: time,
                end: None,
                peak: value,
                duration: 0.0,
            });
        } else if self.value > 0.0 && value <= 0.0 {
            self.result.extinctions += 1;
            self.result.last_extinction = Some(time);
            self.close(time, true, keep_episodes);
        }
        if let Some(episode) = &mut self.current {
            episode.peak = episode.peak.max(value);
        }
        self.result.minimum = self.result.minimum.min(value);
        self.result.maximum = self.result.maximum.max(value);
        self.value = value;
        self.last_time = time;
    }

    fn close(&mut self, time: f64, closed: bool, keep: usize) {
        if let Some(mut episode) = self.current.take() {
            episode.end = closed.then_some(time);
            episode.duration = time - episode.start;
            if self.result.episodes.len() < keep {
                self.result.episodes.push(episode);
            } else {
                self.result.omitted_episodes += 1;
                self.result.omitted_duration += episode.duration;
                self.result.omitted_peak = Some(
                    self.result
                        .omitted_peak
                        .map_or(episode.peak, |peak| peak.max(episode.peak)),
                );
            }
        }
    }
}

pub(crate) struct Monitor<'a> {
    spec: &'a Spec,
    trackers: Vec<Tracker>,
    output: Observation,
}

impl<'a> Monitor<'a> {
    pub(crate) fn new(spec: &'a Spec, id: u64) -> Self {
        Self {
            spec,
            trackers: Vec::new(),
            output: Observation {
                id,
                ..Observation::default()
            },
        }
    }

    pub(crate) fn start(&mut self, params: &[f64], values: &[f64]) {
        let env = Env {
            params,
            observables: values,
            attributes: [&[], &[]],
        };
        for (spec, index, predicate) in &self.spec.trackers {
            let value = values[*index];
            self.trackers.push(Tracker {
                result: TrackerResult {
                    observable: spec.observable.clone(),
                    appearances: 0,
                    extinctions: 0,
                    last_appearance: None,
                    last_extinction: None,
                    first_passage: predicate
                        .as_ref()
                        .filter(|p| truth(p.eval(&env)))
                        .map(|_| 0.0),
                    minimum: value,
                    maximum: value,
                    integral: 0.0,
                    episodes: Vec::new(),
                    omitted_episodes: 0,
                    omitted_duration: 0.0,
                    omitted_peak: None,
                },
                value,
                last_time: 0.0,
                current: (value > 0.0).then_some(Episode {
                    start: 0.0,
                    end: None,
                    peak: value,
                    duration: 0.0,
                }),
            });
        }
    }

    pub(crate) fn keeps_exact_event(&self, rule: usize) -> bool {
        self.spec.log.as_ref().is_some_and(|(cap, rules)| {
            self.output.log.len() < *cap && (rules.is_empty() || rules.contains(&rule))
        })
    }

    pub(crate) fn observe(&mut self, event: &Event<'_>, params: &[f64]) {
        let env = Env {
            params,
            observables: event.observables,
            attributes: [&[], &[]],
        };
        for (tracker, (spec, index, predicate)) in self.trackers.iter_mut().zip(&self.spec.trackers)
        {
            let value = event.observables[*index];
            tracker.update(value, event.time, spec.keep_episodes);
            if tracker.result.first_passage.is_none()
                && predicate.as_ref().is_some_and(|p| truth(p.eval(&env)))
            {
                tracker.result.first_passage = Some(event.time);
            }
        }
        if let Some((max, rules)) = &self.spec.log {
            if rules.is_empty()
                || event.rule.is_some_and(|r| rules.contains(&r))
                || event.firings.iter().any(|(r, _)| rules.contains(r))
            {
                if self.output.log.len() < *max {
                    self.output.log.push(log_entry(event));
                } else {
                    self.output.omitted_log_entries += 1;
                }
            }
        }
    }

    pub(crate) fn observe_count(&mut self, indices: &[usize], value: f64, time: f64) {
        for &i in indices {
            self.trackers[i].update(value, time, self.spec.trackers[i].0.keep_episodes);
        }
    }

    pub(crate) fn finish(mut self, time: f64) -> Observation {
        self.output.trackers = self
            .trackers
            .into_iter()
            .zip(&self.spec.trackers)
            .map(|(mut tracker, (spec, _, _))| {
                tracker.result.integral += tracker.value * (time - tracker.last_time);
                tracker.close(time, false, spec.keep_episodes);
                tracker.result
            })
            .collect();
        self.output
    }
}

fn truth(value: f64) -> bool {
    !value.is_nan() && value != 0.0
}

pub fn log_entry(event: &Event<'_>) -> LogEntry {
    let records = |items: &[(usize, &[f64])]| {
        items
            .iter()
            .map(|(type_index, attrs)| InstanceRecord {
                type_index: *type_index,
                attributes: attrs.to_vec(),
                count: 1,
            })
            .collect()
    };
    LogEntry {
        time: event.time,
        rule: event.rule,
        reactants: records(event.reactants),
        products: records(event.products),
        observables: event.observables.to_vec(),
        firings: event.firings.to_vec(),
    }
}
