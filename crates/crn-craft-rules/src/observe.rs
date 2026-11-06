//! Bounded, read-only event observation. Observation consumes no RNG words.
//!
//! There are two ways to watch a run's events:
//!
//! - **In the result** (rules-run/0.2 and later): [`TrackerSpec`] entries
//!   summarize one observable per replicate (crossings of zero, extrema,
//!   time integral, first passage and episodes), and a [`LogSpec`] retains a
//!   bounded number of events. Both appear per replicate in
//!   [`crate::RulesResult::observations`] as an [`Observation`]. Their
//!   retention counts enter the admission check of [`crate::prepare`].
//! - **Live**: an [`EventObserver`] passed to
//!   [`crate::execute_with_observer`] receives every committed [`Event`] on
//!   the worker threads and can stop a replicate.
//!
//! Observation never changes the numbers: a run with trackers, a log or an
//! observer that always continues draws the same random words and reports
//! the same observable values as a run without them. Generated kernels and
//! loops support trackers and logs but not observers.
//!
//! ```
//! use crn_craft_rules::{Control, execute, prepare};
//! use serde_json::json;
//!
//! // Molecules of A appear one at a time and die; an episode is a span with
//! // at least one molecule present.
//! let model = serde_json::to_vec(&json!({
//!     "schema": "crn-craft.rules-model/0.1", "name": "flicker",
//!     "parameters": {"k_in": 0.5, "k_out": 2.0},
//!     "types": [{"name": "A"}], "initial": [],
//!     "observables": {"A": {"type": "A"}},
//!     "rules": [
//!         {"name": "birth", "reactants": [], "rate": "k_in", "products": [{"type": "A"}]},
//!         {"name": "death", "reactants": [{"label": "a", "type": "A"}],
//!          "rate": "k_out", "products": []}
//!     ]
//! }))
//! .unwrap();
//! let request = serde_json::to_vec(&json!({
//!     "schema": "crn-craft.rules-run/0.2",
//!     "model": {"path": "flicker.rules.json", "sha256": crn_craft::sha256(&model)},
//!     "time": {"end": 20.0, "report_times": [20.0]},
//!     "ensemble": {"replicates": 2, "root_seed": 3},
//!     "outputs": {"statistics": [], "keep_paths": 0, "keep_instances": 0},
//!     "budget": {"max_events_per_replicate": 10000, "max_sampling_steps_per_replicate": 100000,
//!                "max_instances": 1, "threads": 1, "max_output_bytes": 1048576},
//!     "trackers": [{"observable": "A", "first_passage": "A >= 2", "keep_episodes": 3}],
//!     "event_log": {"max_entries": 5}
//! }))
//! .unwrap();
//! let result = execute(&prepare(&model, &request).unwrap(), &Control::default()).unwrap();
//!
//! for observation in &result.observations {
//!     let a = &observation.trackers[0];
//!     // Every appearance opens an episode; the retained ones come first.
//!     assert_eq!(a.appearances, a.episodes.len() as u64 + a.omitted_episodes);
//!     assert!(a.episodes.iter().all(|e| e.peak >= 1.0 && e.duration >= 0.0));
//!     assert!(a.integral >= 0.0 && a.integral <= 20.0 * a.maximum);
//!     // The log keeps the first five events and counts the rest.
//!     let events = result.replicates[observation.id as usize].events;
//!     assert_eq!(observation.log.len() as u64 + observation.omitted_log_entries, events);
//! }
//! ```
use crate::InstanceRecord;
use crate::expr::{Env, Expr, Predicate, Scope};
use crate::model::{Compiled, Model};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// What the engine does after an [`EventObserver`] saw an event.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObserverAction {
    /// Keep simulating the replicate.
    Continue,
    /// End the replicate now, as [`Termination::Stopped`](crate::Termination::Stopped)
    /// at the event's time. Its state is held to the end of the report grid
    /// and it counts as complete for statistics. Other replicates continue.
    Stop,
}

/// A committed state change. For a leap, rule is None and multiplicities
/// identify the aggregate firings. Reactants and products describe exact steps.
///
/// Indices refer to the result's lists: types to
/// [`RulesResult::types`](crate::RulesResult::types), rules to
/// [`RulesResult::rules`](crate::RulesResult::rules) and observable values
/// to [`RulesResult::observables`](crate::RulesResult::observables). The
/// borrowed slices are valid only during the call; copy what you keep, for
/// example with [`log_entry`].
pub struct Event<'a> {
    /// Replicate index.
    pub replicate: u64,
    /// Simulated time of the event; for a leap, the end of the leap.
    pub time: f64,
    /// Index of the rule that fired in an exact step; `None` for a leap.
    pub rule: Option<usize>,
    /// The consumed instances of an exact step as (type index, attribute
    /// values), one per reactant pattern in rule order. Empty for a leap.
    pub reactants: &'a [(usize, &'a [f64])],
    /// The produced instances of an exact step as (type index, attribute
    /// values), one per product in rule order. A product that copies a
    /// reactant appears in both lists. Empty for a leap.
    pub products: &'a [(usize, &'a [f64])],
    /// Every observable after the event, in result column order.
    pub observables: &'a [f64],
    /// For a leap, (rule index, firings) for each lowered channel that fired
    /// at least once; a rule with several channels can appear more than once.
    /// Empty for an exact step.
    pub firings: &'a [(usize, u64)],
}

/// Receives committed events during [`crate::execute_with_observer`].
///
/// `on_event` runs on the native worker threads, concurrently for different
/// replicates, so implementations must be [`Sync`] and should use atomics or
/// locks for shared state. It is called once per event, after the state
/// update and tracker observation and before the model's stop condition.
/// Keep it cheap: it sits in the event loop. It must not panic; see the
/// panics section of [`crate::execute_with_observer`].
pub trait EventObserver: Sync {
    /// An observer that sets this to false is never called, and the engine
    /// compiles the per-event observation away.
    const ENABLED: bool = true;
    /// Called for every committed event; return
    /// [`ObserverAction::Stop`] to end this replicate.
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

/// A tracker in a run request's `trackers` (rules-run/0.2 and later).
///
/// It follows one observable over each replicate, from time 0 to the
/// replicate's stop time or horizon, and reports a tracker result per
/// replicate. A positive value is "present": crossings of zero are counted,
/// and every present span is an episode.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TrackerSpec {
    /// Name of a model observable; an unknown name is `CRN-R020` at
    /// `/trackers`. Neither the observable nor the predicate may read the
    /// continuous `ATP` column of an ATP model.
    pub observable: String,
    /// Optional predicate over parameters and observables, with the syntax
    /// of a model's `stop`. The first event time at which it holds (0 if it
    /// holds initially) is reported; a NaN value does not count as a hit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub first_passage: Option<String>,
    /// Episodes kept in full, in time order; later episodes are only
    /// summarized. Default 0.
    #[serde(default)]
    pub keep_episodes: usize,
}

/// The event log in a run request's `event_log` (rules-run/0.2 and later).
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LogSpec {
    /// Entries kept per replicate; later matching events are only counted.
    /// Replicates times `max_entries` times a per-entry bound from the
    /// model's sizes (observables, attributes and the largest reactant plus
    /// product count of a rule) must fit `budget.max_output_bytes`; the
    /// `CRN-R020` message names the bytes needed.
    pub max_entries: usize,
    /// An empty list logs all rules. Leaps match any listed firing.
    /// Unknown rule names are `CRN-R020` at `/event_log`.
    #[serde(default)]
    pub rules: Vec<String>,
}

/// A span during which a tracked observable was positive.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Episode {
    /// Time the observable became positive; 0 if it was positive initially.
    pub start: f64,
    /// Open at the end of observation when absent.
    pub end: Option<f64>,
    /// Largest value during the episode.
    pub peak: f64,
    /// Observed duration, including the observed part of an open episode.
    pub duration: f64,
}

/// The summary of one tracker over one replicate.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TrackerResult {
    /// The tracked observable.
    pub observable: String,
    /// Changes from zero or below to a positive value. A value positive at
    /// time 0 opens an episode but is not an appearance.
    pub appearances: u64,
    /// Changes from a positive value to zero or below.
    pub extinctions: u64,
    /// Time of the last appearance, if any.
    pub last_appearance: Option<f64>,
    /// Time of the last extinction, if any.
    pub last_extinction: Option<f64>,
    /// First time the `first_passage` predicate held; absent without a
    /// predicate or if it never held.
    pub first_passage: Option<f64>,
    /// Smallest value, over the initial value and every post-event value.
    pub minimum: f64,
    /// Largest value, over the initial value and every post-event value.
    pub maximum: f64,
    /// Time integral of the observable from 0 to the end of observation
    /// (the replicate's stop time or horizon). A non-finite integral fails
    /// the replicate with `CRN-R030`.
    pub integral: f64,
    /// The first `keep_episodes` episodes, in time order. An episode still
    /// open at the end is included if there is room.
    pub episodes: Vec<Episode>,
    /// Episodes beyond `keep_episodes`.
    pub omitted_episodes: u64,
    /// Total duration of the omitted episodes.
    pub omitted_duration: f64,
    /// Largest peak among the omitted episodes; absent when none was
    /// omitted.
    pub omitted_peak: Option<f64>,
}

/// One retained event of the event log.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LogEntry {
    /// Simulated time of the event; for a leap, the end of the leap.
    pub time: f64,
    /// Index into the result's `rules` for an exact step; absent for a leap.
    pub rule: Option<usize>,
    /// Consumed instances of an exact step, one per reactant pattern, each
    /// with count 1. Empty for a leap.
    pub reactants: Vec<InstanceRecord>,
    /// Produced instances of an exact step, one per product, each with count
    /// 1; a product that copies a reactant appears in both lists. Empty for
    /// a leap.
    pub products: Vec<InstanceRecord>,
    /// Every observable after the event, in the result's column order.
    pub observables: Vec<f64>,
    /// For a leap, `[rule index, firings]` pairs, one per lowered channel
    /// that fired; empty for an exact step.
    pub firings: Vec<(usize, u64)>,
}

/// The trackers and event log of one replicate, in the result's
/// `observations`.
#[derive(Debug, Default, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Observation {
    /// Replicate index.
    pub id: u64,
    /// One result per requested tracker, in request order.
    pub trackers: Vec<TrackerResult>,
    /// Retained events, in event order.
    pub log: Vec<LogEntry>,
    /// Matching events beyond `max_entries`.
    pub omitted_log_entries: u64,
}

pub(crate) struct Spec {
    trackers: Vec<(TrackerSpec, usize, Option<Predicate>)>,
    log: Option<(usize, Vec<usize>)>,
    sparse: Vec<bool>,
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
                            .map(|e| Predicate::compile(&e, compiled.observables.len()))
                    })
                    .transpose()?;
                Ok::<_, String>((t.clone(), index, predicate))
            })
            .collect::<Result<_, String>>()
            .map_err(|e| ("/trackers", e))?;
        let continuous: Vec<(usize, &str)> = compiled
            .atp
            .iter()
            .map(|atp| (atp.index, "ATP"))
            .chain(
                compiled
                    .continuous
                    .iter()
                    .map(|v| (v.index, v.name.as_str())),
            )
            .collect();
        for (_, index, predicate) in &trackers {
            let mut flags = vec![false; compiled.observables.len()];
            if let Some(expr) = predicate {
                expr.reads(&mut flags);
            }
            if let Some((_, name)) = continuous.iter().find(|(i, _)| *index == *i || flags[*i]) {
                return Err((
                    "/trackers",
                    format!("trackers apply to event-driven observables; {name} is continuous"),
                ));
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
        Ok(Self {
            trackers,
            log,
            sparse: Vec::new(),
        })
    }

    /// Whether any tracker integrates between class changes.
    pub(crate) fn sparse_effective(&self) -> bool {
        self.sparse.iter().any(|&v| v)
    }

    pub(crate) fn enable_sparse(&mut self, net: &crate::flat::Network) {
        self.sparse = self
            .trackers
            .iter()
            .map(|(_, index, predicate)| {
                predicate.is_none()
                    && net.observables.iter().any(|(i, reduce, terms)| {
                        *i == *index
                            && matches!(reduce, crate::model::Reduce::Sum)
                            && matches!(terms.as_slice(), [(_, 1.0)])
                    })
            })
            .collect();
    }

    pub(crate) fn sparse_observables(&self) -> Option<impl Iterator<Item = usize> + '_> {
        (self.log.is_none() && !self.sparse.is_empty() && self.sparse.iter().all(|&v| v))
            .then(|| self.trackers.iter().map(|(_, i, _)| *i))
    }

    pub(crate) fn reads(&self, flags: &mut [bool]) {
        if self.log.is_some() {
            flags.fill(true);
        }
        for (_, index, predicate) in &self.trackers {
            flags[*index] = true;
            if let Some(p) = predicate {
                p.reads(flags);
            }
        }
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

    /// The bytes reserved for one replicate's trackers and event log.
    pub(crate) fn reservation(&self, sizes: &Sizes) -> Reservation {
        let trackers = self
            .trackers
            .iter()
            .map(|(t, _, _)| 256 + t.keep_episodes as u128 * 64)
            .sum::<u128>();
        Reservation {
            trackers,
            entries: self.log.as_ref().map_or(0, |(n, _)| *n as u128),
            entry: sizes.entry_bytes(),
        }
    }
}

/// The model and method sizes that bound one retained log entry.
pub(crate) struct Sizes {
    /// Largest attribute count of a type.
    pub attributes: usize,
    /// Observable columns of the result, including continuous ones.
    pub observables: usize,
    /// Rules of the model.
    pub rules: usize,
    /// Largest number of reactant patterns plus products of one rule: the
    /// participants of an exact step.
    pub participants: usize,
    /// Lowered channels for the leap methods, whose entries list firings
    /// instead of participants; `None` for exact steps.
    pub leap_channels: Option<usize>,
}

/// Longest compact JSON of a finite f64 (`-2.2250738585072014e-308`) or a
/// u64, plus a separating comma.
const NUMBER: u128 = 25;

impl Sizes {
    /// An upper bound on one log entry, both as compact JSON and in memory.
    ///
    /// The model-based bound counts every number at its longest JSON form:
    /// 160 bytes for the keys, time and rule; 25 per observable; per
    /// participant 62 bytes plus 25 per attribute (type index, attributes and
    /// count 1); per leap firing 45 bytes. Up to 0.2.0 the bound assumed 2
    /// reactants and 64 products for every rule. The smaller of the two is
    /// used, so a request admitted before stays admitted; both exceed an
    /// entry's in-memory size.
    fn entry_bytes(&self) -> u128 {
        let attributes = self.attributes as u128;
        let legacy = 160
            + self.observables as u128 * 8
            + self.rules as u128 * 16
            + (2 + crate::model::MAX_PRODUCTS as u128) * (48 + 8 * attributes);
        let body = match self.leap_channels {
            Some(channels) => channels as u128 * 45,
            None => self.participants as u128 * (62 + NUMBER * attributes),
        };
        legacy.min(160 + self.observables as u128 * NUMBER + body)
    }
}

/// The observation bytes reserved per replicate: trackers plus
/// `entries × entry` for the event log.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Reservation {
    pub trackers: u128,
    pub entries: u128,
    pub entry: u128,
}

impl Reservation {
    pub(crate) fn per_replicate(&self) -> u128 {
        self.trackers + self.entries * self.entry
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
        for (i, (tracker, (spec, index, predicate))) in self
            .trackers
            .iter_mut()
            .zip(&self.spec.trackers)
            .enumerate()
        {
            let value = event.observables[*index];
            if self.spec.sparse.get(i) == Some(&true) && value == tracker.value {
                continue;
            }
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

    /// The path's observation, or none if it failed before tracking began
    /// (initial state rejected or cancelled), so no partial list is emitted.
    pub(crate) fn finish(mut self, time: f64) -> Option<Observation> {
        if self.trackers.len() < self.spec.trackers.len() {
            return None;
        }
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
        Some(self.output)
    }
}

fn truth(value: f64) -> bool {
    !value.is_nan() && value != 0.0
}

/// Copies an event into the owned form the event log uses: participants as
/// [`InstanceRecord`]s with count 1. Observers use it to keep events beyond
/// the callback, for example to serialize them.
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
