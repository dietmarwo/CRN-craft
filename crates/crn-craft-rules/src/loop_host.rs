//! Run-time services for a generated rules event loop.
//!
//! A [`crate::RulesLoop`] emitted by [`crate::generate_rules_loop`] owns the
//! arithmetic of every event: rates, selection, state updates and
//! projections. [`LoopHost`] supplies what depends on the request or on run
//! control, in the order interpreted Static Direct uses it: polling, waiting
//! times and uniforms, reports and the horizon, budgets and counters, tracker
//! observation and the stop predicate.
use super::*;
use crate::NoObserver;
use crate::flat::{Network, check_instances};
use rand_chacha::ChaCha8Rng;

/// Why a generated loop stopped early. Opaque: return it unchanged.
pub struct LoopHalt(pub(crate) Halt);

/// The interpreter's services to a generated loop for one replicate. Its
/// methods are an implementation detail of the generator's output and may
/// change with it; the kernel fingerprint ties the two together.
pub struct LoopHost<'h, 'm> {
    path: &'h Path<'h>,
    net: &'h Network,
    words: Words<'h, ChaCha8Rng>,
    out: &'h mut PathOutput,
    time: &'h mut f64,
    monitor: &'h mut Option<Monitor<'m>>,
    live: u64,
    ziggurat: bool,
}

impl LoopHost<'_, '_> {
    /// Polls for cancellation and the wall-time limit every 1024 events.
    #[inline]
    pub fn tick(&mut self) -> Result<(), LoopHalt> {
        if self.out.events % POLL_EVENTS == 0 && self.out.events > 0 {
            self.path.poll().map_err(LoopHalt)?;
        }
        Ok(())
    }

    /// The next event time for `total`, with the request's waiting-time
    /// policy, or the clock error if it does not advance.
    #[inline]
    pub fn next_time(&mut self, total: f64) -> Result<f64, LoopHalt> {
        let time = *self.time;
        let next = if self.ziggurat {
            let path = self.path;
            let words = &mut self.words;
            let unit = crn_craft::ssa::exponential_from_words(|retry| {
                if retry {
                    path.poll()?;
                }
                words.word()
            })
            .map_err(LoopHalt)?;
            time + unit / total
        } else {
            // Open uniform: u = 0 would give a zero wait.
            time - libm::log(1.0 - self.words.open().map_err(LoopHalt)?) / total
        };
        if next <= time {
            return Err(LoopHalt(numerical(
                "/time",
                "clock failed to advance",
                "Rescale rates and time.",
            )));
        }
        Ok(next)
    }

    #[inline]
    pub fn uniform(&mut self) -> Result<f64, LoopHalt> {
        self.words.uniform().map_err(LoopHalt)
    }

    /// The selection error for a nonfinite total or an unresolvable rate.
    #[cold]
    pub fn resolution_error(&self, total: f64) -> LoopHalt {
        LoopHalt(if total.is_finite() {
            Halt::Numerical(Diagnostic::error(
                "CRN-R031",
                "/rules",
                "a positive channel is below 53-bit resolution of the total",
                "Use the tree method or rescale rates.",
            ))
        } else {
            numerical(
                "/rules",
                "unresolvable total propensity",
                "Rescale rates or populations.",
            )
        })
    }

    #[cold]
    pub fn rate_error(&self, channel: usize) -> LoopHalt {
        LoopHalt(self.net.rate_error(channel))
    }

    #[cold]
    pub fn count_error(&self) -> LoopHalt {
        LoopHalt(numerical(
            "/rules",
            "fixed count outside 0..2^53",
            "Reduce population.",
        ))
    }

    #[cold]
    pub fn overflow_error(&self) -> LoopHalt {
        LoopHalt(numerical(
            "/observables",
            "observable overflows",
            "Rescale weights.",
        ))
    }

    /// No channel can fire: the path ends at the horizon.
    #[cold]
    pub fn exhausted(
        &mut self,
        counts: &[u64],
        values: &mut [f64],
    ) -> Result<Option<f64>, LoopHalt> {
        self.net.refresh_cold(counts, values).map_err(LoopHalt)?;
        Ok(None)
    }

    /// Saves the reports due before `next` and returns false, after the
    /// final refresh, when `next` lies beyond the horizon.
    #[inline]
    pub fn reach(
        &mut self,
        next: f64,
        counts: &[u64],
        values: &mut [f64],
    ) -> Result<bool, LoopHalt> {
        let path = self.path;
        if self.out.saved_points < path.report_times.len()
            && path.report_times[self.out.saved_points] < next
        {
            self.net.refresh_cold(counts, values).map_err(LoopHalt)?;
            path.save_before(self.out, values, next);
        }
        if next > path.end {
            self.net.refresh_cold(counts, values).map_err(LoopHalt)?;
            return Ok(false);
        }
        Ok(true)
    }

    /// Admits one more event under the event limit.
    #[inline]
    pub fn begin_event(&self) -> Result<(), LoopHalt> {
        self.path.check_events(self.out, 1).map_err(LoopHalt)
    }

    /// Distinct species alive before the current event.
    #[inline]
    pub fn live(&self) -> u64 {
        self.live
    }

    /// Admits `live` distinct species under the instance limit.
    #[inline]
    pub fn admit(&self, live: u64) -> Result<(), LoopHalt> {
        check_instances(live, self.path.limits).map_err(LoopHalt)
    }

    /// Records that `channel` fired at `next`, after its state update.
    #[inline]
    pub fn fired(&mut self, channel: usize, created: u64, live: u64, next: f64) {
        self.live = live;
        let out = &mut *self.out;
        out.events += 1;
        out.rule_events[self.net.channels[channel].rule] += 1;
        out.created_instances += created;
        out.peak_instances = out.peak_instances.max(live);
        *self.time = next;
    }

    /// Whether trackers observe events of `channel`.
    #[inline]
    pub fn observes(&self, channel: usize) -> bool {
        self.net.observe_affected[channel]
    }

    /// Passes the event to the trackers. Generated loops have no observer,
    /// so this never stops the path.
    #[inline]
    pub fn observe(&mut self, channel: usize, values: &[f64], time: f64) {
        self.path
            .observe_flat(self.net, channel, values, time, self.monitor, &NoObserver);
    }

    /// Finishes a path whose emitted stop predicate fired.
    #[cold]
    pub fn stopped(&mut self, counts: &[u64], values: &mut [f64]) -> Result<(), LoopHalt> {
        self.net.refresh_cold(counts, values).map_err(LoopHalt)
    }

    /// Evaluates the stop predicate as interpreted Static Direct does.
    #[inline]
    pub fn stop(&mut self, counts: &[u64], values: &mut [f64]) -> Result<bool, LoopHalt> {
        match self.path.stop_values(values) {
            Ok(false) => Ok(false),
            Ok(true) => {
                self.net.refresh_cold(counts, values).map_err(LoopHalt)?;
                Ok(true)
            }
            Err(_) => {
                // The diagnostic includes all projection values.
                self.net.refresh_cold(counts, values).map_err(LoopHalt)?;
                self.path.stop_values(values).map_err(LoopHalt)?;
                Ok(false)
            }
        }
    }
}

impl Path<'_> {
    /// Runs a generated loop with the interpreter's random stream and
    /// accounting.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn simulate_generated_loop<R: RngCore + 'static>(
        &self,
        kernel: &dyn crate::RulesLoop,
        net: &Network,
        counts: &mut [u64],
        values: &mut [f64],
        words: &mut Words<'_, R>,
        out: &mut PathOutput,
        time: &mut f64,
        monitor: &mut Option<Monitor<'_>>,
    ) -> Result<Option<f64>, Halt> {
        let rng = (&mut *words.rng as &mut dyn std::any::Any)
            .downcast_mut::<ChaCha8Rng>()
            .expect("rules replicates use the ChaCha8 stream");
        let live = out.peak_instances;
        let mut host = LoopHost {
            path: self,
            net,
            words: Words {
                rng,
                used: words.used,
                max: words.max,
            },
            out,
            time,
            monitor,
            live,
            ziggurat: matches!(
                self.method,
                crate::Method::StaticDirect {
                    waiting_time: Some(crate::RulesWaitingTime::Ziggurat),
                    ..
                }
            ),
        };
        let result = kernel.run(&mut host, counts, values);
        let used = host.words.used;
        words.used = used;
        result.map_err(|halt| halt.0)
    }
}
