//! Fused exact Direct loop for the recognized fixed-population topology.
use super::*;
use crate::flat::{Network, check_instances};
use crate::mtdna::Kernel;

impl Path<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn simulate_mtdna<R: RngCore, O: EventObserver>(
        &self,
        kernel: &Kernel,
        net: &Network,
        counts: &mut [u64],
        values: &mut [f64],
        words: &mut Words<'_, R>,
        out: &mut PathOutput,
        time: &mut f64,
        monitor: &mut Option<Monitor<'_>>,
        observer: &O,
    ) -> Result<Option<f64>, Halt> {
        let mut total: u64 = counts.iter().sum();
        let mut live = out.peak_instances;
        let [wt, neutral, advantaged] = kernel.species;
        let [ww, nw, aw] = kernel.weights;
        let [wb, nb, ab] = kernel.birth;
        loop {
            if out.events % POLL_EVENTS == 0 && out.events > 0 {
                self.poll()?;
            }
            if total == 0 {
                return Ok(None);
            }
            let w = counts[wt] as f64;
            let n = counts[neutral] as f64;
            let a = counts[advantaged] as f64;
            let scale = (kernel.capacity - total) as f64 / (w * ww + n * nw + a * aw);
            let death = total as f64 * kernel.delta;
            let wt_scale = scale * w;
            let wt_end = death + wt_scale * wb;
            let neutral_end = wt_end + scale * n * nb;
            let rate = neutral_end + scale * a * ab;
            if rate == 0.0 {
                return Ok(None);
            }
            let wait = crn_craft::ssa::exponential_from_words(|retry| {
                if retry {
                    self.poll()?;
                }
                words.word()
            })? / rate;
            let next = *time + wait;
            if !wait.is_finite() || next <= *time {
                return Err(numerical(
                    "/time",
                    "clock failed to advance",
                    "Rescale rates or time.",
                ));
            }
            // Match Static Direct's sampling-budget behavior at the horizon.
            let pick = words.uniform()? * rate;
            if out.saved_points < self.report_times.len()
                && self.report_times[out.saved_points] < next
            {
                net.refresh(counts, values)?;
                self.save_before(out, values, next);
            }
            if next > self.end {
                return Ok(None);
            }
            let (slot, species, birth) = if pick < death {
                if pick < w * kernel.delta {
                    (0, wt, false)
                } else if a == 0.0 || pick < (w + n) * kernel.delta {
                    (4, neutral, false)
                } else {
                    (6, advantaged, false)
                }
            } else if wt_end > death && (pick < wt_end || wt_end == rate) {
                let within = pick - death;
                // If rounding lands at a group's upper boundary, select its
                // last positive outcome, including probability endpoints.
                if within < wt_scale * kernel.wt_unchanged || kernel.wt_unchanged == wb {
                    (1, wt, true)
                } else if within < wt_scale * kernel.wt_through_neutral
                    || kernel.wt_through_neutral == wb
                {
                    (2, neutral, true)
                } else {
                    (3, advantaged, true)
                }
            } else if neutral_end > wt_end && (pick < neutral_end || neutral_end == rate) {
                (5, neutral, true)
            } else if rate > neutral_end {
                (7, advantaged, true)
            } else if a > 0.0 {
                (6, advantaged, false)
            } else if n > 0.0 {
                (4, neutral, false)
            } else {
                (0, wt, false)
            };
            self.check_events(out, 1)?;
            if birth {
                if counts[species] == 0 {
                    check_instances(live + 1, self.limits)?;
                    live += 1;
                    out.created_instances += 1;
                    out.peak_instances = out.peak_instances.max(live);
                }
                counts[species] += 1;
                total += 1;
            } else {
                counts[species] -= 1;
                total -= 1;
                live -= u64::from(counts[species] == 0);
            }
            out.events += 1;
            let channel = kernel.channels[slot];
            out.rule_events[net.channels[channel].rule] += 1;
            *time = next;
            let mut stopped = false;
            if O::ENABLED || (monitor.is_some() && kernel.count_trackers.is_none()) {
                net.refresh(counts, values)?;
                stopped = self.observe_flat(net, channel, values, next, monitor, observer);
            } else if let (Some(m), Some(indices)) = (&mut *monitor, &kernel.count_trackers) {
                m.observe_count(&indices[species], counts[species] as f64, next);
            }
            if stopped
                || kernel.threshold.is_some_and(|threshold| {
                    total > 0 && (total - counts[wt]) as f64 >= threshold * total as f64
                })
            {
                return Ok(Some(next));
            }
        }
    }
}
