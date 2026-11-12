//! Exact thinning for the monotone ATP flow on each fixed-count interval.
use super::*;
use crate::rules::atp::{Atp, Flow};
use crate::rules::flat::Network;

impl Path<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn simulate_thinning<R: RngCore, O: EventObserver>(
        &self,
        atp: &Atp,
        net: &Network,
        counts: &mut [u64],
        values: &mut [f64],
        words: &mut Words<'_, R>,
        out: &mut PathOutput,
        time: &mut f64,
        monitor: &mut Option<Monitor<'_>>,
        observer: &O,
    ) -> Result<Option<f64>, Halt> {
        let mut live = out.peak_instances;
        let mut cache = net.rate_cache();
        let mut rates = vec![0.0; net.channels.len()];
        loop {
            self.poll()?;
            let flow = atp.flow(values);
            if !flow.source.is_finite() {
                return Err(atp_error("nonfinite ATP forcing"));
            }
            let zero = flow.zero_time();
            let remaining = self.end - *time;
            let uncapped = remaining.min(zero);
            // The span cap only tightens the envelope: on a monotone interval
            // the endpoint minimum bounds the intensity for any length. Drop
            // it when it is below the clock resolution at this time.
            let mut limit = uncapped.min(atp.envelope_span(&flow));
            if limit > 0.0 && *time + limit <= *time {
                limit = uncapped;
            }
            if limit <= 0.0 {
                if zero == 0.0 {
                    values[atp.index] = 0.0;
                    return Ok(Some(*time));
                }
                if remaining <= 0.0 {
                    return Ok(None);
                }
                return Err(atp_error("ATP envelope interval cannot advance"));
            }
            let endpoint = flow.at(limit);
            if !endpoint.is_finite() {
                return Err(atp_error("nonfinite ATP endpoint"));
            }
            let minimum = if limit == zero {
                0.0
            } else {
                flow.atp.min(endpoint).max(0.0)
            };
            let upper = atp.multiplier(minimum);
            let env = Env {
                params: &self.model.params,
                observables: values,
                attributes: [&[], &[]],
            };
            cache.reset();
            let mut bound = 0.0;
            for (i, a) in rates.iter_mut().enumerate() {
                *a = cache.channel(net, i, counts, &env)?;
                bound += *a
                    * if atp.coupled[net.channels[i].rule] {
                        upper
                    } else {
                        1.0
                    };
            }
            bound *= 1.0 + 16.0 * f64::EPSILON;
            if !bound.is_finite() {
                return Err(atp_error("nonfinite thinning intensity"));
            }
            let dt = if bound > 0.0 {
                -libm::log(words.open()?) / bound
            } else {
                f64::INFINITY
            };
            if dt >= limit {
                let until = if limit == remaining {
                    self.end
                } else {
                    (*time + limit).min(self.end)
                };
                if until <= *time && limit != zero {
                    return Err(atp_error("ATP envelope interval cannot advance"));
                }
                self.save_flow(out, values, atp.index, &flow, *time, until)?;
                values[atp.index] = if limit == zero {
                    0.0
                } else {
                    endpoint.max(0.0)
                };
                *time = until;
                if limit == zero {
                    return Ok(Some(until));
                }
                if until == self.end {
                    return Ok(None);
                }
                continue;
            }
            let next = *time + dt;
            if next <= *time {
                return Err(atp_error(format!(
                    "ATP proposal clock cannot advance (dt={dt}, time={time})"
                )));
            }
            self.save_flow(out, values, atp.index, &flow, *time, next)?;
            let level = flow.at(dt);
            if !level.is_finite() {
                return Err(atp_error("nonfinite ATP at proposal"));
            }
            values[atp.index] = level.max(0.0);
            *time = next;
            let multiplier = atp.multiplier(values[atp.index]);
            let mut total = 0.0;
            let target = words.uniform()? * bound;
            let mut selected = None;
            for (i, &a) in rates.iter().enumerate() {
                let increment = if atp.coupled[net.channels[i].rule] {
                    a * multiplier
                } else {
                    a
                };
                total += increment;
                if selected.is_none() && target < total {
                    selected = Some(i);
                }
            }
            if total > bound {
                return Err(atp_error(format!(
                    "ATP intensity {total} exceeds envelope {bound}"
                )));
            }
            if let Some(i) = selected {
                self.apply_flat(net, i, counts, &mut live, out)?;
                net.refresh_channel(i, counts, values, true)?;
                if self.observe_flat(net, i, values, next, monitor, observer)
                    || self.stop_values(values)?
                {
                    return Ok(Some(next));
                }
            }
        }
    }

    fn save_flow(
        &self,
        out: &mut PathOutput,
        values: &[f64],
        index: usize,
        flow: &Flow,
        start: f64,
        until: f64,
    ) -> Result<(), Halt> {
        let n = values.len();
        while out.saved_points < self.report_times.len()
            && self.report_times[out.saved_points] < until
        {
            let row = out.saved_points * n;
            out.values[row..row + n].copy_from_slice(values);
            let value = flow.at(self.report_times[out.saved_points] - start);
            if !value.is_finite() {
                return Err(atp_error("nonfinite ATP at report time"));
            }
            out.values[row + index] = value.max(0.0);
            out.saved_points += 1;
        }
        Ok(())
    }
}

fn atp_error(message: impl Into<String>) -> Halt {
    numerical("/atp", message, "Rescale ATP coefficients, rates and time.")
}
