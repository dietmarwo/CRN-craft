//! ATP coupling with analytic flow between stochastic events.
//!
//! A `crn-craft.rules-model/0.2` model may carry an `atp` block
//! ([`AtpSpec`]). ATP is a continuous quantity A that follows
//!
//! dA/dt = f × wt − v1 × A − v2 × N
//!
//! between events, where `wt` and `N` are two of the model's observables
//! (wild-type molecules and all molecules). Selected replication rules are
//! multiplied by c / (c + A). The `thinning` method integrates the flow
//! analytically and samples the coupled rules exactly with a thinning
//! envelope; it is the only method for such a model, and refused for others.
//! A replicate stops (`stopped`) when ATP reaches zero while the forcing
//! f × wt − v2 × N is not positive; with positive forcing ATP rises again.
//! Results carry ATP as an extra observable column named `ATP`, after the
//! model's observables, and report it as zero after an ATP stop. Trackers
//! cannot read it.
//!
//! The [method reference](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference/methods.md#the-atp-block-and-thinning)
//! describes the fixture `models/rules/mtdna-three-class-atp.rules.json`.
//! rules-model/0.3 can write the same flow as an ordinary
//! [`crate::rules::continuous`] variable
//! (`models/rules/mtdna-three-class-atp-continuous.rules.json`), simulated by
//! the approximate `hybrid` method instead of exact thinning.
use crate::rules::model::{ObservableInfo, RuleInfo};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// The ATP block of a rules-model/0.2 model.
///
/// Every coefficient must be finite and nonnegative, with a positive
/// `scale`; the observable and rule names must exist and `ATP` must not be
/// an observable name. Violations are `CRN-R015` at `/atp`.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AtpSpec {
    /// Name of the observable counting wild-type molecules, wt.
    pub wild_type: String,
    /// Name of the observable counting all molecules, N.
    pub total: String,
    /// ATP at time 0.
    pub initial: f64,
    /// Production per wild type, f.
    pub production: f64,
    /// First-order ATP consumption, v1.
    pub decay: f64,
    /// Maintenance per mtDNA molecule, v2.
    pub maintenance: f64,
    /// Replication multiplier is c / (c + ATP).
    pub scale: f64,
    /// Names of the rules whose rate is multiplied by c / (c + ATP), each at
    /// most once.
    pub replication_rules: Vec<String>,
}

#[derive(Debug)]
pub(crate) struct Atp {
    pub index: usize,
    pub initial: f64,
    wild_type: usize,
    total: usize,
    production: f64,
    decay: f64,
    maintenance: f64,
    scale: f64,
    pub coupled: Vec<bool>,
}

impl Atp {
    pub fn compile(
        spec: &AtpSpec,
        obs: &[ObservableInfo],
        rules: &[RuleInfo],
    ) -> Result<Self, String> {
        if [spec.initial, spec.production, spec.decay, spec.maintenance]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
            || !spec.scale.is_finite()
            || spec.scale <= 0.0
        {
            return Err(
                "ATP coefficients must be finite and nonnegative, with positive scale".into(),
            );
        }
        if obs.iter().any(|o| o.name == "ATP") {
            return Err("ATP is a reserved output name".into());
        }
        let index = |name: &str| {
            obs.iter()
                .position(|o| o.name == name)
                .ok_or_else(|| format!("unknown ATP input observable {name}"))
        };
        let mut coupled = vec![false; rules.len()];
        for name in &spec.replication_rules {
            let i = rules
                .iter()
                .position(|r| r.name == *name)
                .ok_or_else(|| format!("unknown ATP replication rule {name}"))?;
            if std::mem::replace(&mut coupled[i], true) {
                return Err(format!("duplicate ATP rule {name}"));
            }
        }
        Ok(Self {
            index: obs.len(),
            initial: spec.initial,
            wild_type: index(&spec.wild_type)?,
            total: index(&spec.total)?,
            production: spec.production,
            decay: spec.decay,
            maintenance: spec.maintenance,
            scale: spec.scale,
            coupled,
        })
    }

    pub fn flow(&self, values: &[f64]) -> Flow {
        Flow {
            atp: values[self.index],
            source: self.production * values[self.wild_type]
                - self.maintenance * values[self.total],
            decay: self.decay,
        }
    }

    /// Bound the ATP change to half of the current denominator. The
    /// magnitude of the flow derivative decreases on a fixed-count interval.
    pub fn envelope_span(&self, flow: &Flow) -> f64 {
        let slope = (flow.source - flow.decay * flow.atp).abs();
        if slope > 0.0 {
            (0.5 * flow.atp + 0.5 * self.scale) / slope
        } else {
            f64::INFINITY
        }
    }

    pub fn multiplier(&self, atp: f64) -> f64 {
        1.0 / (1.0 + atp / self.scale)
    }
}

pub(crate) struct Flow {
    pub atp: f64,
    pub source: f64,
    pub decay: f64,
}

impl Flow {
    pub fn at(&self, dt: f64) -> f64 {
        if self.decay == 0.0 {
            self.atp + self.source * dt
        } else {
            let fraction = -libm::expm1(-self.decay * dt);
            let next =
                self.atp * libm::exp(-self.decay * dt) + self.source * (fraction / self.decay);
            if next == 0.0 && self.source >= 0.0 && self.atp > 0.0 {
                f64::from_bits(1)
            } else {
                next
            }
        }
    }

    pub fn zero_time(&self) -> f64 {
        if self.atp <= 0.0 {
            // At zero, the derivative is the source. Positive forcing lifts
            // ATP; otherwise it stays depleted.
            return if self.source > 0.0 {
                f64::INFINITY
            } else {
                0.0
            };
        }
        if self.source >= 0.0 {
            return f64::INFINITY;
        }
        let linear = self.atp / -self.source;
        if self.decay == 0.0 {
            return linear;
        }
        let scaled = linear * self.decay;
        if scaled == 0.0 {
            linear
        } else if scaled.is_finite() {
            libm::log1p(scaled) / self.decay
        } else {
            // log1p(x) = ln(x) to rounding for a ratio beyond f64::MAX.
            (libm::log(self.atp) + libm::log(self.decay) - libm::log(-self.source)) / self.decay
        }
    }
}
