//! ATP coupling with analytic flow between stochastic events.
use crate::model::{ObservableInfo, RuleInfo};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AtpSpec {
    pub wild_type: String,
    pub total: String,
    pub initial: f64,
    /// Production per wild type, f.
    pub production: f64,
    /// First-order ATP consumption, v1.
    pub decay: f64,
    /// Maintenance per mtDNA molecule, v2.
    pub maintenance: f64,
    /// Replication multiplier is c / (c + ATP).
    pub scale: f64,
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
            return 0.0;
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
        } else {
            libm::log1p(scaled) / self.decay
        }
    }
}
