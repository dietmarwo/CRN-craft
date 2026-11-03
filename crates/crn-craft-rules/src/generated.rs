//! Ahead-of-time kernels for a lowered Static Direct model.
use crate::expr::RateKernel;
use crate::flat::Network;
use crate::model::Reduce;
use crate::{Control, Diagnostic, Method, Prepared, RulesResult, request_error};
use std::collections::BTreeMap;
use std::fmt::Write;

/// Trusted, compiled rate and projection code for a single bound rules model.
/// The scheduler retains state changes, clocks, sampling and event accounting.
pub trait RulesKernel: Sync {
    fn fingerprint(&self) -> &str;
    /// Recompute all rates, or just those affected by the last event.
    /// Return the first invalid channel in declaration order.
    fn rates(
        &self,
        previous: Option<usize>,
        counts: &[u64],
        values: &[f64],
        rates: &mut [f64],
    ) -> Result<(), usize>;
    fn projection(&self, index: usize, counts: &[u64]) -> f64;
}

fn network(prepared: &Prepared) -> Result<&Network, Diagnostic> {
    if matches!(prepared.request.method, Some(Method::StaticDirect { .. })) {
        return Ok(prepared.flat.as_ref().expect("lowered Static Direct"));
    }
    Err(request_error(
        "/method",
        "generated rules kernels require Static Direct",
        "Prepare a closed Static Direct model.",
    ))
}

fn fingerprint(prepared: &Prepared, net: &Network) -> String {
    let mut key = format!("rules-kernel/1\n{}\n", prepared.model_sha256);
    for c in &net.channels {
        writeln!(key, "{}:{:?}:{:?}", c.rule, c.reactants, c.changes).unwrap();
    }
    crn_craft::sha256(key.as_bytes())
}

pub fn execute_with_rules_kernel(
    prepared: &Prepared,
    control: &Control,
    kernel: &impl RulesKernel,
) -> Result<RulesResult, Diagnostic> {
    let net = network(prepared)?;
    if kernel.fingerprint() != fingerprint(prepared, net) {
        return Err(request_error(
            "/model",
            "generated rules kernel fingerprint mismatch",
            "Generate and rebuild the kernel for these exact model bytes.",
        ));
    }
    crate::execute_inner(prepared, control, &crate::NoObserver, Some(kernel))
}

fn bits(value: f64) -> String {
    format!("f64::from_bits({:#018x})", value.to_bits())
}

/// Emit a Rust module defining `GeneratedRulesKernel`. The caller compiles
/// this source into the application. The generic solver remains available for
/// other rate and projection families.
pub fn generate_rules_kernel(prepared: &Prepared) -> Result<String, Diagnostic> {
    let net = network(prepared)?;
    let unsupported = || {
        request_error(
            "/method",
            "generated rules kernels support constant or vacancy rates and sum observables",
            "Use interpreted Static Direct for other expressions or reductions.",
        )
    };
    let mut bases = BTreeMap::new();
    let mut definitions = String::new();
    let mut channels = String::new();
    for (i, c) in net.channels.iter().enumerate() {
        let global = match &net.rates[c.rule] {
            RateKernel::Constant(k) => bits(*k),
            RateKernel::Vacancy {
                limit,
                count,
                denominator,
                factors,
            } => {
                let key = (limit.to_bits(), *count, *denominator);
                let next = bases.len();
                let slot = bases.entry(key).or_insert(next);
                if *slot == next {
                    writeln!(definitions, "let base_{slot} = ({} - values[{count}]).max(0.0) / values[{denominator}];", bits(*limit)).unwrap();
                }
                let mut expression = format!("base_{slot}");
                for &factor in factors {
                    write!(expression, " * {}", bits(factor)).unwrap();
                }
                expression
            }
            RateKernel::Expression(_) => return Err(unsupported()),
        };
        let weight = match *c.reactants.as_slice() {
            [] => "1.0".into(),
            [s] => format!("counts[{s}] as f64"),
            [s, t] if s == t => {
                format!("counts[{s}] as f64 * counts[{s}].saturating_sub(1) as f64 * 0.5")
            }
            [s, t] => format!("counts[{s}] as f64 * counts[{t}] as f64"),
            _ => unreachable!("validated arity"),
        };
        writeln!(channels, "{i} => {{ let weight = {weight}; if weight == 0.0 {{ 0.0 }} else {{ ({global}) * (weight * {}) }} }},", bits(c.factor)).unwrap();
    }
    let mut projections = String::new();
    for (index, (_, reduce, terms)) in net.observables.iter().enumerate() {
        if !matches!(reduce, Reduce::Sum) {
            return Err(unsupported());
        }
        write!(
            projections,
            "{index} => {{ let terms: [f64; {}] = [",
            terms.len()
        )
        .unwrap();
        for &(species, weight) in terms {
            write!(projections, "counts[{species}] as f64 * {},", bits(weight)).unwrap();
        }
        writeln!(projections, "]; terms.into_iter().sum() }},").unwrap();
    }
    let mut affected = String::new();
    for row in &net.affected {
        writeln!(affected, "&{row:?},").unwrap();
    }
    let hash = fingerprint(prepared, net);
    Ok(format!(
        r#"// Generated by CRN-craft. Build with the same floating-point policy.
pub struct GeneratedRulesKernel;
impl crn_craft_rules::RulesKernel for GeneratedRulesKernel {{
    fn fingerprint(&self) -> &str {{ "{hash}" }}
    fn rates(&self, previous: Option<usize>, counts: &[u64], values: &[f64],
             rates: &mut [f64]) -> Result<(), usize> {{
        const AFFECTED: &[&[usize]] = &[{affected}];
        let _ = (counts, values);
        {definitions}
        let rate = |i| match i {{
            {channels}
            _ => unreachable!("validated channel"),
        }};
        let mut update = |i| {{
            let a: f64 = rate(i);
            if !a.is_finite() || a < 0.0 {{ return Err(i); }}
            rates[i] = a;
            Ok(())
        }};
        if let Some(previous) = previous {{
            for &i in AFFECTED[previous] {{ update(i)?; }}
        }} else {{
            for i in 0..{} {{ update(i)?; }}
        }}
        Ok(())
    }}
    fn projection(&self, index: usize, counts: &[u64]) -> f64 {{
        let _ = counts;
        match index {{
            {projections}
            _ => unreachable!("validated projection"),
        }}
    }}
}}
"#,
        net.channels.len()
    ))
}
