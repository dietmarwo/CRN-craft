//! Ahead-of-time kernels for a lowered Static Direct model.
//!
//! [`generate_rules_kernel`] emits the rate and projection arithmetic of one
//! prepared model as Rust source; [`execute_with_rules_kernel`] runs a
//! compiled kernel under the Static Direct scheduler. See the
//! [rules kernel checklist](https://github.com/dietmarwo/CRN-craft/blob/main/docs/developer/specialized-kernels.md#ahead-of-time-generated-rules-kernels).
use crate::rules::expr::RateKernel;
use crate::rules::flat::Network;
use crate::rules::model::Reduce;
use crate::rules::{Control, Diagnostic, Method, Prepared, RulesResult, request_error};
use std::collections::BTreeMap;
use std::fmt::Write;

/// Compiled rate and projection code for one lowered rules model, as emitted
/// by [`generate_rules_kernel`]. The Static Direct scheduler keeps state
/// changes, clocks, sampling, budgets, reports and observation.
///
/// The indices are internal: channels and projections are numbered in the
/// lowered network's order, and `values` holds the model's observables in
/// report order. Implement the trait only through the generator; a
/// handwritten implementation must reproduce the interpreter bit for bit.
pub trait RulesKernel: Sync {
    /// SHA-256 of the emitted source with this value left out, keyed by the
    /// model bytes. A kernel runs only if it equals the current generator's
    /// output for the prepared model.
    fn fingerprint(&self) -> &str;
    /// Recompute all rates when `previous` is `None`, otherwise only the
    /// channels that the last fired channel affects. Return the first invalid
    /// channel in declaration order.
    fn rates(
        &self,
        previous: Option<usize>,
        counts: &[u64],
        values: &[f64],
        rates: &mut [f64],
    ) -> Result<(), usize>;
    /// The value of lowered projection `index` for these counts.
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

const PLACEHOLDER: &str = "@FINGERPRINT@";

/// Covers the whole emitted source, so a kernel left stale by a change in
/// lowering, constants, projections or dependency rows is refused.
fn fingerprint(prepared: &Prepared, template: &str) -> String {
    crate::sha256(format!("rules-kernel/2\n{}\n{template}", prepared.model_sha256).as_bytes())
}

/// Run a prepared Static Direct request with a generated kernel. Results are
/// byte-identical to [`crate::rules::execute`]. Fails with `CRN-R020` at `/model`
/// unless the kernel was generated for these exact model bytes by this
/// version of the generator. Observers are not supported on this path.
///
/// The fingerprint check regenerates the source, so it also refuses a kernel
/// after parameter rebinding or a generator change. Waiting-time and tracker
/// options of the request apply as usual.
///
/// # Errors
///
/// `CRN-R020` at `/method` if the request does not select `static_direct`,
/// `CRN-R020` at `/model` on a fingerprint mismatch, and otherwise the
/// errors of [`crate::rules::execute`].
pub fn execute_with_rules_kernel(
    prepared: &Prepared,
    control: &Control,
    kernel: &impl RulesKernel,
) -> Result<RulesResult, Diagnostic> {
    let net = network(prepared)?;
    crate::rules::linear_selector(prepared, "generated rules kernels")?;
    let expected = template(prepared, net)
        .ok()
        .map(|t| fingerprint(prepared, &t));
    if expected.as_deref() != Some(kernel.fingerprint()) {
        return Err(request_error(
            "/model",
            "generated rules kernel fingerprint mismatch",
            "Generate and rebuild the kernel for these exact model bytes.",
        ));
    }
    crate::rules::execute_inner(
        prepared,
        control,
        &crate::rules::NoObserver,
        Some(kernel),
        None,
    )
}

fn bits(value: f64) -> String {
    format!("f64::from_bits({:#018x})", value.to_bits())
}

/// Emit a Rust module defining `GeneratedRulesKernel` for a prepared Static
/// Direct request. The caller compiles this source into the application and
/// runs it with [`execute_with_rules_kernel`]. Supports constant and vacancy
/// rates and sum observables; parameters are compiled in, so another
/// parameter set needs another kernel.
///
/// A constant rate is an expression over parameters only. A vacancy rate has
/// the form `max(0, K - N) / D * c1 * c2 ...`, with `K` and the factors over
/// parameters and `N`, `D` observables. The example
/// `generate_rules_kernel` in the repository writes the source for a model
/// and request file.
///
/// # Errors
///
/// `CRN-R020` at `/method` if the request does not select `static_direct`,
/// or if a rate or an observable is outside the supported forms, including
/// whole-rule propensities and higher-order rules (`crn-craft.model/1.0`).
///
/// ```
/// use crn_craft::rules::{generate_rules_kernel, prepare};
/// # use serde_json::json;
/// # let model = serde_json::to_vec(&json!({
/// #     "schema": "crn-craft.rules-model/0.1", "name": "immigration_death",
/// #     "parameters": {"k_in": 10.0, "k_out": 1.0},
/// #     "types": [{"name": "X"}], "initial": [],
/// #     "observables": {"X": {"type": "X"}},
/// #     "rules": [
/// #         {"name": "immigration", "reactants": [], "rate": "k_in", "products": [{"type": "X"}]},
/// #         {"name": "death", "reactants": [{"label": "x", "type": "X"}], "rate": "k_out", "products": []}
/// #     ]
/// # })).unwrap();
/// # let request = |method: serde_json::Value| serde_json::to_vec(&json!({
/// #     "schema": "crn-craft.rules-run/0.2",
/// #     "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(&model)},
/// #     "time": {"end": 5.0, "report_times": [0.0, 5.0]},
/// #     "ensemble": {"replicates": 4, "root_seed": 1},
/// #     "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
/// #     "budget": {"max_events_per_replicate": 100000, "max_sampling_steps_per_replicate": 1000000,
/// #                "max_instances": 1, "threads": 1, "max_output_bytes": 1048576},
/// #     "method": method
/// # })).unwrap();
/// let prepared = prepare(&model, &request(json!({"kind": "static_direct"}))).unwrap();
/// let source = generate_rules_kernel(&prepared).unwrap();
/// assert!(source.contains("impl crn_craft::rules::RulesKernel for GeneratedRulesKernel"));
/// // The source is deterministic; its fingerprint ties it to the model bytes.
/// assert_eq!(source, generate_rules_kernel(&prepared).unwrap());
///
/// // Other methods have no lowered Static Direct network to compile.
/// let tree = prepare(&model, &request(json!({"kind": "tree"}))).unwrap();
/// assert_eq!(generate_rules_kernel(&tree).unwrap_err().code, "CRN-R020");
/// ```
pub fn generate_rules_kernel(prepared: &Prepared) -> Result<String, Diagnostic> {
    let net = network(prepared)?;
    let template = template(prepared, net)?;
    let hash = fingerprint(prepared, &template);
    Ok(template.replacen(PLACEHOLDER, &hash, 1))
}

/// The emitted source with a placeholder for its fingerprint.
fn template(prepared: &Prepared, net: &Network) -> Result<String, Diagnostic> {
    let _ = prepared;
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
    if net.channels.iter().any(|c| c.whole) {
        return Err(request_error(
            "/method",
            "generated rules kernels do not support whole-rule propensities",
            "Use interpreted Static Direct for models with propensity rules.",
        ));
    }
    if net.channels.iter().any(|c| c.multi.is_some()) {
        return Err(request_error(
            "/method",
            "generated rules kernels do not support higher-order rules",
            "Use interpreted Static Direct or a generated event loop for rules of more than two reactant molecules.",
        ));
    }
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
            RateKernel::Expression { .. } => return Err(unsupported()),
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
    let hash = PLACEHOLDER;
    Ok(format!(
        r#"// Generated by CRN-craft. Build with the same floating-point policy.
pub struct GeneratedRulesKernel;
impl crn_craft::rules::RulesKernel for GeneratedRulesKernel {{
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

#[cfg(test)]
mod tests {
    use super::*;

    fn fingerprint_of(source: &str) -> &str {
        source
            .split("fn fingerprint(&self) -> &str { \"")
            .nth(1)
            .and_then(|rest| rest.split('"').next())
            .unwrap()
    }

    #[test]
    fn the_fingerprint_covers_the_emitted_arithmetic() {
        // Same model bytes and channel shapes; a projection weight changes as
        // it would after a change in lowering. A stale kernel must not pass.
        let model = include_bytes!("../../../../models/rules/mtdna-three-class.rules.json");
        let request = serde_json::json!({
            "schema": "crn-craft.rules-run/0.2",
            "model": {"path": "m.rules.json", "sha256": crate::sha256(model)},
            "time": {"end": 1.0, "report_times": [0.0, 1.0]},
            "ensemble": {"replicates": 1, "root_seed": 1},
            "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
            "budget": {"threads": 1, "max_events_per_replicate": 10u64,
                       "max_sampling_steps_per_replicate": 100u64, "max_instances": 10u64,
                       "max_output_bytes": 1u64 << 20},
            "method": {"kind": "static_direct"}
        });
        let mut prepared =
            crate::rules::prepare(model, &serde_json::to_vec(&request).unwrap()).unwrap();
        let before = generate_rules_kernel(&prepared).unwrap();
        prepared.flat.as_mut().unwrap().observables[0].2[0].1 *= 2.0;
        let after = generate_rules_kernel(&prepared).unwrap();
        assert_ne!(before, after);
        assert_ne!(fingerprint_of(&before), fingerprint_of(&after));
    }
}
