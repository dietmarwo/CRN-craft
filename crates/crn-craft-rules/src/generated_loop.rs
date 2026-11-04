//! Ahead-of-time event loops for a lowered Static Direct model.
use crate::expr::RateKernel;
use crate::flat::Network;
use crate::model::{MAX_COUNT, Reduce};
use crate::{
    Control, Diagnostic, LoopHalt, LoopHost, Method, Prepared, RulesResult, request_error,
};
use std::collections::BTreeMap;
use std::fmt::Write;

/// The whole Static Direct event for one lowered rules model, as emitted by
/// [`generate_rules_loop`]: rates, dependency updates, the total and its
/// resolution check, selection, the state update and projection refresh.
/// [`LoopHost`] supplies waiting times, uniforms, reports, budgets, tracker
/// observation and the stop predicate.
///
/// Implement the trait only through the generator; a handwritten
/// implementation must reproduce the interpreter bit for bit.
pub trait RulesLoop: Sync {
    /// SHA-256 of the emitted source with this value left out, keyed by the
    /// model bytes. A loop runs only if it equals the current generator's
    /// output for the prepared model.
    fn fingerprint(&self) -> &str;
    /// Simulates one replicate from the prepared initial state. Returns the
    /// stop time if the model's stop predicate fired.
    fn run(
        &self,
        host: &mut LoopHost<'_, '_>,
        counts: &mut [u64],
        values: &mut [f64],
    ) -> Result<Option<f64>, LoopHalt>;
}

const PLACEHOLDER: &str = "@FINGERPRINT@";

/// Rate code is emitted once per dependency row, so the source grows with
/// the square of the channel count in densely coupled networks. Beyond this
/// size the loop is refused before it can exhaust the compiler; the rates
/// kernel and the interpreter still apply.
const MAX_SOURCE: usize = 1 << 20;

fn network(prepared: &Prepared) -> Result<&Network, Diagnostic> {
    if matches!(prepared.request.method, Some(Method::StaticDirect { .. })) {
        return Ok(prepared.flat.as_ref().expect("lowered Static Direct"));
    }
    Err(request_error(
        "/method",
        "generated rules loops require Static Direct",
        "Prepare a closed Static Direct model.",
    ))
}

/// Covers the emitted source. Bump the version when `LoopHost` changes what
/// a host call does without changing the emitted text.
fn fingerprint(prepared: &Prepared, template: &str) -> String {
    crn_craft::sha256(format!("rules-loop/1\n{}\n{template}", prepared.model_sha256).as_bytes())
}

/// Run a prepared Static Direct request with a generated event loop. Results
/// are byte-identical to [`crate::execute`] for the same request. Fails with
/// `CRN-R020` at `/model` unless the loop was generated for these exact model
/// bytes by this version of the generator. Observers are not supported on
/// this path.
pub fn execute_with_rules_loop(
    prepared: &Prepared,
    control: &Control,
    kernel: &impl RulesLoop,
) -> Result<RulesResult, Diagnostic> {
    let net = network(prepared)?;
    let expected = fingerprint(prepared, &template(prepared, net)?);
    if expected != kernel.fingerprint() {
        return Err(request_error(
            "/model",
            "generated rules loop fingerprint mismatch",
            "Generate and rebuild the loop for these exact model bytes.",
        ));
    }
    crate::execute_inner(prepared, control, &crate::NoObserver, None, Some(kernel))
}

/// Emit a Rust module defining `GeneratedRulesLoop` for a prepared Static
/// Direct request. The caller compiles this source into the application and
/// runs it with [`execute_with_rules_loop`]. Supports constant and vacancy
/// rates and sum observables, up to 1 MiB of emitted source. Parameters are
/// compiled in, so another parameter set needs another loop. Request
/// settings are not: one loop serves every Static Direct request for the
/// model bytes.
pub fn generate_rules_loop(prepared: &Prepared) -> Result<String, Diagnostic> {
    let net = network(prepared)?;
    let template = template(prepared, net)?;
    let hash = fingerprint(prepared, &template);
    Ok(template.replacen(PLACEHOLDER, &hash, 1))
}

fn bits(value: f64) -> String {
    format!("f64::from_bits({:#018x})", value.to_bits())
}

struct Rate {
    /// Vacancy base slot, if any.
    base: Option<usize>,
    expression: String,
}

/// The emitted source with a placeholder for its fingerprint.
fn template(prepared: &Prepared, net: &Network) -> Result<String, Diagnostic> {
    let unsupported = || {
        request_error(
            "/method",
            "generated rules loops support constant or vacancy rates and sum observables",
            "Use interpreted Static Direct for other expressions or reductions.",
        )
    };
    let mut bases: BTreeMap<(u64, usize, usize), usize> = BTreeMap::new();
    let mut base_code = Vec::new();
    let mut rates = Vec::with_capacity(net.channels.len());
    for (i, c) in net.channels.iter().enumerate() {
        let (base, global) = match &net.rates[c.rule] {
            RateKernel::Constant(k) => (None, bits(*k)),
            RateKernel::Vacancy {
                limit,
                count,
                denominator,
                factors,
            } => {
                let key = (limit.to_bits(), *count, *denominator);
                let next = bases.len();
                let slot = *bases.entry(key).or_insert(next);
                if slot == next {
                    base_code.push(format!(
                        "let base_{slot} = ({} - values[{count}]).max(0.0) / values[{denominator}];",
                        bits(*limit)
                    ));
                }
                let mut expression = format!("base_{slot}");
                for &factor in factors {
                    write!(expression, " * {}", bits(factor)).unwrap();
                }
                (Some(slot), expression)
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
        rates.push(Rate {
            base,
            expression: format!(
                "{{ let weight = {weight}; let a = if weight == 0.0 {{ 0.0 }} else {{ ({global}) * (weight * {}) }}; \
                 if !(a.is_finite() && a >= 0.0) {{ return Err(host.rate_error({i})); }} r[{i}] = a; }}",
                bits(c.factor)
            ),
        });
    }
    // Every changed projection is refreshed after each event, including the
    // ones the interpreter defers to reports. Deferred projections cannot
    // overflow: lowering allows at most 256 species, and `plan_observables`
    // defers only weights bounded for 1024 terms.
    if net
        .observables
        .iter()
        .any(|(_, reduce, _)| !matches!(reduce, Reduce::Sum))
    {
        return Err(unsupported());
    }
    let update = |channels: &mut dyn Iterator<Item = usize>| {
        let channels: Vec<_> = channels.collect();
        let mut code = String::new();
        let mut needed: Vec<usize> = channels.iter().filter_map(|&i| rates[i].base).collect();
        needed.sort_unstable();
        needed.dedup();
        for slot in needed {
            writeln!(code, "{}", base_code[slot]).unwrap();
        }
        for i in channels {
            writeln!(code, "{}", rates[i].expression).unwrap();
        }
        code
    };
    let n = net.channels.len();
    let mut dependencies = String::new();
    for (channel, row) in net.affected.iter().enumerate() {
        writeln!(
            dependencies,
            "{channel} => {{\n{}}}",
            update(&mut row.iter().copied())
        )
        .unwrap();
    }
    writeln!(dependencies, "_ => {{\n{}}}", update(&mut (0..n))).unwrap();
    let mut total = String::new();
    let mut unresolved = String::from("!total.is_finite()");
    let mut pick = String::new();
    for i in 0..n {
        writeln!(total, "total += r[{i}];").unwrap();
        write!(unresolved, " || (r[{i}] > 0.0 && r[{i}] < floor)").unwrap();
        writeln!(
            pick,
            "if r[{i}] > 0.0 {{ channel = {i}; if target < r[{i}] {{ break 'pick; }} target -= r[{i}]; }}"
        )
        .unwrap();
    }
    // The stop predicate runs after events that change a projection it
    // reads. The fused guarded fraction is emitted with the interpreter's
    // fused operations.
    let model = &prepared.compiled;
    let mut stop_reads = vec![false; model.observables.len()];
    if let Some(program) = &model.stop_program {
        program.reads(&mut stop_reads);
    }
    let emitted = model
        .stop
        .as_ref()
        .and_then(|e| e.rust_predicate(&model.params));
    let finish = "host.stopped(&counts[..], &mut values[..])?; return Ok(Some(next));";
    let stop = match (
        model.stop_program.as_ref().and_then(|p| p.fraction()),
        emitted,
    ) {
        (Some((total, part, fraction)), _) => format!(
            "if values[{total}] > 0.0 && values[{part}] >= {} * values[{total}] {{ {finish} }}",
            bits(fraction)
        ),
        // Other predicates are emitted with the bytecode's semantics; an
        // undefined value leaves the diagnostic to the host's bytecode.
        (None, Some(value)) => format!(
            "{{ let v: f64 = {value};\nif v.is_nan() {{ if host.stop(&counts[..], &mut values[..])? {{ return Ok(Some(next)); }} }} else if v != 0.0 {{ {finish} }} }}"
        ),
        (None, None) => {
            "if host.stop(&counts[..], &mut values[..])? { return Ok(Some(next)); }".into()
        }
    };
    let mut events = String::new();
    for (k, c) in net.channels.iter().enumerate() {
        let mut code = String::from("let mut live = host.live(); let mut created = 0u64;\n");
        for (j, &(s, d)) in c.changes.iter().enumerate() {
            writeln!(
                code,
                "let n{j} = match counts[{s}].checked_add_signed({d}) {{ Some(n) if n <= {MAX_COUNT} => n, _ => return Err(host.count_error()) }};\n\
                 live = live + u64::from(n{j} > 0) - u64::from(counts[{s}] > 0);\n\
                 created += u64::from(n{j} > 0 && counts[{s}] == 0);"
            )
            .unwrap();
        }
        writeln!(code, "host.admit(live)?;").unwrap();
        for (j, &(s, _)) in c.changes.iter().enumerate() {
            writeln!(code, "counts[{s}] = n{j};").unwrap();
        }
        writeln!(code, "host.fired({k}, created, live, next);").unwrap();
        for &index in &net.changed[k] {
            let (slot, _, terms) = &net.observables[index];
            write!(code, "{{ let v: f64 = [").unwrap();
            for &(species, weight) in terms {
                write!(code, "counts[{species}] as f64 * {},", bits(weight)).unwrap();
            }
            writeln!(
                code,
                "].into_iter().sum(); if !v.is_finite() {{ return Err(host.overflow_error()); }} values[{slot}] = v; }}"
            )
            .unwrap();
        }
        let stops = net.changed[k]
            .iter()
            .any(|&index| stop_reads[net.observables[index].0]);
        writeln!(
            events,
            "{k} => {{\n{code}if host.observes({k}) {{ host.observe({k}, &values[..], next); }}\n{}}}",
            if stops { stop.as_str() } else { "" }
        )
        .unwrap();
    }
    let species = net.species.len();
    let observables = prepared.compiled.observables.len();
    let hash = PLACEHOLDER;
    let resolution = bits(crate::engine::RESOLUTION);
    let source = format!(
        r#"// Generated by CRN-craft. Build with the same floating-point policy.
pub struct GeneratedRulesLoop;
impl crn_craft_rules::RulesLoop for GeneratedRulesLoop {{
    fn fingerprint(&self) -> &str {{ "{hash}" }}
    #[allow(clippy::all, unused_mut, unused_variables, unused_assignments, unreachable_code, unused_labels, unused_parens, unused_braces)]
    fn run(&self, host: &mut crn_craft_rules::LoopHost<'_, '_>, counts: &mut [u64],
           values: &mut [f64]) -> Result<Option<f64>, crn_craft_rules::LoopHalt> {{
        let counts: &mut [u64; {species}] = counts.try_into().expect("lowered species");
        let values: &mut [f64; {observables}] = values.try_into().expect("observables");
        let mut r = [0.0f64; {n}];
        let mut previous = usize::MAX;
        loop {{
            host.tick()?;
            match previous {{
{dependencies}            }}
            let mut total = -0.0f64;
{total}            if total == 0.0 {{ return host.exhausted(&counts[..], &mut values[..]); }}
            let floor = total * {resolution};
            if {unresolved} {{ return Err(host.resolution_error(total)); }}
            let next = host.next_time(total)?;
            let mut target = host.uniform()? * total;
            let mut channel = usize::MAX;
            'pick: {{
{pick}            }}
            if channel == usize::MAX {{ return host.exhausted(&counts[..], &mut values[..]); }}
            previous = channel;
            if !host.reach(next, &counts[..], &mut values[..])? {{ return Ok(None); }}
            host.begin_event()?;
            match channel {{
{events}                _ => unreachable!("validated channel"),
            }}
        }}
    }}
}}
"#
    );
    if source.len() > MAX_SOURCE {
        return Err(request_error(
            "/method",
            "the generated rules loop would exceed 1 MiB of source",
            "Use the generated rates kernel or interpreted Static Direct for densely coupled networks.",
        ));
    }
    Ok(source)
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
        // it would after a change in lowering. A stale loop must not pass.
        let model = include_bytes!("../../../models/rules/mtdna-three-class.rules.json");
        let request = serde_json::json!({
            "schema": "crn-craft.rules-run/0.2",
            "model": {"path": "m.rules.json", "sha256": crn_craft::sha256(model)},
            "time": {"end": 1.0, "report_times": [0.0, 1.0]},
            "ensemble": {"replicates": 1, "root_seed": 1},
            "outputs": {"statistics": ["mean"], "keep_paths": 0, "keep_instances": 0},
            "budget": {"threads": 1, "max_events_per_replicate": 10u64,
                       "max_sampling_steps_per_replicate": 100u64, "max_instances": 10u64,
                       "max_output_bytes": 1u64 << 20},
            "method": {"kind": "static_direct"}
        });
        let mut prepared = crate::prepare(model, &serde_json::to_vec(&request).unwrap()).unwrap();
        let before = generate_rules_loop(&prepared).unwrap();
        prepared.flat.as_mut().unwrap().observables[0].2[0].1 *= 2.0;
        let after = generate_rules_loop(&prepared).unwrap();
        assert_ne!(before, after);
        assert_ne!(fingerprint_of(&before), fingerprint_of(&after));
    }
}
