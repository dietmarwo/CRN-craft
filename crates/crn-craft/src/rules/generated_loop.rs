//! Ahead-of-time event loops for a lowered Static Direct model.
//!
//! [`generate_rules_loop`] emits the whole Static Direct event of one
//! prepared model as Rust source; [`execute_with_rules_loop`] runs a
//! compiled loop with a [`LoopHost`] per replicate. Benchmark it before
//! adopting it: on networks with many constant-rate channels it can be
//! slower than the interpreter. See the
//! [rules kernel checklist](https://github.com/dietmarwo/CRN-craft/blob/main/docs/specialized-kernels.md#generated-rules-event-loops).
use crate::rules::expr::RateKernel;
use crate::rules::flat::Network;
use crate::rules::model::{MAX_COUNT, Reduce};
use crate::rules::{
    Control, Diagnostic, LoopHalt, LoopHost, Method, Prepared, RulesResult, request_error,
};
use std::collections::BTreeMap;
use std::collections::btree_map::Entry;
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
    crate::sha256(format!("rules-loop/1\n{}\n{template}", prepared.model_sha256).as_bytes())
}

/// Run a prepared Static Direct request with a generated event loop. Results
/// are byte-identical to [`crate::rules::execute`] for the same request. Fails with
/// `CRN-R020` at `/model` unless the loop was generated for these exact model
/// bytes by this version of the generator. Observers are not supported on
/// this path.
///
/// Trackers, event logs, both waiting-time policies, any report grid, budget
/// and worker count work with the same compiled loop.
///
/// # Errors
///
/// `CRN-R020` at `/method` if the request does not select `static_direct`
/// or the source would exceed 1 MiB, `CRN-R020` at `/model` on a
/// fingerprint mismatch, and otherwise the errors of [`crate::rules::execute`].
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
    crate::rules::execute_inner(
        prepared,
        control,
        &crate::rules::NoObserver,
        None,
        Some(kernel),
    )
}

/// Emit a Rust module defining `GeneratedRulesLoop` for a prepared Static
/// Direct request. The caller compiles this source into the application and
/// runs it with [`execute_with_rules_loop`]. Supports every model that
/// Static Direct runs, up to 1 MiB of emitted source: constant, vacancy and
/// expression rates, whole-rule `propensity` expressions, higher-order rules
/// (through [`crate::rules::combination_rate`], the interpreter's
/// arithmetic), and `sum`,
/// `distinct`, `max` and `min` observables. Rate expressions are emitted
/// with the interpreter's operations in its order, so results stay
/// byte-identical: NaN propagation through logic and `if`, short-circuit
/// `&&` and `||`, `min` and `max` that keep NaN, IEEE division and the same
/// `libm` functions. Parameters are compiled in, so another parameter set
/// needs another loop. Request settings are not: one loop serves every
/// Static Direct request for the model bytes.
///
/// The model's stop predicate is compiled in with the interpreter's
/// semantics; an undefined value is passed back to the host for the
/// interpreter's diagnostic. The example `generate_rules_loop` in the
/// repository writes the source for a model and request file.
///
/// # Errors
///
/// `CRN-R020` at `/method` if the request does not select `static_direct`
/// (models with ATP, continuous variables or instances that change identity
/// need another method and cannot select it), or if the source would
/// exceed 1 MiB, which densely coupled networks reach because rate code is
/// repeated in every dependency row.
///
/// ```
/// use crn_craft::rules::{generate_rules_loop, prepare};
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
/// let source = generate_rules_loop(&prepared).unwrap();
/// assert!(source.contains("impl crn_craft_rules::RulesLoop for GeneratedRulesLoop"));
///
/// // Request settings are not compiled in: Ziggurat waits share the loop.
/// let mut ziggurat: serde_json::Value = serde_json::from_slice(
///     &request(json!({"kind": "static_direct", "waiting_time": "ziggurat"}))).unwrap();
/// ziggurat["schema"] = json!("crn-craft.rules-run/0.4");
/// let other = prepare(&model, &serde_json::to_vec(&ziggurat).unwrap()).unwrap();
/// assert_eq!(generate_rules_loop(&other).unwrap(), source);
/// ```
pub fn generate_rules_loop(prepared: &Prepared) -> Result<String, Diagnostic> {
    let net = network(prepared)?;
    let template = template(prepared, net)?;
    let hash = fingerprint(prepared, &template);
    Ok(template.replacen(PLACEHOLDER, &hash, 1))
}

fn bits(value: f64) -> String {
    format!("f64::from_bits({:#018x})", value.to_bits())
}

/// The global factor of one channel's rate.
enum Global {
    /// Emitted inline: a constant, or a vacancy base times its factors.
    Inline(String),
    /// The rate or propensity expression of this rule, shared by every
    /// channel of the rule in a dependency row.
    Rule(usize),
}

struct Rate {
    /// Vacancy base slot, if any.
    base: Option<usize>,
    global: Global,
    /// The reactant weight, as the interpreter computes it.
    weight: String,
    /// The lowered channel factor.
    factor: f64,
    /// A higher-order channel with a `rate`: its (species, count) list as
    /// Rust source, for `combination_rate`.
    higher_order: Option<String>,
}

/// The emitted source with a placeholder for its fingerprint.
#[allow(clippy::too_many_lines)]
fn template(prepared: &Prepared, net: &Network) -> Result<String, Diagnostic> {
    let model = &prepared.compiled;
    let mut bases: BTreeMap<(u64, usize, usize), usize> = BTreeMap::new();
    let mut base_code = Vec::new();
    // Rate and propensity expressions with the interpreter's evaluation
    // order and NaN, infinity and division semantics, per rule.
    let mut rule_code: BTreeMap<usize, String> = BTreeMap::new();
    let mut rates = Vec::with_capacity(net.channels.len());
    for c in &net.channels {
        let (base, global) = match &net.rates[c.rule] {
            RateKernel::Constant(k) => (None, Global::Inline(bits(*k))),
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
                (Some(slot), Global::Inline(expression))
            }
            RateKernel::Expression(expr) => {
                if let Entry::Vacant(entry) = rule_code.entry(c.rule) {
                    // Validation lets rates and propensities read only
                    // parameters and observables, and never draw.
                    entry.insert(expr.rust_value(&model.params).ok_or_else(|| {
                        request_error(
                            "/method",
                            "generated rules loops support deterministic rate expressions over parameters and observables",
                            "Use interpreted Static Direct for this rule.",
                        )
                    })?);
                }
                (None, Global::Rule(c.rule))
            }
        };
        // The interpreter's reactant weight: the combinatorial count, or for
        // a whole-rule propensity 1 while every pattern can take a distinct
        // molecule and 0 otherwise. A higher-order channel's weight is its
        // availability; the interpreter's `combination_rate` then gives the
        // rate, emitted as a call.
        if let Some(list) = &c.multi {
            let available: Vec<String> = list
                .iter()
                .map(|(s, k)| format!("counts[{s}] >= {k}"))
                .collect();
            let reactants: Vec<String> = list.iter().map(|(s, k)| format!("({s}, {k})")).collect();
            rates.push(Rate {
                base,
                global,
                weight: format!("f64::from(u8::from({}))", available.join(" && ")),
                factor: c.factor,
                higher_order: (!c.whole).then(|| reactants.join(", ")),
            });
            continue;
        }
        let weight = match (c.reactants.as_slice(), c.whole) {
            ([], _) => "1.0".into(),
            ([s], false) => format!("counts[{s}] as f64"),
            ([s, t], false) if s == t => {
                format!("counts[{s}] as f64 * counts[{s}].saturating_sub(1) as f64 * 0.5")
            }
            ([s, t], false) => format!("counts[{s}] as f64 * counts[{t}] as f64"),
            ([s], true) => format!("f64::from(u8::from(counts[{s}] > 0))"),
            ([s, t], true) if s == t => format!("f64::from(u8::from(counts[{s}] > 1))"),
            ([s, t], true) => {
                format!("f64::from(u8::from(counts[{s}] > 0 && counts[{t}] > 0))")
            }
            _ => unreachable!("validated arity"),
        };
        rates.push(Rate {
            base,
            global,
            weight,
            factor: c.factor,
            higher_order: None,
        });
    }
    // Every changed projection is refreshed after each event, including the
    // ones the interpreter defers to reports. Deferred projections cannot
    // overflow: `plan_observables` defers only weights bounded for the
    // projection's number of terms (at least 1024).
    let update = |channels: &mut dyn Iterator<Item = usize>| {
        let channels: Vec<_> = channels.collect();
        let mut code = String::new();
        let mut needed: Vec<usize> = channels.iter().filter_map(|&i| rates[i].base).collect();
        needed.sort_unstable();
        needed.dedup();
        for slot in needed {
            writeln!(code, "{}", base_code[slot]).unwrap();
        }
        // An expression shared by several channels of the row is evaluated
        // once, as the interpreter's per-event rate cache does. Evaluation
        // has no side effects, so evaluating it when every weight is 0
        // changes no value.
        let mut uses: BTreeMap<usize, usize> = BTreeMap::new();
        for &i in &channels {
            if let Global::Rule(rule) = rates[i].global {
                *uses.entry(rule).or_insert(0) += 1;
            }
        }
        for (&rule, &n) in &uses {
            if n > 1 {
                writeln!(code, "let g_{rule}: f64 = {};", rule_code[&rule]).unwrap();
            }
        }
        for i in channels {
            let rate = &rates[i];
            let global = match rate.global {
                Global::Inline(ref code) => code.clone(),
                Global::Rule(rule) if uses[&rule] > 1 => format!("g_{rule}"),
                Global::Rule(rule) => rule_code[&rule].clone(),
            };
            let value = match &rate.higher_order {
                None => format!("({global}) * (weight * {})", bits(rate.factor)),
                Some(list) => format!(
                    "crn_craft_rules::combination_rate({global}, {}, &counts[..], &[{list}])",
                    bits(rate.factor)
                ),
            };
            writeln!(
                code,
                "{{ let weight = {}; let a = if weight == 0.0 {{ 0.0 }} else {{ {value} }}; \
                 if !(a.is_finite() && a >= 0.0) {{ return Err(host.rate_error({i})); }} r[{i}] = a; }}",
                rate.weight,
            )
            .unwrap();
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
    let mut stop_reads = vec![false; model.observables.len()];
    if let Some(program) = &model.stop_program {
        program.reads(&mut stop_reads);
    }
    let emitted = model
        .stop
        .as_ref()
        .and_then(|e| e.rust_value(&model.params));
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
            // The interpreter's reductions, with the same iterator adapters
            // over the terms in the same order.
            let (slot, reduce, terms) = &net.observables[index];
            write!(code, "{{ let v: f64 = [").unwrap();
            for &(species, weight) in terms {
                if *reduce == Reduce::Sum {
                    write!(code, "counts[{species}] as f64 * {},", bits(weight)).unwrap();
                } else {
                    write!(code, "(counts[{species}] > 0).then_some({}),", bits(weight)).unwrap();
                }
            }
            let fold = match reduce {
                Reduce::Sum => "].into_iter().sum()",
                Reduce::Distinct => "].into_iter().flatten().sum()",
                Reduce::Max => "].into_iter().flatten().reduce(f64::max).unwrap_or(0.0)",
                Reduce::Min => "].into_iter().flatten().reduce(f64::min).unwrap_or(0.0)",
            };
            writeln!(
                code,
                "{fold}; if !v.is_finite() {{ return Err(host.overflow_error()); }} values[{slot}] = v; }}"
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
    let resolution = bits(crate::rules::engine::RESOLUTION);
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
        let before = generate_rules_loop(&prepared).unwrap();
        prepared.flat.as_mut().unwrap().observables[0].2[0].1 *= 2.0;
        let after = generate_rules_loop(&prepared).unwrap();
        assert_ne!(before, after);
        assert_ne!(fingerprint_of(&before), fingerprint_of(&after));
    }
}
