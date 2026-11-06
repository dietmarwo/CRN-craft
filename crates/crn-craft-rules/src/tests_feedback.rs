//! Regression tests from the first user feedback on the 0.2.0 interface.
use super::*;
use serde_json::{Value, json};

const THREE_CLASS: &str = include_str!("../../../models/rules/mtdna-three-class.rules.json");

/// The default `max_output_bytes` of the Python helper `rules.ensemble`.
const PYTHON_DEFAULT_OUTPUT: u64 = 64 << 20;

/// The request `crn_craft.rules.ensemble` builds for `model` with an event
/// log, as in the user manual's chapter on logs: 25,920 hours, two report
/// times, two threads and the Python default budgets.
fn python_request(model: &[u8], method: &str, replicates: u64, log: &Value) -> Value {
    json!({
        "schema": if method == "mtdna_direct" {
            "crn-craft.rules-run/0.3"
        } else {
            "crn-craft.rules-run/0.2"
        },
        "model": {"path": "model.rules.json", "sha256": crn_craft::sha256(model)},
        "time": {"end": 25_920.0, "report_times": [0.0, 25_920.0]},
        "ensemble": {"replicates": replicates, "root_seed": 4},
        "outputs": {"statistics": ["mean", "variance"], "keep_paths": 0, "keep_instances": 0},
        "budget": {"threads": 2, "max_events_per_replicate": 1_000_000u64,
                   "max_sampling_steps_per_replicate": 8_000_000u64,
                   "max_instances": 100_000u64, "max_output_bytes": PYTHON_DEFAULT_OUTPUT},
        "method": {"kind": method},
        "event_log": log
    })
}

fn prepare_json(model: &[u8], request: &Value) -> Result<Prepared, Vec<Diagnostic>> {
    prepare(model, &serde_json::to_vec(request).unwrap())
}

/// The reported request: the three-class model with a bounded log of the two
/// mutation rules, at the replicate counts the user manual documents (the
/// log chapter's 100, the first ensemble's 200, the sweep's 400 and the
/// tracker chapter's 500), with the Python default output budget.
#[test]
fn reported_event_log_fits_the_python_default_budget() {
    let model = THREE_CLASS.as_bytes();
    let log = json!({"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]});
    for method in ["static_direct", "mtdna_direct"] {
        for replicates in [100, 200, 400, 500] {
            let request = python_request(model, method, replicates, &log);
            if let Err(diagnostics) = prepare_json(model, &request) {
                panic!("{method} with {replicates} replicates: {diagnostics:?}");
            }
        }
    }
}

/// The reservation's bound and the measured compact JSON of the retained log
/// of every replicate, for logs filled to `max_entries`.
fn log_bound_and_measured(model: &[u8], method: Value, entries: usize) -> (u128, Vec<usize>) {
    let mut request = python_request(model, "tree", 3, &json!({"max_entries": entries}));
    request["method"] = method;
    request["time"] = json!({"end": 2_000.0, "report_times": [0.0, 2_000.0]});
    let prepared = prepare_json(model, &request).unwrap();
    let sizes = observation_sizes(
        &prepared.compiled,
        prepared.request.method,
        prepared.flat.as_ref(),
    );
    let reserved = prepared.observations.reservation(&sizes);
    let result = execute(&prepared, &Control::default()).unwrap();
    assert_eq!(result.termination, Termination::Completed);
    let measured = result
        .observations
        .iter()
        .map(|o| {
            assert_eq!(o.log.len(), entries, "log not filled");
            serde_json::to_vec(&o.log).unwrap().len()
        })
        .collect();
    (reserved.entries * reserved.entry, measured)
}

/// The tighter bound still covers what a full log holds: exact steps of the
/// three-class model, leaps, and the clone model on `tree`, whose new
/// identities and drawn weights are attributes, plus a type whose drawn
/// attributes print at full precision.
#[test]
fn reservation_covers_a_full_log() {
    let three = THREE_CLASS.as_bytes();
    let clones = include_bytes!("../../../models/rules/mtdna-fixed.rules.json");
    let noisy = serde_json::to_vec(&json!({
        "schema": "crn-craft.rules-model/0.1", "name": "noisy",
        "parameters": {"k": 1.0},
        "types": [{"name": "A", "attributes": ["u", "v", "w"]}],
        "initial": [{"type": "A", "attributes": {"u": 0.1, "v": -0.3, "w": 1e-300}, "count": 1}],
        "observables": {"A": {"type": "A"}, "spread": {"type": "A", "weight": "x.u * 3.3"}},
        "rules": [{"name": "redraw", "reactants": [{"label": "a", "type": "A"}], "rate": "k",
                   "products": [{"type": "A", "attributes": {
                       "u": "normal(0, 1) * 1e-7", "v": "-uniform(0, 1)", "w": "normal(0, 1) * 1e300"}}]}]
    }))
    .unwrap();
    let cases: [(&[u8], Value); 6] = [
        (three, json!({"kind": "static_direct"})),
        (three, json!({"kind": "next_reaction"})),
        (three, json!({"kind": "tree"})),
        (three, json!({"kind": "tau", "epsilon": 0.03})),
        (clones, json!({"kind": "tree"})),
        (&noisy, json!({"kind": "tree"})),
    ];
    for (model, method) in cases {
        let (bound, measured) = log_bound_and_measured(model, method.clone(), 50);
        for bytes in measured {
            assert!(
                bytes as u128 <= bound,
                "{method}: {bytes} log bytes above the bound {bound}"
            );
        }
    }
}

/// The diagnostic itemizes the reservation and names the budget that fits.
#[test]
fn reservation_error_names_bytes_and_fix() {
    let model = THREE_CLASS.as_bytes();
    let log = json!({"max_entries": 50, "rules": ["mutation_mt"]});
    let mut request = python_request(model, "static_direct", 20_000, &log);
    let diagnostics = prepare_json(model, &request).err().unwrap();
    let error = &diagnostics[0];
    assert_eq!(
        (error.code.as_str(), error.path.as_str()),
        ("CRN-R020", "/budget")
    );
    // Three participants of three attributes and six observables:
    // 160 + 6 * 25 + 3 * (62 + 3 * 25) = 721 bytes per entry.
    let total = 20_000u128 * 50 * 721;
    assert!(
        error
            .message
            .contains("20000 replicates x 36050 bytes (50 log entries x 721 bytes)")
            && error.message.contains(&format!("= {total} bytes"))
            && error
                .message
                .contains(&format!("max_output_bytes = {}", 64u64 << 20)),
        "{}",
        error.message
    );
    assert!(
        error.hint.contains(&format!("at least {total}")),
        "{}",
        error.hint
    );
    // That budget is enough.
    request["budget"]["max_output_bytes"] = json!(u64::try_from(total).unwrap());
    assert!(prepare_json(model, &request).is_ok());
    request["budget"]["max_output_bytes"] = json!(u64::try_from(total).unwrap() - 1);
    assert!(prepare_json(model, &request).is_err());
}
