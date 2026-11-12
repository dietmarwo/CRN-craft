//! A rules result as `crn_craft.rules.Result` reads it: JSON metadata plus
//! little-endian record columns.
//!
//! The extension compiles this file, and `crn_craft.rules.compile_loop`
//! copies it unchanged into every compiled event loop project, so both
//! paths encode results with the same code. It may use only `crn_craft`
//! and `serde_json`.

use crn_craft::rules::{Record, RulesResult, Termination};

/// Splits a result into its metadata document (the result without its
/// records, plus the record diagnostics by replicate id) and the record
/// columns, in a fixed order.
pub fn encode(mut result: RulesResult) -> (String, Vec<(&'static str, Vec<u8>)>) {
    let records = std::mem::take(&mut result.replicates);
    let metadata = serde_json::json!({
        "result": result,
        "record_diagnostics": records.iter().filter_map(|r| {
            r.diagnostic.as_ref().map(|d| (r.id, d))
        }).collect::<Vec<_>>(),
    });
    (metadata.to_string(), columns(&records))
}

fn columns(records: &[Record]) -> Vec<(&'static str, Vec<u8>)> {
    fn column<const N: usize>(records: &[Record], value: impl Fn(&Record) -> [u8; N]) -> Vec<u8> {
        records.iter().flat_map(value).collect()
    }
    vec![
        ("id", column(records, |r| r.id.to_le_bytes())),
        ("started", column(records, |r| [u8::from(r.started)])),
        ("events", column(records, |r| r.events.to_le_bytes())),
        ("stop_time", column(records, |r| r.stop_time.to_le_bytes())),
        (
            "saved_points",
            column(records, |r| (r.saved_points as u64).to_le_bytes()),
        ),
        ("rng_words", column(records, |r| r.rng_words.to_le_bytes())),
        (
            "peak_instances",
            column(records, |r| r.peak_instances.to_le_bytes()),
        ),
        (
            "created_instances",
            column(records, |r| r.created_instances.to_le_bytes()),
        ),
        (
            "termination",
            column(records, |r| {
                [match r.termination {
                    Termination::Completed => 0,
                    Termination::Stopped => 1,
                    Termination::EventLimit => 2,
                    Termination::SamplingLimit => 3,
                    Termination::InstanceLimit => 4,
                    Termination::WallTimeLimit => 5,
                    Termination::Cancelled => 6,
                    Termination::Failed => 7,
                    Termination::OutputLimit => 8,
                }]
            }),
        ),
        (
            "rule_events",
            records
                .iter()
                .flat_map(|r| &r.rule_events)
                .flat_map(|n| n.to_le_bytes())
                .collect(),
        ),
    ]
}
