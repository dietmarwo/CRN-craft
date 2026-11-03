use crate::{failure, run_interruptible};
use crn_craft_rules::{Control, Prepared, Record, Termination};
use pyo3::{prelude::*, types::{PyBytes, PyDict}};
use std::{path::PathBuf, sync::Arc};

#[pyclass(frozen)]
struct PreparedRules(Arc<Prepared>);

fn columns<'py>(py: Python<'py>, records: &[Record]) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);
    macro_rules! column {
        ($name:literal, $value:expr) => {{
            let bytes: Vec<u8> = records.iter().flat_map($value).collect();
            out.set_item($name, PyBytes::new(py, &bytes))?;
        }};
    }
    column!("id", |r| r.id.to_le_bytes());
    column!("started", |r| [u8::from(r.started)]);
    column!("events", |r| r.events.to_le_bytes());
    column!("stop_time", |r| r.stop_time.to_le_bytes());
    column!("saved_points", |r| (r.saved_points as u64).to_le_bytes());
    column!("rng_words", |r| r.rng_words.to_le_bytes());
    column!("peak_instances", |r| r.peak_instances.to_le_bytes());
    column!("created_instances", |r| r.created_instances.to_le_bytes());
    column!("termination", |r| [match r.termination {
        Termination::Completed => 0, Termination::Stopped => 1,
        Termination::EventLimit => 2, Termination::SamplingLimit => 3,
        Termination::InstanceLimit => 4, Termination::WallTimeLimit => 5,
        Termination::Cancelled => 6, Termination::Failed => 7,
        Termination::OutputLimit => 8,
    }]);
    let counts: Vec<u8> = records.iter().flat_map(|r| &r.rule_events)
        .flat_map(|n| n.to_le_bytes()).collect();
    out.set_item("rule_events", PyBytes::new(py, &counts))?;
    Ok(out)
}

#[pymethods]
impl PreparedRules {
    #[pyo3(signature = (cancel_file=None))]
    fn simulate<'py>(
        &self, py: Python<'py>, cancel_file: Option<String>,
    ) -> PyResult<(String, Bound<'py, PyDict>)> {
        let prepared = Arc::clone(&self.0);
        let mut control = Control::default();
        control.cancel_file = cancel_file.map(PathBuf::from);
        let mut result = run_interruptible(py, control, Control::cancel, move |c| {
            crn_craft_rules::execute(&prepared, c)
        })?.map_err(|d| failure(&d))?;
        let records = std::mem::take(&mut result.replicates);
        let metadata = serde_json::json!({
            "result": result,
            "record_diagnostics": records.iter().filter_map(|r| {
                r.diagnostic.as_ref().map(|d| (r.id, d))
            }).collect::<Vec<_>>(),
        });
        Ok((metadata.to_string(), columns(py, &records)?))
    }
}

#[pyfunction]
fn check_rules_json(model: &[u8]) -> String {
    serde_json::to_string(&crn_craft_rules::check_model(model)).expect("serializable check")
}

#[pyfunction]
fn prepare_rules(py: Python<'_>, model: Vec<u8>, request: Vec<u8>) -> PyResult<PreparedRules> {
    py.detach(move || crn_craft_rules::prepare(&model, &request))
        .map(|p| PreparedRules(Arc::new(p))).map_err(|d| failure(&d))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PreparedRules>()?;
    module.add_function(wrap_pyfunction!(check_rules_json, module)?)?;
    module.add_function(wrap_pyfunction!(prepare_rules, module)?)?;
    Ok(())
}
