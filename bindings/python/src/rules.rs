use crate::{failure, rules_payload, run_interruptible};
use crn_craft_rules::observe::{Event, EventObserver, ObserverAction, log_entry};
use crn_craft_rules::{Control, Prepared};
use pyo3::{
    prelude::*,
    types::{PyBytes, PyDict},
};
use std::{
    path::PathBuf,
    sync::{Arc, Mutex},
};

#[pyclass(frozen)]
struct PreparedRules(Arc<Prepared>);

struct PythonObserver<'a> {
    callback: Py<PyAny>,
    error: Mutex<Option<PyErr>>,
    control: &'a Control,
}

impl EventObserver for PythonObserver<'_> {
    fn on_event(&self, event: &Event<'_>) -> ObserverAction {
        // The flag, not the cancel file: the engine polls the file itself.
        if self.control.cancel_flagged() {
            return ObserverAction::Continue;
        }
        let result = Python::attach(|py| -> PyResult<bool> {
            let text = serde_json::to_string(&log_entry(event))
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
            let answer = self.callback.call1(py, (event.replicate, text))?;
            if answer.is_none(py) {
                Ok(true)
            } else {
                answer.extract::<bool>(py)
            }
        });
        match result {
            Ok(true) => ObserverAction::Continue,
            Ok(false) => ObserverAction::Stop,
            Err(error) => {
                let mut first = self.error.lock().expect("observer lock");
                if first.is_none() {
                    *first = Some(error);
                }
                self.control.cancel();
                ObserverAction::Stop
            }
        }
    }
}

#[pymethods]
impl PreparedRules {
    fn bind_parameters(
        &self,
        py: Python<'_>,
        overrides: std::collections::BTreeMap<String, f64>,
    ) -> PyResult<Self> {
        let prepared = Arc::clone(&self.0);
        py.detach(move || prepared.bind_parameters(&overrides))
            .map(|p| Self(Arc::new(p)))
            .map_err(|d| failure(&d))
    }

    #[pyo3(signature = (cancel_file=None, callback=None))]
    fn simulate<'py>(
        &self,
        py: Python<'py>,
        cancel_file: Option<String>,
        callback: Option<Py<PyAny>>,
    ) -> PyResult<(String, Bound<'py, PyDict>)> {
        let prepared = Arc::clone(&self.0);
        let mut control = Control::default();
        control.cancel_file = cancel_file.map(PathBuf::from);
        let result = run_interruptible(py, control, Control::cancel, move |c| {
            if let Some(callback) = callback {
                let observer = PythonObserver {
                    callback,
                    error: Mutex::new(None),
                    control: c,
                };
                let result = crn_craft_rules::execute_with_observer(&prepared, c, &observer);
                if let Some(error) = observer.error.into_inner().expect("observer lock") {
                    return Err(error);
                }
                result.map_err(|d| failure(&d))
            } else {
                crn_craft_rules::execute(&prepared, c).map_err(|d| failure(&d))
            }
        })??;
        let (metadata, columns) = rules_payload::encode(result);
        let out = PyDict::new(py);
        for (name, data) in columns {
            out.set_item(name, PyBytes::new(py, &data))?;
        }
        Ok((metadata, out))
    }

    /// The generated Static Direct event loop of this prepared run, or the
    /// generator's `CRN-R020` refusal.
    fn loop_source(&self, py: Python<'_>) -> PyResult<String> {
        let prepared = Arc::clone(&self.0);
        py.detach(move || crn_craft_rules::generate_rules_loop(&prepared))
            .map_err(|d| failure(&d))
    }
}

/// The driver sources of a compiled-loop project and the minimum Rust
/// version. `compile_loop` adds the generated loop, the model and the
/// manifest.
#[pyfunction]
fn rules_loop_driver() -> (&'static str, Vec<(&'static str, &'static str)>) {
    (
        env!("CARGO_PKG_RUST_VERSION"),
        vec![
            ("src/main.rs", include_str!("rules_loop_main.rs")),
            ("src/rules_payload.rs", include_str!("rules_payload.rs")),
        ],
    )
}

#[pyfunction]
fn describe_rules_json() -> String {
    crn_craft_rules::describe().to_string()
}

#[pyfunction]
fn check_rules_json(model: &[u8]) -> String {
    serde_json::to_string(&crn_craft_rules::check_model(model)).expect("serializable check")
}

/// The rules model converted from a network model's bytes; a refusal raises
/// `NativeError` with the diagnostics.
#[pyfunction]
fn convert_network<'py>(py: Python<'py>, model: &[u8]) -> PyResult<Bound<'py, PyBytes>> {
    crn_craft::rules::convert::from_network(model)
        .map(|bytes| PyBytes::new(py, &bytes))
        .map_err(|d| failure(&d))
}

#[pyfunction]
fn prepare_rules(py: Python<'_>, model: Vec<u8>, request: Vec<u8>) -> PyResult<PreparedRules> {
    py.detach(move || crn_craft_rules::prepare(&model, &request))
        .map(|p| PreparedRules(Arc::new(p)))
        .map_err(|d| failure(&d))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PreparedRules>()?;
    module.add_function(wrap_pyfunction!(describe_rules_json, module)?)?;
    module.add_function(wrap_pyfunction!(check_rules_json, module)?)?;
    module.add_function(wrap_pyfunction!(prepare_rules, module)?)?;
    module.add_function(wrap_pyfunction!(convert_network, module)?)?;
    module.add_function(wrap_pyfunction!(rules_loop_driver, module)?)?;
    Ok(())
}
