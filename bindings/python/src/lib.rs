//! Minimal JSON bridge; numerical semantics and validators live in the core.
mod rules;
mod rules_payload;

use pyo3::{
    create_exception,
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use std::{
    path::PathBuf,
    sync::{Arc, mpsc},
    time::Duration,
};

create_exception!(
    _native,
    NativeError,
    PyValueError,
    "Rejected request; the message is the JSON diagnostic list or object."
);

fn failure<T: serde::Serialize>(value: &T) -> PyErr {
    NativeError::new_err(
        serde_json::to_string(value).unwrap_or_else(|_| "serialization failure".into()),
    )
}

/// How often the waiting thread checks for Python signals.
const SIGNAL_INTERVAL: Duration = Duration::from_millis(50);

/// Runs `work` on a helper thread. The calling thread waits without the GIL
/// and checks for Python signals every 50 ms; a signal handler that raises
/// (Ctrl-C raises KeyboardInterrupt) cancels the run through the shared
/// control (`cancel`), and its exception is raised once the helper has
/// stopped. Signals are only delivered to the main thread, so calls from
/// other Python threads are not interruptible. A cancel file is polled as in
/// the CLI.
fn run_interruptible<C: Send + Sync + 'static, T: Send + 'static>(
    py: Python<'_>,
    control: C,
    cancel: fn(&C),
    work: impl FnOnce(&C) -> T + Send + 'static,
) -> PyResult<T> {
    let control = Arc::new(control);
    let shared = Arc::clone(&control);
    let (sender, mut receiver) = mpsc::channel();
    let helper = std::thread::Builder::new()
        .name("crn-craft".into())
        .spawn(move || {
            let _ = sender.send(work(&shared));
        })
        .map_err(|e| PyRuntimeError::new_err(format!("cannot start worker: {e}")))?;
    let mut signal = None;
    let value = loop {
        let (back, received) = py.detach(move || {
            let received = receiver.recv_timeout(SIGNAL_INTERVAL);
            (receiver, received)
        });
        receiver = back;
        match received {
            Ok(value) => break Some(value),
            Err(mpsc::RecvTimeoutError::Disconnected) => break None,
            Err(mpsc::RecvTimeoutError::Timeout) => {
                if signal.is_none() {
                    if let Err(error) = py.check_signals() {
                        cancel(&control);
                        signal = Some(error);
                    }
                }
            }
        }
    };
    let joined = py.detach(move || helper.join());
    if let Some(error) = signal {
        return Err(error);
    }
    match (value, joined) {
        (Some(value), Ok(())) => Ok(value),
        _ => Err(PyRuntimeError::new_err("CRN-craft worker thread panicked")),
    }
}

/// The unified discovery document, `crn-craft.describe/1.0`.
#[pyfunction]
fn describe_json() -> String {
    crn_craft::describe().to_string()
}
/// The report of `crn-craft check` on bytes (`crn_craft.check`).
#[pyfunction]
#[pyo3(signature = (model, request=None))]
fn check_any_json(py: Python<'_>, model: Vec<u8>, request: Option<Vec<u8>>) -> String {
    py.detach(move || {
        serde_json::to_string(&crn_craft::check_any(&model, request.as_deref()))
            .expect("finite check report")
    })
}
/// The report of `crn-craft check MODEL [REQUEST]` on files.
#[pyfunction]
#[pyo3(signature = (model, request=None))]
fn check_any_files_json(py: Python<'_>, model: String, request: Option<String>) -> String {
    py.detach(move || {
        serde_json::to_string(&crn_craft::check_any_files(
            std::path::Path::new(&model),
            request.as_deref().map(std::path::Path::new),
        ))
        .expect("finite check report")
    })
}
/// `crn-craft run MODEL REQUEST [--archive NEW_DIRECTORY]` for
/// `python -m crn_craft run`: the text the native CLI prints (without the
/// final newline) and whether it exits 0. A result is compact JSON with
/// sorted keys; an archive report and a refusal are pretty-printed.
#[pyfunction]
#[pyo3(signature = (model, request, archive=None, cancel_file=None))]
fn run_files_json(
    py: Python<'_>,
    model: String,
    request: String,
    archive: Option<String>,
    cancel_file: Option<String>,
) -> PyResult<(String, bool)> {
    let mut control = crn_craft::rules::Control::default();
    control.cancel_file = cancel_file.map(PathBuf::from);
    run_interruptible(py, control, crn_craft::rules::Control::cancel, move |c| {
        let (model, request) = (std::path::Path::new(&model), std::path::Path::new(&request));
        let pretty = |value: &serde_json::Value| {
            serde_json::to_string_pretty(value).expect("finite JSON")
        };
        let compact = |value: serde_json::Value| value.to_string();
        let refusal = |errors: Vec<crn_craft::diagnostic::Diagnostic>| {
            (
                pretty(&serde_json::json!({"success": false, "diagnostics": errors})),
                false,
            )
        };
        if let Some(directory) = archive {
            let report = crn_craft::archive::v1::run_files(
                model,
                request,
                std::path::Path::new(&directory),
                c,
            );
            let success = report.valid && report.simulation_complete;
            (
                pretty(&serde_json::to_value(report).expect("finite archive report")),
                success,
            )
        } else if crn_craft::fsp::v1::is_fsp_request_file(request) {
            match crn_craft::fsp::v1::solve_files(model, request, c) {
                Ok(result) => (compact(serde_json::to_value(result).expect("finite")), true),
                Err(errors) => refusal(errors),
            }
        } else {
            match crn_craft::rules::solve_files(model, request, c) {
                Ok(result) => {
                    let complete = result.termination == crn_craft::rules::Termination::Completed;
                    (compact(serde_json::to_value(result).expect("finite")), complete)
                }
                Err(errors) => refusal(errors),
            }
        }
    })
}
#[pyfunction]
fn analyze_json(model: &[u8]) -> PyResult<String> {
    serde_json::to_string(&crn_craft::analyze_bytes(model, None))
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    rules::register(module)?;
    module.add_function(wrap_pyfunction!(describe_json, module)?)?;
    module.add_function(wrap_pyfunction!(check_any_json, module)?)?;
    module.add_function(wrap_pyfunction!(check_any_files_json, module)?)?;
    module.add_function(wrap_pyfunction!(run_files_json, module)?)?;
    module.add_function(wrap_pyfunction!(analyze_json, module)?)?;
    module.add_function(wrap_pyfunction!(verify_json, module)?)?;
    module.add_function(wrap_pyfunction!(schema_json, module)?)?;
    module.add_function(wrap_pyfunction!(export_json, module)?)?;
    module.add_function(wrap_pyfunction!(analyze_state_json, module)?)?;
    module.add("NativeError", module.py().get_type::<NativeError>())?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

/// The report of `crn-craft verify DIRECTORY`: `crn-craft.archive-check/1.0`.
#[pyfunction]
fn verify_json(py: Python<'_>, directory: String) -> String {
    py.detach(move || {
        serde_json::to_string(&crn_craft::archive::v1::verify(std::path::Path::new(
            &directory,
        )))
        .unwrap()
    })
}
#[pyfunction]
fn schema_json(kind: String) -> PyResult<String> {
    crn_craft::schema(&kind)
        .or_else(|| crn_craft::rules::schema(&kind))
        .map(|v| v.to_string())
        .ok_or_else(|| PyValueError::new_err("unknown schema"))
}
#[pyfunction]
fn export_json(model: Vec<u8>) -> String {
    serde_json::to_string(&crn_craft::export::export_gillespy2(&model)).unwrap()
}
#[pyfunction]
#[pyo3(signature = (model, state=None))]
fn analyze_state_json(model: Vec<u8>, state: Option<Vec<u8>>) -> String {
    serde_json::to_string(&crn_craft::analyze_bytes(&model, state.as_deref())).unwrap()
}
