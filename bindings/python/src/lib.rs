//! Minimal JSON bridge; numerical semantics and validators live in the core.
use crn_craft::simulation::Control;
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
/// control, and its exception is raised once the helper has stopped. Signals
/// are only delivered to the main thread, so calls from other Python threads
/// are not interruptible, as before. A cancel file is polled as in the CLI.
fn interruptible<T: Send + 'static>(
    py: Python<'_>,
    cancel_file: Option<String>,
    work: impl FnOnce(&Control) -> T + Send + 'static,
) -> PyResult<T> {
    let mut control = Control::default();
    control.cancel_file = cancel_file.map(PathBuf::from);
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
                        control.cancel();
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

#[pyfunction]
fn describe_json() -> String {
    crn_craft::describe().to_string()
}
#[pyfunction]
fn analyze_json(model: &[u8]) -> PyResult<String> {
    serde_json::to_string(&crn_craft::analyze_bytes(model, None))
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
#[pyfunction]
fn check_json(model: &[u8], run: Option<&[u8]>) -> String {
    serde_json::to_string(&crn_craft::check_bytes(model, run)).unwrap()
}
#[pyfunction]
#[pyo3(signature = (model, run, cancel_file=None))]
fn simulate_json(
    py: Python<'_>,
    model: Vec<u8>,
    run: Vec<u8>,
    cancel_file: Option<String>,
) -> PyResult<String> {
    simulation_json(py, model, run, cancel_file, false)
}
fn simulation_json(
    py: Python<'_>,
    model: Vec<u8>,
    run: Vec<u8>,
    cancel_file: Option<String>,
    cached: bool,
) -> PyResult<String> {
    let prepared = py
        .detach(move || crn_craft::simulation::prepare(&model, &run))
        .map_err(|e| failure(&e))?;
    let result = interruptible(py, cancel_file, move |control| {
        if cached {
            crn_craft::simulation::execute_cached(&prepared, control)
        } else {
            crn_craft::simulation::execute(&prepared, control)
        }
    })?
    .map_err(|e| failure(&e))?;
    serde_json::to_string(&result.result).map_err(|e| PyValueError::new_err(e.to_string()))
}
#[pyfunction]
#[pyo3(signature = (model, run, cancel_file=None))]
fn fsp_json(
    py: Python<'_>,
    model: Vec<u8>,
    run: Vec<u8>,
    cancel_file: Option<String>,
) -> PyResult<String> {
    let result = interruptible(py, cancel_file, move |control| {
        crn_craft::fsp::solve(&model, &run, control)
    })?
    .map_err(|e| failure(&e))?;
    serde_json::to_string(&result).map_err(|e| PyValueError::new_err(e.to_string()))
}
#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(describe_json, module)?)?;
    module.add_function(wrap_pyfunction!(analyze_json, module)?)?;
    module.add_function(wrap_pyfunction!(check_json, module)?)?;
    module.add_function(wrap_pyfunction!(simulate_json, module)?)?;
    module.add_function(wrap_pyfunction!(fsp_json, module)?)?;
    module.add_function(wrap_pyfunction!(simulate_cached_json, module)?)?;
    module.add_function(wrap_pyfunction!(tau_json, module)?)?;
    module.add_function(wrap_pyfunction!(archive_json, module)?)?;
    module.add_function(wrap_pyfunction!(verify_json, module)?)?;
    module.add_function(wrap_pyfunction!(check_files_json, module)?)?;
    module.add_function(wrap_pyfunction!(schema_json, module)?)?;
    module.add_function(wrap_pyfunction!(export_json, module)?)?;
    module.add_function(wrap_pyfunction!(analyze_state_json, module)?)?;
    module.add_function(wrap_pyfunction!(tau_files_json, module)?)?;
    module.add("NativeError", module.py().get_type::<NativeError>())?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (model, run, cancel_file=None))]
fn simulate_cached_json(
    py: Python<'_>,
    model: Vec<u8>,
    run: Vec<u8>,
    cancel_file: Option<String>,
) -> PyResult<String> {
    simulation_json(py, model, run, cancel_file, true)
}
#[pyfunction]
#[pyo3(signature = (model, request, cancel_file=None))]
fn tau_json(
    py: Python<'_>,
    model: Vec<u8>,
    request: Vec<u8>,
    cancel_file: Option<String>,
) -> PyResult<String> {
    let result = interruptible(py, cancel_file, move |control| {
        crn_craft::tau::solve(&model, &request, control)
    })?
    .map_err(|e| failure(&e))?;
    serde_json::to_string(&result).map_err(|e| PyValueError::new_err(e.to_string()))
}
#[pyfunction]
#[pyo3(signature = (model, run, output, cached, cancel_file=None))]
fn archive_json(
    py: Python<'_>,
    model: String,
    run: String,
    output: String,
    cached: bool,
    cancel_file: Option<String>,
) -> PyResult<String> {
    interruptible(py, cancel_file, move |control| {
        let f = if cached {
            crn_craft::archive::simulate_cached_files
        } else {
            crn_craft::archive::simulate_files
        };
        serde_json::to_string(&f(
            std::path::Path::new(&model),
            std::path::Path::new(&run),
            std::path::Path::new(&output),
            control,
        ))
        .unwrap()
    })
}
#[pyfunction]
fn verify_json(py: Python<'_>, directory: String) -> String {
    py.detach(move || {
        serde_json::to_string(&crn_craft::archive::verify_archive(std::path::Path::new(
            &directory,
        )))
        .unwrap()
    })
}
#[pyfunction]
fn check_files_json(model: String, run: Option<String>) -> String {
    serde_json::to_string(&crn_craft::check_files(
        std::path::Path::new(&model),
        run.as_deref().map(std::path::Path::new),
    ))
    .unwrap()
}

#[pyfunction]
fn schema_json(kind: String) -> PyResult<String> {
    crn_craft::schema(&kind)
        .map(|v| v.to_string())
        .ok_or_else(|| PyValueError::new_err("unknown schema"))
}
#[pyfunction]
fn export_json(model: Vec<u8>) -> String {
    serde_json::to_string(&crn_craft::export::export_gillespy2(&model)).unwrap()
}
#[pyfunction]
fn analyze_state_json(model: Vec<u8>, state: Option<Vec<u8>>, legacy: bool) -> String {
    if legacy {
        serde_json::to_string(&crn_craft::analysis_v1::analyze_bytes(
            &model,
            state.as_deref(),
        ))
        .unwrap()
    } else {
        serde_json::to_string(&crn_craft::analyze_bytes(&model, state.as_deref())).unwrap()
    }
}
#[pyfunction]
#[pyo3(signature = (model, request, cancel_file=None))]
fn tau_files_json(
    py: Python<'_>,
    model: String,
    request: String,
    cancel_file: Option<String>,
) -> PyResult<String> {
    interruptible(
        py,
        cancel_file,
        move |control| match crn_craft::tau::solve_files(
            std::path::Path::new(&model),
            std::path::Path::new(&request),
            control,
        ) {
            Ok(result) => serde_json::to_string(&result).unwrap(),
            Err(errors) => serde_json::json!({"success":false,"diagnostics":errors}).to_string(),
        },
    )
}
