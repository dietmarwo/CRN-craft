//! Minimal JSON bridge; numerical semantics and validators live in the core.
use pyo3::{exceptions::PyValueError, prelude::*};
fn failure<T: serde::Serialize>(value: &T) -> PyErr {
    PyValueError::new_err(
        serde_json::to_string(value).unwrap_or_else(|_| "serialization failure".into()),
    )
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
fn simulate_json(py: Python<'_>, model: Vec<u8>, run: Vec<u8>) -> PyResult<String> {
    py.detach(move || {
        let prepared = crn_craft::simulation::prepare(&model, &run).map_err(|e| failure(&e))?;
        let result =
            crn_craft::simulation::execute(&prepared, &crn_craft::simulation::Control::default())
                .map_err(|e| failure(&e))?;
        serde_json::to_string(&result.result).map_err(|e| PyValueError::new_err(e.to_string()))
    })
}
#[pyfunction]
fn fsp_json(py: Python<'_>, model: Vec<u8>, run: Vec<u8>) -> PyResult<String> {
    py.detach(move || {
        let result =
            crn_craft::fsp::solve(&model, &run, &crn_craft::simulation::Control::default())
                .map_err(|e| failure(&e))?;
        serde_json::to_string(&result).map_err(|e| PyValueError::new_err(e.to_string()))
    })
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
    Ok(())
}

#[pyfunction]
fn simulate_cached_json(py: Python<'_>, model: Vec<u8>, run: Vec<u8>) -> PyResult<String> {
    py.detach(move || {
        let prepared = crn_craft::simulation::prepare(&model, &run).map_err(|e| failure(&e))?;
        let result = crn_craft::simulation::execute_cached(
            &prepared,
            &crn_craft::simulation::Control::default(),
        )
        .map_err(|e| failure(&e))?;
        serde_json::to_string(&result.result).map_err(|e| PyValueError::new_err(e.to_string()))
    })
}
#[pyfunction]
fn tau_json(py: Python<'_>, model: Vec<u8>, request: Vec<u8>) -> PyResult<String> {
    py.detach(move || {
        let result =
            crn_craft::tau::solve(&model, &request, &crn_craft::simulation::Control::default())
                .map_err(|e| failure(&e))?;
        serde_json::to_string(&result).map_err(|e| PyValueError::new_err(e.to_string()))
    })
}
#[pyfunction]
fn archive_json(
    py: Python<'_>,
    model: String,
    run: String,
    output: String,
    cached: bool,
    cancel_file: Option<String>,
) -> String {
    py.detach(move || {
        let mut control = crn_craft::simulation::Control::default();
        control.cancel_file = cancel_file.map(std::path::PathBuf::from);
        let f = if cached {
            crn_craft::archive::simulate_cached_files
        } else {
            crn_craft::archive::simulate_files
        };
        serde_json::to_string(&f(
            std::path::Path::new(&model),
            std::path::Path::new(&run),
            std::path::Path::new(&output),
            &control,
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
fn tau_files_json(py: Python<'_>, model: String, request: String) -> String {
    py.detach(move || {
        match crn_craft::tau::solve_files(
            std::path::Path::new(&model),
            std::path::Path::new(&request),
            &crn_craft::simulation::Control::default(),
        ) {
            Ok(result) => serde_json::to_string(&result).unwrap(),
            Err(errors) => serde_json::json!({"success":false,"diagnostics":errors}).to_string(),
        }
    })
}
