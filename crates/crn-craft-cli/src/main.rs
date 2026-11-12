#![forbid(unsafe_code)]
use crn_craft::{
    CheckReport, analyze_files, check_files, describe, diagnostic::Diagnostic,
    export_gillespy2_file, schema,
};
use serde_json::Value;
use std::{
    env,
    io::{self, Write},
    path::Path,
    process::ExitCode,
};

const HELP: &str = "crn-craft — CRN solvers and network analysis\n\nUsage:\n  crn-craft describe [model|run|result|result-v1|archive|archive-check|check|state|analysis-1.0|export|fsp-result|tau-request|tau-result|ssa-run|ssa-result|rules-model|rules-check|rules-run|rules-result|model-1.0|check-1.0|run-1.0|result-1.0|fsp-result-1.0|archive-1.0|archive-check-1.0]\n  crn-craft check MODEL [RUN]\n  crn-craft analyze MODEL [STATE]\n  crn-craft export MODEL gillespy2\n  crn-craft simulate MODEL RUN NEW_DIRECTORY [--cancel-file PATH]\n  crn-craft verify DIRECTORY\n  crn-craft simulate-cached MODEL RUN NEW_DIRECTORY [--cancel-file PATH]\n  crn-craft fsp MODEL RUN [--cancel-file PATH]\n  crn-craft tau MODEL REQUEST [--cancel-file PATH]\n  crn-craft ssa MODEL REQUEST [--cancel-file PATH]\n  crn-craft rules-check RULES_MODEL\n  crn-craft rules RULES_MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]\n  crn-craft convert MODEL\n  crn-craft --version\n\nanalyze and export read crn-craft.model/0.1 (through convert), rules-model/0.1-0.4\nand model/1.0; analyze prints analysis/1.0. The other network commands read\ncrn-craft.model/0.1 and refuse crn-craft.model/1.0 (CRN-E032); rules-check and\nrules read rules-model/0.1-0.4 and model/1.0.\nrules solves a run/1.0 request whose method is fsp and prints fsp-result/1.0.\nrules --archive writes the model, request and result with a manifest\n(archive/1.0) and prints an archive-check/1.0 report; it exits 0 only for a\ncompleted run. verify reads archive/0.1 and archive/1.0 directories.\nAll data output is JSON; convert prints the converted rules model. simulate exits 0 only for a completed ensemble.\nverify exits 0 for a valid archive, including an incomplete simulation.\nOther errors and incomplete simulations exit 1. Optional ODE uses crn-diffsol.\n";

mod sorted;

fn emit(value: &Value) -> io::Result<()> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    serde_json::to_writer_pretty(&mut out, value)?;
    out.write_all(b"\n")?;
    out.flush()
}

/// FSP and tau results are checked against max_output_bytes in compact form,
/// so they are written compactly.
fn emit_compact(value: &Value) -> io::Result<()> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    serde_json::to_writer(&mut out, value)?;
    out.write_all(b"\n")?;
    out.flush()
}

/// Prints the bytes [`emit_compact`] prints for `serde_json::to_value(result)`
/// (compact JSON with sorted keys) without building the value.
fn emit_sorted(result: &impl serde::Serialize) -> io::Result<()> {
    let bytes = sorted::to_vec(result).map_err(io::Error::other)?;
    let mut out = io::stdout().lock();
    out.write_all(&bytes)?;
    out.write_all(b"\n")
}

fn main() -> ExitCode {
    let args: Vec<_> = env::args_os().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        return match io::stdout().write_all(HELP.as_bytes()) {
            Ok(()) => ExitCode::SUCCESS,
            Err(_) => ExitCode::FAILURE,
        };
    }
    if args.len() == 1 && args[0] == "--version" {
        return match emit(
            &serde_json::json!({"version": env!("CARGO_PKG_VERSION"), "milestone": "SSA-acceleration"}),
        ) {
            Ok(()) => ExitCode::SUCCESS,
            Err(_) => ExitCode::FAILURE,
        };
    }
    if args.len() == 2 && args[0] == "convert" {
        return convert(Path::new(&args[1]));
    }
    if args.first().is_some_and(|arg| arg == "rules")
        && args.len() >= 3
        && rules_options(&args[3..]).is_some()
    {
        return rules(&args);
    }
    let mut compact = false;
    let (value, success) = match args.first().and_then(|arg| arg.to_str()) {
        Some("describe") if args.len() == 1 => {
            // The native CLI adds the rules engine to the core description.
            let mut value = describe();
            value["schema"] = "crn-craft.describe/0.15".into();
            if let Some(commands) = value["commands"].as_array_mut() {
                commands.extend(["rules-check".into(), "rules".into(), "convert".into()]);
            }
            value["dynamic_species"] = crn_craft_rules::describe();
            (value, true)
        }
        Some("describe") if args.len() == 2 => match args[1]
            .to_str()
            .and_then(|kind| schema(kind).or_else(|| crn_craft_rules::schema(kind)))
        {
            Some(value) => (value, true),
            None => usage_error(),
        },
        Some("check") if args.len() == 2 || args.len() == 3 => {
            let report = check_files(Path::new(&args[1]), args.get(2).map(Path::new));
            let success = report.valid;
            (
                serde_json::to_value(report).expect("finite check report"),
                success,
            )
        }
        Some("analyze") if args.len() == 2 || args.len() == 3 => {
            let report = analyze_files(Path::new(&args[1]), args.get(2).map(Path::new));
            let success = report.success;
            (
                serde_json::to_value(report).expect("finite analysis report"),
                success,
            )
        }
        Some("export") if args.len() == 3 && args[2] == "gillespy2" => {
            let report = export_gillespy2_file(Path::new(&args[1]));
            let success = report.success;
            (
                serde_json::to_value(report).expect("finite export report"),
                success,
            )
        }
        Some("simulate" | "simulate-cached")
            if args.len() == 4 || (args.len() == 6 && args[4] == "--cancel-file") =>
        {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(5).map(std::path::PathBuf::from);
            let runner = if args[0] == "simulate-cached" {
                crn_craft::archive::simulate_cached_files
            } else {
                crn_craft::archive::simulate_files
            };
            let report = runner(
                Path::new(&args[1]),
                Path::new(&args[2]),
                Path::new(&args[3]),
                &control,
            );
            let success = report.valid && report.simulation_complete;
            (
                serde_json::to_value(report).expect("finite archive report"),
                success,
            )
        }
        Some("ssa") if args.len() == 3 || (args.len() == 5 && args[3] == "--cancel-file") => {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(4).map(std::path::PathBuf::from);
            match crn_craft::ssa::solve_files(Path::new(&args[1]), Path::new(&args[2]), &control) {
                Ok(result) => {
                    let complete = result.termination == crn_craft::ssa::Termination::Completed;
                    compact = true;
                    (serde_json::to_value(result).unwrap(), complete)
                }
                Err(errors) => (
                    serde_json::json!({"success":false,"diagnostics":errors}),
                    false,
                ),
            }
        }
        Some("rules-check") if args.len() == 2 => {
            let report = crn_craft_rules::check_model_file(Path::new(&args[1]));
            let success = report.valid;
            (
                serde_json::to_value(report).expect("finite rules check report"),
                success,
            )
        }
        Some("tau") if args.len() == 3 || (args.len() == 5 && args[3] == "--cancel-file") => {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(4).map(std::path::PathBuf::from);
            match crn_craft::tau::solve_files(Path::new(&args[1]), Path::new(&args[2]), &control) {
                Ok(result) => {
                    let complete = result.termination == crn_craft::simulation::Stop::Completed;
                    compact = true;
                    (serde_json::to_value(result).unwrap(), complete)
                }
                Err(errors) => (
                    serde_json::json!({"success":false,"diagnostics":errors}),
                    false,
                ),
            }
        }
        Some("fsp") if args.len() == 3 || (args.len() == 5 && args[3] == "--cancel-file") => {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(4).map(std::path::PathBuf::from);
            let check = check_files(Path::new(&args[1]), Some(Path::new(&args[2])));
            if !check.valid {
                (serde_json::to_value(check).unwrap(), false)
            } else {
                match (std::fs::read(&args[1]), std::fs::read(&args[2])) {
                    (Ok(m), Ok(r)) => match crn_craft::fsp::solve(&m, &r, &control) {
                        Ok(result) => {
                            compact = true;
                            (serde_json::to_value(result).unwrap(), true)
                        }
                        Err(errors) => (
                            serde_json::json!({"success":false,"diagnostics":errors}),
                            false,
                        ),
                    },
                    _ => usage_error(),
                }
            }
        }
        Some("verify") if args.len() == 2 => {
            // archive/1.0 or the network engine's archive/0.1.
            let report = crn_craft::archive::verify_any(Path::new(&args[1]));
            let success = report.valid();
            (
                serde_json::to_value(report).expect("finite archive report"),
                success,
            )
        }
        _ => usage_error(),
    };
    let written = if compact {
        emit_compact(&value)
    } else {
        emit(&value)
    };
    match written {
        Ok(()) if success => ExitCode::SUCCESS,
        Ok(()) => ExitCode::FAILURE,
        Err(error) => {
            eprintln!("CRN-E030: cannot write JSON output: {error}");
            ExitCode::FAILURE
        }
    }
}

/// `rules MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]`: a
/// rules or run/1.0 simulation, or a run/1.0 `fsp` request, printed as compact
/// JSON with sorted keys. With `--archive` the archive report goes to standard
/// output and the result into the new directory.
fn rules(args: &[std::ffi::OsString]) -> ExitCode {
    let (archive, cancel_file) = rules_options(&args[3..]).expect("checked options");
    let mut control = crn_craft_rules::Control::default();
    control.cancel_file = cancel_file.map(std::path::PathBuf::from);
    let (model, request) = (Path::new(&args[1]), Path::new(&args[2]));
    let (written, success) = if let Some(directory) = archive {
        let report =
            crn_craft::archive::v1::run_files(model, request, Path::new(directory), &control);
        let success = report.valid && report.simulation_complete;
        (
            emit(&serde_json::to_value(report).expect("finite archive report")),
            success,
        )
    } else if crn_craft::fsp::v1::is_fsp_request_file(request) {
        // A run/1.0 request whose method is `fsp`.
        match crn_craft::fsp::v1::solve_files(model, request, &control) {
            Ok(result) => (emit_sorted(&result), true),
            Err(errors) => (
                emit(&serde_json::json!({"success":false,"diagnostics":errors})),
                false,
            ),
        }
    } else {
        match crn_craft_rules::solve_files(model, request, &control) {
            Ok(result) => (
                emit_sorted(&result),
                result.termination == crn_craft_rules::Termination::Completed,
            ),
            Err(errors) => (
                emit(&serde_json::json!({"success":false,"diagnostics":errors})),
                false,
            ),
        }
    };
    match written {
        Ok(()) if success => ExitCode::SUCCESS,
        Ok(()) => ExitCode::FAILURE,
        Err(error) => {
            eprintln!("CRN-E030: cannot write JSON output: {error}");
            ExitCode::FAILURE
        }
    }
}

/// Prints the converted model's exact bytes, or the refusal as JSON.
fn convert(path: &Path) -> ExitCode {
    let (written, success) = match crn_craft::rules::convert::from_network_file(path) {
        Ok(bytes) => (io::stdout().lock().write_all(&bytes), true),
        Err(errors) => (
            emit(&serde_json::json!({"success": false, "diagnostics": errors})),
            false,
        ),
    };
    match written {
        Ok(()) if success => ExitCode::SUCCESS,
        Ok(()) => ExitCode::FAILURE,
        Err(error) => {
            eprintln!("CRN-E030: cannot write JSON output: {error}");
            ExitCode::FAILURE
        }
    }
}

/// The options after `rules MODEL REQUEST`: `--archive NEW_DIRECTORY` and
/// `--cancel-file PATH`, each at most once, in any order. `None` for anything
/// else.
fn rules_options(
    options: &[std::ffi::OsString],
) -> Option<(Option<&std::ffi::OsStr>, Option<&std::ffi::OsStr>)> {
    let (mut archive, mut cancel_file) = (None, None);
    for pair in options.chunks(2) {
        let [name, value] = pair else {
            return None;
        };
        let slot = match name.to_str() {
            Some("--archive") => &mut archive,
            Some("--cancel-file") => &mut cancel_file,
            _ => return None,
        };
        if slot.replace(value.as_os_str()).is_some() {
            return None;
        }
    }
    Some((archive, cancel_file))
}

fn usage_error() -> (Value, bool) {
    let report = CheckReport::from_error(Diagnostic::error(
        "CRN-E031",
        "",
        "unknown command or wrong argument count",
        "Use --help. Supported: describe, check, analyze, export, simulate, simulate-cached, fsp, tau, ssa, rules-check, rules, convert and verify.",
    ));
    (
        serde_json::to_value(report).expect("finite usage report"),
        false,
    )
}
