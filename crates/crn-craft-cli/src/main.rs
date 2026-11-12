#![forbid(unsafe_code)]
use crn_craft::{
    CheckReport, analyze_files, check_any_files, describe, diagnostic::Diagnostic,
    export_gillespy2_file, schema, unified::MILESTONE,
};
use serde_json::Value;
use std::{
    env,
    io::{self, Write},
    path::Path,
    process::ExitCode,
};

const HELP: &str =
    "crn-craft — reaction networks: exact and approximate simulation, distributions, analysis

Usage:
  crn-craft describe [KIND]
  crn-craft check MODEL [REQUEST]
  crn-craft run MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]
  crn-craft verify DIRECTORY
  crn-craft analyze MODEL [STATE]
  crn-craft convert MODEL
  crn-craft export MODEL gillespy2
  crn-craft --version

Aliases: rules = run, rules-check = check.

Models are crn-craft.model/1.0 or crn-craft.rules-model/0.1-0.4; check, analyze,
convert and export also read a network crn-craft.model/0.1 (through the
converter; convert prints it as crn-craft.model/1.0). Requests are
crn-craft.run/1.0 (a method fsp prints crn-craft.fsp-result/1.0) or
crn-craft.rules-run/0.1-0.6; run prints the result as compact JSON.
check reports whether run would accept MODEL and REQUEST, without running.
run --archive writes the model, request and result with a manifest
(crn-craft.archive/1.0) and prints a crn-craft.archive-check/1.0 report.
verify reads crn-craft.archive/1.0 and 0.1 directories.

describe KIND: model-1.0 run-1.0 result-1.0 check-1.0 fsp-result-1.0
  analysis-1.0 archive-1.0 archive-check-1.0 export state; still read:
  rules-model rules-run rules-result rules-check model archive archive-check.

All output is JSON. run exits 0 only for a completed ensemble (with
--archive, a valid archive of one); verify exits 0 for a valid archive,
also of an incomplete run; check, analyze and export exit 0 for valid input.
Other errors exit 1.
";

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
            &serde_json::json!({"version": env!("CARGO_PKG_VERSION"), "milestone": MILESTONE}),
        ) {
            Ok(()) => ExitCode::SUCCESS,
            Err(_) => ExitCode::FAILURE,
        };
    }
    if args.len() == 2 && args[0] == "convert" {
        return convert(Path::new(&args[1]));
    }
    // `rules` is the 0.6.0 name of `run`.
    if args
        .first()
        .is_some_and(|arg| arg == "run" || arg == "rules")
        && args.len() >= 3
        && run_options(&args[3..]).is_some()
    {
        return run(&args);
    }
    let mut compact = false;
    let (value, success) = match args.first().and_then(|arg| arg.to_str()) {
        Some("describe") if args.len() == 1 => (describe(), true),
        Some("describe") if args.len() == 2 => match args[1]
            .to_str()
            .and_then(|kind| schema(kind).or_else(|| crn_craft::rules::schema(kind)))
        {
            Some(value) => (value, true),
            None => usage_error(),
        },
        // `rules-check` is the 0.6.0 name of `check`.
        Some("check" | "rules-check") if args.len() == 2 || args.len() == 3 => {
            let report = check_any_files(Path::new(&args[1]), args.get(2).map(Path::new));
            let success = report.valid();
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
        Some("verify") if args.len() == 2 => {
            // archive/1.0 or the network engine's archive/0.1.
            let report = crn_craft::archive::verify_any(Path::new(&args[1]));
            let success = report.valid();
            (
                serde_json::to_value(report).expect("finite archive report"),
                success,
            )
        }
        Some(command) if LEGACY.contains(&command) => match legacy(command, &args) {
            Some((value, success, is_compact)) => {
                compact = is_compact;
                (value, success)
            }
            None => usage_error(),
        },
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

/// The network engine's commands of release 0.6.0. They are hidden (absent
/// from `--help`, the usage message and `describe`) and keep working until
/// task M8b deletes them with the old engine. LEGACY: delete in M8b.
const LEGACY: [&str; 5] = ["simulate", "simulate-cached", "ssa", "tau", "fsp"];

/// Runs one hidden network command: `(output, success, compact)`, or `None`
/// for wrong arguments. LEGACY: delete in M8b.
fn legacy(command: &str, args: &[std::ffi::OsString]) -> Option<(Value, bool, bool)> {
    match command {
        "simulate" | "simulate-cached"
            if args.len() == 4 || (args.len() == 6 && args[4] == "--cancel-file") =>
        {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(5).map(std::path::PathBuf::from);
            let runner = if command == "simulate-cached" {
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
            Some((
                serde_json::to_value(report).expect("finite archive report"),
                success,
                false,
            ))
        }
        "ssa" if args.len() == 3 || (args.len() == 5 && args[3] == "--cancel-file") => {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(4).map(std::path::PathBuf::from);
            Some(
                match crn_craft::ssa::solve_files(
                    Path::new(&args[1]),
                    Path::new(&args[2]),
                    &control,
                ) {
                    Ok(result) => {
                        let complete = result.termination == crn_craft::ssa::Termination::Completed;
                        (serde_json::to_value(result).unwrap(), complete, true)
                    }
                    Err(errors) => (
                        serde_json::json!({"success":false,"diagnostics":errors}),
                        false,
                        false,
                    ),
                },
            )
        }
        "tau" if args.len() == 3 || (args.len() == 5 && args[3] == "--cancel-file") => {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(4).map(std::path::PathBuf::from);
            Some(
                match crn_craft::tau::solve_files(
                    Path::new(&args[1]),
                    Path::new(&args[2]),
                    &control,
                ) {
                    Ok(result) => {
                        let complete = result.termination == crn_craft::simulation::Stop::Completed;
                        (serde_json::to_value(result).unwrap(), complete, true)
                    }
                    Err(errors) => (
                        serde_json::json!({"success":false,"diagnostics":errors}),
                        false,
                        false,
                    ),
                },
            )
        }
        "fsp" if args.len() == 3 || (args.len() == 5 && args[3] == "--cancel-file") => {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(4).map(std::path::PathBuf::from);
            let check = crn_craft::check_files(Path::new(&args[1]), Some(Path::new(&args[2])));
            if !check.valid {
                return Some((serde_json::to_value(check).unwrap(), false, false));
            }
            let (Ok(model), Ok(run)) = (std::fs::read(&args[1]), std::fs::read(&args[2])) else {
                return None;
            };
            Some(match crn_craft::fsp::solve(&model, &run, &control) {
                Ok(result) => (serde_json::to_value(result).unwrap(), true, true),
                Err(errors) => (
                    serde_json::json!({"success":false,"diagnostics":errors}),
                    false,
                    false,
                ),
            })
        }
        _ => None,
    }
}

/// `run MODEL REQUEST [--archive NEW_DIRECTORY] [--cancel-file PATH]` (alias
/// `rules`): a rules or run/1.0 simulation, or a run/1.0 `fsp` request,
/// printed as compact JSON with sorted keys. With `--archive` the archive
/// report goes to standard output and the result into the new directory.
fn run(args: &[std::ffi::OsString]) -> ExitCode {
    let (archive, cancel_file) = run_options(&args[3..]).expect("checked options");
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

/// The options after `run MODEL REQUEST`: `--archive NEW_DIRECTORY` and
/// `--cancel-file PATH`, each at most once, in any order. `None` for anything
/// else.
fn run_options(
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
        "Use --help. Commands: describe, check, run, verify, analyze, convert and export (rules and rules-check are aliases of run and check).",
    ));
    (
        serde_json::to_value(report).expect("finite usage report"),
        false,
    )
}
