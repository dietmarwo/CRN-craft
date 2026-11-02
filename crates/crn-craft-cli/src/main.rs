#![forbid(unsafe_code)]
use crn_craft::{
    CheckReport, analyze_files, analyze_legacy_files, check_files, describe,
    diagnostic::Diagnostic, export_gillespy2_file, schema,
};
use serde_json::Value;
use std::{
    env,
    io::{self, Write},
    path::Path,
    process::ExitCode,
};

const HELP: &str = "crn-craft — CRN solvers and network analysis\n\nUsage:\n  crn-craft describe [model|run|result|result-v1|archive|archive-check|check|state|analysis|analysis-v1|export|fsp-result|tau-request|tau-result|ssa-run|ssa-result|rules-model|rules-check|rules-run|rules-result]\n  crn-craft check MODEL [RUN]\n  crn-craft analyze MODEL [STATE]\n  crn-craft analyze-v1 MODEL [STATE]\n  crn-craft export MODEL gillespy2\n  crn-craft simulate MODEL RUN NEW_DIRECTORY [--cancel-file PATH]\n  crn-craft verify DIRECTORY\n  crn-craft simulate-cached MODEL RUN NEW_DIRECTORY [--cancel-file PATH]\n  crn-craft fsp MODEL RUN [--cancel-file PATH]\n  crn-craft tau MODEL REQUEST [--cancel-file PATH]\n  crn-craft ssa MODEL REQUEST [--cancel-file PATH]\n  crn-craft rules-check RULES_MODEL\n  crn-craft rules RULES_MODEL REQUEST [--cancel-file PATH]\n  crn-craft --version\n\nAll data output is JSON. simulate exits 0 only for a completed ensemble.\nverify exits 0 for a valid archive, including an incomplete simulation.\nOther errors and incomplete simulations exit 1. Optional ODE uses crn-diffsol.\n";

fn emit(value: &Value) -> io::Result<()> {
    let mut out = io::stdout().lock();
    serde_json::to_writer_pretty(&mut out, value)?;
    out.write_all(b"\n")
}

/// FSP and tau results are checked against max_output_bytes in compact form,
/// so they are written compactly.
fn emit_compact(value: &Value) -> io::Result<()> {
    let mut out = io::stdout().lock();
    serde_json::to_writer(&mut out, value)?;
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
    let mut compact = false;
    let (value, success) = match args.first().and_then(|arg| arg.to_str()) {
        Some("describe") if args.len() == 1 => (describe(), true),
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
        Some("analyze-v1") if args.len() == 2 || args.len() == 3 => {
            let report = analyze_legacy_files(Path::new(&args[1]), args.get(2).map(Path::new));
            let success = report.success;
            (
                serde_json::to_value(report).expect("finite legacy report"),
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
        Some("rules") if args.len() == 3 || (args.len() == 5 && args[3] == "--cancel-file") => {
            let mut control = crn_craft_rules::Control::default();
            control.cancel_file = args.get(4).map(std::path::PathBuf::from);
            match crn_craft_rules::solve_files(Path::new(&args[1]), Path::new(&args[2]), &control) {
                Ok(result) => {
                    let complete = result.termination == crn_craft_rules::Termination::Completed;
                    compact = true;
                    (serde_json::to_value(result).unwrap(), complete)
                }
                Err(errors) => (
                    serde_json::json!({"success":false,"diagnostics":errors}),
                    false,
                ),
            }
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
            let report = crn_craft::archive::verify_archive(Path::new(&args[1]));
            let success = report.valid;
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

fn usage_error() -> (Value, bool) {
    let report = CheckReport::from_error(Diagnostic::error(
        "CRN-E031",
        "",
        "unknown command or wrong argument count",
        "Use --help. Supported: describe, check, analyze, analyze-v1, export, simulate, simulate-cached, fsp, tau and verify.",
    ));
    (
        serde_json::to_value(report).expect("finite usage report"),
        false,
    )
}
