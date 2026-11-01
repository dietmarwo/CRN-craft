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

const HELP: &str = "crn-craft — M2 Direct SSA and network analysis\n\nUsage:\n  crn-craft describe [model|run|result|result-v1|archive|archive-check|check|state|analysis|analysis-v1|export]\n  crn-craft check MODEL [RUN]\n  crn-craft analyze MODEL [STATE]\n  crn-craft analyze-v1 MODEL [STATE]\n  crn-craft export MODEL gillespy2\n  crn-craft simulate MODEL RUN NEW_DIRECTORY [--cancel-file PATH]\n  crn-craft verify DIRECTORY\n  crn-craft --version\n\nAll data output is JSON. simulate exits 0 only for a completed ensemble.\nverify exits 0 for a valid archive, including an incomplete simulation.\nOther errors and incomplete simulations exit 1. FSP is not implemented.\n";

fn emit(value: &Value) -> io::Result<()> {
    let mut out = io::stdout().lock();
    serde_json::to_writer_pretty(&mut out, value)?;
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
            &serde_json::json!({"version": env!("CARGO_PKG_VERSION"), "milestone": "M2"}),
        ) {
            Ok(()) => ExitCode::SUCCESS,
            Err(_) => ExitCode::FAILURE,
        };
    }
    let (value, success) = match args.first().and_then(|arg| arg.to_str()) {
        Some("describe") if args.len() == 1 => (describe(), true),
        Some("describe") if args.len() == 2 => match args[1].to_str().and_then(schema) {
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
        Some("simulate") if args.len() == 4 || (args.len() == 6 && args[4] == "--cancel-file") => {
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = args.get(5).map(std::path::PathBuf::from);
            let report = crn_craft::archive::simulate_files(
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
    match emit(&value) {
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
        "Use --help. M2 supports describe, check, analyze, analyze-v1, export, simulate and verify; FSP is not implemented.",
    ));
    (
        serde_json::to_value(report).expect("finite usage report"),
        false,
    )
}
