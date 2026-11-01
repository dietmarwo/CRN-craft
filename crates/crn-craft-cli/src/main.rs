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

const HELP: &str = "crn-craft — M1b network analysis and propensities\n\nUsage:\n  crn-craft describe [model|run|result|check|state|analysis|analysis-v1|export]\n  crn-craft check MODEL [RUN]\n  crn-craft analyze MODEL [STATE]\n  crn-craft analyze-v1 MODEL [STATE]\n  crn-craft export MODEL gillespy2\n  crn-craft --version\n\nAll data output is JSON. check/analyze exit 0 on success and 1 for errors.\nAnalysis reports conservation limits explicitly. No simulation is implemented.\n";

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
            &serde_json::json!({"version": env!("CARGO_PKG_VERSION"), "milestone": "M1b"}),
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
        "Use --help. M1b supports describe, check, analyze, analyze-v1 and export; simulation is not implemented.",
    ));
    (
        serde_json::to_value(report).expect("finite usage report"),
        false,
    )
}
