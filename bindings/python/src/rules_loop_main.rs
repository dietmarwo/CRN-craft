//! Runs one rules model with its generated Static Direct event loop.
//!
//! `crn_craft.rules.compile_loop` writes this file, `rules_payload.rs`, the
//! generated loop and the exact model bytes into a small Cargo project and
//! builds it against the `crn-craft` crate. The extension does not load the
//! program; Python runs it as a subprocess.
//!
//! Usage: `crn-craft-loop [--cancel-file PATH] < request.json`
//!
//! The program reads a rules-run request for the embedded model on standard
//! input and runs it with `execute_with_rules_loop`. On success it exits
//! with status 0 and writes the result as the extension returns it: the
//! header line `crn-craft-rules-loop/1`, then for the metadata and each
//! record column a frame of name length (u32, little-endian), name, data
//! length (u64, little-endian) and data. A rejected request exits with
//! status 2 and writes its diagnostics as JSON. Usage and I/O errors exit
//! with status 64 or 1 and a message on standard error.

mod generated_loop;
mod rules_payload;

use crn_craft::rules::{Control, execute_with_rules_loop, prepare};
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::ExitCode;

/// The exact model bytes the loop was generated for.
const MODEL: &[u8] = include_bytes!("../model.rules.json");
const HEADER: &[u8] = b"crn-craft-rules-loop/1\n";

fn rejected(diagnostics: serde_json::Result<String>) -> ExitCode {
    let text = diagnostics.unwrap_or_else(|_| "serialization failure".into());
    let mut out = std::io::stdout().lock();
    match out.write_all(text.as_bytes()).and_then(|()| out.flush()) {
        Ok(()) => ExitCode::from(2),
        Err(error) => {
            eprintln!("crn-craft-loop: cannot write diagnostics: {error}");
            ExitCode::FAILURE
        }
    }
}

fn write_frames(metadata: &str, columns: &[(&str, Vec<u8>)]) -> std::io::Result<()> {
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    out.write_all(HEADER)?;
    let frames = std::iter::once(("metadata", metadata.as_bytes()))
        .chain(columns.iter().map(|(name, data)| (*name, data.as_slice())));
    for (name, data) in frames {
        let size = u32::try_from(name.len()).expect("short column name");
        out.write_all(&size.to_le_bytes())?;
        out.write_all(name.as_bytes())?;
        out.write_all(&(data.len() as u64).to_le_bytes())?;
        out.write_all(data)?;
    }
    out.flush()
}

fn main() -> ExitCode {
    let mut control = Control::default();
    let mut args = std::env::args_os().skip(1);
    while let Some(arg) = args.next() {
        match (arg.to_str(), args.next()) {
            (Some("--cancel-file"), Some(path)) => control.cancel_file = Some(PathBuf::from(path)),
            _ => {
                eprintln!("usage: crn-craft-loop [--cancel-file PATH] < request.json");
                return ExitCode::from(64);
            }
        }
    }
    let mut request = Vec::new();
    if let Err(error) = std::io::stdin().read_to_end(&mut request) {
        eprintln!("crn-craft-loop: cannot read the request: {error}");
        return ExitCode::FAILURE;
    }
    let prepared = match prepare(MODEL, &request) {
        Ok(prepared) => prepared,
        Err(diagnostics) => return rejected(serde_json::to_string(&diagnostics)),
    };
    let kernel = generated_loop::GeneratedRulesLoop;
    let result = match execute_with_rules_loop(&prepared, &control, &kernel) {
        Ok(result) => result,
        Err(diagnostic) => return rejected(serde_json::to_string(&diagnostic)),
    };
    let (metadata, columns) = rules_payload::encode(result);
    match write_frames(&metadata, &columns) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("crn-craft-loop: cannot write the result: {error}");
            ExitCode::FAILURE
        }
    }
}
