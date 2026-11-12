//! Sealed archives: model, request and result with a checksummed manifest
//! and a completion marker written last.
//!
//! Durable-write pattern inspected in optimization-jobs storage.rs at cc6b03d.
//! [`v1`] writes and verifies `crn-craft.archive/1.0` directories of rules
//! and `crn-craft.run/1.0` runs; `crn-craft run MODEL REQUEST --archive
//! NEW_DIRECTORY` and `crn-craft verify DIRECTORY` call it. Archive failures
//! are `CRN-E065`. The network engine's `crn-craft.archive/0.1` was removed
//! in 0.7.0 with that engine.
//!
//! A complete archive can hold an incomplete simulation, so read
//! `archive_complete` and `simulation_complete` together. Integrity checks
//! are not signatures and do not establish biological validity. See
//! [docs/archive-contract.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/archive-contract.md).
use crate::diagnostic::Diagnostic;
use serde::Serialize;
use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::Path,
};

pub mod v1;

fn error(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E065",
        "/archive",
        message,
        "Use a new output directory; preserve all files and verify the completion marker and checksums.",
    )
}
fn io_error(e: std::io::Error) -> Diagnostic {
    error(e.to_string())
}

struct LimitedBuffer {
    bytes: Vec<u8>,
    limit: u64,
}
impl Write for LimitedBuffer {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if self.bytes.len() as u64 + bytes.len() as u64 > self.limit {
            return Err(std::io::Error::other(
                "serialized output byte budget exhausted",
            ));
        }
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encoded<T: Serialize>(value: &T, limit: u64) -> Result<Vec<u8>, Diagnostic> {
    let mut buffer = LimitedBuffer {
        bytes: Vec::new(),
        limit,
    };
    serde_json::to_writer(&mut buffer, value).map_err(|e| error(e.to_string()))?;
    buffer.write_all(b"\n").map_err(io_error)?;
    Ok(buffer.bytes)
}
fn write_new(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(io_error)?;
    file.write_all(bytes).map_err(io_error)?;
    file.sync_all().map_err(io_error)
}
fn read_regular(path: &Path, limit: u64) -> Result<Vec<u8>, Diagnostic> {
    let metadata = fs::symlink_metadata(path).map_err(io_error)?;
    if !metadata.is_file() || metadata.len() > limit {
        return Err(error("not a bounded regular archive file"));
    }
    let file = File::open(path).map_err(io_error)?;
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(io_error)?;
    if bytes.len() as u64 > limit {
        return Err(error("archive file grew beyond its bound"));
    }
    Ok(bytes)
}
