//! Structured diagnostics shared by every contract, solver and the rules engine.
//!
//! A [`Diagnostic`] carries a stable code (`CRN-E…` for errors, `CRN-W…` for
//! warnings), an RFC 6901 JSON Pointer into the offending document, a message
//! and a repair hint. Codes and pointers are the machine contract; messages
//! may be reworded. The full code table is in
//! [docs/reference/diagnostics.md](https://github.com/dietmarwo/CRN-craft/blob/main/docs/reference/diagnostics.md).
//!
//! ```
//! use crn_craft::diagnostic::{Diagnostic, Severity, pointer};
//!
//! let path = pointer("/parameters", "k/on");
//! assert_eq!(path, "/parameters/k~1on");
//! let error = Diagnostic::error("CRN-E015", &path, "non-finite parameter", "Use a finite number.");
//! assert_eq!(error.severity, Severity::Error);
//! assert_eq!(error.code, "CRN-E015");
//! ```
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

// Docs written as `#[doc = concat!(..)]` reach rustdoc but not the generated
// JSON Schema; see the note at the top of lib.rs.

/// Whether a diagnostic blocks use of the input or result: `error` does, `warning` is advisory.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    #[doc = concat!("`error`: the input cannot be used as given, or a computation failed.")]
    #[doc = concat!("Any error makes a check report invalid.")]
    Error,
    #[doc = concat!("`warning`: advisory, for example the structural growth warning `CRN-W001`.")]
    Warning,
}

/// One structured finding with a stable code, the JSON Pointer of the input it concerns, a
/// message and a repair hint.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Diagnostic {
    #[doc = concat!("Whether the finding blocks use of the input or result.")]
    pub severity: Severity,
    /// Stable code such as `CRN-E012` or `CRN-W001`, part of the machine contract.
    pub code: String,
    /// RFC 6901 JSON Pointer. Empty string denotes the document root.
    pub path: String,
    /// Human-readable detail; its wording is not part of the contract.
    pub message: String,
    /// Suggested repair for the finding.
    pub hint: String,
}

impl Diagnostic {
    /// Builds an error diagnostic. `path` is an RFC 6901 pointer, usually built with [`pointer()`].
    pub fn error(code: &str, path: &str, message: impl Into<String>, hint: &str) -> Self {
        Self {
            severity: Severity::Error,
            code: code.into(),
            path: path.into(),
            message: message.into(),
            hint: hint.into(),
        }
    }

    /// Builds a warning diagnostic; otherwise the same as [`Diagnostic::error`].
    pub fn warning(code: &str, path: &str, message: impl Into<String>, hint: &str) -> Self {
        let mut diagnostic = Self::error(code, path, message, hint);
        diagnostic.severity = Severity::Warning;
        diagnostic
    }
}

/// Appends one reference token to the JSON Pointer `parent`, escaping `~` as
/// `~0` and `/` as `~1` (RFC 6901). Pass `""` for the document root.
pub fn pointer(parent: &str, key: impl std::fmt::Display) -> String {
    format!(
        "{}/{}",
        parent,
        key.to_string().replace('~', "~0").replace('/', "~1")
    )
}
