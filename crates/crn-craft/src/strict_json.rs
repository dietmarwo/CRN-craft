//! Reject duplicate keys before serde can overwrite them in maps.
//!
//! [`parse`] is the strict JSON reader behind every contract in this crate and
//! in `crn-craft-rules`. It reads one UTF-8 document of at most
//! [`MAX_DOCUMENT_BYTES`] and [`MAX_JSON_DEPTH`] levels, rejects duplicate keys
//! at any depth, then deserializes into `T`, so unknown fields are rejected
//! wherever `T` denies them. Failures are a single [`Diagnostic`] whose path is
//! an RFC 6901 pointer:
//!
//! | Code | Cause |
//! |---|---|
//! | `CRN-E001` | Invalid JSON (`NaN` and `Infinity` included), invalid UTF-8 or trailing content |
//! | `CRN-E002` | Duplicate key, also when two escaped spellings decode to the same key |
//! | `CRN-E003` | Wrong type, missing or unknown field, unknown version or variant |
//! | `CRN-E004` | Document larger than 16 MiB or nested deeper than 48 levels |
//! | `CRN-E015` | Numeric literal out of range, such as `1e400` |
//!
//! Parsing is not semantic validation; see [`crate::check_any`].
//!
//! ```
//! use crn_craft::strict_json::parse;
//! use std::collections::BTreeMap;
//!
//! let duplicate = parse::<BTreeMap<String, f64>>(br#"{"k": 1, "\u006b": 2}"#).unwrap_err();
//! assert_eq!((duplicate.code.as_str(), duplicate.path.as_str()), ("CRN-E002", "/k"));
//! let ok: BTreeMap<String, f64> = parse(br#"{"k": 1.5}"#).unwrap();
//! assert_eq!(ok["k"], 1.5);
//! ```
use crate::diagnostic::{Diagnostic, pointer};
use serde::de::{self, DeserializeOwned, DeserializeSeed, MapAccess, SeqAccess, Visitor};
use serde_json::{Map, Number, Value};
use std::fmt;

/// Largest accepted input document: 16 MiB, whitespace included (`CRN-E004`).
pub const MAX_DOCUMENT_BYTES: usize = 16 * 1024 * 1024;
/// Deepest accepted nesting level; the root value is level 0 (`CRN-E004`).
pub const MAX_JSON_DEPTH: usize = 48;

struct Seed<'a> {
    path: String,
    depth: usize,
    failure: &'a mut Option<Diagnostic>,
}

impl<'de> DeserializeSeed<'de> for Seed<'_> {
    type Value = Value;
    fn deserialize<D: de::Deserializer<'de>>(self, de: D) -> Result<Value, D::Error> {
        let path = self.path.clone();
        if self.depth > MAX_JSON_DEPTH {
            *self.failure = Some(Diagnostic::error(
                "CRN-E004",
                &path,
                "JSON nesting limit exceeded",
                "Keep every value within 48 levels below the root, which is level 0.",
            ));
            return Err(de::Error::custom("JSON nesting limit exceeded"));
        }
        let failure = &mut *self.failure;
        let result = de.deserialize_any(JsonVisitor(Seed {
            path: self.path,
            depth: self.depth,
            failure,
        }));
        if let Err(ref error) = result {
            if failure.is_none() {
                let message = error.to_string();
                let code = if message.contains("number out of range") {
                    "CRN-E015"
                } else {
                    "CRN-E001"
                };
                *failure = Some(Diagnostic::error(
                    code,
                    &path,
                    message,
                    "Use valid JSON with finite numbers.",
                ));
            }
        }
        result
    }
}

struct JsonVisitor<'a>(Seed<'a>);

impl<'de> Visitor<'de> for JsonVisitor<'_> {
    type Value = Value;
    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a JSON value")
    }
    fn visit_bool<E: de::Error>(self, value: bool) -> Result<Value, E> {
        Ok(Value::Bool(value))
    }
    fn visit_i64<E: de::Error>(self, value: i64) -> Result<Value, E> {
        Ok(value.into())
    }
    fn visit_u64<E: de::Error>(self, value: u64) -> Result<Value, E> {
        Ok(value.into())
    }
    fn visit_f64<E: de::Error>(self, value: f64) -> Result<Value, E> {
        Number::from_f64(value)
            .map(Value::Number)
            .ok_or_else(|| E::custom("non-finite number"))
    }
    fn visit_str<E: de::Error>(self, value: &str) -> Result<Value, E> {
        Ok(Value::String(value.into()))
    }
    fn visit_string<E: de::Error>(self, value: String) -> Result<Value, E> {
        Ok(Value::String(value))
    }
    fn visit_unit<E: de::Error>(self) -> Result<Value, E> {
        Ok(Value::Null)
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut access: A) -> Result<Value, A::Error> {
        let mut values = Vec::new();
        while let Some(value) = access.next_element_seed(Seed {
            path: pointer(&self.0.path, values.len()),
            depth: self.0.depth + 1,
            failure: self.0.failure,
        })? {
            values.push(value);
        }
        Ok(Value::Array(values))
    }
    fn visit_map<A: MapAccess<'de>>(self, mut access: A) -> Result<Value, A::Error> {
        let mut values = Map::new();
        while let Some(key) = access.next_key::<String>()? {
            let path = pointer(&self.0.path, &key);
            if values.contains_key(&key) {
                *self.0.failure = Some(Diagnostic::error(
                    "CRN-E002",
                    &path,
                    format!("duplicate JSON key {key:?}"),
                    "Keep one value per key, including in stoichiometry and parameter maps.",
                ));
                return Err(de::Error::custom("duplicate JSON key"));
            }
            let value = access.next_value_seed(Seed {
                path,
                depth: self.0.depth + 1,
                failure: self.0.failure,
            })?;
            values.insert(key, value);
        }
        Ok(Value::Object(values))
    }
}

/// Parses one strict JSON document into `T`; see the [module docs](self) for
/// the limits and diagnostic codes. Stops at the first failure.
pub fn parse<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, Diagnostic> {
    parse_bounded(bytes, MAX_DOCUMENT_BYTES)
}

// Archive results have a separate encoded-size cap; input contracts keep 16 MiB.
pub(crate) fn parse_bounded<T: DeserializeOwned>(
    bytes: &[u8],
    max_bytes: usize,
) -> Result<T, Diagnostic> {
    if bytes.len() > max_bytes {
        return Err(Diagnostic::error(
            "CRN-E004",
            "",
            if max_bytes == MAX_DOCUMENT_BYTES {
                "document exceeds 16 MiB".into()
            } else {
                format!("document exceeds {max_bytes} bytes")
            },
            "Split or reduce the model.",
        ));
    }
    let mut deserializer = serde_json::Deserializer::from_slice(bytes);
    let mut failure = None;
    let value = Seed {
        path: String::new(),
        depth: 0,
        failure: &mut failure,
    }
    .deserialize(&mut deserializer)
    .map_err(|error| {
        failure.unwrap_or_else(|| {
            Diagnostic::error(
                "CRN-E001",
                "",
                error.to_string(),
                "Use one valid UTF-8 JSON document.",
            )
        })
    })?;
    deserializer.end().map_err(|error| {
        Diagnostic::error("CRN-E001", "", error.to_string(), "Remove trailing data.")
    })?;
    serde_path_to_error::deserialize(value).map_err(|error| {
        let mut path = String::new();
        for segment in error.path() {
            use serde_path_to_error::Segment;
            match segment {
                Segment::Seq { index } => path = pointer(&path, index),
                Segment::Map { key } => path = pointer(&path, key),
                Segment::Enum { .. } | Segment::Unknown => {}
            }
        }
        Diagnostic::error(
            "CRN-E003",
            &path,
            error.inner().to_string(),
            "Compare this field with the schema emitted by describe.",
        )
    })
}
