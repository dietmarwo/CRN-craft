//! Compact JSON with object keys in sorted order, written without building a
//! `serde_json::Value`.
//!
//! [`to_vec`] gives the bytes of `serde_json::to_writer(&serde_json::to_value(x)?)`:
//! `Value` keeps objects in a `BTreeMap<String, Value>`, so its keys come out
//! in byte order of the strings, and a key inserted twice keeps the later
//! value. Scalars are written by `serde_json`'s own serializer, after the
//! conversions `to_value` applies (an `f32` widened to `f64`, a non-finite
//! number as `null`, map keys as strings). Each object's members are
//! serialized into one buffer and copied out in key order, so the cost is
//! one buffer per object instead of a tree of maps and strings.
use serde::ser::{self, Serialize};
use std::borrow::Cow;
use std::fmt::Display;

/// Serialization failure (unsupported map keys or a value's own error).
#[derive(Debug)]
pub struct Error(String);

impl Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

impl ser::Error for Error {
    fn custom<T: Display>(message: T) -> Self {
        Self(message.to_string())
    }
}

impl From<serde_json::Error> for Error {
    fn from(error: serde_json::Error) -> Self {
        Self(error.to_string())
    }
}

/// The compact JSON of `value` with sorted object keys.
pub fn to_vec<T: Serialize + ?Sized>(value: &T) -> Result<Vec<u8>, Error> {
    let mut out = Vec::with_capacity(1 << 16);
    value.serialize(Writer { out: &mut out })?;
    Ok(out)
}

/// Writes one scalar with `serde_json`'s compact serializer.
fn scalar<T: Serialize + ?Sized>(out: &mut Vec<u8>, value: &T) -> Result<(), Error> {
    value.serialize(&mut serde_json::Serializer::new(out))?;
    Ok(())
}

struct Writer<'a> {
    out: &'a mut Vec<u8>,
}

impl<'a> ser::Serializer for Writer<'a> {
    type Ok = ();
    type Error = Error;
    type SerializeSeq = Seq<'a>;
    type SerializeTuple = Seq<'a>;
    type SerializeTupleStruct = Seq<'a>;
    type SerializeTupleVariant = Seq<'a>;
    type SerializeMap = Object<'a>;
    type SerializeStruct = Object<'a>;
    type SerializeStructVariant = Object<'a>;

    fn serialize_bool(self, v: bool) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_i8(self, v: i8) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_i16(self, v: i16) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_i32(self, v: i32) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_i64(self, v: i64) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_i128(self, v: i128) -> Result<(), Error> {
        // `Value` holds 128-bit integers only within the 64-bit range.
        if let Ok(v) = i64::try_from(v) {
            scalar(self.out, &v)
        } else if let Ok(v) = u64::try_from(v) {
            scalar(self.out, &v)
        } else {
            Err(ser::Error::custom("number out of range"))
        }
    }
    fn serialize_u8(self, v: u8) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_u16(self, v: u16) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_u32(self, v: u32) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_u64(self, v: u64) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_u128(self, v: u128) -> Result<(), Error> {
        u64::try_from(v)
            .map_err(|_| ser::Error::custom("number out of range"))
            .and_then(|v| scalar(self.out, &v))
    }
    fn serialize_f32(self, v: f32) -> Result<(), Error> {
        // `Value::from(f32)` widens to f64.
        self.serialize_f64(f64::from(v))
    }
    fn serialize_f64(self, v: f64) -> Result<(), Error> {
        if v.is_finite() {
            scalar(self.out, &v)
        } else {
            self.out.extend_from_slice(b"null");
            Ok(())
        }
    }
    fn serialize_char(self, v: char) -> Result<(), Error> {
        scalar(self.out, &v)
    }
    fn serialize_str(self, v: &str) -> Result<(), Error> {
        scalar(self.out, v)
    }
    fn serialize_bytes(self, v: &[u8]) -> Result<(), Error> {
        // `Value` holds bytes as an array of numbers.
        scalar(
            self.out,
            &v.iter().map(|&b| u64::from(b)).collect::<Vec<_>>(),
        )
    }
    fn serialize_none(self) -> Result<(), Error> {
        self.serialize_unit()
    }
    fn serialize_some<T: Serialize + ?Sized>(self, value: &T) -> Result<(), Error> {
        value.serialize(self)
    }
    fn serialize_unit(self) -> Result<(), Error> {
        self.out.extend_from_slice(b"null");
        Ok(())
    }
    fn serialize_unit_struct(self, _name: &'static str) -> Result<(), Error> {
        self.serialize_unit()
    }
    fn serialize_unit_variant(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
    ) -> Result<(), Error> {
        scalar(self.out, variant)
    }
    fn serialize_newtype_struct<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        value: &T,
    ) -> Result<(), Error> {
        value.serialize(self)
    }
    fn serialize_newtype_variant<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
        value: &T,
    ) -> Result<(), Error> {
        // `{variant: value}`, an object with one member.
        self.out.push(b'{');
        scalar(self.out, variant)?;
        self.out.push(b':');
        value.serialize(Writer {
            out: &mut *self.out,
        })?;
        self.out.push(b'}');
        Ok(())
    }
    fn serialize_seq(self, _len: Option<usize>) -> Result<Seq<'a>, Error> {
        self.out.push(b'[');
        Ok(Seq {
            out: self.out,
            first: true,
            variant: false,
        })
    }
    fn serialize_tuple(self, len: usize) -> Result<Seq<'a>, Error> {
        self.serialize_seq(Some(len))
    }
    fn serialize_tuple_struct(self, _name: &'static str, len: usize) -> Result<Seq<'a>, Error> {
        self.serialize_seq(Some(len))
    }
    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
        _len: usize,
    ) -> Result<Seq<'a>, Error> {
        self.out.push(b'{');
        scalar(self.out, variant)?;
        self.out.extend_from_slice(b":[");
        Ok(Seq {
            out: self.out,
            first: true,
            variant: true,
        })
    }
    fn serialize_map(self, _len: Option<usize>) -> Result<Object<'a>, Error> {
        Ok(Object::new(self.out, None))
    }
    fn serialize_struct(self, _name: &'static str, _len: usize) -> Result<Object<'a>, Error> {
        Ok(Object::new(self.out, None))
    }
    fn serialize_struct_variant(
        self,
        _name: &'static str,
        _index: u32,
        variant: &'static str,
        _len: usize,
    ) -> Result<Object<'a>, Error> {
        Ok(Object::new(self.out, Some(variant)))
    }
}

struct Seq<'a> {
    out: &'a mut Vec<u8>,
    first: bool,
    /// A tuple variant: `{variant: [...]}`.
    variant: bool,
}

impl Seq<'_> {
    fn element<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        if !self.first {
            self.out.push(b',');
        }
        self.first = false;
        value.serialize(Writer {
            out: &mut *self.out,
        })
    }

    fn close(self) {
        self.out.push(b']');
        if self.variant {
            self.out.push(b'}');
        }
    }
}

impl ser::SerializeSeq for Seq<'_> {
    type Ok = ();
    type Error = Error;
    fn serialize_element<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.element(value)
    }
    fn end(self) -> Result<(), Error> {
        self.close();
        Ok(())
    }
}

impl ser::SerializeTuple for Seq<'_> {
    type Ok = ();
    type Error = Error;
    fn serialize_element<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.element(value)
    }
    fn end(self) -> Result<(), Error> {
        self.close();
        Ok(())
    }
}

impl ser::SerializeTupleStruct for Seq<'_> {
    type Ok = ();
    type Error = Error;
    fn serialize_field<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.element(value)
    }
    fn end(self) -> Result<(), Error> {
        self.close();
        Ok(())
    }
}

impl ser::SerializeTupleVariant for Seq<'_> {
    type Ok = ();
    type Error = Error;
    fn serialize_field<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        self.element(value)
    }
    fn end(self) -> Result<(), Error> {
        self.close();
        Ok(())
    }
}

/// An object whose members are serialized into `buffer` (each as
/// `"key":value`) and copied to `out` in key order at the end.
struct Object<'a> {
    out: &'a mut Vec<u8>,
    /// A struct variant: `{variant: {...}}`.
    variant: Option<&'static str>,
    buffer: Vec<u8>,
    /// Each member's key and its `"key":value` range in `buffer`, in
    /// insertion order.
    members: Vec<(Cow<'static, str>, std::ops::Range<usize>)>,
    /// The key of the map member whose value comes next.
    pending: Option<Cow<'static, str>>,
}

impl<'a> Object<'a> {
    fn new(out: &'a mut Vec<u8>, variant: Option<&'static str>) -> Self {
        Self {
            out,
            variant,
            buffer: Vec::with_capacity(256),
            members: Vec::with_capacity(16),
            pending: None,
        }
    }

    fn member<T: Serialize + ?Sized>(
        &mut self,
        key: Cow<'static, str>,
        value: &T,
    ) -> Result<(), Error> {
        let start = self.buffer.len();
        scalar(&mut self.buffer, key.as_ref())?;
        self.buffer.push(b':');
        value.serialize(Writer {
            out: &mut self.buffer,
        })?;
        self.members.push((key, start..self.buffer.len()));
        Ok(())
    }

    fn close(self) -> Result<(), Error> {
        let Self {
            out,
            variant,
            buffer,
            mut members,
            ..
        } = self;
        if let Some(variant) = variant {
            out.push(b'{');
            scalar(out, variant)?;
            out.push(b':');
        }
        // Stable: of two members with one key, the later one follows.
        if !members.is_sorted_by(|a, b| a.0 <= b.0) {
            members.sort_by(|a, b| a.0.cmp(&b.0));
        }
        out.reserve(buffer.len() + members.len() + 2);
        out.push(b'{');
        let mut first = true;
        for (k, (key, range)) in members.iter().enumerate() {
            if members.get(k + 1).is_some_and(|next| next.0 == *key) {
                // A later insert of the same key replaces this one.
                continue;
            }
            if !first {
                out.push(b',');
            }
            first = false;
            out.extend_from_slice(&buffer[range.clone()]);
        }
        out.push(b'}');
        if variant.is_some() {
            out.push(b'}');
        }
        Ok(())
    }
}

impl ser::SerializeStruct for Object<'_> {
    type Ok = ();
    type Error = Error;
    fn serialize_field<T: Serialize + ?Sized>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Error> {
        self.member(Cow::Borrowed(key), value)
    }
    fn end(self) -> Result<(), Error> {
        self.close()
    }
}

impl ser::SerializeStructVariant for Object<'_> {
    type Ok = ();
    type Error = Error;
    fn serialize_field<T: Serialize + ?Sized>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Error> {
        self.member(Cow::Borrowed(key), value)
    }
    fn end(self) -> Result<(), Error> {
        self.close()
    }
}

impl ser::SerializeMap for Object<'_> {
    type Ok = ();
    type Error = Error;
    fn serialize_key<T: Serialize + ?Sized>(&mut self, key: &T) -> Result<(), Error> {
        // `Value`'s map keys: strings, or numbers, booleans and chars as
        // their text.
        let text = match serde_json::to_value(key)? {
            serde_json::Value::String(text) => text,
            serde_json::Value::Number(n) => n.to_string(),
            serde_json::Value::Bool(b) => b.to_string(),
            _ => return Err(ser::Error::custom("key must be a string")),
        };
        self.pending = Some(Cow::Owned(text));
        Ok(())
    }
    fn serialize_value<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Error> {
        let key = self.pending.take().expect("a key before its value");
        self.member(key, value)
    }
    fn end(self) -> Result<(), Error> {
        self.close()
    }
}

#[cfg(test)]
mod tests {
    use super::to_vec;
    use serde::Serialize;
    use serde_json::json;

    fn reference<T: Serialize>(value: &T) -> Vec<u8> {
        serde_json::to_vec(&serde_json::to_value(value).unwrap()).unwrap()
    }

    #[derive(Serialize)]
    enum Shape {
        Unit,
        New(f64),
        Tuple(u8, f32),
        Record { zeta: i64, alpha: Option<String> },
    }

    #[derive(Serialize)]
    struct Mixed {
        zulu: Vec<Shape>,
        alpha: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        skipped: Option<u8>,
        mike: std::collections::BTreeMap<u32, (bool, char)>,
        bravo: [f64; 3],
        unit: (),
        #[serde(flatten)]
        rest: std::collections::HashMap<String, i32>,
        bytes: serde_bytes_like::Bytes,
    }

    mod serde_bytes_like {
        /// Serialized with `serialize_bytes`.
        pub struct Bytes(pub Vec<u8>);
        impl serde::Serialize for Bytes {
            fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                s.serialize_bytes(&self.0)
            }
        }
    }

    #[test]
    fn sorted_bytes_equal_value_bytes() {
        let mixed = Mixed {
            zulu: vec![
                Shape::Unit,
                Shape::New(-0.0),
                Shape::New(f64::NAN),
                Shape::Tuple(7, 0.1),
                Shape::Record {
                    zeta: -3,
                    alpha: Some("é \"q\" \n \u{1}".into()),
                },
                Shape::Record {
                    zeta: i64::MIN,
                    alpha: None,
                },
            ],
            alpha: 1e-310,
            skipped: None,
            mike: [(10, (true, 'x')), (9, (false, '\u{7f}'))].into(),
            bravo: [f64::INFINITY, 2.5e300, 1.0],
            unit: (),
            rest: [("aaa".into(), 1), ("zzz".into(), 2)].into(),
            bytes: serde_bytes_like::Bytes(vec![0, 255]),
        };
        assert_eq!(to_vec(&mixed).unwrap(), reference(&mixed));
        for value in [
            json!({"b": [1, {"d": null, "c": [true, "x"]}], "a": -0.0, "": {}, "A": []}),
            json!([[], {}, 0.30000000000000004, u64::MAX, "\u{2028}"]),
        ] {
            assert_eq!(to_vec(&value).unwrap(), reference(&value));
        }
    }

    #[test]
    fn a_repeated_key_keeps_the_later_value() {
        struct Twice;
        impl Serialize for Twice {
            fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                use serde::ser::SerializeMap;
                let mut map = s.serialize_map(None)?;
                map.serialize_entry("k", &1)?;
                map.serialize_entry("a", &2)?;
                map.serialize_entry("k", &3)?;
                map.serialize_entry(&5u8, &4)?;
                map.end()
            }
        }
        assert_eq!(to_vec(&Twice).unwrap(), reference(&Twice));
        assert_eq!(to_vec(&Twice).unwrap(), br#"{"5":4,"a":2,"k":3}"#);
    }
}
