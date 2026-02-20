//! BinaryCIF MessagePack schema types
//!
//! These types map directly to the bCIF MessagePack structure.
//! See: https://github.com/molstar/BinaryCIF

use serde::Deserialize;

/// Top-level bCIF file
#[derive(Deserialize)]
#[allow(dead_code)]
pub struct BcifFile {
    #[serde(deserialize_with = "string_or_bin::deserialize")]
    pub version: String,
    #[serde(deserialize_with = "string_or_bin::deserialize")]
    pub encoder: String,
    #[serde(rename = "dataBlocks")]
    pub data_blocks: Vec<BcifDataBlock>,
}

/// A data block (typically one per PDB entry)
#[derive(Deserialize)]
pub struct BcifDataBlock {
    #[serde(deserialize_with = "string_or_bin::deserialize")]
    pub header: String,
    pub categories: Vec<BcifCategory>,
}

/// A CIF category (e.g., _atom_site, _cell, _struct_conf)
#[derive(Deserialize)]
pub struct BcifCategory {
    #[serde(deserialize_with = "string_or_bin::deserialize")]
    pub name: String,
    #[serde(rename = "rowCount")]
    pub row_count: u32,
    pub columns: Vec<BcifColumn>,
}

/// A CIF column within a category
#[derive(Deserialize)]
pub struct BcifColumn {
    #[serde(deserialize_with = "string_or_bin::deserialize")]
    pub name: String,
    pub data: BcifData,
    pub mask: Option<BcifData>,
}

/// Encoded binary data blob
#[derive(Deserialize)]
pub struct BcifData {
    #[serde(with = "serde_bytes")]
    pub data: Vec<u8>,
    pub encoding: Vec<BcifEncoding>,
}

/// An encoding step in the encoding chain.
///
/// Encodings are applied in reverse order during decoding.
#[derive(Clone)]
#[allow(dead_code)]
pub enum BcifEncoding {
    ByteArray {
        data_type: i32,
    },
    FixedPoint {
        factor: f64,
        src_type: i32,
    },
    IntervalQuantization {
        min: f64,
        max: f64,
        num_steps: i32,
        src_type: i32,
    },
    RunLength {
        src_type: i32,
        src_size: i32,
    },
    Delta {
        origin: i64,
        src_type: i32,
    },
    IntegerPacking {
        byte_count: i32,
        src_size: i32,
        is_unsigned: bool,
    },
    StringArray {
        data_encoding: Vec<BcifEncoding>,
        string_data: String,
        offset_encoding: Vec<BcifEncoding>,
        offsets: Vec<u8>,
    },
}

/// Manual Deserialize for BcifEncoding.
///
/// RCSB's bCIF encoder uses MessagePack `bin` format for map keys and the
/// `kind` tag value. `rmp-serde`'s `#[serde(tag = "kind")]` can't match
/// bin-encoded tags, so we deserialize the map manually and dispatch on `kind`.
impl<'de> Deserialize<'de> for BcifEncoding {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        deserializer.deserialize_map(EncodingVisitor)
    }
}

/// Visitor that reads a map and dispatches based on the `kind` field.
///
/// Handles map keys encoded as either MessagePack `str` or `bin`.
struct EncodingVisitor;

impl<'de> serde::de::Visitor<'de> for EncodingVisitor {
    type Value = BcifEncoding;

    fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("a bCIF encoding map with a 'kind' field")
    }

    fn visit_map<A>(self, mut map: A) -> Result<BcifEncoding, A::Error>
    where
        A: serde::de::MapAccess<'de>,
    {
        use serde::de::Error;

        let mut kind: Option<String> = None;
        // Numeric fields
        let mut data_type: Option<i32> = None;
        let mut factor: Option<f64> = None;
        let mut src_type: Option<i32> = None;
        let mut min: Option<f64> = None;
        let mut max: Option<f64> = None;
        let mut num_steps: Option<i32> = None;
        let mut src_size: Option<i32> = None;
        let mut origin: Option<i64> = None;
        let mut byte_count: Option<i32> = None;
        let mut is_unsigned: Option<bool> = None;
        // StringArray fields
        let mut data_encoding: Option<Vec<BcifEncoding>> = None;
        let mut string_data: Option<String> = None;
        let mut offset_encoding: Option<Vec<BcifEncoding>> = None;
        let mut offsets: Option<Vec<u8>> = None;

        while let Some(key) = map.next_key::<StringOrBin>()? {
            match key.0.as_str() {
                "kind" => kind = Some(map.next_value::<StringOrBin>()?.0),
                "type" => data_type = Some(map.next_value()?),
                "factor" => factor = Some(map.next_value()?),
                "srcType" => src_type = Some(map.next_value()?),
                "min" => min = Some(map.next_value()?),
                "max" => max = Some(map.next_value()?),
                "numSteps" => num_steps = Some(map.next_value()?),
                "srcSize" => src_size = Some(map.next_value()?),
                "origin" => origin = Some(map.next_value()?),
                "byteCount" => byte_count = Some(map.next_value()?),
                "isUnsigned" => is_unsigned = Some(map.next_value()?),
                "dataEncoding" => data_encoding = Some(map.next_value()?),
                "stringData" => string_data = Some(map.next_value::<StringOrBin>()?.0),
                "offsetEncoding" => offset_encoding = Some(map.next_value()?),
                "offsets" => offsets = Some(map.next_value()?),
                _ => {
                    // Skip unknown fields
                    map.next_value::<serde::de::IgnoredAny>()?;
                }
            }
        }

        let kind = kind.ok_or_else(|| A::Error::missing_field("kind"))?;

        match kind.as_str() {
            "ByteArray" => Ok(BcifEncoding::ByteArray {
                data_type: data_type.ok_or_else(|| A::Error::missing_field("type"))?,
            }),
            "FixedPoint" => Ok(BcifEncoding::FixedPoint {
                factor: factor.ok_or_else(|| A::Error::missing_field("factor"))?,
                src_type: src_type.ok_or_else(|| A::Error::missing_field("srcType"))?,
            }),
            "IntervalQuantization" => Ok(BcifEncoding::IntervalQuantization {
                min: min.ok_or_else(|| A::Error::missing_field("min"))?,
                max: max.ok_or_else(|| A::Error::missing_field("max"))?,
                num_steps: num_steps.ok_or_else(|| A::Error::missing_field("numSteps"))?,
                src_type: src_type.ok_or_else(|| A::Error::missing_field("srcType"))?,
            }),
            "RunLength" => Ok(BcifEncoding::RunLength {
                src_type: src_type.ok_or_else(|| A::Error::missing_field("srcType"))?,
                src_size: src_size.ok_or_else(|| A::Error::missing_field("srcSize"))?,
            }),
            "Delta" => Ok(BcifEncoding::Delta {
                origin: origin.ok_or_else(|| A::Error::missing_field("origin"))?,
                src_type: src_type.ok_or_else(|| A::Error::missing_field("srcType"))?,
            }),
            "IntegerPacking" => Ok(BcifEncoding::IntegerPacking {
                byte_count: byte_count.ok_or_else(|| A::Error::missing_field("byteCount"))?,
                src_size: src_size.ok_or_else(|| A::Error::missing_field("srcSize"))?,
                is_unsigned: is_unsigned
                    .ok_or_else(|| A::Error::missing_field("isUnsigned"))?,
            }),
            "StringArray" => Ok(BcifEncoding::StringArray {
                data_encoding: data_encoding
                    .ok_or_else(|| A::Error::missing_field("dataEncoding"))?,
                string_data: string_data
                    .ok_or_else(|| A::Error::missing_field("stringData"))?,
                offset_encoding: offset_encoding
                    .ok_or_else(|| A::Error::missing_field("offsetEncoding"))?,
                offsets: offsets.ok_or_else(|| A::Error::missing_field("offsets"))?,
            }),
            other => Err(A::Error::unknown_variant(
                other,
                &[
                    "ByteArray",
                    "FixedPoint",
                    "IntervalQuantization",
                    "RunLength",
                    "Delta",
                    "IntegerPacking",
                    "StringArray",
                ],
            )),
        }
    }
}

/// Newtype for deserializing a String from either MessagePack `str` or `bin`.
struct StringOrBin(String);

impl<'de> Deserialize<'de> for StringOrBin {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        string_or_bin::deserialize(deserializer).map(StringOrBin)
    }
}

/// Helper module for deserializing MessagePack str or bin as Rust String.
///
/// Some bCIF encoders (e.g., RCSB's `python-mmcif library`) encode string
/// fields using MessagePack binary (bin) format instead of string (str) format.
/// This visitor accepts both.
mod string_or_bin {
    use serde::de;

    pub fn deserialize<'de, D>(deserializer: D) -> Result<String, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        struct StringOrBytes;

        impl<'de> de::Visitor<'de> for StringOrBytes {
            type Value = String;

            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("a string or byte array")
            }

            fn visit_str<E: de::Error>(self, v: &str) -> Result<String, E> {
                Ok(v.to_owned())
            }

            fn visit_string<E: de::Error>(self, v: String) -> Result<String, E> {
                Ok(v)
            }

            fn visit_bytes<E: de::Error>(self, v: &[u8]) -> Result<String, E> {
                String::from_utf8(v.to_vec()).map_err(E::custom)
            }

            fn visit_byte_buf<E: de::Error>(self, v: Vec<u8>) -> Result<String, E> {
                String::from_utf8(v).map_err(E::custom)
            }
        }

        deserializer.deserialize_any(StringOrBytes)
    }
}

/// Helper module for deserializing MessagePack binary data as Vec<u8>.
///
/// MessagePack binary (bin) format must be deserialized with serde_bytes
/// to avoid being interpreted as a sequence of integers.
mod serde_bytes {
    use serde::Deserialize;

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Vec<u8>, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        // rmp-serde deserializes MessagePack bin as bytes when target is Vec<u8>
        Vec::<u8>::deserialize(deserializer)
    }
}
