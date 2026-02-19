mod reader;

pub use reader::read;

/// A dynamically-typed value produced by deserializing a Python pickle stream.
#[derive(Debug, Clone, PartialEq)]
pub enum PickleValue {
    None,
    Bool(bool),
    Int(i64),
    Float(f64),
    Bytes(Vec<u8>),
    String(String),
    List(Vec<PickleValue>),
    Tuple(Vec<PickleValue>),
    /// Key-value pairs preserving insertion order.
    Dict(Vec<(PickleValue, PickleValue)>),
}

impl PickleValue {
    pub fn as_int(&self) -> Option<i64> {
        match self {
            PickleValue::Int(v) => Some(*v),
            PickleValue::Bool(b) => Some(if *b { 1 } else { 0 }),
            PickleValue::Float(f) => Some(*f as i64),
            _ => None,
        }
    }

    pub fn as_float(&self) -> Option<f64> {
        match self {
            PickleValue::Float(v) => Some(*v),
            PickleValue::Int(v) => Some(*v as f64),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            PickleValue::String(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_bytes(&self) -> Option<&[u8]> {
        match self {
            PickleValue::Bytes(b) => Some(b),
            _ => None,
        }
    }

    pub fn as_list(&self) -> Option<&[PickleValue]> {
        match self {
            PickleValue::List(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_tuple(&self) -> Option<&[PickleValue]> {
        match self {
            PickleValue::Tuple(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_dict(&self) -> Option<&[(PickleValue, PickleValue)]> {
        match self {
            PickleValue::Dict(v) => Some(v),
            _ => None,
        }
    }

    /// Look up a key in a Dict by string key.
    pub fn get(&self, key: &str) -> Option<&PickleValue> {
        self.as_dict().and_then(|pairs| {
            pairs.iter().find_map(|(k, v)| {
                if k.as_str() == Some(key) {
                    Some(v)
                } else {
                    None
                }
            })
        })
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            PickleValue::Bool(b) => Some(*b),
            PickleValue::Int(i) => Some(*i != 0),
            _ => None,
        }
    }
}
