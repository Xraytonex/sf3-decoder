//! Scaffold's block-based object files (`BBProteinSequence`,
//! `SpectrumCollation`, ...): a Java-serialization-like dialect.
//!
//! The descriptor stream names the classes (`cm.classPos`) and where each
//! record starts (`idToLocation`: part file and byte offset). In a record:
//! `0x73` object, `0x75` array, `0x70` null, `0x74` string with a
//! little-endian u16 length (`0x7c`: i64 length), `0x7e` an enum constant
//! (class, then its name as a string); a class reference is `0x71` + a little-endian
//! i32 (1-based index into the class list). Object fields are big-endian in
//! descriptor order, superclass first; array lengths and primitive arrays
//! are little-endian. Classes that write custom data do so little-endian
//! with no end marker, so each needs its own reader here.

use std::collections::HashMap;
use std::sync::Arc;

use super::java::{self, ArrayData, SC_WRITE_METHOD, Stream};
use super::records::Container;
use super::{Sf3Error, layout};

const MAX_DEPTH: usize = 256;
/// Values one file may decode to; far above any real project's, far below
/// what exhausts memory.
const MAX_VALUES: usize = 5_000_000;

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum BbValue {
    Null,
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(String),
    Object(BbObject),
    Array(Vec<BbValue>),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct BbObject {
    pub class: String,
    pub fields: Vec<(Arc<str>, BbValue)>,
    /// Custom data, flattened: map entries as key, value, key, value...
    pub custom: Vec<BbValue>,
}

impl BbObject {
    pub fn field(&self, name: &str) -> Option<&BbValue> {
        self.fields
            .iter()
            .find(|(field, _)| field.as_ref() == name)
            .map(|(_, value)| value)
    }

    pub fn int(&self, name: &str) -> Option<i64> {
        match self.field(name)? {
            BbValue::Int(value) => Some(*value),
            _ => None,
        }
    }

    pub fn str(&self, name: &str) -> Option<&str> {
        match self.field(name)? {
            BbValue::Str(value) => Some(value),
            _ => None,
        }
    }
}

/// Every record of a block-based file, parsed. A record that fails to parse
/// fails the whole file: a partial read would be a guess.
pub(crate) fn read(container: &Container, name: &str) -> Result<Vec<BbObject>, Sf3Error> {
    let descriptor = container.stream(name)?;
    let (classes, locations) = describe(&descriptor)
        .ok_or_else(|| layout(name, "descriptor lacks its class list or record index"))?;
    let mut parts: HashMap<i32, Vec<u8>> = HashMap::new();
    let mut records = Vec::with_capacity(locations.len());
    let mut budget = MAX_VALUES;
    let mut seen = std::collections::HashSet::new();
    for (file, offset) in locations {
        if (file, offset) == (-1, -1) {
            continue; // an unused id
        }
        if file < 0 || offset < 0 || !seen.insert((file, offset)) {
            return Err(layout(
                name,
                &format!("bad record location {file}:{offset}"),
            ));
        }
        let data = match parts.entry(file) {
            std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
            std::collections::hash_map::Entry::Vacant(entry) => entry.insert(
                container
                    .file(&format!("{name}{file}"))?
                    .ok_or_else(|| layout(name, &format!("part {file} is missing")))?,
            ),
        };
        let data: &[u8] = data;
        let offset = usize::try_from(offset).unwrap_or(usize::MAX);
        // Some files put the record id and a continuation before the object.
        let header = (data.get(offset) != Some(&0x73))
            .then(|| super::records::read_i32(data, offset))
            .flatten();
        let mut reader = Reader {
            data,
            pos: offset.saturating_add(if header.is_some() { 8 } else { 0 }),
            stream: &descriptor,
            classes: &classes,
            depth: 0,
            budget,
        };
        let result = reader.value();
        budget = reader.budget;
        match result {
            Ok(BbValue::Object(object))
                if header.is_none_or(|id| object.int("id") == Some(i64::from(id))) =>
            {
                records.push(object);
            }
            Ok(_) => {
                return Err(layout(
                    name,
                    &format!("record at {offset} is not an object with the id its header gives"),
                ));
            }
            Err(problem) => {
                return Err(layout(name, &format!("record at {offset}: {problem}")));
            }
        }
    }
    Ok(records)
}

/// The class list and each record's (part, offset).
/// Class descriptor indices, then each record's (part, byte offset).
type Description = (Vec<usize>, Vec<(i32, i32)>);

fn describe(stream: &Stream) -> Option<Description> {
    let serialization = *stream.objects_of("BlockBasedSerization").first()?;
    let class_map = stream.object(serialization, "cm")?;
    let positions = stream.object(class_map, "classPos")?;
    let classes: Vec<usize> = stream.objects[positions]
        .annotations
        .iter()
        .flat_map(|(_, values)| values)
        .filter_map(|value| match value {
            java::Value::ClassDesc(index) => Some(*index),
            _ => None,
        })
        .collect();
    let index = stream.object(serialization, "idToLocation")?;
    let (ArrayData::Ints(files), ArrayData::Ints(offsets)) = (
        stream.array(index, "file")?,
        stream.array(index, "location")?,
    ) else {
        return None;
    };
    (files.len() == offsets.len()).then(|| {
        (
            classes,
            files.iter().copied().zip(offsets.iter().copied()).collect(),
        )
    })
}

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
    stream: &'a Stream,
    classes: &'a [usize],
    depth: usize,
    /// Values still allowed for this file.
    budget: usize,
}

impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8], String> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.data.len())
            .ok_or_else(|| format!("ends early at byte {}", self.pos))?;
        let bytes = &self.data[self.pos..end];
        self.pos = end;
        Ok(bytes)
    }

    fn le_i32(&mut self) -> Result<i32, String> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().expect("four")))
    }

    fn length(&mut self, element_size: usize) -> Result<usize, String> {
        let length = usize::try_from(self.le_i32()?).map_err(|_| "negative length".to_string())?;
        if length.saturating_mul(element_size) > self.data.len() - self.pos {
            return Err(format!("length {length} runs past the record"));
        }
        self.spend(length)?;
        Ok(length)
    }

    fn spend(&mut self, values: usize) -> Result<(), String> {
        self.budget = self
            .budget
            .checked_sub(values)
            .ok_or("more values than any real project holds")?;
        Ok(())
    }

    fn class(&mut self) -> Result<Option<usize>, String> {
        match self.take(1)?[0] {
            0x70 => Ok(None),
            0x71 => {
                let id = self.le_i32()?;
                usize::try_from(id)
                    .ok()
                    .and_then(|id| id.checked_sub(1))
                    .and_then(|index| self.classes.get(index))
                    .map(|&class| Some(class))
                    .ok_or_else(|| format!("unknown class {id}"))
            }
            other => Err(format!("expected a class reference, got {other:#x}")),
        }
    }

    fn value(&mut self) -> Result<BbValue, String> {
        self.spend(1)?;
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return Err("nesting too deep".into());
        }
        let value = match self.take(1)?[0] {
            0x70 => BbValue::Null,
            0x74 => {
                let length =
                    usize::from(u16::from_le_bytes(self.take(2)?.try_into().expect("two")));
                BbValue::Str(String::from_utf8_lossy(self.take(length)?).into_owned())
            }
            0x7c => {
                // A long string: little-endian i64 length (e.g. titin's sequence).
                let length = i64::from_le_bytes(self.take(8)?.try_into().expect("eight"));
                let length = usize::try_from(length)
                    .ok()
                    .filter(|&length| length <= self.data.len() - self.pos)
                    .ok_or_else(|| format!("long string of {length} bytes runs past the record"))?;
                BbValue::Str(String::from_utf8_lossy(self.take(length)?).into_owned())
            }
            0x75 => self.array()?,
            0x73 => self.object()?,
            // An enum constant: its class, then its name (`GOAnnotation`'s
            // evidence type, "IMP", in 5.3.3 protein records).
            0x7e => {
                self.class()?.ok_or("enum without a class")?;
                match self.value()? {
                    name @ BbValue::Str(_) => name,
                    _ => return Err("enum constant without a name".into()),
                }
            }
            other => return Err(format!("unknown tag {other:#x} at byte {}", self.pos - 1)),
        };
        self.depth -= 1;
        Ok(value)
    }

    fn array(&mut self) -> Result<BbValue, String> {
        let class = self.class()?.ok_or("array without a class")?;
        let element = self.stream.classes[class].name.as_bytes().get(1).copied();
        let primitive = |code: u8| -> Option<usize> {
            Some(match code {
                b'B' | b'Z' => 1,
                b'S' | b'C' => 2,
                b'I' | b'F' => 4,
                b'J' | b'D' => 8,
                _ => return None,
            })
        };
        match element {
            Some(code) if primitive(code).is_some() => {
                let size = primitive(code).expect("checked");
                let length = self.length(size)?;
                let bytes = self.take(length * size)?.to_vec();
                Ok(BbValue::Array(
                    bytes
                        .chunks_exact(size)
                        .map(|chunk| le_primitive(code, chunk))
                        .collect(),
                ))
            }
            _ => {
                let length = self.length(1)?;
                (0..length)
                    .map(|_| self.value())
                    .collect::<Result<_, _>>()
                    .map(BbValue::Array)
            }
        }
    }

    fn object(&mut self) -> Result<BbValue, String> {
        let class = self.class()?.ok_or("object without a class")?;
        let mut chain = vec![class];
        while let Some(parent) = self.stream.classes[*chain.last().expect("non-empty")].super_class
        {
            if chain.len() > 64 {
                return Err("class hierarchy too deep".into());
            }
            chain.push(parent);
        }
        let mut object = BbObject {
            class: self.stream.classes[class].name.clone(),
            fields: Vec::new(),
            custom: Vec::new(),
        };
        for desc in chain.into_iter().rev() {
            let desc = &self.stream.classes[desc];
            for field in &desc.fields {
                let value = self.field(field.code)?;
                object.fields.push((field.name.clone(), value));
            }
            if desc.flags & SC_WRITE_METHOD != 0 {
                let name = desc.name.clone();
                object.custom.extend(self.custom(&name)?);
            }
        }
        Ok(BbValue::Object(object))
    }

    fn field(&mut self, code: u8) -> Result<BbValue, String> {
        let mut be = |n: usize| -> Result<Vec<u8>, String> { Ok(self.take(n)?.to_vec()) };
        Ok(match code {
            b'B' => BbValue::Int(i64::from(be(1)?[0] as i8)),
            b'Z' => BbValue::Bool(be(1)?[0] != 0),
            b'C' | b'S' => BbValue::Int(i64::from(i16::from_be_bytes(
                be(2)?.try_into().expect("two"),
            ))),
            b'I' => BbValue::Int(i64::from(i32::from_be_bytes(
                be(4)?.try_into().expect("four"),
            ))),
            b'J' => BbValue::Int(i64::from_be_bytes(be(8)?.try_into().expect("eight"))),
            b'F' => BbValue::Float(f64::from(f32::from_be_bytes(
                be(4)?.try_into().expect("four"),
            ))),
            b'D' => BbValue::Float(f64::from_be_bytes(be(8)?.try_into().expect("eight"))),
            b'L' | b'[' => self.value()?,
            other => return Err(format!("unknown field type {:?}", other as char)),
        })
    }

    fn custom(&mut self, class: &str) -> Result<Vec<BbValue>, String> {
        match class {
            "java.util.TreeMap" => {
                let size = self.length(2)?;
                let mut entries = Vec::new();
                for _ in 0..size {
                    entries.push(self.value()?);
                    entries.push(self.value()?);
                }
                Ok(entries)
            }
            "gnu.trove.TIntArrayList" => {
                let size = self.length(4)?;
                Ok((0..size)
                    .map(|_| self.le_i32().map(|value| BbValue::Int(i64::from(value))))
                    .collect::<Result<_, _>>()?)
            }
            "java.util.HashSet" => {
                self.take(8)?; // capacity, load factor
                let size = self.length(1)?;
                (0..size).map(|_| self.value()).collect()
            }
            other => Err(format!("no reader for custom data of {other}")),
        }
    }
}

fn le_primitive(code: u8, chunk: &[u8]) -> BbValue {
    match code {
        b'B' => BbValue::Int(i64::from(chunk[0] as i8)),
        b'Z' => BbValue::Bool(chunk[0] != 0),
        b'S' | b'C' => BbValue::Int(i64::from(i16::from_le_bytes(
            chunk.try_into().expect("two"),
        ))),
        b'I' => BbValue::Int(i64::from(i32::from_le_bytes(
            chunk.try_into().expect("four"),
        ))),
        b'J' => BbValue::Int(i64::from_le_bytes(chunk.try_into().expect("eight"))),
        b'F' => BbValue::Float(f64::from(f32::from_le_bytes(
            chunk.try_into().expect("four"),
        ))),
        _ => BbValue::Float(f64::from_le_bytes(chunk.try_into().expect("eight"))),
    }
}

#[cfg(test)]
mod tests {
    use super::super::java::{ClassDesc, FieldDesc};
    use super::*;

    fn class(name: &str, flags: u8, fields: &[(u8, &str)]) -> ClassDesc {
        ClassDesc {
            name: name.into(),
            flags,
            fields: fields
                .iter()
                .map(|(code, name)| FieldDesc {
                    code: *code,
                    name: (*name).into(),
                })
                .collect(),
            super_class: None,
        }
    }

    #[test]
    fn objects_mix_big_endian_fields_with_little_endian_lengths() {
        let stream = Stream {
            classes: vec![
                class(
                    "Protein",
                    0x02,
                    &[(b'I', "id"), (b'L', "name"), (b'L', "ids")],
                ),
                class("gnu.trove.TIntArrayList", 0x03, &[]),
            ],
            ..Stream::default()
        };
        let mut bytes = vec![0x73, 0x71, 1, 0, 0, 0];
        bytes.extend_from_slice(&7_i32.to_be_bytes());
        bytes.extend_from_slice(&[0x74, 3, 0]);
        bytes.extend_from_slice(b"ABC");
        bytes.extend_from_slice(&[0x73, 0x71, 2, 0, 0, 0]);
        bytes.extend_from_slice(&2_i32.to_le_bytes());
        bytes.extend_from_slice(&5_i32.to_le_bytes());
        bytes.extend_from_slice(&6_i32.to_le_bytes());
        let mut reader = Reader {
            data: &bytes,
            pos: 0,
            stream: &stream,
            classes: &[0, 1],
            depth: 0,
            budget: MAX_VALUES,
        };
        let BbValue::Object(protein) = reader.value().unwrap() else {
            panic!("not an object")
        };
        assert_eq!(reader.pos, bytes.len());
        assert_eq!(protein.int("id"), Some(7));
        assert_eq!(protein.str("name"), Some("ABC"));
        let Some(BbValue::Object(list)) = protein.field("ids") else {
            panic!("no list")
        };
        assert_eq!(list.custom, vec![BbValue::Int(5), BbValue::Int(6)]);
    }

    #[test]
    fn enum_constants_read_as_their_names() {
        // A 5.3.3 protein record's GO annotations (BBProteinSequence13, record
        // at 80000), cut to one: a HashSet (capacity 16, load factor 0.75 =
        // 00 00 40 3f LE, size 1) of a GOAnnotation whose goboID is
        // 00 00 41 f7 big-endian = 16887 (GO:0016887) and whose evidence type
        // is the enum constant 0x7e, class 3, string "IMP".
        let stream = Stream {
            classes: vec![
                class("java.util.HashSet", 0x03, &[]),
                class(
                    "com.proteomesoftware.scaffold.datamodel.GOAnnotation",
                    0x02,
                    &[(b'I', "goboID"), (b'L', "evidenceType")],
                ),
                class(
                    "com.proteomesoftware.scaffold.datamodel.GOAnnotation$EvidenceType",
                    0x12,
                    &[],
                ),
            ],
            ..Stream::default()
        };
        let bytes = [
            0x73, 0x71, 1, 0, 0, 0, 0x10, 0, 0, 0, 0, 0, 0x40, 0x3f, 1, 0, 0, 0, 0x73, 0x71, 2, 0,
            0, 0, 0, 0, 0x41, 0xf7, 0x7e, 0x71, 3, 0, 0, 0, 0x74, 3, 0, b'I', b'M', b'P',
        ];
        let mut reader = Reader {
            data: &bytes,
            pos: 0,
            stream: &stream,
            classes: &[0, 1, 2],
            depth: 0,
            budget: MAX_VALUES,
        };
        let BbValue::Object(set) = reader.value().unwrap() else {
            panic!("not an object")
        };
        assert_eq!(reader.pos, bytes.len());
        let [BbValue::Object(annotation)] = set.custom.as_slice() else {
            panic!("{:?}", set.custom)
        };
        assert_eq!(annotation.int("goboID"), Some(16_887));
        assert_eq!(annotation.str("evidenceType"), Some("IMP"));
    }

    #[test]
    fn unknown_custom_data_is_an_error_not_a_guess() {
        let stream = Stream {
            classes: vec![class("java.util.Vector", 0x03, &[])],
            ..Stream::default()
        };
        let bytes = [0x73, 0x71, 1, 0, 0, 0, 0, 0, 0, 0];
        let mut reader = Reader {
            data: &bytes,
            pos: 0,
            stream: &stream,
            classes: &[0],
            depth: 0,
            budget: MAX_VALUES,
        };
        assert!(reader.value().unwrap_err().contains("java.util.Vector"));
    }
}
