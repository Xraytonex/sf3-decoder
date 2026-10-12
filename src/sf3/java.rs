//! A reader for Java Object Serialization streams (protocol 2), from the
//! public specification. Scaffold stores its project graph, descriptors,
//! and lookup tables this way. Nothing is executed: the stream is parsed
//! into plain values held in arenas and referenced by index.

use std::sync::Arc;

const MAGIC: u16 = 0xACED;
const VERSION: u16 = 5;
const BASE_HANDLE: u32 = 0x7E_0000;
/// Nesting deeper than this is refused rather than risking the stack.
const MAX_DEPTH: usize = 1_000;
/// Values one stream may decode to: about 3× the largest real project seen
/// (a 5.3.3 Experiment of 128 MB holding 856,290 score objects decodes to
/// ~27 million; decoding it peaks at 2.2 GB), far below what exhausts
/// memory. Primitive array elements count one each.
const MAX_VALUES: usize = 80_000_000;

const TC_NULL: u8 = 0x70;
const TC_REFERENCE: u8 = 0x71;
const TC_CLASSDESC: u8 = 0x72;
const TC_OBJECT: u8 = 0x73;
const TC_STRING: u8 = 0x74;
const TC_ARRAY: u8 = 0x75;
const TC_CLASS: u8 = 0x76;
const TC_BLOCKDATA: u8 = 0x77;
const TC_ENDBLOCKDATA: u8 = 0x78;
const TC_RESET: u8 = 0x79;
const TC_BLOCKDATALONG: u8 = 0x7A;
const TC_LONGSTRING: u8 = 0x7C;
const TC_PROXYCLASSDESC: u8 = 0x7D;
const TC_ENUM: u8 = 0x7E;

pub const SC_WRITE_METHOD: u8 = 0x01;
const SC_EXTERNALIZABLE: u8 = 0x04;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("byte {offset}: {problem}")]
pub struct JavaError {
    pub offset: usize,
    pub problem: String,
}

/// A field or element value. Objects, arrays, strings, class descriptors,
/// and block data live in the [`Stream`]'s arenas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    Char(u16),
    Str(usize),
    Object(usize),
    Array(usize),
    ClassDesc(usize),
    Class(usize),
    Enum(usize, usize),
    Block(usize),
}

#[derive(Debug, Clone)]
pub struct FieldDesc {
    pub code: u8,
    /// Shared with every object of the class (a 5.3.3 project held ~27
    /// million values; a copy per object doubled the memory).
    pub name: Arc<str>,
}

#[derive(Debug, Clone)]
pub struct ClassDesc {
    pub name: String,
    pub flags: u8,
    pub fields: Vec<FieldDesc>,
    pub super_class: Option<usize>,
}

#[derive(Debug, Clone, Default)]
pub struct Object {
    pub class: usize,
    pub fields: Vec<(Arc<str>, Value)>,
    /// Data written by `writeObject` / `writeExternal`, per class in the
    /// hierarchy that wrote any.
    pub annotations: Vec<(usize, Vec<Value>)>,
}

/// Every element type is decoded (to stay in step with the stream), though
/// the decoder only reads some.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum ArrayData {
    Bytes(Vec<u8>),
    Ints(Vec<i32>),
    Longs(Vec<i64>),
    Floats(Vec<f32>),
    Doubles(Vec<f64>),
    Shorts(Vec<i16>),
    Chars(String),
    Values(Vec<Value>),
}

#[derive(Debug, Clone)]
pub struct Array {
    pub data: ArrayData,
}

/// A parsed stream: its top-level contents plus the arenas they point into.
#[derive(Debug, Default)]
pub struct Stream {
    pub roots: Vec<Value>,
    pub classes: Vec<ClassDesc>,
    pub objects: Vec<Object>,
    pub arrays: Vec<Array>,
    pub strings: Vec<String>,
    pub blocks: Vec<Vec<u8>>,
}

impl Stream {
    pub fn parse(data: &[u8]) -> Result<Self, JavaError> {
        let mut reader = Reader {
            data,
            pos: 0,
            depth: 0,
            budget: MAX_VALUES,
            handles: Vec::new(),
            stream: Stream::default(),
        };
        let magic = reader.u16()?;
        let version = reader.u16()?;
        if magic != MAGIC || version != VERSION {
            return Err(reader.error("not a Java serialization stream (protocol 2)"));
        }
        while reader.pos < data.len() {
            let value = reader.content()?;
            reader.stream.roots.push(value);
        }
        Ok(reader.stream)
    }

    pub fn class_name(&self, object: usize) -> &str {
        &self.classes[self.objects[object].class].name
    }

    /// Whether the object's class name ends with `.{simple}` (or is `simple`).
    pub fn is(&self, object: usize, simple: &str) -> bool {
        let name = self.class_name(object);
        name == simple
            || name
                .strip_suffix(simple)
                .is_some_and(|rest| rest.ends_with('.') || rest.ends_with('$'))
    }

    pub fn field(&self, object: usize, name: &str) -> Option<Value> {
        self.objects[object]
            .fields
            .iter()
            .find(|(field, _)| field.as_ref() == name)
            .map(|(_, value)| *value)
    }

    pub fn int(&self, object: usize, name: &str) -> Option<i64> {
        match self.field(object, name)? {
            Value::Int(value) => Some(value),
            _ => None,
        }
    }

    pub fn float(&self, object: usize, name: &str) -> Option<f64> {
        match self.field(object, name)? {
            Value::Float(value) => Some(value),
            _ => None,
        }
    }

    pub fn str(&self, object: usize, name: &str) -> Option<&str> {
        match self.field(object, name)? {
            Value::Str(index) => Some(&self.strings[index]),
            _ => None,
        }
    }

    pub fn object(&self, object: usize, name: &str) -> Option<usize> {
        match self.field(object, name)? {
            Value::Object(index) => Some(index),
            _ => None,
        }
    }

    pub fn array(&self, object: usize, name: &str) -> Option<&ArrayData> {
        match self.field(object, name)? {
            Value::Array(index) => Some(&self.arrays[index].data),
            _ => None,
        }
    }

    /// Every object reachable from the roots, each once, in a stable order.
    pub fn reachable_objects(&self) -> Vec<usize> {
        let mut seen = vec![false; self.objects.len()];
        let mut seen_arrays = vec![false; self.arrays.len()];
        let mut order = Vec::new();
        let mut stack: Vec<Value> = self.roots.iter().rev().copied().collect();
        while let Some(value) = stack.pop() {
            match value {
                Value::Object(index) if !seen[index] => {
                    seen[index] = true;
                    order.push(index);
                    let object = &self.objects[index];
                    for (_, values) in object.annotations.iter().rev() {
                        stack.extend(values.iter().rev());
                    }
                    stack.extend(object.fields.iter().rev().map(|(_, value)| *value));
                }
                Value::Array(index) if !seen_arrays[index] => {
                    seen_arrays[index] = true;
                    if let ArrayData::Values(values) = &self.arrays[index].data {
                        stack.extend(values.iter().rev());
                    }
                }
                _ => {}
            }
        }
        order
    }

    /// Reachable objects whose class name ends with `simple`.
    pub fn objects_of(&self, simple: &str) -> Vec<usize> {
        self.reachable_objects()
            .into_iter()
            .filter(|&object| self.is(object, simple))
            .collect()
    }

    /// The `int` values a `gnu.trove.TIntArrayList` wrote: its block data
    /// as big-endian ints, a count first.
    pub fn trove_ints(&self, list: usize) -> Option<Vec<i32>> {
        let mut raw = Vec::new();
        for (_, values) in &self.objects[list].annotations {
            for value in values {
                if let Value::Block(index) = value {
                    raw.extend_from_slice(&self.blocks[*index]);
                }
            }
        }
        let ints: Vec<i32> = raw
            .as_chunks::<4>()
            .0
            .iter()
            .map(|chunk| i32::from_be_bytes(*chunk))
            .collect();
        let count = usize::try_from(*ints.first()?).ok()?;
        ints.get(1..=count).map(<[i32]>::to_vec)
    }
}

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
    depth: usize,
    /// Values still allowed for this stream.
    budget: usize,
    handles: Vec<Value>,
    stream: Stream,
}

impl Reader<'_> {
    fn error(&self, problem: impl Into<String>) -> JavaError {
        JavaError {
            offset: self.pos,
            problem: problem.into(),
        }
    }

    fn take(&mut self, n: usize) -> Result<&[u8], JavaError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.data.len())
            .ok_or_else(|| self.error(format!("stream ends early (wanted {n} bytes)")))?;
        let bytes = &self.data[self.pos..end];
        self.pos = end;
        Ok(bytes)
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N], JavaError> {
        let bytes = self.take(N)?;
        Ok(bytes.try_into().expect("length checked"))
    }

    fn u8(&mut self) -> Result<u8, JavaError> {
        Ok(self.array::<1>()?[0])
    }

    fn u16(&mut self) -> Result<u16, JavaError> {
        Ok(u16::from_be_bytes(self.array()?))
    }

    fn i32(&mut self) -> Result<i32, JavaError> {
        Ok(i32::from_be_bytes(self.array()?))
    }

    fn i64(&mut self) -> Result<i64, JavaError> {
        Ok(i64::from_be_bytes(self.array()?))
    }

    fn peek(&self) -> Result<u8, JavaError> {
        self.data
            .get(self.pos)
            .copied()
            .ok_or_else(|| self.error("stream ends early"))
    }

    fn utf(&mut self) -> Result<String, JavaError> {
        let length = usize::from(self.u16()?);
        Ok(String::from_utf8_lossy(self.take(length)?).into_owned())
    }

    fn length(&mut self, element_size: usize) -> Result<usize, JavaError> {
        let length = self.i32()?;
        let length = usize::try_from(length).map_err(|_| self.error("negative length"))?;
        if length.saturating_mul(element_size) > self.data.len() - self.pos {
            return Err(self.error(format!("length {length} runs past the stream")));
        }
        self.spend(length)?;
        Ok(length)
    }

    fn spend(&mut self, values: usize) -> Result<(), JavaError> {
        self.budget = self
            .budget
            .checked_sub(values)
            .ok_or_else(|| self.error("more values than any real project holds"))?;
        Ok(())
    }

    fn handle(&self, offset: usize) -> Result<Value, JavaError> {
        let raw = u32::try_from(self.i32_at(offset)?).unwrap_or(0);
        raw.checked_sub(BASE_HANDLE)
            .and_then(|index| self.handles.get(index as usize).copied())
            .ok_or_else(|| self.error(format!("unknown handle {raw:#x}")))
    }

    fn i32_at(&self, offset: usize) -> Result<i32, JavaError> {
        self.data
            .get(offset..offset + 4)
            .map(|bytes| i32::from_be_bytes(bytes.try_into().expect("four bytes")))
            .ok_or_else(|| self.error("stream ends early"))
    }

    fn reference(&mut self) -> Result<Value, JavaError> {
        let offset = self.pos;
        self.take(4)?;
        self.handle(offset)
    }

    fn enter(&mut self) -> Result<(), JavaError> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return Err(self.error("nesting too deep"));
        }
        Ok(())
    }

    fn content(&mut self) -> Result<Value, JavaError> {
        match self.peek()? {
            TC_BLOCKDATA => {
                self.pos += 1;
                let length = usize::from(self.u8()?);
                self.block(length)
            }
            TC_BLOCKDATALONG => {
                self.pos += 1;
                let length = self.length(1)?;
                self.block(length)
            }
            _ => self.value(),
        }
    }

    fn block(&mut self, length: usize) -> Result<Value, JavaError> {
        let bytes = self.take(length)?.to_vec();
        self.stream.blocks.push(bytes);
        Ok(Value::Block(self.stream.blocks.len() - 1))
    }

    fn value(&mut self) -> Result<Value, JavaError> {
        self.spend(1)?;
        self.enter()?;
        let tc = self.u8()?;
        let value = match tc {
            TC_NULL => Value::Null,
            TC_REFERENCE => self.reference()?,
            TC_STRING | TC_LONGSTRING => {
                let text = if tc == TC_STRING {
                    self.utf()?
                } else {
                    let length = usize::try_from(self.i64()?)
                        .map_err(|_| self.error("negative string length"))?;
                    String::from_utf8_lossy(self.take(length)?).into_owned()
                };
                self.stream.strings.push(text);
                let value = Value::Str(self.stream.strings.len() - 1);
                self.handles.push(value);
                value
            }
            TC_CLASSDESC | TC_PROXYCLASSDESC => {
                self.pos -= 1;
                self.class_desc()?.map_or(Value::Null, Value::ClassDesc)
            }
            TC_CLASS => {
                let class = self.class_desc()?.ok_or_else(|| self.error("null class"))?;
                let value = Value::Class(class);
                self.handles.push(value);
                value
            }
            TC_OBJECT => self.new_object()?,
            TC_ARRAY => self.new_array()?,
            TC_ENUM => {
                let class = self
                    .class_desc()?
                    .ok_or_else(|| self.error("null enum class"))?;
                let slot = self.handles.len();
                self.handles.push(Value::Null);
                let name = match self.value()? {
                    Value::Str(index) => index,
                    _ => return Err(self.error("enum constant without a name")),
                };
                let value = Value::Enum(class, name);
                self.handles[slot] = value;
                value
            }
            TC_RESET => {
                self.handles.clear();
                self.value()?
            }
            other => return Err(self.error(format!("unsupported type code {other:#x}"))),
        };
        self.depth -= 1;
        Ok(value)
    }

    fn class_desc(&mut self) -> Result<Option<usize>, JavaError> {
        self.enter()?;
        let result = match self.u8()? {
            TC_NULL => None,
            TC_REFERENCE => match self.reference()? {
                Value::ClassDesc(index) => Some(index),
                _ => return Err(self.error("reference is not a class descriptor")),
            },
            TC_PROXYCLASSDESC => {
                let count = self.length(2)?;
                let names = (0..count)
                    .map(|_| self.utf())
                    .collect::<Result<Vec<_>, _>>()?;
                let index = self.push_class(ClassDesc {
                    name: format!("Proxy[{}]", names.join(",")),
                    flags: 0,
                    fields: Vec::new(),
                    super_class: None,
                });
                self.annotation()?;
                self.stream.classes[index].super_class = self.class_desc()?;
                Some(index)
            }
            TC_CLASSDESC => {
                let name = self.utf()?;
                self.i64()?; // serialVersionUID
                let index = self.push_class(ClassDesc {
                    name,
                    flags: 0,
                    fields: Vec::new(),
                    super_class: None,
                });
                let flags = self.u8()?;
                let count = usize::from(self.u16()?);
                let mut fields = Vec::with_capacity(count);
                for _ in 0..count {
                    let code = self.u8()?;
                    let name = self.utf()?;
                    if matches!(code, b'[' | b'L') {
                        self.value()?; // the field's type name
                    }
                    fields.push(FieldDesc {
                        code,
                        name: name.into(),
                    });
                }
                self.stream.classes[index].flags = flags;
                self.stream.classes[index].fields = fields;
                self.annotation()?;
                self.stream.classes[index].super_class = self.class_desc()?;
                Some(index)
            }
            other => return Err(self.error(format!("expected a class descriptor, got {other:#x}"))),
        };
        self.depth -= 1;
        Ok(result)
    }

    fn push_class(&mut self, class: ClassDesc) -> usize {
        self.stream.classes.push(class);
        let index = self.stream.classes.len() - 1;
        self.handles.push(Value::ClassDesc(index));
        index
    }

    fn annotation(&mut self) -> Result<Vec<Value>, JavaError> {
        let mut values = Vec::new();
        loop {
            if self.peek()? == TC_ENDBLOCKDATA {
                self.pos += 1;
                return Ok(values);
            }
            values.push(self.content()?);
        }
    }

    /// Class descriptors from the root of the hierarchy down.
    fn chain(&self, class: usize) -> Result<Vec<usize>, JavaError> {
        let mut chain = vec![class];
        let mut current = class;
        while let Some(parent) = self.stream.classes[current].super_class {
            if chain.len() > 64 || chain.contains(&parent) {
                return Err(self.error("class hierarchy loops"));
            }
            chain.push(parent);
            current = parent;
        }
        chain.reverse();
        Ok(chain)
    }

    fn new_object(&mut self) -> Result<Value, JavaError> {
        let class = self
            .class_desc()?
            .ok_or_else(|| self.error("object without a class"))?;
        self.stream.objects.push(Object {
            class,
            ..Object::default()
        });
        let index = self.stream.objects.len() - 1;
        self.handles.push(Value::Object(index));
        for desc in self.chain(class)? {
            let flags = self.stream.classes[desc].flags;
            if flags & SC_EXTERNALIZABLE != 0 {
                let data = self.annotation()?;
                self.stream.objects[index].annotations.push((desc, data));
                continue;
            }
            let fields = self.stream.classes[desc].fields.clone();
            for field in fields {
                let value = self.primitive_or_value(field.code)?;
                self.stream.objects[index].fields.push((field.name, value));
            }
            if flags & SC_WRITE_METHOD != 0 {
                let data = self.annotation()?;
                self.stream.objects[index].annotations.push((desc, data));
            }
        }
        Ok(Value::Object(index))
    }

    fn primitive_or_value(&mut self, code: u8) -> Result<Value, JavaError> {
        Ok(match code {
            b'B' => Value::Int(i64::from(self.u8()? as i8)),
            b'C' => Value::Char(self.u16()?),
            b'D' => Value::Float(f64::from_be_bytes(self.array()?)),
            b'F' => Value::Float(f64::from(f32::from_be_bytes(self.array()?))),
            b'I' => Value::Int(i64::from(self.i32()?)),
            b'J' => Value::Int(self.i64()?),
            b'S' => Value::Int(i64::from(i16::from_be_bytes(self.array()?))),
            b'Z' => Value::Bool(self.u8()? != 0),
            b'L' | b'[' => self.value()?,
            other => return Err(self.error(format!("unknown field type {:?}", other as char))),
        })
    }

    fn new_array(&mut self) -> Result<Value, JavaError> {
        let class = self
            .class_desc()?
            .ok_or_else(|| self.error("array without a class"))?;
        let element = self.stream.classes[class]
            .name
            .as_bytes()
            .get(1)
            .copied()
            .ok_or_else(|| self.error("array class without an element type"))?;
        self.stream.arrays.push(Array {
            data: ArrayData::Values(Vec::new()),
        });
        let index = self.stream.arrays.len() - 1;
        self.handles.push(Value::Array(index));
        let data = match element {
            b'B' | b'Z' => {
                let length = self.length(1)?;
                ArrayData::Bytes(self.take(length)?.to_vec())
            }
            b'I' => {
                let length = self.length(4)?;
                ArrayData::Ints(
                    self.take(length * 4)?
                        .as_chunks::<4>()
                        .0
                        .iter()
                        .map(|c| i32::from_be_bytes(*c))
                        .collect(),
                )
            }
            b'J' => {
                let length = self.length(8)?;
                ArrayData::Longs(
                    self.take(length * 8)?
                        .as_chunks::<8>()
                        .0
                        .iter()
                        .map(|c| i64::from_be_bytes(*c))
                        .collect(),
                )
            }
            b'F' => {
                let length = self.length(4)?;
                ArrayData::Floats(
                    self.take(length * 4)?
                        .as_chunks::<4>()
                        .0
                        .iter()
                        .map(|c| f32::from_be_bytes(*c))
                        .collect(),
                )
            }
            b'D' => {
                let length = self.length(8)?;
                ArrayData::Doubles(
                    self.take(length * 8)?
                        .as_chunks::<8>()
                        .0
                        .iter()
                        .map(|c| f64::from_be_bytes(*c))
                        .collect(),
                )
            }
            b'S' => {
                let length = self.length(2)?;
                ArrayData::Shorts(
                    self.take(length * 2)?
                        .as_chunks::<2>()
                        .0
                        .iter()
                        .map(|c| i16::from_be_bytes(*c))
                        .collect(),
                )
            }
            b'C' => {
                let length = self.length(2)?;
                let units: Vec<u16> = self
                    .take(length * 2)?
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .map(|c| u16::from_be_bytes(*c))
                    .collect();
                ArrayData::Chars(String::from_utf16_lossy(&units))
            }
            b'L' | b'[' => {
                let length = self.length(1)?;
                let mut values = Vec::with_capacity(length);
                for _ in 0..length {
                    values.push(self.value()?);
                }
                ArrayData::Values(values)
            }
            other => {
                return Err(self.error(format!("unknown array element type {:?}", other as char)));
            }
        };
        self.stream.arrays[index].data = data;
        Ok(Value::Array(index))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `new Point(3, 4)` for `class Point implements Serializable { int x; int y; }`,
    /// followed by a back-reference to it and a string, hand-assembled from
    /// the specification.
    fn point_stream() -> Vec<u8> {
        let mut bytes = vec![0xAC, 0xED, 0x00, 0x05, TC_OBJECT, TC_CLASSDESC];
        bytes.extend_from_slice(&[0x00, 0x05]);
        bytes.extend_from_slice(b"Point");
        bytes.extend_from_slice(&7_i64.to_be_bytes()); // serialVersionUID
        bytes.push(0x02); // SC_SERIALIZABLE
        bytes.extend_from_slice(&[0x00, 0x02]);
        bytes.extend_from_slice(&[b'I', 0x00, 0x01, b'x']);
        bytes.extend_from_slice(&[b'I', 0x00, 0x01, b'y']);
        bytes.push(TC_ENDBLOCKDATA);
        bytes.push(TC_NULL); // no superclass
        bytes.extend_from_slice(&3_i32.to_be_bytes());
        bytes.extend_from_slice(&4_i32.to_be_bytes());
        // Handle 0 is the class descriptor, handle 1 the object.
        bytes.push(TC_REFERENCE);
        bytes.extend_from_slice(&(BASE_HANDLE + 1).to_be_bytes());
        bytes.extend_from_slice(&[TC_STRING, 0x00, 0x02, b'h', b'i']);
        bytes
    }

    #[test]
    fn objects_fields_references_and_strings_parse() {
        let stream = Stream::parse(&point_stream()).unwrap();
        assert_eq!(stream.roots.len(), 3);
        let Value::Object(point) = stream.roots[0] else {
            panic!("{:?}", stream.roots[0])
        };
        assert!(stream.is(point, "Point"));
        assert_eq!(stream.int(point, "x"), Some(3));
        assert_eq!(stream.int(point, "y"), Some(4));
        assert_eq!(stream.roots[1], Value::Object(point));
        let Value::Str(text) = stream.roots[2] else {
            panic!("{:?}", stream.roots[2])
        };
        assert_eq!(stream.strings[text], "hi");
        assert_eq!(stream.objects_of("Point"), vec![point]);
    }

    #[test]
    fn truncated_and_foreign_input_is_refused() {
        let bytes = point_stream();
        let error = Stream::parse(&bytes[..bytes.len() - 3]).unwrap_err();
        assert!(error.problem.contains("ends early"), "{error}");
        assert!(Stream::parse(b"PK\x03\x04").is_err());
    }

    #[test]
    fn a_huge_declared_length_is_refused_before_allocating() {
        let mut bytes = vec![0xAC, 0xED, 0x00, 0x05, TC_ARRAY, TC_CLASSDESC, 0x00, 0x02];
        bytes.extend_from_slice(b"[I");
        bytes.extend_from_slice(&0_i64.to_be_bytes());
        bytes.extend_from_slice(&[0x02, 0x00, 0x00, TC_ENDBLOCKDATA, TC_NULL]);
        bytes.extend_from_slice(&i32::MAX.to_be_bytes());
        let error = Stream::parse(&bytes).unwrap_err();
        assert!(error.problem.contains("runs past"), "{error}");
    }
}
