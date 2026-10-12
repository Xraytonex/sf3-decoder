//! The `.sf3` container (a SQLite file of gzip-compressed inner files) and
//! its fixed-size record files.
//!
//! A record type `Name` has a descriptor stream `Name` holding a
//! `MathmaticTranslationTable`, and parts `Name0`, `Name1`, ... of
//! fixed-size blocks. Record id → part `id >> FILE_SHIFT`, block
//! `id & (2^FILE_SHIFT − 1)`; the table lists the last block used in each
//! part, and every block starts with its own id (checked on every read).

use std::io::Read;
use std::path::Path;

use std::cell::Cell;

use flate2::read::MultiGzDecoder;
use rusqlite::{Connection, OpenFlags, OptionalExtension};

use super::java::Stream;
use super::{Sf3Error, layout};

/// No inner file of a real project comes near this; a larger one is
/// refused rather than inflated.
const MAX_INNER_FILE: u64 = 2 << 30;
/// Everything inflated from one project, so a small file can't expand into
/// more memory than a real project ever needs.
const MAX_INFLATED: u64 = 8 << 30;
/// `FILE_SHIFT` values that keep every record id inside an `i32`.
const SHIFTS: std::ops::RangeInclusive<u32> = 1..=30;
pub(crate) const STRING_BLOCK: usize = 68;
/// A string block: its chain's id, 60 bytes of text, the next block.
const STRING_TEXT: std::ops::Range<usize> = 4..64;
/// `PeptideIdentification` and other records pad text with these.
const STRING_PADDING: u8 = 0x02;
const INT_LIST_END: i32 = i32::MIN;

pub(crate) struct Container {
    conn: Connection,
    inflated: Cell<u64>,
}

impl Container {
    pub fn open(path: &Path) -> Result<Self, Sf3Error> {
        let conn = Connection::open_with_flags(
            path,
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        )
        .map_err(|error| Sf3Error::NotSf3(format!("cannot open as SQLite: {error}")))?;
        let has_files: bool = conn
            .query_row(
                "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name = 'Files'",
                [],
                |row| row.get::<_, i64>(0),
            )
            .map(|count| count > 0)
            .map_err(|error| Sf3Error::NotSf3(format!("cannot read as SQLite: {error}")))?;
        if !has_files {
            return Err(Sf3Error::NotSf3("no Files table".into()));
        }
        Ok(Self {
            conn,
            inflated: Cell::new(0),
        })
    }

    /// An inner file, decompressed; `None` when the project has no such file.
    pub fn file(&self, name: &str) -> Result<Option<Vec<u8>>, Sf3Error> {
        let blob: Option<Option<Vec<u8>>> = self
            .conn
            .query_row(
                "SELECT file FROM Files WHERE filename = ?1",
                [name],
                |row| row.get(0),
            )
            .optional()
            .map_err(|error| Sf3Error::Container(format!("{name}: {error}")))?;
        let Some(blob) = blob.flatten() else {
            return Ok(None);
        };
        if !blob.starts_with(&[0x1f, 0x8b]) {
            return Ok(Some(blob));
        }
        let budget = MAX_INNER_FILE.min(MAX_INFLATED.saturating_sub(self.inflated.get()));
        let mut out = Vec::new();
        MultiGzDecoder::new(blob.as_slice())
            .take(budget + 1)
            .read_to_end(&mut out)
            .map_err(|error| Sf3Error::Container(format!("{name}: {error}")))?;
        let size = out.len() as u64;
        if size > budget {
            return Err(Sf3Error::Container(format!(
                "{name}: inflates beyond what a real project holds"
            )));
        }
        self.inflated.set(self.inflated.get() + size);
        Ok(Some(out))
    }

    /// The Scaffold version that wrote the file (`Version.version_name`).
    pub fn version(&self) -> Option<String> {
        self.conn
            .query_row("SELECT version_name FROM Version", [], |row| row.get(0))
            .ok()
    }

    /// The `ProteinGroup` table: (protein group, MS sample, probability).
    /// Empty when the project has none.
    pub fn protein_group_rows(&self) -> Result<Vec<(i64, i64, f64)>, Sf3Error> {
        let table = |error: rusqlite::Error| Sf3Error::Container(format!("ProteinGroup: {error}"));
        let exists: i64 = self
            .conn
            .query_row(
                "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name = 'ProteinGroup'",
                [],
                |row| row.get(0),
            )
            .map_err(table)?;
        if exists == 0 {
            return Ok(Vec::new());
        }
        let mut statement = self
            .conn
            .prepare("SELECT peg, msSampleId, matchProbability FROM ProteinGroup ORDER BY id")
            .map_err(table)?;
        let rows = statement
            .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))
            .map_err(table)?;
        rows.collect::<Result<_, _>>().map_err(table)
    }

    pub fn stream(&self, name: &str) -> Result<Stream, Sf3Error> {
        let data = self
            .file(name)?
            .ok_or_else(|| Sf3Error::Missing(name.to_string()))?;
        Stream::parse(&data).map_err(|error| Sf3Error::Java {
            file: name.to_string(),
            error,
        })
    }

    /// Like [`Container::stream`], but a missing file is `None`.
    pub fn optional_stream(&self, name: &str) -> Result<Option<Stream>, Sf3Error> {
        match self.stream(name) {
            Ok(stream) => Ok(Some(stream)),
            Err(Sf3Error::Missing(_)) => Ok(None),
            Err(error) => Err(error),
        }
    }

    /// Fixed-size records; every block in use must carry its own id.
    pub fn records(&self, name: &str, size: usize) -> Result<Records, Sf3Error> {
        let records = self.load(name, size, false)?;
        records.check_ids()?;
        Ok(records)
    }

    /// String records (68-byte blocks chained by a next pointer): a block
    /// carries the id of its chain, which starts at the block with that id.
    pub fn strings(&self, name: &str) -> Result<Records, Sf3Error> {
        self.chains(name, STRING_BLOCK)
    }

    /// Chained records: fixed-size blocks, each carrying its chain's id
    /// first and the next block's index in its part last (−1 ends the
    /// chain). Strings and peak lists are stored this way.
    pub fn chains(&self, name: &str, size: usize) -> Result<Records, Sf3Error> {
        let records = self.load(name, size, true)?;
        records.check_chain_ids()?;
        Ok(records)
    }

    /// Chained parts may run past `2^FILE_SHIFT` blocks: a chain continues
    /// in its own part, wherever there is room.
    fn load(&self, name: &str, size: usize, chained: bool) -> Result<Records, Sf3Error> {
        let descriptor = self.stream(name)?;
        let (shift, last) = translation_table(&descriptor)
            .ok_or_else(|| layout(name, "no translation table in its descriptor"))?;
        let ids_fit = (last.len() as u64) << shift <= i32::MAX as u64 + 1;
        if !SHIFTS.contains(&shift) || !ids_fit {
            return Err(layout(
                name,
                &format!(
                    "{} parts of 2^{shift} records don't fit record ids",
                    last.len()
                ),
            ));
        }
        let mut parts = Vec::with_capacity(last.len());
        for (part, &last_block) in last.iter().enumerate() {
            let data = self.file(&format!("{name}{part}"))?.unwrap_or_default();
            if data.len() % size != 0 {
                return Err(layout(
                    name,
                    &format!(
                        "part {part} is {} bytes, not a multiple of {size}",
                        data.len()
                    ),
                ));
            }
            let addressable = chained || i64::from(last_block) < 1_i64 << shift;
            if last_block < -1 || !addressable {
                return Err(layout(
                    name,
                    &format!("part {part} lists block {last_block}"),
                ));
            }
            let used = usize::try_from(last_block).map_or(0, |last| last + 1);
            if used * size > data.len() {
                return Err(layout(
                    name,
                    &format!(
                        "part {part} lists {used} blocks but holds {}",
                        data.len() / size
                    ),
                ));
            }
            parts.push((data, used));
        }
        let count = parts.iter().map(|(_, used)| *used).sum();
        Ok(Records {
            name: name.to_string(),
            size,
            shift,
            parts,
            count,
        })
    }
}

/// `FILE_SHIFT` and the last used block of each part.
fn translation_table(stream: &Stream) -> Option<(u32, Vec<i32>)> {
    let table = *stream.objects_of("MathmaticTranslationTable").first()?;
    let shift = u32::try_from(stream.int(table, "FILE_SHIFT")?).ok()?;
    let files = stream.object(table, "files")?;
    Some((shift, stream.trove_ints(files)?))
}

pub(crate) struct Records {
    name: String,
    size: usize,
    shift: u32,
    /// Each part's bytes and how many of its blocks are in use.
    parts: Vec<(Vec<u8>, usize)>,
    /// Blocks in use, which bounds any chain of blocks.
    count: usize,
}

impl Records {
    /// Every block in use must carry its own id; id 0 is never used.
    fn check_ids(&self) -> Result<(), Sf3Error> {
        for id in self.ids() {
            if self.block(id).and_then(|block| read_i32(block, 0)) != Some(id) {
                return Err(layout(
                    &self.name,
                    &format!("block of record {id} carries another id"),
                ));
            }
        }
        Ok(())
    }

    /// A string block carries the id of the chain it belongs to, which
    /// starts in the same part at or before it.
    fn check_chain_ids(&self) -> Result<(), Sf3Error> {
        let mask = (1_usize << self.shift) - 1;
        for (part, (data, used)) in self.parts.iter().enumerate() {
            for block in usize::from(part == 0)..*used {
                let id =
                    read_i32(&data[block * self.size..], 0).and_then(|id| usize::try_from(id).ok());
                let fits =
                    id.is_some_and(|id| id > 0 && id >> self.shift == part && id & mask <= block);
                if !fits {
                    return Err(layout(
                        &self.name,
                        &format!("string block {block} of part {part} has a bad id"),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Block `block` of `part`, if in use.
    fn block_at(&self, part: usize, block: usize) -> Option<&[u8]> {
        let (data, used) = self.parts.get(part)?;
        (block < *used).then(|| &data[block * self.size..(block + 1) * self.size])
    }

    /// Ids of every record in use, ascending.
    pub fn ids(&self) -> impl Iterator<Item = i32> + '_ {
        self.parts
            .iter()
            .enumerate()
            .flat_map(move |(part, (_, used))| {
                let first = usize::from(part == 0);
                (first..*used)
                    .filter_map(move |block| i32::try_from((part << self.shift) | block).ok())
            })
    }

    fn block(&self, id: i32) -> Option<&[u8]> {
        let id = usize::try_from(id).ok().filter(|&id| id > 0)?;
        let (data, used) = self.parts.get(id >> self.shift)?;
        let block = id & ((1 << self.shift) - 1);
        (block < *used).then(|| &data[block * self.size..(block + 1) * self.size])
    }

    /// The record with this id, if stored.
    pub fn get(&self, id: i32) -> Option<&[u8]> {
        self.block(id)
    }

    /// The blocks of a chained record, in order: each block's last four
    /// bytes give the next block's index in the same part, −1 ends it.
    /// `None` when a block is missing, belongs to another chain, or the
    /// chain never ends.
    pub fn chain(&self, id: i32) -> Option<Vec<&[u8]>> {
        let start = usize::try_from(id).ok().filter(|&id| id > 0)?;
        let part = start >> self.shift;
        let mut next = start & ((1 << self.shift) - 1);
        let mut blocks = Vec::new();
        for _ in 0..=self.count {
            // Chains continue within their part: `next` is a block index there.
            let block = self.block_at(part, next)?;
            if read_i32(block, 0) != Some(id) {
                return None; // the chain wandered into another record
            }
            blocks.push(block);
            match read_i32(block, self.size - 4)? {
                -1 => return Some(blocks),
                following => next = usize::try_from(following).ok()?,
            }
        }
        None // a chain that never ends
    }

    /// A string record: the text of a chain of blocks.
    pub fn string(&self, id: i32) -> Option<String> {
        let mut bytes = Vec::new();
        for block in self.chain(id)? {
            let text = &block[STRING_TEXT];
            let end = text
                .iter()
                .position(|&byte| byte == STRING_PADDING)
                .unwrap_or(text.len());
            bytes.extend_from_slice(&text[..end]);
        }
        Some(String::from_utf8_lossy(&bytes).into_owned())
    }

    /// An int-list record: values after the id, up to `INT_MIN`.
    /// `None` when the list has no terminator in its block: a longer list
    /// would continue somewhere unknown.
    pub fn int_list(&self, id: i32) -> Option<Vec<i32>> {
        let block = self.block(id)?;
        let values: Vec<i32> = (4..self.size)
            .step_by(4)
            .filter_map(|offset| read_i32(block, offset))
            .collect();
        let end = values.iter().position(|&value| value == INT_LIST_END)?;
        Some(values[..end].to_vec())
    }
}

pub(crate) fn read_i32(block: &[u8], offset: usize) -> Option<i32> {
    block
        .get(offset..offset + 4)
        .map(|bytes| i32::from_le_bytes(bytes.try_into().expect("four bytes")))
}

pub(crate) fn read_f32(block: &[u8], offset: usize) -> Option<f32> {
    block
        .get(offset..offset + 4)
        .map(|bytes| f32::from_le_bytes(bytes.try_into().expect("four bytes")))
}

pub(crate) fn read_f64(block: &[u8], offset: usize) -> Option<f64> {
    block
        .get(offset..offset + 8)
        .map(|bytes| f64::from_le_bytes(bytes.try_into().expect("eight bytes")))
}

/// Records over one part of `blocks` (shift 15), for tests.
#[cfg(test)]
pub(crate) fn test_records(size: usize, blocks: &[Vec<u8>]) -> Records {
    let data: Vec<u8> = blocks.iter().flatten().copied().collect();
    Records {
        name: "Test".into(),
        size,
        shift: 15,
        parts: vec![(data, blocks.len())],
        count: blocks.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn records(size: usize, blocks: &[Vec<u8>]) -> Records {
        let data: Vec<u8> = blocks.iter().flatten().copied().collect();
        Records {
            name: "Test".into(),
            size,
            shift: 15,
            parts: vec![(data, blocks.len())],
            count: blocks.len(),
        }
    }

    fn string_block(id: i32, text: &[u8], next: i32) -> Vec<u8> {
        let mut block = id.to_le_bytes().to_vec();
        let mut padded = text.to_vec();
        padded.resize(60, STRING_PADDING);
        block.extend_from_slice(&padded);
        block.extend_from_slice(&next.to_le_bytes());
        block
    }

    #[test]
    fn strings_follow_their_chain_of_blocks() {
        let long = [b'A'; 60];
        let store = records(
            68,
            &[
                vec![0; 68],
                string_block(1, &long, 2),
                string_block(1, b"BC", -1),
                string_block(3, b"short", -1),
            ],
        );
        store.check_ids().unwrap_err(); // block 2 carries its chain's id
        store.check_chain_ids().unwrap();
        assert_eq!(store.string(1).unwrap(), format!("{}BC", "A".repeat(60)));
        assert_eq!(store.string(3).as_deref(), Some("short"));
        assert_eq!(store.string(9), None);
    }

    #[test]
    fn chains_continue_inside_their_own_part() {
        // Two blocks per part (shift 1); part 0 runs past its two ids with a
        // continuation block, and part 1's chain points at its own block 1.
        let part0: Vec<u8> = [
            vec![0; 68],
            string_block(1, &[b'A'; 60], 2),
            string_block(1, b"B", -1),
        ]
        .concat();
        let part1: Vec<u8> = [string_block(2, &[b'C'; 60], 1), string_block(2, b"D", -1)].concat();
        let store = Records {
            name: "Test".into(),
            size: 68,
            shift: 1,
            parts: vec![(part0, 3), (part1, 2)],
            count: 5,
        };
        store.check_chain_ids().unwrap();
        assert_eq!(store.string(1).unwrap(), format!("{}B", "A".repeat(60)));
        assert_eq!(store.string(2).unwrap(), format!("{}D", "C".repeat(60)));
    }

    #[test]
    fn int_lists_end_at_int_min() {
        let mut block = 1_i32.to_le_bytes().to_vec();
        for value in [10, 9, INT_LIST_END, 77] {
            block.extend_from_slice(&value.to_le_bytes());
        }
        block.resize(256, 0);
        let store = records(256, &[vec![0; 256], block]);
        store.check_ids().unwrap();
        assert_eq!(store.int_list(1), Some(vec![10, 9]));
        assert_eq!(store.ids().collect::<Vec<_>>(), vec![1]);
    }
}
