//! Stored MS/MS peaks.
//!
//! `Spectrum` 4 names a `PeakList` record (180 bytes): its id, 20
//! (intensity, m/z) pairs, and at 164 the head of a chain of
//! `PeakListMassIntensities` blocks (232 bytes: the chain's id, 28 more
//! pairs, the next block's index in its part; NaN after the last pair, 0 at
//! 164 for no chain). The 20 are Scaffold's picks and the chain everything
//! else, so a spectrum is both together. Scaffold keeps at most 300 peaks.
//! Unused slots of the 20 hold (0, 0), or (−1, −1) in 5.2.2. All f32 LE.

use super::Sf3Error;
use super::records::{Container, Records, read_f32, read_i32};

const SUMMARY_BLOCK: usize = 180;
const PEAK_BLOCK: usize = 232;
const SUMMARY_PAIRS: usize = 20;
const CHAIN_HEAD: usize = 164;
const EMPTY_SLOTS: [(f32, f32); 2] = [(0.0, 0.0), (-1.0, -1.0)];

pub(crate) struct PeakStore {
    summaries: Records,
    chains: Records,
}

impl PeakStore {
    /// `None` when the project stores no peaks.
    pub fn load(container: &Container) -> Result<Option<Self>, Sf3Error> {
        if container.optional_stream("PeakList")?.is_none() {
            return Ok(None);
        }
        Ok(Some(Self {
            summaries: container.records("PeakList", SUMMARY_BLOCK)?,
            chains: container.chains("PeakListMassIntensities", PEAK_BLOCK)?,
        }))
    }

    /// (m/z, intensity) pairs sorted by m/z.
    pub fn read(&self, peak_list: i32) -> Result<Vec<(f64, f32)>, &'static str> {
        let summary = self
            .summaries
            .get(peak_list)
            .ok_or("its PeakList record is missing")?;
        let mut pairs: Vec<(f32, f32)> = (0..SUMMARY_PAIRS)
            .filter_map(|k| pair(summary, 4 + 8 * k))
            .filter(|pair| !EMPTY_SLOTS.contains(pair))
            .collect();
        let head = read_i32(summary, CHAIN_HEAD).unwrap_or(0);
        if head != 0 {
            let blocks = self
                .chains
                .chain(head)
                .ok_or("its chain of further peaks is broken")?;
            'blocks: for block in blocks {
                for k in 0..(PEAK_BLOCK - 8) / 8 {
                    match pair(block, 4 + 8 * k) {
                        Some(found) => pairs.push(found),
                        None => break 'blocks,
                    }
                }
            }
        }
        if pairs
            .iter()
            .any(|&(intensity, mz)| !(mz.is_finite() && mz > 0.0 && intensity.is_finite()))
        {
            return Err("a peak has an invalid m/z or intensity");
        }
        let mut peaks: Vec<(f64, f32)> = pairs
            .into_iter()
            .map(|(intensity, mz)| (f64::from(mz), intensity))
            .collect();
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        Ok(peaks)
    }
}

/// An (intensity, m/z) pair at `offset`; `None` at a NaN, which ends a list.
fn pair(block: &[u8], offset: usize) -> Option<(f32, f32)> {
    let intensity = read_f32(block, offset)?;
    let mz = read_f32(block, offset + 4)?;
    (!intensity.is_nan() && !mz.is_nan()).then_some((intensity, mz))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sf3::records::test_records;

    fn summary(id: i32, pairs: &[(f32, f32)], head: i32) -> Vec<u8> {
        let mut block = id.to_le_bytes().to_vec();
        for k in 0..SUMMARY_PAIRS {
            let (intensity, mz) = pairs.get(k).copied().unwrap_or((0.0, 0.0));
            block.extend_from_slice(&intensity.to_le_bytes());
            block.extend_from_slice(&mz.to_le_bytes());
        }
        block.extend_from_slice(&head.to_le_bytes());
        block.resize(SUMMARY_BLOCK, 0);
        block
    }

    fn chain_block(id: i32, pairs: &[(f32, f32)], next: i32) -> Vec<u8> {
        let mut block = id.to_le_bytes().to_vec();
        for k in 0..28 {
            let (intensity, mz) = pairs.get(k).copied().unwrap_or((f32::NAN, f32::NAN));
            block.extend_from_slice(&intensity.to_le_bytes());
            block.extend_from_slice(&mz.to_le_bytes());
        }
        block.extend_from_slice(&next.to_le_bytes());
        block
    }

    #[test]
    fn a_spectrum_is_its_twenty_peaks_and_its_chain_sorted_by_mz() {
        // Record 1: two picked peaks, then a chain of two blocks (1 → 2)
        // with 28 + 1 more; record 2: picked peaks only, (−1, −1) padding.
        let full: Vec<(f32, f32)> = (0..28).map(|k| (10.0, 300.0 + k as f32)).collect();
        let store = PeakStore {
            summaries: test_records(
                SUMMARY_BLOCK,
                &[
                    vec![0; SUMMARY_BLOCK],
                    summary(1, &[(900.0, 750.25), (500.0, 120.5)], 1),
                    summary(2, &[(7.0, 410.0), (-1.0, -1.0), (-1.0, -1.0)], 0),
                ],
            ),
            chains: test_records(
                PEAK_BLOCK,
                &[
                    vec![0; PEAK_BLOCK],
                    chain_block(1, &full, 2),
                    chain_block(1, &[(3.0, 99.0)], -1),
                ],
            ),
        };
        let one = store.read(1).unwrap();
        assert_eq!(one.len(), 2 + 28 + 1);
        assert_eq!(one[..3], [(99.0, 3.0), (120.5, 500.0), (300.0, 10.0)]);
        assert_eq!(one.last(), Some(&(750.25, 900.0)));
        assert_eq!(store.read(2).unwrap(), vec![(410.0, 7.0)]);
        assert_eq!(store.read(3), Err("its PeakList record is missing"));
    }

    #[test]
    fn a_broken_chain_or_bad_peak_is_reported_not_guessed() {
        let store = PeakStore {
            summaries: test_records(
                SUMMARY_BLOCK,
                &[
                    vec![0; SUMMARY_BLOCK],
                    summary(1, &[(9.0, 100.0)], 5),
                    summary(2, &[(9.0, -3.0)], 0),
                ],
            ),
            chains: test_records(PEAK_BLOCK, &[vec![0; PEAK_BLOCK]]),
        };
        assert_eq!(store.read(1), Err("its chain of further peaks is broken"));
        assert_eq!(store.read(2), Err("a peak has an invalid m/z or intensity"));
    }
}
