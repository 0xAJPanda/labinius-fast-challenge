//! Build-time public parameters, independent of the input table and its seed.
use crate::{ntt::N, sms::{self, CHUNK, DIGITS, STREAMS}};
use std::borrow::Cow;

const ELEMENTS: usize = 4096;
const CHUNKS: usize = ELEMENTS / CHUNK;
const DATA_PER_CHUNK: usize = CHUNK * N * STREAMS;
const SUMS_PER_CHUNK: usize = DIGITS * N;
const _: () = assert!(ELEMENTS % CHUNK == 0 && CHUNK % 8 == 0);

#[repr(C, align(64))]
struct Aligned<const L: usize>([u8; L]);

static DATA: Aligned<{2 * ELEMENTS * N * STREAMS}> =
    Aligned(*include_bytes!(concat!(env!("OUT_DIR"), "/sms-data.bin")));
static SUMS: Aligned<{4 * CHUNKS * SUMS_PER_CHUNK}> =
    Aligned(*include_bytes!(concat!(env!("OUT_DIR"), "/sms-sums.bin")));
static MARGINS: [u8; 8 * CHUNKS] =
    *include_bytes!(concat!(env!("OUT_DIR"), "/sms-margins.bin"));

pub fn matrix() -> sms::Matrix {
    #[cfg(target_endian = "little")]
    // SAFETY: the arrays have 64-byte alignment and exact, compile-time lengths;
    // every i16/i32 bit pattern is valid, and the storage lives for the process.
    let (data, sums) = unsafe {
        (std::slice::from_raw_parts(DATA.0.as_ptr().cast::<i16>(), DATA.0.len() / 2),
         std::slice::from_raw_parts(SUMS.0.as_ptr().cast::<i32>(), SUMS.0.len() / 4))
    };
    let chunks = (0..CHUNKS).map(|c| {
        #[cfg(target_endian = "little")]
        let (data, sums) = (
            Cow::Borrowed(&data[c * DATA_PER_CHUNK..][..DATA_PER_CHUNK]),
            Cow::Borrowed(&sums[c * SUMS_PER_CHUNK..][..SUMS_PER_CHUNK]));
        #[cfg(target_endian = "big")]
        let (data, sums) = (
            Cow::Owned(DATA.0[c * DATA_PER_CHUNK * 2..][..DATA_PER_CHUNK * 2]
                .chunks_exact(2).map(|v| i16::from_le_bytes(v.try_into().unwrap())).collect()),
            Cow::Owned(SUMS.0[c * SUMS_PER_CHUNK * 4..][..SUMS_PER_CHUNK * 4]
                .chunks_exact(4).map(|v| i32::from_le_bytes(v.try_into().unwrap())).collect()));
        let margin = f64::from_le_bytes(MARGINS[c * 8..][..8].try_into().unwrap());
        assert!(margin >= 8.5);
        sms::Chunk { elements: CHUNK, data, sums, margin }
    }).collect();
    sms::Matrix { elements: ELEMENTS, chunks }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_parameters_match_public_expansion() {
        let embedded = matrix();
        let generated = sms::matrix(ELEMENTS, 4);
        assert_eq!(embedded.elements, generated.elements);
        assert_eq!(embedded.chunks.len(), generated.chunks.len());
        for (a, b) in embedded.chunks.iter().zip(&generated.chunks) {
            assert_eq!(a.elements, b.elements);
            assert_eq!(a.data, b.data);
            assert_eq!(a.sums, b.sums);
            assert_eq!(a.margin.to_bits(), b.margin.to_bits());
        }
    }
}
