//! The fold with register-resident output blocks and nibble lanes.
//!
//! Every challenge term adds a column's bit mask, shifted by `4 * position` coefficients, into
//! the coefficient sums. A column's mask holds one coefficient per nibble (`0` or `1`), so a
//! shift by `4 * position` coefficients is a shift by `2 * position` bytes and a 512-bit witness
//! element is 16 vectors. Terms are sorted by position; a block of 128 coefficients is built in
//! four registers by walking the terms that reach it, one unaligned load and one add per vector.
//! Nibble sums are widened into `i16` lanes every 15 terms.

use crate::WITNESS_COEFFICIENTS;
use labinius::challenge::ShortChallenge;
use labinius::fields::scalar::B128;

const MASK_BYTES: usize = WITNESS_COEFFICIENTS / 2;
pub const PAD: usize = 64;
pub const STRIDE: usize = MASK_BYTES + 2 * PAD;
pub const BLOCKS: usize = 10;
const BLOCK_COEFFICIENTS: usize = 128;
pub const SUMS: usize = BLOCKS * BLOCK_COEFFICIENTS;
const POSITIONS: usize = 162;
/// Terms a nibble lane takes before it must be widened.
const SEGMENT: usize = 15;
/// Columns whose masks are walked together, so that they stay in the first-level cache.
const TILE: usize = 64;

/// Terms of all challenges ordered by position and sign.
pub struct Plan {
    /// `starts[s][p]..starts[s][p + 1]` index `columns[s]`: the terms at position `p` whose mask
    /// is added (`s = 0`) or subtracted (`s = 1`).
    starts: [Vec<u32>; 2],
    columns: [Vec<u32>; 2],
    tiles: usize,
}

impl Plan {
    pub fn new(challenges: &[ShortChallenge]) -> Plan {
        let tiles = challenges.len().div_ceil(TILE);
        let mut lists: [Vec<Vec<u32>>; 2] =
            [vec![Vec::new(); tiles * POSITIONS], vec![Vec::new(); tiles * POSITIONS]];
        for (j, c) in challenges.iter().enumerate() {
            for t in 0..c.weight {
                let p = c.positions[t] as usize;
                let negative = ((c.signs >> t) & 1 == 1) ^ (p % 2 == 1);
                lists[negative as usize][(j / TILE) * POSITIONS + p].push((j * STRIDE + PAD) as u32);
            }
        }
        let mut starts = [vec![0u32; tiles * POSITIONS + 1], vec![0u32; tiles * POSITIONS + 1]];
        let mut columns = [Vec::new(), Vec::new()];
        for s in 0..2 {
            for i in 0..tiles * POSITIONS {
                columns[s].extend_from_slice(&lists[s][i]);
                starts[s][i + 1] = columns[s].len() as u32;
            }
        }
        Plan { starts, columns, tiles }
    }
}

/// Per-worker scratch: the expanded masks of every column, zero-padded on both sides.
pub fn scratch(columns: usize) -> Vec<u8> {
    vec![0u8; columns * STRIDE]
}

/// Mask byte `2p + h` holds coefficient `4p + 2h` in its low nibble and `4p + 2h + 1` in its
/// high one: bit `p` of word `2h` and of word `2h + 1`.
#[cfg(target_arch = "aarch64")]
#[inline(always)]
fn expand(words: &[B128], mask: &mut [u8; MASK_BYTES]) {
    use core::arch::aarch64::*;
    assert_eq!(words.len(), 4);
    unsafe {
        let bits: [u8; 16] = [1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128];
        let bits = vld1q_u8(bits.as_ptr());
        let w: [uint8x16_t; 4] = core::array::from_fn(|k| vld1q_u8(words.as_ptr().add(k) as *const u8));
        for c in 0..8u8 {
            let idx = vcombine_u8(vdup_n_u8(2 * c), vdup_n_u8(2 * c + 1));
            let bit = |k: usize| vshrq_n_u8::<7>(vtstq_u8(vqtbl1q_u8(w[k], idx), bits));
            let even = vsliq_n_u8::<4>(bit(0), bit(1));
            let odd = vsliq_n_u8::<4>(bit(2), bit(3));
            let out = mask.as_mut_ptr().add(32 * c as usize);
            vst1q_u8(out, vzip1q_u8(even, odd));
            vst1q_u8(out.add(16), vzip2q_u8(even, odd));
        }
    }
}

/// The `i16` coefficient sums of ring element `m` before the ring reduction.
#[cfg(target_arch = "aarch64")]
pub fn sums(
    values: &[B128],
    column_len: usize,
    plan: &Plan,
    m: usize,
    masks: &mut [u8],
    out: &mut [i16; SUMS],
) {
    use core::arch::aarch64::*;
    let columns = masks.len() / STRIDE;
    for j in 0..columns {
        let mask: &mut [u8; MASK_BYTES] =
            (&mut masks[j * STRIDE + PAD..][..MASK_BYTES]).try_into().unwrap();
        expand(&values[j * column_len + 4 * m..][..4], mask);
    }
    let base = masks.as_ptr();
    for b in 0..BLOCKS {
        let first = (32 * b).saturating_sub(127);
        let last = (32 * b + 31).min(POSITIONS - 1);
        unsafe {
            // even[v] and odd[v]: sixteen byte lanes of vector v as two i16 vectors each.
            let mut even = [[vdupq_n_s16(0); 2]; 4];
            let mut odd = [[vdupq_n_s16(0); 2]; 4];
            for s in 0..2 {
                let mut acc = [vdupq_n_u8(0); 4];
                let mut pending = 0usize;
                macro_rules! flush {
                    () => {
                        let low = vdupq_n_u8(0x0f);
                        for v in 0..4 {
                            let e = vandq_u8(acc[v], low);
                            let o = vshrq_n_u8::<4>(acc[v]);
                            if s == 0 {
                                even[v][0] = vreinterpretq_s16_u16(vaddw_u8(vreinterpretq_u16_s16(even[v][0]), vget_low_u8(e)));
                                even[v][1] = vreinterpretq_s16_u16(vaddw_high_u8(vreinterpretq_u16_s16(even[v][1]), e));
                                odd[v][0] = vreinterpretq_s16_u16(vaddw_u8(vreinterpretq_u16_s16(odd[v][0]), vget_low_u8(o)));
                                odd[v][1] = vreinterpretq_s16_u16(vaddw_high_u8(vreinterpretq_u16_s16(odd[v][1]), o));
                            } else {
                                even[v][0] = vreinterpretq_s16_u16(vsubw_u8(vreinterpretq_u16_s16(even[v][0]), vget_low_u8(e)));
                                even[v][1] = vreinterpretq_s16_u16(vsubw_high_u8(vreinterpretq_u16_s16(even[v][1]), e));
                                odd[v][0] = vreinterpretq_s16_u16(vsubw_u8(vreinterpretq_u16_s16(odd[v][0]), vget_low_u8(o)));
                                odd[v][1] = vreinterpretq_s16_u16(vsubw_high_u8(vreinterpretq_u16_s16(odd[v][1]), o));
                            }
                            acc[v] = vdupq_n_u8(0);
                        }
                        pending = 0;
                    };
                }
                for tile in 0..plan.tiles {
                  for p in first..=last {
                    let shift = 64 * b as isize - 2 * p as isize;
                    let at = tile * POSITIONS + p;
                    let from = plan.starts[s][at] as usize;
                    let to = plan.starts[s][at + 1] as usize;
                    for &col in &plan.columns[s][from..to] {
                        let src = base.offset(col as isize + shift);
                        acc[0] = vaddq_u8(acc[0], vld1q_u8(src));
                        acc[1] = vaddq_u8(acc[1], vld1q_u8(src.add(16)));
                        acc[2] = vaddq_u8(acc[2], vld1q_u8(src.add(32)));
                        acc[3] = vaddq_u8(acc[3], vld1q_u8(src.add(48)));
                        pending += 1;
                        if pending == SEGMENT {
                            flush!();
                        }
                    }
                  }
                }
                if pending > 0 {
                    flush!();
                }
            }
            let dst = out.as_mut_ptr().add(b * BLOCK_COEFFICIENTS);
            for v in 0..4 {
                vst1q_s16(dst.add(32 * v), vzip1q_s16(even[v][0], odd[v][0]));
                vst1q_s16(dst.add(32 * v + 8), vzip2q_s16(even[v][0], odd[v][0]));
                vst1q_s16(dst.add(32 * v + 16), vzip1q_s16(even[v][1], odd[v][1]));
                vst1q_s16(dst.add(32 * v + 24), vzip2q_s16(even[v][1], odd[v][1]));
            }
        }
    }
}
