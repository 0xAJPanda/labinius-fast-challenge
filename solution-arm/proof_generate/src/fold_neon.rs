//! The fold `v_m = sum_j c_j W_{j,m}` bit-sliced on NEON, bucketed by challenge position.
//!
//! A challenge term at position `pos` with sign `s` adds `s W_{j,m}` shifted by `4 pos`, and
//! coefficient `4 p + k` of `W_{j,m}` is bit `p` of word `k`, so the shift moves `p` by `pos` and
//! keeps `k`. Group the terms by position: per position the contribution to `raw[4 (pos + p) + k]`
//! is `D_pos[k][p] = sum_{+} bit_p(w_{j,k}) - sum_{-} bit_p(w_{j,k})`, a count over a list of
//! columns that is the same for every ring element `m`.
//!
//! Each word is one 128-lane NEON register (lane `p`), so the count is a column of carry-save
//! (Harley–Seal) adders over registers: sixteen inputs per step into ones/twos/fours/eights, the
//! sixteens rippled into seven more planes (counts up to `2^11 - 1`). A negative term counts the
//! complemented word, `1 - b`, so the count exceeds `D` by the number of negative terms at that
//! position, a per-position bias subtracted at read-out. The eleven planes are turned back into
//! integers with an 8x8 bit transpose and byte zips, and added into an `i16` raw row in
//! coefficient order with de-interleaving loads (`ld4`/`st4`). The ring reduction is the
//! scalar one of the reference fold. Every step is exact.
use crate::{rq, COLUMN_LEN, RING_ELEMENTS};
use core::arch::aarch64::*;
use labinius::challenge::ShortChallenge;
use labinius::fields::scalar::B128;
use labinius::ring::Representation;
use labinius::RingElement648;

/// Positions whose shifted window stays inside the raw row: `4 pos + 512 <= 2 N`.
const POSITIONS: usize = (2 * rq::N - 4 * 128) / 4 + 1;
/// Planes of a position's counter: counts below `2^PLANES`.
const PLANES: usize = 11;

/// Per position, the columns of its terms as byte offsets into a gathered word plane
/// ([`Scratch`]), the sign in bit 31.
pub struct Buckets {
    columns: usize,
    entries: Vec<u32>,
    starts: Vec<usize>,
    bias: Vec<i16>,
}

impl Buckets {
    /// `None` when a term falls outside the window or a position has too many terms.
    pub fn of(challenges: &[ShortChallenge]) -> Option<Buckets> {
        let mut lists: Vec<Vec<u32>> = vec![Vec::new(); POSITIONS];
        if challenges.len() >= 1 << 26 {
            return None;
        }
        for (j, c) in challenges.iter().enumerate() {
            for t in 0..c.weight {
                let position = c.positions[t] as usize;
                if position >= POSITIONS {
                    return None;
                }
                let negative = ((c.signs >> t) & 1 == 1) ^ (position % 2 == 1);
                lists[position].push((j * 16) as u32 | (negative as u32) << 31);
            }
        }
        if lists.iter().any(|l| l.len() >= 1 << PLANES) {
            return None;
        }
        let mut entries = Vec::new();
        let mut starts = vec![0];
        let mut bias = Vec::with_capacity(POSITIONS);
        for l in &lists {
            bias.push(l.iter().filter(|&&e| e >> 31 != 0).count() as i16);
            entries.extend_from_slice(l);
            starts.push(entries.len());
        }
        Some(Buckets {
            columns: challenges.len(),
            entries,
            starts,
            bias,
        })
    }
}

#[inline(always)]
unsafe fn csa(a: uint8x16_t, b: uint8x16_t, c: uint8x16_t) -> (uint8x16_t, uint8x16_t) {
    let u = veorq_u8(a, b);
    (vbslq_u8(u, c, a), veorq_u8(u, c))
}

/// Word `k` of the entry `e`'s column, complemented for a negative term.
#[inline(always)]
unsafe fn word(base: *const u8, e: u32) -> uint8x16_t {
    let x = vld1q_u8(base.add((e & 0x7fff_ffff) as usize));
    veorq_u8(x, vdupq_n_u8(((e as i32) >> 31) as u8))
}

/// The eleven count planes of the words `k` of the columns in `list`, `base` at word `k` of
/// element `m` of column 0.
#[inline(always)]
unsafe fn count(base: *const u8, list: &[u32]) -> [uint8x16_t; PLANES] {
    let z = vdupq_n_u8(0);
    let (mut o, mut t2, mut f4, mut e8) = (z, z, z, z);
    let mut high = [z; PLANES - 4];
    let mut chunks = list.chunks_exact(16);
    let mut step = |x: [uint8x16_t; 16], o: &mut uint8x16_t, t2: &mut uint8x16_t, f4: &mut uint8x16_t, e8: &mut uint8x16_t| {
        let (ta, no) = csa(*o, x[0], x[1]);
        let (tb, no) = csa(no, x[2], x[3]);
        let (fa, nt) = csa(*t2, ta, tb);
        let (ta, no) = csa(no, x[4], x[5]);
        let (tb, no) = csa(no, x[6], x[7]);
        let (fb, nt) = csa(nt, ta, tb);
        let (ea, nf) = csa(*f4, fa, fb);
        let (ta, no) = csa(no, x[8], x[9]);
        let (tb, no) = csa(no, x[10], x[11]);
        let (fa, nt) = csa(nt, ta, tb);
        let (ta, no) = csa(no, x[12], x[13]);
        let (tb, no) = csa(no, x[14], x[15]);
        let (fb, nt) = csa(nt, ta, tb);
        let (eb, nf) = csa(nf, fa, fb);
        let (s16, ne) = csa(*e8, ea, eb);
        *o = no;
        *t2 = nt;
        *f4 = nf;
        *e8 = ne;
        let mut carry = s16;
        for h in high.iter_mut() {
            let c = vandq_u8(*h, carry);
            *h = veorq_u8(*h, carry);
            carry = c;
        }
    };
    for ch in &mut chunks {
        let x: [uint8x16_t; 16] = core::array::from_fn(|i| word(base, *ch.get_unchecked(i)));
        step(x, &mut o, &mut t2, &mut f4, &mut e8);
    }
    let rest = chunks.remainder();
    if !rest.is_empty() {
        let x: [uint8x16_t; 16] =
            core::array::from_fn(|i| if i < rest.len() { word(base, rest[i]) } else { z });
        step(x, &mut o, &mut t2, &mut f4, &mut e8);
    }
    let mut planes = [z; PLANES];
    planes[0] = o;
    planes[1] = t2;
    planes[2] = f4;
    planes[3] = e8;
    planes[4..].copy_from_slice(&high);
    planes
}

/// `out[t].byte[b].bit[i] = v[i].byte[b].bit[t]`: the 8x8 bit transpose per byte lane.
#[inline(always)]
unsafe fn transpose_bits(mut v: [uint8x16_t; 8]) -> [uint8x16_t; 8] {
    let m1 = vdupq_n_u8(0x55);
    for i in [0usize, 2, 4, 6] {
        let (x, y) = (v[i], v[i + 1]);
        v[i] = vbslq_u8(m1, x, vshlq_n_u8::<1>(y));
        v[i + 1] = vbslq_u8(m1, vshrq_n_u8::<1>(x), y);
    }
    let m2 = vdupq_n_u8(0x33);
    for i in [0usize, 1, 4, 5] {
        let (x, y) = (v[i], v[i + 2]);
        v[i] = vbslq_u8(m2, x, vshlq_n_u8::<2>(y));
        v[i + 2] = vbslq_u8(m2, vshrq_n_u8::<2>(x), y);
    }
    for i in 0..4 {
        let (x, y) = (v[i], v[i + 4]);
        v[i] = vsliq_n_u8::<4>(x, y);
        v[i + 4] = vsriq_n_u8::<4>(y, x);
    }
    v
}

/// The bytes `byte[p]` for `p` in lane order from the transposed planes, where `t[q].byte[b]`
/// holds lane `8 b + q`: `out[c].byte[i]` is lane `16 c + i`.
#[inline(always)]
unsafe fn lane_order(t: [uint8x16_t; 8]) -> [uint8x16_t; 8] {
    let z: [uint8x16_t; 8] = [
        vzip1q_u8(t[0], t[1]),
        vzip2q_u8(t[0], t[1]),
        vzip1q_u8(t[2], t[3]),
        vzip2q_u8(t[2], t[3]),
        vzip1q_u8(t[4], t[5]),
        vzip2q_u8(t[4], t[5]),
        vzip1q_u8(t[6], t[7]),
        vzip2q_u8(t[6], t[7]),
    ];
    let h = |x: uint8x16_t| vreinterpretq_u16_u8(x);
    // y[g][r]: lanes b in 4 r .. 4 r + 4 of the four rows 4 g .. 4 g + 4.
    let y = [
        [
            vzip1q_u16(h(z[0]), h(z[2])),
            vzip2q_u16(h(z[0]), h(z[2])),
            vzip1q_u16(h(z[1]), h(z[3])),
            vzip2q_u16(h(z[1]), h(z[3])),
        ],
        [
            vzip1q_u16(h(z[4]), h(z[6])),
            vzip2q_u16(h(z[4]), h(z[6])),
            vzip1q_u16(h(z[5]), h(z[7])),
            vzip2q_u16(h(z[5]), h(z[7])),
        ],
    ];
    let w = |x: uint16x8_t| vreinterpretq_u32_u16(x);
    let b = |x: uint32x4_t| vreinterpretq_u8_u32(x);
    let mut out = [vdupq_n_u8(0); 8];
    for r in 0..4 {
        out[2 * r] = b(vzip1q_u32(w(y[0][r]), w(y[1][r])));
        out[2 * r + 1] = b(vzip2q_u32(w(y[0][r]), w(y[1][r])));
    }
    out
}

/// `raw[4 (pos + p) + k] += count[k][p] - bias` for `p < 128`, `k < 4`.
#[inline(always)]
unsafe fn add_counts(planes: &[[uint8x16_t; PLANES]; 4], bias: i16, raw: *mut i16) {
    let z = vdupq_n_u8(0);
    let mut lo = [[z; 8]; 4];
    let mut hi = [[z; 8]; 4];
    for k in 0..4 {
        let p = &planes[k];
        lo[k] = lane_order(transpose_bits([p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]));
        hi[k] = lane_order(transpose_bits([p[8], p[9], p[10], z, z, z, z, z]));
    }
    let b = vdupq_n_s16(bias);
    for c in 0..8 {
        for half in 0..2 {
            let at = raw.add(4 * (16 * c + 8 * half));
            let mut r = vld4q_s16(at);
            let value = |k: usize| {
                let (l, h) = if half == 0 {
                    (vget_low_u8(lo[k][c]), vget_low_u8(hi[k][c]))
                } else {
                    (vget_high_u8(lo[k][c]), vget_high_u8(hi[k][c]))
                };
                let v = vorrq_u16(vmovl_u8(l), vshll_n_u8::<8>(h));
                vsubq_s16(vreinterpretq_s16_u16(v), b)
            };
            r.0 = vaddq_s16(r.0, value(0));
            r.1 = vaddq_s16(r.1, value(1));
            r.2 = vaddq_s16(r.2, value(2));
            r.3 = vaddq_s16(r.3, value(3));
            vst4q_s16(at, r);
        }
    }
}

/// A worker's copy of one ring element's words across the columns: word `k` of column `j` at
/// `16 (k columns + j)`. The columns sit `16 COLUMN_LEN` bytes apart in the table, so in place
/// all of one element's words share a cache set; gathered, the counts read them from L1.
pub struct Scratch(Vec<u128>);

impl Scratch {
    pub fn new(buckets: &Buckets) -> Self {
        Scratch(vec![0u128; 4 * buckets.columns])
    }
}

/// `v_m` for ring element `m`.
pub fn fold_element(
    values: &[B128],
    buckets: &Buckets,
    scratch: &mut Scratch,
    m: usize,
) -> RingElement648 {
    let columns = buckets.columns;
    assert!(m < RING_ELEMENTS && values.len() >= columns * COLUMN_LEN);
    assert_eq!(scratch.0.len(), 4 * columns);
    for j in 0..columns {
        let w = &values[j * COLUMN_LEN + 4 * m..][..4];
        for k in 0..4 {
            scratch.0[k * columns + j] = w[k].0;
        }
    }
    let mut raw = [0i16; 2 * rq::N];
    unsafe {
        let base = scratch.0.as_ptr() as *const u8;
        for pos in 0..POSITIONS {
            let list = &buckets.entries[buckets.starts[pos]..buckets.starts[pos + 1]];
            if list.is_empty() {
                continue;
            }
            let planes: [[uint8x16_t; PLANES]; 4] =
                core::array::from_fn(|k| count(base.add(16 * k * columns), list));
            add_counts(&planes, buckets.bias[pos], raw.as_mut_ptr().add(4 * pos));
        }
    }
    let mut wide = [0i32; 2 * rq::N];
    for (w, &r) in wide.iter_mut().zip(&raw) {
        *w = r as i32;
    }
    for k in (rq::N..2 * rq::N - 1).rev() {
        let value = wide[k];
        wide[k - rq::SPLIT] += value;
        wide[k - rq::N] -= value;
    }
    let mut element = RingElement648::zero(Representation::Coefficients);
    for (slot, &value) in element.v.iter_mut().zip(&wide[..rq::N]) {
        *slot = i16::try_from(value).expect("a fold coefficient does not fit the wire codec");
    }
    element
}
