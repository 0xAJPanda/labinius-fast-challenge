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
unsafe fn csa<const L: usize>(
    a: [uint8x16_t; L],
    b: [uint8x16_t; L],
    c: [uint8x16_t; L],
) -> ([uint8x16_t; L], [uint8x16_t; L]) {
    let mut hi = a;
    let mut lo = a;
    for i in 0..L {
        let u = veorq_u8(a[i], b[i]);
        hi[i] = vbslq_u8(u, c[i], a[i]);
        lo[i] = veorq_u8(u, c[i]);
    }
    (hi, lo)
}

/// Words `bases[l]` of the entry `e`'s column, complemented for a negative term.
#[inline(always)]
unsafe fn words<const L: usize>(bases: &[*const u8; L], e: u32) -> [uint8x16_t; L] {
    let off = (e & 0x7fff_ffff) as usize;
    let mask = vdupq_n_u8(((e as i32) >> 31) as u8);
    core::array::from_fn(|l| veorq_u8(vld1q_u8(bases[l].add(off)), mask))
}

/// The eleven count planes of the words at `bases[l]` (one word plane each) of the columns in
/// `list`. The `L` counters are independent dependency chains in one pass over the list.
#[inline(always)]
unsafe fn count<const L: usize, const P: usize, const H: usize>(
    bases: [*const u8; L],
    list: &[u32],
) -> [[uint8x16_t; P]; L] {
    const { assert!(H + 4 == P) };
    let z = [vdupq_n_u8(0); L];
    let (mut o, mut t2, mut f4, mut e8) = (z, z, z, z);
    let mut high = [z; H];
    let n = list.len();
    let mut i = 0;
    while i < n {
        let ld = |k: usize| {
            if i + k < n {
                words(&bases, *list.get_unchecked(i + k))
            } else {
                z
            }
        };
        let (ta, no) = csa(o, ld(0), ld(1));
        let (tb, no) = csa(no, ld(2), ld(3));
        let (fa, nt) = csa(t2, ta, tb);
        let (ta, no) = csa(no, ld(4), ld(5));
        let (tb, no) = csa(no, ld(6), ld(7));
        let (fb, nt) = csa(nt, ta, tb);
        let (ea, nf) = csa(f4, fa, fb);
        let (ta, no) = csa(no, ld(8), ld(9));
        let (tb, no) = csa(no, ld(10), ld(11));
        let (fa, nt) = csa(nt, ta, tb);
        let (ta, no) = csa(no, ld(12), ld(13));
        let (tb, no) = csa(no, ld(14), ld(15));
        let (fb, nt) = csa(nt, ta, tb);
        let (eb, nf) = csa(nf, fa, fb);
        let (s16, ne) = csa(e8, ea, eb);
        o = no;
        t2 = nt;
        f4 = nf;
        e8 = ne;
        let mut carry = s16;
        for h in high.iter_mut() {
            for l in 0..L {
                let c = vandq_u8(h[l], carry[l]);
                h[l] = veorq_u8(h[l], carry[l]);
                carry[l] = c;
            }
        }
        i += 16;
    }
    core::array::from_fn(|l| {
        let mut planes = [vdupq_n_u8(0); P];
        planes[0] = o[l];
        planes[1] = t2[l];
        planes[2] = f4[l];
        planes[3] = e8[l];
        for (p, h) in planes[4..].iter_mut().zip(&high) {
            *p = h[l];
        }
        planes
    })
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
unsafe fn add_counts<const P: usize>(planes: &[[uint8x16_t; P]; 4], bias: i16, raw: *mut i16) {
    let z = vdupq_n_u8(0);
    let mut lo = [[z; 8]; 4];
    let mut hi = [[z; 8]; 4];
    for k in 0..4 {
        let p = &planes[k];
        lo[k] = lane_order(transpose_bits([p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]));
        if P > 8 {
            let q = |i: usize| if i < P { p[i] } else { z };
            hi[k] = lane_order(transpose_bits([q(8), q(9), q(10), z, z, z, z, z]));
        }
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

/// Ring elements gathered together: each column then gives `64 BLOCK` contiguous bytes, not one
/// line per element at a 256 KiB stride.
pub const BLOCK: usize = 8;

impl Scratch {
    pub fn new(buckets: &Buckets) -> Self {
        Scratch(vec![0u128; BLOCK * 4 * buckets.columns])
    }
}

/// `v_m` for ring element `m`.
#[allow(dead_code)]
pub fn fold_element(
    values: &[B128],
    buckets: &Buckets,
    scratch: &mut Scratch,
    m: usize,
) -> RingElement648 {
    let mut out = [RingElement648::zero(Representation::Coefficients)];
    fold_elements(values, buckets, scratch, m, &mut out);
    let [e] = out;
    e
}

/// `v_m` for the ring elements `m0 .. m0 + out.len()`, at most [`BLOCK`] of them.
pub fn fold_elements(
    values: &[B128],
    buckets: &Buckets,
    scratch: &mut Scratch,
    m0: usize,
    out: &mut [RingElement648],
) {
    let columns = buckets.columns;
    let count = out.len();
    assert!(count <= BLOCK && m0 + count <= RING_ELEMENTS);
    assert!(values.len() >= columns * COLUMN_LEN);
    assert_eq!(scratch.0.len(), BLOCK * 4 * columns);
    for j in 0..columns {
        let w = &values[j * COLUMN_LEN + 4 * m0..][..4 * count];
        for (i, x) in w.iter().enumerate() {
            scratch.0[i * columns + j] = x.0;
        }
    }
    for (e, slot) in out.iter_mut().enumerate() {
        let plane0 = scratch.0[4 * e * columns..].as_ptr() as *const u8;
        *slot = element(plane0, columns, buckets);
    }
}

/// One element from its four gathered word planes at `base`.
fn element(base: *const u8, columns: usize, buckets: &Buckets) -> RingElement648 {
    let mut raw = [0i16; 2 * rq::N];
    unsafe {
        for pos in 0..POSITIONS {
            let list = &buckets.entries[buckets.starts[pos]..buckets.starts[pos + 1]];
            if list.is_empty() {
                continue;
            }
            let plane = |k: usize| base.add(16 * k * columns);
            let at = raw.as_mut_ptr().add(4 * pos);
            let bias = buckets.bias[pos];
            if list.len() < 1 << 8 {
                // Counts below 256: eight planes, no high byte.
                let [a, b] = count::<2, 8, 4>([plane(0), plane(1)], list);
                let [c, d] = count::<2, 8, 4>([plane(2), plane(3)], list);
                add_counts(&[a, b, c, d], bias, at);
            } else {
                let [a, b] = count::<2, PLANES, { PLANES - 4 }>([plane(0), plane(1)], list);
                let [c, d] = count::<2, PLANES, { PLANES - 4 }>([plane(2), plane(3)], list);
                add_counts(&[a, b, c, d], bias, at);
            }
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
