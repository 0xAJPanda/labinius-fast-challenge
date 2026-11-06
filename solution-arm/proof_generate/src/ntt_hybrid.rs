use super::*;
use crate::ntt::{centered, for_blocks, primes, ROW};

const PA: usize = 0;
const PB: usize = 1;
const P0: u32 = PRIMES[PA];
const P1: u32 = PRIMES[PB];
pub const TAIL: u32 = 1024;
pub const MODULUS: u128 = P0 as u128 * P1 as u128 * TAIL as u128;
const _: () = assert!(MODULUS == 38_712_288_256);
const _: () = assert!(MODULUS >= 34_816_965_049, "the module-BKZ floor");
pub const CHUNK: usize = 44;
pub const TERMS: usize = 836;
const _: () = assert!(
    (CHUNK as u64) * (TERMS as u64) * (TAIL as u64 / 2) < (P0 as u64 * P1 as u64) / 2,
    "a chunk's tail sums are below half of p0 p1"
);

pub const ROWS: usize = 4;

pub struct HybridMatrix {
    pub elements: usize,
    pub data: Vec<i16>,
}

impl HybridMatrix {
    #[inline]
    pub fn row(&self, m: usize, r: usize) -> &[i16] {
        &self.data[(ROWS * m + r) * ROW..][..ROW]
    }
}

#[inline]
pub fn tail_centered(v: u128) -> i64 {
    let r = (v % TAIL as u128) as i64;
    if r >= TAIL as i64 / 2 {
        r - TAIL as i64
    } else {
        r
    }
}

pub fn matrix(elements: usize, workers: usize) -> HybridMatrix {
    let plans = primes();
    let mut data = vec![0i16; ROWS * elements * ROW];
    for_blocks(&mut data, ROWS * ROW, 8, workers, |first, block| {
        let mut coefficients = vec![[0u128; N]; 8];
        let mut lanes = vec![[[0i16; 8]; N]; ROWS];
        for (g, chunk) in block.chunks_mut(8 * ROWS * ROW).enumerate() {
            let base = first + 8 * g;
            let here = chunk.len() / (ROWS * ROW);
            for (e, c) in coefficients.iter_mut().enumerate().take(here) {
                crate::rq::ajtai_element(MODULUS, elements, base + e, c);
            }
            for r in 0..ROWS {
                let i = if r < 2 { PA } else { PB };
                let p = plans[i].p;
                for (e, c) in coefficients.iter().enumerate().take(here) {
                    for (slot, &v) in lanes[r].iter_mut().zip(c) {
                        slot[e] = if r % 2 == 0 {
                            centered((v % p as u128) as i64, p)
                        } else {
                            centered(tail_centered(v), p)
                        };
                    }
                }
                super::transform_general(i, &mut lanes[r]);
                for (e, row) in chunk.chunks_exact_mut(ROWS * ROW).enumerate().take(here) {
                    for (o, l) in row[r * ROW..][..N].iter_mut().zip(&lanes[r]) {
                        *o = centered(l[e] as i64, p);
                    }
                }
            }
        }
    });
    HybridMatrix { elements, data }
}

#[inline(always)]
unsafe fn pass_d2<const S: usize>(
    x: &mut [int16x8_t; N],
    kern: &Kernel,
    main_row: &[i16],
    tail_row: &[i16],
    acc: &mut [int32x4_t],
    tail: &mut [int32x4_t],
) {
    let s: Schedule = schedule(S);
    let k = vld1q_s16(kern.k.as_ptr());
    debug_assert!(main_row.len() == ROW && tail_row.len() == ROW);
    debug_assert!(acc.len() == 2 * N && tail.len() == 2 * N);
    for q in 0..N / 6 {
        let t = vld1q_s16(kern.d[q].as_ptr());
        let base = 6 * q;
        let mut v: [int16x8_t; 6] = core::array::from_fn(|i| load(x, base + i));
        reduce_inputs(&mut v, s.before[6], s.cheap, &[0, 3], k);
        let (y0, y1, y2) = r3!(s, 0, 1, 2, 3, v[0], v[1], v[2], t, k);
        let (y3, y4, y5) = r3!(s, 4, 5, 6, 7, v[3], v[4], v[5], t, k);
        let y = [y0, y1, y2, y3, y4, y5];
        mac6(acc, base, &y, vld1q_s16(main_row.as_ptr().add(base)));
        mac6(tail, base, &y, vld1q_s16(tail_row.as_ptr().add(base)));
    }
}

#[inline(always)]
unsafe fn narrow(acc: &[int32x4_t], out: &mut [int16x8_t; N]) {
    for (s, o) in out.iter_mut().enumerate() {
        *o = vcombine_s16(vmovn_s32(acc[2 * s]), vmovn_s32(acc[2 * s + 1]));
    }
}

#[inline(always)]
unsafe fn inverse_tail<const P: usize>(x: &mut [int16x8_t; N]) {
    super::inverse::<P>(x);
}

#[inline(always)]
unsafe fn canonical(x: int16x8_t, p: i16) -> int16x8_t {
    let h = vdupq_n_s16((p - 1) / 2);
    let pv = vdupq_n_s16(p);
    let over = vcgtq_s16(x, h);
    let under = vcltq_s16(x, vnegq_s16(h));
    let x = vsubq_s16(x, vandq_s16(vreinterpretq_s16_u16(over), pv));
    vaddq_s16(x, vandq_s16(vreinterpretq_s16_u16(under), pv))
}

fn inv_p0_pair() -> [i16; 8] {
    let inv = crate::ntt::pow_mod(P0 as u64 % P1 as u64, P1 as u64 - 2, P1 as u64);
    let t = super::twiddle(P1, inv);
    [t[0], t[1], 0, 0, 0, 0, 0, 0]
}

#[inline(always)]
unsafe fn recombine(r0: &[int16x8_t; N], r1: &[int16x8_t; N], tail16: &mut [int16x8_t; N], inv: int16x8_t, k1: int16x8_t) {
    let p0 = vdupq_n_s16(P0 as i16);
    for s in 0..N {
        let d = vsubq_s16(r1[s], r0[s]);
        let x1 = canonical(reduce(mulc::<0, 1>(d, inv, k1), k1), P1 as i16);
        let v = vmlaq_s16(r0[s], x1, p0);
        tail16[s] = vaddq_s16(tail16[s], v);
    }
}

#[inline(never)]
unsafe fn sums<const P: usize>(
    words: &[(*const u8, usize); LANES],
    a: &HybridMatrix,
    m: &mut [int16x8_t; BITS],
    x: &mut [int16x8_t; N],
    acc: &mut [int32x4_t],
    tail: &mut [int32x4_t],
    y: &mut [int16x8_t; N],
    mut chunk: impl FnMut(usize, &mut [int16x8_t; N]),
) {
    const { assert!(P < 2) };
    let kern = &kernels()[P];
    let s: Schedule = SCHEDULES[P];
    let p = PRIMES[P] as i32;
    let (rm, rt) = (2 * P, 2 * P + 1);
    acc.fill(vdupq_n_s32(0));
    tail.fill(vdupq_n_s32(0));
    let mut since = 0usize;
    for e in 0..a.elements {
        let at: [*const u8; LANES] = core::array::from_fn(|c| words[c].0.add(e * words[c].1));
        masks(&at, m);
        pass_a(m, x, kern);
        pass_b::<P>(x, kern);
        pass_c::<P>(x, kern);
        pass_d2::<P>(x, kern, a.row(e, rm), a.row(e, rt), acc, tail);
        if (e + 1) % s.accumulate == 0 {
            reduce_sums(acc, p, kern.v32);
        }
        since += 1;
        let end = (e + 1) % CHUNK == 0 || e + 1 == a.elements;
        if end || since == s.accumulate {
            reduce_sums(tail, p, kern.v32);
            since = 0;
        }
        if end {
            narrow(tail, y);
            inverse_tail::<P>(y);
            chunk(e / CHUNK, y);
            tail.fill(vdupq_n_s32(0));
        }
    }
    reduce_sums(acc, p, kern.v32);
}

#[inline]
fn crt2(r0: u32, r1: u32) -> u64 {
    const INV: u64 = crate::ntt::pow_mod(P0 as u64, P1 as u64 - 2, P1 as u64);
    let (r0, r1) = (r0 as u64, r1 as u64);
    let x1 = (r1 + P1 as u64 - r0 % P1 as u64) % P1 as u64 * INV % P1 as u64;
    r0 + P0 as u64 * x1
}

#[inline]
pub fn crt_tail(r01: u64, t: u32) -> u128 {
    const PP: u64 = P0 as u64 * P1 as u64;
    const INV: u64 = {
        let mut x = PP;
        let mut i = 0;
        while i < 5 {
            x = x.wrapping_mul(2u64.wrapping_sub(PP.wrapping_mul(x)));
            i += 1;
        }
        x % TAIL as u64
    };
    const _: () = assert!(PP.wrapping_mul(INV) % TAIL as u64 == 1);
    let k = ((t as u64 + TAIL as u64 - r01 % TAIL as u64) % TAIL as u64) * INV % TAIL as u64;
    (r01 + PP * k) as u128
}

pub fn commit(values: &[B128], a: &HybridMatrix, workers: usize) -> Vec<u128> {
    let len = 4 * a.elements;
    assert_eq!(values.len() % len, 0);
    let columns = values.len() / len;
    let chunks = a.elements.div_ceil(CHUNK);
    let mut out = vec![0u128; columns * N];
    for_blocks(&mut out, N, LANES, workers, |first, block| unsafe {
        let mut m = Box::new([vdupq_n_s16(0); BITS]);
        let mut x = Box::new([vdupq_n_s16(0); N]);
        let mut y = Box::new([vdupq_n_s16(0); N]);
        let mut acc = vec![vdupq_n_s32(0); 2 * 2 * N];
        let mut tail = vec![vdupq_n_s32(0); 2 * N];
        let mut held: Vec<[int16x8_t; N]> = vec![[vdupq_n_s16(0); N]; chunks];
        let mut tail16 = Box::new([vdupq_n_s16(0); N]);
        let inv = vld1q_s16(inv_p0_pair().as_ptr());
        let k1 = vld1q_s16(kernels()[PB].k.as_ptr());
        let count = block.len() / N;
        for (b, target) in block.chunks_mut(LANES * N).enumerate() {
            let col0 = first + LANES * b;
            let here = (count - LANES * b).min(LANES);
            let words: [(*const u8, usize); LANES] = core::array::from_fn(|c| {
                if c < here {
                    (values[(col0 + c) * len..].as_ptr() as *const u8, 64)
                } else {
                    (ZERO_WORDS.as_ptr() as *const u8, 0)
                }
            });
            let (a0, a1) = acc.split_at_mut(2 * N);
            sums::<PA>(&words, a, &mut m, &mut x, a0, &mut tail, &mut y, |c, r| {
                held[c] = *r;
            });
            tail16.fill(vdupq_n_s16(0));
            sums::<PB>(&words, a, &mut m, &mut x, a1, &mut tail, &mut y, |c, r| {
                recombine(&held[c], r, &mut tail16, inv, k1);
            });
            let lanes: Vec<[i32; 4]> =
                acc.iter().map(|&v| core::mem::transmute::<int32x4_t, [i32; 4]>(v)).collect();
            let t16: Vec<[i16; 8]> =
                tail16.iter().map(|&v| core::mem::transmute::<int16x8_t, [i16; 8]>(v)).collect();
            let plans = primes();
            for (c, column) in target.chunks_exact_mut(N).enumerate() {
                let mut residues = [[0u32; N]; 2];
                for (i, r) in residues.iter_mut().enumerate() {
                    let p = PRIMES[i] as i32;
                    for (slot, value) in r.iter_mut().enumerate() {
                        let v = lanes[i * 2 * N + 2 * slot + c / 4][c % 4];
                        *value = v.rem_euclid(p) as u32;
                    }
                }
                plans[PA].inverse(&mut residues[0]);
                plans[PB].inverse(&mut residues[1]);
                for (t, o) in column.iter_mut().enumerate() {
                    let r01 = crt2(residues[0][t], residues[1][t]);
                    let tl = (t16[t][c] as u16 as u32) % TAIL;
                    *o = crt_tail(r01, tl);
                }
            }
        }
    });
    out
}

#[cfg(test)]
mod bound_tests {
    use super::*;

    #[test]
    fn terms_bound_is_exact() {
        let mut weight = vec![0usize; N];
        for d in 0..N + 511 {
            let mut poly = vec![0i64; 2 * N];
            poly[d] = 1;
            for k in (N..2 * N).rev() {
                let c = poly[k];
                if c != 0 {
                    poly[k] = 0;
                    poly[k - 324] += c;
                    poly[k - N] -= c;
                }
            }
            let pairs = (0..N).filter(|&j| d >= j && d - j < 512).count();
            for (i, &c) in poly[..N].iter().enumerate() {
                weight[i] += pairs * c.unsigned_abs() as usize;
            }
        }
        assert_eq!(weight.iter().copied().max().unwrap(), TERMS);
    }
}
