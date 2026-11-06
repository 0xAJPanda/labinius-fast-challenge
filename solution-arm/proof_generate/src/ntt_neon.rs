use super::{for_blocks, primes, schedule, Matrix, Schedule, ALL, FIRST, N, PRIMES, ROW, SCHEDULES};
use core::arch::aarch64::*;
use labinius::fields::scalar::B128;
use std::sync::OnceLock;

const LANES: usize = 8;
const BITS: usize = 512;

pub struct Kernel {
    k: [i16; 8],
    a: [i16; 8],
    g: [i16; 8],
    b: [[i16; 16]; 4],
    c: [[i16; 16]; 24],
    d: [[i16; 8]; 108],
    v32: i32,
}

fn twiddle(p: u32, w: u64) -> [i16; 2] {
    let c = super::centered(w as i64, p) as i64;
    let twisted = (2 * c * 32768 + p as i64).div_euclid(2 * p as i64);
    assert!(twisted.abs() <= 1 << 14);
    [c as i16, twisted as i16]
}

impl Kernel {
    fn new(index: usize) -> Kernel {
        let plan = &primes()[index];
        let p = plan.p;
        let tw = |e: usize| twiddle(p, plan.root(e));
        let r3 = |e: u32| -> [i16; 4] {
            let t = e as usize / 3;
            let (x, y) = (tw(t), tw(2 * t));
            [x[0], x[1], y[0], y[1]]
        };
        let eps = tw(648);
        let v = (((1u64 << 27) / p as u64 + 1) / 2) as i16;
        let k = [p as i16, v, eps[0], eps[1], 0, 0, 0, 0];
        let c = |e: usize| super::centered(plan.root(e) as i64, p);
        let a = [c(324), c(1620), c(162), c(486), c(810), c(486), 0, 0];
        let mut g = [0i16; 8];
        for (i, e) in [324usize, 1620, 162, 810].into_iter().enumerate() {
            g[2 * i..2 * i + 2].copy_from_slice(&tw(e));
        }
        let mut b = [[0i16; 16]; 4];
        for (n, slot) in b.iter_mut().enumerate() {
            slot[..4].copy_from_slice(&r3(plan.exps[2][2 * n]));
            slot[4..8].copy_from_slice(&r3(plan.exps[2][2 * n + 1]));
            slot[8..10].copy_from_slice(&tw(plan.exps[1][n] as usize / 2));
        }
        let mut cc = [[0i16; 16]; 24];
        for (n, slot) in cc.iter_mut().enumerate() {
            slot[..4].copy_from_slice(&r3(plan.exps[3][n]));
            for k in 0..3 {
                slot[4 + 4 * k..8 + 4 * k].copy_from_slice(&r3(plan.exps[4][3 * n + k]));
            }
        }
        let mut d = [[0i16; 8]; 108];
        for (q, slot) in d.iter_mut().enumerate() {
            slot[..4].copy_from_slice(&r3(plan.exps[5][2 * q]));
            slot[4..].copy_from_slice(&r3(plan.exps[5][2 * q + 1]));
        }
        Kernel {
            k,
            a,
            g,
            b,
            c: cc,
            d,
            v32: ((1u64 << 32) / p as u64).div_ceil(2) as i32,
        }
    }
}

pub fn kernels() -> &'static [Kernel; 3] {
    static K: OnceLock<[Kernel; 3]> = OnceLock::new();
    K.get_or_init(|| [Kernel::new(0), Kernel::new(1), Kernel::new(2)])
}

#[inline(always)]
unsafe fn mulc<const W: i32, const T: i32>(a: int16x8_t, t: int16x8_t, k: int16x8_t) -> int16x8_t {
    let mut lo = vmulq_laneq_s16::<W>(a, t);
    let hi = vqrdmulhq_laneq_s16::<T>(a, t);
    core::arch::asm!(
        "mls {lo:v}.8h, {hi:v}.8h, {k:v}.h[0]",
        lo = inout(vreg) lo,
        hi = in(vreg) hi,
        k = in(vreg_low16) k,
        options(pure, nomem, nostack, preserves_flags)
    );
    lo
}

#[inline(always)]
pub unsafe fn reduce(a: int16x8_t, k: int16x8_t) -> int16x8_t {
    let t = vqdmulhq_laneq_s16::<1>(a, k);
    let t = vrshrq_n_s16::<11>(t);
    vmlsq_laneq_s16::<0>(a, t, k)
}

#[inline(always)]
pub unsafe fn reduce_cheap(a: int16x8_t, k: int16x8_t) -> int16x8_t {
    vmlsq_laneq_s16::<0>(a, vrshrq_n_s16::<14>(a), k)
}

#[inline(always)]
unsafe fn reduce_inputs<const L: usize>(
    v: &mut [int16x8_t; L],
    mode: u8,
    cheap: bool,
    first: &[usize],
    k: int16x8_t,
) {
    if mode == ALL {
        for (i, y) in v.iter_mut().enumerate() {
            *y = if cheap && !first.contains(&i) { reduce_cheap(*y, k) } else { reduce(*y, k) };
        }
    } else if mode == FIRST {
        for &i in first {
            v[i] = reduce(v[i], k);
        }
    }
}

macro_rules! r3 {
    ($s:expr, $t0:literal, $t1:literal, $t2:literal, $t3:literal, $a:expr, $b:expr, $c:expr, $t:expr, $k:expr) => {
        match ($s.reduce_d, $s.cheap) {
            (false, _) => radix3::<$t0, $t1, $t2, $t3, 0>($a, $b, $c, $t, $k),
            (true, false) => radix3::<$t0, $t1, $t2, $t3, 1>($a, $b, $c, $t, $k),
            (true, true) => radix3::<$t0, $t1, $t2, $t3, 2>($a, $b, $c, $t, $k),
        }
    };
}

#[inline(always)]
unsafe fn radix3<const T0: i32, const T1: i32, const T2: i32, const T3: i32, const RD: u8>(
    a: int16x8_t,
    b: int16x8_t,
    c: int16x8_t,
    t: int16x8_t,
    k: int16x8_t,
) -> (int16x8_t, int16x8_t, int16x8_t) {
    let b1 = mulc::<T0, T1>(b, t, k);
    let c1 = mulc::<T2, T3>(c, t, k);
    let s = vaddq_s16(b1, c1);
    let mut d = vsubq_s16(b1, c1);
    if RD == 1 {
        d = reduce(d, k);
    } else if RD == 2 {
        d = reduce_cheap(d, k);
    }
    let e = mulc::<2, 3>(d, k, k);
    (
        vaddq_s16(a, s),
        vaddq_s16(vsubq_s16(a, c1), e),
        vsubq_s16(vsubq_s16(a, b1), e),
    )
}

#[inline(always)]
unsafe fn masks(words: &[*const u8; LANES], out: &mut [int16x8_t; BITS]) {
    let bit: [uint16x8_t; 8] = core::array::from_fn(|i| vdupq_n_u16(1 << i));
    for k in 0..4 {
        let r: [uint8x16_t; 8] = core::array::from_fn(|c| vld1q_u8(words[c].add(16 * k)));
        let a: [uint8x16_t; 8] = core::array::from_fn(|i| {
            let (x, y) = (r[2 * (i % 4)], r[2 * (i % 4) + 1]);
            if i < 4 { vzip1q_u8(x, y) } else { vzip2q_u8(x, y) }
        });
        let b: [uint16x8_t; 8] = core::array::from_fn(|i| {
            let h = i / 4;
            let (x, y) = (
                vreinterpretq_u16_u8(a[4 * h + 2 * ((i % 4) / 2)]),
                vreinterpretq_u16_u8(a[4 * h + 2 * ((i % 4) / 2) + 1]),
            );
            if i % 2 == 0 { vzip1q_u16(x, y) } else { vzip2q_u16(x, y) }
        });
        let t: [uint8x16_t; 8] = core::array::from_fn(|v| {
            let h = v / 4;
            let z = (v % 4) / 2;
            let (x, y) = (vreinterpretq_u32_u16(b[4 * h + z]), vreinterpretq_u32_u16(b[4 * h + 2 + z]));
            vreinterpretq_u8_u32(if v % 2 == 0 { vzip1q_u32(x, y) } else { vzip2q_u32(x, y) })
        });
        for (v, &pair) in t.iter().enumerate() {
            for (half, column_bytes) in [vmovl_u8(vget_low_u8(pair)), vmovl_high_u8(pair)]
                .into_iter()
                .enumerate()
            {
                let byte = 2 * v + half;
                for i in 0..8 {
                    out[32 * byte + 4 * i + k] = vreinterpretq_s16_u16(vtstq_u16(column_bytes, bit[i]));
                }
            }
        }
    }
}

#[inline(always)]
unsafe fn load(x: &[int16x8_t; N], i: usize) -> int16x8_t {
    *x.get_unchecked(i)
}

#[inline(always)]
unsafe fn store(x: &mut [int16x8_t; N], i: usize, v: int16x8_t) {
    *x.get_unchecked_mut(i) = v;
}

#[inline(always)]
unsafe fn pass_a(m: &[int16x8_t; BITS], x: &mut [int16x8_t; N], kern: &Kernel) {
    let a = vld1q_s16(kern.a.as_ptr());
    let (cw, cw5) = (vdupq_laneq_s16::<0>(a), vdupq_laneq_s16::<1>(a));
    let (ca, caw) = (vdupq_laneq_s16::<2>(a), vdupq_laneq_s16::<3>(a));
    let (cb, cbw) = (vdupq_laneq_s16::<4>(a), vdupq_laneq_s16::<5>(a));
    for j in 0..162 {
        let (m0, m1, m2) = (m[j], m[j + 162], m[j + 324]);
        let u = vsubq_s16(m0, vandq_s16(m2, cw));
        let v = vsubq_s16(m0, vandq_s16(m2, cw5));
        let mut w = vandq_s16(m1, ca);
        let mut w2 = vandq_s16(m1, cb);
        if j < BITS - 486 {
            let m3 = m[j + 486];
            w = vaddq_s16(w, vandq_s16(m3, caw));
            w2 = vaddq_s16(w2, vandq_s16(m3, cbw));
        }
        store(x, j, vsubq_s16(u, w));
        store(x, j + 162, vaddq_s16(u, w));
        store(x, j + 324, vsubq_s16(v, w2));
        store(x, j + 486, vaddq_s16(v, w2));
    }
}

#[inline(always)]
unsafe fn pass_b<const S: usize>(x: &mut [int16x8_t; N], kern: &Kernel) {
    let s: Schedule = schedule(S);
    let k = vld1q_s16(kern.k.as_ptr());
    for n in 0..4 {
        let t0 = vld1q_s16(kern.b[n].as_ptr());
        let t1 = vld1q_s16(kern.b[n].as_ptr().add(8));
        let base = 162 * n;
        for j in 0..27 {
            let mut v: [int16x8_t; 6] = core::array::from_fn(|i| load(x, base + j + 27 * i));
            reduce_inputs(&mut v, s.before[2], s.cheap, &[0, 1, 2], k);
            for i in 0..3 {
                let b1 = mulc::<0, 1>(v[i + 3], t1, k);
                let a = v[i];
                v[i] = vaddq_s16(a, b1);
                v[i + 3] = vsubq_s16(a, b1);
            }
            reduce_inputs(&mut v, s.before[3], s.cheap, &[0, 3], k);
            let (y0, y1, y2) = r3!(s, 0, 1, 2, 3, v[0], v[1], v[2], t0, k);
            let (y3, y4, y5) = r3!(s, 4, 5, 6, 7, v[3], v[4], v[5], t0, k);
            for (i, y) in [y0, y1, y2, y3, y4, y5].into_iter().enumerate() {
                store(x, base + j + 27 * i, y);
            }
        }
    }
}

#[inline(always)]
unsafe fn pass_c<const S: usize>(x: &mut [int16x8_t; N], kern: &Kernel) {
    let s: Schedule = schedule(S);
    let k = vld1q_s16(kern.k.as_ptr());
    for n in 0..24 {
        let t0 = vld1q_s16(kern.c[n].as_ptr());
        let t1 = vld1q_s16(kern.c[n].as_ptr().add(8));
        let base = 27 * n;
        for j in 0..3 {
            let mut v: [int16x8_t; 9] = core::array::from_fn(|i| load(x, base + j + 3 * i));
            reduce_inputs(&mut v, s.before[4], s.cheap, &[0, 1, 2], k);
            for i in 0..3 {
                let (a, b, c) = r3!(s, 0, 1, 2, 3, v[i], v[i + 3], v[i + 6], t0, k);
                (v[i], v[i + 3], v[i + 6]) = (a, b, c);
            }
            reduce_inputs(&mut v, s.before[5], s.cheap, &[0, 3, 6], k);
            let r = [
                r3!(s, 4, 5, 6, 7, v[0], v[1], v[2], t0, k),
                r3!(s, 0, 1, 2, 3, v[3], v[4], v[5], t1, k),
                r3!(s, 4, 5, 6, 7, v[6], v[7], v[8], t1, k),
            ];
            for (g, (a, b, c)) in r.into_iter().enumerate() {
                store(x, base + j + 9 * g, a);
                store(x, base + j + 9 * g + 3, b);
                store(x, base + j + 9 * g + 6, c);
            }
        }
    }
}

#[inline(always)]
unsafe fn mac6(acc: &mut [int32x4_t], base: usize, y: &[int16x8_t; 6], av: int16x8_t) {
    macro_rules! slot {
        ($i:literal) => {{
            let lo = acc.get_unchecked_mut(2 * (base + $i));
            *lo = vmlsl_laneq_s16::<$i>(*lo, vget_low_s16(y[$i]), av);
            let hi = acc.get_unchecked_mut(2 * (base + $i) + 1);
            *hi = vmlsl_high_laneq_s16::<$i>(*hi, y[$i], av);
        }};
    }
    slot!(0);
    slot!(1);
    slot!(2);
    slot!(3);
    slot!(4);
    slot!(5);
}

#[inline(always)]
unsafe fn pass_d<const S: usize, const MAC: bool>(
    x: &mut [int16x8_t; N],
    kern: &Kernel,
    row: &[i16],
    acc: &mut [int32x4_t],
) {
    let s: Schedule = schedule(S);
    let k = vld1q_s16(kern.k.as_ptr());
    debug_assert!(!MAC || (row.len() == ROW && acc.len() == 2 * N));
    for q in 0..N / 6 {
        let t = vld1q_s16(kern.d[q].as_ptr());
        let base = 6 * q;
        let mut v: [int16x8_t; 6] = core::array::from_fn(|i| load(x, base + i));
        reduce_inputs(&mut v, s.before[6], s.cheap, &[0, 3], k);
        let (y0, y1, y2) = r3!(s, 0, 1, 2, 3, v[0], v[1], v[2], t, k);
        let (y3, y4, y5) = r3!(s, 4, 5, 6, 7, v[3], v[4], v[5], t, k);
        let y = [y0, y1, y2, y3, y4, y5];
        if MAC {
            let av = vld1q_s16(row.as_ptr().add(base));
            mac6(acc, base, &y, av);
        } else {
            for (i, &v) in y.iter().enumerate() {
                store(x, base + i, v);
            }
        }
    }
}

#[inline(never)]
#[cfg_attr(not(test), allow(dead_code))]
pub unsafe fn transform<const P: usize>(m: &[int16x8_t; BITS], x: &mut [int16x8_t; N]) {
    let kern = &kernels()[P];
    pass_a(m, x, kern);
    pass_b::<P>(x, kern);
    pass_c::<P>(x, kern);
    pass_d::<P, false>(x, kern, &[], &mut []);
}

#[inline(always)]
unsafe fn pass_a_general(x: &mut [int16x8_t; N], kern: &Kernel) {
    let g = vld1q_s16(kern.g.as_ptr());
    let k = vld1q_s16(kern.k.as_ptr());
    for j in 0..N / 2 {
        let (a, b) = (load(x, j), load(x, j + N / 2));
        store(x, j, reduce(vaddq_s16(a, mulc::<0, 1>(b, g, k)), k));
        store(x, j + N / 2, reduce(vaddq_s16(a, mulc::<2, 3>(b, g, k)), k));
    }
    for (base, j) in (0..162).map(|j| (0, j)).chain((0..162).map(|j| (324, j))) {
        let a = load(x, base + j);
        let b = load(x, base + j + 162);
        let t = if base == 0 { mulc::<4, 5>(b, g, k) } else { mulc::<6, 7>(b, g, k) };
        store(x, base + j, vaddq_s16(a, t));
        store(x, base + j + 162, vsubq_s16(a, t));
    }
}

pub fn transform_general(prime: usize, x: &mut [[i16; 8]; N]) {
    unsafe {
        let mut v: [int16x8_t; N] = core::array::from_fn(|j| vld1q_s16(x[j].as_ptr()));
        let kern = &kernels()[prime];
        pass_a_general(&mut v, kern);
        pass_b::<3>(&mut v, kern);
        pass_c::<3>(&mut v, kern);
        pass_d::<3, false>(&mut v, kern, &[], &mut []);
        for (o, y) in x.iter_mut().zip(&v) {
            vst1q_s16(o.as_mut_ptr(), *y);
        }
    }
}

#[inline(always)]
unsafe fn reduce_sums(acc: &mut [int32x4_t], p: i32, v: i32) {
    for a in acc.iter_mut() {
        let t = vqrdmulhq_n_s32(*a, v);
        *a = vmlsq_n_s32(*a, t, p);
    }
}

#[inline(never)]
unsafe fn sums<const P: usize>(
    words: &[(*const u8, usize); LANES],
    a: &Matrix,
    m: &mut [int16x8_t; BITS],
    x: &mut [int16x8_t; N],
    acc: &mut [int32x4_t],
) {
    const { assert!(P < 3) };
    let kern = &kernels()[P];
    let s: Schedule = SCHEDULES[P];
    let p = PRIMES[P] as i32;
    acc.fill(vdupq_n_s32(0));
    for e in 0..a.elements {
        let at: [*const u8; LANES] = core::array::from_fn(|c| words[c].0.add(e * words[c].1));
        masks(&at, m);
        pass_a(m, x, kern);
        pass_b::<P>(x, kern);
        pass_c::<P>(x, kern);
        pass_d::<P, true>(x, kern, a.row(e, P), acc);
        if (e + 1) % s.accumulate == 0 {
            reduce_sums(acc, p, kern.v32);
        }
    }
    reduce_sums(acc, p, kern.v32);
}

static ZERO_WORDS: [u128; 4] = [0; 4];

pub fn commit(values: &[B128], a: &Matrix, workers: usize) -> Vec<u128> {
    let len = 4 * a.elements;
    assert_eq!(values.len() % len, 0);
    let columns = values.len() / len;
    let mut out = vec![0u128; columns * N];
    for_blocks(&mut out, N, LANES, workers, |first, block| {
        let mut m = Box::new([unsafe { vdupq_n_s16(0) }; BITS]);
        let mut x = Box::new([unsafe { vdupq_n_s16(0) }; N]);
        let mut acc = vec![unsafe { vdupq_n_s32(0) }; 3 * 2 * N];
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
            let (a0, rest) = acc.split_at_mut(2 * N);
            let (a1, a2) = rest.split_at_mut(2 * N);
            unsafe {
                sums::<0>(&words, a, &mut m, &mut x, a0);
                sums::<1>(&words, a, &mut m, &mut x, a1);
                sums::<2>(&words, a, &mut m, &mut x, a2);
            }
            let lanes: Vec<[i32; 4]> =
                acc.iter().map(|&v| unsafe { core::mem::transmute::<int32x4_t, [i32; 4]>(v) }).collect();
            for (c, column) in target.chunks_exact_mut(N).enumerate() {
                let mut residues = [[0u32; N]; 3];
                for (i, r) in residues.iter_mut().enumerate() {
                    let p = PRIMES[i] as i32;
                    for (slot, value) in r.iter_mut().enumerate() {
                        let v = lanes[i * 2 * N + 2 * slot + c / 4][c % 4];
                        *value = v.rem_euclid(p) as u32;
                    }
                }
                super::finish_column(&residues, column);
            }
        }
    });
    out
}

#[cfg(test)]
mod tests {
    use super::super::{product_bound, reduced_bound, tests::xorshift};
    use super::*;

    fn lanes(v: int16x8_t) -> [i16; 8] {
        unsafe { core::mem::transmute(v) }
    }

    fn every_i16() -> impl Iterator<Item = int16x8_t> {
        (0..1 << 13).map(|i| unsafe {
            let a: [i16; 8] = core::array::from_fn(|l| (i * 8 + l as i32 - 32768) as i16);
            vld1q_s16(a.as_ptr())
        })
    }

    #[test]
    fn reduction_is_exact_and_bounded() {
        for (i, kern) in kernels().iter().enumerate() {
            let p = PRIMES[i] as i32;
            let k = unsafe { vld1q_s16(kern.k.as_ptr()) };
            let mut worst = 0;
            for a in every_i16() {
                let r = lanes(unsafe { reduce(a, k) });
                for (x, y) in lanes(a).into_iter().zip(r) {
                    assert_eq!((x as i32 - y as i32).rem_euclid(p), 0);
                    worst = worst.max((y as i32).abs());
                }
            }
            assert!(worst as u32 <= reduced_bound(p as u32), "prime {p}: {worst}");
        }
    }

    #[test]
    fn cheap_reduction_is_exact_and_bounded() {
        for (i, kern) in kernels().iter().enumerate() {
            if !SCHEDULES[i].cheap {
                continue;
            }
            let p = PRIMES[i] as i32;
            let k = unsafe { vld1q_s16(kern.k.as_ptr()) };
            for a in every_i16() {
                let r = lanes(unsafe { reduce_cheap(a, k) });
                for (x, y) in lanes(a).into_iter().zip(r) {
                    assert_eq!((x as i32 - y as i32).rem_euclid(p), 0);
                    let bound = super::super::cheap_bound(p as u32, x.unsigned_abs() as u32);
                    assert!(y.unsigned_abs() as u32 <= bound, "prime {p}: {x} -> {y}");
                }
            }
        }
    }

    #[test]
    fn twiddle_products_are_exact_and_bounded() {
        for (i, kern) in kernels().iter().enumerate() {
            let p = PRIMES[i] as i64;
            let k = unsafe { vld1q_s16(kern.k.as_ptr()) };
            let mut pairs: Vec<[i16; 2]> = vec![[kern.k[2], kern.k[3]]];
            for row in kern.b.iter().chain(kern.c.iter()) {
                pairs.extend(row.chunks_exact(2).map(|c| [c[0], c[1]]));
            }
            for row in kern.d.iter() {
                pairs.extend(row.chunks_exact(2).map(|c| [c[0], c[1]]));
            }
            for [w, wt] in pairs {
                let t = unsafe { vld1q_s16([w, wt, 0, 0, 0, 0, 0, 0].as_ptr()) };
                for a in every_i16() {
                    let r = lanes(unsafe { mulc::<0, 1>(a, t, k) });
                    for (x, y) in lanes(a).into_iter().zip(r) {
                        assert_eq!((x as i64 * w as i64 - y as i64).rem_euclid(p), 0);
                        assert!(y.unsigned_abs() as u32 <= product_bound(p as u32, x.unsigned_abs() as u32));
                    }
                }
            }
        }
    }

    #[test]
    fn transform_matches_scalar() {
        let mut state = 7u64;
        for case in 0..12 {
            let words: [[u128; 4]; 8] = core::array::from_fn(|c| {
                core::array::from_fn(|k| match (case, c) {
                    (0, _) => u128::MAX,
                    (1, _) => 0,
                    (2, 0) => (k == 0) as u128,
                    (2, 1) => ((k == 3) as u128) << 127,
                    (2, 2) => if k == 1 { 1 << 64 } else { 0 },
                    _ => (xorshift(&mut state) as u128) << 64 | xorshift(&mut state) as u128,
                })
            });
            let ptrs: [*const u8; 8] = core::array::from_fn(|c| words[c].as_ptr() as *const u8);
            let mut m = [unsafe { vdupq_n_s16(0) }; BITS];
            unsafe { masks(&ptrs, &mut m) };
            for (i, plan) in primes().iter().enumerate() {
                let mut x = [unsafe { vdupq_n_s16(0) }; N];
                unsafe {
                    match i {
                        0 => transform::<0>(&m, &mut x),
                        1 => transform::<1>(&m, &mut x),
                        _ => transform::<2>(&m, &mut x),
                    }
                }
                for c in 0..8 {
                    let b: Vec<B128> = words[c].iter().map(|&w| B128(w)).collect();
                    let mut s = [0u32; N];
                    super::super::element_bits(&b, &mut s);
                    plan.forward(&mut s);
                    for (slot, &want) in s.iter().enumerate() {
                        let got = lanes(x[slot])[c] as i64;
                        let bound = super::super::transform_bound(plan.p, SCHEDULES[i]);
                        assert!(got.unsigned_abs() <= bound as u64);
                        assert_eq!(
                            (got + want as i64).rem_euclid(plan.p as i64),
                            0,
                            "case {case} prime {} column {c} slot {slot}",
                            plan.p
                        );
                    }
                }
            }
        }
    }
}
