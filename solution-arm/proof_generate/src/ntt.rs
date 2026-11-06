use labinius::fields::scalar::B128;
use std::sync::OnceLock;

pub const N: usize = 648;
const HALF: usize = N / 2;
const ORDER: usize = 1944;
const RADICES: [usize; 6] = [2, 2, 3, 3, 3, 3];
pub const PRIMES: [u32; 3] = [3889, 9721, 17497];
pub const MODULUS: u128 = PRIMES[0] as u128 * PRIMES[1] as u128 * PRIMES[2] as u128;
const _: () = assert!(MODULUS == 661_473_542_593);
const _: () = {
    let mut i = 0;
    while i < 3 {
        assert!((PRIMES[i] - 1) % ORDER as u32 == 0, "every prime is 1 mod 1944");
        assert!(is_prime(PRIMES[i]));
        assert!(PRIMES[i] < 1 << 15, "residues are i16 lanes");
        i += 1;
    }
};

const fn is_prime(p: u32) -> bool {
    let mut d = 2;
    while d * d <= p {
        if p % d == 0 {
            return false;
        }
        d += 1;
    }
    p > 1
}

const fn pow_mod(base: u64, mut e: u64, p: u64) -> u64 {
    let mut b = base % p;
    let mut r = 1 % p;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

const fn primitive_root(p: u32) -> u32 {
    let p = p as u64;
    let mut g = 2;
    loop {
        let z = pow_mod(g, (p - 1) / ORDER as u64, p);
        if pow_mod(z, 972, p) != 1 && pow_mod(z, 648, p) != 1 {
            return z as u32;
        }
        g += 1;
    }
}

pub const fn centered(x: i64, p: u32) -> i16 {
    let p = p as i64;
    let r = x.rem_euclid(p);
    (if r > p / 2 { r - p } else { r }) as i16
}

pub struct Prime {
    pub p: u32,
    pub pow: Vec<u32>,
    pub exps: Vec<Vec<u32>>,
    inv_radix: [u32; 4],
    inv_split: u32,
}

impl Prime {
    fn new(p: u32) -> Prime {
        let zeta = primitive_root(p) as u64;
        let pow: Vec<u32> = (0..ORDER as u64).map(|i| pow_mod(zeta, i, p as u64) as u32).collect();
        let mut exps = vec![vec![HALF as u32, (5 * HALF) as u32]];
        for &r in &RADICES {
            let next: Vec<u32> = exps
                .last()
                .unwrap()
                .iter()
                .flat_map(|&e| {
                    assert_eq!(e as usize % r, 0, "every split has its roots");
                    (0..r).map(move |k| (e as usize / r + k * ORDER / r) as u32)
                })
                .collect();
            exps.push(next);
        }
        let inv = |x: u64| pow_mod(x, p as u64 - 2, p as u64) as u32;
        let w = pow[HALF] as u64;
        let w5 = pow[5 * HALF] as u64;
        Prime {
            p,
            inv_radix: [0, 1, inv(2), inv(3)],
            inv_split: inv((w + p as u64 - w5) % p as u64),
            pow,
            exps,
        }
    }

    #[inline]
    pub fn root(&self, e: usize) -> u64 {
        self.pow[e % ORDER] as u64
    }

    #[cfg_attr(not(test), allow(dead_code))]
    pub fn leaves(&self) -> &[u32] {
        &self.exps[RADICES.len()]
    }

    pub fn forward(&self, x: &mut [u32; N]) {
        match self.p {
            3889 => self.forward_mod(x, 3889),
            9721 => self.forward_mod(x, 9721),
            17497 => self.forward_mod(x, 17497),
            p => self.forward_mod(x, p as u64),
        }
    }

    #[inline(always)]
    fn forward_mod(&self, x: &mut [u32; N], p: u64) {
        debug_assert_eq!(p, self.p as u64);
        let (w, w5) = (self.root(HALF), self.root(5 * HALF));
        for j in 0..HALF {
            let (a, b) = (x[j] as u64, x[j + HALF] as u64);
            x[j] = ((a + w * b) % p) as u32;
            x[j + HALF] = ((a + w5 * b) % p) as u32;
        }
        let mut len = HALF;
        for (level, &r) in RADICES.iter().enumerate() {
            let sub = len / r;
            for (node, &e) in self.exps[level].iter().enumerate() {
                let base = node * len;
                let tau = e as usize / r;
                for j in 0..sub {
                    let mut a = [0u64; 3];
                    for (i, slot) in a.iter_mut().enumerate().take(r) {
                        *slot = x[base + j + i * sub] as u64 * self.root(tau * i) % p;
                    }
                    for k in 0..r {
                        let mut acc = 0u64;
                        for (i, &v) in a.iter().enumerate().take(r) {
                            acc += v * self.root(k * i * (ORDER / r));
                        }
                        x[base + j + k * sub] = (acc % p) as u32;
                    }
                }
            }
            len = sub;
        }
    }

    pub fn inverse(&self, x: &mut [u32; N]) {
        match self.p {
            3889 => self.inverse_mod(x, 3889),
            9721 => self.inverse_mod(x, 9721),
            17497 => self.inverse_mod(x, 17497),
            p => self.inverse_mod(x, p as u64),
        }
    }

    #[inline(always)]
    fn inverse_mod(&self, x: &mut [u32; N], p: u64) {
        debug_assert_eq!(p, self.p as u64);
        let mut len = 1;
        for level in (0..RADICES.len()).rev() {
            let r = RADICES[level];
            let parent = len * r;
            for (node, &e) in self.exps[level].iter().enumerate() {
                let base = node * parent;
                let tau = e as usize / r;
                for j in 0..len {
                    let mut y = [0u64; 3];
                    for (k, slot) in y.iter_mut().enumerate().take(r) {
                        *slot = x[base + j + k * len] as u64;
                    }
                    for i in 0..r {
                        let mut acc = 0u64;
                        for (k, &v) in y.iter().enumerate().take(r) {
                            acc += v * self.root(ORDER - k * i * (ORDER / r) % ORDER);
                        }
                        let scale = self.root(ORDER - tau * i % ORDER) * self.inv_radix[r] as u64 % p;
                        x[base + j + i * len] = (acc % p * scale % p) as u32;
                    }
                }
            }
            len = parent;
        }
        let w = self.root(HALF);
        for j in 0..HALF {
            let (u, v) = (x[j] as u64, x[j + HALF] as u64);
            let b = (u + p - v) * self.inv_split as u64 % p;
            x[j] = ((u + p * p - w * b) % p) as u32;
            x[j + HALF] = b as u32;
        }
    }
}

pub fn primes() -> &'static [Prime; 3] {
    static PLANS: OnceLock<[Prime; 3]> = OnceLock::new();
    PLANS.get_or_init(|| PRIMES.map(Prime::new))
}

#[inline]
pub fn crt(r: [u32; 3]) -> u128 {
    const P0: u64 = PRIMES[0] as u64;
    const P1: u64 = PRIMES[1] as u64;
    const P2: u64 = PRIMES[2] as u64;
    const INV_P0_MOD_P1: u64 = pow_mod(P0, P1 - 2, P1);
    const INV_P0P1_MOD_P2: u64 = pow_mod(P0 * P1 % P2, P2 - 2, P2);
    let (r0, r1, r2) = (r[0] as u64, r[1] as u64, r[2] as u64);
    debug_assert!(r0 < P0 && r1 < P1 && r2 < P2);
    let x1 = (r1 + P1 - r0 % P1) % P1 * INV_P0_MOD_P1 % P1;
    let low = r0 + P0 * x1;
    let x2 = (r2 + P2 - low % P2) % P2 * INV_P0P1_MOD_P2 % P2;
    (low + P0 * P1 * x2) as u128
}

#[derive(Clone, Copy)]
pub struct Schedule {
    pub before: [u8; 7],
    pub reduce_d: bool,
    pub cheap: bool,
    pub accumulate: usize,
}

pub const NONE: u8 = 0;
pub const ALL: u8 = 1;
pub const FIRST: u8 = 2;

pub const SCHEDULES: [Schedule; 3] = [
    Schedule { before: [NONE; 7], reduce_d: false, cheap: false, accumulate: 35 },
    Schedule {
        before: [NONE, NONE, NONE, FIRST, FIRST, FIRST, FIRST],
        reduce_d: false,
        cheap: false,
        accumulate: 21,
    },
    Schedule {
        before: [NONE, NONE, FIRST, ALL, ALL, ALL, ALL],
        reduce_d: true,
        cheap: true,
        accumulate: 7,
    },
];

pub const GENERAL: Schedule = Schedule {
    before: [ALL; 7],
    reduce_d: true,
    cheap: false,
    accumulate: 0,
};

pub const fn schedule(index: usize) -> Schedule {
    if index < 3 {
        SCHEDULES[index]
    } else {
        GENERAL
    }
}

pub const fn general_transform_bound(p: u32) -> u32 {
    const MAX: u32 = i16::MAX as u32;
    let h = (p - 1) / 2;
    assert!(h + product_bound(p, h) <= MAX, "level 0: a + w b");
    let r = reduced_bound(p);
    assert!(r + product_bound(p, r) <= MAX, "level 1: a +- tau b");
    transform_bound(p, GENERAL)
}

pub const fn product_bound(p: u32, b: u32) -> u32 {
    ((p as u64 * (b as u64 + 32768)).div_ceil(65536)) as u32
}

pub const fn reduced_bound(p: u32) -> u32 {
    p.div_ceil(2) + (3 * p).div_ceil(4096)
}

pub const fn cheap_bound(p: u32, b: u32) -> u32 {
    let off = (p as i64 - (1 << 14)).unsigned_abs();
    ((b as u64 * off + (p as u64) * (1 << 13)).div_ceil(1 << 14)) as u32
}

const fn level_inputs(p: u32, b: u32, mode: u8, cheap: bool) -> (u32, u32) {
    match mode {
        ALL if cheap => (reduced_bound(p), cheap_bound(p, b)),
        ALL => (reduced_bound(p), reduced_bound(p)),
        FIRST => (reduced_bound(p), b),
        _ => (b, b),
    }
}

pub const fn transform_bound(p: u32, s: Schedule) -> u32 {
    const MAX: u32 = i16::MAX as u32;
    let h = (p - 1) / 2;
    let b = 1 + 3 * h;
    assert!(b <= MAX);
    let (a, x) = level_inputs(p, b, s.before[2], s.cheap);
    let mut b = a + product_bound(p, x);
    assert!(b <= MAX);
    let mut l = 3;
    while l <= 6 {
        let (a, x) = level_inputs(p, b, s.before[l], s.cheap);
        let m = product_bound(p, x);
        assert!(2 * m <= MAX, "b' + c' and b' - c' fit");
        let e = if s.reduce_d && s.cheap {
            product_bound(p, cheap_bound(p, 2 * m))
        } else if s.reduce_d {
            product_bound(p, reduced_bound(p))
        } else {
            product_bound(p, 2 * m)
        };
        let y0 = a + 2 * m;
        let y12 = a + m + e;
        b = if y0 > y12 { y0 } else { y12 };
        assert!(b <= MAX, "a level's outputs fit i16");
        l += 1;
    }
    b
}

pub const fn accumulator_bound(p: u32, s: Schedule) -> u64 {
    p as u64 + s.accumulate as u64 * ((p as u64 - 1) / 2) * transform_bound(p, s) as u64
}

const _: () = {
    let mut i = 0;
    while i < 3 {
        let p = PRIMES[i];
        transform_bound(p, SCHEDULES[i]);
        general_transform_bound(p);
        assert!(accumulator_bound(p, SCHEDULES[i]) <= i32::MAX as u64);
        assert!(reduced_bound(p) < (1 << 14));
        i += 1;
    }
};

pub const ROW: usize = 656;
const _: () = assert!(ROW >= 6 * (N / 6 - 1) + 8 && ROW % 8 == 0);

pub struct Matrix {
    pub elements: usize,
    pub data: Vec<i16>,
}

impl Matrix {
    #[inline]
    pub fn row(&self, m: usize, prime: usize) -> &[i16] {
        &self.data[(3 * m + prime) * ROW..][..ROW]
    }

    #[inline]
    #[cfg_attr(target_arch = "aarch64", allow(dead_code))]
    pub fn slots(&self, m: usize, prime: usize) -> &[i16] {
        &self.row(m, prime)[..N]
    }
}

pub fn for_blocks<T: Send>(
    out: &mut [T],
    unit: usize,
    align: usize,
    workers: usize,
    f: impl Fn(usize, &mut [T]) + Sync,
) {
    let groups = out.len() / unit;
    let workers = workers.clamp(1, groups.div_ceil(align).max(1));
    let edges: Vec<(usize, usize)> = (0..workers)
        .map(|t| {
            let (s, e) = labinius::par::block(groups, workers, t, align);
            (s * unit, e * unit)
        })
        .collect();
    let blocks = labinius::par::split_mut(out, &edges);
    let f = &f;
    std::thread::scope(|scope| {
        for (block, &(s, _)) in blocks.into_iter().zip(&edges) {
            scope.spawn(move || f(s / unit, block));
        }
    });
}

#[cfg_attr(not(test), allow(dead_code))]
pub fn transform_matrix(matrix: &[u128], workers: usize) -> Matrix {
    assert_eq!(matrix.len() % N, 0);
    let elements = matrix.len() / N;
    let plans = primes();
    let mut data = vec![0i16; 3 * elements * ROW];
    for_blocks(&mut data, 3 * ROW, 1, workers, |first, block| {
        for (i, element) in block.chunks_exact_mut(3 * ROW).enumerate() {
            let a = &matrix[(first + i) * N..][..N];
            for (plan, out) in plans.iter().zip(element.chunks_exact_mut(ROW)) {
                let mut x = [0u32; N];
                for (slot, &c) in x.iter_mut().zip(a) {
                    *slot = (c % plan.p as u128) as u32;
                }
                plan.forward(&mut x);
                for (o, &v) in out.iter_mut().zip(&x) {
                    *o = centered(v as i64, plan.p);
                }
            }
        }
    });
    Matrix { elements, data }
}

pub fn matrix(elements: usize, workers: usize) -> Matrix {
    let plans = primes();
    let mut data = vec![0i16; 3 * elements * ROW];
    for_blocks(&mut data, 3 * ROW, 8, workers, |first, block| {
        let mut coefficients = vec![[0u128; N]; 8];
        let mut lanes = vec![[[0i16; 8]; N]; 3];
        for (g, chunk) in block.chunks_mut(8 * 3 * ROW).enumerate() {
            let base = first + 8 * g;
            let here = chunk.len() / (3 * ROW);
            for (e, c) in coefficients.iter_mut().enumerate().take(here) {
                crate::rq::ajtai_element(MODULUS, elements, base + e, c);
            }
            for (i, plan) in plans.iter().enumerate() {
                for (e, c) in coefficients.iter().enumerate().take(here) {
                    for (slot, &v) in lanes[i].iter_mut().zip(c) {
                        slot[e] = centered((v as u64 % plan.p as u64) as i64, plan.p);
                    }
                }
                transform_lanes(i, &mut lanes[i], here);
                for (e, row) in chunk.chunks_exact_mut(3 * ROW).enumerate() {
                    for (o, l) in row[i * ROW..][..N].iter_mut().zip(&lanes[i]) {
                        *o = centered(l[e] as i64, plan.p);
                    }
                }
            }
        }
    });
    Matrix { elements, data }
}

fn transform_lanes(i: usize, x: &mut [[i16; 8]; N], count: usize) {
    #[cfg(target_arch = "aarch64")]
    {
        let _ = count;
        neon::transform_general(i, x);
    }
    #[cfg(not(target_arch = "aarch64"))]
    transform_lanes_scalar(i, x, count);
}

#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
fn transform_lanes_scalar(i: usize, x: &mut [[i16; 8]; N], count: usize) {
    let plan = &primes()[i];
    for e in 0..count {
        let mut v: [u32; N] =
            core::array::from_fn(|j| (x[j][e] as i64).rem_euclid(plan.p as i64) as u32);
        plan.forward(&mut v);
        for (slot, &y) in x.iter_mut().zip(&v) {
            slot[e] = y as i16;
        }
    }
}

#[inline]
#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
pub fn element_bits(words: &[B128], out: &mut [u32; N]) {
    out.fill(0);
    for (k, w) in words.iter().enumerate() {
        for q in 0..128 {
            out[4 * q + k] = ((w.0 >> q) & 1) as u32;
        }
    }
}

pub fn finish_column(acc: &[[u32; N]; 3], out: &mut [u128]) {
    let mut residues = [[0u32; N]; 3];
    for ((plan, r), a) in primes().iter().zip(residues.iter_mut()).zip(acc) {
        *r = *a;
        plan.inverse(r);
    }
    for (t, o) in out.iter_mut().enumerate() {
        *o = crt([residues[0][t], residues[1][t], residues[2][t]]);
    }
}

#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
pub fn scalar_commit(values: &[B128], a: &Matrix, workers: usize) -> Vec<u128> {
    let len = 4 * a.elements;
    assert_eq!(values.len() % len, 0);
    assert!((a.elements as u64) << 29 < u64::MAX >> 1);
    let columns = values.len() / len;
    let plans = primes();
    let mut out = vec![0u128; columns * N];
    for_blocks(&mut out, N, 1, workers, |first, block| {
        for (c, target) in block.chunks_exact_mut(N).enumerate() {
            let column = &values[(first + c) * len..][..len];
            let mut acc = [[0u64; N]; 3];
            let mut bits = [0u32; N];
            for m in 0..a.elements {
                element_bits(&column[4 * m..4 * m + 4], &mut bits);
                for (i, plan) in plans.iter().enumerate() {
                    let mut x = bits;
                    plan.forward(&mut x);
                    let p = plan.p as i64;
                    for ((s, &v), &w) in acc[i].iter_mut().zip(&x).zip(a.slots(m, i)) {
                        *s += v as u64 * (w as i64).rem_euclid(p) as u64;
                    }
                }
            }
            let mut reduced = [[0u32; N]; 3];
            for (i, plan) in plans.iter().enumerate() {
                for (r, &s) in reduced[i].iter_mut().zip(&acc[i]) {
                    *r = (s % plan.p as u64) as u32;
                }
            }
            finish_column(&reduced, target);
        }
    });
    out
}

#[cfg(target_arch = "aarch64")]
#[path = "ntt_neon.rs"]
pub mod neon;

pub fn commit(values: &[B128], a: &Matrix, workers: usize) -> Vec<u128> {
    #[cfg(target_arch = "aarch64")]
    {
        neon::commit(values, a, workers)
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        scalar_commit(values, a, workers)
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    pub fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn random_residues(seed: u64, p: u32) -> [u32; N] {
        let mut s = seed | 1;
        core::array::from_fn(|_| (xorshift(&mut s) % p as u64) as u32)
    }

    #[test]
    fn roots_and_leaves() {
        for plan in primes() {
            let p = plan.p as u64;
            assert_eq!(plan.root(1) * plan.root(ORDER - 1) % p, 1);
            assert_eq!(pow_mod(plan.root(1), ORDER as u64, p), 1);
            for d in [972u64, 648] {
                assert_ne!(pow_mod(plan.root(1), d, p), 1);
            }
            let mut seen = vec![false; ORDER];
            for &e in plan.leaves() {
                assert!(e % 2 == 1 && e % 3 != 0, "leaf {e} is not primitive");
                assert!(!seen[e as usize]);
                seen[e as usize] = true;
            }
            assert_eq!(plan.leaves().len(), N);
            assert_eq!((plan.root(HALF) + plan.root(5 * HALF)) % p, 1);
            assert_eq!(plan.root(HALF) * plan.root(5 * HALF) % p, 1);
        }
    }

    #[test]
    fn forward_is_evaluation() {
        for (n, plan) in primes().iter().enumerate() {
            let p = plan.p as u64;
            let x = random_residues(17 + n as u64, plan.p);
            let mut y = x;
            plan.forward(&mut y);
            for (i, &e) in plan.leaves().iter().enumerate() {
                let mut want = 0u64;
                for (j, &c) in x.iter().enumerate() {
                    want = (want + c as u64 * plan.root(e as usize * j)) % p;
                }
                assert_eq!(y[i] as u64, want, "prime {} slot {i}", plan.p);
            }
        }
    }

    #[test]
    fn round_trip() {
        for (n, plan) in primes().iter().enumerate() {
            let mut inputs = vec![[0u32; N], [plan.p - 1; N], [1u32; N]];
            for t in [0usize, 1, 323, 324, 511, 647] {
                let mut x = [0u32; N];
                x[t] = 1;
                inputs.push(x);
            }
            for s in 0..4 {
                inputs.push(random_residues(100 * n as u64 + s, plan.p));
            }
            for x in inputs {
                let mut y = x;
                plan.forward(&mut y);
                assert!(y.iter().all(|&v| v < plan.p));
                plan.inverse(&mut y);
                assert_eq!(y, x, "prime {}", plan.p);
            }
        }
    }

    #[test]
    fn product_matches_schoolbook() {
        for (n, plan) in primes().iter().enumerate() {
            for s in 0..3u64 {
                let a = random_residues(1000 + 10 * n as u64 + s, plan.p);
                let mut b = random_residues(2000 + 10 * n as u64 + s, plan.p);
                if s == 2 {
                    b = core::array::from_fn(|i| (i < 512 && i % 3 != 1) as u32);
                }
                let want = labinius::scalar::mul_mod_phi(&a, &b, plan.p as u16);
                let (mut x, mut y) = (a, b);
                plan.forward(&mut x);
                plan.forward(&mut y);
                let mut z: [u32; N] =
                    core::array::from_fn(|i| (x[i] as u64 * y[i] as u64 % plan.p as u64) as u32);
                plan.inverse(&mut z);
                assert_eq!(z, want, "prime {} case {s}", plan.p);
            }
        }
    }

    #[test]
    fn lane_transforms_agree() {
        let mut s = 5u64;
        for (i, &p) in PRIMES.iter().enumerate() {
            let h = (p as i64 - 1) / 2;
            for case in 0..3 {
                let x: [[i16; 8]; N] = core::array::from_fn(|_| {
                    core::array::from_fn(|l| match (case, l % 2) {
                        (0, 0) => h as i16,
                        (0, _) => -h as i16,
                        _ => ((xorshift(&mut s) % (2 * h + 1) as u64) as i64 - h) as i16,
                    })
                });
                let (mut a, mut b) = (x, x);
                transform_lanes(i, &mut a, 8);
                transform_lanes_scalar(i, &mut b, 8);
                for (u, v) in a.iter().flatten().zip(b.iter().flatten()) {
                    assert_eq!((*u as i64 - *v as i64).rem_euclid(p as i64), 0, "prime {p} case {case} u {u} v {v}");
                }
            }
        }
    }

    #[test]
    fn crt_is_exact() {
        let mut s = 99u64;
        let mut cases: Vec<u128> = vec![0, 1, MODULUS - 1, MODULUS - 2, 3888, 3889, 37_804_969];
        for _ in 0..10_000 {
            cases.push(((xorshift(&mut s) as u128) << 64 | xorshift(&mut s) as u128) % MODULUS);
        }
        for x in cases {
            let r = PRIMES.map(|p| (x % p as u128) as u32);
            assert_eq!(crt(r), x);
        }
    }
}
