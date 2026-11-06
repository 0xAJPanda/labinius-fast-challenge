//! Precomputes the transformed public matrix of the default (hybrid) modulus.
//!
//! The matrix `A` depends on the fixed matrix seed and on `Q` alone, never on the instance seed
//! (`rq::MATRIX_SEED`), so it is public parameter data a prover may precompute. This script
//! expands it exactly as `rq::ajtai_element` does, lifts its residues modulo 1024 exactly as
//! `hybrid::matrix` does, transforms the four rows of every element modulo 3889 and 9721 with
//! the transform tree of `ntt::Prime`, checks the tail certificate of every chunk, and writes
//! the rows in `hybrid::Matrix` layout (little-endian `i16`, centered residues). The prover
//! embeds the image; `hybrid::tests::embedded_matrix_is_the_runtime_matrix` checks it against
//! the run-time construction. Ported from the x86 track's precomputed key (submission 68fe4b84).
use std::io::Write;

const N: usize = 648;
const HALF: usize = N / 2;
const ROW: usize = 656;
const ORDER: usize = 1944;
const RADICES: [usize; 6] = [2, 2, 3, 3, 3, 3];
const ELEMENTS: usize = 4096;
const PRODUCT: i64 = 3889 * 9721;
const MODULUS: u128 = PRODUCT as u128 * 1024;
const CHUNK: usize = 256;
const MATRIX_SEED: [u8; 32] = [0x5a; 32];
const MATRIX_DOMAIN: &[u8] = b"labinius-fast/ajtai-matrix/v1";

fn pow_mod(base: u64, mut e: u64, p: u64) -> u64 {
    let mut b = base % p;
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

fn primitive_root(p: u64) -> u64 {
    let mut g = 2;
    loop {
        let z = pow_mod(g, (p - 1) / ORDER as u64, p);
        if pow_mod(z, 972, p) != 1 && pow_mod(z, 648, p) != 1 {
            return z;
        }
        g += 1;
    }
}

struct Plan {
    p: u64,
    pow: Vec<u64>,
    exps: Vec<Vec<usize>>,
}

impl Plan {
    fn new(p: u64) -> Plan {
        let zeta = primitive_root(p);
        let pow = (0..ORDER as u64).map(|i| pow_mod(zeta, i, p)).collect();
        let mut exps = vec![vec![HALF, 5 * HALF]];
        for &r in &RADICES {
            let next = exps
                .last()
                .unwrap()
                .iter()
                .flat_map(|&e: &usize| {
                    assert_eq!(e % r, 0);
                    (0..r).map(move |k| e / r + k * ORDER / r)
                })
                .collect();
            exps.push(next);
        }
        Plan { p, pow, exps }
    }

    fn root(&self, e: usize) -> u64 {
        self.pow[e % ORDER]
    }

    /// `ntt::Prime::forward`, residues in `[0, p)`.
    fn forward(&self, x: &mut [u64; N]) {
        let p = self.p;
        let (w, w5) = (self.root(HALF), self.root(5 * HALF));
        for j in 0..HALF {
            let (a, b) = (x[j], x[j + HALF]);
            x[j] = (a + w * b) % p;
            x[j + HALF] = (a + w5 * b) % p;
        }
        let mut len = HALF;
        for (level, &r) in RADICES.iter().enumerate() {
            let sub = len / r;
            for (node, &e) in self.exps[level].iter().enumerate() {
                let base = node * len;
                let tau = e / r;
                for j in 0..sub {
                    let mut a = [0u64; 3];
                    for (i, slot) in a.iter_mut().enumerate().take(r) {
                        *slot = x[base + j + i * sub] * self.root(tau * i) % p;
                    }
                    for k in 0..r {
                        let mut acc = 0u64;
                        for (i, &v) in a.iter().enumerate().take(r) {
                            acc += v * self.root(k * i * (ORDER / r));
                        }
                        x[base + j + k * sub] = acc % p;
                    }
                }
            }
            len = sub;
        }
    }
}

fn centered(x: i64, p: u64) -> i16 {
    let p = p as i64;
    let r = x.rem_euclid(p);
    (if r > p / 2 { r - p } else { r }) as i16
}

/// `rq::ajtai_element(MODULUS, ELEMENTS, index)`.
fn element(index: usize) -> [u128; N] {
    let bits = 128 - (MODULUS - 1).leading_zeros();
    let width = bits.div_ceil(8) as usize;
    let mask = (1u128 << bits) - 1;
    let mut hasher = blake3::Hasher::new_keyed(&MATRIX_SEED);
    hasher.update(MATRIX_DOMAIN);
    hasher.update(&MODULUS.to_le_bytes());
    hasher.update(&(ELEMENTS as u64).to_le_bytes());
    hasher.update(&(index as u64).to_le_bytes());
    let mut reader = hasher.finalize_xof();
    let mut out = [0u128; N];
    for slot in out.iter_mut() {
        loop {
            let mut bytes = [0u8; 16];
            reader.fill(&mut bytes[..width]);
            let candidate = u128::from_le_bytes(bytes) & mask;
            if candidate < MODULUS {
                *slot = candidate;
                break;
            }
        }
    }
    out
}

/// `hybrid::lift_pairs` of the centered residues modulo 1024.
fn tail(a: &[u128; N]) -> [i16; N] {
    let mut t: [i32; N] = core::array::from_fn(|j| ((a[j] as i64 + 512) & 1023) as i32 - 512);
    let norm = |a: i32, b: i32| a.abs() + b.abs() + (a + b).abs();
    for j in 0..HALF {
        let (a, b) = (t[j], t[j + HALF]);
        let mut best = (a, b);
        for u in -1..=1 {
            for v in -1..=1 {
                let (x, y) = (a + 1024 * u, b + 1024 * v);
                if norm(x, y) < norm(best.0, best.1) {
                    best = (x, y);
                }
            }
        }
        (t[j], t[j + HALF]) = best;
    }
    t.map(|v| v as i16)
}

/// `hybrid::certificate` of one chunk's tails: the worst exact extremum of a chunk sum.
fn certificate(raw: &[[i16; N]]) -> i64 {
    let mut pos = [0i64; N];
    let mut neg = pos;
    for a in raw {
        let mut ap = [0i64; 649];
        let mut an = ap;
        let mut sp = [0i64; 325];
        let mut sn = sp;
        for j in 0..648 {
            let v = a[j] as i64;
            ap[j + 1] = ap[j] + v.max(0);
            an[j + 1] = an[j] + (-v).max(0);
        }
        for j in 0..324 {
            let v = a[j] as i64 + a[j + 324] as i64;
            sp[j + 1] = sp[j] + v.max(0);
            sn[j + 1] = sn[j] + (-v).max(0);
        }
        let range = |p: &[i64], n: &[i64], lo: usize, hi: usize| {
            if lo > hi { (0, 0) } else { (p[hi + 1] - p[lo], n[hi + 1] - n[lo]) }
        };
        for j in 0..648 {
            let mut add = |(p, n): (i64, i64), minus: bool| {
                if minus {
                    pos[j] += n;
                    neg[j] += p;
                } else {
                    pos[j] += p;
                    neg[j] += n;
                }
            };
            if j < 324 {
                add(range(&ap, &an, 0, j), false);
                if j <= 186 {
                    add(range(&sp, &sn, j + 137, 323), true);
                    add(range(&ap, &an, 324, 647), true);
                } else {
                    add(range(&ap, &an, j + 137, 647), true);
                }
            } else {
                add(range(&sp, &sn, j.saturating_sub(511), 323), false);
                add(range(&ap, &an, 324, j), false);
                if j < 511 {
                    add(range(&ap, &an, j - 187, 323), false);
                }
            }
        }
    }
    pos.into_iter().chain(neg).max().unwrap()
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let out = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("hybrid-key.bin");
    let plans = [Plan::new(3889), Plan::new(9721)];
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get()).min(16);
    let per = ELEMENTS.div_ceil(workers).div_ceil(CHUNK) * CHUNK;
    let mut data = vec![0i16; 4 * ELEMENTS * ROW];
    std::thread::scope(|scope| {
        for (block, first) in data.chunks_mut(per * 4 * ROW).zip((0..).step_by(per)) {
            let plans = &plans;
            scope.spawn(move || {
                let mut tails = Vec::with_capacity(per);
                for (i, rows) in block.chunks_exact_mut(4 * ROW).enumerate() {
                    let a = element(first + i);
                    let t = tail(&a);
                    for (p, plan) in plans.iter().enumerate() {
                        let mut x: [u64; N] = core::array::from_fn(|j| (a[j] % plan.p as u128) as u64);
                        plan.forward(&mut x);
                        for (o, &v) in rows[(2 * p) * ROW..][..N].iter_mut().zip(&x) {
                            *o = centered(v as i64, plan.p);
                        }
                        let mut x: [u64; N] =
                            core::array::from_fn(|j| (t[j] as i64).rem_euclid(plan.p as i64) as u64);
                        plan.forward(&mut x);
                        for (o, &v) in rows[(2 * p + 1) * ROW..][..N].iter_mut().zip(&x) {
                            *o = centered(v as i64, plan.p);
                        }
                    }
                    tails.push(t);
                }
                for chunk in tails.chunks(CHUNK) {
                    let worst = certificate(chunk);
                    assert!(2 * worst < PRODUCT, "exact CRT capacity exceeded: {worst}");
                }
            });
        }
    });
    let mut bytes = Vec::with_capacity(2 * data.len());
    for v in &data {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::File::create(out).unwrap().write_all(&bytes).unwrap();
}
