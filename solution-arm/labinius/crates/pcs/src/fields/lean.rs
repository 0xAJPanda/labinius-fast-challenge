//! The field switch without the full-size sumcheck polynomials.
//!
//! [`SwitchProver`](super::crossfield::SwitchProver) materialises `a[j] = L(eq[j])` and
//! `p[j] = lift(pi0[j])` at full length (2^24 entries of 24 bytes each) and halves both every
//! round. [`LeanSwitch`] sends the same messages from much less:
//!
//! - `eq` is a tensor product, `eq[q B + j] = pre[q] E[j]` with `E` the expansion of the last
//!   `l - d` coordinates (`B = 2^(l - d)` entries) and `pre` that of the first `d`. One pass over
//!   the words computes, per block `q < 2^d`, the 128 partial sums `V_q[i] = sum_j bit_i(pi0[q B +
//!   j]) E[j]` (the four-Russians kernel against the small `E`). Collapsing the blocks level by
//!   level with the coordinate factors `f_m(0) = 1 + r_m`, `f_m(1) = r_m` gives the sets of every
//!   prefix, and the root set is the partials the proof carries.
//! - `a_k[j] = L_k(E_k[j])` for a GF(2)-linear `L_k : B128 -> F162` kept as its 128 basis images,
//!   `E_k` the expansion of coordinates `k..l`. The fold by `s_k` gives
//!   `L_{k+1}(x) = L_k((1 + r_k) x) + s_k L_k(x)`.
//! - `p_k[j] = sum_z w_k(z) lift(pi0[z, j])` over the bound prefixes `z`, with `w_k(z)` the
//!   `F162` eq weights of the challenges.
//! - For `k < d` both message terms are bilinear in the bits of `pi0` and of `E_{k+1}`, so they
//!   are read off the prefix sets: `sum_j M(E[j]) lift(P[j]) = sum_i X^i M(U[i])` with
//!   `U[i] = sum_j bit_i(P[j]) E[j]`.
//! - After `d` rounds `p_d` is built in one more pass over the words (PMULL, unreduced, reduced
//!   once per entry) and `a_d` from `E`; the remaining `l - d` rounds run on `2^(l - d)` entries.
//!
//! The messages, the final evaluation and the partials are the field elements the full prover
//! computes (the sums are only reassociated), so the transcript and the proof are unchanged.
use super::crossfield::PACK;
use super::scalar::{B128, F162};
use super::sumcheck::Poly;
use crate::par;

/// Prefix coordinates the partial sets cover.
pub const PREFIX: usize = 6;

/// `x * X` in `F162` (`X^162 = X^81 + 1`).
#[inline]
fn mul_x(a: F162) -> F162 {
    let [l0, l1, l2] = a.0;
    let top = (l2 >> 33) & 1;
    let n0 = l0 << 1;
    let n1 = (l1 << 1) | (l0 >> 63);
    let n2 = ((l2 << 1) | (l1 >> 63)) & ((1u64 << 34) - 1);
    let mut r = F162([n0, n1, n2]);
    if top == 1 {
        r.0[0] ^= 1;
        r.0[1] ^= 1 << (81 - 64);
    }
    r
}

/// A linear map `B128 -> F162` as per-byte tables of its basis images.
struct Tables(Vec<F162>);

impl Tables {
    fn new(images: &[F162]) -> Self {
        assert_eq!(images.len(), PACK);
        let mut t = vec![F162::ZERO; 16 * 256];
        for b in 0..16 {
            for s in 1..256usize {
                let i = s.trailing_zeros() as usize;
                t[256 * b + s] = t[256 * b + (s & (s - 1))] + images[8 * b + i];
            }
        }
        Self(t)
    }

    #[inline]
    fn apply(&self, x: B128) -> F162 {
        let mut acc = F162::ZERO;
        for b in 0..16 {
            acc += self.0[256 * b + ((x.0 >> (8 * b)) & 0xff) as usize];
        }
        acc
    }

    /// `sum_i X^i M(u[i])`, by Horner from the top.
    fn phi(&self, u: &[B128]) -> F162 {
        let mut acc = F162::ZERO;
        for &x in u.iter().rev() {
            acc = mul_x(acc) + self.apply(x);
        }
        acc
    }
}

/// The images of the basis under `x -> L(c x)`.
fn compose_scale(l: &Tables, c: B128) -> Vec<F162> {
    (0..PACK).map(|b| l.apply(B128(1u128 << b) * c)).collect()
}

/// The factor of coordinate value `bit`: `r` for one, `1 + r` for zero.
#[inline]
fn factor(r: B128, bit: usize) -> B128 {
    if bit == 1 {
        r
    } else {
        B128::ONE + r
    }
}

pub struct LeanSwitch<'a> {
    pi0: &'a [B128],
    r_hi: &'a [B128],
    /// The expansion of `r_hi[PREFIX..]`.
    e: Vec<B128>,
    /// `sets[m]`: the `2^m` prefix sets of length `m`, 128 sums each, prefix MSB first.
    sets: Vec<Vec<B128>>,
    /// `L_k`, as basis images.
    l: Vec<F162>,
    /// `w_k(z)` for the `2^k` bound prefixes.
    w: Vec<F162>,
    a: Option<Poly>,
    p: Option<Poly>,
    round: usize,
    half: usize,
    threads: usize,
}

impl<'a> LeanSwitch<'a> {
    /// The expansion, the prefix sets and the partials; [`partials`](Self::partials) is the
    /// vector [`Workspace::partial_evals`](super::crossfield::Workspace::partial_evals) returns.
    pub fn new(pi0: &'a [B128], r_hi: &'a [B128], threads: usize) -> Self {
        let l = r_hi.len();
        let n = pi0.len();
        assert!(l > PREFIX && n == 1 << l, "2^l words, more coordinates than the prefix");
        let blocks = 1usize << PREFIX;
        let b = n >> PREFIX;
        let mut e = vec![B128::ZERO; b];
        super::crossfield::eq_expand_b128_into(&r_hi[PREFIX..], &mut e);
        let threads = threads.max(1);
        let parts = threads.min(b / 64).max(1);
        let ev = &e;
        let per: Vec<Vec<B128>> = par::map_blocks(b, parts, 64, |_, s, t| {
            let mut out = vec![B128::ZERO; blocks * PACK];
            for q in 0..blocks {
                let v = rows(&pi0[q * b + s..q * b + t], &ev[s..t]);
                out[q * PACK..][..PACK].copy_from_slice(&v);
            }
            out
        });
        let mut top = vec![B128::ZERO; blocks * PACK];
        for part in per {
            for (x, y) in top.iter_mut().zip(part) {
                *x = *x + y;
            }
        }
        let mut sets = vec![Vec::new(); PREFIX + 1];
        sets[PREFIX] = top;
        for m in (0..PREFIX).rev() {
            let (f0, f1) = (factor(r_hi[m], 0), factor(r_hi[m], 1));
            let child = &sets[m + 1];
            let mut cur = vec![B128::ZERO; (1 << m) * PACK];
            for p in 0..1usize << m {
                for i in 0..PACK {
                    cur[p * PACK + i] =
                        f0 * child[2 * p * PACK + i] + f1 * child[(2 * p + 1) * PACK + i];
                }
            }
            sets[m] = cur;
        }
        Self {
            pi0,
            r_hi,
            e,
            sets,
            l: Vec::new(),
            w: vec![F162::ONE],
            a: None,
            p: None,
            round: 0,
            half: n,
            threads,
        }
    }

    pub fn partials(&self) -> Vec<B128> {
        self.sets[0].clone()
    }

    /// Fix the batch weights `L_0(e_b) = batch[b]` before the first message.
    pub fn start(&mut self, batch: &[F162]) {
        assert_eq!(batch.len(), PACK);
        self.l = batch.to_vec();
    }

    pub fn msg(&mut self) -> [F162; 2] {
        self.half /= 2;
        let k = self.round;
        if k < PREFIX {
            let lt = Tables::new(&self.l);
            let m0 = Tables::new(&compose_scale(&lt, factor(self.r_hi[k], 0)));
            let sets = &self.sets[k + 1];
            let mut e0 = F162::ZERO;
            let mut einf = F162::ZERO;
            let mut sum = [B128::ZERO; PACK];
            for (z, &wz) in self.w.iter().enumerate() {
                let u0 = &sets[2 * z * PACK..][..PACK];
                let u1 = &sets[(2 * z + 1) * PACK..][..PACK];
                for i in 0..PACK {
                    sum[i] = u0[i] + u1[i];
                }
                e0 += wz * m0.phi(u0);
                einf += wz * lt.phi(&sum);
            }
            [e0, einf]
        } else {
            super::sumcheck::msg(self.a.as_ref().unwrap(), self.p.as_ref().unwrap(), self.half)
        }
    }

    pub fn fold(&mut self, s: F162) {
        let k = self.round;
        self.round += 1;
        if k < PREFIX {
            let lt = Tables::new(&self.l);
            let m0 = compose_scale(&lt, factor(self.r_hi[k], 0));
            self.l = m0.iter().zip(&self.l).map(|(&m, &x)| m + s * x).collect();
            self.w = self
                .w
                .iter()
                .flat_map(|&wz| [wz * (F162::ONE + s), wz * s])
                .collect();
            if k + 1 == PREFIX {
                self.materialize();
            }
        } else {
            super::sumcheck::fold(self.a.as_mut().unwrap(), self.p.as_mut().unwrap(), self.half, s);
        }
    }

    pub fn final_eval(&self) -> F162 {
        assert_eq!(self.half, 1);
        self.p.as_ref().unwrap().get(0)
    }

    /// `a_d = L_d(E)` and `p_d[j] = sum_z w_d(z) lift(pi0[z B + j])`.
    fn materialize(&mut self) {
        let b = self.e.len();
        let mut a = Poly::zero(b);
        super::crossfield::poly_from_eq_pub(&[&self.e], &[self.l.clone()], &mut a);
        let mut p = Poly::zero(b);
        let parts = self.threads.min(b / 64).max(1);
        let edges: Vec<(usize, usize)> = (0..parts).map(|t| par::block(b, parts, t, 8)).collect();
        let [w0, w1, w2] = &mut p.w;
        let (c0, c1, c2) = (
            par::split_mut(w0, &edges),
            par::split_mut(w1, &edges),
            par::split_mut(w2, &edges),
        );
        let (pi0, w) = (self.pi0, &self.w);
        std::thread::scope(|scope| {
            for (((x0, x1), x2), &(s, t)) in c0.into_iter().zip(c1).zip(c2).zip(&edges) {
                scope.spawn(move || weighted(pi0, b, w, s, t, [x0, x1, x2]));
            }
        });
        self.a = Some(a);
        self.p = Some(p);
    }
}

/// `v[k] = sum_j bit_k(pi0[j]) eq[j]`.
#[cfg(target_arch = "aarch64")]
fn rows(pi0: &[B128], eq: &[B128]) -> Vec<B128> {
    super::neon::partial_evals_rows(pi0, eq)
}

#[cfg(not(target_arch = "aarch64"))]
fn rows(pi0: &[B128], eq: &[B128]) -> Vec<B128> {
    let mut v = vec![B128::ZERO; PACK];
    for (&p, &e) in pi0.iter().zip(eq) {
        for (k, vk) in v.iter_mut().enumerate() {
            if p.bit(k) {
                *vk = *vk + e;
            }
        }
    }
    v
}

/// `out[j - s] = sum_z w[z] lift(pi0[z b + j])` for `j` in `[s, t)`.
fn weighted(pi0: &[B128], b: usize, w: &[F162], s: usize, t: usize, out: [&mut [u64]; 3]) {
    #[cfg(target_arch = "aarch64")]
    if super::neon::pmull() {
        // SAFETY: the CPU has PMULL.
        unsafe { super::neon::weighted_lift_rows(pi0, b, w, s, t, out) };
        return;
    }
    let [o0, o1, o2] = out;
    for j in s..t {
        let mut acc = F162::ZERO;
        for (z, &wz) in w.iter().enumerate() {
            acc += wz * F162::from_b128(pi0[z * b + j]);
        }
        o0[j - s] = acc.0[0];
        o1[j - s] = acc.0[1];
        o2[j - s] = acc.0[2];
    }
}

#[cfg(test)]
mod tests {
    use super::super::crossfield::{eq_expand_f162, SwitchProver, Workspace};
    use super::*;

    fn rng(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn mul_x_matches_product() {
        let mut st = 77u64;
        let x = F162([1 << 1, 0, 0]);
        for _ in 0..100 {
            let a = F162([rng(&mut st), rng(&mut st), rng(&mut st) & ((1 << 34) - 1)]);
            assert_eq!(mul_x(a), a * x);
        }
    }

    #[test]
    fn lean_matches_full() {
        for l in [PREFIX + 1, PREFIX + 3, 14] {
            let n = 1usize << l;
            let mut st = 0x1234_5678_9abc_def1u64 ^ l as u64;
            let pi0: Vec<B128> = (0..n)
                .map(|_| B128((rng(&mut st) as u128) << 64 | rng(&mut st) as u128))
                .collect();
            let r_hi: Vec<B128> = (0..l)
                .map(|_| B128((rng(&mut st) as u128) << 64 | rng(&mut st) as u128))
                .collect();
            let rp: Vec<F162> = (0..7)
                .map(|_| F162([rng(&mut st), rng(&mut st), rng(&mut st) & ((1 << 34) - 1)]))
                .collect();
            let batch = eq_expand_f162(&rp);
            let rs: Vec<F162> = (0..l)
                .map(|_| F162([rng(&mut st), rng(&mut st), rng(&mut st) & ((1 << 34) - 1)]))
                .collect();
            let mut ws = Workspace::new(l, 1);
            let want_partials = ws.partial_evals(0, &pi0, &r_hi);
            let mut full = SwitchProver::from_workspace(ws, &pi0, &[F162::ONE], &batch);
            let mut want = Vec::new();
            for &r in &rs {
                want.push(full.msg());
                full.fold(r);
            }
            let want_final = full.final_eval();
            for threads in [1, 3, 8] {
                let mut lean = LeanSwitch::new(&pi0, &r_hi, threads);
                assert_eq!(lean.partials(), want_partials, "l {l}");
                lean.start(&batch);
                for (k, &r) in rs.iter().enumerate() {
                    assert_eq!(lean.msg(), want[k], "l {l} threads {threads} round {k}");
                    lean.fold(r);
                }
                assert_eq!(lean.final_eval(), want_final, "l {l} threads {threads}");
            }
        }
    }
}
