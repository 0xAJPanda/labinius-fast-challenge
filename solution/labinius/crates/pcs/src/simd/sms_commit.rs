//! Commitment modulo 9721 * 2^22 with one witness NTT.
//!
//! Stream zero computes the ordinary residue modulo 9721. Six auxiliary streams
//! encode paired radix-13 digits of the public matrix modulo 2^22. Pairing
//! coefficients j and j+324 uses Z=X^324, Z^2-Z+1=0, and the quadratic norm
//! a^2+a*b+b^2 to select short residue pairs. Integer carries preserve both
//! coefficients exactly modulo 2^22; public position resolves digit ties.
//!
//! Each digit sum is reconstructed around its public mean from its residue modulo
//! 9721. The public matrix determines the exact second-moment bounds, and the build
//! script refuses margins below its unchanged floor. Larger certified chunks
//! reduce inverse-transform work. Every invocation still consumes its fresh table.
//!
//! Only public key data and reconstruction bounds are precomputed.
use super::*;
use crate::simd::ntt::gen_small;

pub const P: u16 = 9721;
pub const TAIL_BITS: u32 = 22;
pub const MODULUS: u128 = (P as u128) << TAIL_BITS;
/// Accumulation streams: the base residue and six paired digits.
pub const STREAMS: usize = 7;
const DIGIT_RADIX: i64 = 13;
const DIGIT_HALF: i64 = DIGIT_RADIX / 2;
const DIGIT_BOUND: i8 = (2 * DIGIT_RADIX / 3) as i8;
const PAIR_RESIDUES: usize = (DIGIT_RADIX * DIGIT_RADIX) as usize;
const TOP_DIGIT_BOUND: i8 = {
    let mut bound = 3 * (1i64 << (TAIL_BITS - 1));
    let mut d = 0;
    while d < DIGITS - 1 {
        bound = (bound + DIGIT_BOUND as i64) / DIGIT_RADIX;
        d += 1;
    }
    bound as i8
};
pub const DIGITS: usize = STREAMS - 1;
/// Common radix of the auxiliary digit streams.
pub const fn digit_radix(d: usize) -> i64 {
    let _ = d; DIGIT_RADIX
}
/// Weight of digit `d`, `DIGIT_RADIX^d`.
pub const fn digit_weight(d: usize) -> u32 {
    let mut w = 1u32;
    let mut i = 0;
    while i < d {
        w *= DIGIT_RADIX as u32;
        i += 1;
    }
    w
}
/// The balanced representatives cover `[-REP_LOW, 2^22 - REP_LOW)`, all of `Z / 2^22`.
const REP_LOW: i64 = {
    // Canonical scalar digits at their least value.
    let mut low = 0i64;
    let mut d = 0;
    while d < DIGITS {
        low += DIGIT_HALF * digit_weight(d) as i64;
        d += 1;
    }
    low
};
const _: () = {
    // The digits' range is contiguous and at least 2^22 wide.
    let high = REP_LOW;
    assert!(REP_LOW + high + 1 >= 1 << TAIL_BITS);
};
/// Pair batches per reconstruction chunk; the build checks the public matrix's margin.
pub const CHUNK_PAIRS: usize = 104;
/// Chunks of each half column: [`CHUNK_PAIRS`] pair batches each and a shorter last one, so
/// that every chunk lies in one half (the commitment's work units are halves).
pub const fn chunks_per_half(np: usize) -> usize {
    (np / 2).div_ceil(CHUNK_PAIRS)
}
/// Chunks of a column of `elements` ring elements.
pub const fn chunk_count(elements: usize) -> usize {
    2 * chunks_per_half(elements / 2)
}
/// The pair batches of chunk `c` of a column of `np` pair batches.
pub fn chunk_range(np: usize, c: usize) -> core::ops::Range<usize> {
    let (h, l) = (c / chunks_per_half(np), c % chunks_per_half(np));
    let start = h * (np / 2) + l * CHUNK_PAIRS;
    start..(start + CHUNK_PAIRS).min((h + 1) * (np / 2))
}
/// Digits reconstructed per chunk: all of them.
pub const CHUNKED: usize = DIGITS;
const _: () = assert!(CHUNKED % 2 == 0);
/// Entries of [`SmsKeys::mur`] for `elements` ring elements per column.
pub const fn mur_len(elements: usize) -> usize {
    CHUNKED * chunk_count(elements) * N
}
/// Pair batches transformed in lockstep per accumulator read-modify-write.
pub const LOCK: usize = 8;
const TAIL_MASK: u32 = (1u32 << TAIL_BITS) - 1;

/// Centered `2^16 mod q`, the weight [`reduce_vec_c`] gives the high half.
const fn r16c(q: u16) -> i32 {
    (65536 % q as i32) - q as i32
}
const fn acc_after_reduce_c(q: u16) -> i64 {
    32768 * (1 + r16c(q).unsigned_abs() as i64)
}
/// One dword-product step: two products of a lazily reduced transform lane and a centered key.
const fn pair_per_batch(q: u16) -> i64 {
    2 * w_bound(q) * a_bound(q)
}
/// Pair batches between two fold-backs: the largest power of two `P` with
/// `acc_after_reduce_c + P * pair_per_batch <= i32::MAX`.
const fn pair_period(q: u16) -> usize {
    let mut p = 1usize;
    while acc_after_reduce_c(q) + 2 * (p as i64) * pair_per_batch(q) <= i32::MAX as i64 {
        p *= 2;
    }
    p
}
/// [`red16`]'s bias: the largest multiple of q below 2^31.
const fn bias(q: u16) -> i64 {
    ((1i64 << 31) / q as i64) * q as i64
}
pub const PP: usize = pair_period(P);
const _: () = {
    assert!(acc_after_reduce_c(P) + (PP as i64) * pair_per_batch(P) <= bias(P));
    assert!(CHUNK_PAIRS % PP == 0 && (PP % LOCK == 0 || LOCK % PP == 0) && CHUNK_PAIRS % LOCK == 0 && LOCK <= 2 * PP);
    assert!((P as i32 - 1) <= gen_small::TwI::<P, true>::IN_BOUND);
};

const fn qinv16(q: u16) -> u32 {
    let q = q as u32;
    let mut x = q;
    let mut i = 0;
    while i < 4 {
        x = x.wrapping_mul(2u32.wrapping_sub(q.wrapping_mul(x))) & 0xffff;
        i += 1;
    }
    x
}
const fn red16_w(q: u16) -> (i32, i32) {
    let wl = (1u32 << 16) % q as u32;
    let wh = ((1u64 << 32) % q as u64) as u32;
    let qi = qinv16(q);
    (
        (wl | (wh << 16)) as i32,
        (((wl * qi) & 0xffff) | (((wh * qi) & 0xffff) << 16)) as i32,
    )
}
/// A representative of `x mod q` with `|r| < 2q`, lane-wise, for `|x| <= bias(q)`.
#[inline(always)]
unsafe fn red16<const Q: u16>(x: __m512i) -> __m512i {
    let (w, wp) = red16_w(Q);
    let y = _mm512_add_epi32(x, _mm512_set1_epi32(bias(Q) as i32));
    let m = _mm512_mullo_epi16(y, _mm512_set1_epi32(wp));
    let r = _mm512_sub_epi16(
        _mm512_mulhi_epu16(y, _mm512_set1_epi32(w)),
        _mm512_mulhi_epu16(m, _mm512_set1_epi16(Q as i16)),
    );
    _mm512_madd_epi16(r, _mm512_set1_epi16(1))
}
// finish_stream only reads accumulators after a fold. Their signed dword
// range is much smaller than red16's general contract. For s=floor(x/4096),
// round(s*13807/32768) approximates x/P closely enough that x-qhat*P is in
// (-P,P). Multiplication and subtraction still give the exact residue.
const FOLDED_RECIP: i16 = (((1i64 << 27) + P as i64 / 2) / P as i64) as i16;
const _: () = {
    assert!(LOCK == PP);
    let a = acc_after_reduce_c(P);
    let max_s = (a + 4095) / 4096;
    let delta = (FOLDED_RECIP as i64 * P as i64 - (1i64 << 27)).abs();
    let error = 4095 + (max_s * delta + 32767) / 32768 + (P as i64 + 1) / 2;
    assert!(max_s <= i16::MAX as i64 && FOLDED_RECIP > 0);
    assert!((max_s * FOLDED_RECIP as i64 + 16384) / 32768 < 32768);
    assert!(error < P as i64);
};

#[inline(always)]
unsafe fn red_folded(x: __m512i) -> __m512i {
    let s = _mm512_srai_epi32::<12>(x);
    let qhat = _mm512_mulhrs_epi16(s, _mm512_set1_epi16(FOLDED_RECIP));
    // The zero high-word weight discards the high words of s and qhat.
    let product = _mm512_madd_epi16(qhat, _mm512_set1_epi32(P as i32));
    _mm512_sub_epi32(x, product)
}

/// Canonical representative in `[0, Q)` via [`red16`].
#[inline(always)]
unsafe fn canon_red16<const Q: u16>(x: __m512i) -> __m512i {
    let r = red16::<Q>(x);
    let q = _mm512_set1_epi32(Q as i32);
    let q2 = _mm512_set1_epi32(2 * Q as i32);
    let r = _mm512_mask_add_epi32(r, _mm512_cmp_epi32_mask(r, _mm512_setzero_si512(), _MM_CMPINT_LT), r, q2);
    _mm512_mask_sub_epi32(r, _mm512_cmp_epi32_mask(r, q, _MM_CMPINT_NLT), r, q)
}
fn zeroed_slice<T>(n: usize) -> Box<[T]> {
    unsafe {
        let l = std::alloc::Layout::array::<T>(n).unwrap();
        let p = std::alloc::alloc_zeroed(l) as *mut T;
        if p.is_null() {
            std::alloc::handle_alloc_error(l)
        }
        Box::from_raw(core::ptr::slice_from_raw_parts_mut(p, n))
    }
}

// =============================================================================================
// public parameters
// =============================================================================================

/// Words of the pair key per lockstep step: one dword per pair batch per slot and stream.
pub const LOCK_WORDS: usize = N * STREAMS * LOCK;

/// The public data the commitment reads: the keys of all streams in the lockstep layout, the
/// digit chunks' centered mean residues, and the total of the means.
pub struct SmsKeys<'a> {
    /// `pairs[(l * N + j) * STREAMS * LOCK + LOCK s + h]`: for lockstep step `l` (ring elements
    /// `2 LOCK l .. 2 LOCK (l + 1)` of every column), slot `j` and stream `s`, the dword of the
    /// stream's NTT-domain keys of ring elements `2 LOCK l + 2h` (low half) and `+ 1` (high half),
    /// centered. Auxiliary streams include the public per-Phi_6-half scale
    /// consumed by `intt_sms_scaled_batch32_9721`; the base stream is unscaled.
    pub pairs: &'a [u32],
    /// Paired means: `mur[2*((d/2*chunks+c)*N+m)+d%2]` is the centered mean
    /// residue of digit `d`, chunk `c`, coefficient `m`. Adjacent words hold
    /// the two digits consumed by one inverse batch and one dword broadcast.
    pub mur: &'a [i16],
    /// `mu[m] = sum_{d,c} digit_weight(d) * mu'_{d,c}[m]`, wrapping.
    pub mu: &'a [u32],
    /// Ring elements per column.
    pub elements: usize,
}

/// Owned [`SmsKeys`] together with the statistical margin they were certified at.
pub struct SmsKeyData {
    pub pairs: Vec<u32>,
    pub mur: Vec<i16>,
    pub mu: Vec<u32>,
    pub elements: usize,
    /// The smallest `((P - 1) / 2) / sigma` over every chunk, digit and coefficient.
    pub margin: f64,
    /// The same per digit, and where each digit's minimum is (chunk, coefficient, |mean|).
    pub margins: [(f64, usize, usize, i64); DIGITS],
}

impl SmsKeyData {
    pub fn keys(&self) -> SmsKeys<'_> {
        SmsKeys { pairs: &self.pairs, mur: &self.mur, mu: &self.mu, elements: self.elements }
    }
}

/// Canonical scalar digits, retained as a regression reference for the radix identity.
/// Public matrix construction uses the paired encoding below.
pub fn digits(x: u128) -> [i8; DIGITS] {
    let mut r = (x & TAIL_MASK as u128) as i64;
    if r >= (1 << TAIL_BITS) - REP_LOW {
        r -= 1 << TAIL_BITS;
    }
    let mut out = [0i8; DIGITS];
    for d in 0..DIGITS {
        let b = digit_radix(d);
        let v = (r + DIGIT_HALF).rem_euclid(b) - DIGIT_HALF;
        out[d] = v as i8;
        r = (r - v) / b;
    }
    debug_assert_eq!(r, 0);
    out
}

/// Squared norm in the quadratic pair basis: Z^2 - Z + 1 = 0, Z = X^324.
const fn pair_norm(a: i64, b: i64) -> i64 { a * a + a * b + b * b }

#[derive(Clone, Copy)]
struct PairChoices { values: [[i8; 2]; 9], count: usize }

const fn pair_digit_table() -> [PairChoices; PAIR_RESIDUES] {
    let mut table = [PairChoices { values: [[0; 2]; 9], count: 0 }; PAIR_RESIDUES];
    let mut i = 0;
    while i < PAIR_RESIDUES {
        let a = (i / DIGIT_RADIX as usize) as i64 - DIGIT_HALF;
        let b = (i % DIGIT_RADIX as usize) as i64 - DIGIT_HALF;
        let mut best = i64::MAX;
        let mut u = -1;
        while u <= 1 {
            let mut v = -1;
            while v <= 1 {
                let (p, q) = (a + DIGIT_RADIX * u, b + DIGIT_RADIX * v);
                let cost = pair_norm(p, q);
                if cost < best { best = cost; table[i].count = 0; }
                if cost == best {
                    let k = table[i].count;
                    table[i].values[k] = [p as i8, q as i8];
                    table[i].count += 1;
                }
                v += 1;
            }
            u += 1;
        }
        let mut k = 0;
        while k < table[i].count {
            assert!(table[i].values[k][0] >= -DIGIT_BOUND && table[i].values[k][0] <= DIGIT_BOUND);
            assert!(table[i].values[k][1] >= -DIGIT_BOUND && table[i].values[k][1] <= DIGIT_BOUND);
            k += 1;
        }
        i += 1;
    }
    table
}
const PAIR_DIGITS: [PairChoices; PAIR_RESIDUES] = pair_digit_table();

/// Choose a short representative of a coefficient pair modulo 2^22. Only public
/// key coefficients and their position select among ties; no instance data enters.
fn pair_representative(x: u128, y: u128, salt: usize) -> (i64, i64) {
    let modulus = 1i64 << TAIL_BITS;
    let center = |x: u128| {
        let r = (x & TAIL_MASK as u128) as i64;
        if r >= modulus / 2 { r - modulus } else { r }
    };
    let (a, b) = (center(x), center(y));
    let mut choices = [(0i64, 0i64); 9];
    let mut count = 0;
    let mut best = i64::MAX;
    for u in -1..=1 {
        for v in -1..=1 {
            let (p, q) = (a + u * modulus, b + v * modulus);
            let cost = pair_norm(p, q);
            if cost < best { best = cost; count = 0; }
            if cost == best { choices[count] = (p, q); count += 1; }
        }
    }
    choices[salt % count]
}

/// Encode a selected integer pair with minimum-norm radix residues and exact carries.
/// The final carry is retained directly, preserving the two integer identities.
fn encode_pair(mut a: i64, mut b: i64, salt: usize) -> ([i8; DIGITS], [i8; DIGITS]) {
    let (mut da, mut db) = ([0i8; DIGITS], [0i8; DIGITS]);
    for d in 0..DIGITS - 1 {
        let i = (a + DIGIT_HALF).rem_euclid(DIGIT_RADIX) as usize * DIGIT_RADIX as usize
            + (b + DIGIT_HALF).rem_euclid(DIGIT_RADIX) as usize;
        let choices = PAIR_DIGITS[i];
        let [p, q] = choices.values[salt.wrapping_add(d) % choices.count];
        da[d] = p;
        db[d] = q;
        a = (a - p as i64) / DIGIT_RADIX;
        b = (b - q as i64) / DIGIT_RADIX;
    }
    debug_assert!((-(TOP_DIGIT_BOUND as i64)..=TOP_DIGIT_BOUND as i64).contains(&a)
        && (-(TOP_DIGIT_BOUND as i64)..=TOP_DIGIT_BOUND as i64).contains(&b));
    da[DIGITS - 1] = a as i8;
    db[DIGITS - 1] = b as i8;
    (da, db)
}

/// Public-key representative search: redistribute the digit energy while
/// keeping the represented integer pair no longer than ordinary centering.
/// The fixed weights balance the tight lower-digit reconstruction margins.
pub fn digits_pair(x: u128, y: u128, salt: usize) -> ([i8; DIGITS], [i8; DIGITS]) {
    const WEIGHTS: [i64; DIGITS] = [1, 1, 3, 2, 2, 1];
    let score = |digits: &([i8; DIGITS], [i8; DIGITS])| {
        (0..DIGITS).map(|d| WEIGHTS[d] * pair_norm(digits.0[d] as i64, digits.1[d] as i64)).sum::<i64>()
    };
    let modulus = 1i64 << TAIL_BITS;
    let center = |v: u128| {
        let r = (v & TAIL_MASK as u128) as i64;
        if r >= modulus / 2 { r - modulus } else { r }
    };
    let (a, b) = (center(x), center(y));
    let limit = pair_norm(a, b);
    let (p, q) = pair_representative(x, y, salt);
    let mut best = encode_pair(p, q, salt);
    let mut cost = score(&best);
    for u in -1..=1 {
        for v in -1..=1 {
            let (p, q) = (a + u * modulus, b + v * modulus);
            if pair_norm(p, q) > limit { continue; }
            let candidate = encode_pair(p, q, salt);
            let next = score(&candidate);
            if next < cost { cost = next; best = candidate; }
        }
    }
    best
}

/// `sum_i M(L)[m][i]^k` over the 512 witness coefficients `i` for `k = 1` and `k = 2`, where
/// `M(L)[m][i]` is coefficient `m` of `X^i L mod X^648 - X^324 + 1`. For `m >= 324` it is
/// `L[m - i] + L[m + 324 - i]`, for `m < 324` it is `L[m - i] - L[m + 648 - i] - L[m + 972 - i]`
/// (`L` zero outside `0..648`): each a function of `x = m - i` alone, so both sums are windows of
/// a prefix sum.
fn moment_rows(l: &[i8; N], s1: &mut [i64; N], s2: &mut [i64; N]) {
    let at = |x: i64| -> i64 { if (0..N as i64).contains(&x) { l[x as usize] as i64 } else { 0 } };
    // x runs over -511 ..= 647; index x + 511.
    const W: usize = 512 + N - 1;
    let mut fh = [0i64; W + 1];
    let mut fh2 = [0i64; W + 1];
    let mut fl = [0i64; W + 1];
    let mut fl2 = [0i64; W + 1];
    for k in 0..W {
        let x = k as i64 - 511;
        let h = at(x) + at(x + 324);
        let g = at(x) - at(x + 648) - at(x + 972);
        fh[k + 1] = fh[k] + h;
        fh2[k + 1] = fh2[k] + h * h;
        fl[k + 1] = fl[k] + g;
        fl2[k + 1] = fl2[k] + g * g;
    }
    for m in 0..N {
        // x from m - 511 to m: indices m .. m + 512
        let (a, b) = (m, m + 512);
        if m >= 324 {
            s1[m] += fh[b] - fh[a];
            s2[m] += fh2[b] - fh2[a];
        } else {
            s1[m] += fl[b] - fl[a];
            s2[m] += fl2[b] - fl2[a];
        }
    }
}

impl SmsKeyData {
    /// The keys of the matrix `a` (`elements * 648` coefficients modulo [`MODULUS`]).
    pub fn build(a: &[u128], elements: usize) -> SmsKeyData {
        use crate::key::CommitmentKey;
        use crate::ring::Modulus;
        assert_eq!(a.len(), elements * N);
        assert!((elements / 4) % LOCK == 0 && ((elements / 4) % CHUNK_PAIRS) % LOCK == 0, "whole lockstep steps per chunk");
        let base = Modulus::ALL.into_iter().find(|m| m.prime() == P).unwrap();
        let len_f162 = 4 * elements;
        let mut dig = vec![[0i8; DIGITS]; a.len()];
        for (ring_index, (coefficients, out)) in
            a.chunks_exact(N).zip(dig.chunks_exact_mut(N)).enumerate()
        {
            for j in 0..N / 2 {
                let (lo, hi) = digits_pair(coefficients[j], coefficients[j + N / 2], ring_index * (N / 2) + j);
                out[j] = lo;
                out[j + N / 2] = hi;
            }
        }
        let mut keys = Vec::with_capacity(STREAMS);
        keys.push(CommitmentKey::from_coefficients(len_f162, a, MODULUS, base, &[]));
        for d in 0..DIGITS {
            let c: Vec<u128> = dig.iter().map(|v| (v[d] as i64).rem_euclid(P as i64) as u128).collect();
            keys.push(CommitmentKey::from_coefficients(len_f162, &c, P as u128, base, &[]));
        }
        let steps = elements / (2 * LOCK);
        let mut pairs = vec![0u32; steps * LOCK_WORDS];
        for l in 0..steps {
            let j0 = 2 * LOCK * l;
            let (b, lane) = (j0 / 32, j0 % 32);
            for (s, key) in keys.iter().enumerate() {
                let batch = &key.row(0)[b];
                for j in 0..N {
                    for h in 0..LOCK {
                        // Only auxiliary streams use the preconditioned inverse.
                        // The base stream retains its original commitment encoding.
                        let scale = gen_small::sms_key_scales_9721()[usize::from(j >= N / 2)];
                        let scaled = |v: i16| -> u32 {
                            if s == 0 { return v as u16 as u32; }
                            let x = (v as i64 * scale as i64).rem_euclid(P as i64);
                            let centered = if x > (P / 2) as i64 { x - P as i64 } else { x };
                            centered as i16 as u16 as u32
                        };
                        let lo = scaled(batch.v[j][lane + 2 * h]);
                        let hi = scaled(batch.v[j][lane + 2 * h + 1]);
                        pairs[(l * N + j) * STREAMS * LOCK + LOCK * s + h] = lo | (hi << 16);
                    }
                }
            }
        }
        // Chunk moments: the ring elements of chunk `c` ([`chunk_range`]) of every column.
        let chunks = chunk_count(elements);
        let mut mur = vec![0i16; mur_len(elements)];
        let mut mu = vec![0u32; N];
        let mut margin = f64::INFINITY;
        let mut margins = [(f64::INFINITY, 0usize, 0usize, 0i64); DIGITS];
        let half = ((P - 1) / 2) as f64;
        for d in 0..DIGITS {
            let nc = chunks;
            for c in 0..nc {
                let mut s1 = [0i64; N];
                let mut s2 = [0i64; N];
                let ring = {
                    let r = chunk_range(elements / 2, c);
                    2 * r.start..2 * r.end
                };
                for j in ring {
                    let l: [i8; N] = core::array::from_fn(|m| dig[j * N + m][d]);
                    moment_rows(&l, &mut s1, &mut s2);
                }
                for m in 0..N {
                    // mean s1 / 2, variance s2 / 4
                    let mean = s1[m].div_euclid(2);
                    let sigma = (s2[m] as f64 / 4.0).sqrt();
                    if sigma > 0.0 {
                        let z = (half - 0.5) / sigma;
                        margin = margin.min(z);
                        if z < margins[d].0 {
                            margins[d] = (z, c, m, mean.abs());
                        }
                    }
                    let r = mean.rem_euclid(P as i64);
                    mur[2 * (((d / 2) * chunks + c) * N + m) + d % 2] = if r > (P as i64 - 1) / 2 { (r - P as i64) as i16 } else { r as i16 };
                    mu[m] = mu[m].wrapping_add((mean as u32).wrapping_mul(digit_weight(d)));
                }
            }
        }
        SmsKeyData { pairs, mur, mu, elements, margin, margins }
    }
}

// =============================================================================================
// the kernel's sink
// =============================================================================================

/// `v[j][s]`: slot `j`, stream `s`, with sixteen column accumulators.
#[repr(C, align(64))]
pub struct SAcc {
    pub v: [[[i32; 16]; STREAMS]; N],
}

#[repr(C, align(64))]
struct Blk27([i16; 27 * 32]);

/// Fold-back constants as memory operands: the `vpmaddwd` weights `(1, r16c)`, the low-word sign
/// bit, and the old masked correction (unused by the two-uop fold).
///
/// The fold is `vpxord` with `0x8000` and `vpmaddwd` with `(1, r16c)`: flipping bit 15 makes the
/// low word, read as signed, `l - 32768` for the unsigned low word `l`, so the result is
/// `x - 32768 (mod q)` with the same bound as [`reduce_vec_c`] (`|r| <= 32768 (1 + |r16c|)`).
/// Every fold therefore subtracts 32768 from every lane; [`fold_bias`] pre-loads each stream with
/// the total of the folds it will take before its next reset, so the residues come out exact.
#[repr(C, align(64))]
struct Fold([i32; 48]);
static FOLD: Fold = Fold({
    let mut f = [0i32; 48];
    let mut i = 0;
    while i < 16 {
        f[i] = 1 | ((r16c(P) as u16 as i32) << 16);
        f[16 + i] = 0x8000;
        f[32 + i] = r16c(P);
        i += 1;
    }
    f
});

// One slot's stream accumulators, round by round, with plain and folded variants.
include!("sms_slot_gen.rs");
const _: () = assert!(STREAMS == 7 && LOCK == 8);

/// Whether a group step prefetches the next step's witness words into L2.
const INPUT_PREFETCH: bool = false;

/// Slots ahead the sink prefetches its key lines. The accumulator lines are L1-resident already,
/// so prefetching them only spends issue slots.
const PREFETCH_SLOTS: usize = 3;

/// Pair batches between two folds of the sink: every [`PP`] pair batches, at the end of a call
/// (every `PP / LOCK` calls) or, for `LOCK > PP`, inside it as well (the slot macro folds after
/// every eight lockstep products).
const FOLD_EVERY: usize = PP;
const _: () = assert!(LOCK <= PP || PP == 8);

/// The value stream `s` starts from so that its two-uop folds cancel by its next reset: `32768`
/// per fold, modulo `P`. The base runs over the commitment range; each digit
/// runs over its reconstruction chunk. `span` counts that range's pair batches.
fn fold_bias(span: usize) -> i32 {
    debug_assert!(span % FOLD_EVERY == 0);
    ((span / FOLD_EVERY) as u64 * 32768 % P as u64) as i32
}

/// Every slot of stream `s` set to `v`.
#[target_feature(enable = "avx512f")]
unsafe fn fill_stream(a: &mut SAcc, s: usize, v: i32) {
    let p = (a.v.as_mut_ptr() as *mut __m512i).add(s);
    let x = _mm512_set1_epi32(v);
    for j in 0..N {
        _mm512_store_si512(p.add(STREAMS * j), x);
    }
}

/// The 4-batch lockstep sink: the batches' 27-slot blocks are staged in L1 and every slot's nine
/// accumulators take all four products in one read-modify-write, folded back when `RED`.
///
/// `FIRST` marks a chunk's first lockstep step: the digit streams then start from their fold bias
/// broadcast out of `b` instead of being loaded from memory, because [`finish_stream`] left that
/// value there and no longer writes it into all 648 accumulator rows.
struct Sink<const RED: bool, const FIRST: bool = false> {
    k: *const u32,
    acc: *mut i32,
    b: *const i32,
}
impl<const RED: bool, const FIRST: bool> vb::SliceSink<LOCK> for Sink<RED, FIRST> {
    #[inline(always)]
    unsafe fn rows(&mut self, blk: usize, w: [*const i16; LOCK], r0: usize, r1: usize) {
        sink_rows::<RED, FIRST>(self.k, self.acc, self.b, blk, w, r0, r1);
    }
}

/// [`Sink::rows`], one copy for every call site: the transform's hot loop and this body together
/// then fit the half of the decoded-uop cache a hyperthread gets.
#[inline(never)]
#[target_feature(enable = "avx512f,avx512bw,avx512vnni")]
unsafe fn sink_rows<const RED: bool, const FIRST: bool>(kp: *const u32, accp: *mut i32, bias: *const i32, blk: usize, w: [*const i16; LOCK], r0: usize, r1: usize) {
    let s = Sink::<RED, FIRST> { k: kp, acc: accp, b: bias };
    s.rows_inline(blk, w, r0, r1);
}

impl<const RED: bool, const FIRST: bool> Sink<RED, FIRST> {
    #[inline(always)]
    unsafe fn rows_inline(&self, blk: usize, w: [*const i16; LOCK], r0: usize, r1: usize) {
        for i in r0..r1 {
            let j = 27 * blk + i;
            let x: [__m512i; LOCK] = core::array::from_fn(|h| _mm512_load_si512((w[h] as *const __m512i).add(i)));
            let k = self.k.add(j * STREAMS * LOCK);
            let a = self.acc.add(j * STREAMS * 16);
            if PREFETCH_SLOTS > 0 && j + PREFETCH_SLOTS < N {
                let kf = k.add(PREFETCH_SLOTS * STREAMS * LOCK) as *const i8;
                for l in 0..(STREAMS * LOCK * 4).div_ceil(64) {
                    _mm_prefetch::<_MM_HINT_T0>(kf.add(64 * l));
                }
            }
            macro_rules! call {
                ($m:ident, $mode:ident, $($xn:ident = $xi:expr),*) => {
                    core::arch::asm!(concat!($m!($mode), "/* {f} */"), a = in(reg) a, k = in(reg) k, f = in(reg) FOLD.0.as_ptr(),
                        $($xn = in(zmm_reg) x[$xi],)*
                        out("zmm20") _, out("zmm21") _, out("zmm22") _, out("zmm23") _, out("zmm24") _, out("zmm25") _,
                        out("zmm26") _, out("zmm27") _, out("zmm28") _, out("zmm29") _, out("k1") _, options(nostack, preserves_flags))
                };
            }
            macro_rules! call_first {
                ($m:ident, $mode:ident, $($xn:ident = $xi:expr),*) => {
                    core::arch::asm!(concat!($m!($mode), "/* {f} */"), a = in(reg) a, k = in(reg) k, f = in(reg) FOLD.0.as_ptr(),
                        b = in(reg) self.b,
                        $($xn = in(zmm_reg) x[$xi],)*
                        out("zmm20") _, out("zmm21") _, out("zmm22") _, out("zmm23") _, out("zmm24") _, out("zmm25") _,
                        out("zmm26") _, out("zmm27") _, out("zmm28") _, out("zmm29") _, out("k1") _, options(nostack, preserves_flags))
                };
            }
            if LOCK == 8 {
                if RED && FIRST { call_first!(sms_slot8, reduce_first, x0 = 0, x1 = 1, x2 = 2, x3 = 3, x4 = 4 % LOCK, x5 = 5 % LOCK, x6 = 6 % LOCK, x7 = 7 % LOCK) }
                else if FIRST { call_first!(sms_slot8, plain_first, x0 = 0, x1 = 1, x2 = 2, x3 = 3, x4 = 4 % LOCK, x5 = 5 % LOCK, x6 = 6 % LOCK, x7 = 7 % LOCK) }
                else if RED { call!(sms_slot8, reduce, x0 = 0, x1 = 1, x2 = 2, x3 = 3, x4 = 4 % LOCK, x5 = 5 % LOCK, x6 = 6 % LOCK, x7 = 7 % LOCK) }
                else { call!(sms_slot8, plain, x0 = 0, x1 = 1, x2 = 2, x3 = 3, x4 = 4 % LOCK, x5 = 5 % LOCK, x6 = 6 % LOCK, x7 = 7 % LOCK) }
            } else if LOCK == 16 {
                call!(sms_slot16, reduce, x0 = 0, x1 = 1, x2 = 2, x3 = 3, x4 = 4 % LOCK, x5 = 5 % LOCK, x6 = 6 % LOCK, x7 = 7 % LOCK,
                    x8 = 8 % LOCK, x9 = 9 % LOCK, x10 = 10 % LOCK, x11 = 11 % LOCK, x12 = 12 % LOCK, x13 = 13 % LOCK, x14 = 14 % LOCK, x15 = 15 % LOCK)
            } else {
                if RED { call!(sms_slot4, reduce, x0 = 0, x1 = 1, x2 = 2, x3 = 3) }
                else { call!(sms_slot4, plain, x0 = 0, x1 = 1, x2 = 2, x3 = 3) }
            }
        }
    }
}

const PAIR_PACK_PERM: [i64; 8] = [0, 2, 4, 6, 1, 3, 5, 7];

/// Stream `s` of the accumulators modulo 9721 as `i16` representatives (`|r| < 2q`) into lanes
/// `16 half .. 16 half + 16` of `out`. The stream's next chunk starts from the worker's bias
/// instead ([`Sink`] with `FIRST`), so nothing is written back into the accumulator here.
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn finish_stream(a: &mut SAcc, s: usize, out: &mut Batch32, half: usize) {
    let p = (a.v.as_mut_ptr() as *mut __m512i).add(s);
    let perm = _mm512_loadu_si512(PAIR_PACK_PERM.as_ptr() as *const __m512i);
    let l0 = 16 * half;
    for k in 0..N / 2 {
        let j = 2 * k;
        let x0 = _mm512_load_si512(p.add(STREAMS * j));
        let x1 = _mm512_load_si512(p.add(STREAMS * (j + 1)));
        let packed = _mm512_packs_epi32(red_folded(x0), red_folded(x1));
        let u = _mm512_permutexvar_epi64(perm, packed);
        _mm256_storeu_si256(out.v[j].as_mut_ptr().add(l0) as *mut __m256i, _mm512_castsi512_si256(u));
        _mm256_storeu_si256(out.v[j + 1].as_mut_ptr().add(l0) as *mut __m256i, _mm512_extracti64x4_epi64::<1>(u));
    }
}

/// One sixteen-column total; each inverse pair is weighted and combined before this store.
#[repr(C, align(64))]
struct Total {
    t: [[i32; 16]; N],
}

/// Digits `d` and `d+1` of one chunk share the two 16-column lane banks.
/// The paired mean words are broadcast together after interleaving the lanes.
/// Each bank uses its own radix weight before both contribute to one total.
#[target_feature(enable = "avx512f,avx512bw,avx512vl")]
unsafe fn reconstruct(b: &Batch32, d: usize, keys: &SmsKeys, c0: usize, total: &mut Total) {
    let chunks = chunk_count(keys.elements);
    let means = keys.mur.as_ptr().add(2 * (((d / 2) * chunks + c0) * N)).cast::<i32>();
    let q = _mm512_set1_epi16(P as i16);
    let hq = _mm512_set1_epi16(((P - 1) / 2) as i16);
    let nhq = _mm512_set1_epi16(-(((P - 1) / 2) as i16));
    let narrow = digit_weight(d + 1) <= i16::MAX as u32;
    let weights = if narrow { digit_weight(d) | (digit_weight(d + 1) << 16) }
        else { 1 | ((DIGIT_RADIX as u32) << 16) };
    let weights = _mm512_set1_epi32(weights as i32);
    let wd = _mm512_set1_epi32(digit_weight(d) as i32);
    let interleave: [i16; 32] = core::array::from_fn(|i| (i / 2 + 16 * (i % 2)) as i16);
    let interleave = _mm512_loadu_si512(interleave.as_ptr().cast());
    for m in 0..N {
        let x = _mm512_load_si512(b.v[m].as_ptr() as *const __m512i);
        let x = _mm512_permutexvar_epi16(interleave, x);
        // The public image stores (low mean, high mean) as one little-endian
        // dword. Unaligned reads also support caller-owned i16 slices.
        let mu = _mm512_set1_epi32(core::ptr::read_unaligned(means.add(m)));
        let mut y = _mm512_sub_epi16(x, mu);
        y = _mm512_mask_sub_epi16(y, _mm512_cmpgt_epi16_mask(y, hq), y, q);
        y = _mm512_mask_add_epi16(y, _mm512_cmplt_epi16_mask(y, nhq), y, q);
        // Adjacent words are the two digits of the same column. Both full
        // weights fit signed words for pairs 0/1 and 2/3; pair 4/5 first uses
        // weights (1,13), then scales its dword result by 13^4.
        let mut v = _mm512_madd_epi16(y, weights);
        if !narrow { v = _mm512_mullo_epi32(v, wd); }
        let t = total.t[m].as_mut_ptr() as *mut __m512i;
        _mm512_store_si512(t, _mm512_add_epi32(_mm512_load_si512(t), v));
    }
}

/// Per-worker scratch, allocated once and reused by every group the worker takes.
struct SmsWorker {
    acc: Box<[SAcc]>,
    pend: Box<[Batch32]>,
    total: Box<[Total]>,
    idx: Box<[BinaryIndex32]>,
    w: Box<[Blk27]>,
    /// The digit streams' start value for the next chunk, broadcast by the first lockstep step.
    bias: Box<[i32; 1]>,
}
impl SmsWorker {
    fn new() -> Self {
        SmsWorker {
            acc: zeroed_slice(1),
            pend: zeroed_slice(1),
            total: zeroed_slice(1),
            idx: zeroed_slice(LOCK),
            w: zeroed_slice(2 * LOCK),
            bias: Box::new([0i32; 1]),
        }
    }
}

/// 16 x 16 u32 transpose: row `s` holds columns `0..16` of slot `j0 + s`; output `c` holds slots
/// `j0 .. j0 + 16` of column `c`.
#[inline(always)]
unsafe fn transpose16x16_u32(r: &[__m512i; 16]) -> [__m512i; 16] {
    let mut u = [[_mm512_setzero_si512(); 4]; 4];
    for q in 0..4 {
        let t0 = _mm512_unpacklo_epi32(r[4 * q], r[4 * q + 1]);
        let t1 = _mm512_unpackhi_epi32(r[4 * q], r[4 * q + 1]);
        let t2 = _mm512_unpacklo_epi32(r[4 * q + 2], r[4 * q + 3]);
        let t3 = _mm512_unpackhi_epi32(r[4 * q + 2], r[4 * q + 3]);
        u[q] = [
            _mm512_unpacklo_epi64(t0, t2),
            _mm512_unpackhi_epi64(t0, t2),
            _mm512_unpacklo_epi64(t1, t3),
            _mm512_unpackhi_epi64(t1, t3),
        ];
    }
    let mut out = [_mm512_setzero_si512(); 16];
    for m in 0..4 {
        let t0 = _mm512_shuffle_i64x2::<0x88>(u[0][m], u[1][m]);
        let t1 = _mm512_shuffle_i64x2::<0xDD>(u[0][m], u[1][m]);
        let t2 = _mm512_shuffle_i64x2::<0x88>(u[2][m], u[3][m]);
        let t3 = _mm512_shuffle_i64x2::<0xDD>(u[2][m], u[3][m]);
        let v = [
            _mm512_shuffle_i64x2::<0x88>(t0, t2),
            _mm512_shuffle_i64x2::<0x88>(t1, t3),
            _mm512_shuffle_i64x2::<0xDD>(t0, t2),
            _mm512_shuffle_i64x2::<0xDD>(t1, t3),
        ];
        for l in 0..4 {
            out[4 * l + m] = v[l];
        }
    }
    out
}

/// What [`write_rows`] writes: a whole column's rows, or one of the two halves of its
/// ring-element range ([`HALVES`]). The first half to finish leaves its raw partial rows (the
/// base residue canonical in `[0, 9721)`, the tail total without the means and unmasked); the
/// second adds its own to them and writes the final rows.
#[derive(Clone, Copy, PartialEq, Eq)]
enum RowsMode {
    Full,
    First,
    Second,
}

/// The group's two rows per column: the base stream canonical in `[0, 9721)` (NTT domain) and the
/// tail `t mod 2^22` (coefficients); accumulators and total cleared for the next group.
#[target_feature(enable = "avx512f,avx512bw,avx512vl")]
unsafe fn write_rows(acc: &mut SAcc, total: &mut Total, keys: &SmsKeys, gc: usize, dst: *mut [u32; N], mode: RowsMode) {
    let zero = _mm512_setzero_si512();
    let mask = _mm512_set1_epi32(TAIL_MASK as i32);
    let q = _mm512_set1_epi32(P as i32);
    let pa = acc.v.as_mut_ptr() as *mut __m512i;
    let pt = total.t.as_mut_ptr() as *mut __m512i;
    let mut j0 = 0;
    while j0 < N {
        let n = 16.min(N - j0);
        let mut r0 = [zero; 16];
        let mut rt = [zero; 16];
        for s in 0..n {
            let j = j0 + s;
            let v0 = _mm512_load_si512(pa.add(STREAMS * j));
            _mm512_store_si512(pa.add(STREAMS * j), zero);
            r0[s] = canon_red16::<P>(v0);
            let a = _mm512_load_si512(pt.add(j));
            _mm512_store_si512(pt.add(j), zero);
            let mu = if mode == RowsMode::First { 0 } else { keys.mu[j] as i32 };
            rt[s] = _mm512_add_epi32(a, _mm512_set1_epi32(mu));
        }
        let c0 = transpose16x16_u32(&r0);
        let ct = transpose16x16_u32(&rt);
        let m = if n == 16 { 0xFFFF } else { (1u16 << n) - 1 };
        for c in 0..gc {
            let (db, dt) = ((dst.add(2 * c) as *mut i32).add(j0), (dst.add(2 * c + 1) as *mut i32).add(j0));
            let (mut b, mut t) = (c0[c], ct[c]);
            if mode == RowsMode::Second {
                b = _mm512_add_epi32(b, _mm512_maskz_loadu_epi32(m, db));
                b = _mm512_mask_sub_epi32(b, _mm512_cmp_epi32_mask(b, q, _MM_CMPINT_NLT), b, q);
                t = _mm512_add_epi32(t, _mm512_maskz_loadu_epi32(m, dt));
            }
            if mode != RowsMode::First {
                t = _mm512_and_si512(t, mask);
            }
            _mm512_mask_storeu_epi32(db, m, b);
            _mm512_mask_storeu_epi32(dt, m, t);
        }
        j0 += 16;
    }
}

/// Parts of a column's ring-element range committed as separate work units by
/// [`commit_paired`]: the commitment's coarse units otherwise leave cores idle at its end. Each
/// part is whole digit chunks and one half-column chunk of the last digit.
const HALVES: usize = 1;

/// Where a group writes its rows: in full, or as one half of the range under its group's lock
/// (`state` counts the halves already written).
enum RowsOut<'a> {
    Full,
    Half(&'a std::sync::Mutex<u8>),
}

/// Commits a group of up to 16 columns across all `np` pair batches and writes its `2 g` rows to
/// `dst`. `sync(b_done)` is called at every [`PP`] boundary so hyperthread siblings can share the
/// key stream in L2.
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx512vbmi,avx512vbmi2,avx512vnni,gfni")]
unsafe fn group<F: FnMut(usize)>(cols: &[&[B128]], keys: &SmsKeys, np: usize, range: core::ops::Range<usize>, out: RowsOut, ws: &mut SmsWorker, dst: *mut [u32; N], mut sync: F) {
    let g = cols.len();
    assert!(g >= 1 && g <= 16 && (np / 2) % LOCK == 0 && ((np / 2) % CHUNK_PAIRS) % LOCK == 0 && 2 * np == keys.elements);
    let hl = np / 2;
    let span = range.end - range.start;
    assert!(range.start % (np / 2) == 0 && (span == np || span == np / 2), "whole halves of the column");
    static ZERO8: [B128; 8] = [B128::ZERO; 8];
    let kbase = keys.pairs.as_ptr();
    for s in 0..STREAMS {
        // Stream 0 runs over the range; the digits' spans are whole chunks of it.
        let b = if s == 0 { fold_bias(span) } else { fold_bias(CHUNK_PAIRS.min(hl)) };
        fill_stream(&mut ws.acc[0], s, b);
    }
    // The range's first chunk starts from the same value; every later chunk end rewrites it.
    ws.bias[0] = fold_bias(CHUNK_PAIRS.min(hl));
    let acc = ws.acc[0].v.as_mut_ptr() as *mut i32;
    let w = ws.w.as_mut_ptr() as *mut i16;
    for b in range.clone().step_by(LOCK) {
        // The next step's witness words (16 columns x 1 KB, 256 KB apart): the hardware
        // prefetchers only follow some of the sixteen streams, and the table comes from DRAM.
        if INPUT_PREFETCH && b + LOCK < range.end {
            for l in 0..g {
                let p = cols[l].as_ptr().add(8 * (b + LOCK)) as *const i8;
                for line in 0..(8 * LOCK * 16) / 64 {
                    _mm_prefetch::<_MM_HINT_T1>(p.add(64 * line));
                }
            }
        }
        for h in 0..LOCK {
            let cols16: [*const B128; 16] = core::array::from_fn(|l| {
                if l < g { cols[l].as_ptr().add(8 * (b + h)) } else { ZERO8.as_ptr() }
            });
            crate::simd::transpose_f162::slice_b128_cols16_into(&cols16, &mut ws.idx[h]);
        }
        let ins: [&BinaryIndex32; LOCK] = core::array::from_fn(|h| &ws.idx[h]);
        let k = kbase.add((b / LOCK) * LOCK_WORDS);
        let bias = ws.bias.as_ptr();
        // A chunk's first step: the digit streams start from `ws.bias`, the value the last chunk
        // end left there, instead of loading 648 accumulator rows each.
        let first = (b - b / hl * hl) % CHUNK_PAIRS == 0;
        let mode = bench_mode();
        if mode == 4 { continue; }
        if mode == 1 {
            struct Nop;
            impl vb::SliceSink<LOCK> for Nop { unsafe fn rows(&mut self, _: usize, w: [*const i16; LOCK], _: usize, _: usize) { core::hint::black_box(w); } }
            vb::ntt_bin_batch32_pipelined::<P, LOCK, _>(ins, w, &mut Nop);
        } else if mode == 2 {
            let mut sk = Sink::<false> { k, acc, b: bias };
            for kk in 0..24 { vb::SliceSink::<LOCK>::rows(&mut sk, kk, core::array::from_fn(|h| w.add(32 * 27 * h) as *const i16), 0, 27); }
        } else if first && ((b + LOCK) % PP == 0 || LOCK > PP) {
            vb::ntt_bin_batch32_pipelined::<P, LOCK, _>(ins, w, &mut Sink::<true, true> { k, acc, b: bias });
            sync(b + LOCK);
        } else if first {
            vb::ntt_bin_batch32_pipelined::<P, LOCK, _>(ins, w, &mut Sink::<false, true> { k, acc, b: bias });
        } else if (b + LOCK) % PP == 0 || LOCK > PP {
            vb::ntt_bin_batch32_pipelined::<P, LOCK, _>(ins, w, &mut Sink::<true> { k, acc, b: bias });
            sync(b + LOCK);
        } else {
            vb::ntt_bin_batch32_pipelined::<P, LOCK, _>(ins, w, &mut Sink::<false> { k, acc, b: bias });
        }
        // Chunk ends: every CHUNK_PAIRS pair batches into a half, and the half's end.
        let (e, cph) = (b + LOCK, chunks_per_half(np));
        let (eh, local) = ((e - 1) / hl, e - (e - 1) / hl * hl);
        if (local % CHUNK_PAIRS == 0 || local == hl) && mode != 3 {
            let chunk = eh * cph + (local - 1) / CHUNK_PAIRS;
            // The next chunk's length: the rest of this half, or the next half's first chunk.
            let next = if local == hl { CHUNK_PAIRS.min(hl) } else { CHUNK_PAIRS.min(hl - local) };
            ws.bias[0] = fold_bias(next);
            for d in (0..CHUNKED).step_by(2) {
                let pd = &mut ws.pend[0];
                finish_stream(&mut ws.acc[0], 1 + d, pd, 0);
                finish_stream(&mut ws.acc[0], 2 + d, pd, 1);
                pd.representation = Representation::Ntt;
                gen_small::intt_sms_scaled_batch32_9721(pd);
                reconstruct(pd, d, keys, chunk, &mut ws.total[0]);
            }
        }
    }
    match out {
        RowsOut::Full => write_rows(&mut ws.acc[0], &mut ws.total[0], keys, g, dst, RowsMode::Full),
        RowsOut::Half(state) => {
            let mut done = state.lock().unwrap();
            let mode = if *done == 0 { RowsMode::First } else { RowsMode::Second };
            write_rows(&mut ws.acc[0], &mut ws.total[0], keys, g, dst, mode);
            *done += 1;
        }
    }
}

/// Benchmark-only phase isolation (`SMS_MODE`): 0 everything, 1 no accumulation, 2 accumulation
/// only, 3 no chunk reconstruction, 4 slicing only. Results are meaningless unless 0.
#[inline(always)]
fn bench_mode() -> u32 {
    #[cfg(test)]
    {
        static M: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
        *M.get_or_init(|| std::env::var("SMS_MODE").ok().and_then(|v| v.parse().ok()).unwrap_or(0))
    }
    #[cfg(not(test))]
    {
        0
    }
}

/// A batch counter on a cache line of its own.
#[repr(align(128))]
struct Progress(std::sync::atomic::AtomicUsize);
const SIBLING_WAIT: std::time::Duration = std::time::Duration::from_micros(500);
/// Maximum lead (in pair batches) between two hyperthread siblings on one physical core.
const PAIR_LEAD: usize = 2 * PP;

/// Diagnostics of the last [`commit_paired`], for the profile line only: per worker, the
/// microseconds from the commitment's start to the end of its last group, and the 16-column
/// groups it ran; and the number of sibling pairs that gave up their coupling.
pub static WORKER_END_US: [std::sync::atomic::AtomicU32; 16] = [const { std::sync::atomic::AtomicU32::new(0) }; 16];
pub static WORKER_GROUPS: [std::sync::atomic::AtomicU32; 16] = [const { std::sync::atomic::AtomicU32::new(0) }; 16];
pub static DECOUPLED_PAIRS: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
/// When (microseconds into the commitment, shifted by 4) and by which waiting worker (low 4 bits)
/// each pair decoupled.
pub static DECOUPLED_AT: [std::sync::atomic::AtomicU32; 16] = [const { std::sync::atomic::AtomicU32::new(0) }; 16];

/// The diagnostics above as ` sms_end_us=a/b/.. sms_groups=a/b/.. sms_decoupled=n`.
pub fn worker_report(workers: usize) -> String {
    use std::sync::atomic::Ordering::Relaxed;
    let w = workers.min(16);
    let ends: Vec<String> = (0..w).map(|i| WORKER_END_US[i].load(Relaxed).to_string()).collect();
    let groups: Vec<String> = (0..w).map(|i| WORKER_GROUPS[i].load(Relaxed).to_string()).collect();
    let nd = (DECOUPLED_PAIRS.load(Relaxed) as usize).min(16);
    let at: Vec<String> = (0..nd).map(|i| { let v = DECOUPLED_AT[i].load(Relaxed); format!("{}@{}", v & 15, v >> 4) }).collect();
    format!(" sms_end_us={} sms_groups={} sms_decoupled={} sms_decoupled_at={}", ends.join("/"), groups.join("/"), DECOUPLED_PAIRS.load(Relaxed), at.join("/"))
}

/// The commitment on sibling pairs: hyperthreads `2k`, `2k + 1` (placed on core `k`) take
/// 16-column work units of one ring-element part ([`HALVES`]) each and advance through the key
/// within [`PAIR_LEAD`] pair batches of each other, so both read every key line from their core's
/// L2. Units are handed out dynamically, part-major: a coupled pair claims two adjacent units
/// (two column groups, one part, one key range) at a time, the second for the sibling that asks
/// next. A sibling kept waiting longer than [`SIBLING_WAIT`] (descheduled: the host's writeback
/// of the freshly written table, an interrupt) ends the pair's coupling, and from then on each
/// of the two claims units on its own, so a slowed CPU does less of the work instead of holding
/// up the commitment's end. The two parts of a column group may run on different workers; the
/// second to finish writes its rows.
fn commit_paired(keys: &SmsKeys, columns: &[&[B128]], np: usize, workers: usize, placement: Option<&crate::par::CorePlacement>, dst: usize) {
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    let g = 16usize;
    let pairs = workers / 2;
    let groups = columns.len().div_ceil(g);
    let units = groups * HALVES;
    let locks: Vec<std::sync::Mutex<u8>> = (0..groups).map(|_| std::sync::Mutex::new(0)).collect();
    let next = AtomicUsize::new(0);
    // Per pair: the unit its last pair claim left for the other sibling.
    let pending: Vec<std::sync::Mutex<Option<usize>>> = (0..pairs).map(|_| std::sync::Mutex::new(None)).collect();
    let done: Vec<Progress> = (0..2 * pairs).map(|_| Progress(AtomicUsize::new(0))).collect();
    let coupled: Vec<AtomicBool> = (0..pairs).map(|_| AtomicBool::new(true)).collect();
    let t0 = std::time::Instant::now();
    crate::par::pool().run(2 * pairs, |worker| {
        let _guard = placement.and_then(|p| p.pin(worker));
        let pair = worker / 2;
        let mine = &done[worker].0;
        let other = &done[worker ^ 1].0;
        let coupled = &coupled[pair];
        let mut ws: Option<SmsWorker> = None;
        let mut ran = 0u32;
        // Pair batches of finished units: the siblings' common clock while they are coupled.
        let mut progressed = 0usize;
        loop {
            let unit = {
                let mut left = pending[pair].lock().unwrap();
                match left.take() {
                    Some(u) => Some(u),
                    None => {
                        let two = coupled.load(Ordering::Relaxed);
                        let u = next.fetch_add(if two { 2 } else { 1 }, Ordering::Relaxed);
                        if two && u + 1 < units {
                            *left = Some(u + 1);
                        }
                        (u < units).then_some(u)
                    }
                }
            };
            let Some(u) = unit else { break };
            let (part, gi) = (u / groups, u % groups);
            let range = part * np / HALVES..(part + 1) * np / HALVES;
            let lo = gi * g;
            let n = g.min(columns.len() - lo);
            let ws = ws.get_or_insert_with(SmsWorker::new);
            ran += 1;
            let (from, base) = (range.start, progressed);
            unsafe {
                let out = if HALVES == 1 { RowsOut::Full } else { RowsOut::Half(&locks[gi]) };
                group(&columns[lo..lo + n], keys, np, range.clone(), out, ws, (dst as *mut [u32; N]).add(2 * lo), |b_done| {
                    let d = base + (b_done - from);
                    mine.store(d, Ordering::Release);
                    if !coupled.load(Ordering::Relaxed) {
                        return;
                    }
                    if other.load(Ordering::Acquire) + PAIR_LEAD >= d {
                        return;
                    }
                    let start = std::time::Instant::now();
                    let mut spins = 0u32;
                    while other.load(Ordering::Acquire) + PAIR_LEAD < d {
                        core::hint::spin_loop();
                        spins += 1;
                        if spins % 64 == 0 && start.elapsed() > SIBLING_WAIT {
                            coupled.store(false, Ordering::Relaxed);
                            let k = DECOUPLED_PAIRS.fetch_add(1, Ordering::Relaxed) as usize;
                            if k < 16 {
                                DECOUPLED_AT[k].store(((t0.elapsed().as_micros() as u32) << 4) | worker as u32, Ordering::Relaxed);
                            }
                            return;
                        }
                    }
                })
            };
            progressed += range.len();
            mine.store(progressed, Ordering::Release);
        }
        mine.store(usize::MAX / 2, Ordering::Release);
        if worker < 16 {
            WORKER_END_US[worker].store(t0.elapsed().as_micros() as u32, Ordering::Relaxed);
            WORKER_GROUPS[worker].store(ran, Ordering::Relaxed);
        }
    });
}

/// The commitment's rows: `rows[2c]` is column `c` modulo 9721 in the NTT domain (canonical),
/// `rows[2c + 1]` its coefficients modulo `2^22`.
pub fn commit_columns_rows(keys: &SmsKeys, columns: &[&[B128]], workers: usize) -> Vec<[u32; N]> {
    assert!(!columns.is_empty(), "at least one witness column");
    let len = columns[0].len();
    for c in columns {
        assert_eq!(c.len(), len);
    }
    assert_eq!(len % 8, 0);
    let np = len / 8;
    assert_eq!(keys.elements, 2 * np);
    assert_eq!(keys.pairs.len(), (np / LOCK) * LOCK_WORDS);
    let g = 16usize;
    let groups = columns.len().div_ceil(g);
    let mut rows: Vec<[u32; N]> = Vec::with_capacity(columns.len() * 2);
    let dst = rows.as_mut_ptr() as usize;
    let paired = if workers >= 2 && workers % 2 == 0 && groups > 1 { crate::par::CorePlacement::for_sibling_pairs(workers) } else { None };
    if let Some(placement) = paired {
        commit_paired(keys, columns, np, workers, Some(&placement), dst);
    } else {
        let workers = workers.clamp(1, groups);
        let next = std::sync::atomic::AtomicUsize::new(0);
        let placement = crate::par::CorePlacement::for_logical(workers).or_else(|| crate::par::CorePlacement::for_workers(workers));
        crate::par::pool().run(workers, |worker| {
            let _guard = placement.as_ref().and_then(|p| p.pin(worker));
            let mut ws: Option<SmsWorker> = None;
            loop {
                let gi = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if gi >= groups {
                    break;
                }
                let c0 = gi * g;
                let gc = g.min(columns.len() - c0);
                let ws = ws.get_or_insert_with(SmsWorker::new);
                unsafe { group(&columns[c0..c0 + gc], keys, np, 0..np, RowsOut::Full, ws, (dst as *mut [u32; N]).add(2 * c0), |_| {}) };
            }
        });
    }
    unsafe { rows.set_len(columns.len() * 2) };
    rows
}

/// Column `c`'s coefficients modulo [`MODULUS`] from its two rows (reference recombination).
pub fn column_coefficients(base: &[u32; N], tail: &[u32; N]) -> Vec<u64> {
    let mut out = vec![0i16; N];
    crate::bd::coefficients_from_rows(P, &[base], &mut out);
    let inv = inv_tail();
    (0..N)
        .map(|i| {
            let r = (out[i] as i64).rem_euclid(P as i64) as u64;
            let t = tail[i] as u64;
            // v = t + 2^22 k with v ≡ r (mod P)
            let k = ((r + P as u64 - t % P as u64) % P as u64) * inv % P as u64;
            t + (k << TAIL_BITS)
        })
        .collect()
}
/// `2^-22 mod 9721`.
pub const fn inv_tail() -> u64 {
    let t = (1u64 << TAIL_BITS) % P as u64;
    let mut x = 1u64;
    while (x * t) % P as u64 != 1 {
        x += 1;
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    /// `t = sum_j A_j s_j mod Q` straight from the coefficients.
    fn reference(a: &[u128], col: &[B128]) -> Vec<u64> {
        let q = MODULUS;
        let elements = col.len() / 4;
        let mut acc = vec![0u128; 2 * N];
        for j in 0..elements {
            let aj = &a[j * N..(j + 1) * N];
            for c in 0..512 {
                let (bit, k) = (c / 4, c % 4);
                if (col[4 * j + k].0 >> bit) & 1 == 1 {
                    for m in 0..N {
                        acc[c + m] = (acc[c + m] + aj[m]) % q;
                    }
                }
            }
        }
        for k in (N..2 * N).rev() {
            let v = acc[k];
            acc[k - 324] = (acc[k - 324] + v) % q;
            acc[k - 648] = (acc[k - 648] + q - v) % q;
        }
        acc[..N].iter().map(|&x| x as u64).collect()
    }

    #[test]
    fn digits_reconstruct() {
        let mut s = 0x1234u64;
        for _ in 0..10000 {
            let x = (xorshift(&mut s) as u128) % MODULUS;
            let d = digits(x);
            let mut v = 0i64;
            for k in (0..DIGITS).rev() {
                v = v * digit_radix(k) + d[k] as i64;
            }
            assert_eq!((v.rem_euclid(1 << TAIL_BITS)) as u128, x % (1 << TAIL_BITS));
            for k in 0..DIGITS - 1 {
                assert!((-(DIGIT_HALF as i8)..=DIGIT_HALF as i8).contains(&d[k]));
            }
            assert!((-(DIGIT_HALF as i8)..=DIGIT_HALF as i8).contains(&d[DIGITS - 1]));
        }
    }

    #[test]
    fn shaped_pairs_preserve_both_coefficients() {
        let modulus = 1i64 << TAIL_BITS;
        let check = |x: u128, y: u128, salt: usize| {
            let (a, b) = digits_pair(x, y, salt);
            let value = |d: [i8; DIGITS]| d.iter().enumerate()
                .map(|(i, &v)| v as i64 * digit_weight(i) as i64).sum::<i64>();
            assert!(a[..DIGITS - 1].iter().chain(&b[..DIGITS - 1])
                .all(|v| (-DIGIT_BOUND..=DIGIT_BOUND).contains(v)));
            assert!((-TOP_DIGIT_BOUND..=TOP_DIGIT_BOUND).contains(&a[DIGITS - 1]));
            assert!((-TOP_DIGIT_BOUND..=TOP_DIGIT_BOUND).contains(&b[DIGITS - 1]));
            let (vx, vy) = (value(a), value(b));
            assert_eq!(vx.rem_euclid(modulus) as u128, x & TAIL_MASK as u128);
            assert_eq!(vy.rem_euclid(modulus) as u128, y & TAIL_MASK as u128);
            let center = |x: u128| {
                let r = (x & TAIL_MASK as u128) as i64;
                if r >= modulus / 2 { r - modulus } else { r }
            };
            assert!(pair_norm(vx, vy) <= pair_norm(center(x), center(y)));
        };
        let edges = [0, 1, (modulus / 3) as u128, (modulus / 2) as u128,
            (2 * modulus / 3) as u128, (modulus - 1) as u128];
        for x in edges { for y in edges { for salt in 0..3 { check(x, y, salt); } } }
        let mut state = 0x628d54a7c13e90bfu64;
        for salt in 0..10000 {
            let x = ((xorshift(&mut state) as u128) << 64) | xorshift(&mut state) as u128;
            let y = ((xorshift(&mut state) as u128) << 64) | xorshift(&mut state) as u128;
            check(x, y, salt);
        }
        for a in -DIGIT_HALF..=DIGIT_HALF { for b in -DIGIT_HALF..=DIGIT_HALF {
            let choices = PAIR_DIGITS[((a + DIGIT_HALF) * DIGIT_RADIX + b + DIGIT_HALF) as usize];
            for &[p, q] in &choices.values[..choices.count] {
                assert_eq!((p as i64 - a).rem_euclid(DIGIT_RADIX), 0);
                assert_eq!((q as i64 - b).rem_euclid(DIGIT_RADIX), 0);
                for u in -2..=2 { for v in -2..=2 {
                    assert!(pair_norm(p as i64, q as i64) <= pair_norm(a + DIGIT_RADIX * u, b + DIGIT_RADIX * v));
                } }
            }
        } }
    }

    #[test]
    fn moments_match_direct() {
        let mut s = 0x77u64;
        let l: [i8; N] = core::array::from_fn(|_| ((xorshift(&mut s) % (2 * TOP_DIGIT_BOUND as u64 + 1)) as i8) - TOP_DIGIT_BOUND);
        let mut s1 = [0i64; N];
        let mut s2 = [0i64; N];
        moment_rows(&l, &mut s1, &mut s2);
        let mut d1 = [0i64; N];
        let mut d2 = [0i64; N];
        for i in 0..512 {
            let mut k = vec![0i64; 2 * N];
            for m in 0..N {
                k[m + i] = l[m] as i64;
            }
            for j in (N..2 * N).rev() {
                let v = k[j];
                k[j - 324] += v;
                k[j - 648] -= v;
            }
            for m in 0..N {
                d1[m] += k[m];
                d2[m] += k[m] * k[m];
            }
        }
        assert_eq!(s1, d1);
        assert_eq!(s2, d2);
    }

    #[test]
    fn folded_reduction_exhaustive_domain() {
        if !std::is_x86_feature_detected!("avx512bw") { return; }
        let limit = acc_after_reduce_c(P) as i32;
        let mut base = -limit;
        while base <= limit {
            let words: [i32; 16] = core::array::from_fn(|lane| (base + lane as i32).min(limit));
            let mut got = [0i32; 16];
            unsafe {
                let x = _mm512_loadu_si512(words.as_ptr().cast());
                _mm512_storeu_si512(got.as_mut_ptr().cast(), red_folded(x));
            }
            for lane in 0..16 {
                assert!(got[lane].abs() < P as i32, "x={}, r={}", words[lane], got[lane]);
                assert_eq!((words[lane] - got[lane]) % P as i32, 0,
                    "x={}, r={}", words[lane], got[lane]);
            }
            base += 16;
        }
    }

    #[test]
    fn paired_inverse_banks_match_scalar_reconstruction() {
        let elements = 512;
        let chunks = chunk_count(elements);
        let mut seed = 0xa2e78139ca040d65;
        let mut mur = vec![0i16; mur_len(elements)];
        for m in &mut mur { *m = (xorshift(&mut seed) % P as u64) as i16 - (P / 2) as i16; }
        let means = [0u32; N];
        let keys = SmsKeys { pairs: &[], mur: &mur, mu: &means, elements };
        let mut batch = zeroed_slice::<Batch32>(1);
        let mut total = zeroed_slice::<Total>(1);
        for row in &mut total[0].t { for v in row { *v = xorshift(&mut seed) as i32; } }
        for d in (0..DIGITS).step_by(2) {
            for c in 0..chunks {
                for row in &mut batch[0].v { for v in row {
                    *v = (xorshift(&mut seed) % P as u64) as i16 - (P / 2) as i16;
                } }
                batch[0].v[0][0] = (P / 2) as i16;
                batch[0].v[0][31] = -(P as i16 / 2);
                let before = total[0].t.to_vec();
                unsafe { reconstruct(&batch[0], d, &keys, c, &mut total[0]) };
                for m in 0..N { for column in 0..16 {
                    let mut expected = before[m][column];
                    for bank in 0..2 {
                        let digit = d + bank;
                        let lane = column + 16 * bank;
                        let mean = mur[2 * (((digit / 2) * chunks + c) * N + m) + digit % 2] as i32;
                        let delta = (batch[0].v[m][lane] as i32 - mean + P as i32 / 2)
                            .rem_euclid(P as i32) - P as i32 / 2;
                        expected = expected.wrapping_add(delta.wrapping_mul(digit_weight(digit) as i32));
                    }
                    assert_eq!(total[0].t[m][column], expected, "pair {d}, chunk {c}, row {m}, column {column}");
                } }
            }
        }
    }

    #[test]
    fn commitment_matches_reference() {
        commitment_matches_reference_at(512, 7.5); // full and short chunks in each half
        commitment_matches_reference_at(1024, 7.5); // an inverse batch crosses the half boundary
    }

    fn commitment_matches_reference_at(elements: usize, floor: f64) {
        let ncols = 37usize;
        let mut s = 0x9e3779b97f4a7c15u64;
        let a: Vec<u128> = (0..elements * N).map(|_| ((xorshift(&mut s) as u128) << 64 | xorshift(&mut s) as u128) % MODULUS).collect();
        let data = SmsKeyData::build(&a, elements);
        eprintln!("margin {}", data.margin);
        assert!(data.margin > floor);
        let len = 4 * elements;
        let words: Vec<B128> = (0..len * ncols).map(|_| B128(((xorshift(&mut s) as u128) << 64) | xorshift(&mut s) as u128)).collect();
        let cols: Vec<&[B128]> = words.chunks_exact(len).collect();
        let keys = data.keys();
        for workers in [1usize, 3, 4] {
            let rows = commit_columns_rows(&keys, &cols, workers);
            for c in [0usize, 1, 15, 16, 31, 36] {
                let got = column_coefficients(&rows[2 * c], &rows[2 * c + 1]);
                let want = reference(&a, cols[c]);
                assert!(got == want, "column {c} with {workers} workers differs");
            }
        }
    }
}

#[cfg(test)]
mod bench {
    use super::*;
    /// Full-shape timing against the hybrid commitment (`cargo test --release -- --ignored sms_bench --nocapture`).
    #[test]
    #[ignore]
    fn sms_bench() {
        let elements = 4096usize;
        let ncols = 1024usize;
        let mut s = 0x9e3779b97f4a7c15u64;
        let mut next = || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
        let a: Vec<u128> = (0..elements * N).map(|_| ((next() as u128) << 64 | next() as u128) % MODULUS).collect();
        let t = std::time::Instant::now();
        let data = SmsKeyData::build(&a, elements);
        eprintln!("keys {:?} margin {:.2}", t.elapsed(), data.margin);
        let len = 4 * elements;
        let words: Vec<B128> = (0..len * ncols).map(|_| B128(((next() as u128) << 64) | next() as u128)).collect();
        let cols: Vec<&[B128]> = words.chunks_exact(len).collect();
        let keys = data.keys();
        let workers = std::env::var("SMS_WORKERS").ok().and_then(|v| v.parse().ok()).unwrap_or(crate::par::threads());
        let _ = crate::par::pool();
        let reps: usize = std::env::var("SMS_REPS").ok().and_then(|v| v.parse().ok()).unwrap_or(7);
        let mut ts = Vec::new();
        let mut chk = 0u32;
        for _ in 0..reps {
            let t = std::time::Instant::now();
            let rows = commit_columns_rows(&keys, &cols, workers);
            ts.push(t.elapsed().as_secs_f64() * 1e3);
            chk ^= rows[1][0];
        }
        ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eprintln!("sms commit min {:.1} med {:.1} ms ({chk})", ts[0], ts[ts.len() / 2]);
        if std::env::var_os("SMS_ONLY").is_some() {return;}
        let base = crate::key::CommitmentKey::random(len, 11, crate::ring::Modulus::Q3889_FS_S, &[crate::ring::Modulus::Q9721_FS_S]);
        let tail = crate::key::CommitmentKey::random(len, 12, crate::ring::Modulus::Q3889_FS_S, &[crate::ring::Modulus::Q9721_FS_S]);
        let pk = super::super::hybrid::pair_keys(&base, &tail, workers);
        for _ in 0..3 {
            let t = std::time::Instant::now();
            let rows = super::super::hybrid::commit_columns_rows_pairs(&pk, &cols, workers);
            eprintln!("hybrid commit {:?} ({})", t.elapsed(), rows[2][0]);
        }
    }
}

#[cfg(test)]
mod bench_parts {
    use super::*;
    #[test]
    #[ignore]
    fn sms_parts_bench() {
        let mut b = zeroed_slice::<Batch32>(1);
        let mut s = 7u64;
        for j in 0..N { for l in 0..32 { s ^= s << 13; s ^= s >> 7; s ^= s << 17; b[0].v[j][l] = ((s % 19441) as i32 - 9720) as i16; } }
        let n = 20000;
        let t = std::time::Instant::now();
        for _ in 0..n {
            b[0].representation = Representation::Ntt;
            unsafe { gen_small::intt_gen_batch32::<P>(&mut b[0]) };
            for j in 0..N { b[0].v[j][0] = b[0].v[j][0].wrapping_mul(3) % 9721; }
        }
        let dt = t.elapsed().as_secs_f64() / n as f64;
        eprintln!("intt_gen_batch32::<9721>: {:.2} us = {:.0} cycles at 3.5 GHz", dt * 1e6, dt * 3.5e9);
        let mut acc = zeroed_slice::<SAcc>(1);
        let t = std::time::Instant::now();
        for _ in 0..n {
            unsafe { finish_stream(&mut acc[0], 3, &mut b[0], 1) };
        }
        let dt = t.elapsed().as_secs_f64() / n as f64;
        eprintln!("finish_stream: {:.2} us", dt * 1e6);
    }
}
