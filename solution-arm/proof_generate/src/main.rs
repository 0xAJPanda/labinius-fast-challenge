//! The portable track's prover: the same CLI, transcript, and proof bytes as the x86 track's
//! `labinius-fast`, computed without AVX-512.
//!
//! The commitment `C_j = A W_j` is computed directly in `R_Q` rather than through the library's
//! per-prime NTT kernels. A witness ring element is binary, so `A_m W_{j,m}` is a sum of shifted
//! copies of `A_m`; grouping eight consecutive coefficients turns each group into one lookup of a
//! 256-entry table of precomputed partial sums (the "four Russians" trick). The field switch,
//! row evaluation, and fold follow the checker's reference prover.
#[allow(dead_code)]
mod rq;

use bincode::Options;
use labinius::challenge::{sample_short_challenge, ShortChallenge, DEFAULT_BOUND, DEFAULT_WEIGHT};
use labinius::eval::eq_table;
use labinius::fields::{
    crossfield as cf,
    scalar::{B128, F162},
};
use labinius::ring::Representation;
use labinius::scheme::FoldingSource;
use labinius::{wire, FoldedWitness, RingElement648, RowEvaluation, Transcript};
use serde::Serialize;

const LOG_VARS: usize = 31;
const LOG_PACKED: usize = LOG_VARS - cf::LOG_PACK;
const COLUMN_LOG: u32 = 10;
const COLUMNS: usize = 1 << COLUMN_LOG;
const ROW_LOG: usize = LOG_PACKED - COLUMN_LOG as usize;
const COLUMN_LEN: usize = 1 << ROW_LOG;
const RING_ELEMENTS: usize = COLUMN_LEN / 4;
const LIMIT: u64 = 64 * 1024 * 1024;
const MIN_MODULUS: u128 = 34_816_965_049;
/// The x86 track's default, so both tracks write the same proof for the same seed.
const MODULUS: u128 = 40_772_829_184;

/// Coefficients of a witness ring element per table lookup in [`commit`].
const GROUP: usize = 8;
/// Coefficients a binary witness ring element can have set: F162 lifts B128, so bit `p < 128`.
const WITNESS_COEFFICIENTS: usize = 4 * 128;
const GROUPS: usize = WITNESS_COEFFICIENTS / GROUP;
/// One table row: `A_m` shifted by up to `GROUP - 1`, padded to an even length.
const ROW: usize = (GROUP - 1 + rq::N + 1) & !1;
/// The unreduced product of a witness column with `A`, padded like [`ROW`].
const RAW: usize = (GROUPS - 1) * GROUP + ROW;
/// Largest modulus [`commit`]'s `u64` accumulator holds: a raw coefficient sums at most
/// `WITNESS_COEFFICIENTS * RING_ELEMENTS = 2^21` residues below `Q`.
const FAST_MODULUS: u128 = 1 << 43;

#[derive(Serialize)]
struct Proof {
    version: u32,
    modulus: u128,
    commitment: Vec<u8>,
    partials: Vec<u128>,
    rounds: Vec<[[u64; 3]; 2]>,
    opened: [u64; 3],
    row: Vec<u8>,
    folded: Vec<u8>,
}

fn codec() -> impl Options {
    bincode::DefaultOptions::new()
        .with_fixint_encoding()
        .with_limit(LIMIT)
        .reject_trailing_bytes()
}

fn seed() -> [u8; 32] {
    let s = std::env::var("LABINIUS_SEED").expect("LABINIUS_SEED is required");
    assert_eq!(s.len(), 64);
    core::array::from_fn(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).expect("invalid seed"))
}

fn parse_modulus(text: &str) -> Result<u128, String> {
    let q: u128 = text
        .parse()
        .map_err(|_| format!("not a decimal integer: {text}"))?;
    if q < MIN_MODULUS {
        return Err("commitment modulus is below the security bound".into());
    }
    if q >= rq::MAX_MODULUS {
        return Err("commitment modulus is above what the checker's arithmetic covers".into());
    }
    Ok(q)
}

fn threads() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

fn read_table(path: &str) -> Vec<B128> {
    let bytes = std::fs::read(path).expect("read witness table");
    assert_eq!(
        bytes.len(),
        (1 << LOG_VARS) / 8,
        "witness table has the wrong length"
    );
    bytes
        .chunks_exact(16)
        .map(|chunk| B128(u128::from_le_bytes(chunk.try_into().unwrap())))
        .collect()
}

fn transcript(seed: &[u8; 32], modulus: u128, commitment: &[u8]) -> Transcript {
    let mut t = Transcript::new(b"labinius-fast/v1");
    t.absorb_bytes(seed);
    t.absorb_u64(LOG_VARS as u64);
    t.absorb_u64(COLUMN_LOG as u64);
    t.absorb_bytes(&modulus.to_le_bytes());
    t.absorb_bytes(commitment);
    t
}

/// The table of group `pattern`s for ring element `element` of `A`: row `s` is the sum of
/// `A_m X^r` over the set bits `r` of `s`. Each row is a previous row plus one shifted copy.
fn build_table(table: &mut [u64], element: &[u64]) {
    table[..ROW].fill(0);
    for s in 1usize..1 << GROUP {
        let r = s.ilog2() as usize;
        let (done, rest) = table.split_at_mut(s * ROW);
        let row = &mut rest[..ROW];
        row.copy_from_slice(&done[(s & !(1 << r)) * ROW..][..ROW]);
        for (slot, &a) in row[r..r + rq::N].iter_mut().zip(element) {
            *slot += a;
        }
    }
}

/// Group `g` of ring element `m` of a column as a table index. Coefficient `4 p + k` is bit `p`
/// of word `k`, so the group's eight coefficients are bits `2g` and `2g + 1` of the four words.
#[inline(always)]
fn pattern(words: &[u128; 4], g: usize) -> usize {
    let mut s = 0;
    for (k, &w) in words.iter().enumerate() {
        let pair = (w >> (2 * g)) as usize & 3;
        s |= (pair & 1) << k | (pair >> 1) << (4 + k);
    }
    s
}

/// A raw product reduced modulo `q` and folded into `R_Q` with `X^648 = X^324 - 1`.
fn reduce(raw: &[u64], q: u64) -> Vec<u128> {
    let mut r = vec![0u64; 2 * rq::N - 1];
    for (slot, &x) in r.iter_mut().zip(raw) {
        *slot = x % q;
    }
    for k in (rq::N..2 * rq::N - 1).rev() {
        let v = r[k];
        let up = &mut r[k - rq::SPLIT];
        *up += v;
        if *up >= q {
            *up -= q;
        }
        let down = &mut r[k - rq::N];
        *down += q - v;
        if *down >= q {
            *down -= q;
        }
    }
    r[..rq::N].iter().map(|&x| x as u128).collect()
}

/// Every column's commitment, flattened. Columns go out over the cores; each worker walks the
/// ring elements once, building that element's table and applying it to all of its columns.
fn commit(values: &[B128], modulus: u128) -> Vec<u128> {
    let matrix = rq::ajtai_matrix(modulus, RING_ELEMENTS);
    if modulus >= FAST_MODULUS {
        return commit_reference(values, &matrix, modulus);
    }
    let matrix: Vec<u64> = matrix.iter().map(|&a| a as u64).collect();
    let q = modulus as u64;
    let per = COLUMNS.div_ceil(threads());
    let mut out = vec![0u128; COLUMNS * rq::N];
    std::thread::scope(|scope| {
        for (block, first) in out.chunks_mut(per * rq::N).zip((0..).step_by(per)) {
            let matrix = &matrix;
            scope.spawn(move || {
                let count = block.len() / rq::N;
                let mut table = vec![0u64; ROW << GROUP];
                let mut raw = vec![0u64; count * RAW];
                for (m, element) in matrix.chunks_exact(rq::N).enumerate() {
                    build_table(&mut table, element);
                    for (c, acc) in raw.chunks_exact_mut(RAW).enumerate() {
                        let base = (first + c) * COLUMN_LEN + 4 * m;
                        let words: [u128; 4] = core::array::from_fn(|k| values[base + k].0);
                        for g in 0..GROUPS {
                            let s = pattern(&words, g);
                            if s == 0 {
                                continue;
                            }
                            let row = &table[s * ROW..][..ROW];
                            for (slot, &t) in acc[g * GROUP..][..ROW].iter_mut().zip(row) {
                                *slot += t;
                            }
                        }
                    }
                }
                for (target, acc) in block.chunks_exact_mut(rq::N).zip(raw.chunks_exact(RAW)) {
                    target.copy_from_slice(&reduce(acc, q));
                }
            });
        }
    });
    out
}

/// The checker's schoolbook commitment, for a modulus too wide for [`commit`]'s `u64` lanes.
fn commit_reference(values: &[B128], matrix: &[u128], modulus: u128) -> Vec<u128> {
    let per = COLUMNS.div_ceil(threads());
    let mut out = vec![0u128; COLUMNS * rq::N];
    std::thread::scope(|scope| {
        for (block, columns) in out
            .chunks_mut(per * rq::N)
            .zip(values.chunks(per * COLUMN_LEN))
        {
            scope.spawn(move || {
                for (target, column) in block
                    .chunks_mut(rq::N)
                    .zip(columns.chunks_exact(COLUMN_LEN))
                {
                    let witness: Vec<F162> = column.iter().map(|&x| F162::from_b128(x)).collect();
                    target.copy_from_slice(&rq::commit_binary_column(matrix, &witness, modulus));
                }
            });
        }
    });
    out
}

fn dot_lifted(eq: &[F162], column: &[B128]) -> F162 {
    eq.iter()
        .zip(column)
        .fold(F162::ZERO, |acc, (&e, &x)| acc + e * F162::from_b128(x))
}

/// `u_j = sum_i eq(p0, i) W[i + COLUMN_LEN j]`, one field element per column.
fn row_evaluate(values: &[B128], p0: &[F162]) -> Vec<F162> {
    let eq = eq_table(p0);
    let per = COLUMNS.div_ceil(threads());
    let mut out = vec![F162::ZERO; COLUMNS];
    std::thread::scope(|scope| {
        for (block, columns) in out.chunks_mut(per).zip(values.chunks(per * COLUMN_LEN)) {
            let eq = &eq;
            scope.spawn(move || {
                for (slot, column) in block.iter_mut().zip(columns.chunks_exact(COLUMN_LEN)) {
                    *slot = dot_lifted(eq, column);
                }
            });
        }
    });
    out
}

/// The sumcheck coordinates split into the row point `p0` and the column point `p1`.
fn split_point(coords: &[F162]) -> (Vec<F162>, Vec<F162>) {
    assert_eq!(
        coords.len(),
        LOG_PACKED,
        "one coordinate per packed variable"
    );
    let (high, low) = coords.split_at(COLUMN_LOG as usize);
    (
        low.iter().rev().copied().collect(),
        high.iter().rev().copied().collect(),
    )
}

fn folding_challenges(t: &mut Transcript, row: &RowEvaluation) -> Vec<ShortChallenge> {
    row.absorb(t);
    (0..COLUMNS)
        .map(|_| sample_short_challenge(t, DEFAULT_WEIGHT, DEFAULT_BOUND).0)
        .collect()
}

/// `v_m = sum_j c_j W_{j,m}` over the integers. A challenge sits in `R_648` as `c(-X^4)`, so each
/// of its terms adds or subtracts the column's binary element shifted by a multiple of four.
fn fold_element(values: &[B128], challenges: &[ShortChallenge], m: usize) -> RingElement648 {
    let mut raw = [0i32; 2 * rq::N];
    let mut bits = [0i32; rq::N];
    for (j, challenge) in challenges.iter().enumerate() {
        let packed = &values[j * COLUMN_LEN + 4 * m..][..4];
        for (k, x) in packed.iter().enumerate() {
            for p in 0..128 {
                bits[4 * p + k] = ((x.0 >> p) & 1) as i32;
            }
        }
        for t in 0..challenge.weight {
            let position = challenge.positions[t] as usize;
            let negative = ((challenge.signs >> t) & 1 == 1) ^ (position % 2 == 1);
            let window = &mut raw[4 * position..][..WITNESS_COEFFICIENTS];
            if negative {
                for (slot, &b) in window.iter_mut().zip(&bits) {
                    *slot -= b;
                }
            } else {
                for (slot, &b) in window.iter_mut().zip(&bits) {
                    *slot += b;
                }
            }
        }
    }
    for k in (rq::N..2 * rq::N - 1).rev() {
        let value = raw[k];
        raw[k - rq::SPLIT] += value;
        raw[k - rq::N] -= value;
    }
    let mut element = RingElement648::zero(Representation::Coefficients);
    for (slot, &value) in element.v.iter_mut().zip(&raw[..rq::N]) {
        *slot = i16::try_from(value).expect("a fold coefficient does not fit the wire codec");
    }
    element
}

fn fold(values: &[B128], challenges: &[ShortChallenge]) -> Vec<RingElement648> {
    let per = RING_ELEMENTS.div_ceil(threads());
    let mut out: Vec<RingElement648> = (0..RING_ELEMENTS)
        .map(|_| RingElement648::zero(Representation::Coefficients))
        .collect();
    std::thread::scope(|scope| {
        for (block, base) in out.chunks_mut(per).zip((0..).step_by(per)) {
            scope.spawn(move || {
                for (i, slot) in block.iter_mut().enumerate() {
                    *slot = fold_element(values, challenges, base + i);
                }
            });
        }
    });
    out
}

fn pack_folded(folded: &FoldedWitness) -> Vec<u8> {
    let mut out = Vec::with_capacity(4 * rq::N * folded.elements().len());
    for e in folded.elements() {
        for &c in &e.v {
            out.extend_from_slice(&(c as i32).to_le_bytes());
        }
    }
    out
}

fn prove(seed: &[u8; 32], modulus: u128, values: &[B128], timer: &mut Timer) -> Proof {
    let commitment = rq::pack_ring_elements(&commit(values, modulus), modulus);
    timer.lap("commit");
    let mut t = transcript(seed, modulus, &commitment);
    let lo = t.sample_b128(b"r-lo", 7);
    let hi = t.sample_b128(b"r-hi", LOG_PACKED);
    let mut ws = cf::Workspace::new(LOG_PACKED, 1);
    let partials = ws.partial_evals(0, values, &hi);
    let value = partials
        .iter()
        .zip(cf::eq_expand_b128(&lo))
        .fold(B128::ZERO, |a, (&v, e)| a + v * e);
    t.absorb_b128(&[value]);
    t.absorb_b128(&partials);
    let batch = cf::eq_expand_f162(&t.sample_f162(b"switch-batch", 7));
    let mut switch = cf::SwitchProver::from_workspace(ws, values, &[F162::ONE], &batch);
    let mut rounds = Vec::with_capacity(LOG_PACKED);
    let mut coords = Vec::with_capacity(LOG_PACKED);
    for _ in 0..LOG_PACKED {
        let msg = switch.msg();
        t.absorb_f162(&msg);
        let r = t.sample_f162(b"switch-round", 1)[0];
        rounds.push(msg.map(|v| v.0));
        coords.push(r);
        switch.fold(r);
    }
    let opened = switch.final_eval();
    t.absorb_f162(&[opened]);
    timer.lap("switch");
    let (p0, _) = split_point(&coords);
    let row = RowEvaluation::of(row_evaluate(values, &p0));
    timer.lap("row");
    let challenges = folding_challenges(&mut t, &row);
    let folded = FoldedWitness::of(fold(values, &challenges));
    timer.lap("fold");
    Proof {
        version: 1,
        modulus,
        commitment,
        partials: partials.iter().map(|v| v.0).collect(),
        rounds,
        opened: opened.0,
        row: wire::pack_row_evaluation(&row),
        folded: pack_folded(&folded),
    }
}

/// Phase times on stderr, diagnostic only.
struct Timer(std::time::Instant, Vec<String>);

impl Timer {
    fn lap(&mut self, phase: &str) {
        self.1
            .push(format!("{phase}_ms={}", self.0.elapsed().as_millis()));
        self.0 = std::time::Instant::now();
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(
        args.len() == 5 || args.len() == 7,
        "usage: labinius-fast --prove FILE --table FILE [--modulus Q]"
    );
    assert_eq!(args[1], "--prove", "unknown command");
    assert_eq!(args[3], "--table", "witness table path is required");
    let seed = seed();
    let modulus = if args.len() == 7 {
        assert_eq!(args[5], "--modulus", "unknown prover option");
        parse_modulus(&args[6]).expect("invalid commitment modulus")
    } else {
        MODULUS
    };
    let mut timer = Timer(std::time::Instant::now(), Vec::new());
    let values = read_table(&args[4]);
    timer.lap("read");
    let proof = prove(&seed, modulus, &values, &mut timer);
    let bytes = codec().serialize(&proof).expect("serialize proof");
    std::fs::write(&args[2], bytes).expect("write proof");
    timer.lap("write");
    eprintln!("labinius-arm {}", timer.1.join(" "));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn column(seed: u64) -> Vec<B128> {
        let mut state = seed;
        (0..COLUMN_LEN)
            .map(|_| {
                let mut word = 0u128;
                for _ in 0..2 {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    word = word << 64 | state as u128;
                }
                B128(word)
            })
            .collect()
    }

    /// The table lookup agrees with the checker's coefficient-by-coefficient commitment.
    #[test]
    fn grouped_commit_matches_reference() {
        let modulus = MODULUS;
        let values: Vec<B128> = (0..COLUMNS as u64)
            .flat_map(|j| {
                if j < 2 {
                    column(j + 1)
                } else {
                    vec![B128::ZERO; COLUMN_LEN]
                }
            })
            .collect();
        let fast = commit(&values, modulus);
        let matrix = rq::ajtai_matrix(modulus, RING_ELEMENTS);
        for j in 0..2 {
            let witness: Vec<F162> = values[j * COLUMN_LEN..][..COLUMN_LEN]
                .iter()
                .map(|&x| F162::from_b128(x))
                .collect();
            assert_eq!(
                fast[j * rq::N..][..rq::N],
                rq::commit_binary_column(&matrix, &witness, modulus)[..],
                "column {j}"
            );
        }
        assert!(fast[2 * rq::N..].iter().all(|&x| x == 0));
    }
}
