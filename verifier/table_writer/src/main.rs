//! Writes the witness table for `LABINIUS_SEED`, byte for byte what
//! `labinius-fast-verify --write-table` writes: the BLAKE3 XOF keyed by the seed, updated with
//! [`TABLE_CONTEXT`], filled to `(1 << LOG_VARS) / 8` bytes. Verification never reads this file;
//! it rebuilds the table from the seed.
const LOG_VARS: usize = 31;
const TABLE_CONTEXT: &[u8] = b"labinius-fast/binary-table/v1";

fn seed() -> [u8; 32] {
    let s = std::env::var("LABINIUS_SEED").expect("LABINIUS_SEED is required");
    assert_eq!(s.len(), 64);
    core::array::from_fn(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).expect("invalid seed"))
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(args.len() == 2, "usage: labinius-table-writer FILE");
    let mut hasher = blake3::Hasher::new_keyed(&seed());
    hasher.update(TABLE_CONTEXT);
    let mut bytes = vec![0; (1 << LOG_VARS) / 8];
    hasher.finalize_xof().fill(&mut bytes);
    std::fs::write(&args[1], bytes).expect("write witness table");
}
