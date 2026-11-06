//! Precomputes the public commitment key of the default modulus: the hybrid matrix of
//! `hybrid.rs` (the BLAKE3 expansion of `A` for `Q = 3889 * 9721 * 1024`, its residues and
//! centered tail lifts transformed, with the worst-case CRT certificate of every chunk checked).
//! It depends only on the modulus and the instance shape, never on a seed or a table, and is
//! exactly the array `hybrid::matrix` computes at run time; `main.rs` embeds it and computes
//! the matrix as before when the image is empty (non-AArch64 build hosts).
#![allow(dead_code, unused_imports)]

#[path = "src/hybrid.rs"]
mod hybrid;
#[path = "src/ntt.rs"]
mod ntt;
#[path = "src/rq.rs"]
mod rq;

/// `main.rs`'s `RING_ELEMENTS`: `2^(31 - 7 - 10) / 4`.
const RING_ELEMENTS: usize = 4096;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    for f in ["src/ntt.rs", "src/ntt_neon.rs", "src/hybrid.rs", "src/rq.rs"] {
        println!("cargo:rerun-if-changed={f}");
    }
    let out = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("hybrid-key.bin");
    std::fs::write(out, image()).unwrap();
}

#[cfg(target_arch = "aarch64")]
fn image() -> Vec<u8> {
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get());
    let m = hybrid::matrix(RING_ELEMENTS, workers);
    m.data.iter().flat_map(|v| v.to_le_bytes()).collect()
}

#[cfg(not(target_arch = "aarch64"))]
fn image() -> Vec<u8> {
    Vec::new()
}
