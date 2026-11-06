//! Precomputes the public commitment key of the default modulus `Q = 9721 * 2^22`: the matrix of
//! `sms.rs` (the BLAKE3 expansion of `A` for `Q`, its residues modulo 9721 and its six balanced
//! radix-16 digit polynomials, all transformed, plus every chunk's exact digit means and its
//! checked reconstruction margin). It depends only on the modulus and the instance shape, never
//! on a seed, a table or the transcript, and is exactly what `sms::matrix` computes at run time;
//! `main.rs` embeds it and computes the matrix as before when the image is absent.
#![allow(dead_code, unused_imports, unused_variables)]

#[path = "src/hybrid.rs"]
mod hybrid;
#[path = "src/ntt.rs"]
mod ntt;
#[path = "src/rq.rs"]
mod rq;
#[path = "src/sms.rs"]
mod sms;

/// `main.rs`'s `RING_ELEMENTS`: `2^(31 - 7 - 10) / 4`.
const RING_ELEMENTS: usize = 4096;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    for f in ["src/ntt.rs", "src/ntt_neon.rs", "src/sms_dot.rs", "src/sms.rs", "src/hybrid.rs", "src/rq.rs"] {
        println!("cargo:rerun-if-changed={f}");
    }
    let dir = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let (data, means, margins) = image();
    std::fs::write(dir.join("sms-data.bin"), data).unwrap();
    std::fs::write(dir.join("sms-means.bin"), means).unwrap();
    std::fs::write(dir.join("sms-margins.bin"), margins).unwrap();
}

#[cfg(target_arch = "aarch64")]
fn image() -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get());
    let m = sms::matrix(RING_ELEMENTS, workers);
    let data = m.chunks.iter().flat_map(|c| c.data.iter().flat_map(|v| v.to_le_bytes())).collect();
    let means = m.chunks.iter().flat_map(|c| c.means.iter().flat_map(|v| v.to_le_bytes())).collect();
    let margins = m.chunks.iter().flat_map(|c| c.margin.to_le_bytes()).collect();
    (data, means, margins)
}

#[cfg(not(target_arch = "aarch64"))]
fn image() -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    (Vec::new(), Vec::new(), Vec::new())
}
