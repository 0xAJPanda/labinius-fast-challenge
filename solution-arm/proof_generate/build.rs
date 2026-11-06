//! Precompute the commitment matrix of the default (hybrid) modulus in the transform domain.
//! `A` is public parameter data: a function of the fixed matrix seed, the modulus and the element
//! count alone, never of the instance seed. So it is computed here, once, at build time, with the
//! prover's own NEON kernels, and embedded in the binary (the x86 track embeds its key the same
//! way). On a build host that is not AArch64, an empty image is written and the prover builds the
//! matrix when it runs.
#![allow(dead_code, unused_imports)]
use std::env;
use std::path::PathBuf;

#[cfg(target_arch = "aarch64")]
#[path = "src/ntt.rs"]
mod ntt;
#[cfg(target_arch = "aarch64")]
#[path = "src/rq.rs"]
mod rq;

/// Ring elements per column (`parameters`: 2^24 words / 1024 columns / 4); the prover checks the
/// image's length before it uses it.
const RING_ELEMENTS: usize = 4096;

fn main() {
    let dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    for f in ["build.rs", "src/ntt.rs", "src/ntt_neon.rs", "src/ntt_hybrid.rs", "src/rq.rs"] {
        println!("cargo:rerun-if-changed={f}");
    }
    println!("cargo:rerun-if-env-changed=CARGO_CFG_TARGET_ARCH");
    std::fs::write(dir.join("hybrid-matrix.bin"), image()).expect("write the matrix image");
}

#[cfg(target_arch = "aarch64")]
fn image() -> Vec<u8> {
    // The build script runs on the host; embed only for an AArch64 target as well.
    if env::var("CARGO_CFG_TARGET_ARCH").as_deref() != Ok("aarch64") {
        return Vec::new();
    }
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get());
    let m = ntt::neon::hybrid::matrix(RING_ELEMENTS, workers);
    m.data.iter().flat_map(|v| v.to_le_bytes()).collect()
}

#[cfg(not(target_arch = "aarch64"))]
fn image() -> Vec<u8> {
    Vec::new()
}
