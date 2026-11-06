#![allow(dead_code)]
#[path = "src/ntt.rs"]
mod ntt;
#[path = "src/rq.rs"]
mod rq;
#[path = "src/hybrid.rs"]
mod hybrid;

fn main() {
    for path in ["build.rs", "src/ntt.rs", "src/ntt_neon.rs", "src/rq.rs", "src/hybrid.rs"] {
        println!("cargo:rerun-if-changed={path}");
    }
    let matrix = hybrid::matrix(4096, 8);
    let mut image = Vec::with_capacity(matrix.data.len() * 2);
    for &coefficient in matrix.data.iter() {
        image.extend_from_slice(&coefficient.to_le_bytes());
    }
    let directory = std::path::PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    std::fs::write(directory.join("hybrid-matrix.bin"), image).unwrap();
}
