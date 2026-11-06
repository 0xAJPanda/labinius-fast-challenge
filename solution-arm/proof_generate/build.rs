#![allow(dead_code)]
#[path = "src/ntt.rs"]
mod ntt;
#[path = "src/rq.rs"]
mod rq;
#[path = "src/hybrid.rs"]
mod hybrid;
#[path = "src/sms.rs"]
mod sms;

fn main() {
    for path in ["build.rs", "src/ntt.rs", "src/ntt_neon.rs", "src/rq.rs", "src/hybrid.rs", "src/sms.rs", "src/sms_dot.rs"] {
        println!("cargo:rerun-if-changed={path}");
    }
    let matrix = sms::matrix(4096, 8);
    let mut data = Vec::new();
    let mut sums = Vec::new();
    let mut margins = Vec::new();
    for chunk in &matrix.chunks {
        for &x in chunk.data.iter() { data.extend_from_slice(&x.to_le_bytes()); }
        for &x in chunk.sums.iter() { sums.extend_from_slice(&x.to_le_bytes()); }
        margins.extend_from_slice(&chunk.margin.to_le_bytes());
    }
    let out = std::path::PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    std::fs::write(out.join("sms-data.bin"), data).unwrap();
    std::fs::write(out.join("sms-sums.bin"), sums).unwrap();
    std::fs::write(out.join("sms-margins.bin"), margins).unwrap();
}
