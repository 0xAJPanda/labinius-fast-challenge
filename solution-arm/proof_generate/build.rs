//! Build only the instance-independent public key, using the runtime generator verbatim.
#![allow(dead_code)]
#[path = "src/ntt.rs"] mod ntt;
#[path = "src/rq.rs"] mod rq;
#[path = "src/hybrid.rs"] mod hybrid;
#[path = "src/sms.rs"] mod sms;

fn main() {
    for path in ["build.rs", "src/ntt.rs", "src/ntt_neon.rs", "src/rq.rs", "src/hybrid.rs", "src/sms.rs", "src/sms_dot.rs"] {
        println!("cargo:rerun-if-changed={path}");
    }
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("linux") {
        println!("cargo:rustc-link-arg-bins=-Wl,-z,max-page-size=0x200000");
    }
    let elements = 1 << 12;
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get());
    let matrix = sms::matrix(elements, workers);
    let mut bytes = Vec::new();
    let mut layout = Vec::new();
    for chunk in matrix.chunks {
        let data_offset = bytes.len();
        bytes.extend(chunk.data.iter().flat_map(|x| x.to_le_bytes()));
        while bytes.len() % 4 != 0 { bytes.push(0); }
        let means_offset = bytes.len();
        bytes.extend(chunk.means.iter().flat_map(|x| x.to_le_bytes()));
        layout.push((chunk.elements, data_offset, chunk.data.len(), means_offset, chunk.means.len(), chunk.margin.to_bits()));
    }
    let dir = std::path::PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    std::fs::write(dir.join("matrix.bin"), &bytes).expect("write public matrix");
    std::fs::write(dir.join("matrix_params.rs"), format!(
        "const KEY_BYTES: usize = {};\nconst KEY_MODULUS: u128 = {};\nconst KEY_PRIMES: [u32; 3] = {:?};\nconst KEY_ELEMENTS: usize = {};\nconst KEY_N: usize = {};\nconst KEY_CHUNK: usize = {};\nconst KEY_DIGITS: usize = {};\nconst KEY_STREAMS: usize = {};\nconst KEY_LAYOUT: [(usize, usize, usize, usize, usize, u64); {}] = {:?};\n",
        bytes.len(), sms::MODULUS, ntt::PRIMES, elements, ntt::N, sms::CHUNK, sms::DIGITS, sms::STREAMS, layout.len(), layout,
    )).expect("write public matrix parameters");
}
