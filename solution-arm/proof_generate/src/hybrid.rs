//! Two witness transforms, with an exactly reconstructed residue modulo 1024.
//! The certificate bounds every possible binary witness, not sampled inputs.
use crate::ntt::{self, N, ROW, for_blocks, centered, primes};
use labinius::fields::scalar::B128;
pub const PRODUCT: i64 = 3889 * 9721;
pub const MODULUS: u128 = PRODUCT as u128 * 1024;
pub const CHUNK: usize = 224;
const _: () = assert!(MODULUS >= 34_816_965_049 && MODULUS < (1u128 << 96));

pub struct Matrix { pub elements: usize, pub data: Vec<i16> }
impl Matrix {
    pub fn row(&self, e: usize, p: usize, tail: usize) -> &[i16] {
        &self.data[(4*e+2*p+tail)*ROW..][..ROW]
    }
}

/// Exact extrema of each output coefficient over independent witness bits.
/// Prefix sums combine overlapping monomial reductions before taking signs.
pub fn certificate(raw: &[[i16; N]]) -> i64 {
    let mut pos = [0i64; N]; let mut neg = pos;
    for a in raw {
        let mut ap = [0i64; 649]; let mut an = ap;
        let mut sp = [0i64; 325]; let mut sn = sp;
        for j in 0..648 { let v=a[j] as i64; ap[j+1]=ap[j]+v.max(0); an[j+1]=an[j]+(-v).max(0); }
        for j in 0..324 { let v=a[j] as i64+a[j+324] as i64; sp[j+1]=sp[j]+v.max(0); sn[j+1]=sn[j]+(-v).max(0); }
        let range = |p: &[i64], n: &[i64], lo: usize, hi: usize| {
            if lo>hi {(0,0)} else {(p[hi+1]-p[lo], n[hi+1]-n[lo])}
        };
        for j in 0..648 {
            let mut add = |(p,n): (i64,i64), minus: bool| {
                if minus {pos[j]+=n; neg[j]+=p;} else {pos[j]+=p; neg[j]+=n;}
            };
            if j<324 {
                add(range(&ap,&an,0,j),false);
                if j<=186 { add(range(&sp,&sn,j+137,323),true); add(range(&ap,&an,324,647),true); }
                else { add(range(&ap,&an,j+137,647),true); }
            } else {
                add(range(&sp,&sn,j.saturating_sub(511),323),false);
                add(range(&ap,&an,324,j),false);
                if j<511 {add(range(&ap,&an,j-187,323),false);}
            }
        }
    }
    let worst = pos.into_iter().chain(neg).max().unwrap();
    assert!(2*worst < PRODUCT, "exact CRT capacity exceeded: {worst}");
    worst
}

pub fn matrix(elements: usize, workers: usize) -> Matrix {
    let mut data = vec![0i16; 4*elements*ROW];
    // Block boundaries coincide with the certified reconstruction boundaries.
    for_blocks(&mut data, 4*ROW, CHUNK, workers, |first, block| {
        for (g, chunk) in block.chunks_mut(CHUNK*4*ROW).enumerate() {
            let start = first+g*CHUNK;
            let here = chunk.len()/(4*ROW);
            let mut raw = vec![[0i16; N]; here];
            let mut coefficients = vec![[0u128; N]; 8];
            for base in (0..here).step_by(8) {
                let count = (here-base).min(8);
                for e in 0..count {
                    crate::rq::ajtai_element(MODULUS,elements,start+base+e,&mut coefficients[e]);
                    // Center residues: map 512..1023 to -512..-1, leave 0..511.
                    for j in 0..N {raw[base+e][j] = (((coefficients[e][j] as i64+512)&1023)-512) as i16;}
                }
                for p in 0..2 { for tail in 0..2 {
                    let prime=primes()[p].p;
                    let mut lanes = vec![[0i16;8];N];
                    for j in 0..N {for e in 0..count {
                        lanes[j][e]=centered(if tail==0 {(coefficients[e][j] % prime as u128) as i64} else {raw[base+e][j] as i64},prime);
                    }}
                    ntt::transform_lanes(p,lanes.as_mut_slice().try_into().unwrap(),count);
                    for e in 0..count {for j in 0..N {
                        chunk[(4*(base+e)+2*p+tail)*ROW+j]=centered(lanes[j][e] as i64,prime);
                    }}
                }}
            }
            certificate(&raw);
        }
    });
    Matrix {elements,data}
}

#[inline]
pub fn crt_pair(a: u32,b: u32) -> i64 {
    // 3889^-1 modulo 9721 = 3242.
    a as i64+3889*((b as i64-a as i64).rem_euclid(9721)*3242%9721)
}
#[inline]
pub fn combine(a: u32,b: u32,tail: i64) -> u128 {
    let x=crt_pair(a,b);
    (x+PRODUCT*((tail-x)*153).rem_euclid(1024)) as u128
}

pub fn commit(values: &[B128], a: &Matrix, workers: usize) -> Vec<u128> {
    #[cfg(target_arch="aarch64")]
    {ntt::neon::hybrid_commit(values,a,workers)}
    #[cfg(not(target_arch="aarch64"))]
    {let _=(values,a,workers); panic!("hybrid commitment requires AArch64")}
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn crt_boundaries() {
        assert_eq!(3889*3242%9721,1); assert_eq!(PRODUCT*153%1024,1);
        for x in [0,1,3888,3889,PRODUCT/2,PRODUCT-1,MODULUS as i64-1] {
            assert_eq!(combine((x%3889) as u32,(x%9721) as u32,x&1023),x as u128);
        }
    }
    #[test] fn certificate_matches_monomials() {
        let a:[i16;N]=core::array::from_fn(|j| ((j*79+17)%1024) as i16-512);
        let mut pos=[0i64;N]; let mut neg=pos;
        for b in 0..512 {
            let mut shifted=[0i64;1160];
            for j in 0..N {shifted[j+b]=a[j] as i64;}
            for j in (N..1160).rev() {let v=shifted[j]; shifted[j-324]+=v; shifted[j-648]-=v;}
            for j in 0..N {pos[j]+=shifted[j].max(0); neg[j]+=(-shifted[j]).max(0);}
        }
        assert_eq!(certificate(&[a]),pos.into_iter().chain(neg).max().unwrap());
    }
    #[test] fn commitment_matches_schoolbook() {
        for elements in [1,3,225] {
            let a=matrix(elements,2);
            let mut s=19u64;
            let values:Vec<B128>=(0..4*elements*9).map(|_| {
                s^=s<<13;s^=s>>7;s^=s<<17;
                B128(s as u128 | ((s.rotate_left(17) as u128)<<64))
            }).collect();
            let raw=crate::rq::ajtai_matrix(MODULUS,elements);
            let want=crate::commit_reference(&values,&raw,MODULUS);
            assert_eq!(commit(&values,&a,2),want,"elements {elements}");
        }
    }
}
