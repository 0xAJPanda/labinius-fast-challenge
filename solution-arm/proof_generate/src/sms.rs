//! One binary NTT with balanced radix-sixteen public-key digits.
//! Reconstruction about the public mean has a certified Hoeffding margin;
//! unlike hybrid.rs, this is a statistical guarantee for fresh uniform tables.
use crate::ntt::{self, N, centered, for_blocks};
use labinius::fields::scalar::B128;
pub const P:i64=9721;
pub const TAIL:i64=1<<22;
pub const MODULUS:u128=P as u128*TAIL as u128;
pub const DIGITS:usize=6;
pub const STREAMS:usize=DIGITS+1;
pub const CHUNK:usize=64;
pub struct Chunk {pub elements:usize,pub data:std::borrow::Cow<'static,[i16]>,pub means:std::borrow::Cow<'static,[i32]>,pub margin:f64}
impl Chunk {pub fn key(&self,b:usize,j:usize,s:usize)->&[i16] {&self.data[((b*N+j)*STREAMS+s)*8..][..8]}}
pub struct Matrix {pub elements:usize,pub chunks:Vec<Chunk>}
pub fn digits(mut x:i64)->[i8;DIGITS] {
    let mut d=[0i8;DIGITS];
    x &= TAIL-1;
    for a in &mut d[..5] {*a=((x+8)&15) as i8-8;x=(x+8)>>4;}
    d[5]=(x&3) as i8;d
}
/// Exact first and second moments of the coefficient's binary linear form.
fn moments(l:&[i8;N],s1:&mut [i64;N],s2:&mut [i64;N]) {
    let at=|x:i64|if (0..N as i64).contains(&x){l[x as usize] as i64}else{0};
    const W:usize=512+N-1;
    let mut h=[0i64;W+1];let mut h2=h;let mut g=h;let mut g2=h;
    for k in 0..W {
        let x=k as i64-511;let a=at(x)+at(x+324);let b=at(x)-at(x+648)-at(x+972);
        h[k+1]=h[k]+a;h2[k+1]=h2[k]+a*a;g[k+1]=g[k]+b;g2[k+1]=g2[k]+b*b;
    }
    for j in 0..N {
        if j>=324 {s1[j]+=h[j+512]-h[j];s2[j]+=h2[j+512]-h2[j];}
        else {s1[j]+=g[j+512]-g[j];s2[j]+=g2[j+512]-g2[j];}
    }
}
pub fn matrix(elements:usize,workers:usize)->Matrix {
    let mut chunks:Vec<Chunk>=(0..elements.div_ceil(CHUNK)).map(|c|Chunk {elements:(elements-c*CHUNK).min(CHUNK),data:Vec::new().into(),means:Vec::new().into(),margin:f64::INFINITY}).collect();
    for_blocks(&mut chunks,1,1,workers,|first,blocks| {
        for (ci,chunk) in blocks.iter_mut().enumerate() {
            let start=(first+ci)*CHUNK;
            let mut data=vec![0i16;chunk.elements.div_ceil(8)*N*STREAMS*8];
            let mut s1=vec![[0i64;N];DIGITS];let mut s2=s1.clone();
            let mut coeff=vec![[0u128;N];8];let mut ds=vec![[[0i8;DIGITS];N];8];
            for base in (0..chunk.elements).step_by(8) {
                let count=(chunk.elements-base).min(8);
                for e in 0..count {
                    crate::rq::ajtai_element(MODULUS,elements,start+base+e,&mut coeff[e]);
                    for j in 0..N {ds[e][j]=digits(coeff[e][j] as i64);}
                    for d in 0..DIGITS {
                        let l=core::array::from_fn(|j|ds[e][j][d]);moments(&l,&mut s1[d],&mut s2[d]);
                    }
                }
                for stream in 0..STREAMS {
                    let mut lanes=vec![[0i16;8];N];
                    for j in 0..N {for e in 0..count {
                        lanes[j][e]=if stream==0 {centered((coeff[e][j]%P as u128) as i64,P as u32)} else {ds[e][j][stream-1] as i16};
                    }}
                    ntt::transform_lanes(1,lanes.as_mut_slice().try_into().unwrap(),count);
                    for e in 0..count {for j in 0..N {data[((base/8*N+j)*STREAMS+stream)*8+e]=centered(lanes[j][e] as i64,P as u32);}}
                }
            }
            chunk.data=data.into();
            let mut means=vec![0i32;DIGITS*N];
            for d in 0..DIGITS {for j in 0..N {
                let mean=s1[d][j].div_euclid(2);means[d*N+j]=mean as i32;
                if s2[d][j]!=0 {chunk.margin=chunk.margin.min(4859.0/(s2[d][j] as f64/4.0).sqrt());}
            }}
            chunk.means=means.into();
            assert!(chunk.margin>=8.5,"insufficient reconstruction margin: {}",chunk.margin);
        }
    });
    eprintln!("sms digits={DIGITS} chunk={CHUNK} min_margin={:.3}",chunks.iter().map(|c|c.margin).fold(f64::INFINITY,f64::min));
    Matrix {elements,chunks}
}
/// [`matrix`] over `elements` elements from its precomputed image (`build.rs`): the chunks' key
/// data and means concatenated, and their margins; `None` if the image has another shape.
pub fn from_image(data:&'static [i16],means:&'static [i32],margins:&'static [f64],elements:usize)->Option<Matrix> {
    let count=elements.div_ceil(CHUNK);
    let sizes:Vec<usize>=(0..count).map(|c|(elements-c*CHUNK).min(CHUNK)).collect();
    let total:usize=sizes.iter().map(|e|e.div_ceil(8)*N*STREAMS*8).sum();
    if data.len()!=total || means.len()!=count*DIGITS*N || margins.len()!=count {return None;}
    let (mut d,mut m)=(data,means);
    let chunks=sizes.iter().zip(margins).map(|(&e,&margin)| {
        let (a,rest)=d.split_at(e.div_ceil(8)*N*STREAMS*8);d=rest;
        let (b,rest)=m.split_at(DIGITS*N);m=rest;
        assert!(margin>=8.5,"insufficient reconstruction margin: {margin}");
        Chunk {elements:e,data:a.into(),means:b.into(),margin}
    }).collect();
    Some(Matrix {elements,chunks})
}
#[inline]
pub fn combine(base:i64,tail:i64)->u128 {
    const INV:i64=980; // (2^22)^-1 mod 9721, checked below.
    (tail+TAIL*((base-tail).rem_euclid(P)*INV%P)) as u128
}
pub fn commit(values:&[B128],a:&Matrix,workers:usize)->Vec<u128> {
    #[cfg(target_arch="aarch64")] {ntt::neon::sms_commit(values,a,workers)}
    #[cfg(not(target_arch="aarch64"))] {let _=(values,a,workers);panic!("AArch64 required")}
}
#[cfg(test)] mod tests {
    use super::*;
    #[test] fn digit_encoding() {
        let mut s=7;
        for _ in 0..100000 {
            s=ntt::tests::xorshift(&mut s);let x=s as i64&(TAIL-1);let ds=digits(x);
            let y=ds.iter().enumerate().map(|(d,&v)|(v as i64)<<(4*d)).sum::<i64>();
            assert_eq!((y-x).rem_euclid(TAIL),0);assert!(ds[..5].iter().all(|&x|(-8..8).contains(&x)));assert!((0..4).contains(&ds[5]));
            let z=(s as u128%MODULUS) as i64;assert_eq!(combine(z%P,z&(TAIL-1)),z as u128);
        }
    }
    #[test] fn moment_reference() {
        let a:[i8;N]=core::array::from_fn(|j|(j%8) as i8-4);
        let (mut s1,mut s2)=([0i64;N],[0i64;N]);moments(&a,&mut s1,&mut s2);
        let (mut r1,mut r2)=([0i64;N],[0i64;N]);
        for b in 0..512 {let mut x=[0i64;1160];for j in 0..N {x[j+b]=a[j] as i64;}
            for j in (N..1160).rev() {let v=x[j];x[j-324]+=v;x[j-648]-=v;}
            for j in 0..N {r1[j]+=x[j];r2[j]+=x[j]*x[j];}
        }assert_eq!(s1,r1);assert_eq!(s2,r2);
    }
    #[test] fn commitment_reference() {
        for elements in [1,3,129] {
            let a=matrix(elements,2);let mut s=17;
            let values:Vec<B128>=(0..4*elements*9).map(|_| {
                let lo=ntt::tests::xorshift(&mut s);let hi=ntt::tests::xorshift(&mut s);B128(lo as u128|((hi as u128)<<64))
            }).collect();
            let raw=crate::rq::ajtai_matrix(MODULUS,elements);
            assert_eq!(commit(&values,&a,2),crate::commit_reference(&values,&raw,MODULUS));
        }
    }
}
