//! SVE halfword dot products, four ring elements per i64 accumulator lane.
use super::*;
const _: () = {
    let bound=super::super::transform_bound(PRIMES[1],SCHEDULES[1]) as u64;
    assert!(8*bound*4860 < i32::MAX as u64);
    assert!(4096*bound*4860 < (1u64<<41));
};
#[inline(always)]
unsafe fn pack(a:int16x8_t,b:int16x8_t,c:int16x8_t,d:int16x8_t)->[int16x8_t;4] {
    let ab0=vreinterpretq_s32_s16(vzip1q_s16(a,b));
    let ab1=vreinterpretq_s32_s16(vzip2q_s16(a,b));
    let cd0=vreinterpretq_s32_s16(vzip1q_s16(c,d));
    let cd1=vreinterpretq_s32_s16(vzip2q_s16(c,d));
    [vreinterpretq_s16_s32(vzip1q_s32(ab0,cd0)),vreinterpretq_s16_s32(vzip2q_s32(ab0,cd0)),
     vreinterpretq_s16_s32(vzip1q_s32(ab1,cd1)),vreinterpretq_s16_s32(vzip2q_s32(ab1,cd1))]
}
#[inline(always)]
unsafe fn dot(acc:&mut [int64x2_t],a:[int16x8_t;4],b:[int16x8_t;4],key:int16x8_t) {
    let key=vnegq_s16(key);
    let k0=vreinterpretq_s16_s64(vdupq_laneq_s64::<0>(vreinterpretq_s64_s16(key)));
    let k1=vreinterpretq_s16_s64(vdupq_laneq_s64::<1>(vreinterpretq_s64_s16(key)));
    let(mut s0,mut s1,mut s2,mut s3)=(acc[0],acc[1],acc[2],acc[3]);
    #[cfg(target_feature="sve")]
    core::arch::asm!(
        "sdot z0.d, z4.h, z12.h", "sdot z1.d, z5.h, z12.h",
        "sdot z2.d, z6.h, z12.h", "sdot z3.d, z7.h, z12.h",
        "sdot z0.d, z8.h, z13.h", "sdot z1.d, z9.h, z13.h",
        "sdot z2.d, z10.h, z13.h", "sdot z3.d, z11.h, z13.h",
        inout("v0") s0, inout("v1") s1, inout("v2") s2, inout("v3") s3,
        in("v4") a[0], in("v5") a[1], in("v6") a[2], in("v7") a[3],
        in("v8") b[0], in("v9") b[1], in("v10") b[2], in("v11") b[3],
        in("v12") k0, in("v13") k1,
        options(pure,nomem,nostack,preserves_flags)
    );
    #[cfg(not(target_feature="sve"))]
    {
        let term=|x:int16x8_t,k:int16x8_t| {
            let lo=vmull_s16(vget_low_s16(x),vget_low_s16(k));
            let hi=vmull_high_s16(x,k);
            vpaddlq_s32(vpaddq_s32(lo,hi))
        };
        s0=vaddq_s64(s0,vaddq_s64(term(a[0],k0),term(b[0],k1)));
        s1=vaddq_s64(s1,vaddq_s64(term(a[1],k0),term(b[1],k1)));
        s2=vaddq_s64(s2,vaddq_s64(term(a[2],k0),term(b[2],k1)));
        s3=vaddq_s64(s3,vaddq_s64(term(a[3],k0),term(b[3],k1)));
    }
    acc[..4].copy_from_slice(&[s0,s1,s2,s3]);
}
pub(super) unsafe fn tile(chunk:&crate::sms::Chunk,batch:usize,x:&[[int16x8_t;N]],acc:&mut [int64x2_t]) {
    const S:usize=crate::sms::STREAMS;
    for j in 0..N {
        let a=pack(x[0][j],x[1][j],x[2][j],x[3][j]);
        let b=pack(x[4][j],x[5][j],x[6][j],x[7][j]);
        for s in 0..S {dot(&mut acc[4*(j*S+s)..],a,b,vld1q_s16(chunk.key(batch,j,s).as_ptr()));}
    }
}
/// The largest full main-stream sum is below 2^41. Double precision gives
/// enough quotient accuracy that rounding has error < 1/(2p); p is odd.
#[inline(always)]
unsafe fn narrow(v:int64x2_t)->int32x2_t {
    let q=vcvtnq_s64_f64(vmulq_n_f64(vcvtq_f64_s64(v),1.0/9721.0));
    vmovn_s64(vsubq_s64(v,vmull_n_s32(vmovn_s64(q),9721)))
}
pub(super) unsafe fn inverse(stream:usize,acc:&[int64x2_t],x:&mut [int16x8_t;N]) {
    for j in 0..N {
        let at=4*(j*crate::sms::STREAMS+stream);
        let lo=vcombine_s32(narrow(acc[at]),narrow(acc[at+1]));
        let hi=vcombine_s32(narrow(acc[at+2]),narrow(acc[at+3]));
        x[j]=vcombine_s16(vmovn_s32(lo),vmovn_s32(hi));
    }
    inverse_lanes(1,x);
}
#[cfg(test)] mod tests {
    use super::*;
    #[test] fn kernels() {unsafe {
        use std::hint::black_box;
        let s=crate::sms::STREAMS;
        let key:Vec<i16>=(0..N*s*8).map(|i|((i*1337+7)%9721)as i16-4860).collect();
        let chunk=crate::sms::Chunk{elements:8,data:key.into(),sums:vec![].into(),margin:0.0};
        let x:Vec<[int16x8_t;N]>=(0..8).map(|e|core::array::from_fn(|j| {
            let v:[i16;8]=core::array::from_fn(|c|((j*199+c*117+e*31)%9721)as i16-4860);vld1q_s16(v.as_ptr())
        })).collect();
        let mut a=vec![vdupq_n_s32(0);2*N*s];let mut b=vec![vdupq_n_s64(0);4*N*s];
        sms_tile(&chunk,0,&x,&mut a);tile(&chunk,0,&x,&mut b);
        for j in 0..N*s {let av:[i32;8]=core::mem::transmute([a[2*j],a[2*j+1]]);let bv:[i64;8]=core::mem::transmute([b[4*j],b[4*j+1],b[4*j+2],b[4*j+3]]);for c in 0..8 {assert_eq!(av[c] as i64,bv[c]);}}
        let t=std::time::Instant::now();for _ in 0..1000 {sms_tile(black_box(&chunk),0,black_box(&x),black_box(&mut a));}let neon=t.elapsed();
        let t=std::time::Instant::now();for _ in 0..1000 {tile(black_box(&chunk),0,black_box(&x),black_box(&mut b));}let sve=t.elapsed();
        eprintln!("MAC 1000 tiles: smlsl={neon:?} sve-sdot={sve:?}");
        for v in [-((1i64)<<41),-9722,-4861,-4860,-1,0,4860,4861,9722,1<<41] {
            let r:[i32;2]=core::mem::transmute(narrow(vdupq_n_s64(v)));assert_eq!((r[0] as i64-v).rem_euclid(9721),0);assert!((r[0] as i64).abs()<=4860);
        }
    }}
}
