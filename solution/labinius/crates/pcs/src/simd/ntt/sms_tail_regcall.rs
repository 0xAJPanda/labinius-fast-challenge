// Private calling contract for the SMS binary NTT tail, q = 9721.
//
// The straight-line body is the expanded-twiddle `tail27_mem_b` schedule, with
// its modulus reload removed and its 27 Montgomery-root memory operands read
// from zmm31. There is no change to the butterfly arithmetic or output layout.
// Keep the direct comparison with `tail27_mem_b` when editing either schedule.
//
// Inputs: rdi = prefix, rsi = output, rdx/rcx/r8 = expanded twiddle blocks,
// zmm27 = 9721, zmm31 = Montgomery omega. The body writes rax, zmm0..26,
// zmm28..30, and the output buffer. Other GPRs, vector registers, mask
// registers, and flags survive. CALL/RET use the stack, so the call site must
// not declare `nostack`. This is a private assembly contract, not a C ABI call.

core::arch::global_asm!(r#".text
.p2align 4
.globl labinius_sms_tail27_regcall
.hidden labinius_sms_tail27_regcall
.type labinius_sms_tail27_regcall,@function
labinius_sms_tail27_regcall:
lea rax,[rip+{table}+{cv_offset}]
vmovdqa64 zmm0,ZMMWORD PTR [rdi]
vmovdqa64 zmm9,ZMMWORD PTR [rdi+0x240]
vmovdqa64 zmm18,ZMMWORD PTR [rdi+0x480]
vpmullw zmm28,zmm9,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm9,zmm9,ZMMWORD PTR [rdx]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm9,zmm9,zmm28
vpmullw zmm28,zmm18,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm18,zmm18,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm18,zmm18,zmm28
vpsubw zmm28,zmm9,zmm18
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpsrlw zmm30,zmm0,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm0,zmm0,zmm30
vmovdqa64 zmm3,ZMMWORD PTR [rdi+0xc0]
vmovdqa64 zmm12,ZMMWORD PTR [rdi+0x300]
vmovdqa64 zmm21,ZMMWORD PTR [rdi+0x540]
vpmullw zmm29,zmm12,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm12,zmm12,ZMMWORD PTR [rdx]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm12,zmm12,zmm29
vpmullw zmm29,zmm21,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm21,zmm21,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm21,zmm21,zmm29
vpsubw zmm29,zmm12,zmm21
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpsrlw zmm30,zmm3,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm3,zmm3,zmm30
vpaddw zmm30,zmm9,zmm18
vpsubw zmm18,zmm0,zmm18
vpaddw zmm18,zmm18,zmm28
vpsubw zmm9,zmm0,zmm9
vpsubw zmm9,zmm9,zmm28
vpaddw zmm0,zmm0,zmm30
vmovdqa64 zmm6,ZMMWORD PTR [rdi+0x180]
vmovdqa64 zmm15,ZMMWORD PTR [rdi+0x3c0]
vmovdqa64 zmm24,ZMMWORD PTR [rdi+0x600]
vpmullw zmm28,zmm15,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm15,zmm15,ZMMWORD PTR [rdx]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm15,zmm15,zmm28
vpmullw zmm28,zmm24,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm24,zmm24,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm24,zmm24,zmm28
vpsubw zmm28,zmm15,zmm24
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpsrlw zmm30,zmm6,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm6,zmm6,zmm30
vpaddw zmm30,zmm12,zmm21
vpsubw zmm21,zmm3,zmm21
vpaddw zmm21,zmm21,zmm29
vpsubw zmm12,zmm3,zmm12
vpsubw zmm12,zmm12,zmm29
vpaddw zmm3,zmm3,zmm30
vmovdqa64 zmm1,ZMMWORD PTR [rdi+0x40]
vmovdqa64 zmm10,ZMMWORD PTR [rdi+0x280]
vmovdqa64 zmm19,ZMMWORD PTR [rdi+0x4c0]
vpmullw zmm29,zmm10,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm10,zmm10,ZMMWORD PTR [rdx]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm10,zmm10,zmm29
vpmullw zmm29,zmm19,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm19,zmm19,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm19,zmm19,zmm29
vpsubw zmm29,zmm10,zmm19
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpsrlw zmm30,zmm1,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm1,zmm1,zmm30
vpaddw zmm30,zmm15,zmm24
vpsubw zmm24,zmm6,zmm24
vpaddw zmm24,zmm24,zmm28
vpsubw zmm15,zmm6,zmm15
vpsubw zmm15,zmm15,zmm28
vpaddw zmm6,zmm6,zmm30
vmovdqa64 zmm4,ZMMWORD PTR [rdi+0x100]
vmovdqa64 zmm13,ZMMWORD PTR [rdi+0x340]
vmovdqa64 zmm22,ZMMWORD PTR [rdi+0x580]
vpmullw zmm28,zmm13,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm13,zmm13,ZMMWORD PTR [rdx]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm13,zmm13,zmm28
vpmullw zmm28,zmm22,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm22,zmm22,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm22,zmm22,zmm28
vpsubw zmm28,zmm13,zmm22
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpsrlw zmm30,zmm4,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm4,zmm4,zmm30
vpaddw zmm30,zmm10,zmm19
vpsubw zmm19,zmm1,zmm19
vpaddw zmm19,zmm19,zmm29
vpsubw zmm10,zmm1,zmm10
vpsubw zmm10,zmm10,zmm29
vpaddw zmm1,zmm1,zmm30
vmovdqa64 zmm7,ZMMWORD PTR [rdi+0x1c0]
vmovdqa64 zmm16,ZMMWORD PTR [rdi+0x400]
vmovdqa64 zmm25,ZMMWORD PTR [rdi+0x640]
vpmullw zmm29,zmm16,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm16,zmm16,ZMMWORD PTR [rdx]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm16,zmm16,zmm29
vpmullw zmm29,zmm25,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm25,zmm25,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm25,zmm25,zmm29
vpsubw zmm29,zmm16,zmm25
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpsrlw zmm30,zmm7,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm7,zmm7,zmm30
vpaddw zmm30,zmm13,zmm22
vpsubw zmm22,zmm4,zmm22
vpaddw zmm22,zmm22,zmm28
vpsubw zmm13,zmm4,zmm13
vpsubw zmm13,zmm13,zmm28
vpaddw zmm4,zmm4,zmm30
vmovdqa64 zmm2,ZMMWORD PTR [rdi+0x80]
vmovdqa64 zmm11,ZMMWORD PTR [rdi+0x2c0]
vmovdqa64 zmm20,ZMMWORD PTR [rdi+0x500]
vpmullw zmm28,zmm11,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm11,zmm11,ZMMWORD PTR [rdx]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm11,zmm11,zmm28
vpmullw zmm28,zmm20,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm20,zmm20,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm20,zmm20,zmm28
vpsubw zmm28,zmm11,zmm20
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpsrlw zmm30,zmm2,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm2,zmm2,zmm30
vpaddw zmm30,zmm16,zmm25
vpsubw zmm25,zmm7,zmm25
vpaddw zmm25,zmm25,zmm29
vpsubw zmm16,zmm7,zmm16
vpsubw zmm16,zmm16,zmm29
vpaddw zmm7,zmm7,zmm30
vmovdqa64 zmm5,ZMMWORD PTR [rdi+0x140]
vmovdqa64 zmm14,ZMMWORD PTR [rdi+0x380]
vmovdqa64 zmm23,ZMMWORD PTR [rdi+0x5c0]
vpmullw zmm29,zmm14,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm14,zmm14,ZMMWORD PTR [rdx]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm14,zmm14,zmm29
vpmullw zmm29,zmm23,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm23,zmm23,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm23,zmm23,zmm29
vpsubw zmm29,zmm14,zmm23
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpsrlw zmm30,zmm5,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm5,zmm5,zmm30
vpaddw zmm30,zmm11,zmm20
vpsubw zmm20,zmm2,zmm20
vpaddw zmm20,zmm20,zmm28
vpsubw zmm11,zmm2,zmm11
vpsubw zmm11,zmm11,zmm28
vpaddw zmm2,zmm2,zmm30
vmovdqa64 zmm8,ZMMWORD PTR [rdi+0x200]
vmovdqa64 zmm17,ZMMWORD PTR [rdi+0x440]
vmovdqa64 zmm26,ZMMWORD PTR [rdi+0x680]
vpmullw zmm28,zmm17,ZMMWORD PTR [rdx+0x40]
vpmulhw zmm17,zmm17,ZMMWORD PTR [rdx]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm17,zmm17,zmm28
vpmullw zmm28,zmm26,ZMMWORD PTR [rdx+0xc0]
vpmulhw zmm26,zmm26,ZMMWORD PTR [rdx+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm26,zmm26,zmm28
vpsubw zmm28,zmm17,zmm26
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpsrlw zmm30,zmm8,0xb
vpermw zmm30,zmm30,ZMMWORD PTR [rax+0x200]
vpaddw zmm8,zmm8,zmm30
vpaddw zmm30,zmm14,zmm23
vpsubw zmm23,zmm5,zmm23
vpaddw zmm23,zmm23,zmm29
vpsubw zmm14,zmm5,zmm14
vpsubw zmm14,zmm14,zmm29
vpaddw zmm5,zmm5,zmm30
vpaddw zmm30,zmm17,zmm26
vpsubw zmm26,zmm8,zmm26
vpaddw zmm26,zmm26,zmm28
vpsubw zmm17,zmm8,zmm17
vpsubw zmm17,zmm17,zmm28
vpaddw zmm8,zmm8,zmm30
vpmullw zmm28,zmm3,ZMMWORD PTR [rcx+0x40]
vpmulhw zmm3,zmm3,ZMMWORD PTR [rcx]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm3,zmm3,zmm28
vpmullw zmm28,zmm6,ZMMWORD PTR [rcx+0xc0]
vpmulhw zmm6,zmm6,ZMMWORD PTR [rcx+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm6,zmm6,zmm28
vpsubw zmm28,zmm3,zmm6
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpmullw zmm29,zmm21,ZMMWORD PTR [rcx+0x140]
vpmulhw zmm21,zmm21,ZMMWORD PTR [rcx+0x100]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm21,zmm21,zmm29
vpmullw zmm29,zmm24,ZMMWORD PTR [rcx+0x1c0]
vpmulhw zmm24,zmm24,ZMMWORD PTR [rcx+0x180]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm24,zmm24,zmm29
vpsubw zmm29,zmm21,zmm24
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm3,zmm6
vpsubw zmm6,zmm0,zmm6
vpaddw zmm6,zmm6,zmm28
vpsubw zmm3,zmm0,zmm3
vpsubw zmm3,zmm3,zmm28
vpaddw zmm0,zmm0,zmm30
vpmullw zmm28,zmm12,ZMMWORD PTR [rcx+0x240]
vpmulhw zmm12,zmm12,ZMMWORD PTR [rcx+0x200]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm12,zmm12,zmm28
vpmullw zmm28,zmm15,ZMMWORD PTR [rcx+0x2c0]
vpmulhw zmm15,zmm15,ZMMWORD PTR [rcx+0x280]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm15,zmm15,zmm28
vpsubw zmm28,zmm12,zmm15
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm21,zmm24
vpsubw zmm24,zmm18,zmm24
vpaddw zmm24,zmm24,zmm29
vpsubw zmm21,zmm18,zmm21
vpsubw zmm21,zmm21,zmm29
vpaddw zmm18,zmm18,zmm30
vpmullw zmm29,zmm4,ZMMWORD PTR [rcx+0x40]
vpmulhw zmm4,zmm4,ZMMWORD PTR [rcx]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm4,zmm4,zmm29
vpmullw zmm29,zmm7,ZMMWORD PTR [rcx+0xc0]
vpmulhw zmm7,zmm7,ZMMWORD PTR [rcx+0x80]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm7,zmm7,zmm29
vpsubw zmm29,zmm4,zmm7
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm12,zmm15
vpsubw zmm15,zmm9,zmm15
vpaddw zmm15,zmm15,zmm28
vpsubw zmm12,zmm9,zmm12
vpsubw zmm12,zmm12,zmm28
vpaddw zmm9,zmm9,zmm30
vpmullw zmm28,zmm22,ZMMWORD PTR [rcx+0x140]
vpmulhw zmm22,zmm22,ZMMWORD PTR [rcx+0x100]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm22,zmm22,zmm28
vpmullw zmm28,zmm25,ZMMWORD PTR [rcx+0x1c0]
vpmulhw zmm25,zmm25,ZMMWORD PTR [rcx+0x180]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm25,zmm25,zmm28
vpsubw zmm28,zmm22,zmm25
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm4,zmm7
vpsubw zmm7,zmm1,zmm7
vpaddw zmm7,zmm7,zmm29
vpsubw zmm4,zmm1,zmm4
vpsubw zmm4,zmm4,zmm29
vpaddw zmm1,zmm1,zmm30
vpmullw zmm29,zmm13,ZMMWORD PTR [rcx+0x240]
vpmulhw zmm13,zmm13,ZMMWORD PTR [rcx+0x200]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm13,zmm13,zmm29
vpmullw zmm29,zmm16,ZMMWORD PTR [rcx+0x2c0]
vpmulhw zmm16,zmm16,ZMMWORD PTR [rcx+0x280]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm16,zmm16,zmm29
vpsubw zmm29,zmm13,zmm16
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm22,zmm25
vpsubw zmm25,zmm19,zmm25
vpaddw zmm25,zmm25,zmm28
vpsubw zmm22,zmm19,zmm22
vpsubw zmm22,zmm22,zmm28
vpaddw zmm19,zmm19,zmm30
vpmullw zmm28,zmm5,ZMMWORD PTR [rcx+0x40]
vpmulhw zmm5,zmm5,ZMMWORD PTR [rcx]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm5,zmm5,zmm28
vpmullw zmm28,zmm8,ZMMWORD PTR [rcx+0xc0]
vpmulhw zmm8,zmm8,ZMMWORD PTR [rcx+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm8,zmm8,zmm28
vpsubw zmm28,zmm5,zmm8
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm13,zmm16
vpsubw zmm16,zmm10,zmm16
vpaddw zmm16,zmm16,zmm29
vpsubw zmm13,zmm10,zmm13
vpsubw zmm13,zmm13,zmm29
vpaddw zmm10,zmm10,zmm30
vpmullw zmm29,zmm23,ZMMWORD PTR [rcx+0x140]
vpmulhw zmm23,zmm23,ZMMWORD PTR [rcx+0x100]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm23,zmm23,zmm29
vpmullw zmm29,zmm26,ZMMWORD PTR [rcx+0x1c0]
vpmulhw zmm26,zmm26,ZMMWORD PTR [rcx+0x180]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm26,zmm26,zmm29
vpsubw zmm29,zmm23,zmm26
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm5,zmm8
vpsubw zmm8,zmm2,zmm8
vpaddw zmm8,zmm8,zmm28
vpsubw zmm5,zmm2,zmm5
vpsubw zmm5,zmm5,zmm28
vpaddw zmm2,zmm2,zmm30
vpmullw zmm28,zmm14,ZMMWORD PTR [rcx+0x240]
vpmulhw zmm14,zmm14,ZMMWORD PTR [rcx+0x200]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm14,zmm14,zmm28
vpmullw zmm28,zmm17,ZMMWORD PTR [rcx+0x2c0]
vpmulhw zmm17,zmm17,ZMMWORD PTR [rcx+0x280]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm17,zmm17,zmm28
vpsubw zmm28,zmm14,zmm17
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm23,zmm26
vpsubw zmm26,zmm20,zmm26
vpaddw zmm26,zmm26,zmm29
vpsubw zmm23,zmm20,zmm23
vpsubw zmm23,zmm23,zmm29
vpaddw zmm20,zmm20,zmm30
vpaddw zmm30,zmm14,zmm17
vpsubw zmm17,zmm11,zmm17
vpaddw zmm17,zmm17,zmm28
vpsubw zmm14,zmm11,zmm14
vpsubw zmm14,zmm14,zmm28
vpaddw zmm11,zmm11,zmm30
vpmullw zmm28,zmm1,ZMMWORD PTR [r8+0x40]
vpmulhrsw zmm30,zmm0,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm0,zmm0,zmm30
vpmulhw zmm1,zmm1,ZMMWORD PTR [r8]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm1,zmm1,zmm28
vpmullw zmm28,zmm2,ZMMWORD PTR [r8+0xc0]
vpmulhw zmm2,zmm2,ZMMWORD PTR [r8+0x80]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm2,zmm2,zmm28
vpsubw zmm28,zmm1,zmm2
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpmullw zmm29,zmm7,ZMMWORD PTR [r8+0x140]
vpmulhrsw zmm30,zmm6,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm6,zmm6,zmm30
vpmulhw zmm7,zmm7,ZMMWORD PTR [r8+0x100]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm7,zmm7,zmm29
vpmullw zmm29,zmm8,ZMMWORD PTR [r8+0x1c0]
vpmulhw zmm8,zmm8,ZMMWORD PTR [r8+0x180]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm8,zmm8,zmm29
vpsubw zmm29,zmm7,zmm8
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm1,zmm2
vpsubw zmm2,zmm0,zmm2
vpaddw zmm2,zmm2,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x40],zmm2
vpsubw zmm1,zmm0,zmm1
vpsubw zmm1,zmm1,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x80],zmm1
vpaddw zmm0,zmm0,zmm30
vmovdqa64 ZMMWORD PTR [rsi],zmm0
vpmullw zmm28,zmm4,ZMMWORD PTR [r8+0x240]
vpmulhrsw zmm30,zmm3,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm3,zmm3,zmm30
vpmulhw zmm4,zmm4,ZMMWORD PTR [r8+0x200]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm4,zmm4,zmm28
vpmullw zmm28,zmm5,ZMMWORD PTR [r8+0x2c0]
vpmulhw zmm5,zmm5,ZMMWORD PTR [r8+0x280]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm5,zmm5,zmm28
vpsubw zmm28,zmm4,zmm5
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm7,zmm8
vpsubw zmm8,zmm6,zmm8
vpaddw zmm8,zmm8,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x100],zmm8
vpsubw zmm7,zmm6,zmm7
vpsubw zmm7,zmm7,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x140],zmm7
vpaddw zmm6,zmm6,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0xc0],zmm6
vpmullw zmm29,zmm19,ZMMWORD PTR [r8+0x340]
vpmulhrsw zmm30,zmm18,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm18,zmm18,zmm30
vpmulhw zmm19,zmm19,ZMMWORD PTR [r8+0x300]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm19,zmm19,zmm29
vpmullw zmm29,zmm20,ZMMWORD PTR [r8+0x3c0]
vpmulhw zmm20,zmm20,ZMMWORD PTR [r8+0x380]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm20,zmm20,zmm29
vpsubw zmm29,zmm19,zmm20
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm4,zmm5
vpsubw zmm5,zmm3,zmm5
vpaddw zmm5,zmm5,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x1c0],zmm5
vpsubw zmm4,zmm3,zmm4
vpsubw zmm4,zmm4,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x200],zmm4
vpaddw zmm3,zmm3,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0x180],zmm3
vpmullw zmm28,zmm25,ZMMWORD PTR [r8+0x440]
vpmulhrsw zmm30,zmm24,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm24,zmm24,zmm30
vpmulhw zmm25,zmm25,ZMMWORD PTR [r8+0x400]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm25,zmm25,zmm28
vpmullw zmm28,zmm26,ZMMWORD PTR [r8+0x4c0]
vpmulhw zmm26,zmm26,ZMMWORD PTR [r8+0x480]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm26,zmm26,zmm28
vpsubw zmm28,zmm25,zmm26
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm19,zmm20
vpsubw zmm20,zmm18,zmm20
vpaddw zmm20,zmm20,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x280],zmm20
vpsubw zmm19,zmm18,zmm19
vpsubw zmm19,zmm19,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x2c0],zmm19
vpaddw zmm18,zmm18,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0x240],zmm18
vpmullw zmm29,zmm22,ZMMWORD PTR [r8+0x540]
vpmulhrsw zmm30,zmm21,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm21,zmm21,zmm30
vpmulhw zmm22,zmm22,ZMMWORD PTR [r8+0x500]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm22,zmm22,zmm29
vpmullw zmm29,zmm23,ZMMWORD PTR [r8+0x5c0]
vpmulhw zmm23,zmm23,ZMMWORD PTR [r8+0x580]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm23,zmm23,zmm29
vpsubw zmm29,zmm22,zmm23
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm25,zmm26
vpsubw zmm26,zmm24,zmm26
vpaddw zmm26,zmm26,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x340],zmm26
vpsubw zmm25,zmm24,zmm25
vpsubw zmm25,zmm25,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x380],zmm25
vpaddw zmm24,zmm24,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0x300],zmm24
vpmullw zmm28,zmm10,ZMMWORD PTR [r8+0x640]
vpmulhrsw zmm30,zmm9,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm9,zmm9,zmm30
vpmulhw zmm10,zmm10,ZMMWORD PTR [r8+0x600]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm10,zmm10,zmm28
vpmullw zmm28,zmm11,ZMMWORD PTR [r8+0x6c0]
vpmulhw zmm11,zmm11,ZMMWORD PTR [r8+0x680]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm11,zmm11,zmm28
vpsubw zmm28,zmm10,zmm11
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm22,zmm23
vpsubw zmm23,zmm21,zmm23
vpaddw zmm23,zmm23,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x400],zmm23
vpsubw zmm22,zmm21,zmm22
vpsubw zmm22,zmm22,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x440],zmm22
vpaddw zmm21,zmm21,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0x3c0],zmm21
vpmullw zmm29,zmm16,ZMMWORD PTR [r8+0x740]
vpmulhrsw zmm30,zmm15,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm15,zmm15,zmm30
vpmulhw zmm16,zmm16,ZMMWORD PTR [r8+0x700]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm16,zmm16,zmm29
vpmullw zmm29,zmm17,ZMMWORD PTR [r8+0x7c0]
vpmulhw zmm17,zmm17,ZMMWORD PTR [r8+0x780]
vpmulhw zmm29,zmm29,zmm27
vpsubw zmm17,zmm17,zmm29
vpsubw zmm29,zmm16,zmm17
vpmullw zmm30,zmm29,ZMMWORD PTR [rax+0x80]
vpmulhw zmm29,zmm29,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm29,zmm29,zmm30
vpaddw zmm30,zmm10,zmm11
vpsubw zmm11,zmm9,zmm11
vpaddw zmm11,zmm11,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x4c0],zmm11
vpsubw zmm10,zmm9,zmm10
vpsubw zmm10,zmm10,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x500],zmm10
vpaddw zmm9,zmm9,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0x480],zmm9
vpmullw zmm28,zmm13,ZMMWORD PTR [r8+0x840]
vpmulhrsw zmm30,zmm12,ZMMWORD PTR [rax+0xc0]
vpmullw zmm30,zmm30,zmm27
vpsubw zmm12,zmm12,zmm30
vpmulhw zmm13,zmm13,ZMMWORD PTR [r8+0x800]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm13,zmm13,zmm28
vpmullw zmm28,zmm14,ZMMWORD PTR [r8+0x8c0]
vpmulhw zmm14,zmm14,ZMMWORD PTR [r8+0x880]
vpmulhw zmm28,zmm28,zmm27
vpsubw zmm14,zmm14,zmm28
vpsubw zmm28,zmm13,zmm14
vpmullw zmm30,zmm28,ZMMWORD PTR [rax+0x80]
vpmulhw zmm28,zmm28,zmm31
vpmulhw zmm30,zmm30,zmm27
vpsubw zmm28,zmm28,zmm30
vpaddw zmm30,zmm16,zmm17
vpsubw zmm17,zmm15,zmm17
vpaddw zmm17,zmm17,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x580],zmm17
vpsubw zmm16,zmm15,zmm16
vpsubw zmm16,zmm16,zmm29
vmovdqa64 ZMMWORD PTR [rsi+0x5c0],zmm16
vpaddw zmm15,zmm15,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0x540],zmm15
vpaddw zmm30,zmm13,zmm14
vpsubw zmm14,zmm12,zmm14
vpaddw zmm14,zmm14,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x640],zmm14
vpsubw zmm13,zmm12,zmm13
vpsubw zmm13,zmm13,zmm28
vmovdqa64 ZMMWORD PTR [rsi+0x680],zmm13
vpaddw zmm12,zmm12,zmm30
vmovdqa64 ZMMWORD PTR [rsi+0x600],zmm12
ret
.size labinius_sms_tail27_regcall,.-labinius_sms_tail27_regcall"#,table=sym T9721,cv_offset=const core::mem::offset_of!(Tables,cv));
unsafe extern "C"{fn labinius_sms_tail27_regcall();}

#[inline(always)]
unsafe fn tail27_sms_regcall(
    bp: *const i16, op: *mut i16, t4: *const u32, t5: *const u32, t6: *const u32,
) {
    core::arch::asm!(
        "call {tail}",
        tail = sym labinius_sms_tail27_regcall,
        in("rdi") bp, in("rsi") op, in("rdx") t4, in("rcx") t5, in("r8") t6,
        in("zmm27") _mm512_set1_epi16(9721),
        in("zmm31") _mm512_set1_epi16(Params::<9721>::to_mont(Params::<9721>::OMEGA)),
        out("rax") _,
        out("zmm0") _,out("zmm1") _,out("zmm2") _,out("zmm3") _,out("zmm4") _,out("zmm5") _,out("zmm6") _,out("zmm7") _,out("zmm8") _,out("zmm9") _,out("zmm10") _,out("zmm11") _,out("zmm12") _,out("zmm13") _,out("zmm14") _,out("zmm15") _,out("zmm16") _,out("zmm17") _,out("zmm18") _,out("zmm19") _,out("zmm20") _,out("zmm21") _,out("zmm22") _,out("zmm23") _,out("zmm24") _,out("zmm25") _,out("zmm26") _,out("zmm28") _,out("zmm29") _,out("zmm30") _,
        options(preserves_flags),
    );
}

#[cfg(test)]
mod sms_regcall_tests {
    use super::*;
    #[repr(C, align(64))]
    struct Block([i16; 27 * 32]);

    #[test]
    fn register_call_matches_original_leaf() {
        let mut state = 0x8fa6_e352_4c91_d7b0u64;
        let mut input = Block([0; 27 * 32]);
        let mut expected = Block([0; 27 * 32]);
        let mut got = Block([0; 27 * 32]);
        let bound = 3 * 9721;
        for fixture in 0..6 {
            for (i, x) in input.0.iter_mut().enumerate() {
                state ^= state << 13; state ^= state >> 7; state ^= state << 17;
                *x = match fixture {
                    0 => 0,
                    1 => bound as i16,
                    2 => -bound as i16,
                    3 => if i % 2 == 0 { bound as i16 } else { -bound as i16 },
                    _ => ((state % (2 * bound + 1) as u64) as i32 - bound) as i16,
                };
            }
            for block in 0..24 {
                unsafe {
                    let t4 = SMS_TWIDDLES.t4.as_ptr().add(4 * block).cast::<u32>();
                    let t5 = SMS_TWIDDLES.t5.as_ptr().add(12 * block).cast::<u32>();
                    let t6 = SMS_TWIDDLES.t6.as_ptr().add(36 * block).cast::<u32>();
                    let cv = T9721.cv.as_ptr().cast::<i16>();
                    tail27_mem_b(input.0.as_ptr(), expected.0.as_mut_ptr(), t4, t5, t6, cv);
                    tail27_sms_regcall(input.0.as_ptr(), got.0.as_mut_ptr(), t4, t5, t6);
                }
                assert_eq!(got.0, expected.0, "fixture {fixture}, block {block}");
            }
        }
    }
}
