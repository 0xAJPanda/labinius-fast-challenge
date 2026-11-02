#!/usr/bin/env python3
# Regenerates `sms_slot8` in src/simd/sms_slot_gen.rs: python3 tools/gen_sms_slot8.py 7
# For 8 streams it reproduces the original eight-stream macro byte for byte.
"""Generate the sms_slot8 macro (LOCK = 8) for a given number of accumulation streams."""
import sys

def gen(streams: int) -> str:
    L = 8
    out = ["macro_rules! sms_slot8 {"]
    def loads(first):
        r = ['        "vmovdqa64 zmm20, zmmword ptr [{a} + 0]\\n",']
        for s in range(1, streams):
            if first:
                r.append(f'        "vpbroadcastd zmm{20+s}, dword ptr [{{b}}]\\n",')
            else:
                r.append(f'        "vmovdqa64 zmm{20+s}, zmmword ptr [{{a}} + {64*s}]\\n",')
        return r
    def macs():
        r = []
        for h in range(L):
            for s in range(streams):
                r.append(f'        "vpdpwssd zmm{20+s}, {{x{h}}}, dword ptr [{{k}} + {4*(L*s+h)}]{{{{1to16}}}}\\n",')
        return r
    def stores(red):
        r = []
        for s in range(streams):
            if red:
                r.append(f'        "vpxord zmm29, zmm{20+s}, dword ptr [{{f}} + 64]{{{{1to16}}}}\\n",')
                r.append('        "vpmaddwd zmm29, zmm29, zmmword ptr [{f}]\\n",')
                r.append(f'        "vmovdqa64 zmmword ptr [{{a}} + {64*s}], zmm29\\n",')
            else:
                r.append(f'        "vmovdqa64 zmmword ptr [{{a}} + {64*s}], zmm{20+s}\\n",')
        return r
    for name, first, red in [("plain", False, False), ("reduce", False, True), ("reduce_first", True, True), ("plain_first", True, False)]:
        out.append(f"    ({name}) => {{ concat!(")
        out += loads(first) + macs() + stores(red)
        out.append("    ) };")
    out.append("}")
    return "\n".join(out) + "\n"

if __name__ == "__main__":
    sys.stdout.write(gen(int(sys.argv[1])))
