# Controller list (one RunHeader per acquisition device) in both files; identify the MS controller.
using Mmap
files = [(expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw"), 2866, [1320253552]),
         (expanduser("~/Projects/pioneer-gui-test/rawtest/HELA_MOUSESIL_METHTEST_07312024_06.raw"), 4264, [463444680, 505572910, 505935206, 638409553])]
for (path, first_ptr, rhs) in files
    b = Mmap.mmap(path); rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
    println("== ", basename(path))
    for o in first_ptr-32:4:first_ptr+16*length(rhs)+8
        println("  @$o (ptr$(o >= first_ptr ? "+" : "")$(o - first_ptr)): u32 ", rd(UInt32, o), "  i32 ", rd(Int32, o))
    end
    for rh in rhs
        println("  RunHeader @$rh: +0..+16 Int32 ", [rd(Int32, rh + o) for o in 0:4:16],
                "; index/data/events addr ", [rd(UInt64, rh + o) for o in (7408, 7416, 7448)])
    end
end
