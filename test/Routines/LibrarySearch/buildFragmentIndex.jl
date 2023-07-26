function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end

@testset "buildFragmentIndex.jl" begin

frag_ions = [
    #Same fragment bin. Split into two RT bins
    FragmentIon(100.0, UInt32(1), 200.0, Ref(1e6), 10.0, UInt8(2)),
    FragmentIon(100.0, UInt32(1), 201.0, Ref(2e6), 10.1, UInt8(2)),
    FragmentIon(100.0, UInt32(1), 203.0, Ref(3e6), 10.1, UInt8(2)),
    FragmentIon(100.0, UInt32(1), 204.0, Ref(4e6), 100.0, UInt8(2)),
    FragmentIon(100.0, UInt32(1), 205.0, Ref(5e6), 100.0, UInt8(2)),
    FragmentIon(100.01, UInt32(1), 206.0, Ref(6e6), 102.0, UInt8(2)),

    #New framgent bin but otherwise the same
    FragmentIon(200.0, UInt32(1), 207.0, Ref(7e6), 10.0, UInt8(2)),
    FragmentIon(200.0, UInt32(1), 208.0, Ref(8e6), 10.1, UInt8(2)),
    FragmentIon(200.0, UInt32(1), 209.0, Ref(9e6), 10.1, UInt8(2)),
    FragmentIon(200.0, UInt32(1), 210.0, Ref(10e6), 100.0, UInt8(2)),
    FragmentIon(200.0, UInt32(1), 211.0, Ref(11e6), 100.0, UInt8(2)),
    FragmentIon(200.01, UInt32(1), 212.0, Ref(12e6), 102.0, UInt8(2)),


]
#Build a toy framgnent index. 
f_index = buildFragmentIndex!(frag_ions, 20.0, 5.0, low_frag_mz = 50.0, low_prec_mz = 50.0);

@test length(f_index.fragment_bins) == 4
@test [length(x) for x in f_index.rt_bins] ==[2, 1, 2, 1]
@test [length(x.precs) for x in f_index.precursor_bins] ==[3, 2, 1, 3, 2, 1]
@test Tol(getLowMZ(f_index.fragment_bins[1]), 100.0)
@test Tol(getLowMZ(f_index.fragment_bins[2]), 100.01)
@test Tol(getLowMZ(f_index.fragment_bins[3]), 200.0)
@test Tol(getLowMZ(f_index.fragment_bins[4]), 200.01)

#Make sure precursor bins are sorted
@test all([issorted(prec_bin.precs, by = x->getPrecMZ(x)) for prec_bin in f_index.precursor_bins])
#Reverse frag ions and try against
frag_ions = frag_ions[[x for x in length(frag_ions):-1:1]]

f_index = buildFragmentIndex!(frag_ions, 20.0, 5.0, low_frag_mz = 50.0, low_prec_mz = 50.0);

@test length(f_index.fragment_bins) == 4
@test [length(x) for x in f_index.rt_bins] ==[2, 1, 2, 1]
@test [length(x.precs) for x in f_index.precursor_bins] ==[3, 2, 1, 3, 2, 1]
@test Tol(getLowMZ(f_index.fragment_bins[1]), 100.0)
@test Tol(getLowMZ(f_index.fragment_bins[2]), 100.01)
@test Tol(getLowMZ(f_index.fragment_bins[3]), 200.0)
@test Tol(getLowMZ(f_index.fragment_bins[4]), 200.01)

@test all([issorted(prec_bin.precs, by = x->getPrecMZ(x)) for prec_bin in f_index.precursor_bins])

#Test low and high frag mz 
f_index = buildFragmentIndex!(frag_ions, 20.0, 5.0, low_frag_mz = 150.0, low_prec_mz = 50.0);
@test length(f_index.precursor_bins) == 3
@test Tol(f_index.precursor_bins[1].precs[1].prec_mz, 207.0)

f_index = buildFragmentIndex!(frag_ions, 20.0, 5.0, low_frag_mz = 50.0, high_frag_mz = 150.0, low_prec_mz = 50.0);
@test length(f_index.precursor_bins) == 3
@test Tol(f_index.precursor_bins[1].precs[1].prec_mz, 200.0)
@test Tol(f_index.precursor_bins[3].precs[1].prec_mz, 206.0)




end