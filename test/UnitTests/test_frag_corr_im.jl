# n_scans_in_window of _add_fragment_chromatogram_features! on ion-mobility slice data (timsTOF):
# PSMs of one precursor spread over mobility slices and cycles.
using Test
using DataFrames
using Pioneer

@testset "n_scans_in_window on mobility slices" begin
    # Precursor 5: cycle 1 one slice, cycle 2 two slices, cycles 3-4 one slice; precursor 9: a singleton.
    psms = DataFrame(precursor_idx = UInt32[5, 5, 5, 5, 5, 9], scan_idx = UInt32[1, 2, 3, 4, 5, 6],
                     cycle_idx = Int32[1, 2, 2, 3, 4, 7], weight = Float32[1, 4, 3, 2, 1, 5],
                     irt_obs = Float32[1, 2, 2, 3, 4, 7])
    for r in 1:8; psms[!, Symbol("frag$(r)_int")] = Float32[1, 4, 3, 2, 1, 5] .* Float32(r); end
    Pioneer._add_fragment_chromatogram_features!(psms)
    @test psms.n_scans_in_window == UInt32[1, 2, 2, 1, 1, 1]
    @test psms.n_scans == UInt32[5, 5, 5, 5, 5, 1]
    # Three slices in one cycle.
    psms2 = DataFrame(precursor_idx = fill(UInt32(3), 4), scan_idx = UInt32.(1:4), cycle_idx = Int32[1, 1, 1, 2],
                      weight = Float32[1, 5, 1, 2], irt_obs = Float32[1, 1, 1, 2])
    for r in 1:8; psms2[!, Symbol("frag$(r)_int")] = Float32[1, 5, 1, 2] .* Float32(r); end
    Pioneer._add_fragment_chromatogram_features!(psms2)
    @test psms2.n_scans_in_window == UInt32[3, 3, 3, 1]
end

@testset "n_scans_in_window without a cycle column stays 1" begin
    psms = DataFrame(precursor_idx = UInt32[1, 1, 1], scan_idx = UInt32[1, 2, 3],
                     weight = Float32[1, 2, 1], irt_obs = Float32[1, 2, 3])
    for r in 1:8; psms[!, Symbol("frag$(r)_int")] = Float32[1, 2, 1] .* Float32(r); end
    Pioneer._add_fragment_chromatogram_features!(psms)
    @test all(psms.n_scans_in_window .== 1)
    @test psms.n_scans == UInt32[3, 3, 3]
end
