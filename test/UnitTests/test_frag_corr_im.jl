# n_scans_in_window and frag_corr_effective_n_im of _add_fragment_chromatogram_features! on ion-mobility
# slice data (timsTOF): PSMs of one precursor spread over mobility slices and cycles.
using Test
using DataFrames
using Pioneer

@testset "n_scans_in_window and frag_corr_effective_n_im on mobility slices" begin
    # Precursor 5: slice 100 in cycles 1-4 (fragment ranks 1-4 track the weight), slice 108 in cycle 2
    # only; precursor 9: a singleton. scan s -> im_scans[s].
    im_scans = UInt16[100, 100, 108, 100, 100, 50]
    cycle    = Int32[1, 2, 2, 3, 4, 7]
    w        = Float32[1, 4, 3, 2, 1, 5]
    psms = DataFrame(precursor_idx = UInt32[5, 5, 5, 5, 5, 9], scan_idx = UInt32[1, 2, 3, 4, 5, 6],
                     cycle_idx = cycle, weight = w, irt_obs = Float32[1, 2, 2, 3, 4, 7])
    for r in 1:8
        psms[!, Symbol("frag$(r)_int")] = r <= 4 ? w .* Float32(r) : zeros(Float32, 6)
    end
    Pioneer._add_fragment_chromatogram_features!(psms; im_scans = im_scans)
    @test psms.n_scans_in_window == UInt32[1, 2, 2, 1, 1, 1]
    @test psms.n_scans[1] == 5
    # Slice 100 spans four cycles with ranks 1-4 perfectly correlated to the weight.
    rw = Pioneer._fragment_rank_weights(8)
    expected = Pioneer._positive_corr_summary(Float32[1, 1, 1, 1, 0, 0, 0, 0], rw)[2]
    @test all(psms.frag_corr_effective_n_im[[1, 2, 4, 5]] .≈ expected)
    @test psms.frag_corr_effective_n_im[3] == 0f0      # slice 108 holds a single PSM
    @test psms.frag_corr_effective_n_im[6] == 0f0      # singleton precursor
    # The global feature is unchanged: all five PSMs of precursor 5, still perfectly correlated.
    @test psms.frag_corr_effective_n[1] ≈ expected
    @test psms.frag_corr_effective_n[3] ≈ expected
end

@testset "frag_corr_effective_n_im: uncorrelated slice gets a lower value than the correlated one" begin
    im_scans = UInt16[100, 100, 100, 116, 116, 116]
    w        = Float32[1, 5, 1, 1, 5, 1]
    psms = DataFrame(precursor_idx = UInt32[2, 2, 2, 2, 2, 2], scan_idx = UInt32[1, 2, 3, 4, 5, 6],
                     cycle_idx = Int32[1, 2, 3, 1, 2, 3], weight = w, irt_obs = Float32[1, 2, 3, 1, 2, 3])
    for r in 1:8
        # slice 100: fragments follow the weight; slice 116: fragments anti-follow it
        psms[!, Symbol("frag$(r)_int")] = Float32[1, 5, 1, 5, 1, 5] .* Float32(r)
    end
    Pioneer._add_fragment_chromatogram_features!(psms; im_scans = im_scans)
    @test psms.frag_corr_effective_n_im[1] > 0f0
    @test psms.frag_corr_effective_n_im[4] == 0f0
    @test all(psms.n_scans_in_window .== 2)
end

@testset "no mobility data: n_scans_in_window is 1, effective_n_im is 0" begin
    psms = DataFrame(precursor_idx = UInt32[1, 1, 1], scan_idx = UInt32[1, 2, 3], cycle_idx = Int32[1, 2, 3],
                     weight = Float32[1, 2, 1], irt_obs = Float32[1, 2, 3])
    for r in 1:8; psms[!, Symbol("frag$(r)_int")] = Float32[1, 2, 1] .* Float32(r); end
    Pioneer._add_fragment_chromatogram_features!(psms)
    @test all(psms.n_scans_in_window .== 1)
    @test all(psms.frag_corr_effective_n_im .== 0f0)
    @test all(isfinite, psms.frag_corr_effective_n_im)
end
