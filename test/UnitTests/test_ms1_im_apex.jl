# Ion-mobility apex features of _add_ms1_chromatogram_features! (packet / slice data):
# ms1_apex_offset_im, ms1_apex_offset_im_cycle, ms1_weight_apex_to_m0_apex_im.
using Test
using DataFrames
using Pioneer

@testset "MS1 ion-mobility apex features" begin
    # Precursor 7: six candidate PSMs over two cycles; precursor 9: a singleton.
    # scan s -> im_scans[s], cycle_idxs[s]
    im_scans   = UInt16[100, 108, 116, 100, 108, 300, 50]
    cycle_idxs = Int32[  1,   1,   1,   2,   2,   2,  3]
    psms = DataFrame(
        precursor_idx    = UInt32[7, 7, 7, 7, 7, 7, 9],
        scan_idx         = Int64[1, 2, 3, 4, 5, 6, 7],
        ms1_m0_intensity = Float32[10, 50, 20, 5, 30, 0, 1],
        ms1_m1_intensity = Float32[5, 25, 10, 2, 15, 0, 1],
        weight           = Float32[1, 2, 9, 1, 3, 1, 1],   # weight apex at scan 3 (IM 116)
        irt_obs          = Float32[1, 1, 1, 2, 2, 2, 3],
    )
    Pioneer._add_ms1_chromatogram_features!(psms; im_scans = im_scans, cycle_idxs = cycle_idxs)
    # Global M0 apex is scan 2 (IM 108), whatever the PSM's cycle.
    @test psms.ms1_apex_offset_im[1:6] == Float32[8, 0, 8, 8, 0, 192]
    # Within-cycle apexes: cycle 1 -> scan 2 (IM 108), cycle 2 -> scan 5 (IM 108).
    @test psms.ms1_apex_offset_im_cycle[1:6] == Float32[8, 0, 8, 8, 0, 192]
    # Weight apex (scan 3, IM 116) against M0 apex (IM 108), one value per precursor.
    @test all(psms.ms1_weight_apex_to_m0_apex_im[1:6] .== 8f0)
    # Singleton precursor keeps zeros (fewer than 2 points), like the iRT features.
    @test psms.ms1_apex_offset_im[7] == 0f0
    @test psms.ms1_apex_offset_im_cycle[7] == 0f0
    @test psms.ms1_weight_apex_to_m0_apex_im[7] == 0f0
    # The iRT features are unchanged by the addition.
    @test psms.ms1_apex_offset_irt[1:6] == Float32[0, 0, 0, 1, 1, 1]
end

@testset "MS1 ion-mobility apex features: undefined values and no mobility data" begin
    im_scans   = UInt16[10, 20, 30]
    cycle_idxs = Int32[1, 2, 2]
    # No M0 signal anywhere -> all three undefined (-1).
    psms = DataFrame(precursor_idx = UInt32[1, 1, 1], scan_idx = Int64[1, 2, 3],
                     ms1_m0_intensity = Float32[0, 0, 0], ms1_m1_intensity = Float32[0, 0, 0],
                     weight = Float32[1, 2, 3], irt_obs = Float32[1, 2, 2])
    Pioneer._add_ms1_chromatogram_features!(psms; im_scans = im_scans, cycle_idxs = cycle_idxs)
    @test all(psms.ms1_apex_offset_im .== -1f0)
    @test all(psms.ms1_apex_offset_im_cycle .== -1f0)
    @test all(psms.ms1_weight_apex_to_m0_apex_im .== -1f0)
    # PSM alone in its cycle -> within-cycle offset undefined, global offset defined.
    psms2 = DataFrame(precursor_idx = UInt32[1, 1, 1], scan_idx = Int64[1, 2, 3],
                      ms1_m0_intensity = Float32[4, 9, 3], ms1_m1_intensity = Float32[2, 4, 1],
                      weight = Float32[1, 2, 3], irt_obs = Float32[1, 2, 2])
    Pioneer._add_ms1_chromatogram_features!(psms2; im_scans = im_scans, cycle_idxs = cycle_idxs)
    @test psms2.ms1_apex_offset_im == Float32[10, 0, 10]
    @test psms2.ms1_apex_offset_im_cycle == Float32[-1, 0, 10]
    @test all(psms2.ms1_weight_apex_to_m0_apex_im .== 10f0)     # weight apex scan 3 (IM 30) vs M0 apex scan 2 (IM 20)
    # No mobility data (Thermo / Sciex): the columns exist and are zero.
    psms3 = psms2[:, [:precursor_idx, :scan_idx, :ms1_m0_intensity, :ms1_m1_intensity, :weight, :irt_obs]]
    Pioneer._add_ms1_chromatogram_features!(psms3)
    @test all(psms3.ms1_apex_offset_im .== 0f0)
    @test all(psms3.ms1_apex_offset_im_cycle .== 0f0)
    @test all(psms3.ms1_weight_apex_to_m0_apex_im .== 0f0)
    # Every value finite (feature-finiteness contract).
    for c in (:ms1_apex_offset_im, :ms1_apex_offset_im_cycle, :ms1_weight_apex_to_m0_apex_im)
        @test all(isfinite, psms2[!, c])
    end
end
