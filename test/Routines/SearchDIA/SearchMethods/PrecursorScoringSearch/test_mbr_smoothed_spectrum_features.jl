import Pioneer
using Arrow
using DataFrames
using Test

function _mbr_smoothed_main_table(;
    precursor_idx,
    scan_idx,
    ms_file_idx,
    frag1,
    frag2,
)
    n = length(precursor_idx)
    return DataFrame(
        precursor_idx = UInt32.(precursor_idx),
        scan_idx = UInt32.(scan_idx),
        weight = ones(Float32, n),
        log2_intensity_explained = fill(Float16(1), n),
        irt_pred = zeros(Float32, n),
        irt_obs = zeros(Float32, n),
        log_by_ratio_m0 = zeros(Float16, n),
        rt = zeros(Float32, n),
        ms_file_idx = fill(UInt32(ms_file_idx), n),
        frag1_smoothed_intensity = Float32.(frag1),
        frag2_smoothed_intensity = Float32.(frag2),
        frag3_smoothed_intensity = zeros(Float32, n),
        frag4_smoothed_intensity = zeros(Float32, n),
        frag5_smoothed_intensity = zeros(Float32, n),
        frag6_smoothed_intensity = zeros(Float32, n),
        frag7_smoothed_intensity = zeros(Float32, n),
        frag8_smoothed_intensity = zeros(Float32, n),
    )
end

function _mbr_smoothed_pass1_table(; precursor_idx, scan_idx, score)
    n = length(precursor_idx)
    return DataFrame(
        precursor_idx = UInt32.(precursor_idx),
        scan_idx = UInt32.(scan_idx),
        trace_prob_prepass = Float32.(score),
        trace_prob_infold = Float32.(score),
    )
end

@testset "MBR smoothed spectrum Hellinger features are paired true/false evidence" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        f2 = joinpath(dir, "run2_fold0.arrow")

        Arrow.write(f1, _mbr_smoothed_main_table(
            precursor_idx = [20, 22],
            scan_idx = [1001, 1002],
            ms_file_idx = 1,
            frag1 = [100, 50],
            frag2 = [0, 50],
        ))
        Arrow.write(f1 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_smoothed_pass1_table(
            precursor_idx = [20, 22],
            scan_idx = [1001, 1002],
            score = [0.95, 0.40],
        ))

        Arrow.write(f2, _mbr_smoothed_main_table(
            precursor_idx = [20, 21],
            scan_idx = [2001, 2002],
            ms_file_idx = 2,
            frag1 = [100, 0],
            frag2 = [0, 100],
        ))
        Arrow.write(f2 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_smoothed_pass1_table(
            precursor_idx = [20, 21],
            scan_idx = [2001, 2002],
            score = [0.90, 0.80],
        ))

        partner_col = zeros(UInt32, 22)
        partner_col[20] = UInt32(21)

        donors = Pioneer.build_mbr_donor_dict_streaming_with_pass1([f1, f2])
        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(f1, donors, partner_col)
        mbr = DataFrame(Arrow.Table(f1 * Pioneer.MBR_SIDECAR_SUFFIX))

        @test isapprox(mbr.MBR_smoothed_frag_hellinger_true[1], 0.0f0; atol = 1.0f-6)
        @test isapprox(mbr.MBR_smoothed_frag_hellinger_false[1], 1.0f0; atol = 1.0f-6)
        @test mbr.MBR_smoothed_frag_hellinger_true[2] == -1.0f0
        @test mbr.MBR_smoothed_frag_hellinger_false[2] == -1.0f0
    end
end

@testset "MBR smoothed spectrum FTR features map true to false positionally" begin
    true_idx = findfirst(==(:MBR_smoothed_frag_hellinger_true), Pioneer.FTR_FEATURES_F_TRUE)

    @test true_idx !== nothing
    @test Pioneer.FTR_FEATURES_F_FALSE[true_idx] == :MBR_smoothed_frag_hellinger_false
end
