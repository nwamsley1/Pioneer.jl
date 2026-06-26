import Pioneer
using Arrow
using DataFrames
using Test

struct MBREmpiricalMockPrecursors <: Pioneer.LibraryPrecursors
    mz::Vector{Float32}
    irt::Vector{Float32}
end

Pioneer.getMz(p::MBREmpiricalMockPrecursors) = p.mz
Pioneer.getIrt(p::MBREmpiricalMockPrecursors) = p.irt

function _mbr_empirical_fragment_lookup()
    frags = Pioneer.CompactFrag{Float32}[
        Pioneer.CompactFrag(
            UInt32(1), 500.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(3), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(1), 600.0f0, Float16(1),
            false, true, false, false,
            UInt8(1), UInt8(4), UInt8(2), UInt8(2), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(2), 600.0f0, Float16(1),
            false, true, false, false,
            UInt8(1), UInt8(4), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(2), 500.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(3), UInt8(2), UInt8(2), UInt8(0),
        ),
    ]
    return Pioneer.StandardFragmentLookup{Float32}(frags, UInt64[1, 3, 5])
end

function _mbr_empirical_main_table(;
    precursor_idx,
    scan_idx,
    ms_file_idx,
    target,
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
        target = Bool.(target),
        cv_fold = zeros(UInt8, n),
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

function _mbr_empirical_pass1_table(; precursor_idx, scan_idx, score)
    return DataFrame(
        precursor_idx = UInt32.(precursor_idx),
        scan_idx = UInt32.(scan_idx),
        trace_prob_prepass = Float32.(score),
        trace_prob_infold = Float32.(score),
    )
end

@testset "MBR empirical smoothed Hellinger aligns spectra by library annotation" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        f2 = joinpath(dir, "run2_fold0.arrow")

        Arrow.write(f1, _mbr_empirical_main_table(
            precursor_idx = [1],
            scan_idx = [1001],
            ms_file_idx = 1,
            target = [true],
            frag1 = [100],
            frag2 = [0],
        ))
        Arrow.write(f1 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [1],
            scan_idx = [1001],
            score = [0.40],
        ))

        Arrow.write(f2, _mbr_empirical_main_table(
            precursor_idx = [2],
            scan_idx = [2001],
            ms_file_idx = 2,
            target = [false],
            frag1 = [0],
            frag2 = [100],
        ))
        Arrow.write(f2 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [2],
            scan_idx = [2001],
            score = [0.91],
        ))

        precursors = MBREmpiricalMockPrecursors(
            Float32[500, 501],
            Float32[10, 10.01],
        )
        fragment_keys = Pioneer.build_mbr_fragment_annotation_keys(_mbr_empirical_fragment_lookup())
        partner_pools = Pioneer.build_counterfactual_partner_pools([f1, f2], precursors)
        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1([f1, f2])

        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            f1,
            donor_dict,
            partner_pools,
            fragment_keys,
        )
        mbr = DataFrame(Arrow.Table(f1 * Pioneer.MBR_SIDECAR_SUFFIX))

        @test isapprox(mbr.MBR_smoothed_frag_hellinger_false[1], 0.0f0; atol = 1.0f-6)
        @test mbr.MBR_is_missing_false[1] == false
        @test mbr.MBR_max_pair_prob_false[1] == 0.91f0
        @test :MBR_smoothed_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_smoothed_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE
    end
end
