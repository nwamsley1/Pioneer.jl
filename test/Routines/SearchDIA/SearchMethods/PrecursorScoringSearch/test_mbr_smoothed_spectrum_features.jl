import Pioneer
using Arrow
using DataFrames
using Test

struct MBRMockPrecursors <: Pioneer.LibraryPrecursors
    mz::Vector{Float32}
    irt::Vector{Float32}
end

Pioneer.getMz(p::MBRMockPrecursors) = p.mz
Pioneer.getIrt(p::MBRMockPrecursors) = p.irt

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

function _mbr_test_fragment_lookup(rank_positions_by_pid::Dict{UInt32, Vector{UInt8}})
    max_pid = maximum(keys(rank_positions_by_pid))
    frags = Pioneer.CompactFrag{Float32}[]
    ranges = Vector{UInt64}(undef, Int(max_pid) + 1)
    next_idx = UInt64(1)
    @inbounds for pid in UInt32(1):max_pid
        ranges[Int(pid)] = next_idx
        positions = get(rank_positions_by_pid, pid, UInt8[])
        for (rank, ion_position) in enumerate(positions)
            push!(frags, Pioneer.CompactFrag(
                pid,
                Float32(100 + Int(pid) + rank),
                Float16(1),
                true,
                false,
                false,
                false,
                UInt8(1),
                ion_position,
                UInt8(2),
                UInt8(rank),
                UInt8(0),
            ))
            next_idx += UInt64(1)
        end
    end
    ranges[Int(max_pid) + 1] = next_idx
    return Pioneer.StandardFragmentLookup(frags, ranges)
end

function _mbr_test_annotation_keys(rank_positions_by_pid::Dict{UInt32, Vector{UInt8}})
    lookup = _mbr_test_fragment_lookup(rank_positions_by_pid)
    return Pioneer._mbr_fragment_annotation_key_table(lookup, maximum(keys(rank_positions_by_pid)))
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
        annotation_keys = _mbr_test_annotation_keys(Dict(
            UInt32(20) => UInt8[1, 2],
            UInt32(21) => UInt8[1, 2],
            UInt32(22) => UInt8[1, 2],
        ))

        donors = Pioneer.build_mbr_donor_dict_streaming_with_pass1([f1, f2])
        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            f1,
            donors,
            partner_col,
            annotation_keys,
        )
        mbr = DataFrame(Arrow.Table(f1 * Pioneer.MBR_SIDECAR_SUFFIX))

        @test isapprox(mbr.MBR_smoothed_frag_hellinger_true[1], 0.0f0; atol = 1.0f-6)
        @test isapprox(mbr.MBR_smoothed_frag_hellinger_false[1], 1.0f0; atol = 1.0f-6)
        @test mbr.MBR_smoothed_frag_hellinger_true[2] == -1.0f0
        @test mbr.MBR_smoothed_frag_hellinger_false[2] == -1.0f0
    end
end

@testset "MBR smoothed spectrum Hellinger aligns counterfactual fragments by library annotation key" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        f2 = joinpath(dir, "run2_fold0.arrow")

        Arrow.write(f1, _mbr_smoothed_main_table(
            precursor_idx = [30],
            scan_idx = [3001],
            ms_file_idx = 1,
            frag1 = [100],
            frag2 = [0],
        ))
        Arrow.write(f1 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_smoothed_pass1_table(
            precursor_idx = [30],
            scan_idx = [3001],
            score = [0.95],
        ))

        Arrow.write(f2, _mbr_smoothed_main_table(
            precursor_idx = [31],
            scan_idx = [3002],
            ms_file_idx = 2,
            frag1 = [0],
            frag2 = [100],
        ))
        Arrow.write(f2 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_smoothed_pass1_table(
            precursor_idx = [31],
            scan_idx = [3002],
            score = [0.90],
        ))

        partner_col = zeros(UInt32, 31)
        partner_col[30] = UInt32(31)
        annotation_keys = _mbr_test_annotation_keys(Dict(
            UInt32(30) => UInt8[1, 2],
            UInt32(31) => UInt8[2, 1],
        ))

        donors = Pioneer.build_mbr_donor_dict_streaming_with_pass1([f1, f2])
        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            f1,
            donors,
            partner_col,
            annotation_keys,
        )
        mbr = DataFrame(Arrow.Table(f1 * Pioneer.MBR_SIDECAR_SUFFIX))

        @test isapprox(mbr.MBR_smoothed_frag_hellinger_false[1], 0.0f0; atol = 1.0f-6)
    end
end

@testset "MBR false donor selection uses alternate partner with cross-file donor" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        f2 = joinpath(dir, "run2_fold0.arrow")

        Arrow.write(f1, _mbr_smoothed_main_table(
            precursor_idx = [40, 41],
            scan_idx = [4001, 4002],
            ms_file_idx = 1,
            frag1 = [100, 0],
            frag2 = [0, 100],
        ))
        Arrow.write(f1 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_smoothed_pass1_table(
            precursor_idx = [40, 41],
            scan_idx = [4001, 4002],
            score = [0.40, 0.95],
        ))

        Arrow.write(f2, _mbr_smoothed_main_table(
            precursor_idx = [40, 42],
            scan_idx = [5001, 5002],
            ms_file_idx = 2,
            frag1 = [100, 100],
            frag2 = [0, 0],
        ))
        Arrow.write(f2 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_smoothed_pass1_table(
            precursor_idx = [40, 42],
            scan_idx = [5001, 5002],
            score = [0.90, 0.88],
        ))

        partner_candidates = zeros(UInt32, 2, 42)
        partner_candidates[1, 40] = UInt32(41)
        partner_candidates[2, 40] = UInt32(42)
        annotation_keys = _mbr_test_annotation_keys(Dict(
            UInt32(40) => UInt8[1, 2],
            UInt32(41) => UInt8[1, 2],
            UInt32(42) => UInt8[1, 2],
        ))

        donors = Pioneer.build_mbr_donor_dict_streaming_with_pass1([f1, f2])
        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            f1,
            donors,
            partner_candidates,
            annotation_keys,
        )
        mbr = DataFrame(Arrow.Table(f1 * Pioneer.MBR_SIDECAR_SUFFIX))

        @test mbr.MBR_is_missing_false[1] == false
        @test mbr.MBR_max_pair_prob_false[1] == Float32(0.88)
    end
end

@testset "MBR counterfactual partner candidate matrix ranks nearest iRT partners" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        Arrow.write(f1, DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            target = Bool[true, false, false],
            cv_fold = UInt8[0, 0, 0],
        ))
        precursors = MBRMockPrecursors(
            Float32[500, 500, 500],
            Float32[10, 11, 30],
        )

        candidates = Pioneer.build_counterfactual_partner_candidate_matrix([f1], precursors, 2)

        @test candidates[1, 1] == UInt32(2)
        @test candidates[2, 1] == UInt32(3)
    end
end
