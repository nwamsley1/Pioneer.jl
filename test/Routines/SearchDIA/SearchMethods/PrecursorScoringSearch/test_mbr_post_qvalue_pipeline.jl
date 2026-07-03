import Pioneer
using Arrow
using DataFrames
using Test

struct MBRPostQMockPrecursors <: Pioneer.LibraryPrecursors
    mz::Vector{Float32}
    irt::Vector{Float32}
end

Pioneer.getMz(p::MBRPostQMockPrecursors) = p.mz
Pioneer.getIrt(p::MBRPostQMockPrecursors) = p.irt
Pioneer.getCharge(p::MBRPostQMockPrecursors) = fill(UInt8(2), length(p.mz))
Pioneer.getLength(p::MBRPostQMockPrecursors) = fill(UInt8(9), length(p.mz))

struct MBRPostQMockPrecursorsWithLength <: Pioneer.LibraryPrecursors
    mz::Vector{Float32}
    irt::Vector{Float32}
    length::Vector{UInt8}
end

Pioneer.getMz(p::MBRPostQMockPrecursorsWithLength) = p.mz
Pioneer.getIrt(p::MBRPostQMockPrecursorsWithLength) = p.irt
Pioneer.getCharge(p::MBRPostQMockPrecursorsWithLength) = fill(UInt8(2), length(p.mz))
Pioneer.getLength(p::MBRPostQMockPrecursorsWithLength) = p.length

function _mbr_post_q_fragment_lookup()
    frags = Pioneer.CompactFrag{Float32}[
        Pioneer.CompactFrag(
            UInt32(1), 500.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(3), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(2), 600.0f0, Float16(1),
            false, true, false, false,
            UInt8(1), UInt8(4), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(3), 700.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(5), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(4), 800.0f0, Float16(1),
            false, true, false, false,
            UInt8(1), UInt8(6), UInt8(2), UInt8(1), UInt8(0),
        ),
    ]
    return Pioneer.StandardFragmentLookup{Float32}(frags, UInt64[1, 2, 3, 4, 5])
end

function _mbr_post_q_main_table(; ms_file_idx, scan_offset = 0)
    precursor_idx = UInt32[1, 2]
    scan_idx = UInt32[scan_offset + 101, scan_offset + 102]
    score = Float32[0.95, 0.05]
    n = length(precursor_idx)
    return DataFrame(
        precursor_idx = precursor_idx,
        scan_idx = scan_idx,
        weight = Float32[10, 1],
        log2_intensity_explained = Float16[1, 0.25],
        irt_pred = Float32[10, 20],
        irt_obs = Float32[10.1, 20.1],
        log_by_ratio_m0 = zeros(Float16, n),
        rt = Float32[5, 6],
        n_scans = Float32[8, 4],
        main_pep = Float32[0.01, 0.99],
        smoothed_2d_shadow_hellinger = Float32[0.1, 0.2],
        ms_file_idx = fill(UInt32(ms_file_idx), n),
        target = Bool[true, false],
        cv_fold = UInt8[0, 1],
        trace_prob = score,
        trace_prob_prepass = score,
        trace_prob_infold = score,
        prec_prob = score,
        qval = Float32[0, 1],
        pep = Float32[0, 1],
        global_qval = Float32[0, 1],
        global_pep = Float32[0, 1],
        mbr_recovered = falses(n),
        frag1_smoothed_intensity = Float32[100, 0],
        frag2_smoothed_intensity = Float32[0, 100],
        frag3_smoothed_intensity = zeros(Float32, n),
        frag4_smoothed_intensity = zeros(Float32, n),
        frag5_smoothed_intensity = zeros(Float32, n),
        frag6_smoothed_intensity = zeros(Float32, n),
        frag7_smoothed_intensity = zeros(Float32, n),
        frag8_smoothed_intensity = zeros(Float32, n),
    )
end

function _mbr_post_q_candidate_table(; ms_file_idx, scan_offset = 0, passing = false)
    precursor_idx = UInt32[1, 2, 3, 4]
    scan_idx = UInt32[scan_offset + 101, scan_offset + 102, scan_offset + 103, scan_offset + 104]
    score = passing ? Float32[0.95, 0.90, 0.96, 0.91] :
                      Float32[0.50, 0.40, 0.51, 0.41]
    qval = passing ? fill(0.001f0, 4) : fill(0.05f0, 4)
    n = length(precursor_idx)
    return DataFrame(
        precursor_idx = precursor_idx,
        scan_idx = scan_idx,
        weight = Float32[10, 2, 11, 3],
        log2_intensity_explained = Float16[1, 0.25, 1, 0.25],
        irt_pred = Float32[10, 10.05, 20, 20.05],
        irt_obs = Float32[10.1, 10.1, 20.1, 20.1],
        log_by_ratio_m0 = zeros(Float16, n),
        rt = Float32[5, 5.1, 6, 6.1],
        n_scans = Float32[8, 4, 9, 5],
        main_pep = Float32[0.01, 0.99, 0.01, 0.99],
        smoothed_2d_shadow_hellinger = Float32[0.1, 0.2, 0.1, 0.2],
        ms_file_idx = fill(UInt32(ms_file_idx), n),
        target = Bool[true, false, true, false],
        cv_fold = UInt8[0, 0, 1, 1],
        trace_prob = score,
        trace_prob_prepass = score,
        trace_prob_infold = score,
        prec_prob = score,
        qval = qval,
        pep = qval,
        global_qval = fill(0.001f0, n),
        global_pep = fill(0.001f0, n),
        mbr_recovered = falses(n),
        frag1_smoothed_intensity = Float32[100, 0, 100, 0],
        frag2_smoothed_intensity = Float32[0, 100, 0, 100],
        frag3_smoothed_intensity = zeros(Float32, n),
        frag4_smoothed_intensity = zeros(Float32, n),
        frag5_smoothed_intensity = zeros(Float32, n),
        frag6_smoothed_intensity = zeros(Float32, n),
        frag7_smoothed_intensity = zeros(Float32, n),
        frag8_smoothed_intensity = zeros(Float32, n),
    )
end

@testset "post-qvalue MBR runs on filtered files with pass1 scores in main" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1.arrow")
        f2 = joinpath(dir, "run2.arrow")
        excluded = joinpath(dir, "excluded.arrow")

        Arrow.write(f1, _mbr_post_q_main_table(ms_file_idx = 1, scan_offset = 1000))
        Arrow.write(f2, _mbr_post_q_main_table(ms_file_idx = 2, scan_offset = 2000))
        Arrow.write(excluded, _mbr_post_q_main_table(ms_file_idx = 3, scan_offset = 3000))

        refs = Pioneer.PSMFileReference[
            Pioneer.PSMFileReference(f1),
            Pioneer.PSMFileReference(f2),
        ]
        precursors = MBRPostQMockPrecursors(Float32[500, 600], Float32[10, 20])

        @test !isfile(f1 * Pioneer.PASS1_SIDECAR_SUFFIX)
        summary = Pioneer.run_mbr_after_qvalue_filter!(
            refs,
            precursors,
            _mbr_post_q_fragment_lookup(),
        )

        @test summary.n_files == 2
        @test summary.n_recovery_sidecars == 2
        @test summary.n_merged == 2

        for path in (f1, f2)
            df = DataFrame(Arrow.Table(path))
            @test :mbr_recovered in propertynames(df)
            @test :MBR_transfer_candidate in propertynames(df)
            @test :mbr_target_decoy_prob in propertynames(df)
            @test :ftr_qval_true in propertynames(df)
            @test :ftr_pep_true in propertynames(df)
            @test !isfile(path * Pioneer.PASS1_SIDECAR_SUFFIX)
            @test !isfile(path * Pioneer.MBR_SIDECAR_SUFFIX)
            @test !isfile(path * Pioneer.RECOVERY_SIDECAR_SUFFIX)
        end

        excluded_df = DataFrame(Arrow.Table(excluded))
        @test !(:MBR_transfer_candidate in propertynames(excluded_df))
        @test !isfile(excluded * Pioneer.PASS1_SIDECAR_SUFFIX)
    end
end

@testset "post-qvalue MBR evaluates failed receiver rows with passing cross-run donors" begin
    mktempdir() do dir
        candidate_path = joinpath(dir, "candidate_unfiltered.arrow")
        donor_path = joinpath(dir, "donor_filtered.arrow")

        Arrow.write(candidate_path, _mbr_post_q_candidate_table(ms_file_idx = 1, scan_offset = 1000))
        Arrow.write(donor_path, _mbr_post_q_candidate_table(ms_file_idx = 2, scan_offset = 2000, passing = true))

        candidate_refs = Pioneer.PSMFileReference[Pioneer.PSMFileReference(candidate_path)]
        donor_refs = Pioneer.PSMFileReference[Pioneer.PSMFileReference(donor_path)]
        precursors = MBRPostQMockPrecursors(
            Float32[500, 600, 700, 800],
            Float32[10, 10.05, 20, 20.05],
        )

        summary = Pioneer.run_mbr_after_qvalue_filter!(
            candidate_refs,
            donor_refs,
            precursors,
            _mbr_post_q_fragment_lookup(),
        )

        @test summary.n_candidates > 0

        candidate_df = DataFrame(Arrow.Table(candidate_path))
        @test any(candidate_df.MBR_transfer_candidate)
        @test !any(candidate_df.qval .<= 0.01f0)
        @test !isfile(candidate_path * Pioneer.PASS1_SIDECAR_SUFFIX)
        @test !isfile(donor_path * Pioneer.PASS1_SIDECAR_SUFFIX)
    end
end

@testset "MBR counterfactuals only pair to receiver-run eligible precursors" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        donor_path = joinpath(dir, "donor.arrow")

        receiver_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            scan_idx = UInt32[101, 102, 103],
            weight = Float32[10, 2, 3],
            log2_intensity_explained = Float16[1, 0.5, 0.5],
            irt_pred = Float32[10.0, 10.02, 10.20],
            irt_obs = Float32[10.1, 10.1, 10.1],
            log_by_ratio_m0 = zeros(Float16, 3),
            rt = Float32[5.0, 5.1, 5.2],
            n_scans = Float32[8, 4, 5],
            main_pep = Float32[0.99, 0.001, 0.99],
            smoothed_2d_shadow_hellinger = Float32[0.1, 0.2, 0.3],
            ms_file_idx = fill(UInt32(1), 3),
            target = trues(3),
            cv_fold = zeros(UInt8, 3),
            trace_prob = Float32[0.40, 0.99, 0.50],
            trace_prob_prepass = Float32[0.40, 0.99, 0.50],
            trace_prob_infold = Float32[0.40, 0.99, 0.50],
            prec_prob = Float32[0.40, 0.99, 0.50],
            qval = Float32[0.05, 0.001, 0.05],
            pep = Float32[0.05, 0.001, 0.05],
            global_qval = Float32[0.001, 0.001, 0.001],
            global_pep = Float32[0.001, 0.001, 0.001],
            mbr_recovered = falses(3),
            frag1_smoothed_intensity = Float32[100, 80, 60],
            frag2_smoothed_intensity = Float32[0, 20, 40],
            frag3_smoothed_intensity = zeros(Float32, 3),
            frag4_smoothed_intensity = zeros(Float32, 3),
            frag5_smoothed_intensity = zeros(Float32, 3),
            frag6_smoothed_intensity = zeros(Float32, 3),
            frag7_smoothed_intensity = zeros(Float32, 3),
            frag8_smoothed_intensity = zeros(Float32, 3),
        )
        donor_df = receiver_df[2:3, :]
        donor_df[!, :ms_file_idx] .= UInt32(2)
        donor_df[!, :scan_idx] .= UInt32[202, 203]
        donor_df[!, :trace_prob] .= Float32[0.99, 0.70]
        donor_df[!, :trace_prob_prepass] .= Float32[0.99, 0.70]
        donor_df[!, :trace_prob_infold] .= Float32[0.99, 0.70]
        donor_df[!, :prec_prob] .= Float32[0.99, 0.70]

        Arrow.write(receiver_path, receiver_df)
        Arrow.write(donor_path, donor_df)
        for (path, df) in ((receiver_path, receiver_df), (donor_path, donor_df))
            Arrow.write(path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = df.precursor_idx,
                scan_idx = df.scan_idx,
                trace_prob_prepass = df.trace_prob_prepass,
                trace_prob_infold = df.trace_prob_infold,
            ))
        end

        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1([donor_path])
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [receiver_path],
            MBRPostQMockPrecursors(Float32[500, 501, 502], Float32[10.0, 10.02, 10.20]),
        )
        eligible_by_file = Pioneer.build_counterfactual_receiver_eligibility(
            [receiver_path];
            q_value_threshold = 0.01f0,
        )

        side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            Pioneer.build_mbr_fragment_annotation_keys(_mbr_post_q_fragment_lookup()),
            counterfactual_eligibility_by_file = eligible_by_file,
            passing_score_floor = 0.0f0,
        )
        side = DataFrame(Arrow.Table(side_path))

        @test side.MBR_best_pair_prob_false[1] == 0.70f0
        @test !side.MBR_best_is_missing_false[1]
    end
end

@testset "MBR counterfactuals can pair to dataset-passing precursors absent from receiver run" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        other_path = joinpath(dir, "other_run.arrow")

        receiver_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            scan_idx = UInt32[101, 102],
            weight = Float32[10, 2],
            log2_intensity_explained = Float16[1, 0.5],
            irt_pred = Float32[10.0, 10.02],
            irt_obs = Float32[10.1, 10.1],
            log_by_ratio_m0 = zeros(Float16, 2),
            rt = Float32[5.0, 5.1],
            n_scans = Float32[8, 4],
            main_pep = Float32[0.99, 0.001],
            smoothed_2d_shadow_hellinger = Float32[0.1, 0.2],
            ms_file_idx = fill(UInt32(1), 2),
            target = trues(2),
            cv_fold = zeros(UInt8, 2),
            trace_prob = Float32[0.40, 0.99],
            trace_prob_prepass = Float32[0.40, 0.99],
            trace_prob_infold = Float32[0.40, 0.99],
            prec_prob = Float32[0.40, 0.99],
            qval = Float32[0.05, 0.001],
            pep = Float32[0.05, 0.001],
            global_qval = Float32[0.001, 0.001],
            global_pep = Float32[0.001, 0.001],
            mbr_recovered = falses(2),
            frag1_smoothed_intensity = Float32[100, 80],
            frag2_smoothed_intensity = Float32[0, 20],
            frag3_smoothed_intensity = zeros(Float32, 2),
            frag4_smoothed_intensity = zeros(Float32, 2),
            frag5_smoothed_intensity = zeros(Float32, 2),
            frag6_smoothed_intensity = zeros(Float32, 2),
            frag7_smoothed_intensity = zeros(Float32, 2),
            frag8_smoothed_intensity = zeros(Float32, 2),
        )
        other_df = DataFrame(
            precursor_idx = UInt32[2, 3],
            scan_idx = UInt32[202, 203],
            weight = Float32[2, 3],
            log2_intensity_explained = Float16[0.5, 0.5],
            irt_pred = Float32[10.02, 10.20],
            irt_obs = Float32[10.1, 10.1],
            log_by_ratio_m0 = zeros(Float16, 2),
            rt = Float32[5.1, 5.2],
            n_scans = Float32[4, 5],
            main_pep = Float32[0.001, 0.001],
            smoothed_2d_shadow_hellinger = Float32[0.2, 0.3],
            ms_file_idx = fill(UInt32(2), 2),
            target = trues(2),
            cv_fold = zeros(UInt8, 2),
            trace_prob = Float32[0.99, 0.70],
            trace_prob_prepass = Float32[0.99, 0.70],
            trace_prob_infold = Float32[0.99, 0.70],
            prec_prob = Float32[0.99, 0.70],
            qval = Float32[0.001, 0.001],
            pep = Float32[0.001, 0.001],
            global_qval = Float32[0.001, 0.001],
            global_pep = Float32[0.001, 0.001],
            mbr_recovered = falses(2),
            frag1_smoothed_intensity = Float32[80, 60],
            frag2_smoothed_intensity = Float32[20, 40],
            frag3_smoothed_intensity = zeros(Float32, 2),
            frag4_smoothed_intensity = zeros(Float32, 2),
            frag5_smoothed_intensity = zeros(Float32, 2),
            frag6_smoothed_intensity = zeros(Float32, 2),
            frag7_smoothed_intensity = zeros(Float32, 2),
            frag8_smoothed_intensity = zeros(Float32, 2),
        )

        Arrow.write(receiver_path, receiver_df)
        Arrow.write(other_path, other_df)
        for (path, df) in ((receiver_path, receiver_df), (other_path, other_df))
            Arrow.write(path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = df.precursor_idx,
                scan_idx = df.scan_idx,
                trace_prob_prepass = df.trace_prob_prepass,
                trace_prob_infold = df.trace_prob_infold,
            ))
        end

        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1([other_path])
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [receiver_path, other_path],
            MBRPostQMockPrecursors(Float32[500, 501, 502], Float32[10.0, 10.02, 10.20]),
        )
        eligible_by_file = Pioneer.build_counterfactual_receiver_eligibility(
            [receiver_path, other_path];
            q_value_threshold = 0.01f0,
        )

        side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            Pioneer.build_mbr_fragment_annotation_keys(_mbr_post_q_fragment_lookup()),
            counterfactual_eligibility_by_file = eligible_by_file,
            passing_score_floor = 0.0f0,
        )
        side = DataFrame(Arrow.Table(side_path))

        @test side.MBR_best_pair_prob_false[1] == 0.70f0
        @test !side.MBR_best_is_missing_false[1]
    end
end

@testset "MBR counterfactuals can pair to same-charge different-length precursors" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        other_path = joinpath(dir, "other_run.arrow")

        receiver_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            scan_idx = UInt32[101, 102],
            weight = Float32[10, 2],
            log2_intensity_explained = Float16[1, 0.5],
            irt_pred = Float32[10.0, 10.02],
            irt_obs = Float32[10.1, 10.1],
            log_by_ratio_m0 = zeros(Float16, 2),
            rt = Float32[5.0, 5.1],
            n_scans = Float32[8, 4],
            main_pep = Float32[0.99, 0.001],
            smoothed_2d_shadow_hellinger = Float32[0.1, 0.2],
            ms_file_idx = fill(UInt32(1), 2),
            target = trues(2),
            cv_fold = zeros(UInt8, 2),
            trace_prob = Float32[0.40, 0.99],
            trace_prob_prepass = Float32[0.40, 0.99],
            trace_prob_infold = Float32[0.40, 0.99],
            prec_prob = Float32[0.40, 0.99],
            qval = Float32[0.05, 0.001],
            pep = Float32[0.05, 0.001],
            global_qval = Float32[0.001, 0.001],
            global_pep = Float32[0.001, 0.001],
            mbr_recovered = falses(2),
            frag1_smoothed_intensity = Float32[100, 80],
            frag2_smoothed_intensity = Float32[0, 20],
            frag3_smoothed_intensity = zeros(Float32, 2),
            frag4_smoothed_intensity = zeros(Float32, 2),
            frag5_smoothed_intensity = zeros(Float32, 2),
            frag6_smoothed_intensity = zeros(Float32, 2),
            frag7_smoothed_intensity = zeros(Float32, 2),
            frag8_smoothed_intensity = zeros(Float32, 2),
        )
        other_df = DataFrame(
            precursor_idx = UInt32[3],
            scan_idx = UInt32[203],
            weight = Float32[3],
            log2_intensity_explained = Float16[0.5],
            irt_pred = Float32[10.20],
            irt_obs = Float32[10.1],
            log_by_ratio_m0 = zeros(Float16, 1),
            rt = Float32[5.2],
            n_scans = Float32[5],
            main_pep = Float32[0.001],
            smoothed_2d_shadow_hellinger = Float32[0.3],
            ms_file_idx = fill(UInt32(2), 1),
            target = trues(1),
            cv_fold = zeros(UInt8, 1),
            trace_prob = Float32[0.70],
            trace_prob_prepass = Float32[0.70],
            trace_prob_infold = Float32[0.70],
            prec_prob = Float32[0.70],
            qval = Float32[0.001],
            pep = Float32[0.001],
            global_qval = Float32[0.001],
            global_pep = Float32[0.001],
            mbr_recovered = falses(1),
            frag1_smoothed_intensity = Float32[60],
            frag2_smoothed_intensity = Float32[40],
            frag3_smoothed_intensity = zeros(Float32, 1),
            frag4_smoothed_intensity = zeros(Float32, 1),
            frag5_smoothed_intensity = zeros(Float32, 1),
            frag6_smoothed_intensity = zeros(Float32, 1),
            frag7_smoothed_intensity = zeros(Float32, 1),
            frag8_smoothed_intensity = zeros(Float32, 1),
        )

        Arrow.write(receiver_path, receiver_df)
        Arrow.write(other_path, other_df)
        for (path, df) in ((receiver_path, receiver_df), (other_path, other_df))
            Arrow.write(path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = df.precursor_idx,
                scan_idx = df.scan_idx,
                trace_prob_prepass = df.trace_prob_prepass,
                trace_prob_infold = df.trace_prob_infold,
            ))
        end

        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1([other_path])
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [receiver_path, other_path],
            MBRPostQMockPrecursorsWithLength(
                Float32[500, 501, 502],
                Float32[10.0, 10.02, 10.20],
                UInt8[9, 9, 12],
            ),
        )
        eligible_by_file = Pioneer.build_counterfactual_receiver_eligibility(
            [receiver_path, other_path];
            q_value_threshold = 0.01f0,
        )

        side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            Pioneer.build_mbr_fragment_annotation_keys(_mbr_post_q_fragment_lookup()),
            counterfactual_eligibility_by_file = eligible_by_file,
            passing_score_floor = 0.0f0,
        )
        side = DataFrame(Arrow.Table(side_path))

        @test side.MBR_best_pair_prob_false[1] == 0.70f0
        @test !side.MBR_best_is_missing_false[1]
    end
end

@testset "MBR counterfactuals choose closest m/z within receiver iRT tolerance" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        other_path = joinpath(dir, "other_run.arrow")

        receiver_df = DataFrame(
            precursor_idx = UInt32[1],
            scan_idx = UInt32[101],
            weight = Float32[10],
            log2_intensity_explained = Float16[1],
            irt_pred = Float32[10.0],
            irt_obs = Float32[10.1],
            log_by_ratio_m0 = zeros(Float16, 1),
            rt = Float32[5.0],
            n_scans = Float32[8],
            main_pep = Float32[0.99],
            smoothed_2d_shadow_hellinger = Float32[0.1],
            ms_file_idx = fill(UInt32(1), 1),
            target = trues(1),
            cv_fold = zeros(UInt8, 1),
            trace_prob = Float32[0.40],
            trace_prob_prepass = Float32[0.40],
            trace_prob_infold = Float32[0.40],
            prec_prob = Float32[0.40],
            qval = Float32[0.05],
            pep = Float32[0.05],
            global_qval = Float32[0.001],
            global_pep = Float32[0.001],
            mbr_recovered = falses(1),
            frag1_smoothed_intensity = Float32[100],
            frag2_smoothed_intensity = Float32[0],
            frag3_smoothed_intensity = zeros(Float32, 1),
            frag4_smoothed_intensity = zeros(Float32, 1),
            frag5_smoothed_intensity = zeros(Float32, 1),
            frag6_smoothed_intensity = zeros(Float32, 1),
            frag7_smoothed_intensity = zeros(Float32, 1),
            frag8_smoothed_intensity = zeros(Float32, 1),
        )
        other_df = DataFrame(
            precursor_idx = UInt32[2, 3, 4],
            scan_idx = UInt32[202, 203, 204],
            weight = Float32[2, 3, 4],
            log2_intensity_explained = Float16[0.5, 0.5, 0.5],
            irt_pred = Float32[10.01, 10.50, 12.00],
            irt_obs = Float32[10.1, 10.1, 10.1],
            log_by_ratio_m0 = zeros(Float16, 3),
            rt = Float32[5.1, 5.2, 5.3],
            n_scans = Float32[4, 5, 6],
            main_pep = Float32[0.001, 0.001, 0.001],
            smoothed_2d_shadow_hellinger = Float32[0.2, 0.3, 0.4],
            ms_file_idx = fill(UInt32(2), 3),
            target = trues(3),
            cv_fold = zeros(UInt8, 3),
            trace_prob = Float32[0.60, 0.80, 0.95],
            trace_prob_prepass = Float32[0.60, 0.80, 0.95],
            trace_prob_infold = Float32[0.60, 0.80, 0.95],
            prec_prob = Float32[0.60, 0.80, 0.95],
            qval = Float32[0.001, 0.001, 0.001],
            pep = Float32[0.001, 0.001, 0.001],
            global_qval = Float32[0.001, 0.001, 0.001],
            global_pep = Float32[0.001, 0.001, 0.001],
            mbr_recovered = falses(3),
            frag1_smoothed_intensity = Float32[80, 60, 40],
            frag2_smoothed_intensity = Float32[20, 40, 60],
            frag3_smoothed_intensity = zeros(Float32, 3),
            frag4_smoothed_intensity = zeros(Float32, 3),
            frag5_smoothed_intensity = zeros(Float32, 3),
            frag6_smoothed_intensity = zeros(Float32, 3),
            frag7_smoothed_intensity = zeros(Float32, 3),
            frag8_smoothed_intensity = zeros(Float32, 3),
        )

        Arrow.write(receiver_path, receiver_df)
        Arrow.write(other_path, other_df)
        for (path, df) in ((receiver_path, receiver_df), (other_path, other_df))
            Arrow.write(path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = df.precursor_idx,
                scan_idx = df.scan_idx,
                trace_prob_prepass = df.trace_prob_prepass,
                trace_prob_infold = df.trace_prob_infold,
            ))
        end

        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1([other_path])
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [receiver_path, other_path],
            MBRPostQMockPrecursors(
                Float32[500.0, 520.0, 500.2, 500.1],
                Float32[10.0, 10.01, 10.50, 12.0],
            ),
        )
        eligible_by_file = Pioneer.build_counterfactual_receiver_eligibility(
            [receiver_path, other_path];
            q_value_threshold = 0.01f0,
        )

        side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            Pioneer.build_mbr_fragment_annotation_keys(_mbr_post_q_fragment_lookup()),
            counterfactual_eligibility_by_file = eligible_by_file,
            irt_tolerance_by_file = Dict(UInt32(1) => 0.75f0),
            passing_score_floor = 0.0f0,
        )
        side = DataFrame(Arrow.Table(side_path))

        @test side.MBR_best_pair_prob_false[1] == 0.80f0
        @test !side.MBR_best_is_missing_false[1]
    end
end

@testset "MBR recovered rows remap precursor probabilities from transfer score" begin
    df = DataFrame(
        prec_prob = Float32[0.90, 0.80, 0.70],
        mbr_recovered = Bool[true, true, false],
        mbr_target_decoy_prob = Float32[0.25, 0.75, NaN32],
    )

    qval_spline(score) = score >= 0.60f0 ? 0.005f0 : 0.020f0
    _, op = Pioneer.remap_mbr_recovered_prec_probs(qval_spline, 0.01f0)
    op(df)

    @test df.prec_prob[1] < df.prec_prob[2]
    @test df.prec_prob[2] < 0.60f0
    @test df.prec_prob[3] == 0.70f0
end

@testset "MBR recovery cutoff uses FTR q-value rather than PEP" begin
    qvals = Float32[0.009, 0.02, 0.001]
    peps = Float32[0.5, 0.0, 0.2]

    @test Pioneer._mbr_recovery_mask(qvals, peps, 0.01f0) == Bool[true, false, true]
end

@testset "MBR semi-supervised transfer training uses counterfactual FTR labels" begin
    @test Pioneer.MBR_SEMISUPERVISED_FTR_THRESHOLD == 0.03f0
    @test !isdefined(Pioneer, :MBR_SEMISUPERVISED_FDR_THRESHOLD)

    n_targets = 50
    top_scores = collect(Float32.(range(1.0, 0.9; length = n_targets)))
    top_scores[end] = 0.895f0
    bottom_scores = vcat(Float32[0.897], collect(Float32.(range(0.5, 0.1; length = n_targets - 1))))
    metrics = Pioneer._mbr_transfer_iteration_metrics(
        vcat(top_scores, bottom_scores),
        n_targets,
    )

    @test metrics.qvals_top[end] <= Pioneer.MBR_SEMISUPERVISED_FTR_THRESHOLD
    @test metrics.pep_top[end] > Pioneer.MBR_SEMISUPERVISED_FTR_THRESHOLD
    @test metrics.positive_top[end]

    scores = Float32[0.99, 0.98, 0.97, 0.96, 0.50, 0.40, 0.30, 0.20]
    target_top = Bool[true, false, true, true]
    metrics = Pioneer._mbr_transfer_iteration_metrics(scores, length(target_top))

    @test metrics.qvals_top == zeros(Float32, length(target_top))
    @test metrics.positive_top == trues(length(target_top))

    positive_top = Bool[true, false, false, false]

    labels = Pioneer._mbr_transfer_training_labels(positive_top)
    mask = Pioneer._mbr_transfer_training_mask(positive_top)

    @test labels == Bool[true, false, false, false, false, false, false, false]
    @test mask == Bool[true, false, false, false, true, true, true, true]

    valid_transfer_top = Bool[true, false, true, false]
    valid_labels = Pioneer._mbr_transfer_training_labels(valid_transfer_top)
    valid_mask = Pioneer._mbr_transfer_training_mask(valid_transfer_top)
    valid_metrics = Pioneer._mbr_transfer_iteration_metrics(
        scores,
        length(valid_transfer_top);
        eval_labels = valid_labels,
        eval_mask = valid_mask,
    )

    @test valid_metrics.positive_top == valid_transfer_top
    @test valid_metrics.n_positive == count(valid_transfer_top)
    @test valid_metrics.qvals_top[1] == 0.0f0
    @test isinf(valid_metrics.qvals_top[2])
    @test valid_metrics.qvals_top[3] == 0.0f0
    @test isinf(valid_metrics.qvals_top[4])
end

@testset "MBR FTR model uses iRT difference but not raw RT difference" begin
    expected_true = Symbol[
        :MBR_best_pair_prob_true,
        :MBR_worst_pair_prob_true,
        :MBR_best_log2_weight_ratio_true,
        :MBR_worst_log2_weight_ratio_true,
        :MBR_best_log2_explained_ratio_true,
        :MBR_worst_log2_explained_ratio_true,
        :MBR_best_abs_n_scans_diff_true,
        :MBR_worst_abs_n_scans_diff_true,
        :MBR_best_log2_n_scans_ratio_true,
        :MBR_worst_log2_n_scans_ratio_true,
        :MBR_best_irt_diff_true,
        :MBR_worst_irt_diff_true,
        :MBR_best_observed_irt_diff_true,
        :MBR_worst_observed_irt_diff_true,
        :MBR_single_donor_true,
        :MBR_best_smoothed_frag_hellinger_true,
        :MBR_worst_smoothed_frag_hellinger_true,
        :MBR_best_donor_library_hellinger_true,
        :MBR_worst_donor_library_hellinger_true,
    ]

    for feature in expected_true
        @test feature in Pioneer.FTR_FEATURES_F_TRUE
    end

    @test :MBR_best_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE
    @test !(:MBR_max_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_max_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test :MBR_best_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_observed_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_observed_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_observed_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_observed_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_single_donor_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_single_donor_false in Pioneer.FTR_FEATURES_F_FALSE
    @test !(:MBR_best_log_by_diff_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_best_log_by_diff_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_best_rt_diff_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_best_rt_diff_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_log2_weight_ratio_worst_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_smoothed_frag_hellinger_worst_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_donor_library_hellinger_worst_true in Pioneer.FTR_FEATURES_F_TRUE)
end

@testset "MBR observed iRT difference and single-donor flag use best/worst donors" begin
    mktempdir() do dir
        best_donor_path = joinpath(dir, "best_donor.arrow")
        worst_donor_path = joinpath(dir, "worst_donor.arrow")
        single_donor_path = joinpath(dir, "single_donor.arrow")
        receiver_path = joinpath(dir, "receiver.arrow")
        single_receiver_path = joinpath(dir, "single_receiver.arrow")

        best_donor_df = _mbr_post_q_main_table(ms_file_idx = 1, scan_offset = 1000)
        best_donor_df[!, :irt_obs] = Float32[11.5, 19.0]
        worst_donor_df = _mbr_post_q_main_table(ms_file_idx = 2, scan_offset = 2000)
        worst_donor_df[!, :irt_obs] = Float32[17.0, 24.0]
        worst_donor_df[!, :n_scans] = Float32[5, 4]
        worst_donor_df[!, :trace_prob_prepass] = Float32[0.70, 0.60]
        worst_donor_df[!, :trace_prob_infold] = Float32[0.70, 0.60]
        single_donor_df = _mbr_post_q_main_table(ms_file_idx = 4, scan_offset = 4000)
        single_donor_df[!, :irt_obs] = Float32[30.0, 35.0]
        receiver_df = _mbr_post_q_main_table(ms_file_idx = 3, scan_offset = 3000)
        receiver_df[!, :irt_obs] = Float32[14.25, 22.0]
        single_receiver_df = _mbr_post_q_main_table(ms_file_idx = 5, scan_offset = 5000)
        single_receiver_df[!, :irt_obs] = Float32[31.25, 37.0]

        Arrow.write(best_donor_path, best_donor_df)
        Arrow.write(worst_donor_path, worst_donor_df)
        Arrow.write(single_donor_path, single_donor_df)
        Arrow.write(receiver_path, receiver_df)
        Arrow.write(single_receiver_path, single_receiver_df)
        for (path, df) in (
            (best_donor_path, best_donor_df),
            (worst_donor_path, worst_donor_df),
            (single_donor_path, single_donor_df),
            (receiver_path, receiver_df),
            (single_receiver_path, single_receiver_df),
        )
            Arrow.write(path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = df.precursor_idx,
                scan_idx = df.scan_idx,
                trace_prob_prepass = df.trace_prob_prepass,
                trace_prob_infold = df.trace_prob_infold,
            ))
        end

        precursors = MBRPostQMockPrecursors(Float32[500, 600], Float32[10, 20])
        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1(
            [best_donor_path, worst_donor_path],
        )
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [best_donor_path, worst_donor_path, receiver_path],
            precursors,
        )

        side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            Pioneer.build_mbr_fragment_annotation_keys(_mbr_post_q_fragment_lookup()),
            passing_score_floor = 0.0f0,
        )
        side = DataFrame(Arrow.Table(side_path))

        receiver_residual = receiver_df.irt_pred[1] - receiver_df.irt_obs[1]
        best_donor_residual = best_donor_df.irt_pred[1] - best_donor_df.irt_obs[1]
        worst_donor_residual = worst_donor_df.irt_pred[1] - worst_donor_df.irt_obs[1]

        @test side.MBR_best_observed_irt_diff_true[1] == abs(14.25f0 - 11.5f0)
        @test side.MBR_worst_observed_irt_diff_true[1] == abs(14.25f0 - 17.0f0)
        @test side.MBR_best_irt_diff_true[1] == abs(receiver_residual - best_donor_residual)
        @test side.MBR_worst_irt_diff_true[1] == abs(receiver_residual - worst_donor_residual)
        @test isapprox(
            side.MBR_worst_log2_n_scans_ratio_true[1],
            log2((receiver_df.n_scans[1] + 1.0f0) / (worst_donor_df.n_scans[1] + 1.0f0));
            atol = 1.0f-6,
        )
        @test side.MBR_single_donor_true[1] == 0.0f0
        @test :MBR_worst_irt_diff_true in propertynames(side)
        @test :MBR_worst_irt_diff_false in propertynames(side)
        @test :MBR_worst_log2_n_scans_ratio_true in propertynames(side)
        @test :MBR_worst_log2_n_scans_ratio_false in propertynames(side)
        @test !(:MBR_best_log_by_diff_true in propertynames(side))
        @test !(:MBR_best_log_by_diff_false in propertynames(side))

        single_donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1([single_donor_path])
        single_partner_pools = Pioneer.build_counterfactual_partner_pools(
            [single_donor_path, single_receiver_path],
            precursors,
        )
        single_side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            single_receiver_path,
            single_donor_dict,
            single_partner_pools,
            Pioneer.build_mbr_fragment_annotation_keys(_mbr_post_q_fragment_lookup()),
            passing_score_floor = 0.0f0,
        )
        single_side = DataFrame(Arrow.Table(single_side_path))

        @test single_side.MBR_single_donor_true[1] == 1.0f0
        @test single_side.MBR_worst_observed_irt_diff_true[1] == -1.0f0
        @test single_side.MBR_worst_irt_diff_true[1] == -1.0f0
        @test single_side.MBR_single_donor_false[1] == 1.0f0
        @test single_side.MBR_worst_observed_irt_diff_false[1] == -1.0f0
        @test single_side.MBR_worst_irt_diff_false[1] == -1.0f0
    end
end
