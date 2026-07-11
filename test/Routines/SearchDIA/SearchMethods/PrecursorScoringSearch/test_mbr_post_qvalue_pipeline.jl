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
        Pioneer.CompactFrag(
            UInt32(5), 900.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(7), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(6), 1000.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(8), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(7), 1100.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(9), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(8), 1200.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(10), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(9), 1300.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(11), UInt8(2), UInt8(1), UInt8(0),
        ),
    ]
    return Pioneer.StandardFragmentLookup{Float32}(frags, UInt64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
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

@testset "post-qvalue MBR counterfactuals can pair to receiver-run passing precursors" begin
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
        donor_df = copy(receiver_df)
        donor_df[!, :ms_file_idx] .= UInt32(2)
        donor_df[!, :scan_idx] .= UInt32[201, 202, 203]
        donor_df[!, :trace_prob] .= Float32[0.80, 0.99, 0.70]
        donor_df[!, :trace_prob_prepass] .= Float32[0.80, 0.99, 0.70]
        donor_df[!, :trace_prob_infold] .= Float32[0.80, 0.99, 0.70]
        donor_df[!, :prec_prob] .= Float32[0.80, 0.99, 0.70]

        Arrow.write(receiver_path, receiver_df)
        Arrow.write(donor_path, donor_df)

        old_keep_sidecars = get(ENV, "PIONEER_MBR_KEEP_SIDECARS", nothing)
        old_cf = get(ENV, "PIONEER_MBR_N_COUNTERFACTUALS", nothing)
        ENV["PIONEER_MBR_KEEP_SIDECARS"] = "1"
        ENV["PIONEER_MBR_N_COUNTERFACTUALS"] = "1"
        try
            Pioneer.run_mbr_after_qvalue_filter!(
                [Pioneer.PSMFileReference(receiver_path)],
                [Pioneer.PSMFileReference(donor_path)],
                MBRPostQMockPrecursors(Float32[500, 501, 502], Float32[10.0, 10.02, 10.20]),
                _mbr_post_q_fragment_lookup(),
            )
            side = DataFrame(Arrow.Table(receiver_path * Pioneer.MBR_SIDECAR_SUFFIX))

            @test side.MBR_best_pair_prob_false[1] == 0.99f0
            @test !side.MBR_best_is_missing_false[1]
        finally
            if old_keep_sidecars === nothing
                delete!(ENV, "PIONEER_MBR_KEEP_SIDECARS")
            else
                ENV["PIONEER_MBR_KEEP_SIDECARS"] = old_keep_sidecars
            end
            if old_cf === nothing
                delete!(ENV, "PIONEER_MBR_N_COUNTERFACTUALS")
            else
                ENV["PIONEER_MBR_N_COUNTERFACTUALS"] = old_cf
            end
        end
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
            precursor_idx = UInt32[1, 2, 3],
            scan_idx = UInt32[201, 202, 203],
            weight = Float32[10, 2, 3],
            log2_intensity_explained = Float16[1.0, 0.5, 0.5],
            irt_pred = Float32[10.0, 10.02, 10.20],
            irt_obs = Float32[10.1, 10.1, 10.1],
            log_by_ratio_m0 = zeros(Float16, 3),
            rt = Float32[5.0, 5.1, 5.2],
            n_scans = Float32[8, 4, 5],
            main_pep = Float32[0.001, 0.001, 0.001],
            smoothed_2d_shadow_hellinger = Float32[0.1, 0.2, 0.3],
            ms_file_idx = fill(UInt32(2), 3),
            target = trues(3),
            cv_fold = zeros(UInt8, 3),
            trace_prob = Float32[0.80, 0.99, 0.70],
            trace_prob_prepass = Float32[0.80, 0.99, 0.70],
            trace_prob_infold = Float32[0.80, 0.99, 0.70],
            prec_prob = Float32[0.80, 0.99, 0.70],
            qval = Float32[0.001, 0.001, 0.001],
            pep = Float32[0.001, 0.001, 0.001],
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

@testset "MBR counterfactuals require same-charge same-length precursors" begin
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

        @test side.MBR_best_pair_prob_false[1] == -1.0f0
        @test side.MBR_best_is_missing_false[1]
    end
end

@testset "MBR counterfactuals choose closest library iRT within same length and charge" begin
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
            precursor_idx = UInt32[1, 2, 3, 4, 5, 6, 7, 8, 9],
            scan_idx = UInt32[201, 202, 203, 204, 205, 206, 207, 208, 209],
            weight = Float32[10, 2, 3, 4, 5, 6, 7, 8, 9],
            log2_intensity_explained = Float16[1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            irt_pred = Float32[10.0, 10.01, 10.50, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00],
            irt_obs = fill(10.1f0, 9),
            log_by_ratio_m0 = zeros(Float16, 9),
            rt = Float32[5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8],
            n_scans = Float32[8, 4, 5, 6, 7, 8, 9, 10, 11],
            main_pep = fill(0.001f0, 9),
            smoothed_2d_shadow_hellinger = Float32[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            ms_file_idx = fill(UInt32(2), 9),
            target = trues(9),
            cv_fold = zeros(UInt8, 9),
            trace_prob = Float32[0.90, 0.60, 0.80, 0.95, 0.70, 0.65, 0.55, 0.45, 0.35],
            trace_prob_prepass = Float32[0.90, 0.60, 0.80, 0.95, 0.70, 0.65, 0.55, 0.45, 0.35],
            trace_prob_infold = Float32[0.90, 0.60, 0.80, 0.95, 0.70, 0.65, 0.55, 0.45, 0.35],
            prec_prob = Float32[0.90, 0.60, 0.80, 0.95, 0.70, 0.65, 0.55, 0.45, 0.35],
            qval = fill(0.001f0, 9),
            pep = fill(0.001f0, 9),
            global_qval = fill(0.001f0, 9),
            global_pep = fill(0.001f0, 9),
            mbr_recovered = falses(9),
            frag1_smoothed_intensity = Float32[100, 80, 60, 40, 20, 10, 10, 10, 10],
            frag2_smoothed_intensity = Float32[0, 20, 40, 60, 80, 90, 90, 90, 90],
            frag3_smoothed_intensity = zeros(Float32, 9),
            frag4_smoothed_intensity = zeros(Float32, 9),
            frag5_smoothed_intensity = zeros(Float32, 9),
            frag6_smoothed_intensity = zeros(Float32, 9),
            frag7_smoothed_intensity = zeros(Float32, 9),
            frag8_smoothed_intensity = zeros(Float32, 9),
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
                Float32[500.0, 650.0, 500.2, 500.1, 500.3, 500.4, 500.5, 500.6, 500.7],
                Float32[10.0, 10.01, 10.50, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0],
                UInt8[9, 9, 9, 9, 9, 9, 9, 9, 9],
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

        @test side.MBR_best_pair_prob_false[1] == 0.60f0
        @test side.MBR_best_pair_prob_false2[1] == 0.80f0
        @test side.MBR_best_pair_prob_false3[1] == 0.95f0
        @test side.MBR_best_pair_prob_false4[1] == 0.70f0
        @test side.MBR_best_pair_prob_false5[1] == 0.65f0
        @test side.MBR_best_pair_prob_false6[1] == 0.55f0
        @test side.MBR_best_pair_prob_false7[1] == 0.45f0
        @test side.MBR_best_pair_prob_false8[1] == 0.35f0
        @test !side.MBR_best_is_missing_false[1]
        @test !side.MBR_best_is_missing_false2[1]
        @test !side.MBR_best_is_missing_false3[1]
        @test !side.MBR_best_is_missing_false4[1]
        @test !side.MBR_best_is_missing_false5[1]
        @test !side.MBR_best_is_missing_false6[1]
        @test !side.MBR_best_is_missing_false7[1]
        @test !side.MBR_best_is_missing_false8[1]
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

@testset "MBR combined error cutoff spends unused global FDR budget" begin
    slack_result = Pioneer._mbr_combined_error_recovery(
        Float32[0.95, 0.90, 0.85, 0.80],
        Bool[true, true, true, true],
        Float32[0.10, 0.89, 0.84, 0.80];
        base_targets = 100,
        base_decoys = 0,
        alpha = 0.01f0,
    )

    @test slack_result.recovered == Bool[true, true, true, false]
    @test slack_result.n_recovered == 3
    @test slack_result.mbr_false_transfers == 1
    @test slack_result.mbr_decoys == 0
    @test slack_result.combined_error_rate <= 0.01f0
    @test slack_result.combined_error_qvals[3] <= 0.01f0
    @test slack_result.combined_error_qvals[4] > 0.01f0

    decoy_result = Pioneer._mbr_combined_error_recovery(
        Float32[0.95, 0.94],
        Bool[false, true],
        Float32[0.95, 0.10];
        base_targets = 100,
        base_decoys = 0,
        alpha = 0.01f0,
    )

    @test decoy_result.recovered == Bool[true, true]
    @test decoy_result.mbr_decoys == 1
    @test decoy_result.mbr_false_transfers == 0
    @test decoy_result.total_errors == 1

    no_dilution_result = Pioneer._mbr_combined_error_recovery(
        fill(0.99f0, 200),
        trues(200),
        fill(-Inf32, 200);
        base_targets = 100,
        base_decoys = 2,
        alpha = 0.01f0,
    )

    @test no_dilution_result.n_recovered == 0
    @test all(>(0.01f0), no_dilution_result.combined_error_qvals)
end

@testset "MBR candidates do not require an available counterfactual" begin
    old_cf = get(ENV, "PIONEER_MBR_N_COUNTERFACTUALS", nothing)
    old_disable_hellinger = get(ENV, "PIONEER_MBR_DISABLE_HELLINGER_CONTRAST", nothing)
    ENV["PIONEER_MBR_N_COUNTERFACTUALS"] = "1"
    ENV["PIONEER_MBR_DISABLE_HELLINGER_CONTRAST"] = "1"
    try
        psms = DataFrame(
            precursor_idx = UInt32[10, 20],
            scan_idx = UInt32[100, 200],
            ms_file_idx = UInt32[1, 2],
            cv_fold = UInt8[0, 1],
            target = Bool[true, true],
            qval = Float32[0.02, 0.02],
            global_qval = Float32[0.001, 0.001],
            trace_prob_prepass = Float32[0.95, 0.96],
            trace_prob_infold = Float32[0.80, 0.81],
            MBR_best_pair_prob_true = Float32[0.90, 0.91],
            MBR_best_pair_prob_false = Float32[-1.0, -1.0],
            MBR_best_is_missing_true = Bool[false, false],
            MBR_best_is_missing_false = Bool[true, true],
        )

        summary = Pioneer.apply_mbr_filter_paired!(
            psms;
            alpha = 0.01f0,
            q_thresh = 0.01f0,
            prob_thresh_override = 0.0f0,
        )

        @test summary.n_candidates == 2
        @test psms.MBR_transfer_candidate == Bool[true, true]
        @test all(isfinite, psms.ftr_qval_true)
    finally
        if old_cf === nothing
            delete!(ENV, "PIONEER_MBR_N_COUNTERFACTUALS")
        else
            ENV["PIONEER_MBR_N_COUNTERFACTUALS"] = old_cf
        end
        if old_disable_hellinger === nothing
            delete!(ENV, "PIONEER_MBR_DISABLE_HELLINGER_CONTRAST")
        else
            ENV["PIONEER_MBR_DISABLE_HELLINGER_CONTRAST"] = old_disable_hellinger
        end
    end
end

@testset "MBR recovery alpha can be overridden for threshold sweeps" begin
    @test Pioneer._mbr_recovery_alpha_from_env(Dict{String, String}()) == 0.01f0
    @test Pioneer._mbr_recovery_alpha_from_env(Dict("PIONEER_MBR_FTR_ALPHA" => "0.02")) == 0.02f0
    @test Pioneer._mbr_recovery_alpha_from_env(Dict("PIONEER_MBR_FTR_ALPHA" => "0.05")) == 0.05f0
    @test_throws ErrorException Pioneer._mbr_recovery_alpha_from_env(Dict("PIONEER_MBR_FTR_ALPHA" => "0"))
    @test_throws ErrorException Pioneer._mbr_recovery_alpha_from_env(Dict("PIONEER_MBR_FTR_ALPHA" => "1.1"))
    @test_throws ErrorException Pioneer._mbr_recovery_alpha_from_env(Dict("PIONEER_MBR_FTR_ALPHA" => "nope"))
end

@testset "MBR counterfactual count can be selected for sweeps" begin
    @test Pioneer._mbr_n_counterfactuals_from_env(Dict{String, String}()) == 2
    @test Pioneer._mbr_n_counterfactuals_from_env(Dict("PIONEER_MBR_N_COUNTERFACTUALS" => "3")) == 3
    @test Pioneer._mbr_n_counterfactuals_from_env(Dict("PIONEER_MBR_N_COUNTERFACTUALS" => "4")) == 4
    @test Pioneer._mbr_n_counterfactuals_from_env(Dict("PIONEER_MBR_N_COUNTERFACTUALS" => "8")) == 8
    @test_throws ErrorException Pioneer._mbr_n_counterfactuals_from_env(Dict("PIONEER_MBR_N_COUNTERFACTUALS" => "0"))
    @test_throws ErrorException Pioneer._mbr_n_counterfactuals_from_env(Dict("PIONEER_MBR_N_COUNTERFACTUALS" => "9"))
    @test_throws ErrorException Pioneer._mbr_n_counterfactuals_from_env(Dict("PIONEER_MBR_N_COUNTERFACTUALS" => "nope"))
end

@testset "MBR partial counterfactual coverage keeps candidates and masks missing blocks" begin
    psms = DataFrame(
        MBR_best_is_missing_true = Bool[false, false, false, false],
        MBR_best_is_missing_false = Bool[false, true, false, true],
        MBR_best_is_missing_false2 = Bool[true, false, true, true],
        MBR_best_is_missing_false3 = Bool[true, true, true, true],
    )
    present = Pioneer._mbr_counterfactual_present_matrix(
        psms,
        [1, 2, 3, 4];
        n_counterfactuals = 3,
    )
    @test present == Bool[
        true false false
        false true false
        true false false
        false false false
    ]

    frame_mask = Pioneer._mbr_transfer_frame_mask(present)
    @test frame_mask == Bool[
        true, true, true, true,
        true, false, true, false,
        false, true, false, false,
        false, false, false, false,
    ]
    @test frame_mask[4]

    training_mask = Pioneer._mbr_transfer_training_mask(
        Bool[true, false, true, true];
        n_counterfactuals = 3,
        frame_mask = frame_mask,
    )
    @test training_mask == Bool[
        true, false, true, true,
        true, false, true, false,
        false, true, false, false,
        false, false, false, false,
    ]
    @test training_mask[4]
end

@testset "MBR FTR q-values use the top scoring counterfactual per candidate" begin
    @test Pioneer._mbr_use_top_counterfactual_qvalues_from_env(Dict{String,String}())
    @test !Pioneer._mbr_use_top_counterfactual_qvalues_from_env(Dict(
        "PIONEER_MBR_USE_TOP_COUNTERFACTUAL_QVALUES" => "0",
    ))

    scores = Float32[
        0.50, 0.50, 0.50,
        0.10, 0.90, 0.30,
        0.80, 0.20, 0.40,
    ]
    base_mask = Bool[
        true, true, true,
        true, true, false,
        true, true, true,
    ]
    top_cf_mask = Pioneer._mbr_top_counterfactual_eval_mask(
        scores,
        3;
        n_counterfactuals = 2,
        eval_mask = base_mask,
    )

    @test top_cf_mask == Bool[
        true, true, true,
        false, true, false,
        true, false, true,
    ]

    hard_cf_scores = Float32[
        0.90, 0.80,
        0.95, 0.10,
        0.20, 0.70,
    ]
    metrics = Pioneer._mbr_transfer_iteration_metrics(
        hard_cf_scores,
        2;
        n_counterfactuals = 2,
    )

    @test metrics.qvalue_eval_mask == Bool[
        true, true,
        true, false,
        false, true,
    ]
    @test metrics.qvals_top[1] ≈ 0.5f0
    @test metrics.qvals_top[2] ≈ 0.5f0
    @test isinf(metrics.qvals_double[4])
    @test isinf(metrics.qvals_double[5])

    all_cf_metrics = Pioneer._mbr_transfer_iteration_metrics(
        hard_cf_scores,
        2;
        n_counterfactuals = 2,
        use_top_counterfactual_qvalues = false,
    )
    @test all_cf_metrics.qvalue_eval_mask == trues(6)
    @test all(isfinite, all_cf_metrics.qvals_double)
end

@testset "MBR initial positives use candidate-local receiver q-values" begin
    @test Pioneer.MBR_INITIAL_RECEIVER_Q_THRESHOLD == 0.05f0

    receiver_scores = Float32[
        0.99, 0.98, 0.97, 0.96, 0.95, 0.94,
        0.93, 0.92, 0.91, 0.90, 0.89, 0.88,
    ]
    target_top = Bool[
        true, true, true, true, true, true,
        true, true, true, true, false, true,
    ]

    positive_top, receiver_qvals = Pioneer._mbr_initial_receiver_positive_top(
        receiver_scores,
        target_top,
    )

    @test positive_top == Bool[
        true, true, true, true, true, true,
        true, true, true, true, false, false,
    ]
    @test receiver_qvals[1] == 0.0f0
    @test receiver_qvals[10] <= 0.05f0
    @test receiver_qvals[11] > 0.05f0
    @test receiver_qvals[12] > 0.05f0
    @test !positive_top[11]
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

    two_counterfactual_labels = Pioneer._mbr_transfer_counterfactual_labels(
        length(positive_top);
        n_counterfactuals = 2,
    )
    two_counterfactual_training_labels = Pioneer._mbr_transfer_training_labels(
        positive_top;
        n_counterfactuals = 2,
    )
    two_counterfactual_training_mask = Pioneer._mbr_transfer_training_mask(
        positive_top;
        n_counterfactuals = 2,
    )

    @test two_counterfactual_labels == Bool[
        true, true, true, true,
        false, false, false, false,
        false, false, false, false,
    ]
    @test two_counterfactual_training_labels == Bool[
        true, false, false, false,
        false, false, false, false,
        false, false, false, false,
    ]
    @test two_counterfactual_training_mask == Bool[
        true, false, false, false,
        true, true, true, true,
        true, true, true, true,
    ]

    four_counterfactual_labels = Pioneer._mbr_transfer_counterfactual_labels(
        length(positive_top);
        n_counterfactuals = 4,
    )
    four_counterfactual_training_mask = Pioneer._mbr_transfer_training_mask(
        positive_top;
        n_counterfactuals = 4,
    )
    @test length(four_counterfactual_labels) == 20
    @test four_counterfactual_labels[1:4] == trues(4)
    @test all(.!four_counterfactual_labels[5:20])
    @test four_counterfactual_training_mask == Bool[
        true, false, false, false,
        true, true, true, true,
        true, true, true, true,
        true, true, true, true,
        true, true, true, true,
    ]

    two_cf_scores = Float32[0.97, 0.98, 0.10]
    two_cf_training_metrics = Pioneer._mbr_transfer_iteration_metrics(
        two_cf_scores,
        1;
        n_counterfactuals = 2,
    )

    @test_throws MethodError Pioneer._mbr_transfer_iteration_metrics(
        two_cf_scores,
        1;
        n_counterfactuals = 2,
        scale_qvalues_by_counterfactuals = true,
    )

    @test two_cf_training_metrics.qvals_top[1] ≈ 1.0f0

    partial_cf_mask = Bool[
        true,
        true,
        false,
    ]
    partial_cf_training_metrics = Pioneer._mbr_transfer_iteration_metrics(
        two_cf_scores,
        1;
        n_counterfactuals = 2,
        eval_mask = partial_cf_mask,
    )

    @test partial_cf_training_metrics.qvals_top[1] ≈ 1.0f0

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

@testset "MBR Hellinger rank uses only selected counterfactuals" begin
    psms = DataFrame(
        MBR_best_smoothed_frag_hellinger_true = Float32[0.20, 0.80],
        MBR_best_smoothed_frag_hellinger_false = Float32[0.60, 0.10],
        MBR_best_smoothed_frag_hellinger_false2 = Float32[0.40, 0.90],
        MBR_best_smoothed_frag_hellinger_false3 = Float32[0.05, 0.20],
    )

    Pioneer._mbr_add_best_hellinger_rank_features!(
        psms;
        n_counterfactuals = 2,
    )

    @test psms.MBR_best_smoothed_frag_hellinger_rank_true == Float32[1, 2]
    @test psms.MBR_best_smoothed_frag_hellinger_rank_false == Float32[3, 1]
    @test psms.MBR_best_smoothed_frag_hellinger_rank_false2 == Float32[2, 3]
    @test !(:MBR_best_smoothed_frag_hellinger_rank_false3 in propertynames(psms))
end

@testset "MBR FTR debug export writes selected true and counterfactual blocks" begin
    mktempdir() do dir
        sub = DataFrame(
            precursor_idx = UInt32[10, 20],
            ms_file_idx = UInt32[1, 2],
            scan_idx = UInt32[100, 200],
            cv_fold = UInt8[0, 1],
            target = Bool[true, true],
            qval = Float32[0.02, 0.03],
            global_qval = Float32[0.001, 0.001],
            trace_prob_prepass = Float32[0.50, 0.60],
            trace_prob_infold = Float32[0.55, 0.65],
            MBR_best_pair_prob_true = Float32[0.90, 0.85],
            MBR_best_pair_prob_false = Float32[0.40, 0.30],
        )
        available_true = Symbol[:trace_prob_infold, :MBR_best_pair_prob_true]
        x_blocks = Matrix{Float32}[
            Float32[0.55 0.90; 0.65 0.85],
            Float32[0.55 0.40; 0.65 0.30],
        ]

        Pioneer._mbr_write_ftr_debug_tables!(
            dir,
            sub,
            available_true,
            x_blocks,
            Float32[0.99, 0.80, 0.20, 0.10],
            Float32[0.00, 0.01, 1.00, 1.00],
            Float32[0.00, 0.02, 0.90, 0.95],
            Bool[true, true, false, false],
            Bool[true, true, true, true],
            Bool[true, false];
            n_counterfactuals = 1,
            alpha = 0.01f0,
            q_thresh = 0.01f0,
            prob_thresh = 0.5f0,
            best_iter = 1,
            n_positive = 2,
        )

        frame = DataFrame(Arrow.Table(joinpath(dir, "mbr_ftr_frame.arrow")))
        @test nrow(frame) == 4
        @test frame.mbr_counterfactual_idx == UInt8[0, 0, 1, 1]
        @test frame.mbr_is_true_block == Bool[true, true, false, false]
        @test frame.MBR_best_pair_prob_true == Float32[0.90, 0.85, 0.40, 0.30]
        @test frame.mbr_recovered_top == Bool[true, false, false, false]

        candidates = DataFrame(Arrow.Table(joinpath(dir, "mbr_ftr_candidates.arrow")))
        @test candidates.mbr_candidate_idx == UInt32[1, 2]
        @test candidates.mbr_ftr_qval_true_debug == Float32[0.00, 0.01]
        @test candidates.mbr_recovered_debug == Bool[true, false]

        summary = DataFrame(Arrow.Table(joinpath(dir, "mbr_ftr_summary.arrow")))
        @test only(summary.n_counterfactuals) == 1
        @test only(summary.n_recovered) == 1
        @test read(joinpath(dir, "mbr_ftr_features.txt"), String) ==
            "trace_prob_infold\nMBR_best_pair_prob_true\n"
    end
end

@testset "MBR FTR model uses cross-run, best, and worst donor features while ablating RT, b/y, and non-best Hellinger" begin
    expected_true = Symbol[
        Pioneer.MBR_CROSS_RUN_FTR_FEATURES...,
        :MBR_best_pair_prob_true,
        :MBR_worst_pair_prob_true,
        :MBR_best_run_similarity_true,
        :MBR_log2_weight_lod_ratio,
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
        :MBR_best_smoothed_frag_hellinger_rank_true,
        :MBR_best_corr_frag_hellinger_true,
        :MBR_best_corr_frag_hellinger_rank_true,
        :MBR_best_donor_frag_corr_bitvec_rank_true,
        :MBR_best_receiver_corr_frag_hellinger_true,
        :MBR_best_receiver_corr_frag_hellinger_rank_true,
        :MBR_receiver_frag_corr_bitvec_rank,
        :MBR_best_shared_corr_frag_hellinger_true,
        :MBR_best_shared_corr_frag_hellinger_rank_true,
        :MBR_best_shared_corr_frag_bitvec_rank_true,
    ]
    expected_false = Symbol[
        Pioneer.MBR_CROSS_RUN_FTR_FEATURES...,
        :MBR_best_pair_prob_false,
        :MBR_worst_pair_prob_false,
        :MBR_best_run_similarity_false,
        :MBR_log2_weight_lod_ratio,
        :MBR_best_log2_weight_ratio_false,
        :MBR_worst_log2_weight_ratio_false,
        :MBR_best_log2_explained_ratio_false,
        :MBR_worst_log2_explained_ratio_false,
        :MBR_best_abs_n_scans_diff_false,
        :MBR_worst_abs_n_scans_diff_false,
        :MBR_best_log2_n_scans_ratio_false,
        :MBR_worst_log2_n_scans_ratio_false,
        :MBR_best_irt_diff_false,
        :MBR_worst_irt_diff_false,
        :MBR_best_observed_irt_diff_false,
        :MBR_worst_observed_irt_diff_false,
        :MBR_single_donor_false,
        :MBR_best_smoothed_frag_hellinger_false,
        :MBR_best_smoothed_frag_hellinger_rank_false,
        :MBR_best_corr_frag_hellinger_false,
        :MBR_best_corr_frag_hellinger_rank_false,
        :MBR_best_donor_frag_corr_bitvec_rank_false,
        :MBR_best_receiver_corr_frag_hellinger_false,
        :MBR_best_receiver_corr_frag_hellinger_rank_false,
        :MBR_receiver_frag_corr_bitvec_rank,
        :MBR_best_shared_corr_frag_hellinger_false,
        :MBR_best_shared_corr_frag_hellinger_rank_false,
        :MBR_best_shared_corr_frag_bitvec_rank_false,
    ]

    @test Pioneer.FTR_FEATURES_F_TRUE == expected_true
    @test Pioneer.FTR_FEATURES_F_FALSE == expected_false
    @test Pioneer._mbr_ftr_features_true_from_env(Dict{String, String}()) == expected_true
    no_contrast_features = Pioneer._mbr_ftr_features_true_from_env(
        Dict("PIONEER_MBR_DISABLE_HELLINGER_CONTRAST" => "1"),
    )
    @test :MBR_best_smoothed_frag_hellinger_true in no_contrast_features
    @test :MBR_best_corr_frag_hellinger_true in no_contrast_features
    @test :MBR_best_donor_frag_corr_bitvec_rank_true in no_contrast_features
    @test :MBR_best_receiver_corr_frag_hellinger_true in no_contrast_features
    @test :MBR_receiver_frag_corr_bitvec_rank in no_contrast_features
    @test :MBR_best_shared_corr_frag_hellinger_true in no_contrast_features
    @test :MBR_best_shared_corr_frag_bitvec_rank_true in no_contrast_features
    @test !(:MBR_best_smoothed_frag_hellinger_rank_true in no_contrast_features)
    @test !(:MBR_best_corr_frag_hellinger_rank_true in no_contrast_features)
    @test !(:MBR_best_receiver_corr_frag_hellinger_rank_true in no_contrast_features)
    @test !(:MBR_best_shared_corr_frag_hellinger_rank_true in no_contrast_features)
    @test length(no_contrast_features) == length(expected_true) - 4
    @test Pioneer.FTR_FEATURES_F_FALSE2 == Symbol[
        f === :MBR_log2_weight_lod_ratio ? f :
        Symbol(replace(String(f), "_true" => "_false2"))
        for f in expected_true
    ]

    all_features = Set(vcat(Pioneer.FTR_FEATURES_F_TRUE, Pioneer.FTR_FEATURES_F_FALSE))
    @test !(:main_search_prob in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:main_search_prob in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:trace_prob_infold in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:trace_prob_infold in Pioneer.FTR_FEATURES_F_FALSE)
    @test Pioneer.MBR_CROSS_RUN_FTR_FEATURES == Symbol[Pioneer.ADVANCED_FEATURE_SET...]
    @test all(feature -> feature in Pioneer.FTR_FEATURES_F_TRUE, Pioneer.ADVANCED_FEATURE_SET)
    @test all(feature -> feature in Pioneer.FTR_FEATURES_F_FALSE, Pioneer.ADVANCED_FEATURE_SET)
    @test :MBR_best_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_run_similarity_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_run_similarity_false in Pioneer.FTR_FEATURES_F_FALSE
    @test !(:MBR_worst_run_similarity_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_worst_run_similarity_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_median_run_similarity_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_median_run_similarity_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_max_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_max_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test :MBR_best_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test !(:MBR_best_rt_diff_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_best_rt_diff_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_best_log_by_diff_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_best_log_by_diff_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test :MBR_worst_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_observed_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_observed_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_observed_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_observed_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_single_donor_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_single_donor_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_log2_weight_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_log2_weight_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_log2_explained_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_log2_explained_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_abs_n_scans_diff_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_abs_n_scans_diff_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_worst_log2_n_scans_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_worst_log2_n_scans_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
    @test !(:MBR_log2_weight_ratio_worst_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_smoothed_frag_hellinger_worst_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_donor_library_hellinger_worst_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test :MBR_best_smoothed_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_smoothed_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_smoothed_frag_hellinger_rank_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_smoothed_frag_hellinger_rank_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_corr_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_corr_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_corr_frag_hellinger_rank_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_corr_frag_hellinger_rank_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_donor_frag_corr_bitvec_rank_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_donor_frag_corr_bitvec_rank_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_receiver_corr_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_receiver_corr_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_receiver_corr_frag_hellinger_rank_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_receiver_corr_frag_hellinger_rank_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_receiver_frag_corr_bitvec_rank in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_receiver_frag_corr_bitvec_rank in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_shared_corr_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_shared_corr_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_shared_corr_frag_hellinger_rank_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_shared_corr_frag_hellinger_rank_false in Pioneer.FTR_FEATURES_F_FALSE
    @test :MBR_best_shared_corr_frag_bitvec_rank_true in Pioneer.FTR_FEATURES_F_TRUE
    @test :MBR_best_shared_corr_frag_bitvec_rank_false in Pioneer.FTR_FEATURES_F_FALSE
    @test !(:MBR_best_smoothed_frag_hellinger_margin_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_best_smoothed_frag_hellinger_margin_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_worst_smoothed_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_worst_smoothed_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_best_donor_library_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_best_donor_library_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE)
    @test !(:MBR_worst_donor_library_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE)
    @test !(:MBR_worst_donor_library_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE)
end

@testset "MBR observed iRT difference and single-donor flag use best/worst donors" begin
    mktempdir() do dir
        best_donor_path = joinpath(dir, "best_donor.arrow")
        worst_donor_path = joinpath(dir, "worst_donor.arrow")
        min_weight_donor_path = joinpath(dir, "min_weight_donor.arrow")
        single_donor_path = joinpath(dir, "single_donor.arrow")
        receiver_path = joinpath(dir, "receiver.arrow")
        single_receiver_path = joinpath(dir, "single_receiver.arrow")

        best_donor_df = _mbr_post_q_main_table(ms_file_idx = 1, scan_offset = 1000)
        best_donor_df[!, :irt_obs] = Float32[11.5, 19.0]
        worst_donor_df = _mbr_post_q_main_table(ms_file_idx = 2, scan_offset = 2000)
        worst_donor_df[!, :irt_obs] = Float32[17.0, 24.0]
        worst_donor_df[!, :weight] = Float32[8, 4]
        worst_donor_df[!, :n_scans] = Float32[5, 4]
        worst_donor_df[!, :trace_prob_prepass] = Float32[0.70, 0.60]
        worst_donor_df[!, :trace_prob_infold] = Float32[0.70, 0.60]
        min_weight_donor_df = _mbr_post_q_main_table(ms_file_idx = 6, scan_offset = 6000)
        min_weight_donor_df[!, :irt_obs] = Float32[13.0, 26.0]
        min_weight_donor_df[!, :weight] = Float32[2, 0.5]
        min_weight_donor_df[!, :n_scans] = Float32[3, 2]
        min_weight_donor_df[!, :trace_prob_prepass] = Float32[0.90, 0.50]
        min_weight_donor_df[!, :trace_prob_infold] = Float32[0.90, 0.50]
        single_donor_df = _mbr_post_q_main_table(ms_file_idx = 4, scan_offset = 4000)
        single_donor_df[!, :irt_obs] = Float32[30.0, 35.0]
        receiver_df = _mbr_post_q_main_table(ms_file_idx = 3, scan_offset = 3000)
        receiver_df[!, :irt_obs] = Float32[14.25, 22.0]
        receiver_df[!, :rt] = Float32[6.25, 22.0]
        receiver_df[!, :log_by_ratio_m0] = Float16[1.5, 0.0]
        single_receiver_df = _mbr_post_q_main_table(ms_file_idx = 5, scan_offset = 5000)
        single_receiver_df[!, :irt_obs] = Float32[31.25, 37.0]

        Arrow.write(best_donor_path, best_donor_df)
        Arrow.write(worst_donor_path, worst_donor_df)
        Arrow.write(min_weight_donor_path, min_weight_donor_df)
        Arrow.write(single_donor_path, single_donor_df)
        Arrow.write(receiver_path, receiver_df)
        Arrow.write(single_receiver_path, single_receiver_df)
        for (path, df) in (
            (best_donor_path, best_donor_df),
            (worst_donor_path, worst_donor_df),
            (min_weight_donor_path, min_weight_donor_df),
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
            [best_donor_path, worst_donor_path, min_weight_donor_path];
            passing_score_floor = 0.0f0,
        )
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [best_donor_path, worst_donor_path, min_weight_donor_path, receiver_path],
            precursors,
        )

        side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            Pioneer.build_mbr_fragment_annotation_keys(_mbr_post_q_fragment_lookup()),
            passing_score_floor = 0.0f0,
            run_similarity = Pioneer._MBRRunSimilarity(Dict(
                (UInt32(3), UInt32(1)) => 0.75f0,
                (UInt32(3), UInt32(2)) => 0.50f0,
                (UInt32(3), UInt32(6)) => 0.25f0,
            )),
        )
        side = DataFrame(Arrow.Table(side_path))

        receiver_residual = receiver_df.irt_pred[1] - receiver_df.irt_obs[1]
        best_donor_residual = best_donor_df.irt_pred[1] - best_donor_df.irt_obs[1]
        min_weight_donor_residual = min_weight_donor_df.irt_pred[1] - min_weight_donor_df.irt_obs[1]
        min_weight_false_residual = min_weight_donor_df.irt_pred[2] - min_weight_donor_df.irt_obs[2]

        @test side.MBR_best_pair_prob_true[1] == 0.95f0
        @test side.MBR_best_run_similarity_true[1] == 0.75f0
        @test side.MBR_best_run_similarity_false[1] == 0.75f0
        @test side.MBR_worst_run_similarity_true[1] == 0.25f0
        @test side.MBR_worst_run_similarity_false[1] == 0.25f0
        @test side.MBR_median_run_similarity_true[1] == 0.50f0
        @test side.MBR_median_run_similarity_false[1] == 0.50f0
        @test side.MBR_best_observed_irt_diff_true[1] == abs(14.25f0 - 11.5f0)
        @test side.MBR_best_irt_diff_true[1] == abs(receiver_residual - best_donor_residual)
        @test isapprox(
            side.MBR_best_log2_n_scans_ratio_true[1],
            log2((receiver_df.n_scans[1] + 1.0f0) / (best_donor_df.n_scans[1] + 1.0f0));
            atol = 1.0f-6,
        )
        @test side.MBR_worst_pair_prob_true[1] == 0.90f0
        @test side.MBR_worst_observed_irt_diff_true[1] == abs(14.25f0 - 13.0f0)
        @test side.MBR_worst_irt_diff_true[1] == abs(receiver_residual - min_weight_donor_residual)
        @test side.MBR_best_rt_diff_true[1] == abs(6.25f0 - best_donor_df.rt[1])
        @test side.MBR_best_log_by_diff_true[1] == 1.5f0
        @test isapprox(
            side.MBR_worst_log2_n_scans_ratio_true[1],
            log2((receiver_df.n_scans[1] + 1.0f0) / (min_weight_donor_df.n_scans[1] + 1.0f0));
            atol = 1.0f-6,
        )
        @test side.MBR_single_donor_true[1] == 0.0f0
        @test side.MBR_worst_pair_prob_false[1] == 0.50f0
        @test side.MBR_worst_observed_irt_diff_false[1] == abs(14.25f0 - 26.0f0)
        @test side.MBR_worst_irt_diff_false[1] == abs(receiver_residual - min_weight_false_residual)
        @test side.MBR_single_donor_false[1] == 0.0f0
        @test :MBR_worst_irt_diff_true in propertynames(side)
        @test :MBR_worst_irt_diff_false in propertynames(side)
        @test :MBR_worst_log2_n_scans_ratio_true in propertynames(side)
        @test :MBR_worst_log2_n_scans_ratio_false in propertynames(side)
        @test :MBR_best_log_by_diff_true in propertynames(side)
        @test :MBR_best_log_by_diff_false in propertynames(side)

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
        @test single_side.MBR_worst_pair_prob_true[1] == -1.0f0
        @test single_side.MBR_worst_run_similarity_true[1] == -1.0f0
        @test single_side.MBR_median_run_similarity_true[1] == 0.0f0
        @test single_side.MBR_worst_log2_weight_ratio_true[1] == -1.0f0
        @test single_side.MBR_worst_log2_explained_ratio_true[1] == -1.0f0
        @test single_side.MBR_worst_abs_n_scans_diff_true[1] == -1.0f0
        @test single_side.MBR_worst_log2_n_scans_ratio_true[1] == -1.0f0
        @test single_side.MBR_worst_observed_irt_diff_true[1] == -1.0f0
        @test single_side.MBR_worst_irt_diff_true[1] == -1.0f0
        @test single_side.MBR_single_donor_false[1] == 1.0f0
        @test single_side.MBR_worst_pair_prob_false[1] == -1.0f0
        @test single_side.MBR_worst_run_similarity_false[1] == -1.0f0
        @test single_side.MBR_median_run_similarity_false[1] == 0.0f0
        @test single_side.MBR_worst_log2_weight_ratio_false[1] == -1.0f0
        @test single_side.MBR_worst_log2_explained_ratio_false[1] == -1.0f0
        @test single_side.MBR_worst_abs_n_scans_diff_false[1] == -1.0f0
        @test single_side.MBR_worst_log2_n_scans_ratio_false[1] == -1.0f0
        @test single_side.MBR_worst_observed_irt_diff_false[1] == -1.0f0
        @test single_side.MBR_worst_irt_diff_false[1] == -1.0f0
    end
end
