using Test
using DataFrames
using Pioneer

using Pioneer: MainSearchIrtCorrectionModel, MainSearchIrtRefinement, _compute_phase2_columns!
using Pioneer: _passing_precursor_targets, refine_mainsearch_irt_predictions!
using Pioneer: reapply_psm_classifier_and_select_best!, train_lgbm_and_select_best
using Pioneer: train_psm_classifier_with_fallback

function _mock_scan_cycle_index(pairs::Vector{Tuple{Int32, Int32}})
    n = maximum(first, pairs)
    scan_to_window_key = fill((Int32(0), Int32(0)), n)
    scan_to_window_pos = zeros(Int32, n)
    key = (Int32(500000), Int32(24000))
    for (scan, pos) in pairs
        scan_to_window_key[Int(scan)] = key
        scan_to_window_pos[Int(scan)] = pos
    end
    return (scan_to_window_key = scan_to_window_key, scan_to_window_pos = scan_to_window_pos)
end

struct MainSearchMockPrecursors <: Pioneer.LibraryPrecursors
    sequence::Vector{String}
    structural_mods::Vector{Union{Missing, String}}
    mz::Vector{Float32}
    irt::Vector{Float32}
end

@testset "MainSearch classifier fallback" begin
    @testset "shared LightGBM configs enable unbalanced labels" begin
        @test Pioneer.SHARED_LGBM_HP.is_unbalance === true
        @test Pioneer.SCORING_LGBM_HP.is_unbalance === true
    end

    @testset "low-data model selection accepts non-LightGBM candidates" begin
        psms = DataFrame(
            target = Bool[
                true, false, false, true,
                true, false, false, true,
                true, false, false, true,
            ],
            cv_fold = UInt8[
                1, 0, 1, 0,
                1, 0, 1, 0,
                1, 0, 1, 0,
            ],
            discriminant = Float32[
                0.9, 0.1, 0.2, 0.8,
                0.85, 0.15, 0.25, 0.75,
                0.88, 0.12, 0.22, 0.78,
            ],
        )

        tiny_lgbm_hp = (
            num_iterations = 3,
            learning_rate = 0.2,
            max_depth = 2,
            num_leaves = 3,
            min_data_in_leaf = 1,
            feature_fraction = 1.0,
            bagging_fraction = 1.0,
            bagging_freq = 0,
            is_unbalance = false,
            max_bin = 16,
            lambda_l1 = 0.0,
            lambda_l2 = 0.0,
        )

        scores, infold_scores, _last_classifier, info = train_psm_classifier_with_fallback(
            psms;
            features = [:discriminant],
            lgbm_hp = tiny_lgbm_hp,
        )

        @test length(scores) == nrow(psms)
        @test infold_scores === nothing
        @test info.low_data
        @test haskey(info.candidate_oof, "probit")
        @test all(isfinite, scores)
    end

end

@testset "MainSearch wide-window core bounds" begin
    @testset "core follows cycle-contiguous passing scans around best PSM" begin
        scan_cycle_index = _mock_scan_cycle_index([
            (Int32(100), Int32(10)),
            (Int32(105), Int32(11)),
            (Int32(110), Int32(12)),
            (Int32(120), Int32(14)),
        ])

        bounds = Pioneer._mainsearch_wide_core_bounds(
            UInt32[100, 105, 110, 120],
            Float32[0.70, 0.95, 0.80, 0.90],
            Bool[true, true, true, true],
            scan_cycle_index,
        )

        @test bounds.scan_min == UInt32(100)
        @test bounds.scan_max == UInt32(110)
        @test bounds.n_scans == UInt16(3)
    end

    @testset "non-passing scan breaks core contiguity" begin
        scan_cycle_index = _mock_scan_cycle_index([
            (Int32(100), Int32(10)),
            (Int32(105), Int32(11)),
            (Int32(110), Int32(12)),
        ])

        bounds = Pioneer._mainsearch_wide_core_bounds(
            UInt32[100, 105, 110],
            Float32[0.70, 0.95, 0.80],
            Bool[true, true, false],
            scan_cycle_index,
        )

        @test bounds.scan_min == UInt32(100)
        @test bounds.scan_max == UInt32(105)
        @test bounds.n_scans == UInt16(2)
    end

    @testset "best PSM row carries wide core metadata" begin
        scan_cycle_index = _mock_scan_cycle_index([
            (Int32(100), Int32(10)),
            (Int32(105), Int32(11)),
            (Int32(110), Int32(12)),
            (Int32(120), Int32(14)),
        ])
        psms = DataFrame(
            precursor_idx = UInt32[7, 7, 7, 7],
            scan_idx = UInt32[100, 105, 110, 120],
            lgbm_score = Float32[0.70, 0.95, 0.80, 0.90],
            weight = Float32[10, 20, 15, 12],
            irt_obs = Float32[1, 2, 3, 4],
            rt = Float32[10, 20, 30, 40],
        )

        best = Pioneer.select_best_per_precursor!(
            psms,
            :lgbm_score;
            pass_mask = Bool[true, true, true, true],
            scan_cycle_index = scan_cycle_index,
        )

        @test nrow(best) == 1
        @test best.scan_idx[1] == UInt32(105)
        @test best.wide_core_scan_min[1] == UInt32(100)
        @test best.wide_core_scan_max[1] == UInt32(110)
        @test best.wide_core_n_scans[1] == UInt16(3)
    end

    @testset "best PSM row counts passing scans from other isotope traces" begin
        psms = DataFrame(
            precursor_idx = UInt32[7, 7, 7, 7, 7],
            scan_idx = UInt32[100, 105, 110, 115, 120],
            cycle_idx = UInt32[1, 2, 1, 2, 1],
            lgbm_score = Float32[0.70, 0.95, 0.80, 0.75, 0.60],
            weight = Float32[10, 20, 15, 12, 8],
            irt_obs = Float32[1, 2, 3, 4, 5],
            rt = Float32[10, 20, 30, 40, 50],
            frag1_int = Float16[10, 20, 15, 12, 8],
            frag2_int = Float16[0, 0, 0, 0, 0],
            frag3_int = Float16[0, 0, 0, 0, 0],
            frag4_int = Float16[0, 0, 0, 0, 0],
            frag5_int = Float16[0, 0, 0, 0, 0],
            frag6_int = Float16[0, 0, 0, 0, 0],
            frag7_int = Float16[0, 0, 0, 0, 0],
            frag8_int = Float16[0, 0, 0, 0, 0],
            isotopes_captured = [
                (Int8(0), Int8(1)),
                (Int8(0), Int8(1)),
                (Int8(1), Int8(2)),
                (Int8(1), Int8(2)),
                (Int8(2), Int8(3)),
            ],
            precursor_fraction_transmitted = Float32[0.7, 0.8, 0.4, 0.4, 0.2],
        )

        best = Pioneer.select_best_per_precursor!(
            psms,
            :lgbm_score;
            pass_mask = Bool[true, true, true, false, true],
        )

        @test nrow(best) == 1
        @test best.scan_idx[1] == UInt32(105)
        @test best.isotopes_captured[1] == (Int8(0), Int8(1))
        @test best.precursor_fraction_transmitted[1] == Float32(0.8)
        @test best.n_scans_other_traces[1] == UInt32(2)
    end

    @testset "best PSM row carries selected-trace fragment support union and intersection" begin
        psms = DataFrame(
            precursor_idx = fill(UInt32(7), 6),
            scan_idx = UInt32[100, 105, 110, 115, 120, 125],
            cycle_idx = UInt32[1, 2, 3, 4, 2, 3],
            lgbm_score = Float32[0.70, 0.95, 0.80, 0.75, 0.60, 0.55],
            weight = Float32[10, 20, 15, 12, 5, 6],
            irt_obs = Float32[1, 2, 3, 4, 2, 3],
            rt = Float32[10, 20, 30, 40, 20, 30],
            frag1_int = Float16[10, 20, 0, 5, 99, 99],
            frag2_int = Float16[0, 10, 15, 0, 99, 99],
            frag3_int = Float16[5, 7, 9, 0, 99, 99],
            frag4_int = Float16[0, 0, 2, 0, 99, 99],
            frag5_int = Float16[0, 0, 0, 0, 99, 99],
            frag6_int = Float16[1, 1, 1, 0, 99, 99],
            frag7_int = Float16[7, 7, 7, 0, 99, 99],
            frag8_int = Float16[8, 8, 0, 0, 99, 99],
            isotopes_captured = [
                (Int8(0), Int8(1)),
                (Int8(0), Int8(1)),
                (Int8(0), Int8(1)),
                (Int8(0), Int8(1)),
                (Int8(1), Int8(2)),
                (Int8(1), Int8(2)),
            ],
            precursor_fraction_transmitted = Float32[0.8, 0.8, 0.8, 0.8, 0.4, 0.4],
        )

        bitvec_rank_table = fill(UInt16(256), 256)
        bitvec_rank_table[Int(0xef) + 1] = UInt16(11)
        bitvec_rank_table[Int(0x64) + 1] = UInt16(22)

        best = Pioneer.select_best_per_precursor!(
            psms,
            :lgbm_score;
            pass_mask = Bool[true, true, true, false, true, true],
            bitvec_rank_table = bitvec_rank_table,
            pred_fragment_intensity_provider! = (buf, pid) -> begin
                @test pid == UInt32(7)
                buf[1] = 30f0
                buf[2] = 25f0
                buf[3] = 21f0
                buf[4] = 2f0
                buf[5] = 0f0
                buf[6] = 3f0
                buf[7] = 21f0
                buf[8] = 16f0
                nothing
            end,
        )

        @test nrow(best) == 1
        @test best.scan_idx[1] == UInt32(105)
        @test best.n_frags_detected_union[1] == UInt8(7)
        @test best.n_frags_detected_intersection[1] == UInt8(3)
        @test best.n_frags_detected_union_bitvec_rank[1] == UInt16(11)
        @test best.n_frags_detected_intersection_bitvec_rank[1] == UInt16(22)
        @test best.frag_observed_sum_hellinger[1] ≈ 0.0f0 atol=1f-6
    end

    @testset "best PSM row carries aligned isotope trace agreement features" begin
        psms = DataFrame(
            precursor_idx = fill(UInt32(7), 9),
            scan_idx = UInt32[100, 105, 110, 115, 120, 125, 130, 135, 140],
            cycle_idx = UInt32[1, 2, 3, 1, 2, 3, 2, 3, 4],
            lgbm_score = Float32[0.70, 0.95, 0.80, 0.60, 0.65, 0.62, 0.55, 0.58, 0.57],
            weight = Float32[1, 4, 1, 1, 4, 1, 1, 4, 1],
            irt_obs = Float32[1, 2, 3, 1, 2, 3, 2, 3, 4],
            rt = Float32[10, 20, 30, 10, 20, 30, 20, 30, 40],
            frag1_int = Float16[2, 8, 2, 2, 8, 2, 1, 4, 1],
            frag2_int = Float16[1, 4, 1, 1, 4, 1, 1, 4, 1],
            frag3_int = Float16[0, 0, 0, 0, 0, 0, 0, 0, 0],
            frag4_int = Float16[0, 0, 0, 0, 0, 0, 0, 0, 0],
            frag5_int = Float16[0, 0, 0, 0, 0, 0, 0, 0, 0],
            frag6_int = Float16[0, 0, 0, 0, 0, 0, 0, 0, 0],
            frag7_int = Float16[0, 0, 0, 0, 0, 0, 0, 0, 0],
            frag8_int = Float16[0, 0, 0, 0, 0, 0, 0, 0, 0],
            isotopes_captured = [
                (Int8(0), Int8(1)), (Int8(0), Int8(1)), (Int8(0), Int8(1)),
                (Int8(1), Int8(2)), (Int8(1), Int8(2)), (Int8(1), Int8(2)),
                (Int8(2), Int8(3)), (Int8(2), Int8(3)), (Int8(2), Int8(3)),
            ],
            precursor_fraction_transmitted = Float32[0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2],
        )

        best = Pioneer.select_best_per_precursor!(
            psms,
            :lgbm_score;
            pass_mask = trues(nrow(psms)),
        )

        @test nrow(best) == 1
        @test best.scan_idx[1] == UInt32(105)
        @test best.n_scans_other_traces[1] == UInt32(6)
        @test best.trace_other_weight_corr[1] ≈ 1.0f0 atol=1f-6
        @test best.trace_other_frag_sum_corr[1] ≈ 1.0f0 atol=1f-6
        @test best.trace_other_apex_delta_irt[1] ≈ 0.0f0
    end

    @testset "missing other isotope trace uses agreement sentinels" begin
        psms = DataFrame(
            precursor_idx = UInt32[7, 7, 7],
            scan_idx = UInt32[100, 105, 110],
            cycle_idx = UInt32[1, 2, 3],
            lgbm_score = Float32[0.70, 0.95, 0.80],
            weight = Float32[1, 4, 1],
            irt_obs = Float32[1, 2, 3],
            rt = Float32[10, 20, 30],
            frag1_int = Float16[2, 8, 2],
            frag2_int = Float16[1, 4, 1],
            frag3_int = Float16[0, 0, 0],
            frag4_int = Float16[0, 0, 0],
            frag5_int = Float16[0, 0, 0],
            frag6_int = Float16[0, 0, 0],
            frag7_int = Float16[0, 0, 0],
            frag8_int = Float16[0, 0, 0],
            isotopes_captured = [
                (Int8(0), Int8(1)),
                (Int8(0), Int8(1)),
                (Int8(0), Int8(1)),
            ],
            precursor_fraction_transmitted = Float32[0.8, 0.8, 0.8],
        )

        best = Pioneer.select_best_per_precursor!(
            psms,
            :lgbm_score;
            pass_mask = trues(nrow(psms)),
        )

        @test best.n_scans_other_traces[1] == UInt32(0)
        @test best.trace_other_weight_corr[1] == -1.0f0
        @test best.trace_other_frag_sum_corr[1] == -1.0f0
        @test best.trace_other_apex_delta_irt[1] == 100.0f0
    end
end

Pioneer.getSequence(p::MainSearchMockPrecursors) = p.sequence
Pioneer.getStructuralMods(p::MainSearchMockPrecursors) = p.structural_mods
Pioneer.getMz(p::MainSearchMockPrecursors) = p.mz
Pioneer.getIrt(p::MainSearchMockPrecursors) = p.irt

@testset "MainSearch iRT refinement" begin
    @testset "passing targets use prescore q-value criteria" begin
        precursor_ids, irt_pred_inputs, irt_corrections = _passing_precursor_targets(
            UInt32[11, 22, 33],
            Bool[true, false, true],
            Float32[0.99, 0.98, 0.97],
            Float32[10.0, 20.0, 30.0],
            Float32[12.0, 18.0, 40.0],
            0.01f0,
        )

        @test precursor_ids == UInt32[11]
        @test irt_pred_inputs == Float32[10.0]
        @test irt_corrections == Float32[2.0]
    end

    @testset "refinement updates predictions and errors before classifier reapply" begin
        precursors = MainSearchMockPrecursors(
            ["AAAA", "CCCC", "AAAA", "CCCC", "DDDD", "DDDD"],
            fill(missing, 6),
            fill(500.0f0, 6),
            fill(10.0f0, 6),
        )
        psms = DataFrame(
            target = Bool[true, true, false, true, true, false],
            precursor_idx = UInt32[1, 2, 5, 3, 4, 6],
            irt_pred = fill(10.0f0, 6),
            irt_obs = Float32[4, 8, 6, 4, 8, 6],
            irt_error = fill(0.0f0, 6),
        )
        psms[!, :irt_error] = abs.(psms.irt_obs .- psms.irt_pred)
        scores = Float32[0.99, 0.98, 0.01, 0.99, 0.98, 0.01]
        best_psms = copy(psms)

        result = refine_mainsearch_irt_predictions!(
            psms,
            best_psms,
            scores,
            MainSearchIrtRefinement(precursors; q_value_threshold = 0.01f0, min_precursors = 2),
        )

        @test result.refined
        @test Set(result.training_target_precursors) == Set(UInt32[1, 2, 3, 4])
        @test psms.irt_pred[1] ≈ 4.0f0 atol = 1f-4
        @test psms.irt_pred[2] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_pred[4] ≈ 4.0f0 atol = 1f-4
        @test psms.irt_pred[5] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_error[1] ≈ 0.0f0 atol = 1f-4
        @test psms.irt_error[2] ≈ 0.0f0 atol = 1f-4
    end

    @testset "refinement trains out-of-fold models when cv_fold is present" begin
        precursors = MainSearchMockPrecursors(
            fill("AAAA", 6),
            fill(missing, 6),
            fill(500.0f0, 6),
            fill(10.0f0, 6),
        )
        psms = DataFrame(
            target = Bool[true, true, true, true, false, false],
            precursor_idx = UInt32[1, 2, 3, 4, 5, 6],
            cv_fold = UInt8[0, 0, 1, 1, 0, 1],
            irt_pred = fill(10.0f0, 6),
            irt_obs = Float32[8, 8, 13, 13, 8, 13],
            irt_error = fill(0.0f0, 6),
        )
        best_psms = copy(psms)
        scores = Float32[0.99, 0.98, 0.99, 0.98, 0.01, 0.01]

        result = refine_mainsearch_irt_predictions!(
            psms,
            best_psms,
            scores,
            MainSearchIrtRefinement(precursors; q_value_threshold = 0.01f0, min_precursors = 2),
        )

        @test result.refined
        @test result.model isa Dict{UInt8, MainSearchIrtCorrectionModel}
        @test sort(collect(keys(result.model))) == UInt8[0, 1]
        @test Set(result.training_target_precursors) == Set(UInt32[1, 2, 3, 4])
        @test psms.irt_pred[1] ≈ 13.0f0 atol = 1f-4
        @test psms.irt_pred[2] ≈ 13.0f0 atol = 1f-4
        @test psms.irt_pred[3] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_pred[4] ≈ 8.0f0 atol = 1f-4
        @test best_psms.irt_pred[1] ≈ 13.0f0 atol = 1f-4
        @test best_psms.irt_pred[3] ≈ 8.0f0 atol = 1f-4
    end

    @testset "phase two iRT difference uses refined prediction column" begin
        irt_diff_col = Vector{Float32}(undef, 2)
        prec_mz_col = Vector{Float32}(undef, 2)
        pair_id_col = Vector{UInt32}(undef, 2)
        entrap_col = Vector{UInt8}(undef, 2)

        _compute_phase2_columns!(
            UInt32[1, 2],
            Float32[12.0, 12.0],
            Float32[11.0, 15.0],
            Float32[100.0, 100.0],
            Float32[401.0, 402.0],
            UInt32[101, 102],
            UInt8[0, 1],
            irt_diff_col,
            prec_mz_col,
            pair_id_col,
            entrap_col,
        )

        @test irt_diff_col == Float32[1.0, 3.0]
        @test prec_mz_col == Float32[401.0, 402.0]
        @test pair_id_col == UInt32[101, 102]
        @test entrap_col == UInt8[0, 1]
    end

end
