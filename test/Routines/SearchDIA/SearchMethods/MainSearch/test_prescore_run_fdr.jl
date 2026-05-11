using Test
using DataFrames

@testset "MainSearch prescore per-run FDR filtering" begin
    fold0 = DataFrame(
        precursor_idx = UInt32[1, 2],
        lgbm_prob = Float32[0.99, 0.98],
        target = Bool[false, true],
        cv_fold = UInt8[0, 0],
    )
    fold1 = DataFrame(
        precursor_idx = UInt32[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        lgbm_prob = Float32[0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.88, 0.10],
        target = Bool[true, true, true, true, true, true, true, true, true, true, false],
        cv_fold = UInt8[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    )

    result = Pioneer._filter_prescore_run_qvalues([fold0, fold1], 0.10f0)

    @test result.n_before == 13
    @test result.n_after == 12
    @test result.n_targets_pass == 11
    @test result.n_decoys_pass == 1

    @test result.filtered_tables[1].precursor_idx == UInt32[1, 2]
    @test result.filtered_tables[2].precursor_idx == UInt32[3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
end

@testset "MainSearch prescore rescues cross-run candidates" begin
    fold0 = DataFrame(
        precursor_idx = UInt32[1, 2],
        lgbm_prob = Float32[0.99, 0.98],
        target = Bool[false, true],
        cv_fold = UInt8[0, 0],
    )
    fold1 = DataFrame(
        precursor_idx = UInt32[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        lgbm_prob = Float32[0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.88, 0.10, 0.09],
        target = Bool[true, true, true, true, true, true, true, true, true, true, false, true],
        cv_fold = UInt8[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    )

    strict_runs = Dict{UInt32, Set{UInt32}}(
        UInt32(13) => Set(UInt32[2]),
        UInt32(14) => Set(UInt32[3]),
    )

    result = Pioneer._filter_prescore_run_qvalues(
        [fold0, fold1],
        0.10f0;
        ms_file_idx = UInt32(1),
        strict_runs_by_precursor = strict_runs,
    )

    @test result.n_before == 14
    @test result.n_after == 14
    @test result.n_strict_pass == 12
    @test result.n_rescue_pass == 2
    @test result.n_rescue_targets == 1
    @test result.n_rescue_decoys == 1

    combined = vcat(result.filtered_tables...)
    @test combined.precursor_idx == UInt32[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    @test combined.prescore_strict_pass == Bool[trues(12); false; false]
    @test combined.prescore_rescue_candidate == Bool[falses(12); true; true]
    @test combined.prescore_donor_run_count == UInt16[zeros(UInt16, 12); 1; 1]
    @test all(combined.prescore_qval .>= 0.0f0)
end

@testset "MBR donor features ignore prescore-rescued rows as donors" begin
    psms = DataFrame(
        pair_id = UInt32[10, 10, 10],
        trace_prob_prepass = Float32[0.99, 0.80, 0.20],
        weight = Float32[16, 8, 4],
        log2_intensity_explained = Float32[6, 4, 1],
        ms_file_idx = UInt32[1, 2, 3],
        decoy = Bool[false, false, false],
        prescore_strict_pass = Bool[false, true, true],
    )

    Pioneer.compute_mbr_features!(psms; score_col = :trace_prob_prepass)

    @test psms.MBR_max_pair_prob[3] == 0.80f0
    @test psms.MBR_log2_weight_ratio[3] == -1.0f0
    @test !psms.MBR_is_missing[3]
end
