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
