# Unit tests for per-scan positional-isomer weight competition
# (phospho localization, Idea-1 Phase A). See
# src/Routines/SearchDIA/SearchMethods/MainSearch/features.jl:
#   add_isomer_competition_features!
using Pioneer, DataFrames, Test

@testset "add_isomer_competition_features!" begin
    # precursors 1,2 = isomers (group 1); 3 = an unrelated co-isolated peptide
    # (group 2); 4,5 = isomers (group 3). Vector indexed by precursor_idx.
    gids = UInt32[1, 1, 2, 3, 3]

    # PSMs, contiguous-by-scan (the function's precondition):
    #  scan 10: group1 weights 3,1 (sum 4 -> 0.75, 0.25); the unrelated
    #           precursor 3 at w=6 must NOT dilute group1 (it is a different
    #           group). This is the property that distinguishes this feature
    #           from weight_ratio_at_scan (which would give 3/6 = 0.5).
    #  scan 11: precursor 1 alone -> fraction 1.
    #  scan 12: group3 split 2,2 -> 0.5, 0.5.
    psms = DataFrame(
        scan_idx      = UInt32[10, 10, 10, 11, 12, 12],
        weight        = Float32[3, 1, 6, 5, 2, 2],
        precursor_idx = UInt32[1, 2, 3, 1, 4, 5],
    )

    Pioneer.add_isomer_competition_features!(psms, gids)

    @test psms.iso_weight_fraction_at_scan ≈ Float32[0.75, 0.25, 1.0, 1.0, 0.5, 0.5]
    @test psms.iso_group_size_at_scan == UInt16[2, 2, 1, 1, 2, 2]
    # sibling fractions within a scan+group sum to 1
    @test psms.iso_weight_fraction_at_scan[1] + psms.iso_weight_fraction_at_scan[2] ≈ 1.0f0
end

@testset "add_isomer_competition_features! empty" begin
    e = DataFrame(scan_idx = UInt32[], weight = Float32[], precursor_idx = UInt32[])
    Pioneer.add_isomer_competition_features!(e, UInt32[])
    @test nrow(e) == 0
    @test hasproperty(e, :iso_weight_fraction_at_scan)
    @test hasproperty(e, :iso_group_size_at_scan)
end
