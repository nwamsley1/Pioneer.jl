using Test
using Pioneer

function _test_mbr_donor(
    prob::Float32,
    file_idx::UInt32,
    pid::UInt32 = UInt32(1),
    n_scans::Float32 = 4.0f0,
)
    return Pioneer._MBRDonorEntry(
        prob,
        pid,
        1.0f0,
        0.25f0,
        0.5f0,
        12.0f0,
        0.1f0,
        n_scans,
        Pioneer.MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT,
        file_idx,
        false,
    )
end

function _test_partner_pools_for_receiver()
    target_by_pid = falses(11)
    target_by_pid[1] = true
    fold_by_pid = zeros(UInt8, 11)
    mz_bin_by_pid = ones(UInt8, 11)
    irt_by_pid = zeros(Float32, 11)
    irt_by_pid[1] = 10.0f0

    decoy_pool = (pids = UInt32.(2:11), irts = Float32[10.01, 10.02, 10.03, 10.04, 10.05, 10.06, 10.07, 10.08, 10.09, 10.10])
    empty_pool = (pids = UInt32[], irts = Float32[])

    return Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((0, 1) => decoy_pool),
        Dict{Int, Pioneer._IrtPool}(0 => decoy_pool, 1 => empty_pool),
        decoy_pool,
    )
end

@testset "MBR donor dictionary keeps top two distinct-file donors" begin
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}()

    Pioneer._accumulate_donor_entries!(
        donor_dict,
        UInt32[1, 1, 1, 1, 1],
        UInt32[1, 1, 1, 1, 1],
        UInt32[10, 11, 12, 13, 14],
        UInt32[10, 11, 12, 13, 14],
        Float32[0.95, 0.90, 0.70, 0.60, 0.98],
        Float32[1, 1, 1, 1, 1],
        Float16[0.5, 0.5, 0.5, 0.5, 0.5],
        Float32[12, 12, 12, 12, 12],
        Float32[11, 11, 11, 11, 11],
        nothing,
        Float32[2, 3, 4, 5, 6],
        ntuple(_ -> zeros(Float32, 5), 8),
        UInt32[10, 10, 20, 30, 40],
        "unit-test",
    )

    entries = donor_dict[UInt32(1)]

    @test length(entries) == 2
    @test [entry.ms_file_idx for entry in entries] == UInt32[40, 10]
    @test [entry.trace_prob for entry in entries] == Float32[0.98, 0.95]
    @test [entry.n_scans for entry in entries] == Float32[6, 2]
    @test Pioneer._donor_for_pid(donor_dict, UInt32(1), UInt32(40)).ms_file_idx == UInt32(10)
    @test Pioneer._donor_for_pid(donor_dict, UInt32(1), UInt32(10)).ms_file_idx == UInt32(40)
end

@testset "MBR false donor selection excludes the same precursor" begin
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(1) => [_test_mbr_donor(0.99f0, UInt32(20), UInt32(1))],
        UInt32(2) => [_test_mbr_donor(0.40f0, UInt32(20), UInt32(2))],
    )
    target_by_pid = trues(2)
    fold_by_pid = zeros(UInt8, 2)
    mz_bin_by_pid = ones(UInt8, 2)
    irt_by_pid = Float32[10.0, 10.01]
    pool = (pids = UInt32[1, 2], irts = Float32[10.0, 10.01])
    partner_pools = Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((0, 1) => pool),
        Dict{Int, Pioneer._IrtPool}(0 => pool, 1 => Pioneer._empty_pool()),
        pool,
    )
    hard_indexes = Pioneer.build_mbr_hard_counterfactual_indexes(donor_dict, partner_pools)

    donor = Pioneer._false_donor_for_pid(
        Dict{UInt32, Union{Nothing, Pioneer._MBRDonorEntry}}(),
        donor_dict,
        partner_pools,
        hard_indexes,
        UInt32(1),
        UInt32(10),
    )

    @test donor !== nothing
    @test donor.precursor_idx == UInt32(2)
end

@testset "MBR false donor selection skips same-file partners" begin
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}()
    for pid in UInt32(2):UInt32(10)
        donor_dict[pid] = [_test_mbr_donor(0.20f0 + 0.01f0 * Float32(pid), UInt32(10))]
    end
    donor_dict[UInt32(11)] = [_test_mbr_donor(0.55f0, UInt32(20))]

    cache = Dict{UInt32, Union{Nothing, Pioneer._MBRDonorEntry}}()
    partner_pools = _test_partner_pools_for_receiver()
    hard_indexes = Pioneer.build_mbr_hard_counterfactual_indexes(donor_dict, partner_pools)
    donor = Pioneer._false_donor_for_pid(
        cache,
        donor_dict,
        partner_pools,
        hard_indexes,
        UInt32(1),
        UInt32(10),
    )

    @test donor !== nothing
    @test donor.trace_prob == 0.55f0
    @test donor.ms_file_idx == UInt32(20)
    @test haskey(cache, UInt32(1))
end

@testset "MBR false donor selection prefers highest scoring local counterfactual" begin
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(2) => [_test_mbr_donor(0.30f0, UInt32(20))],
        UInt32(3) => [_test_mbr_donor(0.85f0, UInt32(20))],
        UInt32(4) => [_test_mbr_donor(0.99f0, UInt32(20))],
    )
    target_by_pid = falses(4)
    target_by_pid[1] = true
    fold_by_pid = zeros(UInt8, 4)
    mz_bin_by_pid = ones(UInt8, 4)
    irt_by_pid = zeros(Float32, 4)
    irt_by_pid[1] = 10.0f0
    decoy_pool = (pids = UInt32[2, 3, 4], irts = Float32[10.01, 10.20, 13.0])
    empty_pool = (pids = UInt32[], irts = Float32[])
    partner_pools = Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((0, 1) => decoy_pool),
        Dict{Int, Pioneer._IrtPool}(0 => decoy_pool, 1 => empty_pool),
        decoy_pool,
    )
    hard_indexes = Pioneer.build_mbr_hard_counterfactual_indexes(donor_dict, partner_pools)

    donor = Pioneer._false_donor_for_pid(
        Dict{UInt32, Union{Nothing, Pioneer._MBRDonorEntry}}(),
        donor_dict,
        partner_pools,
        hard_indexes,
        UInt32(1),
        UInt32(10),
    )

    @test donor !== nothing
    @test donor.trace_prob == 0.85f0
end

@testset "MBR false donor fold fallback uses nearest iRT" begin
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(2) => [_test_mbr_donor(0.30f0, UInt32(20))],
        UInt32(3) => [_test_mbr_donor(0.85f0, UInt32(20))],
    )
    target_by_pid = falses(3)
    target_by_pid[1] = true
    fold_by_pid = zeros(UInt8, 3)
    mz_bin_by_pid = ones(UInt8, 3)
    irt_by_pid = zeros(Float32, 3)
    irt_by_pid[1] = 10.0f0
    primary_pool = (pids = UInt32[], irts = Float32[])
    fold_pool = (pids = UInt32[2, 3], irts = Float32[10.01, 10.20])
    partner_pools = Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((0, 1) => primary_pool),
        Dict{Int, Pioneer._IrtPool}(0 => fold_pool, 1 => primary_pool),
        fold_pool,
    )
    hard_indexes = Pioneer.build_mbr_hard_counterfactual_indexes(donor_dict, partner_pools)

    donor = Pioneer._false_donor_for_pid(
        Dict{UInt32, Union{Nothing, Pioneer._MBRDonorEntry}}(),
        donor_dict,
        partner_pools,
        hard_indexes,
        UInt32(1),
        UInt32(10),
    )

    @test donor !== nothing
    @test donor.trace_prob == 0.30f0
end
