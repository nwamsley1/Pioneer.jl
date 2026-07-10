using Test
using Pioneer
using Arrow
using DataFrames

struct FileAwareMockPrecursors <: Pioneer.LibraryPrecursors
    mz::Vector{Float32}
    irt::Vector{Float32}
    charge::Vector{UInt8}
    length::Vector{UInt8}
end

Pioneer.getMz(p::FileAwareMockPrecursors) = p.mz
Pioneer.getIrt(p::FileAwareMockPrecursors) = p.irt
Pioneer.getCharge(p::FileAwareMockPrecursors) = p.charge
Pioneer.getLength(p::FileAwareMockPrecursors) = p.length

function _test_mbr_donor(prob::Float32, file_idx::UInt32; precursor_idx::UInt32 = UInt32(1))
    return Pioneer._MBRDonorEntry(
        prob,
        precursor_idx,
        1.0f0,
        0.25f0,
        0.5f0,
        12.0f0,
        0.1f0,
        24.0f0,
        8.0f0,
        Pioneer.MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT,
        0.25f0,
        UInt8(0x03),
        UInt16(11),
        file_idx,
        false,
    )
end

function _test_partner_pools_for_receiver()
    target_by_pid = falses(11)
    target_by_pid[1] = true
    fold_by_pid = zeros(UInt8, 11)
    mz_bin_by_pid = ones(UInt8, 11)
    charge_by_pid = fill(UInt8(2), 11)
    length_by_pid = fill(UInt8(9), 11)
    irt_by_pid = zeros(Float32, 11)
    irt_by_pid[1] = 10.0f0

    decoy_pool = (pids = UInt32.(2:11), irts = Float32[10.01, 10.02, 10.03, 10.04, 10.05, 10.06, 10.07, 10.08, 10.09, 10.10])
    empty_pool = (pids = UInt32[], irts = Float32[])

    return Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        charge_by_pid,
        length_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int, Int, Int}, Pioneer._IrtPool}((0, 1, 2, 9) => decoy_pool),
        Dict{Tuple{Int, Int, Int}, Pioneer._IrtPool}((0, 2, 9) => decoy_pool),
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((2, 9) => decoy_pool),
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
        nothing,
        Float32[5, 5, 5, 5, 5],
        Float32[0.25, 0.25, 0.25, 0.25, 0.25],
        ntuple(_ -> zeros(Float32, 5), 8),
        UInt32[10, 10, 20, 30, 40],
        "unit-test",
    )

    entries = donor_dict[UInt32(1)]

    @test length(entries) == 2
    @test [entry.ms_file_idx for entry in entries] == UInt32[40, 10]
    @test [entry.trace_prob for entry in entries] == Float32[0.98, 0.95]
    @test Pioneer._donor_for_pid(donor_dict, UInt32(1), UInt32(40)).ms_file_idx == UInt32(10)
    @test Pioneer._donor_for_pid(donor_dict, UInt32(1), UInt32(10)).ms_file_idx == UInt32(40)
end

@testset "MBR counterfactual partner pool allows same-class non-self precursors" begin
    target_by_pid = trues(3)
    fold_by_pid = zeros(UInt8, 3)
    mz_bin_by_pid = ones(UInt8, 3)
    charge_by_pid = fill(UInt8(2), 3)
    length_by_pid = fill(UInt8(9), 3)
    irt_by_pid = Float32[0, 10.0, 10.01]

    target_pool = (pids = UInt32[1, 2], irts = Float32[10.0, 10.01])
    partner_pools = Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        charge_by_pid,
        length_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int, Int, Int}, Pioneer._IrtPool}((0, 1, 2, 9) => target_pool),
        Dict{Tuple{Int, Int, Int}, Pioneer._IrtPool}((0, 2, 9) => target_pool),
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((2, 9) => target_pool),
    )
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(1) => [_test_mbr_donor(0.99f0, UInt32(10))],
        UInt32(2) => [_test_mbr_donor(0.88f0, UInt32(20))],
    )

    donor = Pioneer._false_donor_for_pid(
        Dict{UInt32, Union{Nothing, Pioneer._MBRDonorEntry}}(),
        donor_dict,
        partner_pools,
        UInt32(1),
        UInt32(10),
    )

    @test donor !== nothing
    @test donor.trace_prob == 0.88f0
    @test donor.ms_file_idx == UInt32(20)
end

@testset "MBR counterfactual partner pool requires same charge and length" begin
    target_by_pid = trues(5)
    fold_by_pid = zeros(UInt8, 5)
    mz_bin_by_pid = ones(UInt8, 5)
    charge_by_pid = UInt8[2, 3, 2, 2, 2]
    length_by_pid = UInt8[9, 9, 10, 9, 9]
    irt_by_pid = Float32[10.0, 10.01, 10.02, 10.50, 11.0]

    same_charge_length_pool = (pids = UInt32[1, 4], irts = Float32[10.0, 10.50])
    partner_pools = Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        charge_by_pid,
        length_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int, Int, Int}, Pioneer._IrtPool}(
            (0, 1, 2, 9) => same_charge_length_pool,
            (0, 1, 3, 9) => (pids = UInt32[2], irts = Float32[10.01]),
            (0, 1, 2, 10) => (pids = UInt32[3], irts = Float32[10.02]),
        ),
        Dict{Tuple{Int, Int, Int}, Pioneer._IrtPool}((0, 2, 9) => same_charge_length_pool),
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((2, 9) => same_charge_length_pool),
    )
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(2) => [_test_mbr_donor(0.92f0, UInt32(20))],
        UInt32(3) => [_test_mbr_donor(0.91f0, UInt32(30))],
        UInt32(4) => [_test_mbr_donor(0.70f0, UInt32(40))],
    )

    donor = Pioneer._false_donor_for_pid(
        Dict{UInt32, Union{Nothing, Pioneer._MBRDonorEntry}}(),
        donor_dict,
        partner_pools,
        UInt32(1),
        UInt32(10),
    )

    @test donor !== nothing
    @test donor.trace_prob == 0.70f0
    @test donor.ms_file_idx == UInt32(40)
end

@testset "MBR false donor selection uses receiver row refined iRT" begin
    target_by_pid = trues(3)
    fold_by_pid = zeros(UInt8, 3)
    mz_bin_by_pid = ones(UInt8, 3)
    charge_by_pid = fill(UInt8(2), 3)
    length_by_pid = fill(UInt8(9), 3)
    irt_by_pid = Float32[100.0, 10.01, 99.9]

    same_charge_length_pool = (pids = UInt32[2, 3], irts = Float32[10.01, 99.9])
    partner_pools = Pioneer._CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        charge_by_pid,
        length_by_pid,
        irt_by_pid,
        Dict{Tuple{Int, Int, Int, Int}, Pioneer._IrtPool}((0, 1, 2, 9) => same_charge_length_pool),
        Dict{Tuple{Int, Int, Int}, Pioneer._IrtPool}((0, 2, 9) => same_charge_length_pool),
        Dict{Tuple{Int, Int}, Pioneer._IrtPool}((2, 9) => same_charge_length_pool),
    )
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(2) => [_test_mbr_donor(0.80f0, UInt32(20); precursor_idx = UInt32(2))],
        UInt32(3) => [_test_mbr_donor(0.90f0, UInt32(30); precursor_idx = UInt32(3))],
    )

    donors = Pioneer._resolve_false_donors_for_pid(
        donor_dict,
        partner_pools,
        UInt32(1),
        UInt32(10),
        nothing,
        10.0f0,
    )

    @test donors[1] !== nothing
    @test donors[1].precursor_idx == UInt32(2)
end

@testset "MBR counterfactual pools keep duplicate refined iRT observations but return unique CF precursors" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        donor_far_path = joinpath(dir, "donor_far.arrow")
        donor_near_path = joinpath(dir, "donor_near.arrow")

        base = DataFrame(
            precursor_idx = UInt32[],
            scan_idx = UInt32[],
            irt_pred = Float32[],
            ms_file_idx = UInt32[],
            target = Bool[],
            cv_fold = UInt8[],
        )
        receiver = copy(base)
        push!(receiver, (UInt32(1), UInt32(101), 10.0f0, UInt32(1), true, UInt8(0)))
        donor_far = copy(base)
        push!(donor_far, (UInt32(2), UInt32(201), 13.0f0, UInt32(2), true, UInt8(0)))
        push!(donor_far, (UInt32(3), UInt32(202), 10.10f0, UInt32(2), true, UInt8(0)))
        donor_near = copy(base)
        push!(donor_near, (UInt32(2), UInt32(301), 10.01f0, UInt32(3), true, UInt8(0)))

        Arrow.write(receiver_path, receiver)
        Arrow.write(donor_far_path, donor_far)
        Arrow.write(donor_near_path, donor_near)

        precursors = FileAwareMockPrecursors(
            Float32[500.0, 501.0, 502.0],
            Float32[10.0, 13.0, 10.10],
            UInt8[2, 2, 2],
            UInt8[9, 9, 9],
        )
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [receiver_path, donor_far_path, donor_near_path],
            precursors,
        )

        donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
            UInt32(2) => [_test_mbr_donor(0.80f0, UInt32(2); precursor_idx = UInt32(2))],
            UInt32(3) => [_test_mbr_donor(0.70f0, UInt32(2); precursor_idx = UInt32(3))],
        )

        donors = Pioneer._resolve_false_donors_for_pid(
            donor_dict,
            partner_pools,
            UInt32(1),
            UInt32(1),
            nothing,
            10.0f0,
        )

        @test donors[1] !== nothing
        @test donors[1].precursor_idx == UInt32(2)
        @test donors[2] !== nothing
        @test donors[2].precursor_idx == UInt32(3)
    end
end

@testset "MBR counterfactual donor selection uses the true donor file" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        true_donor_path = joinpath(dir, "true_donor_file.arrow")
        other_donor_path = joinpath(dir, "other_donor_file.arrow")

        base = DataFrame(
            precursor_idx = UInt32[],
            scan_idx = UInt32[],
            irt_pred = Float32[],
            ms_file_idx = UInt32[],
            target = Bool[],
            cv_fold = UInt8[],
        )
        receiver = copy(base)
        push!(receiver, (UInt32(1), UInt32(101), 10.0f0, UInt32(1), true, UInt8(0)))
        true_donor = copy(base)
        push!(true_donor, (UInt32(1), UInt32(201), 10.0f0, UInt32(2), true, UInt8(0)))
        push!(true_donor, (UInt32(3), UInt32(202), 10.05f0, UInt32(2), true, UInt8(0)))
        other_donor = copy(base)
        push!(other_donor, (UInt32(2), UInt32(301), 10.01f0, UInt32(3), true, UInt8(0)))

        Arrow.write(receiver_path, receiver)
        Arrow.write(true_donor_path, true_donor)
        Arrow.write(other_donor_path, other_donor)

        precursors = FileAwareMockPrecursors(
            Float32[500.0, 501.0, 502.0],
            Float32[10.0, 10.01, 10.05],
            UInt8[2, 2, 2],
            UInt8[9, 9, 9],
        )
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [receiver_path, true_donor_path, other_donor_path],
            precursors,
        )
        donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
            UInt32(2) => [_test_mbr_donor(0.95f0, UInt32(3); precursor_idx = UInt32(2))],
            UInt32(3) => [_test_mbr_donor(0.85f0, UInt32(2); precursor_idx = UInt32(3))],
        )

        donors = Pioneer._resolve_false_donors_for_pid(
            donor_dict,
            partner_pools,
            UInt32(1),
            UInt32(1),
            nothing,
            UInt32(2),
            10.0f0,
        )

        @test donors[1] !== nothing
        @test donors[1].precursor_idx == UInt32(3)
        @test donors[1].ms_file_idx == UInt32(2)
    end
end

@testset "MBR false donor selection uses nearest cross-file iRT partner" begin
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}()
    for pid in UInt32(2):UInt32(10)
        donor_dict[pid] = [_test_mbr_donor(0.20f0 + 0.01f0 * Float32(pid), UInt32(10))]
    end
    donor_dict[UInt32(11)] = [_test_mbr_donor(0.55f0, UInt32(20))]

    cache = Dict{UInt32, Union{Nothing, Pioneer._MBRDonorEntry}}()
    donor = Pioneer._false_donor_for_pid(
        cache,
        donor_dict,
        _test_partner_pools_for_receiver(),
        UInt32(1),
        UInt32(10),
    )

    @test donor !== nothing
    @test donor.trace_prob == 0.55f0
    @test donor.ms_file_idx == UInt32(20)
    @test haskey(cache, UInt32(1))
end
