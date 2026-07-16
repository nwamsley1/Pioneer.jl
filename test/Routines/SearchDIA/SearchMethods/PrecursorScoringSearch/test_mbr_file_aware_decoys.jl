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
        nothing,
        UInt32[10, 10, 20, 30, 40],
        nothing,
        nothing,
        "unit-test",
    )

    entries = donor_dict[UInt32(1)]

    @test length(entries) == 2
    @test [entry.ms_file_idx for entry in entries] == UInt32[40, 10]
    @test [entry.trace_prob for entry in entries] == Float32[0.98, 0.95]
    @test Pioneer._donor_for_pid(donor_dict, UInt32(1), UInt32(40)).ms_file_idx == UInt32(10)
    @test Pioneer._donor_for_pid(donor_dict, UInt32(1), UInt32(10)).ms_file_idx == UInt32(40)
end

@testset "MBR donor dictionary can retain all donor files for similarity selection" begin
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}()

    Pioneer._accumulate_donor_entries!(
        donor_dict,
        UInt32[1, 1, 1, 1],
        UInt32[1, 1, 1, 1],
        UInt32[10, 11, 12, 13],
        UInt32[10, 11, 12, 13],
        Float32[0.95, 0.90, 0.70, 0.60],
        Float32[1, 1, 1, 1],
        Float16[0.5, 0.5, 0.5, 0.5],
        Float32[12, 12, 12, 12],
        Float32[11, 11, 11, 11],
        nothing,
        nothing,
        Float32[5, 5, 5, 5],
        Float32[0.25, 0.25, 0.25, 0.25],
        ntuple(_ -> zeros(Float32, 4), 8),
        nothing,
        UInt32[10, 20, 30, 40],
        nothing,
        nothing,
        "unit-test",
        max_donor_files_per_precursor = typemax(Int),
    )

    @test [entry.ms_file_idx for entry in donor_dict[UInt32(1)]] == UInt32[10, 20, 30, 40]
end

@testset "MBR true donor selection uses IDF-weighted receiver containment" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        donor_low_path = joinpath(dir, "donor_low.arrow")
        donor_high_path = joinpath(dir, "donor_high.arrow")
        donor_equal_path = joinpath(dir, "donor_equal.arrow")

        function write_similarity_table(path, file_idx, pids, qvals)
            Arrow.write(path, DataFrame(
                precursor_idx = UInt32.(pids),
                ms_file_idx = fill(UInt32(file_idx), length(pids)),
                target = trues(length(pids)),
                qval = Float32.(qvals),
            ))
        end

        write_similarity_table(
            receiver_path,
            1,
            [1, 1, 2, 3, 10],
            [0.001, 0.001, 0.001, 0.001, 0.05],
        )
        write_similarity_table(donor_low_path, 2, [1, 10], [0.001, 0.001])
        write_similarity_table(donor_high_path, 3, [1, 2, 10], [0.001, 0.001, 0.001])
        write_similarity_table(donor_equal_path, 4, [1, 10], [0.001, 0.001])

        similarity = Pioneer.build_mbr_run_similarity(
            [receiver_path, donor_low_path, donor_high_path, donor_equal_path];
            q_value_threshold = 0.01f0,
        )

        idf_pid2 = log(5 / 3)
        idf_pid3 = log(5 / 2)
        expected_high = Float32(idf_pid2 / (idf_pid2 + idf_pid3))
        @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(2)) == 0.0f0
        @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(3)) ≈ expected_high
        @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(4)) == 0.0f0

        donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
            UInt32(10) => [
                _test_mbr_donor(0.99f0, UInt32(2); precursor_idx = UInt32(10)),
                _test_mbr_donor(0.70f0, UInt32(3); precursor_idx = UInt32(10)),
            ],
            UInt32(11) => [
                _test_mbr_donor(0.90f0, UInt32(2); precursor_idx = UInt32(11)),
                _test_mbr_donor(0.80f0, UInt32(4); precursor_idx = UInt32(11)),
            ],
        )

        selected = Pioneer._donor_for_pid(donor_dict, UInt32(10), UInt32(1), similarity)
        tied = Pioneer._donor_for_pid(donor_dict, UInt32(11), UInt32(1), similarity)

        @test selected !== nothing
        @test selected.ms_file_idx == UInt32(3)
        @test tied !== nothing
        @test tied.ms_file_idx == UInt32(2)
    end
end

@testset "MBR donor selection ignores cluster-conditioned similarity" begin
    active_files = BitSet((1, 2, 3))
    similarity = Pioneer._MBRRunSimilarity(
        Dict(
            (UInt32(1), UInt32(2)) => 0.8f0,
            (UInt32(1), UInt32(3)) => 0.2f0,
        ),
        Dict{Tuple{UInt32, UInt32}, Float32}(),
        Dict{Tuple{UInt32, UInt32}, Float32}(),
        Dict{UInt32, Float32}(),
        Dict{UInt32, Float32}(),
        active_files,
        Dict{UInt32, BitSet}(),
        UInt32[1],
        Float32[1.0],
        Float32[1.0],
        UInt16[2],
        Pioneer._MBRClusterRunSimilarity[
            Pioneer._MBRClusterRunSimilarity(
                Dict(
                    (UInt32(1), UInt32(2)) => 2.0f0,
                    (UInt32(1), UInt32(3)) => 9.0f0,
                ),
                Dict{Tuple{UInt32, UInt32}, Float32}(),
                Dict(UInt32(1) => 10.0f0),
                Dict{UInt32, Float32}(),
                UInt32(10),
            ),
        ],
    )
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(1) => [
            _test_mbr_donor(0.70f0, UInt32(2)),
            _test_mbr_donor(0.99f0, UInt32(3)),
        ],
    )

    @test Pioneer._mbr_cluster_shared_support(
        similarity,
        UInt32(1),
        UInt32(1),
        UInt32(3),
    ) > Pioneer._mbr_cluster_shared_support(
        similarity,
        UInt32(1),
        UInt32(1),
        UInt32(2),
    )
    selected = Pioneer._donor_for_pid(
        donor_dict,
        UInt32(1),
        UInt32(1),
        similarity,
    )
    @test selected !== nothing
    @test selected.ms_file_idx == UInt32(2)
end

@testset "MBR run similarity includes passing decoy postings" begin
    mktempdir() do dir
        paths = [joinpath(dir, "run$(idx).arrow") for idx in 1:3]
        Arrow.write(paths[1], DataFrame(
            precursor_idx = UInt32[1, 2],
            ms_file_idx = UInt32[1, 1],
            target = Bool[true, false],
            qval = Float32[0.001, 0.001],
        ))
        Arrow.write(paths[2], DataFrame(
            precursor_idx = UInt32[1],
            ms_file_idx = UInt32[2],
            target = Bool[true],
            qval = Float32[0.001],
        ))
        Arrow.write(paths[3], DataFrame(
            precursor_idx = UInt32[3],
            ms_file_idx = UInt32[3],
            target = Bool[true],
            qval = Float32[0.001],
        ))

        similarity = Pioneer.build_mbr_run_similarity(paths; q_value_threshold = 0.01f0)
        shared_idf = log(4 / 3)
        decoy_idf = log(4 / 2)
        expected = Float32(shared_idf / (shared_idf + decoy_idf))
        @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(2)) ≈ expected

        means = Pioneer._mbr_mean_run_similarity_by_file(similarity, 3)
        @test means[UInt32(1)] ≈ expected / 2.0f0
        @test means[UInt32(2)] ≈ 0.5f0
        @test means[UInt32(3)] == 0.0f0
    end
end

@testset "MBR run similarity builds intensity clusters from search outputs" begin
    mktempdir() do dir
        paths = [joinpath(dir, "run$(idx).arrow") for idx in 1:2]
        for (file_idx, path) in enumerate(paths)
            Arrow.write(path, DataFrame(
                precursor_idx = UInt32[1, 2, 3, 4],
                ms_file_idx = fill(UInt32(file_idx), 4),
                target = trues(4),
                qval = fill(0.001f0, 4),
                weight = Float32[100, 200, 300, 400],
            ))
        end

        similarity = Pioneer.build_mbr_run_similarity(paths)

        @test !isempty(similarity.cluster_run_similarity)
        @test length(similarity.cluster_by_precursor) == 4
        @test Pioneer._mbr_precursor_present(similarity, UInt32(1), UInt32(1))
        @test Pioneer._mbr_cluster_agreement(
            similarity,
            UInt32(1),
            UInt32(1),
            UInt32(2),
        ) == 1.0f0
        @test Pioneer._mbr_cluster_active_fraction(
            similarity,
            UInt32(1),
            UInt32(1),
            UInt32(2),
        ) == 1.0f0
        @test Pioneer._mbr_cluster_peer_count(similarity, UInt32(1)) == 3.0f0
    end
end

@testset "MBR intensity clustering adapts leaf count to precursor population" begin
    file_ids = UInt32[1, 2]
    identical_postings = Dict(
        precursor_idx => Pioneer._MBRIntensityPosting[
            Pioneer._MBRIntensityPosting(UInt32(1), log2(100.0f0)),
            Pioneer._MBRIntensityPosting(UInt32(2), log2(100.0f0)),
        ]
        for precursor_idx in UInt32(1):UInt32(1024)
    )
    identical_fit = Pioneer._fit_mbr_intensity_clusters(
        identical_postings,
        file_ids,
    )
    identical_sizes = Dict{UInt32, Int}()
    for cluster_idx in identical_fit.cluster_by_precursor
        identical_sizes[cluster_idx] = get(identical_sizes, cluster_idx, 0) + 1
    end

    @test eltype(identical_fit.cluster_by_precursor) == UInt32
    @test identical_fit.n_clusters == 1
    @test only(values(identical_sizes)) == 1024
    @test all(iszero, identical_fit.assignment_margin)

    separated_postings = Dict{UInt32, Vector{Pioneer._MBRIntensityPosting}}()
    for precursor_idx in UInt32(1):UInt32(128)
        separated_postings[precursor_idx] = [
            Pioneer._MBRIntensityPosting(UInt32(1), log2(100.0f0)),
        ]
    end
    for precursor_idx in UInt32(129):UInt32(256)
        separated_postings[precursor_idx] = [
            Pioneer._MBRIntensityPosting(UInt32(2), log2(100.0f0)),
        ]
    end
    separated_fit = Pioneer._fit_mbr_intensity_clusters(
        separated_postings,
        file_ids,
    )

    @test separated_fit.n_clusters == 2
    @test length(unique(separated_fit.cluster_by_precursor)) == 2
    @test all(x -> x ≈ 1.0f0, separated_fit.assignment_similarity)
    @test all(x -> x ≈ 2.0f0, separated_fit.assignment_margin)

end

@testset "MBR intensity profiles use sparse Pearson correlation" begin
    file_ids = UInt32[1, 2, 3, 4, 5, 6]
    postings = Dict(
        UInt32(1) => [
            Pioneer._MBRIntensityPosting(UInt32(2), 5.0f0),
        ],
        UInt32(2) => [
            Pioneer._MBRIntensityPosting(UInt32(1), 5.0f0),
            Pioneer._MBRIntensityPosting(UInt32(2), 5.0f0),
            Pioneer._MBRIntensityPosting(UInt32(3), 5.0f0),
        ],
        UInt32(3) => [
            Pioneer._MBRIntensityPosting(UInt32(2), 8.0f0),
        ],
        UInt32(4) => [
            Pioneer._MBRIntensityPosting(UInt32(1), 1.0f0),
            Pioneer._MBRIntensityPosting(UInt32(2), 2.0f0),
        ],
        UInt32(5) => [
            Pioneer._MBRIntensityPosting(UInt32(1), 10.0f0),
            Pioneer._MBRIntensityPosting(UInt32(2), 11.0f0),
        ],
    )
    profiles = Pioneer._mbr_intensity_profiles(postings, file_ids)

    function dense_profile(profile)
        values = fill(profile.baseline, length(file_ids))
        for idx in eachindex(profile.file_columns)
            values[profile.file_columns[idx]] += profile.values[idx]
        end
        return values
    end

    one_run = dense_profile(profiles[1])
    three_runs = dense_profile(profiles[2])
    @test sum(one_run) ≈ 0.0f0 atol = 1.0f-6
    @test sum(abs2, one_run) ≈ 1.0f0 atol = 1.0f-6
    @test sum(three_runs) ≈ 0.0f0 atol = 1.0f-6
    @test sum(abs2, three_runs) ≈ 1.0f0 atol = 1.0f-6
    @test Pioneer._mbr_profile_dot(profiles[1], profiles[2], length(file_ids)) ≈
          inv(sqrt(5.0f0)) atol = 1.0f-6
    @test Pioneer._mbr_profile_dot(profiles[1], profiles[3], length(file_ids)) ≈
          1.0f0 atol = 1.0f-6

    raw_left = Float32[1, 2, 0, 0, 0, 0]
    raw_right = Float32[10, 11, 0, 0, 0, 0]
    centered_left = raw_left .- sum(raw_left) / length(raw_left)
    centered_right = raw_right .- sum(raw_right) / length(raw_right)
    expected_raw_correlation = sum(centered_left .* centered_right) /
                               sqrt(sum(abs2, centered_left) * sum(abs2, centered_right))
    @test Pioneer._mbr_profile_dot(profiles[4], profiles[5], length(file_ids)) ≈
          expected_raw_correlation atol = 1.0f-6

    cluster = Pioneer._mbr_spherical_cluster(
        profiles,
        Int[1, 2],
        zeros(Float64, length(file_ids)),
        Int32[],
        zeros(Float32, length(file_ids)),
    )
    expected_centroid = one_run + three_runs
    expected_centroid ./= sqrt(sum(abs2, expected_centroid))
    dense_centroid_corrections = zeros(Float32, length(file_ids))
    Pioneer._mbr_set_dense_centroid!(dense_centroid_corrections, cluster.centroid)
    @test Pioneer._mbr_profile_centroid_dot(
        profiles[1],
        dense_centroid_corrections,
        cluster.centroid,
        length(file_ids),
    ) ≈ sum(one_run .* expected_centroid) atol = 1.0f-6
end

@testset "MBR intensity clusters measure peer agreement" begin
    passed_by_file = Dict(
        UInt32(1) => BitSet(vcat(1:15)),
        UInt32(2) => BitSet(vcat(1:9, 13:18)),
    )
    intensity_postings = Dict{UInt32, Vector{Pioneer._MBRIntensityPosting}}()
    for precursor_idx in UInt32(1):UInt32(6)
        intensity_postings[precursor_idx] = [
            Pioneer._MBRIntensityPosting(UInt32(1), log2(100.0f0)),
            Pioneer._MBRIntensityPosting(UInt32(2), log2(100.0f0)),
        ]
    end
    for precursor_idx in UInt32(7):UInt32(9)
        intensity_postings[precursor_idx] = [
            Pioneer._MBRIntensityPosting(UInt32(1), log2(900.0f0)),
            Pioneer._MBRIntensityPosting(UInt32(2), log2(100.0f0)),
        ]
    end
    for precursor_idx in UInt32(10):UInt32(12)
        intensity_postings[precursor_idx] = [
            Pioneer._MBRIntensityPosting(UInt32(1), log2(900.0f0)),
        ]
    end
    for precursor_idx in UInt32(13):UInt32(15)
        intensity_postings[precursor_idx] = [
            Pioneer._MBRIntensityPosting(UInt32(1), log2(100.0f0)),
            Pioneer._MBRIntensityPosting(UInt32(2), log2(900.0f0)),
        ]
    end
    for precursor_idx in UInt32(16):UInt32(18)
        intensity_postings[precursor_idx] = [
            Pioneer._MBRIntensityPosting(UInt32(2), log2(900.0f0)),
        ]
    end

    similarity = Pioneer._build_mbr_run_similarity_from_passed(
        passed_by_file;
        intensity_postings = intensity_postings,
        target_cluster_size = 6,
        min_child_size = 2,
        min_cluster_gain = 0.001f0,
    )

    human_cluster = similarity.cluster_by_precursor[1]
    yeast_cluster = similarity.cluster_by_precursor[7]
    ecoli_cluster = similarity.cluster_by_precursor[13]
    @test length(unique(similarity.cluster_by_precursor)) == 3
    @test all(==(human_cluster), similarity.cluster_by_precursor[1:6])
    @test all(==(yeast_cluster), similarity.cluster_by_precursor[7:12])
    @test all(==(ecoli_cluster), similarity.cluster_by_precursor[13:18])
    @test length(unique((human_cluster, yeast_cluster, ecoli_cluster))) == 3

    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(7),
        UInt32(1),
        UInt32(2),
    ) == 0.4f0
    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(7),
        UInt32(2),
        UInt32(1),
    ) == 0.4f0
    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(13),
        UInt32(1),
        UInt32(2),
    ) == 0.4f0
    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(13),
        UInt32(2),
        UInt32(1),
    ) == 0.4f0
    @test Pioneer._mbr_cluster_active_fraction(
        similarity,
        UInt32(7),
        UInt32(1),
        UInt32(2),
    ) == 1.0f0
    @test Pioneer._mbr_cluster_shared_support(
        similarity,
        UInt32(7),
        UInt32(1),
        UInt32(2),
    ) == 3.0f0
    @test Pioneer._mbr_cluster_shared_support(
        similarity,
        UInt32(7),
        UInt32(2),
        UInt32(1),
    ) == 3.0f0
    @test Pioneer._mbr_cluster_peer_count(similarity, UInt32(7)) == 5.0f0
    @test Pioneer._mbr_cluster_n_runs_observed(similarity, UInt32(10)) == 1.0f0
    @test all(>(0.0f0), similarity.cluster_assignment_similarity)
end

@testset "MBR cluster agreement excludes the candidate precursor" begin
    cluster_files = Dict(
        UInt32(1) => UInt32[2, 4, 6, 7],
        UInt32(2) => UInt32[2, 3, 4, 6],
        UInt32(3) => UInt32[2],
        UInt32(4) => UInt32[2],
    )
    cluster = Pioneer._build_mbr_cluster_run_similarity(
        cluster_files,
        UInt32[1, 2, 3, 4, 5, 6, 7],
    )
    passed_by_file = Dict(
        file_idx => BitSet(
            Int(precursor_idx)
            for (precursor_idx, run_files) in cluster_files
            if file_idx in run_files
        )
        for file_idx in UInt32(1):UInt32(7)
    )
    similarity = Pioneer._MBRRunSimilarity(
        Dict{Tuple{UInt32, UInt32}, Float32}(),
        Dict{Tuple{UInt32, UInt32}, Float32}(),
        Dict{Tuple{UInt32, UInt32}, Float32}(),
        Dict{UInt32, Float32}(),
        Dict{UInt32, Float32}(),
        BitSet(1:7),
        passed_by_file,
        fill(UInt32(1), 4),
        ones(Float32, 4),
        ones(Float32, 4),
        ones(UInt16, 4),
        Pioneer._MBRClusterRunSimilarity[cluster],
    )

    # All three peers are present in the donor and absent in the receiver.
    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(1),
        UInt32(1),
        UInt32(2),
    ) == 0.0f0
    @test Pioneer._mbr_cluster_active_fraction(
        similarity,
        UInt32(1),
        UInt32(1),
        UInt32(2),
    ) == 1.0f0

    # One peer is present in both runs; the other two are absent in both.
    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(1),
        UInt32(3),
        UInt32(4),
    ) == 1.0f0
    @test Pioneer._mbr_cluster_active_fraction(
        similarity,
        UInt32(1),
        UInt32(3),
        UInt32(4),
    ) ≈ 1.0f0 / 3.0f0

    # One peer is donor-only; the other two are absent in both runs.
    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(1),
        UInt32(5),
        UInt32(6),
    ) ≈ 2.0f0 / 3.0f0
    @test Pioneer._mbr_cluster_active_fraction(
        similarity,
        UInt32(1),
        UInt32(5),
        UInt32(6),
    ) ≈ 1.0f0 / 3.0f0

    # Shared absence is agreement, while active fraction shows no peer support.
    @test Pioneer._mbr_cluster_agreement(
        similarity,
        UInt32(1),
        UInt32(1),
        UInt32(7),
    ) == 1.0f0
    @test Pioneer._mbr_cluster_active_fraction(
        similarity,
        UInt32(1),
        UInt32(1),
        UInt32(7),
    ) == 0.0f0
    @test Pioneer._mbr_cluster_peer_count(similarity, UInt32(1)) == 3.0f0
end

@testset "MBR weighted containment does not penalize donor-only IDs during donor selection" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver.arrow")
        donor_broad_path = joinpath(dir, "donor_broad.arrow")
        donor_focused_path = joinpath(dir, "donor_focused.arrow")

        function write_similarity_table(path, file_idx, pids)
            Arrow.write(path, DataFrame(
                precursor_idx = UInt32.(pids),
                ms_file_idx = fill(UInt32(file_idx), length(pids)),
                target = trues(length(pids)),
                qval = fill(0.001f0, length(pids)),
            ))
        end

        write_similarity_table(receiver_path, 1, [1, 2])
        write_similarity_table(donor_broad_path, 2, [1, 2, 3, 4, 5])
        write_similarity_table(donor_focused_path, 3, [1])

        similarity = Pioneer.build_mbr_run_similarity(
            [receiver_path, donor_broad_path, donor_focused_path];
            q_value_threshold = 0.01f0,
        )

        @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(2)) == 1.0f0
        @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(3)) == 0.0f0
        donor_only_weight = 3 * log(4 / 2)
        expected_reverse = Float32(log(4 / 3) / (log(4 / 3) + donor_only_weight))
        @test Pioneer._mbr_run_similarity(similarity, UInt32(2), UInt32(1)) ≈ expected_reverse

        donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
            UInt32(10) => [
                _test_mbr_donor(0.99f0, UInt32(2); precursor_idx = UInt32(10)),
                _test_mbr_donor(0.70f0, UInt32(3); precursor_idx = UInt32(10)),
            ],
        )

        selected = Pioneer._donor_for_pid(donor_dict, UInt32(10), UInt32(1), similarity)

        @test selected !== nothing
        @test selected.ms_file_idx == UInt32(2)
    end
end

@testset "MBR IDF containment uses sparse complements for common precursors" begin
    passed_by_file = Dict{UInt32, BitSet}()
    for file_idx in UInt32(1):UInt32(100)
        ids = BitSet((1, 1000 + Int(file_idx)))
        file_idx < UInt32(100) && push!(ids, 2)
        passed_by_file[file_idx] = ids
    end

    similarity = Pioneer._build_mbr_run_similarity_from_passed(passed_by_file)
    common_idf = log(101 / 100)
    unique_idf = log(101 / 2)
    expected_shared = Float32(common_idf / (common_idf + unique_idf))

    @test isempty(similarity.coverage)
    @test isempty(similarity.shared_weight)
    @test length(similarity.missing_weight) == 99
    @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(2)) ≈ expected_shared
    @test Pioneer._mbr_run_similarity(similarity, UInt32(1), UInt32(100)) == 0.0f0
    @test Pioneer._mbr_run_similarity(similarity, UInt32(100), UInt32(1)) == 0.0f0
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
