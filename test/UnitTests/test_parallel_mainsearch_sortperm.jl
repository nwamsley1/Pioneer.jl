if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
end

@testset "parallel MainSearch permutations" begin
    workspace = Pioneer.Int32SortPermWorkspace()

    @testset "matches Base ordering and reuses Int32 buffers" begin
        precursor_ids = UInt32[9, 2, 4, 2, 9, 1, 4, 2]
        expected = Int32.(sortperm(precursor_ids))
        observed = Pioneer.parallel_sortperm_int32!(workspace, precursor_ids)

        @test observed === workspace.perm
        @test eltype(observed) === Int32
        @test observed == expected
        @test precursor_ids[observed] == sort(precursor_ids)

        scores = Float32[0.5, 0.9, 0.5, 0.2, 0.9, 0.5, 0.2]
        expected_desc = Int32.(sortperm(scores; rev = true, alg = QuickSort))
        observed_desc = Pioneer.parallel_sortperm_int32!(workspace, scores; rev = true)

        @test observed_desc === observed
        @test observed_desc == expected_desc
        @test length(workspace.temp) == length(scores)
    end

    @testset "precomputed Int32 order preserves PEP results" begin
        scores = Float32[0.5, 0.9, 0.5, 0.2, 0.9, 0.5, 0.2]
        targets = Bool[true, false, false, true, true, false, true]
        expected_peps = similar(scores)
        observed_peps = similar(scores)

        Pioneer.get_PEP!(scores, targets, expected_peps)
        order = Pioneer.parallel_sortperm_int32!(workspace, scores; rev = true)
        Pioneer._get_PEP_from_order!(scores, targets, observed_peps, order, 1.0f0)

        @test observed_peps == expected_peps

        mainsearch_peps, pass_mask = Pioneer._mainsearch_peps_and_pass_mask(
            scores,
            targets,
            workspace;
            pep_threshold = 0.5f0,
        )
        @test mainsearch_peps == expected_peps
        @test pass_mask == (expected_peps .<= 0.5f0)
    end

    @testset "precursor permutation keeps every column in lockstep" begin
        psms = DataFrame(
            precursor_idx = UInt32[3, 1, 2, 1, 3, 2, 1],
            scan_idx = UInt32[30, 10, 20, 11, 31, 21, 12],
            score = Float32[3, 1, 2, 1.1, 3.1, 2.1, 1.2],
            target = Bool[true, false, true, true, false, true, false],
            fragments = [[Float32(i), Float32(i + 1)] for i in 1:7],
        )
        expected = psms[sortperm(psms.precursor_idx), :]

        result = Pioneer.permute_psms_by_precursor_idx!(psms, workspace)

        @test result === psms
        @test psms == expected
    end

    @testset "empty input" begin
        observed = Pioneer.parallel_sortperm_int32!(workspace, Float32[])
        @test observed === workspace.perm
        @test isempty(observed)
        @test isempty(workspace.temp)
    end
end
