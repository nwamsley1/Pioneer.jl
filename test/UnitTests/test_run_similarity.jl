# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

@testset "Run similarity atlas" begin
    @testset "IDF-weighted directional containment" begin
        observed_ids_by_run = Dict(
            UInt32(1) => UInt32[1, 2, 3],
            UInt32(2) => UInt32[1, 2],
            UInt32(3) => UInt32[1],
            UInt32(4) => UInt32[4],
        )
        atlas = Pioneer.build_run_similarity(observed_ids_by_run)
        idf1 = log(5.0 / 4.0)
        idf2 = log(5.0 / 3.0)
        idf3 = log(5.0 / 2.0)

        @test Pioneer.run_similarity(
            atlas,
            UInt32(1),
            UInt32(2),
        ) ≈ Float32((idf1 + idf2) / (idf1 + idf2 + idf3))
        @test Pioneer.run_similarity(
            atlas,
            UInt32(2),
            UInt32(1),
        ) ≈ 1.0f0
        @test Pioneer.run_similarity(
            atlas,
            UInt32(4),
            UInt32(1),
        ) == 0.0f0
        @test Pioneer.is_observed_in_run(atlas, 2, UInt32(1))
        @test !Pioneer.is_observed_in_run(atlas, 2, UInt32(3))

        # The atlas owns a snapshot that can safely survive its builder.
        empty!(observed_ids_by_run[UInt32(1)])
        @test Pioneer.is_observed_in_run(atlas, 2, UInt32(1))
    end

    @testset "complement storage and centrality" begin
        atlas = Pioneer.build_run_similarity(Dict(
            UInt32(1) => UInt32[10, 20],
            UInt32(2) => UInt32[10, 20],
            UInt32(3) => UInt32[10, 20],
            UInt32(4) => UInt32[10, 20],
            UInt32(5) => UInt32[20],
        ))

        @test isempty(atlas.shared_weight)
        @test !isempty(atlas.missing_weight)
        @test Pioneer.run_similarity(
            atlas,
            UInt32(1),
            UInt32(2),
        ) == 1.0f0
        @test Pioneer.run_similarity(
            atlas,
            UInt32(1),
            UInt32(5),
        ) == 0.0f0
        @test Pioneer.run_centrality(atlas, UInt32(1)) ≈ 0.75f0
        @test Pioneer.run_centrality(atlas, UInt32(5)) == 0.0f0
        @test !Pioneer.is_observed_in_run(atlas, 20, UInt32(1))
        @test !Pioneer.is_observed_in_run(atlas, 20, UInt32(5))
    end

    @testset "empty and single-run atlases" begin
        empty_atlas = Pioneer.build_run_similarity(
            Dict{UInt32, Vector{UInt32}}(),
        )
        @test isempty(empty_atlas.run_ids)
        @test Pioneer.run_similarity(
            empty_atlas,
            UInt32(1),
            UInt32(2),
        ) == 0.0f0

        single_atlas = Pioneer.build_run_similarity(Dict(
            UInt32(7) => UInt32[1, 2],
        ))
        @test Pioneer.run_centrality(single_atlas, UInt32(7)) == 0.0f0
        @test Pioneer.run_similarity(
            single_atlas,
            UInt32(7),
            UInt32(7),
        ) == 0.0f0
    end
end
