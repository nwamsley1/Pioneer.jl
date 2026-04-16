@testset "MaxLFQ Sparse Graph Handling" begin
    function run_maxlfq_case(
        peptides::Vector{UInt32},
        experiments::Vector{UInt16},
        abundances::Vector{Union{Missing, Float32}};
        scores::Vector{Union{Missing, Float32}} = Union{Missing, Float32}[fill(1.0f0, length(peptides))...]
    )
        n_runs = length(unique(experiments))

        target_out = Vector{Union{Missing, Bool}}(missing, n_runs)
        entrap_id_out = Vector{Union{Missing, UInt8}}(missing, n_runs)
        species_out = Vector{Union{Missing, String}}(missing, n_runs)
        protein_out = Vector{Union{Missing, String}}(missing, n_runs)
        peptides_out = Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}(missing, n_runs)
        experiments_out = Vector{Union{Missing, UInt32}}(missing, n_runs)
        log2_abundance_out = Vector{Union{Missing, Float32}}(missing, n_runs)
        global_qval_out = Vector{Union{Missing, Float32}}(missing, n_runs)
        qval_out = Vector{Union{Missing, Float32}}(missing, n_runs)
        pep_out = Vector{Union{Missing, Float32}}(missing, n_runs)
        pg_score_out = Vector{Union{Missing, Float32}}(missing, n_runs)
        global_pg_score_out = Vector{Union{Missing, Float32}}(missing, n_runs)
        total_peak_area_out = Vector{Union{Missing, Float32}}(missing, n_runs)

        metrics = Union{Missing, Float32}[fill(0.001f0, length(peptides))...]

        Pioneer.getProtAbundance(
            "P",
            1,
            true,
            UInt8(1),
            "species",
            peptides,
            experiments,
            fill(true, length(peptides)),
            abundances,
            metrics,
            metrics,
            metrics,
            scores,
            scores,
            target_out,
            entrap_id_out,
            species_out,
            protein_out,
            peptides_out,
            experiments_out,
            log2_abundance_out,
            global_qval_out,
            qval_out,
            pep_out,
            pg_score_out,
            global_pg_score_out,
            total_peak_area_out
        )

        return (
            experiments = experiments_out,
            log2_abundance = log2_abundance_out,
            peptides = peptides_out,
            total_peak_area = total_peak_area_out
        )
    end

    @testset "Disconnected singleton run is missing" begin
        result = run_maxlfq_case(
            UInt32[1, 1, 2],
            UInt16[1, 2, 3],
            Union{Missing, Float32}[100.0f0, 200.0f0, 50.0f0]
        )

        @test result.experiments == Union{Missing, UInt32}[UInt32(1), UInt32(2), UInt32(3)]
        @test !ismissing(result.log2_abundance[1])
        @test !ismissing(result.log2_abundance[2])
        @test ismissing(result.log2_abundance[3])
        @test result.log2_abundance[2] - result.log2_abundance[1] ≈ 1.0f0 atol = 1f-4
        @test isequal(result.total_peak_area, Union{Missing, Float32}[100.0f0, 200.0f0, 50.0f0])
    end

    @testset "Sparse connected chain preserves pairwise ratios" begin
        result = run_maxlfq_case(
            UInt32[1, 1, 2, 2],
            UInt16[1, 2, 2, 3],
            Union{Missing, Float32}[100.0f0, 200.0f0, 300.0f0, 600.0f0]
        )

        @test all(x -> !ismissing(x), result.log2_abundance)
        @test result.log2_abundance[2] - result.log2_abundance[1] ≈ 1.0f0 atol = 1f-4
        @test result.log2_abundance[3] - result.log2_abundance[2] ≈ 1.0f0 atol = 1f-4
    end

    @testset "Tied largest component uses higher summed pg_score" begin
        result = run_maxlfq_case(
            UInt32[1, 1, 2, 2],
            UInt16[1, 2, 3, 4],
            Union{Missing, Float32}[100.0f0, 200.0f0, 300.0f0, 600.0f0];
            scores = Union{Missing, Float32}[1.0f0, 1.0f0, 5.0f0, 5.0f0]
        )

        @test ismissing(result.log2_abundance[1])
        @test ismissing(result.log2_abundance[2])
        @test !ismissing(result.log2_abundance[3])
        @test !ismissing(result.log2_abundance[4])
        @test result.log2_abundance[4] - result.log2_abundance[3] ≈ 1.0f0 atol = 1f-4
    end

    @testset "Dense graph remains quantified across runs" begin
        result = run_maxlfq_case(
            UInt32[1, 1, 1, 2, 2, 2],
            UInt16[1, 2, 3, 1, 2, 3],
            Union{Missing, Float32}[100.0f0, 200.0f0, 400.0f0, 10.0f0, 20.0f0, 40.0f0]
        )

        @test all(x -> !ismissing(x), result.log2_abundance)
        @test result.log2_abundance[2] - result.log2_abundance[1] ≈ 1.0f0 atol = 1f-4
        @test result.log2_abundance[3] - result.log2_abundance[2] ≈ 1.0f0 atol = 1f-4
    end

    @testset "Nonpositive intensities remain missing" begin
        result = run_maxlfq_case(
            UInt32[1, 1, 2],
            UInt16[1, 2, 3],
            Union{Missing, Float32}[0.0f0, missing, -5.0f0]
        )

        @test all(ismissing, result.log2_abundance)
    end

    @testset "Singleton protein falls back to total peak area" begin
        result = run_maxlfq_case(
            UInt32[1, 2],
            UInt16[7, 7],
            Union{Missing, Float32}[100.0f0, 25.0f0]
        )

        @test result.experiments == Union{Missing, UInt32}[UInt32(7)]
        @test !ismissing(result.log2_abundance[1])
        @test result.log2_abundance[1] ≈ log2(125.0f0) atol = 1f-4
        @test isequal(result.total_peak_area, Union{Missing, Float32}[125.0f0])
    end

    @testset "Arrow append schema stays nullable across LFQ batches" begin
        batch1_abundance = Pioneer.get_linear_abundance(Union{Missing, Float32}[1.0f0, 2.0f0])
        batch2_abundance = Pioneer.get_linear_abundance(Union{Missing, Float32}[1.0f0, missing])

        @test eltype(batch1_abundance) == Union{Missing, Float32}
        @test eltype(batch2_abundance) == Union{Missing, Float32}

        mktempdir() do tmpdir
            output_path = joinpath(tmpdir, "lfq.arrow")

            batch1 = DataFrame(abundance = batch1_abundance)
            batch2 = DataFrame(abundance = batch2_abundance)

            open(output_path, "w") do io
                Arrow.write(io, batch1; file = false)
            end
            Arrow.append(output_path, batch2)

            appended = DataFrame(Arrow.Table(output_path))
            @test isequal(appended.abundance, Union{Missing, Float32}[2.0f0, 4.0f0, 2.0f0, missing])
        end
    end

    @testset "Single-point precursor integrations do not contribute to LFQ" begin
        input_df = DataFrame(
            inferred_protein_group = ["P001", "P001", "P001"],
            target = Bool[true, true, true],
            entrapment_group_id = UInt8[0, 0, 0],
            species = ["human", "human", "human"],
            precursor_idx = UInt32[1, 1, 2],
            ms_file_idx = UInt16[1, 2, 2],
            use_for_protein_quant = Bool[true, true, false],
            peak_area = Union{Missing, Float32}[100.0f0, 200.0f0, 300.0f0],
            pg_qval = Union{Missing, Float32}[0.001f0, 0.001f0, 0.001f0],
            qlobal_pg_qval = Union{Missing, Float32}[0.001f0, 0.001f0, 0.001f0],
            pg_pep = Union{Missing, Float32}[0.001f0, 0.001f0, 0.001f0],
            pg_score = Union{Missing, Float32}[5.0f0, 5.0f0, 5.0f0],
            global_pg_score = Union{Missing, Float32}[5.0f0, 5.0f0, 5.0f0],
            points_integrated = Union{Missing, UInt32}[2, 1, 3]
        )

        filtered = Pioneer.apply_lfq_preprocessing(input_df, 0.01f0)

        @test nrow(filtered) == 1
        @test filtered.ms_file_idx == UInt16[1]
        @test filtered.precursor_idx == UInt32[1]
        @test filtered.points_integrated == Union{Missing, UInt32}[2]
    end
end
