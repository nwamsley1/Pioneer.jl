using Test
using DataFrames
using Statistics
using Arrow
using Distributions
using Pioneer

@testset "Protein Group Weight Coverage Features" begin
    empty_precursor_consensus = (
        relative_weight = Dict{Tuple{String, Bool, UInt8, UInt32}, Float32}(),
        mean_relative_weight = Dict{Tuple{String, Bool, UInt8}, Float32}(),
        profiled_precursor_count = Dict{Tuple{String, Bool, UInt8}, Int32}(),
        precursors_by_protein = Dict{Tuple{String, Bool, UInt8}, Vector{Pair{UInt32, Float32}}}()
    )

    @testset "Quant Necessary Columns Include Weight" begin
        cols_nombr = Pioneer.get_quant_necessary_columns(false)
        cols_mbr = Pioneer.get_quant_necessary_columns(true)

        @test :weight in cols_nombr
        @test :weight in cols_mbr
    end

    @testset "Group PSMs By Protein Uses Max Weight Within Peptide" begin
        psms = DataFrame(
            inferred_protein_group = ["P1", "P1", "P1", "P2"],
            target = Bool[true, true, true, false],
            entrap_id = UInt8[1, 1, 1, 1],
            precursor_idx = UInt32[11, 12, 13, 21],
            base_pep_id = UInt32[101, 101, 102, 201],
            sequence = ["PEP1", "PEP1", "PEP2", "PEP3"],
            structural_mods = ["", "", "", ""],
            isotopic_mods = ["", "", "", ""],
            use_for_protein_quant = Bool[true, true, true, true],
            missed_cleavage = Int64[0, 1, 0, 0],
            Mox = Int64[0, 0, 0, 0],
            prec_prob = Float32[0.8, 0.9, 0.6, 0.7],
            weight = Float32[100.0, 80.0, 50.0, 25.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = empty_precursor_consensus)

        @test hasproperty(grouped, :top_pep_weight)

        p1 = grouped[(grouped.protein_name .== "P1") .& grouped.target, :]
        @test nrow(p1) == 1
        @test p1.n_peptides[1] == 2
        @test p1.top_pep_weight[1] == 100.0f0
    end

    @testset "Protein Roll-Up Combines Precursors To Modified Peptides To Peptides" begin
        psms = DataFrame(
            inferred_protein_group = ["P1", "P1", "P1", "P1"],
            target = Bool[true, true, true, true],
            entrap_id = UInt8[1, 1, 1, 1],
            precursor_idx = UInt32[11, 12, 13, 14],
            base_pep_id = UInt32[101, 101, 101, 102],
            sequence = ["PEP1", "PEP1", "PEP1", "PEP2"],
            structural_mods = ["", "", "(3,M,Oxidation)", ""],
            isotopic_mods = ["", "", "", ""],
            use_for_protein_quant = Bool[true, true, true, true],
            missed_cleavage = Int64[0, 0, 0, 0],
            Mox = Int64[0, 0, 1, 0],
            prec_prob = Float32[0.80, 0.90, 0.50, 0.60],
            weight = Float32[100.0, 80.0, 60.0, 50.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = empty_precursor_consensus)

        @test grouped.n_peptides[1] == 2
        @test grouped.peptide_list[1] == "PEP1;PEP2"
        expected_pg_score = Float32(-log((0.2 + 0.01) * (0.1 + 0.01) * (0.5 + 0.01) * (0.4 + 0.01)))
        @test grouped.pg_score[1] ≈ expected_pg_score atol = 1e-5
    end

    @testset "Protein Roll-Up Counts Peptides By Sequence" begin
        psms = DataFrame(
            inferred_protein_group = ["P1", "P1", "P1"],
            target = Bool[true, true, true],
            entrap_id = UInt8[1, 1, 1],
            precursor_idx = UInt32[11, 12, 13],
            base_pep_id = UInt32[101, 202, 303],
            sequence = ["PEP1", "PEP1", "PEP2"],
            structural_mods = ["", "(3,M,Oxidation)", ""],
            isotopic_mods = ["", "", ""],
            use_for_protein_quant = Bool[true, true, true],
            missed_cleavage = Int64[0, 0, 0],
            Mox = Int64[0, 1, 0],
            prec_prob = Float32[0.80, 0.50, 0.60],
            weight = Float32[100.0, 60.0, 50.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = empty_precursor_consensus)

        @test grouped.n_peptides[1] == 2
        @test grouped.peptide_list[1] == "PEP1;PEP2"
        expected_pg_score = Float32(-log((0.2 + 0.01) * (0.5 + 0.01) * (0.4 + 0.01)))
        @test grouped.pg_score[1] ≈ expected_pg_score atol = 1e-5
    end

    @testset "MaxLFQ Protein Support Counts Distinguish Precursors Modified Peptides And Sequences" begin
        peptides = Union{Missing, Vector{Union{Missing, UInt32}}}[
            Union{Missing, UInt32}[UInt32(1), UInt32(2), UInt32(3), missing],
            missing,
        ]
        precursor_sequences = ["PEP1", "PEP1", "PEP2"]
        precursor_structural_mods = Union{Missing, String}["", "(3,M,Oxidation)", ""]
        precursor_isotopic_mods = Union{Missing, String}["", "", ""]

        n_precursors, n_modified_peptides, n_peptides = Pioneer.countProteinSupport(
            peptides,
            precursor_sequences,
            precursor_structural_mods,
            precursor_isotopic_mods
        )

        @test n_precursors[1] == 3
        @test n_modified_peptides[1] == 3
        @test n_peptides[1] == 2
        @test ismissing(n_precursors[2])
        @test ismissing(n_modified_peptides[2])
        @test ismissing(n_peptides[2])
    end

    @testset "Grouped Protein Catalog Uses Union Peptide Set" begin
        protein_catalog = Dict(
            (protein_name = "A", target = true, entrap_id = UInt8(1)) => Set(["p1", "p2"]),
            (protein_name = "B", target = true, entrap_id = UInt8(1)) => Set(["p2", "p3"])
        )
        df = DataFrame(
            protein_name = ["A;B"],
            target = Bool[true],
            entrap_id = UInt8[1],
            n_peptides = Int64[2]
        )

        (_, op) = Pioneer.add_protein_features(protein_catalog)
        out = op(copy(df))

        @test out.n_possible_peptides[1] == 3
        @test out.peptide_coverage[1] ≈ (2f0 / 3f0)
    end

    @testset "Weight Calibration Uses Max Peptide Weights" begin
        psms = DataFrame(
            use_for_protein_quant = Bool[true, true, true, true, true, false],
            sequence = ["P1A", "P1A", "P1B", "P2A", "P2B", "P3A"],
            inferred_protein_group = Union{Missing, String}["P1", "P1", "P1", "P2", "P2", missing],
            target = Bool[true, true, true, true, true, true],
            entrap_id = UInt8[1, 1, 1, 1, 1, 1],
            weight = Float32[60.0f0, 40.0f0, 40.0f0, 80.0f0, 10.0f0, 500.0f0]
        )

        model = Pioneer.estimate_weight_detection_model(psms)

        expected_threshold = Float32(quantile(log.([60.0, 40.0, 80.0, 10.0]), 0.05))
        @test model.log_threshold ≈ expected_threshold
        @test length(model.rank_scale_profile) == model.profiled_rank_count
        @test all(x -> x > 0.0f0, model.rank_scale_profile)
        @test model.profiled_rank_count == 1
    end

    @testset "Coverage Log Ratio Uses Exceedance Strength As Expected" begin
        calibration = (
            log_threshold = Float32(log(10.0)),
            rank_drop_profile = Float32[0.0f0, -0.1f0, -0.8f0, -1.5f0],
            rank_scale_profile = Float32[1.0f0, 0.25f0, 0.35f0, 0.45f0],
            profiled_rank_count = 4
        )
        pg_df = DataFrame(
            protein_name = ["P_big", "P_small", "P_balanced", "P_invalid", "P_short"],
            target = Bool[true, true, true, true, true],
            entrap_id = UInt8[1, 1, 1, 1, 1],
            n_peptides = Int64[1, 1, 4, 1, 1],
            pg_score = Float32[2.0, 1.5, 0.5, 3.0, 1.0],
            n_possible_peptides = Int64[20, 2, 4, 12, 1],
            top_pep_weight = Float32[100.0, 100.0, 100.0, 0.0, 100.0]
        )

        (_, op) = Pioneer.add_weight_observation_features(calibration)
        out = op(copy(pg_df))

        for col in (
            :expected_excess_from_top,
            :coverage_log_ratio
        )
            @test hasproperty(out, col)
        end

        @test out.expected_excess_from_top[1] > out.expected_excess_from_top[2]
        @test out.coverage_log_ratio[1] < out.coverage_log_ratio[2]
        @test out.coverage_log_ratio[3] > out.coverage_log_ratio[1]
        @test out.coverage_log_ratio[3] > out.coverage_log_ratio[2]

        # Invalid weight -> neutral values
        @test out.expected_excess_from_top[4] == 0.0f0
        @test out.coverage_log_ratio[4] == 0.0f0

        # N <= 1 -> neutral values
        @test out.expected_excess_from_top[5] == 0.0f0
        @test out.coverage_log_ratio[5] == 0.0f0
    end

    @testset "Protein Auto-Pass pg_score Threshold Can Be Disabled" begin
        df = DataFrame(
            target = Bool[true, false, false, false],
            n_peptides = Int64[1, 1, 2, 7],
            pg_score = Float32[0.9, 1.1, 0.7, 2.3],
            coverage_log_ratio = Float32[0.9, 0.8, 0.7, 0.1],
            precursor_consensus_prefix_shape = Float32[0.8, 0.7, 0.6, 0.5],
            any_common_peps = Bool[true, true, false, true]
        )

        pg_score_threshold = Pioneer.compute_protein_autopass_pg_score_threshold(df; context = "test")

        @test isinf(pg_score_threshold)
    end

    @testset "Protein Semi-Supervised Training Set Mines Failed Targets By Shape Or PEP Threshold" begin
        scores = Float32[
            0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90,
            0.89, 0.88, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82, 0.81, 0.80,
            0.79, 0.78, 0.77, 0.76, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15,
            0.14, 0.13, 0.12, 0.11, 0.10, 0.09
        ]
        prefix_shape = Float32[
            0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00, -0.02, -0.05, -0.08,
            -0.09, -0.095, -0.10, -0.21, -0.22, -0.23, -0.24, -0.55, 0.30, 0.30,
            0.30, 0.30, -0.55, 0.20, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30,
            0.30, 0.30, 0.30, 0.30, 0.30, 0.30
        ]
        n_peptides = Int64[
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 1, 2, 2,
            2, 2, 1, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2
        ]
        targets = Bool[
            true, true, true, true, true, true, true, true, true, true,
            true, true, true, true, true, true, true, true, false, false,
            false, false, true, true, false, false, false, false, false, false,
            false, false, false, false, false, false
        ]
        ss = Pioneer.build_protein_semisupervised_training_set(
            scores,
            targets,
            prefix_shape,
            n_peptides;
            q_value_threshold = 0.01f0
        )
        ss_loose_shape = Pioneer.build_protein_semisupervised_training_set(
            scores,
            targets,
            prefix_shape,
            n_peptides;
            q_value_threshold = 0.01f0,
            mined_negative_prefix_shape_threshold = 1.0f0,
            mined_negative_pep_threshold = 1.0f0
        )

        @test any(ss.positive_mask)
        @test any(ss.confident_positive_mask)
        @test any(ss.mined_negative_mask)
        @test ss.mined_negative_prefix_shape_threshold == -0.20f0
        @test ss.mined_negative_pep_threshold == 0.90f0
        @test all(ss.keep_mask[.!targets])
        @test all(ss.keep_mask)
        @test all(ss.confident_positive_mask .<= ss.positive_mask)
        @test sum(ss.labels) == sum(ss.positive_mask[ss.keep_mask])
        @test all(.!ss.labels[.!targets[ss.keep_mask]])
        @test all((prefix_shape[ss.mined_negative_mask] .<= ss.mined_negative_prefix_shape_threshold) .|
                  (ss.peps[ss.mined_negative_mask] .>= ss.mined_negative_pep_threshold))
        @test ss.positive_mask == (targets .& .!ss.mined_negative_mask)
        @test ss.requested_mined_negative_count == sum(ss.mined_negative_mask)
        @test sum(ss.mined_negative_mask) == ss.requested_mined_negative_count
        @test all(ss.peps[ss.confident_positive_mask] .<= 1.0f0)
        @test sum(ss.mined_negative_mask) <= sum(ss_loose_shape.mined_negative_mask)
    end

    @testset "Protein Semi-Supervised Training Set Singleton Low-Shape Targets Are Mined Negative Even If Q-Passing" begin
        scores = Float32[0.99, 0.76, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90]
        targets = Bool[true, true, false, false, false, false, false, false]
        prefix_shape = Float32[-0.25, -0.25, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30]
        n_peptides = Int64[1, 2, 2, 2, 2, 2, 2, 2]

        ss = Pioneer.build_protein_semisupervised_training_set(
            scores,
            targets,
            prefix_shape,
            n_peptides;
            q_value_threshold = 0.01f0
        )

        @test ss.qvals[1] <= 0.01f0
        @test !ss.confident_positive_mask[1]
        @test !ss.positive_mask[1]
        @test ss.mined_negative_mask[1]
        @test ss.qvals[2] > 0.01f0
        @test ss.mined_negative_mask[2]
        @test !ss.confident_positive_mask[2]
    end

    @testset "Protein Semi-Supervised Training Set Keeps Q-Passing Low-Shape Multi-Peptide Targets Positive" begin
        scores = Float32[0.99, 0.76, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90]
        targets = Bool[true, true, false, false, false, false, false, false]
        prefix_shape = Float32[-0.25, -0.25, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30]
        n_peptides = Int64[2, 2, 2, 2, 2, 2, 2, 2]

        ss = Pioneer.build_protein_semisupervised_training_set(
            scores,
            targets,
            prefix_shape,
            n_peptides;
            q_value_threshold = 0.01f0
        )

        @test ss.qvals[1] <= 0.01f0
        @test ss.confident_positive_mask[1]
        @test ss.positive_mask[1]
        @test !ss.mined_negative_mask[1]
        @test ss.qvals[2] > 0.01f0
        @test ss.mined_negative_mask[2]
    end

    @testset "Protein Semi-Supervised Training Set Rescues High-Shape Failing Targets On Initial Labels" begin
        scores = Float32.(collect(reverse(range(0.99, 0.70, length = 30))))
        prefix_shape = Float32.(vcat(fill(0.30, 15), fill(0.05, 10), fill(-0.25, 5)))
        n_peptides = fill(Int64(2), 30)
        targets = trues(30)
        targets[[1, 9, 17]] .= false

        ss = Pioneer.build_protein_semisupervised_training_set(
            scores,
            targets,
            prefix_shape,
            n_peptides;
            q_value_threshold = 0.01f0
        )

        rescued_positive_mask = targets .& .!ss.confident_positive_mask .& .!ss.mined_negative_mask

        @test ss.mined_negative_prefix_shape_threshold == -0.20f0
        @test ss.mined_negative_pep_threshold == 0.90f0
        @test all(prefix_shape[rescued_positive_mask] .> ss.mined_negative_prefix_shape_threshold)
        @test all(ss.peps[rescued_positive_mask] .< ss.mined_negative_pep_threshold)
        @test any(ss.mined_negative_mask)
        @test any(rescued_positive_mask)
        @test all(ss.positive_mask[rescued_positive_mask])
        @test all(ss.keep_mask[rescued_positive_mask])
        @test all(ss.keep_mask[.!targets])
    end

    @testset "Protein Semi-Supervised Training Set Drops High-Shape Failing Targets After Iteration One" begin
        scores = Float32.(collect(reverse(range(0.99, 0.70, length = 30))))
        prefix_shape = Float32.(vcat(fill(0.30, 15), fill(0.05, 10), fill(-0.25, 5)))
        n_peptides = fill(Int64(2), 30)
        targets = trues(30)
        targets[[1, 9, 17]] .= false

        ss = Pioneer.build_protein_semisupervised_training_set(
            scores,
            targets,
            prefix_shape,
            n_peptides;
            q_value_threshold = 0.01f0,
            keep_non_mined_targets_as_positive = false
        )

        dropped_target_mask = targets .& .!ss.confident_positive_mask .& .!ss.mined_negative_mask

        @test ss.mined_negative_prefix_shape_threshold == -0.20f0
        @test any(dropped_target_mask)
        @test all(.!ss.positive_mask[dropped_target_mask])
        @test all(.!ss.keep_mask[dropped_target_mask])
    end

    @testset "Consensus Relative Weight Builder Uses Quant Precursors And Current pg_score Weighting" begin
        temp_dir = mktempdir()

        try
            file1 = DataFrame(
                inferred_protein_group = ["P", "P", "P"],
                target = Bool[true, true, true],
                entrap_id = UInt8[1, 1, 1],
                precursor_idx = UInt32[1, 2, 3],
                base_pep_id = UInt32[101, 102, 103],
                sequence = ["PEP1", "PEP2", "PEP3"],
                structural_mods = ["", "", ""],
                isotopic_mods = ["", "", ""],
                use_for_protein_quant = Bool[true, true, true],
                MBR_candidate = Bool[false, false, true],
                missed_cleavage = Int64[0, 0, 0],
                Mox = Int64[0, 0, 0],
                prec_prob = Float32[0.99, 0.98, 0.10],
                MBR_boosted_prec_prob = Float32[0.99, 0.98, 0.995],
                weight = Float32[100.0, 50.0, 1000.0]
            )

            file2 = DataFrame(
                inferred_protein_group = ["P", "P"],
                target = Bool[true, true],
                entrap_id = UInt8[1, 1],
                precursor_idx = UInt32[1, 2],
                base_pep_id = UInt32[101, 102],
                sequence = ["PEP1", "PEP2"],
                structural_mods = ["", ""],
                isotopic_mods = ["", ""],
                use_for_protein_quant = Bool[true, true],
                MBR_candidate = Bool[false, false],
                missed_cleavage = Int64[0, 0],
                Mox = Int64[0, 0],
                prec_prob = Float32[0.10, 0.20],
                MBR_boosted_prec_prob = Float32[0.10, 0.20],
                weight = Float32[60.0, 120.0]
            )

            path1 = joinpath(temp_dir, "psms_001.arrow")
            path2 = joinpath(temp_dir, "psms_002.arrow")
            Arrow.write(path1, file1)
            Arrow.write(path2, file2)

            consensus = Pioneer.build_precursor_consensus([
                Pioneer.PSMFileReference(path1),
                Pioneer.PSMFileReference(path2)
            ])

            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(3))] >
                  consensus.relative_weight[("P", true, UInt8(1), UInt32(1))]
            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(3))] > 0.9f0
            @test haskey(consensus.relative_weight, ("P", true, UInt8(1), UInt32(3)))
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

    @testset "Consensus Relative Weight Builder Keeps Only Top Five Protein Runs" begin
        temp_dir = mktempdir()

        try
            refs = Pioneer.PSMFileReference[]
            for run_idx in 1:6
                df = if run_idx <= 5
                    DataFrame(
                        inferred_protein_group = ["P"],
                        target = Bool[true],
                        entrap_id = UInt8[1],
                        precursor_idx = UInt32[1],
                        base_pep_id = UInt32[101],
                        sequence = ["PEP1"],
                        structural_mods = [""],
                        isotopic_mods = [""],
                        use_for_protein_quant = Bool[true],
                        MBR_candidate = Bool[false],
                        missed_cleavage = Int64[0],
                        Mox = Int64[0],
                        prec_prob = Float32[0.95f0 - 0.01f0 * run_idx],
                        weight = Float32[100.0]
                    )
                else
                    DataFrame(
                        inferred_protein_group = ["P"],
                        target = Bool[true],
                        entrap_id = UInt8[1],
                        precursor_idx = UInt32[2],
                        base_pep_id = UInt32[102],
                        sequence = ["PEP2"],
                        structural_mods = [""],
                        isotopic_mods = [""],
                        use_for_protein_quant = Bool[true],
                        MBR_candidate = Bool[false],
                        missed_cleavage = Int64[0],
                        Mox = Int64[0],
                        prec_prob = Float32[0.05],
                        weight = Float32[100.0]
                    )
                end

                path = joinpath(temp_dir, "psms_$(lpad(run_idx, 3, '0')).arrow")
                Arrow.write(path, df)
                push!(refs, Pioneer.PSMFileReference(path))
            end

            consensus = Pioneer.build_precursor_consensus(refs)

            @test consensus.profiled_precursor_count[("P", true, UInt8(1))] == 1
            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(1))] == 1.0f0
            @test !haskey(consensus.relative_weight, ("P", true, UInt8(1), UInt32(2)))
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

    @testset "Consensus Relative Weight Builder Applies Exponential Decay Across Selected Runs" begin
        temp_dir = mktempdir()

        try
            run1 = DataFrame(
                inferred_protein_group = ["P", "P"],
                target = Bool[true, true],
                entrap_id = UInt8[1, 1],
                precursor_idx = UInt32[2, 1],
                base_pep_id = UInt32[102, 101],
                sequence = ["PEP2", "PEP1"],
                structural_mods = ["", ""],
                isotopic_mods = ["", ""],
                use_for_protein_quant = Bool[true, true],
                MBR_candidate = Bool[false, false],
                missed_cleavage = Int64[0, 0],
                Mox = Int64[0, 0],
                prec_prob = Float32[0.90, 0.90],
                weight = Float32[100.0, 80.0]
            )
            run2 = DataFrame(
                inferred_protein_group = ["P", "P"],
                target = Bool[true, true],
                entrap_id = UInt8[1, 1],
                precursor_idx = UInt32[1, 2],
                base_pep_id = UInt32[101, 102],
                sequence = ["PEP1", "PEP2"],
                structural_mods = ["", ""],
                isotopic_mods = ["", ""],
                use_for_protein_quant = Bool[true, true],
                MBR_candidate = Bool[false, false],
                missed_cleavage = Int64[0, 0],
                Mox = Int64[0, 0],
                prec_prob = Float32[0.90, 0.90],
                weight = Float32[100.0, 80.0]
            )

            path1 = joinpath(temp_dir, "psms_001.arrow")
            path2 = joinpath(temp_dir, "psms_002.arrow")
            Arrow.write(path1, run1)
            Arrow.write(path2, run2)

            consensus = Pioneer.build_precursor_consensus([
                Pioneer.PSMFileReference(path1),
                Pioneer.PSMFileReference(path2)
            ])

            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(2))] >
                  consensus.relative_weight[("P", true, UInt8(1), UInt32(1))]
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

    @testset "Consensus Relative Weight Support Uses Observed MBR Precursors" begin
        consensus = (
            relative_weight = Dict(
                ("P", true, UInt8(1), UInt32(11)) => 1.0f0,
                ("P", true, UInt8(1), UInt32(12)) => 0.2f0
            ),
            mean_relative_weight = Dict(
                ("P", true, UInt8(1)) => 0.4f0
            ),
            profiled_precursor_count = Dict(
                ("P", true, UInt8(1)) => Int32(5)
            ),
            shape_strength = Dict(
                ("P", true, UInt8(1)) => 1.0f0
            ),
            shape_confidence_scale = 1.0f0,
            precursors_by_protein = Dict(
                ("P", true, UInt8(1)) => Pair{UInt32, Float32}[UInt32(11) => 1.0f0, UInt32(12) => 0.2f0]
            )
        )

        psms = DataFrame(
            inferred_protein_group = ["P", "P"],
            target = Bool[true, true],
            entrap_id = UInt8[1, 1],
            precursor_idx = UInt32[11, 12],
            base_pep_id = UInt32[101, 102],
            sequence = ["PEP1", "PEP2"],
            structural_mods = ["", ""],
            isotopic_mods = ["", ""],
            use_for_protein_quant = Bool[true, true],
            MBR_candidate = Bool[false, true],
            missed_cleavage = Int64[0, 0],
            Mox = Int64[0, 0],
            prec_prob = Float32[0.8, 0.7],
            MBR_boosted_prec_prob = Float32[0.8, 0.85],
            weight = Float32[100.0, 90.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus)
        @test grouped.precursor_consensus_prefix_shape[1] ≈ 0.121988f0 atol = 1e-5
        @test grouped.pg_score_x_precursor_consensus_prefix_shape[1] ≈
              grouped.pg_score[1] * grouped.precursor_consensus_prefix_shape[1] atol = 1e-5
    end

    @testset "Consensus Prefix Features Reward Expected Singleton And Penalize Low-Ranked Singleton" begin
        consensus = (
            relative_weight = Dict(
                ("P_top", true, UInt8(1), UInt32(101)) => 1.0f0,
                ("P_top", true, UInt8(1), UInt32(102)) => 0.2f0,
                ("P_low", true, UInt8(1), UInt32(201)) => 1.0f0,
                ("P_low", true, UInt8(1), UInt32(202)) => 0.2f0
            ),
            mean_relative_weight = Dict(
                ("P_top", true, UInt8(1)) => 0.6f0,
                ("P_low", true, UInt8(1)) => 0.6f0
            ),
            profiled_precursor_count = Dict(
                ("P_top", true, UInt8(1)) => Int32(2),
                ("P_low", true, UInt8(1)) => Int32(2)
            ),
            shape_strength = Dict(
                ("P_top", true, UInt8(1)) => 1.0f0,
                ("P_low", true, UInt8(1)) => 1.0f0
            ),
            shape_confidence_scale = 1.0f0,
            precursors_by_protein = Dict(
                ("P_top", true, UInt8(1)) => Pair{UInt32, Float32}[UInt32(101) => 1.0f0, UInt32(102) => 0.2f0],
                ("P_low", true, UInt8(1)) => Pair{UInt32, Float32}[UInt32(201) => 1.0f0, UInt32(202) => 0.2f0]
            )
        )

        psms = DataFrame(
            inferred_protein_group = ["P_top", "P_low"],
            target = Bool[true, true],
            entrap_id = UInt8[1, 1],
            precursor_idx = UInt32[101, 202],
            base_pep_id = UInt32[101, 202],
            sequence = ["PEP_TOP", "PEP_LOW"],
            structural_mods = ["", ""],
            isotopic_mods = ["", ""],
            use_for_protein_quant = Bool[true, true],
            MBR_candidate = Bool[false, false],
            missed_cleavage = Int64[0, 0],
            Mox = Int64[0, 0],
            prec_prob = Float32[0.8, 0.8],
            weight = Float32[100.0, 100.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus)

        top = grouped[grouped.protein_name .== "P_top", :]
        low = grouped[grouped.protein_name .== "P_low", :]

        @test top.precursor_consensus_prefix_shape[1] ≈ 0.316060f0 atol = 1e-5
        @test low.precursor_consensus_prefix_shape[1] ≈ -0.210707f0 atol = 1e-5
        @test low.precursor_consensus_prefix_shape[1] < 0.0f0
    end

    @testset "Consensus Prefix Features Penalize Off-Consensus Observed Precursors" begin
        consensus = (
            relative_weight = Dict(
                ("P", true, UInt8(1), UInt32(101)) => 1.0f0,
                ("P", true, UInt8(1), UInt32(102)) => 0.2f0
            ),
            mean_relative_weight = Dict(
                ("P", true, UInt8(1)) => 0.6f0
            ),
            profiled_precursor_count = Dict(
                ("P", true, UInt8(1)) => Int32(2)
            ),
            shape_strength = Dict(
                ("P", true, UInt8(1)) => 1.0f0
            ),
            shape_confidence_scale = 1.0f0,
            precursors_by_protein = Dict(
                ("P", true, UInt8(1)) => Pair{UInt32, Float32}[UInt32(101) => 1.0f0, UInt32(102) => 0.2f0]
            )
        )

        psms = DataFrame(
            inferred_protein_group = ["P", "P"],
            target = Bool[true, true],
            entrap_id = UInt8[1, 1],
            precursor_idx = UInt32[101, 999],
            base_pep_id = UInt32[101, 999],
            sequence = ["PEP_TOP", "PEP_OFF"],
            structural_mods = ["", ""],
            isotopic_mods = ["", ""],
            use_for_protein_quant = Bool[true, true],
            MBR_candidate = Bool[false, false],
            missed_cleavage = Int64[0, 0],
            Mox = Int64[0, 0],
            prec_prob = Float32[0.8, 0.8],
            weight = Float32[100.0, 100.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus)

        @test grouped.precursor_consensus_prefix_shape[1] ≈ 0.0f0 atol = 1e-5
        @test grouped.precursor_consensus_prefix_shape[1] < 0.316060f0
    end

    @testset "Consensus Prefix Features Use Leave-One-Out Consensus By Run" begin
        protein_key = ("P", true, UInt8(1))
        consensus = (
            relative_weight = Dict(
                (protein_key[1], protein_key[2], protein_key[3], UInt32(101)) => 0.5f0,
                (protein_key[1], protein_key[2], protein_key[3], UInt32(102)) => 0.5f0
            ),
            mean_relative_weight = Dict(
                protein_key => 0.5f0
            ),
            profiled_precursor_count = Dict(
                protein_key => Int32(2)
            ),
            shape_strength = Dict(
                protein_key => 1.0f0
            ),
            shape_confidence_scale = 1.0f0,
            selected_run_votes = Dict(
                protein_key => [
                    (pg_score = 1.0f0, run_order = Int64(1), normalized_precursors = Pair{UInt32, Float32}[UInt32(101) => 1.0f0]),
                    (pg_score = 1.0f0, run_order = Int64(2), normalized_precursors = Pair{UInt32, Float32}[UInt32(102) => 1.0f0]),
                    (pg_score = 1.0f0, run_order = Int64(3), normalized_precursors = Pair{UInt32, Float32}[UInt32(103) => 1.0f0])
                ]
            ),
            consensus_target_run_count = Dict(
                protein_key => Int32(2)
            ),
            cached_consensus_weight_sums = Dict(
                protein_key => Dict(UInt32(101) => 1.0, UInt32(102) => 1.0, UInt32(103) => 1.0)
            ),
            cached_protein_total_vote = Dict(
                protein_key => 3.0
            ),
            precursors_by_protein = Dict(
                protein_key => Pair{UInt32, Float32}[UInt32(101) => 0.5f0, UInt32(102) => 0.5f0]
            )
        )

        psms = DataFrame(
            inferred_protein_group = ["P"],
            target = Bool[true],
            entrap_id = UInt8[1],
            precursor_idx = UInt32[101],
            base_pep_id = UInt32[101],
            sequence = ["PEP1"],
            structural_mods = [""],
            isotopic_mods = [""],
            use_for_protein_quant = Bool[true],
            MBR_candidate = Bool[false],
            missed_cleavage = Int64[0],
            Mox = Int64[0],
            prec_prob = Float32[0.8],
            weight = Float32[100.0]
        )

        grouped_default = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus)
        grouped_loo = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus, current_run_order = Int64(1))

        @test grouped_default.precursor_consensus_prefix_shape[1] == 0.0f0
        @test grouped_loo.precursor_consensus_prefix_shape[1] < 0.0f0
        @test grouped_loo.precursor_consensus_prefix_shape[1] < grouped_default.precursor_consensus_prefix_shape[1]
    end

    @testset "Decoy Shape Is Scored Against Itself" begin
        protein_name = "P1"
        consensus = (
            relative_weight = Dict(
                (protein_name, true, UInt8(1), UInt32(101)) => 1.0f0
            ),
            mean_relative_weight = Dict(
                (protein_name, true, UInt8(1)) => 1.0f0
            ),
            profiled_precursor_count = Dict(
                (protein_name, true, UInt8(1)) => Int32(1)
            ),
            shape_strength = Dict(
                (protein_name, true, UInt8(1)) => 1.0f0
            ),
            shape_confidence_scale = 1.0f0,
            precursors_by_protein = Dict(
                (protein_name, true, UInt8(1)) => Pair{UInt32, Float32}[UInt32(101) => 1.0f0]
            )
        )

        psms = DataFrame(
            inferred_protein_group = [protein_name],
            target = Bool[false],
            entrap_id = UInt8[1],
            precursor_idx = UInt32[201],
            base_pep_id = UInt32[201],
            sequence = ["DECOY_PEP"],
            structural_mods = [""],
            isotopic_mods = [""],
            use_for_protein_quant = Bool[true],
            MBR_candidate = Bool[false],
            missed_cleavage = Int64[0],
            Mox = Int64[0],
            prec_prob = Float32[0.8],
            weight = Float32[100.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus)

        @test grouped.precursor_consensus_prefix_shape[1] == 0.0f0
        @test !hasproperty(grouped, :training_precursor_consensus_prefix_shape)
    end

    @testset "Optional Probit Feature Columns Drop Cleanly" begin
        feature_names = Pioneer.protein_probit_feature_names()
        df = DataFrame(pg_score = Float32[0.1, 0.2, 0.3])

        Pioneer.remove_zero_variance_columns!(feature_names, df)

        @test feature_names == [:pg_score]
    end

    @testset "Protein Probit Feature Names Include Consensus Support" begin
        @test Pioneer.protein_probit_feature_names() == [
            :pg_score,
            :peptide_coverage_logit,
            :any_common_peps,
            :coverage_log_ratio,
            :precursor_consensus_prefix_shape,
            :pg_score_x_precursor_consensus_prefix_shape
        ]

        @test Pioneer.protein_probit_feature_names(include_n_possible_peptides = true) == [
            :pg_score,
            :peptide_coverage_logit,
            :n_possible_peptides,
            :any_common_peps,
            :coverage_log_ratio,
            :precursor_consensus_prefix_shape,
            :pg_score_x_precursor_consensus_prefix_shape
        ]
    end

    @testset "Protein Probit Fold Label Scatter Writes File" begin
        temp_dir = mktempdir()

        try
            df = DataFrame(
                file_idx = Int64[1, 2, 3],
                species = ["YEAST", "HUMAN", "YEAST"],
                pg_score = Float32[0.2, 0.5, 0.8],
                coverage_log_ratio = Float32[-0.7, -0.1, 0.4],
                precursor_consensus_prefix_shape = Float32[0.1, 0.3, 0.6],
                target = Bool[true, false, true],
                n_peptides = Int64[1, 4, 9]
            )
            ss = (
                positive_mask = Bool[true, false, false],
                mined_negative_mask = Bool[false, false, true]
            )

            scatter_path = Pioneer.write_protein_probit_fold_label_scatter(
                df,
                UInt8(2),
                temp_dir,
                ss;
                stage = "iter_1",
                context = "test",
                file_idx_to_name = Dict(1 => "Run_A", 2 => "Run_B", 3 => "Run_with_Yeast")
            )
            coverage_path = Pioneer.write_protein_probit_fold_coverage_scatter(
                df,
                UInt8(2),
                temp_dir,
                ss;
                stage = "iter_1",
                context = "test",
                file_idx_to_name = Dict(1 => "Run_A", 2 => "Run_B", 3 => "Run_with_Yeast")
            )
            @test scatter_path !== nothing
            @test isfile(scatter_path)
            @test coverage_path !== nothing
            @test isfile(coverage_path)
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

end
