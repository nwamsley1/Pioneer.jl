"""
Test suite for protein-level entrapment analysis functions
"""

using Test
using DataFrames
using Arrow
push!(LOAD_PATH, dirname(dirname(@__DIR__)))
using EntrapmentAnalysis

@testset "Protein Entrapment Analysis" begin
    
    @testset "assign_protein_entrapment_pairs!" begin
        # Test basic pairing with one entrapment group
        df = DataFrame(
            protein = ["PROT1", "PROT1", "PROT2", "PROT2"],
            entrap_id = UInt8[0, 1, 0, 1],
            file_name = ["file1", "file1", "file1", "file1"]
        )
        
        assign_protein_entrapment_pairs!(df)
        
        @test hasproperty(df, :entrap_pair_id)
        @test !any(ismissing.(df.entrap_pair_id))
        
        # Check that paired proteins have the same pair_id
        prot1_pairs = df[df.protein .== "PROT1", :entrap_pair_id]
        @test length(unique(prot1_pairs)) == 1
        
        prot2_pairs = df[df.protein .== "PROT2", :entrap_pair_id]
        @test length(unique(prot2_pairs)) == 1
        
        # Different proteins should have different pair_ids
        @test prot1_pairs[1] != prot2_pairs[1]
    end
    
    @testset "Multiple entrapment groups" begin
        df = DataFrame(
            protein = ["PROT1", "PROT1", "PROT1", "PROT1"],
            entrap_id = UInt8[0, 1, 2, 0],  # Two originals, two entrapments
            file_name = ["file1", "file1", "file1", "file1"]
        )
        
        assign_protein_entrapment_pairs!(df)
        
        # Each original should be paired with one from each entrapment group
        pair_ids = df.entrap_pair_id
        @test length(unique(pair_ids)) == 2  # Two pairs total
        
        # Check round-robin distribution
        orig_indices = findall(df.entrap_id .== 0)
        @test df.entrap_pair_id[orig_indices[1]] != df.entrap_pair_id[orig_indices[2]]
    end
    
    @testset "Proteins without pairs" begin
        df = DataFrame(
            protein = ["PROT1", "PROT1", "PROT2"],
            entrap_id = UInt8[0, 1, 0],  # PROT2 has no entrapment
            file_name = ["file1", "file1", "file1"]
        )
        
        assign_protein_entrapment_pairs!(df)
        
        # PROT1 should be paired
        prot1_df = df[df.protein .== "PROT1", :]
        @test !any(ismissing.(prot1_df.entrap_pair_id))
        @test length(unique(prot1_df.entrap_pair_id)) == 1
        
        # PROT2 should have missing pair_id
        prot2_df = df[df.protein .== "PROT2", :]
        @test all(ismissing.(prot2_df.entrap_pair_id))
    end
    
    @testset "add_original_target_protein_scores!" begin
        # Test with ms_file_idx column
        df = DataFrame(
            protein = ["PROT1", "PROT1", "PROT2", "PROT2"],
            entrap_id = UInt8[0, 1, 0, 1],
            entrap_pair_id = UInt32[1, 1, 2, 2],
            ms_file_idx = [1, 1, 1, 1],
            pg_score = Float32[100.0, 90.0, 80.0, 70.0]
        )
        
        add_original_target_protein_scores!(df; score_col=:pg_score)
        
        @test hasproperty(df, :pg_score_original_target)
        
        # Original proteins should get their own score
        @test df.pg_score_original_target[1] == 100.0f0
        @test df.pg_score_original_target[3] == 80.0f0
        
        # Entrapment proteins should get their paired original's score
        @test df.pg_score_original_target[2] == 100.0f0
        @test df.pg_score_original_target[4] == 80.0f0
    end
    
    @testset "add_original_target_protein_scores! with file_name" begin
        # Test with file_name column instead of ms_file_idx
        df = DataFrame(
            protein = ["PROT1", "PROT1", "PROT2", "PROT2"],
            entrap_id = UInt8[0, 1, 0, 1],
            entrap_pair_id = UInt32[1, 1, 2, 2],
            file_name = ["file1", "file1", "file2", "file2"],
            pg_score = Float32[100.0, 90.0, 80.0, 70.0]
        )
        
        add_original_target_protein_scores!(df; score_col=:pg_score)
        
        @test hasproperty(df, :pg_score_original_target)
        @test all(df.pg_score_original_target .== [100.0, 100.0, 80.0, 80.0])
    end
    
    @testset "Missing original target" begin
        df = DataFrame(
            protein = ["PROT1", "PROT1"],
            entrap_id = UInt8[1, 1],  # Only entrapments, no original
            entrap_pair_id = UInt32[1, 1],
            ms_file_idx = [1, 1],
            pg_score = Float32[90.0, 85.0]
        )
        
        add_original_target_protein_scores!(df; score_col=:pg_score)
        
        # Should get -1.0 when original is missing
        @test all(df.pg_score_original_target .== -1.0f0)
    end
    
    @testset "create_global_protein_results_df" begin
        # Test with multiple files
        df = DataFrame(
            protein = ["PROT1", "PROT1", "PROT2", "PROT2"],
            ms_file_idx = [1, 2, 1, 2],
            global_pg_score = Float32[100.0, 110.0, 80.0, 75.0],
            entrap_id = UInt8[0, 0, 1, 1]
        )
        
        global_df = create_global_protein_results_df(df; score_col=:global_pg_score)
        
        @test nrow(global_df) == 2  # One row per protein
        @test all(global_df.ms_file_idx .== 0)  # Global indicator
        
        # Should keep best scoring row for each protein
        prot1_row = global_df[global_df.protein .== "PROT1", :]
        @test prot1_row.global_pg_score[1] == 110.0f0
        
        prot2_row = global_df[global_df.protein .== "PROT2", :]
        @test prot2_row.global_pg_score[1] == 80.0f0
    end
    
    @testset "create_global_protein_results_df with file_name" begin
        df = DataFrame(
            protein = ["PROT1", "PROT1"],
            file_name = ["file1", "file2"],
            global_pg_score = Float32[100.0, 110.0],
            entrap_id = UInt8[0, 0]
        )
        
        global_df = create_global_protein_results_df(df; score_col=:global_pg_score)
        
        @test nrow(global_df) == 1
        @test global_df.file_name[1] == "global"
        @test global_df.global_pg_score[1] == 110.0f0
    end
    
    @testset "Multiple score columns" begin
        df = DataFrame(
            protein = ["PROT1", "PROT1"],
            entrap_id = UInt8[0, 1],
            entrap_pair_id = UInt32[1, 1],
            ms_file_idx = [1, 1],
            pg_score = Float32[100.0, 90.0],
            pg_pep = Float32[0.01, 0.02]
        )
        
        add_original_target_protein_scores!(df, [:pg_score, :pg_pep])
        
        @test hasproperty(df, :pg_score_original_target)
        @test hasproperty(df, :pg_pep_original_target)
        @test df.pg_score_original_target == [100.0, 100.0]
        @test df.pg_pep_original_target == [0.01f0, 0.01f0]
    end
end

@testset "Protein EFDR Functions" begin
    # Test add_protein_efdr_columns!
    df = DataFrame(
        protein = ["PROT1", "PROT1", "PROT2", "PROT2"],
        entrap_id = UInt8[0, 1, 0, 1],
        entrap_pair_id = UInt32[1, 1, 2, 2],
        pg_score = Float32[100.0, 90.0, 80.0, 70.0],
        pg_score_original_target = Float32[100.0, 100.0, 80.0, 80.0],
        qval = Float32[0.01, 0.02, 0.01, 0.02]
    )
    
    add_protein_efdr_columns!(df; score_qval_pairs=[(:pg_score, :qval)])
    
    @test hasproperty(df, :pg_score_combined_efdr)
    @test hasproperty(df, :pg_score_paired_efdr)
    @test all(.!ismissing.(df.pg_score_combined_efdr))
    @test all(.!ismissing.(df.pg_score_paired_efdr))
end

@testset "Protein EFDR API Integration" begin
    # Create a mock protein dataset
    mock_data = DataFrame(
        file_name = repeat(["file1.raw", "file2.raw"], inner=4),
        target = fill(true, 8),
        entrap_id = repeat(UInt8[0, 1, 0, 2], 2),
        species = fill("YEAST", 8),
        protein = repeat(["sp|P00330|ADH1_YEAST", "sp|P00330|ADH1_YEAST", 
                         "sp|P00331|ADH2_YEAST", "sp|P00331|ADH2_YEAST"], 2),
        peptides = [missing for _ in 1:8],
        n_peptides = UInt32[10, 10, 8, 8, 10, 10, 8, 8],
        global_qval = Float32[0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.02, 0.02],
        qval = Float32[0.01, 0.015, 0.02, 0.025, 0.01, 0.012, 0.02, 0.022],
        pg_pep = Float32[0.001, 0.001, 0.002, 0.002, 0.001, 0.001, 0.002, 0.002],
        pg_score = Float32[100.5, 95.3, 80.2, 75.1, 98.7, 96.1, 82.3, 77.4],
        global_pg_score = Float32[100.5, 95.3, 80.2, 75.1, 100.5, 95.3, 82.3, 77.4],
        abundance = Float32[1e6, 0.9e6, 0.8e6, 0.75e6, 1.1e6, 0.95e6, 0.85e6, 0.78e6]
    )
    
    # Save to temp file
    temp_path = tempname() * ".arrow"
    Arrow.write(temp_path, mock_data)
    
    # Test that API function can handle protein data
    try
        results = run_protein_efdr_analysis(
            temp_path;
            output_dir=tempname(),
            verbose=false,
            plot_formats=[]  # Skip plots for testing
        )
        
        @test haskey(results, :filtered_data)
        @test haskey(results, :comparison_results)
        @test haskey(results, :output_files)
        
        # Check that EFDR columns were added
        if isa(results.filtered_data, DataFrame)
            @test hasproperty(results.filtered_data, :pg_score_combined_efdr) ||
                  hasproperty(results.filtered_data, :global_pg_score_combined_efdr)
        end
        
    finally
        rm(temp_path, force=true)
    end
end