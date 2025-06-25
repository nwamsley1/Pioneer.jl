@testset "PSM Protein Group Score Update Tests" begin
    
    @testset "Basic Protein Group Score Updates" begin
        temp_dir = mktempdir()
        
        # Create PSM data with protein group information
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4],
            prob = Float32[0.8, 0.7, 0.9, 0.6],
            target = Bool[true, true, false, true],
            entrapment_group_id = UInt8[0, 0, 0, 1],
            inferred_protein_group = ["P1", "P1", "P2", "P2"]
        )
        
        # Create protein group scores
        protein_scores_df = DataFrame(
            protein_name = ["P1", "P2"],
            target = Bool[true, false],
            entrapment_group_id = UInt8[0, 0],
            pg_score = Float32[2.5, 1.8],
            global_pg_score = Float32[2.8, 2.1],
            pg_qval = Float32[0.01, 0.05],
            global_pg_qval = Float32[0.008, 0.04]
        )
        
        psm_file = joinpath(temp_dir, "psms.arrow")
        scores_file = joinpath(temp_dir, "protein_scores.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(scores_file, protein_scores_df)
        
        psm_ref = PSMFileReference(psm_file)
        scores_ref = ProteinGroupFileReference(scores_file)
        
        # Apply protein group score updates
        output_file = joinpath(temp_dir, "updated_psms.arrow")
        updated_ref = update_psms_with_scores(psm_ref, scores_ref, output_file)
        
        @test exists(updated_ref)
        result_df = DataFrame(Arrow.Table(file_path(updated_ref)))
        
        # Should have original columns plus protein score columns
        @test hasproperty(result_df, :precursor_idx)
        @test hasproperty(result_df, :prob)
        @test hasproperty(result_df, :pg_score)
        @test hasproperty(result_df, :global_pg_score)
        @test hasproperty(result_df, :pg_qval)
        @test hasproperty(result_df, :global_qval_pg)
        @test nrow(result_df) == 4
        
        # Verify protein scores were applied correctly
        # PSMs 1,2 belong to P1, PSMs 3,4 belong to P2
        @test result_df.pg_score[1] == Float32(2.5)  # P1
        @test result_df.pg_score[2] == Float32(2.5)  # P1
        @test result_df.pg_score[3] == Float32(1.8)  # P2
        @test result_df.pg_score[4] == Float32(1.8)  # P2
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end