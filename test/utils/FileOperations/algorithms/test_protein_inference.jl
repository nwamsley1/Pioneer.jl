@testset "Protein Inference FileOperations Integration Tests" begin
    
    @testset "File Reference Creation and Validation" begin
        temp_dir = mktempdir()
        
        # Test PSM file reference creation
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            prob = Float32[0.9, 0.8, 0.85],
            target = Bool[true, true, false],
            entrapment_group_id = UInt8[0, 0, 0]
        )
        
        psm_file = joinpath(temp_dir, "test_psms.arrow")
        Arrow.write(psm_file, psm_df)
        psm_ref = PSMFileReference(psm_file)
        
        @test exists(psm_ref)
        @test file_path(psm_ref) == psm_file
        @test row_count(psm_ref) == 3
        
        # Test protein group file reference creation
        pg_df = DataFrame(
            protein_name = ["P1", "P2"],
            target = Bool[true, false],
            entrap_id = UInt8[0, 0],
            n_peptides = Int64[2, 1],
            pg_score = Float32[2.5, 1.8]
        )
        
        pg_file = joinpath(temp_dir, "test_proteins.arrow")
        Arrow.write(pg_file, pg_df)
        pg_ref = ProteinGroupFileReference(pg_file)
        
        @test exists(pg_ref)
        @test file_path(pg_ref) == pg_file
        @test row_count(pg_ref) == 2
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Paired File References" begin
        temp_dir = mktempdir()
        
        # Create matching PSM and protein group files
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            prob = Float32[0.9, 0.8],
            target = Bool[true, true]
        )
        
        pg_df = DataFrame(
            protein_name = ["P1"],
            target = Bool[true],
            n_peptides = Int64[2]
        )
        
        psm_file = joinpath(temp_dir, "paired_psms.arrow")
        pg_file = joinpath(temp_dir, "paired_proteins.arrow")
        
        Arrow.write(psm_file, psm_df)
        Arrow.write(pg_file, pg_df)
        
        # Test paired file creation
        paired_files = PairedSearchFiles(psm_file, pg_file, 1)
        
        @test exists(paired_files.psm_ref)
        @test exists(paired_files.protein_ref)
        @test paired_files.ms_file_idx == 1
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "File Schema Validation" begin
        temp_dir = mktempdir()
        
        # Test PSM schema validation
        complete_psm_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            prob = Float32[0.9, 0.8],
            target = Bool[true, false],
            entrapment_group_id = UInt8[0, 0],
            inferred_protein_group = ["P1", "P1"],
            use_for_protein_quant = Bool[true, true]
        )
        
        psm_file = joinpath(temp_dir, "complete_psms.arrow")
        Arrow.write(psm_file, complete_psm_df)
        psm_ref = PSMFileReference(psm_file)
        
        # Test required columns validation
        required_psm_cols = Set([:precursor_idx, :prob, :target, :entrapment_group_id])
        @test validate_schema(psm_ref, required_psm_cols) === nothing
        
        # Test protein group schema validation  
        complete_pg_df = DataFrame(
            protein_name = ["P1"],
            target = Bool[true],
            entrap_id = UInt8[0],
            n_peptides = Int64[2],
            pg_score = Float32[2.5]
        )
        
        pg_file = joinpath(temp_dir, "complete_proteins.arrow")
        Arrow.write(pg_file, complete_pg_df)
        pg_ref = ProteinGroupFileReference(pg_file)
        
        required_pg_cols = Set([:protein_name, :target, :entrap_id, :n_peptides])
        @test validate_schema(pg_ref, required_pg_cols) === nothing
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Missing File Handling" begin
        temp_dir = mktempdir()
        
        # Test non-existent file references
        missing_psm_file = joinpath(temp_dir, "missing_psms.arrow")
        missing_psm_ref = PSMFileReference(missing_psm_file)
        
        @test !exists(missing_psm_ref)
        @test_throws ErrorException validate_exists(missing_psm_ref)
        
        missing_pg_file = joinpath(temp_dir, "missing_proteins.arrow")
        missing_pg_ref = ProteinGroupFileReference(missing_pg_file)
        
        @test !exists(missing_pg_ref)
        @test_throws ErrorException validate_exists(missing_pg_ref)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "File Operations Error Handling" begin
        temp_dir = mktempdir()
        
        # Test schema validation with missing columns
        incomplete_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            prob = Float32[0.9, 0.8]
            # Missing required columns
        )
        
        incomplete_file = joinpath(temp_dir, "incomplete.arrow")
        Arrow.write(incomplete_file, incomplete_df)
        incomplete_ref = PSMFileReference(incomplete_file)
        
        # Should fail when checking for missing required columns
        required_cols = Set([:precursor_idx, :prob, :target, :entrapment_group_id])
        @test_throws ErrorException validate_schema(incomplete_ref, required_cols)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Multiple File Management" begin
        temp_dir = mktempdir()
        
        # Create multiple file references
        psm_refs = PSMFileReference[]
        pg_refs = ProteinGroupFileReference[]
        
        for i in 1:3
            # Create PSM file
            psm_df = DataFrame(
                precursor_idx = UInt32[i, i+10],
                prob = Float32[0.9, 0.8],
                target = Bool[true, false]
            )
            psm_file = joinpath(temp_dir, "psms_$i.arrow")
            Arrow.write(psm_file, psm_df)
            push!(psm_refs, PSMFileReference(psm_file))
            
            # Create protein group file
            pg_df = DataFrame(
                protein_name = ["P$i"],
                target = Bool[true],
                n_peptides = Int64[1]
            )
            pg_file = joinpath(temp_dir, "proteins_$i.arrow")
            Arrow.write(pg_file, pg_df)
            push!(pg_refs, ProteinGroupFileReference(pg_file))
        end
        
        # Verify all files exist
        @test length(psm_refs) == 3
        @test length(pg_refs) == 3
        
        for psm_ref in psm_refs
            @test exists(psm_ref)
        end
        
        for pg_ref in pg_refs
            @test exists(pg_ref)
        end
        
        # Test file mapping
        file_mapping = Dict{String, String}()
        for (psm_ref, pg_ref) in zip(psm_refs, pg_refs)
            file_mapping[file_path(psm_ref)] = file_path(pg_ref)
        end
        
        @test length(file_mapping) == 3
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end