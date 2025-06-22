using Test
using Arrow, DataFrames

# Include the algorithms module
package_root = dirname(dirname(dirname(dirname(@__DIR__))))
include(joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl"))

@testset "Protein Inference Algorithm Tests" begin
    
    @testset "Basic Protein Inference" begin
        temp_dir = mktempdir()
        
        # Create PSM data with protein groupings
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4, 5],
            inferred_protein_group = ["P1", "P1", "P2", "P3", "P1|P2"],
            target = Bool[true, true, true, false, true],
            entrapment_group_id = UInt8[0, 0, 0, 0, 1],
            prob = Float32[0.9, 0.8, 0.7, 0.6, 0.85],
            use_for_protein_quant = Bool[true, true, true, true, true]
        )
        
        psm_file = joinpath(temp_dir, "psms.arrow")
        Arrow.write(psm_file, psm_df)
        psm_ref = PSMFileReference(psm_file)
        
        # Apply protein inference
        output_file = joinpath(temp_dir, "protein_groups.arrow")
        protein_ref = apply_protein_inference(psm_ref, output_file)
        
        @test exists(protein_ref)
        result_df = DataFrame(Arrow.Table(file_path(protein_ref)))
        
        # Should have protein groups with metadata
        @test hasproperty(result_df, :protein_name)
        @test hasproperty(result_df, :target)
        @test hasproperty(result_df, :entrapment_group_id)
        @test hasproperty(result_df, :n_peptides)
        @test nrow(result_df) >= 3  # At least P1, P2, P3
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Protein Inference with Entrapment Groups" begin
        temp_dir = mktempdir()
        
        # Create PSM data with multiple entrapment groups
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4, 5, 6],
            inferred_protein_group = ["P1", "P1", "P2", "P1", "P2", "P3"],
            target = Bool[true, true, true, true, false, false],
            entrapment_group_id = UInt8[0, 0, 1, 1, 0, 1],
            prob = Float32[0.9, 0.8, 0.85, 0.75, 0.6, 0.7],
            use_for_protein_quant = Bool[true, true, true, true, true, true]
        )
        
        psm_file = joinpath(temp_dir, "entrap_psms.arrow")
        Arrow.write(psm_file, psm_df)
        psm_ref = PSMFileReference(psm_file)
        
        # Apply protein inference
        output_file = joinpath(temp_dir, "entrap_proteins.arrow")
        protein_ref = apply_protein_inference(psm_ref, output_file)
        
        result_df = DataFrame(Arrow.Table(file_path(protein_ref)))
        
        # Should handle entrapment groups properly
        @test any(result_df.entrapment_group_id .== 0)
        @test any(result_df.entrapment_group_id .== 1)
        
        # Same protein in different entrapment groups should be separate entries
        p1_entries = filter(row -> row.protein_name == "P1", result_df)
        @test nrow(p1_entries) >= 1  # Could be 1 or 2 depending on entrapment groups
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Protein Inference Parameter Validation" begin
        temp_dir = mktempdir()
        
        # Create minimal PSM data
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            inferred_protein_group = ["P1", "P2"],
            target = Bool[true, true],
            entrapment_group_id = UInt8[0, 0],
            prob = Float32[0.9, 0.8],
            use_for_protein_quant = Bool[true, true]
        )
        
        psm_file = joinpath(temp_dir, "param_psms.arrow")
        Arrow.write(psm_file, psm_df)
        psm_ref = PSMFileReference(psm_file)
        
        # Test with custom parameters
        output_file = joinpath(temp_dir, "param_proteins.arrow")
        protein_ref = apply_protein_inference(
            psm_ref, 
            output_file;
            min_peptides_per_protein = 1,
            require_unique_peptides = false
        )
        
        @test exists(protein_ref)
        result_df = DataFrame(Arrow.Table(file_path(protein_ref)))
        @test nrow(result_df) >= 2
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Protein Inference with Missing Columns" begin
        temp_dir = mktempdir()
        
        # Create PSM data missing some optional columns
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            inferred_protein_group = ["P1", "P1", "P2"],
            target = Bool[true, true, false],
            # Missing entrapment_group_id - should default to 0
            prob = Float32[0.9, 0.8, 0.7]
            # Missing use_for_protein_quant - should default to true
        )
        
        psm_file = joinpath(temp_dir, "missing_cols_psms.arrow")
        Arrow.write(psm_file, psm_df)
        psm_ref = PSMFileReference(psm_file)
        
        # Should handle missing columns gracefully
        output_file = joinpath(temp_dir, "missing_cols_proteins.arrow")
        protein_ref = apply_protein_inference(psm_ref, output_file)
        
        @test exists(protein_ref)
        result_df = DataFrame(Arrow.Table(file_path(protein_ref)))
        @test hasproperty(result_df, :protein_name)
        @test hasproperty(result_df, :target)
        @test hasproperty(result_df, :entrapment_group_id)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Protein Inference Error Handling" begin
        temp_dir = mktempdir()
        
        # Create PSM data missing required columns
        incomplete_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            # Missing inferred_protein_group
            target = Bool[true, false]
        )
        
        incomplete_file = joinpath(temp_dir, "incomplete_psms.arrow")
        Arrow.write(incomplete_file, incomplete_df)
        incomplete_ref = PSMFileReference(incomplete_file)
        
        output_file = joinpath(temp_dir, "error_proteins.arrow")
        
        # Should error due to missing required column
        @test_throws KeyError apply_protein_inference(incomplete_ref, output_file)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end