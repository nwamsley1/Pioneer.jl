using Test
using Arrow, DataFrames

# Include the I/O module
package_root = dirname(dirname(dirname(dirname(@__DIR__))))
include(joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl"))

@testset "Validation Functions Tests" begin
    
    @testset "MaxLFQ Input Validation" begin
        temp_dir = mktempdir()
        
        # Create valid MaxLFQ input data
        valid_maxlfq_df = DataFrame(
            inferred_protein_group = ["P1", "P1", "P2", "P2"],
            target = Bool[true, true, false, true],
            entrapment_group_id = UInt8[0, 0, 0, 1],
            precursor_idx = UInt32[1, 2, 3, 4],
            pg_qval = Float32[0.01, 0.02, 0.05, 0.03],
            global_qval_pg = Float32[0.01, 0.02, 0.05, 0.03],
            use_for_protein_quant = Bool[true, true, true, false],
            peak_area = Float32[1000.0, 2000.0, 1500.0, 800.0]
        )
        valid_file = joinpath(temp_dir, "valid_maxlfq.arrow")
        Arrow.write(valid_file, valid_maxlfq_df)
        valid_ref = PSMFileReference(valid_file)
        
        # Test valid input passes validation
        @test validate_maxlfq_input(valid_ref) === nothing
        
        # Test with vector input
        @test validate_maxlfq_input([valid_ref]) === nothing
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "MaxLFQ Missing Columns Validation" begin
        temp_dir = mktempdir()
        
        # Test missing required columns one by one
        required_columns = [
            :inferred_protein_group,
            :target,
            :entrapment_group_id,
            :precursor_idx,
            :pg_qval,
            :global_qval_pg,
            :use_for_protein_quant,
            :peak_area
        ]
        
        base_df = DataFrame(
            inferred_protein_group = ["P1", "P2"],
            target = Bool[true, false],
            entrapment_group_id = UInt8[0, 0],
            precursor_idx = UInt32[1, 2],
            pg_qval = Float32[0.01, 0.02],
            global_qval_pg = Float32[0.01, 0.02],
            use_for_protein_quant = Bool[true, true],
            peak_area = Float32[1000.0, 2000.0]
        )
        
        for missing_col in required_columns
            # Create DataFrame missing one required column
            test_df = select(base_df, Not(missing_col))
            test_file = joinpath(temp_dir, "missing_$(missing_col).arrow")
            Arrow.write(test_file, test_df)
            test_ref = PSMFileReference(test_file)
            
            @test_throws ErrorException validate_maxlfq_input(test_ref)
            rm(test_file)
        end
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "MaxLFQ Parameter Validation" begin
        # Test valid parameters
        valid_params = Dict(
            :q_value_threshold => 0.05,
            :batch_size => 100000,
            :min_peptides => 1,
            :use_global_scores => true
        )
        @test validate_maxlfq_parameters(valid_params) === nothing
        
        # Test invalid q_value_threshold
        invalid_qval_params = copy(valid_params)
        invalid_qval_params[:q_value_threshold] = 1.5  # > 1.0
        @test_throws ErrorException validate_maxlfq_parameters(invalid_qval_params)
        
        invalid_qval_params[:q_value_threshold] = -0.1  # < 0
        @test_throws ErrorException validate_maxlfq_parameters(invalid_qval_params)
        
        # Test invalid batch_size
        invalid_batch_params = copy(valid_params)
        invalid_batch_params[:batch_size] = -100  # negative
        @test_throws ErrorException validate_maxlfq_parameters(invalid_batch_params)
        
        invalid_batch_params[:batch_size] = 0  # zero
        @test_throws ErrorException validate_maxlfq_parameters(invalid_batch_params)
        
        # Test invalid min_peptides
        invalid_peptides_params = copy(valid_params)
        invalid_peptides_params[:min_peptides] = 0  # < 1
        @test_throws ErrorException validate_maxlfq_parameters(invalid_peptides_params)
        
        # Test missing required parameters
        incomplete_params = Dict(:q_value_threshold => 0.05)
        @test_throws ErrorException validate_maxlfq_parameters(incomplete_params)
    end
    
    @testset "PSM Schema Validation" begin
        temp_dir = mktempdir()
        
        # Create PSM data with complete schema
        complete_psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            prob = Float32[0.9, 0.8, 0.7],
            prec_prob = Float32[0.85, 0.75, 0.65],
            global_prob = Float32[0.88, 0.78, 0.68],
            target = Bool[true, false, true],
            entrapment_group_id = UInt8[0, 0, 1],
            inferred_protein_group = ["P1", "P1", "P2"],
            pg_score = Float32[2.1, 2.1, 1.5],
            global_pg_score = Float32[2.5, 2.5, 1.8],
            pg_qval = Float32[0.01, 0.01, 0.02],
            global_qval_pg = Float32[0.005, 0.005, 0.015],
            use_for_protein_quant = Bool[true, true, true]
        )
        
        complete_file = joinpath(temp_dir, "complete_psm.arrow")
        Arrow.write(complete_file, complete_psm_df)
        complete_ref = PSMFileReference(complete_file)
        
        # Should pass validation
        @test validate_scoring_psm_schema(complete_ref) === nothing
        
        # Test with minimal required columns
        minimal_df = select(complete_psm_df, [
            :precursor_idx, :prob, :target, :entrapment_group_id,
            :inferred_protein_group, :use_for_protein_quant
        ])
        minimal_file = joinpath(temp_dir, "minimal_psm.arrow")
        Arrow.write(minimal_file, minimal_df)
        minimal_ref = PSMFileReference(minimal_file)
        
        @test validate_scoring_psm_schema(minimal_ref) === nothing
        
        # Test missing required column
        incomplete_df = select(complete_psm_df, Not(:precursor_idx))
        incomplete_file = joinpath(temp_dir, "incomplete_psm.arrow")
        Arrow.write(incomplete_file, incomplete_df)
        incomplete_ref = PSMFileReference(incomplete_file)
        
        @test_throws ErrorException validate_scoring_psm_schema(incomplete_ref)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Protein Group Schema Validation" begin
        temp_dir = mktempdir()
        
        # Create protein group data with complete schema
        complete_protein_df = DataFrame(
            protein_name = ["P1", "P2", "P3"],
            target = Bool[true, true, false],
            entrapment_group_id = UInt8[0, 0, 0],
            pg_score = Float32[2.5, 1.8, 1.2],
            global_pg_score = Float32[3.0, 2.1, 1.5],
            n_peptides = Int64[3, 2, 1],
            pg_qval = Float32[0.01, 0.02, 0.05],
            global_pg_qval = Float32[0.005, 0.015, 0.04]
        )
        
        complete_file = joinpath(temp_dir, "complete_protein.arrow")
        Arrow.write(complete_file, complete_protein_df)
        complete_ref = ProteinGroupFileReference(complete_file)
        
        # Should pass validation
        @test validate_protein_group_schema(complete_ref) === nothing
        
        # Test with minimal required columns
        minimal_df = select(complete_protein_df, [
            :protein_name, :target, :entrapment_group_id, :pg_score
        ])
        minimal_file = joinpath(temp_dir, "minimal_protein.arrow")
        Arrow.write(minimal_file, minimal_df)
        minimal_ref = ProteinGroupFileReference(minimal_file)
        
        @test validate_protein_group_schema(minimal_ref) === nothing
        
        # Test missing required column
        incomplete_df = select(complete_protein_df, Not(:protein_name))
        incomplete_file = joinpath(temp_dir, "incomplete_protein.arrow")
        Arrow.write(incomplete_file, incomplete_df)
        incomplete_ref = ProteinGroupFileReference(incomplete_file)
        
        @test_throws ErrorException validate_protein_group_schema(incomplete_ref)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Sort Order Validation" begin
        temp_dir = mktempdir()
        
        # Create properly sorted data
        sorted_df = DataFrame(
            protein_name = ["P1", "P1", "P2", "P2"],
            target = Bool[true, false, true, false],
            score = Float32[0.9, 0.8, 0.7, 0.6]
        )
        # Ensure it's actually sorted
        sort!(sorted_df, [:protein_name, :target, :score], rev=[false, true, true])
        
        sorted_file = joinpath(temp_dir, "sorted.arrow")
        Arrow.write(sorted_file, sorted_df)
        sorted_ref = PSMFileReference(sorted_file)
        mark_sorted!(sorted_ref, :protein_name, :target, :score)
        
        # Should pass validation
        required_sort = [:protein_name, :target, :score]
        @test validate_sort_order(sorted_ref, required_sort) === nothing
        
        # Test with incorrect sort marking
        incorrect_ref = PSMFileReference(sorted_file)
        mark_sorted!(incorrect_ref, :score)  # Wrong sort order
        @test_throws ErrorException validate_sort_order(incorrect_ref, required_sort)
        
        # Test with unsorted file
        unsorted_ref = PSMFileReference(sorted_file)
        # Don't mark as sorted
        @test_throws ErrorException validate_sort_order(unsorted_ref, required_sort)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "File Existence Validation" begin
        temp_dir = mktempdir()
        
        # Create existing file
        existing_df = DataFrame(a = [1, 2, 3])
        existing_file = joinpath(temp_dir, "exists.arrow")
        Arrow.write(existing_file, existing_df)
        existing_ref = PSMFileReference(existing_file)
        
        # Should pass validation
        @test validate_file_exists(existing_ref) === nothing
        
        # Test with non-existent file
        missing_ref = PSMFileReference(joinpath(temp_dir, "missing.arrow"))
        @test_throws ErrorException validate_file_exists(missing_ref)
        
        # Test with vector of references
        @test validate_file_exists([existing_ref]) === nothing
        @test_throws ErrorException validate_file_exists([existing_ref, missing_ref])
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Combined Validation Scenarios" begin
        temp_dir = mktempdir()
        
        # Create comprehensive test data
        test_df = DataFrame(
            inferred_protein_group = ["P1", "P1", "P2"],
            target = Bool[true, true, false],
            entrapment_group_id = UInt8[0, 0, 0],
            precursor_idx = UInt32[1, 2, 3],
            pg_qval = Float32[0.01, 0.02, 0.05],
            global_qval_pg = Float32[0.01, 0.02, 0.05],
            use_for_protein_quant = Bool[true, true, true],
            peak_area = Float32[1000.0, 2000.0, 1500.0]
        )
        
        # Sort for MaxLFQ requirements
        sort!(test_df, [:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx])
        
        test_file = joinpath(temp_dir, "comprehensive.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        mark_sorted!(test_ref, :inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
        
        # Should pass all validations
        @test validate_file_exists(test_ref) === nothing
        @test validate_maxlfq_input(test_ref) === nothing
        @test validate_sort_order(test_ref, [:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx]) === nothing
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end