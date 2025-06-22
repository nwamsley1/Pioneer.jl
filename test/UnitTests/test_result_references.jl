using Test
using Arrow, DataFrames

# Include the necessary files
cd(@__DIR__)
package_root = dirname(dirname(@__DIR__))
include(joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl"))
# SearchResultReferences is still in SearchMethods (not part of FileOperations refactoring)
include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "SearchResultReferences.jl"))

@testset "SearchResultReferences Tests" begin
    
    @testset "ScoringSearchResultRefs" begin
        # Create temporary test files
        temp_dir = mktempdir()
        
        # Create paired files for 2 MS runs
        paired_files = PairedSearchFiles[]
        
        for i in 1:2
            psm_file = joinpath(temp_dir, "psms_$i.arrow")
            protein_file = joinpath(temp_dir, "proteins_$i.arrow")
            
            # Create PSM data with scoring columns
            psm_df = DataFrame(
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
            Arrow.write(psm_file, psm_df)
            
            # Create protein data
            protein_df = DataFrame(
                protein_name = ["P1", "P2"],
                target = Bool[true, true],
                entrapment_group_id = UInt8[0, 1],
                pg_score = Float32[2.1, 1.5],
                global_pg_score = Float32[2.5, 1.8],
                n_peptides = Int64[2, 1],
                pg_qval = Float32[0.01, 0.02],
                global_pg_qval = Float32[0.005, 0.015]
            )
            Arrow.write(protein_file, protein_df)
            
            push!(paired_files, PairedSearchFiles(psm_file, protein_file, i))
        end
        
        # Test creation
        result_refs = ScoringSearchResultRefs(paired_files)
        @test num_files(result_refs) == 2
        @test length(get_psm_refs(result_refs)) == 2
        @test length(get_protein_refs(result_refs)) == 2
        
        # Test validation passes
        @test validate_scoring_output(result_refs)
        
        # Test accessors
        paired = get_paired_file(result_refs, 1)
        @test paired.ms_file_idx == 1
        @test paired.psm_ref.row_count == 3
        @test paired.protein_ref.row_count == 2
        
        # Test with merged scores
        merged_file = joinpath(temp_dir, "merged_scores.arrow")
        merged_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            global_prob = Float32[0.9, 0.8, 0.7],
            prec_prob = Float32[0.85, 0.75, 0.65],
            target = Bool[true, false, true]
        )
        Arrow.write(merged_file, merged_df)
        
        result_refs_with_merged = ScoringSearchResultRefs(
            paired_files,
            PSMFileReference(merged_file)
        )
        @test !isnothing(result_refs_with_merged.merged_scores_ref)
        @test validate_scoring_output(result_refs_with_merged)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "MaxLFQSearchResultRefs" begin
        # Create temporary test files
        temp_dir = mktempdir()
        
        # Create input paired files
        paired_files = PairedSearchFiles[]
        
        psm_file = joinpath(temp_dir, "psms_1.arrow")
        protein_file = joinpath(temp_dir, "proteins_1.arrow")
        
        # Create PSM data ready for MaxLFQ
        psm_df = DataFrame(
            inferred_protein_group = ["P1", "P1", "P2"],
            target = Bool[true, true, true],
            entrapment_group_id = UInt8[0, 0, 0],
            precursor_idx = UInt32[1, 2, 3],
            peak_area = Float32[100.0, 200.0, 150.0],
            use_for_protein_quant = Bool[true, true, true],
            pg_qval = Float32[0.01, 0.01, 0.02],
            global_qval_pg = Float32[0.005, 0.005, 0.015]
        )
        # Sort as required
        sort!(psm_df, [:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx])
        Arrow.write(psm_file, psm_df)
        
        # Mark as sorted
        psm_ref = PSMFileReference(psm_file)
        mark_sorted!(psm_ref, :inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
        
        protein_df = DataFrame(protein_name = ["P1", "P2"])
        Arrow.write(protein_file, protein_df)
        protein_ref = ProteinGroupFileReference(protein_file)
        
        push!(paired_files, PairedSearchFiles(psm_ref, protein_ref, 1))
        
        # Test creation
        maxlfq_refs = MaxLFQSearchResultRefs(paired_files)
        @test length(maxlfq_refs.input_paired_files) == 1
        @test isnothing(maxlfq_refs.precursors_long_ref)
        @test isnothing(maxlfq_refs.proteins_long_ref)
        
        # Test validation
        @test validate_maxlfq_input(get_input_psm_refs(maxlfq_refs))
        
        # Test with output files
        precursors_long = joinpath(temp_dir, "precursors_long.arrow")
        Arrow.write(precursors_long, psm_df)
        
        maxlfq_refs_with_outputs = MaxLFQSearchResultRefs(
            paired_files,
            precursors_long_ref = PSMFileReference(precursors_long)
        )
        @test !isnothing(maxlfq_refs_with_outputs.precursors_long_ref)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Validation Functions" begin
        temp_dir = mktempdir()
        
        # Test missing required columns
        psm_file = joinpath(temp_dir, "incomplete_psms.arrow")
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2],
            prob = Float32[0.9, 0.8]
            # Missing many required columns
        )
        Arrow.write(psm_file, psm_df)
        
        psm_ref = PSMFileReference(psm_file)
        @test_throws ErrorException validate_scoring_psm_schema(psm_ref)
        
        # Test wrong sort order
        psm_file2 = joinpath(temp_dir, "unsorted_psms.arrow")
        psm_df2 = DataFrame(
            inferred_protein_group = ["P2", "P1"],  # Wrong order
            target = Bool[true, true],
            entrapment_group_id = UInt8[0, 0],
            precursor_idx = UInt32[1, 2],
            peak_area = Float32[100.0, 200.0],
            use_for_protein_quant = Bool[true, true],
            pg_qval = Float32[0.01, 0.01],
            global_qval_pg = Float32[0.005, 0.005]
        )
        Arrow.write(psm_file2, psm_df2)
        
        psm_ref2 = PSMFileReference(psm_file2)
        # Not marked as sorted, so should fail validation
        @test_throws ErrorException validate_maxlfq_input([psm_ref2])
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Utility Functions" begin
        temp_dir = mktempdir()
        
        # Create test files
        psm_paths = String[]
        protein_paths = String[]
        
        for i in 1:2
            psm_file = joinpath(temp_dir, "psms_$i.arrow")
            protein_file = joinpath(temp_dir, "proteins_$i.arrow")
            
            Arrow.write(psm_file, DataFrame(col1 = [1, 2]))
            Arrow.write(protein_file, DataFrame(col1 = [3, 4]))
            
            push!(psm_paths, psm_file)
            push!(protein_paths, protein_file)
        end
        
        # Test creation from paths
        result_refs = create_result_refs_from_paths(psm_paths, protein_paths)
        @test num_files(result_refs) == 2
        
        # Test mismatched lengths
        @test_throws ErrorException create_result_refs_from_paths(psm_paths, [protein_paths[1]])
        
        # Test describe functions (just check they don't error)
        @test describe_results(result_refs) === nothing
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end