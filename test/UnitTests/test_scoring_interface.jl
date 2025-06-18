using Test
using Arrow, DataFrames

# Include the necessary files
cd(@__DIR__)
package_root = dirname(dirname(@__DIR__))
include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "ScoringSearch", "scoring_interface.jl"))

@testset "ScoringInterface Tests" begin
    
    @testset "Basic Interface Functions" begin
        # Create temporary test directory
        temp_dir = mktempdir()
        
        # Create test PSM data
        psm_file = joinpath(temp_dir, "test_psms.arrow")
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4, 5],
            prob = Float32[0.9, 0.85, 0.8, 0.75, 0.7],
            inferred_protein_group = ["P1", "P1", "P2", "P2", "P3"],
            target = Bool[true, true, true, false, true],
            entrapment_group_id = UInt8[0, 0, 0, 0, 1],
            use_for_protein_quant = Bool[true, true, true, true, true]
        )
        Arrow.write(psm_file, psm_df)
        
        # Create PSM reference
        psm_ref = PSMFileReference(psm_file)
        
        # Test generate_scoring_summary with empty results
        empty_results = ScoringSearchResultRefs(PairedSearchFiles[])
        summary = generate_scoring_summary(empty_results)
        @test nrow(summary) == 0
        
        # Test with single paired file
        # First need to create a protein file
        protein_file = joinpath(temp_dir, "test_proteins.arrow")
        protein_df = DataFrame(
            protein_name = ["P1", "P2", "P3"],
            target = Bool[true, true, true],
            entrapment_group_id = UInt8[0, 0, 1],
            pg_score = Float32[2.5, 1.8, 1.2],
            n_peptides = Int64[2, 2, 1]
        )
        Arrow.write(protein_file, protein_df)
        
        paired = PairedSearchFiles(psm_file, protein_file, 1)
        results = ScoringSearchResultRefs([paired])
        
        summary = generate_scoring_summary(results)
        @test nrow(summary) == 1
        @test summary.n_psms[1] == 5
        @test summary.n_proteins[1] == 3
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Global Score Calculation" begin
        temp_dir = mktempdir()
        
        # Create multiple protein group files
        pg_refs = ProteinGroupFileReference[]
        
        for i in 1:2
            pg_file = joinpath(temp_dir, "proteins_$i.arrow")
            pg_df = DataFrame(
                protein_name = ["P1", "P2", "P3"],
                target = Bool[true, true, false],
                entrapment_group_id = UInt8[0, 0, 0],
                pg_score = Float32[2.0 + i, 1.5 + i, 1.0 + i]  # Different scores per file
            )
            Arrow.write(pg_file, pg_df)
            push!(pg_refs, ProteinGroupFileReference(pg_file))
        end
        
        # Calculate global scores (should be max across files)
        global_scores = calculate_global_protein_scores(pg_refs)
        
        @test length(global_scores) == 3
        @test global_scores[("P1", true, UInt8(0))] ≈ 4.0f0  # max(3.0, 4.0)
        @test global_scores[("P2", true, UInt8(0))] ≈ 3.5f0  # max(2.5, 3.5)
        @test global_scores[("P3", false, UInt8(0))] ≈ 3.0f0  # max(2.0, 3.0)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Filter Protein Groups" begin
        temp_dir = mktempdir()
        
        # Create protein group file with q-values
        pg_file = joinpath(temp_dir, "proteins.arrow")
        pg_df = DataFrame(
            protein_name = ["P1", "P2", "P3", "P4"],
            target = Bool[true, true, true, false],
            entrapment_group_id = UInt8[0, 0, 0, 0],
            pg_score = Float32[3.0, 2.5, 2.0, 1.5],
            pg_qval = Float32[0.001, 0.005, 0.02, 0.1],
            global_pg_qval = Float32[0.0005, 0.003, 0.015, 0.08]
        )
        Arrow.write(pg_file, pg_df)
        
        pg_ref = ProteinGroupFileReference(pg_file)
        
        # Filter with threshold
        filtered_path = joinpath(temp_dir, "filtered_proteins.arrow")
        filtered_ref = filter_protein_groups_by_qvalue(pg_ref, filtered_path, 0.01f0)
        
        # Check results
        filtered_df = DataFrame(Arrow.Table(file_path(filtered_ref)))
        @test nrow(filtered_df) == 2  # Only P1 and P2 should pass
        @test all(filtered_df.protein_name .∈ Ref(["P1", "P2"]))
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Merge Protein Groups" begin
        temp_dir = mktempdir()
        
        # Create paired files
        paired_files = PairedSearchFiles[]
        
        for i in 1:2
            psm_file = joinpath(temp_dir, "psms_$i.arrow")
            Arrow.write(psm_file, DataFrame(dummy = [1]))
            
            pg_file = joinpath(temp_dir, "proteins_$i.arrow")
            pg_df = DataFrame(
                protein_name = ["P$i", "P$(i+2)"],
                target = Bool[true, true],
                entrapment_group_id = UInt8[0, 0],
                pg_score = Float32[3.0 - i*0.5, 2.0 - i*0.5]
            )
            sort!(pg_df, [:pg_score, :target], rev=[true, true])
            Arrow.write(pg_file, pg_df)
            
            psm_ref = PSMFileReference(psm_file)
            pg_ref = ProteinGroupFileReference(pg_file)
            mark_sorted!(pg_ref, :pg_score, :target)
            
            push!(paired_files, PairedSearchFiles(psm_ref, pg_ref, i))
        end
        
        result_refs = ScoringSearchResultRefs(paired_files)
        
        # Merge protein groups
        merged_path = joinpath(temp_dir, "merged_proteins.arrow")
        merged_ref = merge_protein_groups(result_refs, merged_path)
        
        # Check merged results
        @test exists(merged_ref)
        merged_df = DataFrame(Arrow.Table(file_path(merged_ref)))
        @test nrow(merged_df) == 4  # P1, P2, P3, P4
        @test Set(merged_df.protein_name) == Set(["P1", "P2", "P3", "P4"])
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end