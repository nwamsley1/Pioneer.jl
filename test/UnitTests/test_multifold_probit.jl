# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

using Test
using DataFrames
using Arrow
using Pioneer
using Random

# Mock LibraryPrecursors for testing
struct MockLibraryPrecursors <: Pioneer.LibraryPrecursors
    data::Nothing
    n::Int64
    accession_numbers_to_pid::Dictionary{String, UInt32}
    pid_to_cv_fold::Vector{UInt8}
end

# Implement required interface
Pioneer.getCvFold(lp::MockLibraryPrecursors, precursor_idx::Integer) = lp.pid_to_cv_fold[precursor_idx]
Base.length(lp::MockLibraryPrecursors) = lp.n

"""
Create mock LibraryPrecursors with specified number of CV folds
"""
function create_mock_precursors(n_precursors::Int, n_folds::Int; seed=42)
    Random.seed!(seed)
    
    # Create CV fold assignments (0 to n_folds-1)
    pid_to_cv_fold = UInt8.(rand(0:(n_folds-1), n_precursors))
    
    # Create mock accession mapping
    accession_numbers_to_pid = Dictionary{String, UInt32}()
    
    return MockLibraryPrecursors(
        nothing,
        n_precursors,
        accession_numbers_to_pid,
        pid_to_cv_fold
    )
end

"""
Create mock PSM data for testing
"""
function create_mock_psm_data(n_psms::Int, n_precursors::Int, protein_names::Vector{String}; seed=43)
    Random.seed!(seed)
    
    # Create PSM DataFrame
    psms = DataFrame(
        precursor_idx = rand(1:n_precursors, n_psms),
        prob = rand(Float32, n_psms),
        target = rand(Bool, n_psms),
        use_for_protein_quant = fill(true, n_psms),
        inferred_protein_group = rand(protein_names, n_psms)
    )
    
    return psms
end

"""
Create mock protein group data
"""
function create_mock_protein_groups(protein_names::Vector{String}; seed=44)
    Random.seed!(seed)
    n_proteins = length(protein_names)
    
    # Create protein groups with features
    protein_groups = DataFrame(
        protein_name = protein_names,
        target = rand(Bool, n_proteins),
        entrap_id = zeros(UInt8, n_proteins),
        pg_score = rand(Float32, n_proteins),
        n_peptides = rand(UInt16(1):UInt16(10), n_proteins),
        n_possible_peptides = rand(UInt16(5):UInt16(20), n_proteins),
        peptide_coverage = rand(Float32, n_proteins),
        total_peptide_length = rand(UInt16(10):UInt16(100), n_proteins),
        log_n_possible_peptides = log.(Float32.(rand(5:20, n_proteins))),
        log_binom_coeff = rand(Float32, n_proteins)
    )
    
    return protein_groups
end

"""
Create mock file references for testing
"""
function create_mock_file_refs(temp_dir::String, protein_groups::DataFrame, psms::DataFrame)
    # Write protein group file
    pg_path = joinpath(temp_dir, "passing_proteins_001.arrow")
    Arrow.write(pg_path, protein_groups)
    pg_ref = Pioneer.ProteinGroupFileReference(pg_path)
    
    # Write PSM file
    psm_path = joinpath(temp_dir, "scored_PSMs_001.arrow")
    Arrow.write(psm_path, psms)
    
    return [pg_ref]
end

@testset "Multi-fold Probit Analysis Tests" begin
    
    @testset "CV Fold Detection" begin
        # Test 2-fold
        precursors_2fold = create_mock_precursors(100, 2)
        folds = Pioneer.detect_unique_cv_folds(precursors_2fold)
        @test length(folds) == 2
        @test folds == [0, 1]
        
        # Test 5-fold
        precursors_5fold = create_mock_precursors(100, 5)
        folds = Pioneer.detect_unique_cv_folds(precursors_5fold)
        @test length(folds) == 5
        @test folds == [0, 1, 2, 3, 4]
        
        # Test 3-fold
        precursors_3fold = create_mock_precursors(100, 3)
        folds = Pioneer.detect_unique_cv_folds(precursors_3fold)
        @test length(folds) == 3
        @test folds == [0, 1, 2]
    end
    
    @testset "PSM Path Mapping" begin
        # Create mock reference
        pg_path = "/path/to/results/passing_proteins_001.arrow"
        pg_ref = Pioneer.ProteinGroupFileReference(pg_path)
        
        # Test path mapping
        psm_path = Pioneer.get_corresponding_psm_path(pg_ref)
        @test psm_path == "/path/to/results/scored_PSMs_001.arrow"
    end
    
    @testset "Protein Group CV Fold Assignment" begin
        # Create temporary directory
        temp_dir = mktempdir()
        
        try
            # Create mock data
            protein_names = ["Protein_A", "Protein_B", "Protein_C", "Protein_D"]
            n_precursors = 20
            precursors = create_mock_precursors(n_precursors, 3)
            
            # Create PSMs with known associations
            psms = DataFrame(
                precursor_idx = [1, 2, 3, 4, 5, 6, 7, 8],
                prob = [0.9, 0.8, 0.7, 0.6, 0.95, 0.85, 0.75, 0.65],
                target = fill(true, 8),
                use_for_protein_quant = fill(true, 8),
                inferred_protein_group = ["Protein_A", "Protein_A", "Protein_B", "Protein_B",
                                         "Protein_C", "Protein_C", "Protein_D", "Protein_D"]
            )
            
            # Create protein groups
            protein_groups = create_mock_protein_groups(protein_names)
            
            # Create file references
            pg_refs = create_mock_file_refs(temp_dir, protein_groups, psms)
            
            # Test CV fold assignment
            Pioneer.assign_protein_group_cv_folds!(protein_groups, pg_refs, precursors)
            
            # Check that cv_fold column was added
            @test hasproperty(protein_groups, :cv_fold)
            @test length(protein_groups.cv_fold) == length(protein_names)
            
            # Check that each protein got the cv_fold of its highest scoring peptide
            # Protein_A should get cv_fold of precursor 1 (highest score 0.9)
            # Protein_C should get cv_fold of precursor 5 (highest score 0.95)
            for (i, protein_name) in enumerate(protein_names)
                protein_psms = filter(row -> row.inferred_protein_group == protein_name, psms)
                if nrow(protein_psms) > 0
                    best_idx = argmax(protein_psms.prob)
                    best_precursor_idx = protein_psms.precursor_idx[best_idx]
                    expected_fold = Pioneer.getCvFold(precursors, best_precursor_idx)
                    @test protein_groups.cv_fold[i] == expected_fold
                end
            end
            
        finally
            rm(temp_dir, recursive=true)
        end
    end
    
    @testset "Multi-fold Probit Regression" begin
        # Create temporary directory
        temp_dir = mktempdir()
        qc_folder = joinpath(temp_dir, "qc")
        mkpath(qc_folder)
        
        try
            # Create larger mock dataset
            protein_names = ["Protein_$i" for i in 1:50]
            n_precursors = 100
            n_folds = 3
            precursors = create_mock_precursors(n_precursors, n_folds)
            
            # Create PSMs
            psms = create_mock_psm_data(200, n_precursors, protein_names)
            
            # Create protein groups with more realistic distribution
            protein_groups = create_mock_protein_groups(protein_names)
            # Ensure we have both targets and decoys
            protein_groups.target[1:25] .= true
            protein_groups.target[26:50] .= false
            
            # Create file references
            pg_refs = create_mock_file_refs(temp_dir, protein_groups, psms)
            
            # Run multi-fold probit analysis
            Pioneer.perform_probit_analysis_multifold(
                protein_groups,
                qc_folder,
                pg_refs,
                precursors;
                show_improvement = true
            )
            
            # Check that original scores were saved
            @test hasproperty(protein_groups, :old_pg_score)
            
            # Check that cv_fold column was removed after analysis
            @test !hasproperty(protein_groups, :cv_fold)
            
            # Check that scores were updated (at least some should be different)
            @test any(protein_groups.pg_score .!= protein_groups.old_pg_score)
            
            # Check that updated file exists
            pg_file = Pioneer.file_path(pg_refs[1])
            @test isfile(pg_file)
            
            # Load and check updated file
            updated_pg = DataFrame(Arrow.Table(pg_file))
            @test !hasproperty(updated_pg, :cv_fold)  # Temporary column should be removed
            @test hasproperty(updated_pg, :pg_score)  # Scores should be present
            
        finally
            rm(temp_dir, recursive=true)
        end
    end
    
    @testset "Edge Cases" begin
        # Test with single fold (no cross-validation possible)
        precursors_1fold = MockLibraryPrecursors(
            nothing,
            10,
            Dictionary{String, UInt32}(),
            fill(UInt8(0), 10)  # All precursors in fold 0
        )
        
        folds = Pioneer.detect_unique_cv_folds(precursors_1fold)
        @test length(folds) == 1
        @test folds == [0]
        
        # Test empty protein groups
        empty_protein_groups = DataFrame(
            protein_name = String[],
            target = Bool[],
            pg_score = Float32[],
            peptide_coverage = Float32[],
            n_possible_peptides = UInt16[]
        )
        
        # Should handle gracefully without errors
        temp_dir = mktempdir()
        try
            Pioneer.assign_protein_group_cv_folds!(
                empty_protein_groups,
                Pioneer.ProteinGroupFileReference[],
                precursors_1fold
            )
            @test nrow(empty_protein_groups) == 0
        finally
            rm(temp_dir, recursive=true)
        end
    end
end

