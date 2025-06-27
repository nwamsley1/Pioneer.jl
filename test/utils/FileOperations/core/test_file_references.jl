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


@testset "FileReference Core Tests" begin
    
    @testset "PSMFileReference Basic Tests" begin
        temp_dir = mktempdir()
        
        # Create PSM test data
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4, 5],
            prob = Float32[0.9, 0.85, 0.8, 0.75, 0.7],
            inferred_protein_group = ["P1", "P1", "P2", "P2", "P3"],
            target = Bool[true, true, true, false, true],
            entrapment_group_id = UInt8[0, 0, 0, 0, 1]
        )
        
        psm_file = joinpath(temp_dir, "test_psms.arrow")
        Arrow.write(psm_file, psm_df)
        
        # Test PSMFileReference creation
        psm_ref = PSMFileReference(psm_file)
        @test psm_ref.file_exists
        @test psm_ref.row_count == 5
        @test psm_ref.sorted_by == ()
        
        # Test accessor functions
        @test file_path(psm_ref) == psm_file
        @test exists(psm_ref) == true
        @test row_count(psm_ref) == 5
        @test sorted_by(psm_ref) == ()
        
        # Test schema detection
        expected_cols = Set([:precursor_idx, :prob, :inferred_protein_group, :target, :entrapment_group_id])
        @test psm_ref.schema.columns == expected_cols
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "ProteinGroupFileReference Basic Tests" begin
        temp_dir = mktempdir()
        
        # Create protein group test data
        protein_df = DataFrame(
            protein_name = ["P1", "P2", "P3"],
            target = Bool[true, true, false],
            entrapment_group_id = UInt8[0, 0, 0],
            pg_score = Float32[2.5, 1.8, 1.2],
            n_peptides = Int64[3, 2, 1]
        )
        
        protein_file = joinpath(temp_dir, "test_proteins.arrow")
        Arrow.write(protein_file, protein_df)
        
        # Test ProteinGroupFileReference creation
        protein_ref = ProteinGroupFileReference(protein_file)
        @test protein_ref.file_exists
        @test protein_ref.row_count == 3
        @test protein_ref.sorted_by == ()
        
        # Test accessor functions
        @test file_path(protein_ref) == protein_file
        @test exists(protein_ref) == true
        @test row_count(protein_ref) == 3
        
        # Test schema detection
        expected_cols = Set([:protein_name, :target, :entrapment_group_id, :pg_score, :n_peptides])
        @test protein_ref.schema.columns == expected_cols
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "ProteinQuantFileReference Specialized Tests" begin
        temp_dir = mktempdir()
        
        # Create protein quantification test data
        protein_quant_df = DataFrame(
            file_name = ["file1.raw", "file2.raw", "file1.raw", "file2.raw", "file3.raw"],
            protein = ["P1", "P1", "P2", "P2", "P3"],
            target = Bool[true, true, true, false, true],
            entrap_id = UInt8[0, 0, 0, 0, 1],
            abundance = Float32[5000.0, 4500.0, 3000.0, 2800.0, 2000.0],
            n_peptides = Int64[3, 3, 2, 2, 1],
            experiments = UInt32[1, 2, 1, 2, 3]
        )
        
        quant_file = joinpath(temp_dir, "protein_quant.arrow")
        Arrow.write(quant_file, protein_quant_df)
        
        # Test ProteinQuantFileReference creation
        quant_ref = ProteinQuantFileReference(quant_file)
        @test quant_ref.file_exists
        @test quant_ref.row_count == 5
        
        # Test specialized accessors
        @test n_protein_groups(quant_ref) == 3  # P1, P2, P3
        @test n_experiments(quant_ref) == 3     # experiments 1, 2, 3
        
        # Test with file_name column instead of experiments
        protein_quant_df2 = select(protein_quant_df, Not(:experiments))
        quant_file2 = joinpath(temp_dir, "protein_quant2.arrow")
        Arrow.write(quant_file2, protein_quant_df2)
        
        quant_ref2 = ProteinQuantFileReference(quant_file2)
        @test n_experiments(quant_ref2) == 3  # Should count unique file_name values
        
        # Test with missing protein column
        no_protein_df = select(protein_quant_df, Not(:protein))
        no_protein_file = joinpath(temp_dir, "no_protein.arrow")
        Arrow.write(no_protein_file, no_protein_df)
        
        no_protein_ref = ProteinQuantFileReference(no_protein_file)
        @test n_protein_groups(no_protein_ref) == 0
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "FileReference Edge Cases" begin
        temp_dir = mktempdir()
        
        # Test with empty file
        empty_df = DataFrame()
        empty_file = joinpath(temp_dir, "empty.arrow")
        Arrow.write(empty_file, empty_df)
        
        empty_ref = PSMFileReference(empty_file)
        @test empty_ref.file_exists
        @test empty_ref.row_count == 0
        @test isempty(empty_ref.schema.columns)
        
        # Test with single-row DataFrame
        single_df = DataFrame(only_col = [42])
        single_file = joinpath(temp_dir, "single.arrow")
        Arrow.write(single_file, single_df)
        
        single_ref = PSMFileReference(single_file)
        @test single_ref.row_count == 1
        @test has_column(single_ref.schema, :only_col)
        
        # Test with DataFrame containing missing values
        missing_df = DataFrame(
            id = [1, 2, 3],
            value = [1.0, missing, 3.0],
            text = ["a", "b", missing]
        )
        missing_file = joinpath(temp_dir, "missing.arrow")
        Arrow.write(missing_file, missing_df)
        
        missing_ref = PSMFileReference(missing_file)
        @test missing_ref.row_count == 3
        @test has_column(missing_ref.schema, :value)
        @test has_column(missing_ref.schema, :text)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "PairedSearchFiles Tests" begin
        temp_dir = mktempdir()
        
        # Create matching PSM and protein files
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            inferred_protein_group = ["P1", "P1", "P2"],
            prob = Float32[0.9, 0.8, 0.7]
        )
        protein_df = DataFrame(
            protein_name = ["P1", "P2"],
            n_peptides = Int64[2, 1],
            score = Float32[1.8, 0.7]
        )
        
        psm_file = joinpath(temp_dir, "paired_psms.arrow")
        protein_file = joinpath(temp_dir, "paired_proteins.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(protein_file, protein_df)
        
        # Test successful pairing with file paths
        paired = PairedSearchFiles(psm_file, protein_file, 42)
        @test paired.ms_file_idx == 42
        @test paired.psm_ref.file_exists
        @test paired.protein_ref.file_exists
        @test paired.psm_ref.row_count == 3
        @test paired.protein_ref.row_count == 2
        
        # Test pairing with reference objects
        psm_ref = PSMFileReference(psm_file)
        protein_ref = ProteinGroupFileReference(protein_file)
        paired2 = PairedSearchFiles(psm_ref, protein_ref, 24)
        @test paired2.ms_file_idx == 24
        
        # Test validation failure - remove one file
        rm(protein_file)
        @test_throws ErrorException PairedSearchFiles(psm_file, protein_file, 1)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end