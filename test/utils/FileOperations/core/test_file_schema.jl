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


@testset "FileSchema Tests" begin
    
    @testset "FileSchema Creation and Basic Operations" begin
        # Test basic schema creation
        schema = FileSchema([:col1, :col2, :col3])
        @test schema.columns == Set([:col1, :col2, :col3])
        
        # Test with empty schema
        empty_schema = FileSchema(Symbol[])
        @test isempty(empty_schema.columns)
        
        # Test duplicate columns are handled by Set
        dup_schema = FileSchema([:col1, :col2, :col1])
        @test dup_schema.columns == Set([:col1, :col2])
    end
    
    @testset "FileSchema Column Operations" begin
        schema = FileSchema([:int_col, :float_col, :string_col, :bool_col])
        
        # Test has_column
        @test has_column(schema, :int_col)
        @test has_column(schema, :float_col)
        @test !has_column(schema, :nonexistent_col)
        
        # Test column validation
        required_subset = Set([:int_col, :float_col])
        @test validate_required_columns(schema, required_subset) === nothing
        
        # Test validation failure
        missing_subset = Set([:int_col, :fake_col1, :fake_col2])
        @test_throws ErrorException validate_required_columns(schema, missing_subset)
    end
    
    @testset "get_column_or_default Tests" begin
        test_df = DataFrame(
            int_col = [1, 2, 3],
            float_col = [1.5, 2.5, 3.5],
            string_col = ["a", "b", "c"]
        )
        
        # Test with complete schema
        test_schema = FileSchema([:int_col, :float_col, :string_col])
        @test get_column_or_default(test_df, test_schema, :int_col, 999) == [1, 2, 3]
        
        # Test with limited schema - column not in schema should return default
        limited_schema = FileSchema([:int_col, :float_col])
        @test get_column_or_default(test_df, limited_schema, :string_col, "missing") == ["missing", "missing", "missing"]
        
        # Test with different default types
        @test get_column_or_default(test_df, limited_schema, :nonexistent_col, true) == [true, true, true]
        @test get_column_or_default(test_df, limited_schema, :nonexistent_col, 42) == [42, 42, 42]
    end
    
    @testset "FileSchema with Complex Column Types" begin
        # Test schema creation with various column types
        complex_schema = FileSchema([:int_col, :float_col, :string_col, :bool_col, :missing_col])
        
        @test complex_schema.columns == Set([:int_col, :float_col, :string_col, :bool_col, :missing_col])
        @test has_column(complex_schema, :int_col)
        @test !has_column(complex_schema, :nonexistent_col)
        
        # Test column validation with complex requirements
        required_subset = Set([:int_col, :float_col])
        @test validate_required_columns(complex_schema, required_subset) === nothing
        
        # Test validation failure with detailed error
        missing_subset = Set([:int_col, :fake_col1, :fake_col2])
        @test_throws ErrorException validate_required_columns(complex_schema, missing_subset)
    end
end