
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