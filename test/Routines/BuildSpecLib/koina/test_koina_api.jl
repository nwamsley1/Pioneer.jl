# Core Koina API functionality tests
# Tests the basic API request/response handling with real Koina API calls

@testset "Koina API Core Functions" begin
    
    @testset "Basic API Request Success" begin
        # Create minimal valid request data for UniSpec
        test_request = Dict(
            "id" => "test_request",
            "inputs" => [
                Dict("name" => "peptide_sequences", "data" => ["PEPTIDE"], "shape" => [1, 1], "datatype" => "BYTES"),
                Dict("name" => "precursor_charges", "data" => [2], "shape" => [1, 1], "datatype" => "INT32"),
                Dict("name" => "collision_energies", "data" => [25.0f0], "shape" => [1, 1], "datatype" => "FP32"),
                Dict("name" => "instrument_types", "data" => ["QE"], "shape" => [1, 1], "datatype" => "BYTES")
            ]
        )
        json_data = JSON.json(test_request)
        test_url = KOINA_URLS["unispec"]
        
        response = make_koina_request(json_data, test_url, max_attempts=3)
        @test isa(response, Dict)
        @test haskey(response, "outputs")
        @test length(response["outputs"]) > 0
        
        # Check that we get the expected UniSpec response structure
        outputs = response["outputs"]
        output_names = [out["name"] for out in outputs]
        @test "annotation" in output_names
        @test "mz" in output_names  
        @test "intensities" in output_names
    end
    
    @testset "Request Error Handling" begin
        # Test with malformed JSON should throw error
        malformed_json = "{ invalid json"
        test_url = KOINA_URLS["unispec"]
        @test_throws Exception make_koina_request(malformed_json, test_url, max_attempts=1)
        
        # Test with unknown model URL should throw error
        valid_request = Dict(
            "id" => "test",
            "inputs" => [
                Dict("name" => "peptide_sequences", "data" => ["PEPTIDE"], "shape" => [1, 1], "datatype" => "BYTES")
            ]
        )
        valid_json = JSON.json(valid_request)
        unknown_url = "https://koina.wilhelmlab.org/v2/models/UnknownModel/infer"
        @test_throws Exception make_koina_request(valid_json, unknown_url, max_attempts=1)
    end
    
    @testset "Different Model URLs" begin
        # Test prosit model (doesn't need instrument types)
        prosit_request = Dict(
            "id" => "prosit_test", 
            "inputs" => [
                Dict("name" => "peptide_sequences", "data" => ["PEPTIDE"], "shape" => [1, 1], "datatype" => "BYTES"),
                Dict("name" => "precursor_charges", "data" => [2], "shape" => [1, 1], "datatype" => "INT32"),
                Dict("name" => "collision_energies", "data" => [25.0f0], "shape" => [1, 1], "datatype" => "FP32")
            ]
        )
        json_data = JSON.json(prosit_request)
        
        response = make_koina_request(json_data, KOINA_URLS["prosit_2020_hcd"], max_attempts=3)
        @test haskey(response, "outputs")
        @test length(response["outputs"]) > 0
        
        # Prosit should return different structure than UniSpec
        output_names = [out["name"] for out in response["outputs"]]
        @test "mz" in output_names  # All models should have mz
    end
    
    @testset "Request Timing and Performance" begin
        # Test that requests complete in reasonable time
        test_request = Dict(
            "id" => "timing_test", 
            "inputs" => [
                Dict("name" => "peptide_sequences", "data" => ["PEPTIDE"], "shape" => [1, 1], "datatype" => "BYTES"),
                Dict("name" => "precursor_charges", "data" => [2], "shape" => [1, 1], "datatype" => "INT32"),
                Dict("name" => "collision_energies", "data" => [25.0f0], "shape" => [1, 1], "datatype" => "FP32")
            ]
        )
        json_data = JSON.json(test_request)
        
        # Test that prosit responds quickly
        start_time = time()
        response = make_koina_request(json_data, KOINA_URLS["prosit_2020_hcd"], max_attempts=3)
        elapsed = time() - start_time
        
        @test haskey(response, "outputs")
        @test elapsed < 30.0  # Should complete within 30 seconds
    end
    
    @testset "Response Structure Validation" begin
        # Test prosit response structure in detail
        test_request = Dict(
            "id" => "validation_test", 
            "inputs" => [
                Dict("name" => "peptide_sequences", "data" => ["PEPTIDE"], "shape" => [1, 1], "datatype" => "BYTES"),
                Dict("name" => "precursor_charges", "data" => [2], "shape" => [1, 1], "datatype" => "INT32"),
                Dict("name" => "collision_energies", "data" => [25.0f0], "shape" => [1, 1], "datatype" => "FP32")
            ]
        )
        json_data = JSON.json(test_request)
        
        response = make_koina_request(json_data, KOINA_URLS["prosit_2020_hcd"], max_attempts=3)
        @test haskey(response, "outputs")
        
        outputs = response["outputs"]
        @test length(outputs) >= 3  # Should have annotation, mz, intensities at minimum
        
        # Find each output by name
        annotation_output = nothing
        mz_output = nothing  
        intensity_output = nothing
        
        for output in outputs
            if output["name"] == "annotation"
                annotation_output = output
            elseif output["name"] == "mz"
                mz_output = output
            elseif output["name"] == "intensities"
                intensity_output = output
            end
        end
        
        @test annotation_output !== nothing
        @test mz_output !== nothing
        @test intensity_output !== nothing
        
        # Check data types
        @test annotation_output["datatype"] == "BYTES"
        @test mz_output["datatype"] == "FP32"
        @test intensity_output["datatype"] == "FP32"
        
        # Check that we have actual data
        @test length(annotation_output["data"]) > 0
        @test length(mz_output["data"]) > 0
        @test length(intensity_output["data"]) > 0
        
        # Check that annotations are strings
        @test all(isa(ann, String) for ann in annotation_output["data"])
        # Check that mz and intensities are numbers
        @test all(isa(val, Number) for val in mz_output["data"])
        @test all(isa(val, Number) for val in intensity_output["data"])
    end
end