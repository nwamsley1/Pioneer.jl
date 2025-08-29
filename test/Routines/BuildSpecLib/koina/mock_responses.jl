# Mock HTTP responses for Koina API testing
# This allows testing without making actual network requests

using HTTP, JSON

# Mock responses for each model type based on actual Koina API response format
const MOCK_RESPONSES = Dict{String, Dict{String, Any}}(
    "unispec" => Dict(
        "outputs" => [
            Dict(
                "name" => "annotation",
                "shape" => [1, 20],
                "datatype" => "BYTES",
                "data" => ["b2^1", "b3^1", "y4^1", "y5^1", "b4^1", "y6^1", "b5^1", "y7^1", 
                          "b6^1", "y8^1", "y9^1", "b7^1", "y10^1", "b8^1", "y11^1", 
                          "b9^1", "y12^1", "b10^1", "y13^1", "y14^1"]
            ),
            Dict(
                "name" => "mz",
                "shape" => [1, 20],
                "datatype" => "FP32",
                "data" => [200.1f0, 300.2f0, 400.3f0, 500.4f0, 350.2f0, 600.5f0, 450.3f0, 
                          700.6f0, 550.4f0, 800.7f0, 900.8f0, 650.5f0, 1000.9f0, 750.6f0,
                          1100.0f0, 850.7f0, 1200.1f0, 950.8f0, 1300.2f0, 1400.3f0]
            ),
            Dict(
                "name" => "intensities",
                "shape" => [1, 20],
                "datatype" => "FP32", 
                "data" => [0.1f0, 0.8f0, 0.5f0, 0.9f0, 0.3f0, 0.7f0, 0.6f0, 0.85f0, 
                          0.4f0, 0.75f0, 0.65f0, 0.2f0, 0.95f0, 0.35f0, 0.55f0,
                          0.45f0, 0.25f0, 0.15f0, 0.05f0, 0.02f0]
            )
        ]
    ),
    
    "prosit_2020_hcd" => Dict(
        "outputs" => [
            Dict(
                "name" => "annotation",
                "shape" => [1, 29], 
                "datatype" => "BYTES",
                "data" => ["b1+1", "b2+1", "b3+1", "b4+1", "b5+1", "b6+1", "b7+1",
                          "y1+1", "y2+1", "y3+1", "y4+1", "y5+1", "y6+1", "y7+1",
                          "b1+2", "b2+2", "b3+2", "b4+2", "b5+2", "b6+2", "b7+2",
                          "y1+2", "y2+2", "y3+2", "y4+2", "y5+2", "y6+2", "y7+2", "p+1"]
            ),
            Dict(
                "name" => "mz", 
                "shape" => [1, 29],
                "datatype" => "FP32",
                "data" => [150.0f0, 250.1f0, 350.2f0, 450.3f0, 550.4f0, 650.5f0, 750.6f0,
                          147.1f0, 234.1f0, 347.2f0, 448.3f0, 561.4f0, 648.5f0, 761.6f0,
                          75.5f0, 125.6f0, 175.6f0, 225.7f0, 275.7f0, 325.8f0, 375.8f0,
                          74.0f0, 117.6f0, 174.1f0, 224.7f0, 281.2f0, 324.8f0, 381.3f0, 800.9f0]
            ),
            Dict(
                "name" => "intensities",
                "shape" => [1, 29], 
                "datatype" => "FP32",
                "data" => [0.001f0, 0.5f0, 0.9f0, 0.7f0, 0.8f0, 0.6f0, 0.4f0,
                          0.05f0, 0.6f0, 0.85f0, 0.75f0, 0.95f0, 0.65f0, 0.45f0,
                          0.002f0, 0.01f0, 0.03f0, 0.02f0, 0.015f0, 0.008f0, 0.005f0,
                          0.001f0, 0.02f0, 0.04f0, 0.03f0, 0.025f0, 0.018f0, 0.012f0, 0.1f0]
            )
        ]
    ),

    "altimeter" => Dict(
        "outputs" => [
            Dict(
                "name" => "annotations", 
                "shape" => [1, 4, 20],
                "datatype" => "INT32",
                "data" => collect(1:80)  # Flattened indices: [1,2,3,4, 5,6,7,8, ...]
            ),
            Dict(
                "name" => "mz",
                "shape" => [1, 20],
                "datatype" => "FP32", 
                "data" => [200.1f0, 300.2f0, 400.3f0, 500.4f0, 350.2f0, 600.5f0, 450.3f0,
                          700.6f0, 550.4f0, 800.7f0, 900.8f0, 650.5f0, 1000.9f0, 750.6f0,
                          1100.0f0, 850.7f0, 1200.1f0, 950.8f0, 1300.2f0, 1400.3f0]
            ),
            Dict(
                "name" => "coefficients",
                "shape" => [1, 4, 20],
                "datatype" => "FP32",
                "data" => [0.1f0, 0.2f0, 0.3f0, 0.4f0, 0.15f0, 0.25f0, 0.35f0, 0.45f0,
                          0.11f0, 0.21f0, 0.31f0, 0.41f0, 0.16f0, 0.26f0, 0.36f0, 0.46f0,
                          0.12f0, 0.22f0, 0.32f0, 0.42f0, 0.17f0, 0.27f0, 0.37f0, 0.47f0,
                          0.13f0, 0.23f0, 0.33f0, 0.43f0, 0.18f0, 0.28f0, 0.38f0, 0.48f0,
                          0.14f0, 0.24f0, 0.34f0, 0.44f0, 0.19f0, 0.29f0, 0.39f0, 0.49f0,
                          0.101f0, 0.201f0, 0.301f0, 0.401f0, 0.151f0, 0.251f0, 0.351f0, 0.451f0,
                          0.102f0, 0.202f0, 0.302f0, 0.402f0, 0.152f0, 0.252f0, 0.352f0, 0.452f0,
                          0.103f0, 0.203f0, 0.303f0, 0.403f0, 0.153f0, 0.253f0, 0.353f0, 0.453f0,
                          0.104f0, 0.204f0, 0.304f0, 0.404f0, 0.154f0, 0.254f0, 0.354f0, 0.454f0,
                          0.105f0, 0.205f0, 0.305f0, 0.405f0, 0.155f0, 0.255f0, 0.355f0, 0.455f0]
            ),
            Dict(
                "name" => "knots",
                "shape" => [6],
                "datatype" => "FP32", 
                "data" => [10.0f0, 20.0f0, 30.0f0, 40.0f0, 50.0f0, 60.0f0]
            )
        ]
    ),

    "AlphaPeptDeep" => Dict(
        "outputs" => [
            Dict(
                "name" => "annotation",
                "shape" => [1, 24],
                "datatype" => "BYTES",
                "data" => ["b1+1", "b2+1", "b3+1", "b4+1", "b5+1", "b6+1", 
                          "y1+1", "y2+1", "y3+1", "y4+1", "y5+1", "y6+1",
                          "b1+2", "b2+2", "b3+2", "b4+2", "b5+2", "b6+2",
                          "y1+2", "y2+2", "y3+2", "y4+2", "y5+2", "y6+2"]
            ),
            Dict(
                "name" => "mz",
                "shape" => [1, 24],
                "datatype" => "FP32",
                "data" => [150.0f0, 250.1f0, 350.2f0, 450.3f0, 550.4f0, 650.5f0,
                          147.1f0, 234.1f0, 347.2f0, 448.3f0, 561.4f0, 648.5f0,
                          75.5f0, 125.6f0, 175.6f0, 225.7f0, 275.7f0, 325.8f0,
                          74.0f0, 117.6f0, 174.1f0, 224.7f0, 281.2f0, 324.8f0]
            ),
            Dict(
                "name" => "intensities", 
                "shape" => [1, 24],
                "datatype" => "FP32",
                "data" => [0.1f0, 0.7f0, 0.9f0, 0.6f0, 0.8f0, 0.5f0,
                          0.15f0, 0.75f0, 0.85f0, 0.65f0, 0.95f0, 0.55f0,
                          0.05f0, 0.25f0, 0.35f0, 0.2f0, 0.3f0, 0.15f0,
                          0.08f0, 0.28f0, 0.38f0, 0.23f0, 0.33f0, 0.18f0]
            )
        ]
    )
)

# Mock HTTP module for testing - replaces actual HTTP requests
module MockHTTP
    import ..MOCK_RESPONSES
    using JSON
    
    # Track failure count for testing retry logic
    mutable struct MockState
        failure_count::Int
        max_failures::Int
    end
    const mock_state = MockState(0, 0)
    
    function set_mock_failures(count::Int)
        mock_state.failure_count = 0
        mock_state.max_failures = count
    end
    
    function extract_model_from_url(url::String)
        if occursin("UniSpec", url)
            return "unispec"
        elseif occursin("Prosit_2020_intensity_HCD", url)
            return "prosit_2020_hcd"
        elseif occursin("Altimeter", url)
            return "altimeter"
        elseif occursin("AlphaPeptDeep", url)
            return "AlphaPeptDeep"
        else
            error("Unknown model URL: $url")
        end
    end
    
    function post(url::String; body::String)
        # Simulate failures for retry testing
        if mock_state.failure_count < mock_state.max_failures
            mock_state.failure_count += 1
            error("Simulated network failure")
        end
        
        try
            request = JSON.parse(body)
            model = extract_model_from_url(url)
            
            if !haskey(MOCK_RESPONSES, model)
                error("Unknown model: $model")
            end
            
            response = MOCK_RESPONSES[model]
            return (body = JSON.json(response), status = 200)
        catch e
            # Return error response for malformed requests
            error_response = Dict("error" => "Malformed request: $(e)")
            return (body = JSON.json(error_response), status = 400)
        end
    end
end

# Helper function to enable mock HTTP for testing
function with_mock_http(f)
    # Store original HTTP.post function
    original_post = HTTP.post
    
    try
        # Replace HTTP.post with mock version
        eval(:(HTTP.post = MockHTTP.post))
        f()
    finally
        # Restore original function
        eval(:(HTTP.post = $original_post))
    end
end

# Helper function to test retry behavior
function with_mock_failures(failure_count::Int, f)
    MockHTTP.set_mock_failures(failure_count)
    try
        with_mock_http(f)
    finally
        MockHTTP.set_mock_failures(0)  # Reset
    end
end