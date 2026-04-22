using DataStructures: OrderedDict

@testset "JSON object compatibility" begin
    @testset "Search defaults merge accepts AbstractDict inputs" begin
        defaults = Dict{String, Any}(
            "global" => Dict{String, Any}(
                "scoring" => Dict{String, Any}(
                    "q_value_threshold" => 0.01,
                    "max_pep" => 0.5,
                ),
            ),
            "paths" => Dict{String, Any}(
                "results" => "default-results",
            ),
        )

        user_params = OrderedDict{String, Any}(
            "global" => OrderedDict{String, Any}(
                "scoring" => OrderedDict{String, Any}(
                    "q_value_threshold" => 0.05,
                ),
            ),
        )

        merged = Pioneer.merge_with_defaults(user_params, defaults)

        @test merged["global"]["scoring"]["q_value_threshold"] == 0.05
        @test merged["global"]["scoring"]["max_pep"] == 0.5
        @test merged["paths"]["results"] == "default-results"
    end

    @testset "Build defaults merge accepts AbstractDict inputs" begin
        defaults = Dict{String, Any}(
            "library_params" => Dict{String, Any}(
                "max_frag_rank" => 12,
                "prediction_model" => "unispec",
            ),
            "predict_fragments" => true,
        )

        user_params = OrderedDict{String, Any}(
            "library_params" => OrderedDict{String, Any}(
                "max_frag_rank" => 20,
            ),
            "predict_fragments" => false,
        )

        merged = Pioneer.merge_with_build_defaults(user_params, defaults)

        @test merged["library_params"]["max_frag_rank"] == 20
        @test merged["library_params"]["prediction_model"] == "unispec"
        @test merged["predict_fragments"] == false
    end
end
