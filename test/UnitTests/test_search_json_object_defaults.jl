using JSON
using Test
using Pioneer
using DataStructures: OrderedDict

@testset "Search config parsing accepts JSON.Object inputs" begin
    mktempdir() do dir
        config_path = joinpath(dir, "search_config.json")
        config = Dict(
            "paths" => Dict(
                "library" => joinpath(dir, "library.poin"),
                "ms_data" => joinpath(dir, "raw"),
                "results" => joinpath(dir, "results"),
            ),
        )
        open(config_path, "w") do io
            JSON.print(io, config)
        end

        ordered_user_params = OrderedDict(
            "paths" => OrderedDict(
                "library" => config["paths"]["library"],
                "ms_data" => config["paths"]["ms_data"],
                "results" => config["paths"]["results"],
            ),
        )
        defaults = Pioneer.get_default_parameters(true)
        @test Pioneer.merge_with_defaults(ordered_user_params, defaults)["paths"]["library"] == config["paths"]["library"]

        @test_nowarn Pioneer.checkParams(config_path)

        parsed = Pioneer.parse_pioneer_parameters(config_path)
        @test parsed.paths.library == config["paths"]["library"]
        @test parsed.paths.ms_data == config["paths"]["ms_data"]
        @test parsed.paths.results == config["paths"]["results"]
    end
end
