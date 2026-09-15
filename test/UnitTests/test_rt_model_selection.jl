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

# Retention-time model selection: `library_params.rt_model` picks the Koina
# model that predicts retention times, and the parse hides the difference in
# output tensor names between models.

const RT_SELECTION_PARAMS =
    joinpath(@__DIR__, "..", "..", "data", "test_build_spec_lib",
             "scenario_b_standard", "params_altimeter.json")

_rt_response(name, values) = Dict{String, Any}(
    "id" => "0",
    "outputs" => Any[Dict{String, Any}(
        "name" => name, "datatype" => "FP32", "shape" => Any[length(values), 1],
        "data" => values)])

@testset "retention time model selection" begin

@testset "registry" begin
    @test Pioneer.DEFAULT_RT_MODEL == "chronologer"
    for name in keys(Pioneer.RT_MODEL_CONFIGS)
        # Every RT model must have an endpoint to send the request to.
        @test haskey(Pioneer.KOINA_URLS, name)
    end
    @test Pioneer.RT_MODEL_CONFIGS["chronologer"].output == :rt
    @test Pioneer.RT_MODEL_CONFIGS["prosit_2024_irt_ptm"].output == :irt
end

@testset "the request is the same for every RT model" begin
    table = DataFrame(koina_sequence = ["PEPTIDEK", "AEGC[UNIMOD:4]DSPK"])
    for name in keys(Pioneer.RT_MODEL_CONFIGS)
        batches = Pioneer.prepare_koina_batch(Pioneer.RetentionTimeModel(name), table; batch_size = 10)
        @test length(batches) == 1
        body = Pioneer.JSON.parse(batches[1])
        @test length(body["inputs"]) == 1
        @test body["inputs"][1]["name"] == "peptide_sequences"
        @test body["inputs"][1]["data"] == ["PEPTIDEK", "AEGC[UNIMOD:4]DSPK"]
    end
end

@testset "parse lands each model's output in :rt" begin
    chrono = Pioneer.parse_koina_batch(Pioneer.RetentionTimeModel("chronologer"),
                                       _rt_response("rt", [1.5, 2.5]))
    @test chrono.fragments.rt == Float32[1.5, 2.5]

    prosit = Pioneer.parse_koina_batch(Pioneer.RetentionTimeModel("prosit_2024_irt_ptm"),
                                       _rt_response("irt", [-17.0, 84.9]))
    @test prosit.fragments.rt == Float32[-17.0, 84.9]

    @testset "a response without the model's output tensor is an error, not an empty column" begin
        @test_throws ErrorException Pioneer.parse_koina_batch(
            Pioneer.RetentionTimeModel("prosit_2024_irt_ptm"), _rt_response("rt", [1.0]))
    end
end

@testset "predict_rt_koina rejects an unknown model before any request" begin
    @test_throws ErrorException Pioneer.predict_rt_koina(
        DataFrame(koina_sequence = ["PEPTIDEK"]); rt_model = "no_such_model")
end

@testset "check_params_bsp" begin
    params = Pioneer.JSON.parsefile(RT_SELECTION_PARAMS, dicttype = Dict{String, Any})

    @testset "absent rt_model defaults to chronologer and is recorded" begin
        checked = Pioneer.check_params_bsp(Pioneer.JSON.json(params))
        @test checked["library_params"]["rt_model"] == "chronologer"
    end

    @testset "a registered rt_model is kept" begin
        p = deepcopy(params)
        p["library_params"]["rt_model"] = "prosit_2024_irt_ptm"
        checked = Pioneer.check_params_bsp(Pioneer.JSON.json(p))
        @test checked["library_params"]["rt_model"] == "prosit_2024_irt_ptm"
    end

    @testset "an unknown rt_model is rejected" begin
        p = deepcopy(params)
        p["library_params"]["rt_model"] = "ssrcalc"
        @test_throws Pioneer.InvalidParametersError Pioneer.check_params_bsp(Pioneer.JSON.json(p))
        p["library_params"]["rt_model"] = 3
        @test_throws Pioneer.InvalidParametersError Pioneer.check_params_bsp(Pioneer.JSON.json(p))
    end
end

end # retention time model selection
