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

# A build is refused when its modifications -- or an unmodified cysteine --
# are outside what the fragment model or the retention-time model was trained
# on. See check_model_mod_support.

const MOD_SUPPORT_PARAMS =
    joinpath(@__DIR__, "..", "..", "data", "test_build_spec_lib",
             "scenario_b_standard", "params_altimeter.json")

_ms_mods(names, patterns) = Dict{String, Any}(
    "name" => names, "pattern" => patterns, "mass" => zeros(length(names)))

function _ms_params(; frag = "altimeter", rt = "chronologer",
                     fixed = _ms_mods(["Unimod:4"], ["C"]),
                     variable = _ms_mods(String[], String[]),
                     predict_fragments = true)
    Dict{String, Any}(
        "predict_fragments" => predict_fragments,
        "library_params" => Dict{String, Any}("prediction_model" => frag, "rt_model" => rt),
        "fixed_mods" => fixed,
        "variable_mods" => variable,
    )
end

_ms_error(p) = try
    Pioneer.check_model_mod_support(p); ""
catch e
    e isa Pioneer.InvalidParametersError || rethrow()
    e.message
end

@testset "model modification support" begin

@testset "unimod_id" begin
    @test Pioneer.unimod_id("Unimod:35") == 35
    @test Pioneer.unimod_id("UNIMOD:4") == 4
    @test Pioneer.unimod_id(" unimod:21 ") == 21
    @test Pioneer.unimod_id("21") == 21
    @test Pioneer.unimod_id("Oxidation") === nothing
    @test Pioneer.unimod_id("") === nothing
    @test Pioneer.unimod_id(35) === nothing
end

@testset "mod_pattern_sites" begin
    @test Pioneer.mod_pattern_sites("[STY]") == Set(['S', 'T', 'Y'])
    @test Pioneer.mod_pattern_sites("^.") == Set(['n'])
    @test Pioneer.mod_pattern_sites("^K") == Set(['n'])
end

@testset "every model declares its support" begin
    for (name, c) in Pioneer.MODEL_CONFIGS
        @test c.supported_mods isa Pioneer.ModSupport
        @test c.free_cys isa Bool
    end
    for (name, c) in Pioneer.RT_MODEL_CONFIGS
        @test c.supported_mods isa Pioneer.ModSupport
        @test c.free_cys isa Bool
    end
    @test Pioneer.MODEL_CONFIGS["prosit_2025_40ptm"].free_cys
    @test !Pioneer.MODEL_CONFIGS["altimeter"].free_cys
    @test Pioneer.RT_MODEL_CONFIGS["chronologer"].free_cys
end

@testset "accepted" begin
    @testset "the stock selection on every model pair" begin
        stock = _ms_mods(["Unimod:35"], ["M"])
        for frag in keys(Pioneer.MODEL_CONFIGS), rt in keys(Pioneer.RT_MODEL_CONFIGS)
            @test _ms_error(_ms_params(frag = frag, rt = rt, variable = stock)) == ""
        end
    end

    @testset "phospho with the PTM fragment model and Chronologer" begin
        @test _ms_error(_ms_params(frag = "prosit_2025_40ptm",
                                   variable = _ms_mods(["Unimod:21"], ["[STY]"]))) == ""
    end

    @testset "free cysteine with the 40-PTM model and either RT model" begin
        for rt in keys(Pioneer.RT_MODEL_CONFIGS)
            @test _ms_error(_ms_params(frag = "prosit_2025_40ptm", rt = rt,
                                       fixed = _ms_mods(String[], String[]))) == ""
        end
    end

    @testset "a fixed mod covering C among other residues counts" begin
        # Carbamidomethyl on K is in the Prosit PTM vocabulary but not in
        # Chronologer's, so this needs the Prosit RT model; what matters here is
        # that the [CK] pattern satisfies the cysteine requirement.
        @test _ms_error(_ms_params(frag = "prosit_2025_40ptm", rt = "prosit_2024_irt_ptm",
                                   fixed = _ms_mods(["Unimod:4"], ["[CK]"]))) == ""
        msg = _ms_error(_ms_params(frag = "prosit_2025_40ptm",
                                   fixed = _ms_mods(["Unimod:4"], ["[CK]"])))
        @test occursin("retention-time model chronologer does not support Unimod:4 on K", msg)
        @test !occursin("carbamidomethylated cysteine", msg)
    end

    @testset "resuming from disk predicts nothing, so nothing is checked" begin
        @test _ms_error(_ms_params(fixed = _ms_mods(String[], String[]),
                                   predict_fragments = false)) == ""
    end
end

@testset "refused" begin
    @testset "free cysteine on a model that assumes carbamidomethyl" begin
        msg = _ms_error(_ms_params(fixed = _ms_mods(String[], String[])))
        @test occursin("fragment model altimeter assumes carbamidomethylated cysteine", msg)
        @test !occursin("retention-time model chronologer assumes", msg)   # Chronologer allows it
        @test occursin("Fragment models that would: prosit_2025_40ptm", msg)
    end

    @testset "a PTM the fragment model never saw" begin
        msg = _ms_error(_ms_params(variable = _ms_mods(["Unimod:21"], ["[STY]"])))
        @test occursin("fragment model altimeter does not support Unimod:21 on S, Unimod:21 on T, Unimod:21 on Y", msg)
        @test occursin("prosit_2024_ptm", msg) && occursin("prosit_2025_40ptm", msg)
    end

    @testset "a PTM the fragment model knows but the RT model does not" begin
        # HexNAc on S: in the Prosit PTM vocabulary, not in Chronologer's 17.
        msg = _ms_error(_ms_params(frag = "prosit_2025_40ptm",
                                   variable = _ms_mods(["Unimod:43"], ["[ST]"])))
        @test !occursin("fragment model", msg)
        @test occursin("retention-time model chronologer does not support Unimod:43 on S, Unimod:43 on T", msg)
        @test occursin("Retention-time models that would: prosit_2024_irt_ptm", msg)
        # ...and switching the RT model resolves it.
        @test _ms_error(_ms_params(frag = "prosit_2025_40ptm", rt = "prosit_2024_irt_ptm",
                                   variable = _ms_mods(["Unimod:43"], ["[ST]"]))) == ""
    end

    @testset "a known accession on an unsupported residue" begin
        # Oxidation is supported on M only for the base models.
        msg = _ms_error(_ms_params(variable = _ms_mods(["Unimod:35"], ["[MW]"])))
        @test occursin("does not support Unimod:35 on W", msg)
        @test !occursin("Unimod:35 on M", msg)
    end

    @testset "an N-terminal pattern is checked as the N-terminus" begin
        # Acetyl: Chronologer allows K and N-term; the base fragment models allow neither.
        msg = _ms_error(_ms_params(variable = _ms_mods(["Unimod:1"], ["^."])))
        @test occursin("fragment model altimeter does not support Unimod:1 on N-term", msg)
        @test !occursin("retention-time model chronologer does not support", msg)
    end

    @testset "no model accepts the selection" begin
        msg = _ms_error(_ms_params(variable = _ms_mods(["Unimod:99999"], ["K"])))
        @test occursin("Fragment models that would: none", msg)
        @test occursin("Retention-time models that would: none", msg)
    end

    @testset "a non-UNIMOD name is refused outright" begin
        msg = _ms_error(_ms_params(variable = _ms_mods(["Oxidation"], ["M"])))
        @test occursin("\"Oxidation\" is not a UNIMOD accession", msg)
    end
end

@testset "wired into check_params_bsp" begin
    params = Pioneer.JSON.parsefile(MOD_SUPPORT_PARAMS, dicttype = Dict{String, Any})
    @test Pioneer.check_params_bsp(Pioneer.JSON.json(params)) isa Dict

    unalkylated = deepcopy(params)
    unalkylated["fixed_mods"] = _ms_mods(String[], String[])
    @test_throws Pioneer.InvalidParametersError Pioneer.check_params_bsp(Pioneer.JSON.json(unalkylated))

    unalkylated["library_params"]["prediction_model"] = "prosit_2025_40ptm"
    @test Pioneer.check_params_bsp(Pioneer.JSON.json(unalkylated)) isa Dict
end

end # model modification support
