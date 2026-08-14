# Unit tests for the fixed/variable modification site conflict check
# (src/Routines/BuildSpecLib/utils/check_params.jl).
#
# A fixed modification occupies every matching residue, so a variable
# modification on the same residue describes a peptide that cannot exist.
# Pioneer does not notice on its own: fillVarModStrings! copies the fixed mods
# and pushes the variable ones on top without checking positions, so the residue
# ends up modified twice and the library carries impossible masses. Nothing
# downstream throws, which is why this is checked at the parameter boundary.
#
# Run standalone: julia --project=. test/UnitTests/test_mod_site_conflicts.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

const MOD_CONFLICT_PARAMS =
    joinpath(@__DIR__, "..", "..", "data", "test_build_spec_lib",
             "scenario_b_standard", "params_altimeter.json")

_mods(names, patterns, masses) =
    Dict{String, Any}("name" => names, "pattern" => patterns, "mass" => masses)

_check(fixed, variable) = Pioneer.check_mod_site_conflicts(
    Dict{String, Any}("fixed_mods" => fixed, "variable_mods" => variable))

@testset "modification site conflicts" begin

@testset "mod_pattern_residues" begin
    # Membership is decided by matching each amino acid, so bracket classes and
    # alternations work without parsing regex syntax.
    @test Pioneer.mod_pattern_residues("C") == Set(['C'])
    @test Pioneer.mod_pattern_residues("[KR]") == Set(['K', 'R'])
    @test Pioneer.mod_pattern_residues("[CK]") == Set(['C', 'K'])
    @test Pioneer.mod_pattern_residues("S|T") == Set(['S', 'T'])

    @testset "an unparseable pattern reports nothing rather than throwing" begin
        # Malformed patterns are another check's business; this one must not be
        # the thing that fails on them.
        @test isempty(Pioneer.mod_pattern_residues("["))
    end
end

@testset "accepted configurations" begin
    @testset "disjoint residues" begin
        @test _check(_mods(["Unimod:4"], ["C"], [57.021464]),
                     _mods(["Unimod:35"], ["M"], [15.994915])) === nothing
    end

    @testset "same modification on different residues" begin
        # Carbamidomethyl fixed on C and variable on K is legitimate: the
        # residues do not overlap.
        @test _check(_mods(["Unimod:4"], ["C"], [57.021464]),
                     _mods(["Unimod:4"], ["K"], [57.021464])) === nothing
    end

    @testset "no variable mods at all" begin
        @test _check(_mods(["Unimod:4"], ["C"], [57.021464]),
                     _mods(String[], String[], Float64[])) === nothing
    end

    @testset "missing or malformed sections are left to other checks" begin
        @test Pioneer.check_mod_site_conflicts(Dict{String, Any}()) === nothing
        @test Pioneer.check_mod_site_conflicts(
            Dict{String, Any}("fixed_mods" => "nonsense",
                              "variable_mods" => "nonsense")) === nothing
    end
end

@testset "rejected configurations" begin
    @testset "the same modification on the same residue" begin
        @test_throws Pioneer.InvalidParametersError _check(
            _mods(["Unimod:4"], ["C"], [57.021464]),
            _mods(["Unimod:4"], ["C"], [57.021464]))
    end

    @testset "different modifications competing for one residue" begin
        # The general case: the fixed mod already occupies every C, so any
        # variable mod on C is impossible, not just carbamidomethyl.
        @test_throws Pioneer.InvalidParametersError _check(
            _mods(["Unimod:4"], ["C"], [57.021464]),
            _mods(["Unimod:35"], ["C"], [15.994915]))
    end

    @testset "partial overlap through a residue class" begin
        @test_throws Pioneer.InvalidParametersError _check(
            _mods(["Unimod:4"], ["[CK]"], [57.021464]),
            _mods(["Unimod:4"], ["K"], [57.021464]))
    end

    @testset "a conflict in any position, not only the first" begin
        @test_throws Pioneer.InvalidParametersError _check(
            _mods(["Unimod:4"], ["C"], [57.021464]),
            _mods(["Unimod:35", "Unimod:21"], ["M", "C"], [15.994915, 79.966331]))
    end
end

@testset "the error says what is wrong and how to fix it" begin
    err = try
        _check(_mods(["Unimod:4"], ["C"], [57.021464]),
               _mods(["Unimod:35"], ["C"], [15.994915]))
        nothing
    catch e
        e
    end
    @test err isa Pioneer.InvalidParametersError

    msg = err.message
    @testset "names both modifications and the shared residue" begin
        @test occursin("Unimod:35", msg)      # the variable one
        @test occursin("Unimod:4", msg)       # the fixed one
        @test occursin("C", msg)              # the residue they share
    end

    @testset "explains the consequence and the remedy" begin
        @test occursin("fixed modification already occupies", msg)
        @test occursin(r"Remove|narrow", msg)
    end
end

@testset "wired into check_params_bsp" begin
    # The helper being correct is not enough: the entry point has to call it.
    params = Pioneer.JSON.parsefile(MOD_CONFLICT_PARAMS, dicttype = Dict{String, Any})

    @testset "the stock parameters are accepted" begin
        @test Pioneer.check_params_bsp(Pioneer.JSON.json(params)) isa Dict
    end

    @testset "a conflicting variable mod is rejected" begin
        conflicting = deepcopy(params)
        conflicting["variable_mods"] = _mods(["Unimod:4"], ["C"], [57.021464])
        @test_throws Pioneer.InvalidParametersError Pioneer.check_params_bsp(
            Pioneer.JSON.json(conflicting))
    end
end

end # top-level testset
