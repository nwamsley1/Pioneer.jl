using Arrow
using Test

using Pioneer: SetPrecursors, getNumEnzymaticTermini,
    getNumVariableModifications
using Pioneer: _num_variable_modifications_at,
    _configured_variable_mod_names

function _write_minimal_precursor_table(path; enzymatic_termini = nothing)
    columns = (
        sequence = ["PEPTIDEK", "SEMIPEP"],
        accession_numbers = ["P1", "P2"],
    )
    table = enzymatic_termini === nothing ? columns : merge(
        columns,
        (num_enzymatic_termini = UInt8.(enzymatic_termini),),
    )
    Arrow.write(path, table)
    return nothing
end

function _write_variable_modification_precursor_table(
    path;
    num_variable_modifications = nothing
)
    columns = (
        sequence = ["FIXEDC", "OXM", "PHOSPHO"],
        accession_numbers = ["P1", "P2", "P3"],
        structural_mods = Union{Missing, String}[
            "(5,C,Unimod:4)",
            "(2,M,Unimod:35)",
            "(3,S,Unimod:21)(7,C,Unimod:4)",
        ],
    )
    table = num_variable_modifications === nothing ? columns : merge(
        columns,
        (
            num_variable_modifications =
                UInt8.(num_variable_modifications),
        ),
    )
    Arrow.write(path, table)
    return nothing
end

@testset "enzymatic termini precursor metadata" begin
    mktempdir() do dir
        current_path = joinpath(dir, "current.arrow")
        _write_minimal_precursor_table(current_path; enzymatic_termini = [2, 1])
        current = SetPrecursors(Arrow.Table(current_path))
        current_values = getNumEnzymaticTermini(current)

        @test collect(current_values) == UInt8[2, 1]

        legacy_path = joinpath(dir, "legacy.arrow")
        _write_minimal_precursor_table(legacy_path)
        legacy = SetPrecursors(Arrow.Table(legacy_path))
        legacy_values = getNumEnzymaticTermini(legacy)

        @test legacy_values == UInt8[2, 2]
    end
end

@testset "variable-modification precursor metadata" begin
    mktempdir() do dir
        current_path = joinpath(dir, "current.arrow")
        _write_variable_modification_precursor_table(
            current_path;
            num_variable_modifications = [0, 1, 1]
        )
        current = SetPrecursors(Arrow.Table(current_path))
        @test collect(getNumVariableModifications(current)) == UInt8[0, 1, 1]

        inferred_path = joinpath(dir, "inferred.arrow")
        _write_variable_modification_precursor_table(inferred_path)
        inferred = SetPrecursors(
            Arrow.Table(inferred_path);
            variable_mod_names = Set(["Unimod:35", "Unimod:21"])
        )
        inferred_values = getNumVariableModifications(inferred)
        @test inferred_values == UInt8[0, 1, 1]

        legacy = SetPrecursors(Arrow.Table(inferred_path))
        legacy_values = getNumVariableModifications(legacy)
        structural_mods = Arrow.Table(inferred_path).structural_mods
        @test legacy_values === nothing
        @test _num_variable_modifications_at(
            legacy_values,
            structural_mods,
            1
        ) == UInt8(0)
        @test _num_variable_modifications_at(
            legacy_values,
            structural_mods,
            2
        ) == UInt8(1)
        # Legacy fallback is deliberately M-oxidation-only and silent.
        @test _num_variable_modifications_at(
            legacy_values,
            structural_mods,
            3
        ) == UInt8(0)

        @test _configured_variable_mod_names(Dict(
            "variable_mods" => Dict("name" => ["Unimod:4"]),
            "fixed_mods" => Dict("name" => ["Unimod:4"]),
        )) === nothing
    end
end
