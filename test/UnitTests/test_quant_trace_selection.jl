# Regression tests for separate-trace quantification.
#
# `trace_mode: "separate"` crashed in IntegrateChromatogramsSearch with
# `ArgumentError: column name :isotopes_captured not found in the data frame`.
# `:isotopes_captured` had been retired from the PSM schema — a passing_psms
# table carries 118 columns and that is not among them — but
# `apply_quant_trace_selection!` still wrote the chosen isotope trace into it.
# Combined mode never calls that function, so only separate-trace runs broke,
# on every dataset.
#
# The key is now returned rather than stored. It has exactly one live consumer,
# `chromatogram_index_key(::SeperateTraces, ...)`, which uses it to look rows up
# in the chromatogram index; nothing reads it back off the frame. So a column
# would add a per-PSM field that only one function writes and only one function
# reads, and would make the separate-trace output schema differ from combined
# for no gain.
#
# Run standalone: julia --project=. test/UnitTests/test_quant_trace_selection.jl
if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
end

@testset "separate-trace quant selection returns the key without adding a column" begin
    # Deliberately without :isotopes_captured — the real passing_psms schema.
    psms = DataFrame(
        precursor_idx = UInt32[10, 20, 30],
        precursor_fraction_transmitted = Float32[0.10, 0.10, 0.10],
    )
    before = Set(Symbol.(names(psms)))

    selected = Dict{UInt32, Tuple{Tuple{Int8, Int8}, Float32}}(
        UInt32(10) => ((Int8(0), Int8(1)), 0.80f0),
        UInt32(20) => ((Int8(1), Int8(2)), 0.65f0),
        # 30 deliberately absent: a precursor with no chromatogram row.
    )

    keys_out = Pioneer.apply_quant_trace_selection!(psms, selected)

    # The crash was writing this onto the frame; it must stay off it.
    @test Set(Symbol.(names(psms))) == before
    @test !hasproperty(psms, :isotopes_captured)

    @test keys_out isa Vector{Tuple{Int8, Int8}}
    @test length(keys_out) == nrow(psms)
    @test keys_out[1] == (Int8(0), Int8(1))
    @test keys_out[2] == (Int8(1), Int8(2))
    # Unmatchable, so the precursor finds no chromatogram group and keeps zero
    # area — not quantifiable, which is correct when nothing was extracted.
    @test keys_out[3] == (Int8(-1), Int8(-1))

    # This one IS a real output column and must reflect the selected trace.
    @test psms.precursor_fraction_transmitted[1] ≈ 0.80f0
    @test psms.precursor_fraction_transmitted[2] ≈ 0.65f0
    # Untouched where nothing was selected.
    @test psms.precursor_fraction_transmitted[3] ≈ 0.10f0
end

# The transmission-ranking rule itself is covered by
# test_chromatogram_integration_trace_order.jl, which also pins the no-column
# contract against a fixture matching the real PSM schema.
