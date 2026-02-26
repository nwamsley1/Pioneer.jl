using Test
using Dictionaries
using Pioneer: harmonize_charge_state_irts!

@testset "charge states share best-scoring iRT" begin
    precursor_dict = Dictionary{UInt32, NamedTuple{(:best_prob, :best_ms_file_idx, :best_scan_idx, :best_irt, :mean_irt, :var_irt, :n, :mz), Tuple{Float32, UInt32, UInt32, Float32, Union{Missing, Float32}, Union{Missing, Float32}, Union{Missing, UInt16}, Float32}}}()

    # Same sequence, different charge states (different precursor IDs)
    insert!(precursor_dict, UInt32(1), (
        best_prob = 0.80f0,
        best_ms_file_idx = UInt32(1),
        best_scan_idx = UInt32(100),
        best_irt = 30.0f0,
        mean_irt = missing,
        var_irt = missing,
        n = missing,
        mz = 500.0f0
    ))
    insert!(precursor_dict, UInt32(2), (
        best_prob = 0.95f0,
        best_ms_file_idx = UInt32(1),
        best_scan_idx = UInt32(101),
        best_irt = 32.0f0,
        mean_irt = missing,
        var_irt = missing,
        n = missing,
        mz = 333.0f0
    ))

    # Different sequence should remain unchanged
    insert!(precursor_dict, UInt32(3), (
        best_prob = 0.70f0,
        best_ms_file_idx = UInt32(2),
        best_scan_idx = UInt32(200),
        best_irt = 45.0f0,
        mean_irt = missing,
        var_irt = missing,
        n = missing,
        mz = 700.0f0
    ))

    sequence_lookup = ["", "PEPTIDE", "PEPTIDE", "OTHERSEQ"]

    harmonize_charge_state_irts!(precursor_dict, sequence_lookup)

    @test precursor_dict[UInt32(1)][:best_irt] == 32.0f0
    @test precursor_dict[UInt32(2)][:best_irt] == 32.0f0
    @test precursor_dict[UInt32(3)][:best_irt] == 45.0f0

    # Sanity: only iRT changed for precursor 1
    @test precursor_dict[UInt32(1)][:best_prob] == 0.80f0
    @test precursor_dict[UInt32(1)][:best_scan_idx] == UInt32(100)
end
