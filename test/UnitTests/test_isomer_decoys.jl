# Unit tests for isomer (localization) decoy core logic.
# Hand-verified on VLSAGSPESIK (L=11, S at 3,6,9): move phospho S6 -> P7 decoy.
# See dev_docs/prosit_ptm_integration/localization_decoy_library_plan.md.
using Pioneer, Test
const _P = Pioneer
const PHOS = 79.966331

@testset "fragment_contains (b/y span)" begin
    L = 11
    @test _P.fragment_contains('b', 6, L, 6) == true    # b6 spans 1..6
    @test _P.fragment_contains('b', 6, L, 7) == false
    @test _P.fragment_contains('b', 7, L, 7) == true    # b7 spans 1..7
    @test _P.fragment_contains('y', 5, L, 7) == true    # y5 spans 7..11
    @test _P.fragment_contains('y', 5, L, 6) == false
    @test _P.fragment_contains('p', 0, L, 6) == true    # precursor spans all
end

@testset "permute_fragment_mz S6->P7 (charge 1)" begin
    L, kold, knew = 11, 6, 7
    mz0 = 1000.0f0
    d = Float32(PHOS)
    # b6 contains 6 not 7 -> loses the phospho -> -Δ
    @test _P.permute_fragment_mz(mz0, 'b', 6, 1, L, kold, knew, PHOS) ≈ mz0 - d
    # y5 contains 7 not 6 -> gains the phospho -> +Δ
    @test _P.permute_fragment_mz(mz0, 'y', 5, 1, L, kold, knew, PHOS) ≈ mz0 + d
    # b7 contains both -> unchanged ; b5 contains neither -> unchanged
    @test _P.permute_fragment_mz(mz0, 'b', 7, 1, L, kold, knew, PHOS) == mz0
    @test _P.permute_fragment_mz(mz0, 'b', 5, 1, L, kold, knew, PHOS) == mz0
    # precursor ion unchanged (mod count preserved)
    @test _P.permute_fragment_mz(mz0, 'p', 0, 1, L, kold, knew, PHOS) == mz0
end

@testset "permute_fragment_mz respects charge" begin
    L = 11; mz0 = 500.0f0
    # charge-2 site-determining ion shifts by Δ/2
    @test _P.permute_fragment_mz(mz0, 'b', 6, 2, L, 6, 7, PHOS) ≈ mz0 - Float32(PHOS)/2
end

@testset "choose_decoy_site spacing + sign" begin
    L = 11; acc = UInt8[3, 6, 9]
    @test _P.choose_decoy_site(L, acc, 6, 1, true)  == 7   # 6+1 (P), non-acceptor
    @test _P.choose_decoy_site(L, acc, 6, 1, false) == 5   # 6-1 (G)
    @test _P.choose_decoy_site(L, acc, 3, 1, false) == 2   # 3-1
    @test _P.choose_decoy_site(L, acc, 9, 1, true)  == 10  # 9+1
    # spacing 2 from S6 -> 8 (E) with +; 4 (A) with -
    @test _P.choose_decoy_site(L, acc, 6, 2, true)  == 8
    @test _P.choose_decoy_site(L, acc, 6, 2, false) == 4
    # if preferred sign lands on an acceptor, fall back to the other sign.
    # from S6 spacing 3: +->9 (acceptor, skip) -> -->3 (acceptor, skip) -> 0
    @test _P.choose_decoy_site(L, acc, 6, 3, true)  == 0
end

@testset "decoy sites are covered by retention candidate set" begin
    # parity: every site choose_decoy_site can pick must be in decoy_neighbor_positions
    seq = "VLSAGSPESIK"; acc = UInt8[3, 6, 9]
    for sp in (1, 2)
        cover = Set(_P.decoy_neighbor_positions(seq, acc; spacing = sp))
        for k in acc, pp in (true, false)
            d = _P.choose_decoy_site(length(seq), acc, Int(k), sp, pp)
            d == 0 || @test UInt8(d) in cover
        end
    end
end
