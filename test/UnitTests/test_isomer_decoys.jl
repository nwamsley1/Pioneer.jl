# Unit tests for isomer (localization) decoy core logic.
# Hand-verified on VLSAGSPESIK (L=11, S at 3,6,9): move phospho S6 -> P7 decoy.
# See dev_docs/prosit_ptm_integration/localization_decoy_library_plan.md.
using Pioneer, Test, DataFrames
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

@testset "generate_isomer_decoys! table emission" begin
    # One real mono-phospho precursor VLSAGSPESIK @ S6; row 1 (odd -> +sign -> d=7,P).
    prec = DataFrame(sequence = ["VLSAGSPESIK"], mods = ["(6,S,Unimod:21)"],
                     mz = Float32[600.0], prec_charge = UInt8[2], pair_id = UInt32[5])
    frag = DataFrame(precursor_idx = UInt32[1, 1, 1, 1], annotation = ["b5", "b6", "b7", "y5"],
                     mz = Float32[500, 500, 500, 500], intensities = Float32[10, 20, 30, 40])
    vm = NamedTuple{(:p, :r), Tuple{Regex, String}}[(p = r"[STY]", r = "Unimod:21")]
    masses = Dict("Unimod:21" => 79.966331f0)

    _P.generate_isomer_decoys!(prec, frag, vm, masses; spacing = 1)

    @test nrow(prec) == 2                              # 1 real + 1 decoy
    @test prec.is_loc_decoy == [false, true]
    @test prec.loc_base_prec_id == UInt32[0, 1]
    @test prec.mz[2] == prec.mz[1]                     # precursor mass preserved
    @test prec.mods[2] == "(7,P,Unimod:21)"            # phospho moved S6 -> P7
    @test prec.sequence[2] == "VLSAGSPESIK"
    @test prec.pair_id[2] != prec.pair_id[1]           # loc-decoy gets a FRESH pair_id
    @test prec.pair_id[2] == UInt32(6)                 # max(5)+1

    dec = frag[frag.precursor_idx .== 2, :]
    d = Float32(PHOS)
    m = Dict(dec.annotation[i] => dec.mz[i] for i in 1:nrow(dec))
    @test nrow(dec) == 4
    @test m["b6"] ≈ 500 - d                            # contains 6 not 7 -> -Δ
    @test m["y5"] ≈ 500 + d                            # contains 7 not 6 -> +Δ
    @test m["b7"] == 500 && m["b5"] == 500             # unchanged
    @test all(frag[frag.precursor_idx .== 1, :mz] .== 500)   # reals untouched
end

@testset "loc-decoy pair_id preserves real target/decoy pairing" begin
    # A real target + its ID-decoy share pair_id=10 (2 rows). Both phospho, so both
    # spawn loc-decoys. After add_pair_indices! the real pair must stay partnered and
    # the loc-decoys must be unpaired (missing).
    prec = DataFrame(
        sequence = ["SAAAAAK", "KAAAAAS"],
        mods = ["(1,S,Unimod:21)", "(7,S,Unimod:21)"],
        mz = Float32[400, 400], prec_charge = UInt8[2, 2],
        is_decoy = [false, true], pair_id = UInt32[10, 10])
    frag = DataFrame(precursor_idx = UInt32[1, 2], annotation = ["b2", "b2"],
                     mz = Float32[200, 200], intensities = Float32[1, 1])
    vm = NamedTuple{(:p, :r), Tuple{Regex, String}}[(p = r"[STY]", r = "Unimod:21")]

    _P.generate_isomer_decoys!(prec, frag, vm, Dict("Unimod:21" => 79.966331f0))
    @test nrow(prec) == 4                               # 2 real + 2 loc-decoy
    @test allunique(prec.pair_id[3:4])                  # fresh ids
    @test all(prec.pair_id[3:4] .> UInt32(10))

    _P.add_pair_indices!(prec)
    @test prec.partner_precursor_idx[1] == 2            # real pair intact
    @test prec.partner_precursor_idx[2] == 1
    @test ismissing(prec.partner_precursor_idx[3])      # loc-decoys unpaired
    @test ismissing(prec.partner_precursor_idx[4])
end

@testset "generate_isomer_decoys! no-op without var mods" begin
    prec = DataFrame(sequence = ["PEPTIDEK"], mods = [""], mz = Float32[400.0])
    frag = DataFrame(precursor_idx = UInt32[1], annotation = ["b2"],
                     mz = Float32[200], intensities = Float32[1])
    empty_vm = NamedTuple{(:p, :r), Tuple{Regex, String}}[]
    _P.generate_isomer_decoys!(prec, frag, empty_vm, Dict("Unimod:21" => 79.966331f0))
    @test nrow(prec) == 1 && all(.!prec.is_loc_decoy)  # no decoys, flag col added
end
