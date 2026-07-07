# Unit tests for isomer (localization) decoys — move-one-random onto shadow sites.
# See dev_docs/prosit_ptm_integration/isomer_decoy_shadow_sites.md.
using Pioneer, Test, DataFrames
const _P = Pioneer
const PHOS = 79.966331
const OXI  = 15.994915
const VM = NamedTuple{(:p, :r), Tuple{Regex, String}}

@testset "fragment_contains + permute_fragment_mz" begin
    L = 11
    @test _P.fragment_contains('b', 6, L, 6) == true
    @test _P.fragment_contains('b', 6, L, 7) == false
    @test _P.fragment_contains('y', 5, L, 7) == true
    @test _P.permute_fragment_mz(1000.0f0, 'b', 6, 1, L, 6, 7, PHOS) ≈ 1000 - Float32(PHOS)
    @test _P.permute_fragment_mz(1000.0f0, 'y', 5, 1, L, 6, 7, PHOS) ≈ 1000 + Float32(PHOS)
    @test _P.permute_fragment_mz(1000.0f0, 'b', 7, 1, L, 6, 7, PHOS) == 1000  # contains both
    @test _P.permute_fragment_mz(500.0f0,  'b', 6, 2, L, 6, 7, PHOS) ≈ 500 - Float32(PHOS)/2
end

# parse a decoy mods string -> set of (position, mod_name)
_sites(s) = Set((parse(Int, m.captures[1]), String(m.captures[3]))
                for m in eachmatch(r"\((\d+),([A-Za-z]),([^)]+)\)", s))

@testset "generate_isomer_decoys!: move-one-random, phospho only" begin
    seq = "VLSAGSPESIK"                      # S at 3,6,9
    vm  = VM[(p = r"[STY]", r = "Unimod:21")]
    shadow = _P.shadow_decoy_sites(seq, _P.real_sites_by_mod(seq, vm); seed = _P.seq_seed(seq))
    shad_p = Set(shadow["Unimod:21"])        # the decoy phospho positions
    @test length(shad_p) == 3                # count-matched to S3,S6,S9

    # peptidoforms: 3 mono + 1 di, all sequence VLSAGSPESIK
    mods = ["(3,S,Unimod:21)", "(6,S,Unimod:21)", "(9,S,Unimod:21)",
            "(3,S,Unimod:21)(6,S,Unimod:21)"]
    prec = DataFrame(sequence = fill(seq, 4), mods = mods,
                     mz = Float32[600,600,600,650], prec_charge = fill(UInt8(2), 4),
                     pair_id = UInt32[1,2,3,4])
    frag = DataFrame(precursor_idx = UInt32[1,2,3,4], annotation = fill("b4", 4),
                     mz = fill(500.0f0, 4), intensities = fill(1.0f0, 4))

    prec, frag = _P.generate_isomer_decoys!(prec, frag, vm, Dict("Unimod:21"=>79.966331f0))

    @test nrow(prec) == 8                     # 4 reals + 4 decoys (all modified)
    dec = prec[prec.is_loc_decoy, :]
    @test nrow(dec) == 4
    @test prec.is_loc_decoy == [false,false,false,false, true,true,true,true]

    for i in 1:nrow(dec)
        sites = _sites(dec.mods[i])
        base = Int(dec.loc_base_prec_id[i])
        bsites = _sites(prec.mods[base])
        # same modification composition (count per mod name)
        @test length(sites) == length(bsites)
        # exactly one mod is on a shadow (decoy) position; the rest on real sites
        onshadow = count(sp -> sp[1] in shad_p, sites)
        @test onshadow == 1
        @test dec.mz[i] == prec.mz[base]      # mass preserved
    end
    # decoy peptidoforms unique
    @test allunique(dec.mods)
    # fresh pair_ids
    @test all(dec.pair_id .> UInt32(4))
    @test allunique(dec.pair_id)
    # decoy fragments: b4 permuted iff it straddles the moved site (mz != 500 or ==)
    @test all(frag[frag.precursor_idx .<= 4, :mz] .== 500.0f0)   # reals untouched
end

@testset "generate_isomer_decoys!: two mods, composition + one-shadow rule" begin
    seq = "SGGSGGMG"                          # S1,S4 phospho ; M7 ox
    vm  = VM[(p = r"[STY]", r = "Unimod:21"), (p = r"M", r = "Unimod:35")]
    rb  = _P.real_sites_by_mod(seq, vm)
    shadow = _P.shadow_decoy_sites(seq, rb; seed = _P.seq_seed(seq))
    shad_all = Set(vcat(values(shadow)...))

    # peptidoform 1: 1 phospho on 2 S sites (C(2,1)=2, ambiguous) + 1 ox on the
    #   only M (C(1,1)=1, NOT ambiguous) -> gets a decoy, but only the phospho moves.
    # peptidoform 2: 2 phospho on both S sites + ox on the only M -> fully saturated
    #   (C(2,2)=1, C(1,1)=1) -> unique at its m/z -> NO decoy.
    mods = ["(1,S,Unimod:21)(7,M,Unimod:35)",           # 1 phospho + 1 ox
            "(1,S,Unimod:21)(4,S,Unimod:21)(7,M,Unimod:35)"]  # 2 phospho + 1 ox (saturated)
    prec = DataFrame(sequence = fill(seq, 2), mods = mods,
                     mz = Float32[500, 600], prec_charge = fill(UInt8(2), 2),
                     pair_id = UInt32[1, 2])
    frag = DataFrame(precursor_idx = UInt32[1, 2], annotation = ["b3", "b3"],
                     mz = Float32[300, 300], intensities = Float32[1, 1])
    masses = Dict("Unimod:21" => 79.966331f0, "Unimod:35" => 15.994915f0)

    prec, frag = _P.generate_isomer_decoys!(prec, frag, vm, masses)
    dec = prec[prec.is_loc_decoy, :]
    @test nrow(dec) == 1                       # only the ambiguous peptidoform 1
    @test Int(dec.loc_base_prec_id[1]) == 1    # decoy derived from peptidoform 1
    let s = _sites(dec.mods[1]), b = _sites(prec.mods[1])
        countname(set) = Dict(n => count(x -> x[2] == n, set) for n in unique(x[2] for x in set))
        @test countname(s) == countname(b)     # per-mod-name composition preserved
        @test count(sp -> sp[1] in shad_all, s) == 1     # exactly one mod on a shadow
        # the moved mod is the phospho (ambiguous); the ox (Unimod:35) stays on M7
        @test (7, "Unimod:35") in s
    end
end

@testset "no decoy for unambiguous localization (C(c,n) == 1)" begin
    vm = VM[(p = r"[STY]", r = "Unimod:21")]
    masses = Dict("Unimod:21" => 79.966331f0)
    # single phosphosite: 1 S, 1 phospho -> C(1,1)=1 -> no decoy
    prec1 = DataFrame(sequence = ["AASAAK"], mods = ["(3,S,Unimod:21)"],
                      mz = Float32[400], prec_charge = UInt8[2], pair_id = UInt32[1])
    frag1 = DataFrame(precursor_idx = UInt32[1], annotation = ["b3"],
                      mz = Float32[300], intensities = Float32[1])
    p1, _ = _P.generate_isomer_decoys!(prec1, frag1, vm, masses)
    @test nrow(p1) == 1 && !any(p1.is_loc_decoy)

    # fully saturated: 2 S, 2 phospho -> C(2,2)=1 -> no decoy
    prec2 = DataFrame(sequence = ["ASASAK"], mods = ["(2,S,Unimod:21)(4,S,Unimod:21)"],
                      mz = Float32[500], prec_charge = UInt8[2], pair_id = UInt32[1])
    frag2 = DataFrame(precursor_idx = UInt32[1], annotation = ["b3"],
                      mz = Float32[300], intensities = Float32[1])
    p2, _ = _P.generate_isomer_decoys!(prec2, frag2, vm, masses)
    @test nrow(p2) == 1 && !any(p2.is_loc_decoy)

    # 2 S, 1 phospho -> C(2,1)=2 -> DOES get a decoy (control)
    prec3 = DataFrame(sequence = ["ASASAK"], mods = ["(2,S,Unimod:21)"],
                      mz = Float32[420], prec_charge = UInt8[2], pair_id = UInt32[1])
    frag3 = DataFrame(precursor_idx = UInt32[1], annotation = ["b3"],
                      mz = Float32[300], intensities = Float32[1])
    p3, _ = _P.generate_isomer_decoys!(prec3, frag3, vm, masses)
    @test nrow(p3) == 2 && count(p3.is_loc_decoy) == 1
end

@testset "unmodified peptidoform gets no decoy" begin
    seq = "SGGSGGMG"
    vm = VM[(p = r"[STY]", r = "Unimod:21")]
    prec = DataFrame(sequence = [seq], mods = [""], mz = Float32[400.0])
    frag = DataFrame(precursor_idx = UInt32[1], annotation = ["b2"], mz = Float32[200], intensities = Float32[1])
    prec, frag = _P.generate_isomer_decoys!(prec, frag, vm, Dict("Unimod:21"=>79.966331f0))
    @test nrow(prec) == 1 && all(.!prec.is_loc_decoy)
end

@testset "multi-mod: ambiguity requires c>n for at least ONE mod type" begin
    # phospho (STY) + oxidation (M). Decoy iff C(c_p,n_p)*C(c_o,n_o) > 1.
    vm = VM[(p = r"[STY]", r = "Unimod:21"), (p = r"M", r = "Unimod:35")]
    masses = Dict("Unimod:21" => 79.966331f0, "Unimod:35" => 15.994915f0)
    mk(seq, mods, mz) = (DataFrame(sequence=[seq], mods=[mods], mz=Float32[mz], prec_charge=UInt8[2], pair_id=UInt32[1]),
                         DataFrame(precursor_idx=UInt32[1], annotation=["b3"], mz=Float32[300], intensities=Float32[1]))
    # 1 phospho site (S3) + 1 ox site (M5), 1 phospho + 1 ox -> C=1*1=1 -> NO decoy (the caveat)
    p, f = mk("AASAMK", "(3,S,Unimod:21)(5,M,Unimod:35)", 500); r, _ = _P.generate_isomer_decoys!(p, f, vm, masses)
    @test nrow(r) == 1 && !any(r.is_loc_decoy)
    # 2 phospho sites (S3,S4) + 1 ox site, 1 phospho + 1 ox -> phospho ambiguous -> decoy
    p, f = mk("AASSMK", "(3,S,Unimod:21)(5,M,Unimod:35)", 520); r, _ = _P.generate_isomer_decoys!(p, f, vm, masses)
    @test nrow(r) == 2 && count(r.is_loc_decoy) == 1
    # 1 phospho site + 2 ox sites (M5,M6), 1 phospho + 1 ox -> OX ambiguous -> decoy
    p, f = mk("AASAMM", "(3,S,Unimod:21)(5,M,Unimod:35)", 540); r, _ = _P.generate_isomer_decoys!(p, f, vm, masses)
    @test nrow(r) == 2 && count(r.is_loc_decoy) == 1
end

@testset "loc-decoy pair_id preserves real target/decoy pairing" begin
    seq = "VLSAGSPESIK"
    vm = VM[(p = r"[STY]", r = "Unimod:21")]
    prec = DataFrame(sequence = [seq, seq], mods = ["(3,S,Unimod:21)", "(6,S,Unimod:21)"],
                     mz = Float32[600, 600], prec_charge = UInt8[2, 2],
                     is_decoy = [false, true], pair_id = UInt32[10, 10])
    frag = DataFrame(precursor_idx = UInt32[1, 2], annotation = ["b4", "b4"],
                     mz = Float32[500, 500], intensities = Float32[1, 1])
    prec, frag = _P.generate_isomer_decoys!(prec, frag, vm, Dict("Unimod:21"=>79.966331f0))
    @test nrow(prec) == 4
    @test all(prec.pair_id[3:4] .> UInt32(10)) && allunique(prec.pair_id[3:4])
    _P.add_pair_indices!(prec)
    @test prec.partner_precursor_idx[1] == 2 && prec.partner_precursor_idx[2] == 1
    @test ismissing(prec.partner_precursor_idx[3]) && ismissing(prec.partner_precursor_idx[4])
end
