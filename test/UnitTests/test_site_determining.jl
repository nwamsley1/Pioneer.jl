# Unit tests for site-determining retention + shadow decoy sites.
# See dev_docs/prosit_ptm_integration/isomer_decoy_shadow_sites.md.
using Pioneer, Test, DataFrames
using Polynomials: ImmutablePolynomial

const _P = Pioneer
const VM = NamedTuple{(:p, :r), Tuple{Regex, String}}
_matchvm(seq, vm) = _P.matchVarMods(String(seq), vm)

@testset "acceptor_positions + real_sites_by_mod (two mods)" begin
    seq = "SGGSGGMG"                 # S1 G2 G3 S4 G5 G6 M7 G8
    two = VM[(p = r"[STY]", r = "Unimod:21"), (p = r"M", r = "Unimod:35")]
    @test _P.acceptor_positions(_matchvm(seq, VM[(p=r"[STY]", r="Unimod:21")])) == UInt8[1, 4]
    rb = _P.real_sites_by_mod(seq, two)
    @test rb["Unimod:21"] == [1, 4]
    @test rb["Unimod:35"] == [7]
end

@testset "seq_seed is deterministic + sequence-specific" begin
    @test _P.seq_seed("SGGSGGMG") == _P.seq_seed("SGGSGGMG")
    @test _P.seq_seed("SGGSGGMG") != _P.seq_seed("SGGSGGMA")
end

@testset "shadow_decoy_sites: count-matched, distinct, non-acceptor, deterministic" begin
    seq = "SGGSGGMG"
    rb = Dict("Unimod:21" => [1, 4], "Unimod:35" => [7])
    s1 = _P.shadow_decoy_sites(seq, rb; seed = _P.seq_seed(seq))
    s2 = _P.shadow_decoy_sites(seq, rb; seed = _P.seq_seed(seq))
    @test s1 == s2                                   # deterministic
    @test length(s1["Unimod:21"]) == 2               # count-matched to 2 real S
    @test length(s1["Unimod:35"]) == 1               # count-matched to 1 real M
    allshadow = vcat(values(s1)...)
    @test allunique(allshadow)                       # mutually distinct
    @test all(p -> !(p in (1, 4, 7)), allshadow)     # none on a real acceptor
    @test all(p -> 1 <= p <= length(seq), allshadow)
end

@testset "compute_gap_sites: real-only vs union(real, shadow)" begin
    seq = "SGGSGGMG"
    two = VM[(p = r"[STY]", r = "Unimod:21"), (p = r"M", r = "Unimod:35")]
    @test _P.compute_gap_sites(seq, two, false) == UInt8[1, 4, 7]      # real acceptors only
    g = _P.compute_gap_sites(seq, two, true)                          # union with shadow
    shadow = vcat(values(_P.shadow_decoy_sites(seq, Dict("Unimod:21"=>[1,4],"Unimod:35"=>[7]);
                                               seed = _P.seq_seed(seq)))...)
    @test g == UInt8.(sort(unique(vcat([1,4,7], shadow))))
    @test issubset(UInt8[1,4,7], g)                                    # real sites present
    @test all(UInt8(p) in g for p in shadow)                          # shadow sites present
end

@testset "no-op: no acceptor residues" begin
    seq = "LLEELLK"
    two = VM[(p = r"[STY]", r = "Unimod:21")]
    @test _P.compute_gap_sites(seq, two, true) == UInt8[]
    @test _P.gap_cover_indices(UInt8[], 7, Char[], UInt8[], Int[], 5) == Int[]
end

# --- gap_cover_indices (unchanged logic) on VLSAGSPESIK (L=11) ---
function _frag_layout()
    base = Char[fill('b', 10); fill('y', 10)]
    ord  = UInt8[collect(1:10); collect(1:10)]
    return base, ord
end

@testset "gap_cover_indices picks the top in-gap ion outside top-N" begin
    base, ord = _frag_layout(); L = 11
    total = zeros(Float32, 20)
    for (r, t) in zip([1,2,9,10,11,12,19,20], 100.0f0:-1f0:93f0); total[r] = t; end
    total[4] = 50f0; total[7] = 45f0
    total[3]=40f0; total[5]=39f0; total[16]=38f0; total[17]=37f0; total[18]=36f0
    total[6]=30f0; total[8]=29f0; total[13]=28f0; total[14]=27f0; total[15]=26f0
    block = sortperm(total, rev = true)
    @test _P.gap_cover_indices(UInt8[3, 6, 9], L, base, ord, block, 8) == [4, 7]
end

@testset "filter_fragments! wires gap retention (VLSAGSPESIK, real ctx)" begin
    fb = _P.FragBoundModel(ImmutablePolynomial(Float32[0]), ImmutablePolynomial(Float32[1f6]))
    precdf = DataFrame(sequence = ["VLSAGSPESIK"], mods = [""], mz = Float32[600.0])
    build_ctx(gs) = _P.build_agnostic_frag_filter_ctx(
        fb, precdf, Dict{String,Int8}();
        y_start = UInt8(3), b_start = UInt8(2), include_p = false, include_isotope = false,
        include_immonium = false, include_internal = false, include_neutral_diff = true,
        max_frag_charge = UInt8(3), max_frag_rank = UInt8(5),
        length_to_frag_count_multiple = Float32(1000), min_frag_intensity = 0.0f0, gap_sites = gs)
    anns = String[]; mzs = Float32[]
    for i in 1:10; push!(anns, "b$i"); push!(mzs, Float32(150 + i*7)); end
    for j in 1:10; push!(anns, "y$j"); push!(mzs, Float32(500 + j*7)); end
    intmap = Dict("b2"=>1f6,"b9"=>1f6,"b10"=>1f6,"y9"=>1f6,"y10"=>1f6,
                  "b4"=>500f0,"b7"=>450f0,"b3"=>100f0,"b5"=>90f0,"b6"=>80f0,"b8"=>70f0,
                  "y3"=>60f0,"y4"=>50f0,"y5"=>40f0,"y6"=>30f0,"y7"=>20f0,"y8"=>10f0,
                  "b1"=>5f0,"y1"=>5f0,"y2"=>5f0)
    ints = Float32[intmap[a] for a in anns]
    mkdf() = DataFrame(annotation = copy(anns), mz = copy(mzs),
                       intensities = copy(ints), precursor_idx = fill(UInt32(1), 20))
    model = _P.InstrumentAgnosticModel("test")
    on  = Set(_P.filter_fragments!(mkdf(), model, build_ctx([UInt8[3,6,9]])).annotation)
    off = Set(_P.filter_fragments!(mkdf(), model, build_ctx(Vector{UInt8}[])).annotation)
    @test !("b4" in off) && !("b7" in off)
    @test ("b4" in on) && ("b7" in on)
    @test length(on) == length(off) + 2
end
