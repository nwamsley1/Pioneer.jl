# Unit tests for site-determining fragment retention (phospho localization).
# Hand-verified on the worked example VLSAGSPESIK (S at 3,6,9). See
# dev_docs/prosit_ptm_integration/site_determining_fragment_retention.md.
using Pioneer, Test, DataFrames
using Polynomials: ImmutablePolynomial

const _P = Pioneer
_vm(seq) = _P.matchVarMods(String(seq),
    NamedTuple{(:p, :r), Tuple{Regex, String}}[(p = r"[STY]", r = "Unimod:21")])

@testset "acceptor / decoy / gap sites (VLSAGSPESIK)" begin
    seq = "VLSAGSPESIK"                      # V L S A G S P E S I K ; S at 3,6,9
    vm = _vm(seq)
    @test _P.acceptor_positions(vm) == UInt8[3, 6, 9]
    # neighbours of 3->{2,4}, 6->{5,7}, 9->{8,10}; all non-STY
    @test _P.decoy_neighbor_positions(seq, UInt8[3, 6, 9]) == UInt8[2, 4, 5, 7, 8, 10]
    @test _P.compute_gap_sites(seq, vm, false) == UInt8[3, 6, 9]
    @test _P.compute_gap_sites(seq, vm, true)  == UInt8[2, 3, 4, 5, 6, 7, 8, 9, 10]
end

@testset "no-op: no acceptor residues" begin
    seq = "LLEELLK"
    vm = _vm(seq)
    @test _P.acceptor_positions(vm) == UInt8[]
    @test _P.compute_gap_sites(seq, vm, true) == UInt8[]
    # empty gap_sites -> no extra retained
    @test _P.gap_cover_indices(UInt8[], 7, Char[], UInt8[], Int[], 5) == Int[]
end

# Fragment layout for VLSAGSPESIK (L=11): rows 1..10 = b1..b10, rows 11..20 = y1..y10.
# cleavage: b_i -> i ; y_m -> 11 - m.
# gaps for S={3,6,9}: (3,6) crossers c in {3,4,5} ; (6,9) crossers c in {6,7,8}.
function _frag_layout()
    base = Char[fill('b', 10); fill('y', 10)]
    ord  = UInt8[collect(1:10); collect(1:10)]
    return base, ord
end

@testset "gap_cover_indices picks the top in-gap ion outside top-N" begin
    base, ord = _frag_layout()
    L = 11
    # Assign intensities so the 8 non-crossers dominate (top-N), and among gap
    # crossers b4 (row 4) is the strongest of gap(3,6) and b7 (row 7) the
    # strongest of gap(6,9). topn_cut = 8 -> both gap winners lie outside top-N.
    total = zeros(Float32, 20)
    noncross = [1, 2, 9, 10, 11, 12, 19, 20]        # c in {1,2,9,10} on both series
    for (r, t) in zip(noncross, 100.0f0:-1f0:93f0); total[r] = t; end
    total[4] = 50f0   # b4  (gap 3,6 winner)
    total[7] = 45f0   # b7  (gap 6,9 winner)
    total[3] = 40f0; total[5] = 39f0; total[16] = 38f0; total[17] = 37f0; total[18] = 36f0
    total[6] = 30f0; total[8] = 29f0; total[13] = 28f0; total[14] = 27f0; total[15] = 26f0

    block = sortperm(total, rev = true)             # descending-intensity order
    extra = _P.gap_cover_indices(UInt8[3, 6, 9], L, base, ord, block, 8)
    @test extra == [4, 7]                            # b4 then b7

    # sanity: b4 has cleavage 4 (in gap 3,6), b7 cleavage 7 (in gap 6,9)
    @test ord[4] == 4 && base[4] == 'b'
    @test ord[7] == 7 && base[7] == 'b'
end

@testset "per-peptidoform: different forms keep different in-gap ions" begin
    base, ord = _frag_layout()
    L = 11
    # Form A: y7 (row 17, cleavage 11-7=4, in gap 3,6) is the strongest gap(3,6) crosser.
    tA = fill(1f0, 20); tA[[1,2,9,10,11,12,19,20]] .= 100f0
    tA[17] = 50f0                                   # y7 strongest in gap(3,6)
    tA[7]  = 45f0                                   # b7 strongest in gap(6,9)
    bA = sortperm(tA, rev = true)
    exA = _P.gap_cover_indices(UInt8[3,6,9], L, base, ord, bA, 8)
    @test 17 in exA                                 # form A keeps y7 for gap(3,6)

    # Form B: b4 (row 4) is the strongest gap(3,6) crosser.
    tB = fill(1f0, 20); tB[[1,2,9,10,11,12,19,20]] .= 100f0
    tB[4] = 50f0                                    # b4 strongest in gap(3,6)
    tB[7] = 45f0
    bB = sortperm(tB, rev = true)
    exB = _P.gap_cover_indices(UInt8[3,6,9], L, base, ord, bB, 8)
    @test 4 in exB                                  # form B keeps b4 for the same gap
    @test !(17 in exB)                              # and NOT y7
end

@testset "uncoverable gap contributes nothing (no crosser present)" begin
    base, ord = _frag_layout()
    L = 11
    # Only gap(3,6) crossers carry intensity; every gap(6,9) crosser is zero.
    total = zeros(Float32, 20)
    total[[1,2,9,10,11,12,19,20]] .= 100f0
    total[4] = 50f0                                 # gap(3,6) covered
    # rows 6,7,8,13,14,15 (gap 6,9 crossers) left at 0 -> excluded by a real
    # keep-filter; here we simulate by omitting them from `block`.
    block = [r for r in sortperm(total, rev = true) if !(r in (6,7,8,13,14,15))]
    extra = _P.gap_cover_indices(UInt8[3,6,9], L, base, ord, block, 8)
    @test extra == [4]                              # gap(6,9) uncoverable -> nothing
end

@testset "filter_fragments! wires gap retention (VLSAGSPESIK, real ctx)" begin
    fb = _P.FragBoundModel(ImmutablePolynomial(Float32[0]), ImmutablePolynomial(Float32[1f6]))
    precdf = DataFrame(sequence = ["VLSAGSPESIK"], mods = [""], mz = Float32[600.0])
    build_ctx(gs) = _P.build_agnostic_frag_filter_ctx(
        fb, precdf, Dict{String,Int8}();
        y_start = UInt8(3), b_start = UInt8(2), include_p = false, include_isotope = false,
        include_immonium = false, include_internal = false, include_neutral_diff = true,
        max_frag_charge = UInt8(3), max_frag_rank = UInt8(5),
        length_to_frag_count_multiple = Float32(1000), min_frag_intensity = 0.0f0,
        gap_sites = gs)

    # b1..b10, y1..y10; non-crossers (c in {1,2,9,10}) huge so they fill top-5;
    # b4 (gap 3,6) and b7 (gap 6,9) are the top crossers but below the top-5 cut.
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

    @test !("b4" in off) && !("b7" in off)          # top-5 alone drops the gap ions
    @test ("b4" in on) && ("b7" in on)              # retention keeps them
    @test issubset(off, on)                          # ON == OFF + gap ions
    @test length(on) == length(off) + 2
end
