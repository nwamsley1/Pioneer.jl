using Test
using Arrow
using DataFrames
using Random

using Pioneer
using Pioneer: BasicMassSpecData, getMsOrder, getCycleIdx, getCenterMz,
               map_any_hit_to_center!, filter_to_center_bin!, ZTGeometry

const MSTEP = 1.0221f0
const MSTART = 400.0f0
const MBINS = 24

"""Synthetic ZT arrow: `n_cycles` cycles, each an MS1 head plus `MBINS` contiguous MS2 bins."""
function write_map_arrow(path::String; n_cycles::Int = 6)
    per = MBINS + 1
    n = n_cycles * per
    ms_order  = Vector{UInt8}(undef, n)
    center_mz = Vector{Union{Missing,Float32}}(undef, n)
    isol_w    = Vector{Union{Missing,Float32}}(undef, n)
    cycle_idx = Vector{UInt32}(undef, n)
    for c in 1:n_cycles, b in 0:MBINS
        i = (c - 1) * per + b + 1
        cycle_idx[i] = UInt32(c)
        if b == 0
            ms_order[i] = UInt8(1); center_mz[i] = missing; isol_w[i] = missing
        else
            ms_order[i] = UInt8(2)
            center_mz[i] = MSTART + Float32(b - 1) * MSTEP
            isol_w[i] = MSTEP
        end
    end
    mz  = [Vector{Union{Missing,Float32}}(Float32[100f0]) for _ in 1:n]
    int = [Vector{Union{Missing,Float32}}(Float32[10f0])  for _ in 1:n]
    Arrow.write(path, DataFrame(
        mz_array = mz, intensity_array = int,
        scanHeader = fill("s", n), scanNumber = Int32.(1:n),
        basePeakMz = fill(100f0, n), basePeakIntensity = fill(1000f0, n),
        retentionTime = Float32.(range(0, length = n, step = 0.01f0)),
        lowMz = fill(100f0, n), highMz = fill(1000f0, n), TIC = fill(5f3, n),
        centerMz = center_mz, isolationWidthMz = isol_w,
        msOrder = ms_order, packetType = Int32.(fill(0, n)), cycle_idx = cycle_idx))
    return BasicMassSpecData(path)
end

"""
Independent serial reference: the literal Dict/Set implementation, finding each precursor's bin
by scanning +/-halfbins neighbours. Deliberately the slow obvious way — shares no logic with the
counting-sort implementation under test.
"""
function reference_map(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, pmz, halfbins)
    nspec = length(spectra)
    center_sets = Dict{Int, Set{UInt32}}()
    seen = Set{Tuple{UInt32, UInt32}}()
    for si in all_scan_idxs
        rng = scan_to_prec_idx[si]
        ismissing(rng) && continue
        cyc = getCycleIdx(spectra, si)
        for r in rng
            p = precursors_passed[r]
            key = (UInt32(cyc), p)
            (key in seen) && continue
            push!(seen, key)
            best = 0; bestd = Inf32
            for d in -halfbins:halfbins
                sj = si + d
                (sj < 1 || sj > nspec) && continue
                getMsOrder(spectra, sj) == 2 || continue
                getCycleIdx(spectra, sj) == cyc || continue
                cm = getCenterMz(spectra, sj)
                ismissing(cm) && continue
                dd = abs(Float32(cm) - Float32(pmz[p]))
                if dd < bestd; bestd = dd; best = sj; end
            end
            best == 0 && continue
            push!(get!(center_sets, best, Set{UInt32}()), p)
        end
    end
    return center_sets
end

"""Per-scan candidate sets after an in-place re-anchoring."""
function scan_sets_map(scan_to_prec_idx, precursors_passed, nspec)
    out = Dict{Int, Set{UInt32}}()
    for si in 1:nspec
        rng = scan_to_prec_idx[si]
        ismissing(rng) && continue
        isempty(rng) && continue
        out[si] = Set(precursors_passed[r] for r in rng)
    end
    return out
end

"""Emit each precursor into every MS2 bin within `tol` Da of its m/z — the wide-emit input."""
function wide_emissions(spectra, all_scan_idxs, pmz, tol::Float32)
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    fill!(s2p, missing)
    passed = UInt32[]
    for si in all_scan_idxs
        cm = getCenterMz(spectra, si)
        ismissing(cm) && continue
        start = length(passed) + 1
        for p in 1:length(pmz)
            abs(pmz[p] - Float32(cm)) <= tol && push!(passed, UInt32(p))
        end
        s2p[si] = length(passed) >= start ? (start:length(passed)) : missing
    end
    return s2p, passed
end

@testset "ZT map_any_hit_to_center! matches a serial reference" begin
    spectra = write_map_arrow(joinpath(mktempdir(), "m.arrow"); n_cycles = 6)
    nspec = length(spectra)
    all_scan_idxs = [i for i in 1:nspec if getMsOrder(spectra, i) == 2]
    geom = ZTGeometry(MSTEP, MSTEP, Int32(MBINS), Int32(6), 6.5f0)

    # precursors spread across the ramp, deliberately off bin centres by varying amounts
    rng = MersenneTwister(4242)
    pmz = Float32[MSTART + (MBINS - 1) * MSTEP * rand(rng) for _ in 1:60]

    @testset "tol = $tol" for tol in (0.51f0, 1.0f0, 2.0f0, 4.0f0)
        s2p, passed = wide_emissions(spectra, all_scan_idxs, pmz, tol)
        halfbins = ceil(Int, tol / MSTEP) + 1
        want = reference_map(copy(s2p), copy(passed), spectra, all_scan_idxs, pmz, halfbins)

        new_passed = map_any_hit_to_center!(s2p, passed, spectra, all_scan_idxs, pmz, geom)
        got = scan_sets_map(s2p, new_passed, nspec)

        # 1. identical per-scan candidate sets
        @test got == want

        # 2. every emitted (scan, precursor) landed on the precursor's OWN bin
        for (si, ps) in got, p in ps
            @test abs(pmz[p] - Float32(getCenterMz(spectra, si))) <= MSTEP / 2 + 1f-3
        end

        # 3. sorted and duplicate-free within each scan
        for si in 1:nspec
            r = s2p[si]
            ismissing(r) && continue
            v = new_passed[r]
            @test issorted(v)
            @test length(unique(v)) == length(v)
        end

        # 4. ranges tile new_passed exactly, in scan order
        off = 0
        for si in 1:nspec
            r = s2p[si]
            ismissing(r) && continue
            @test first(r) == off + 1
            off = last(r)
        end
        @test off == length(new_passed)
    end

    @testset "narrow tol reduces to filter_to_center_bin!" begin
        # With emissions confined to a precursor's own bin, re-anchoring must be a no-op
        # relative to the center-bin filter.
        s2p_a, passed_a = wide_emissions(spectra, all_scan_idxs, pmz, MSTEP / 2)
        s2p_b, passed_b = wide_emissions(spectra, all_scan_idxs, pmz, MSTEP / 2)
        out_a = map_any_hit_to_center!(s2p_a, passed_a, spectra, all_scan_idxs, pmz, geom)
        out_b = filter_to_center_bin!(s2p_b, passed_b, spectra, all_scan_idxs, pmz)
        @test scan_sets_map(s2p_a, out_a, nspec) == scan_sets_map(s2p_b, out_b, nspec)
    end

    @testset "determinism" begin
        a2p, ap = wide_emissions(spectra, all_scan_idxs, pmz, 2.0f0)
        b2p, bp = wide_emissions(spectra, all_scan_idxs, pmz, 2.0f0)
        ra = map_any_hit_to_center!(a2p, ap, spectra, all_scan_idxs, pmz, geom)
        rb = map_any_hit_to_center!(b2p, bp, spectra, all_scan_idxs, pmz, geom)
        @test ra == rb
        @test isequal(a2p, b2p)
    end

    @testset "empty input" begin
        s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, nspec); fill!(s2p, missing)
        out = map_any_hit_to_center!(s2p, UInt32[], spectra, all_scan_idxs, pmz, geom)
        @test isempty(out)
    end
end
