using Test
using Arrow
using DataFrames
using Random

using Pioneer: BasicMassSpecData, getMsOrder, getCycleIdx, expand_to_metascans!, _sort_dedup!

"""
Write a synthetic scanning-quad (ZT) arrow: `n_cycles` cycles, each one MS1 head scan followed
by `bins` contiguous MS2 bins whose centerMz ramps by `step`.
"""
function write_zt_ms_arrow(path::String; n_cycles::Int = 5, bins::Int = 20, step::Float64 = 1.0221)
    per_cycle = bins + 1
    n = n_cycles * per_cycle

    ms_order  = Vector{UInt8}(undef, n)
    center_mz = Vector{Union{Missing, Float32}}(undef, n)
    isol_w    = Vector{Union{Missing, Float32}}(undef, n)
    cycle_idx = Vector{UInt32}(undef, n)

    for c in 1:n_cycles, b in 0:bins
        i = (c - 1) * per_cycle + b + 1
        cycle_idx[i] = UInt32(c)
        if b == 0                       # MS1 head of the cycle
            ms_order[i]  = UInt8(1)
            center_mz[i] = missing
            isol_w[i]    = missing
        else
            ms_order[i]  = UInt8(2)
            center_mz[i] = Float32(400.0 + (b - 1) * step)
            isol_w[i]    = Float32(step)
        end
    end

    mz  = [Vector{Union{Missing,Float32}}(Float32[100.0, 200.0]) for _ in 1:n]
    int = [Vector{Union{Missing,Float32}}(Float32[10.0, 20.0])   for _ in 1:n]

    Arrow.write(path, DataFrame(
        mz_array = mz, intensity_array = int,
        scanHeader = ["s$(i)" for i in 1:n], scanNumber = Int32.(1:n),
        basePeakMz = Float32.(fill(100, n)), basePeakIntensity = Float32.(fill(1000, n)),
        retentionTime = Float32.(range(0, length = n, step = 0.01)),
        lowMz = Float32.(fill(100, n)), highMz = Float32.(fill(1000, n)),
        TIC = Float32.(fill(5000, n)),
        centerMz = center_mz, isolationWidthMz = isol_w,
        msOrder = ms_order, packetType = Int32.(fill(0, n)),
        cycle_idx = cycle_idx,
    ))
    return path
end

"""
Independent serial reference: the straightforward `Set`-per-scan union over ±k same-cycle MS2
neighbours. Deliberately written the slow, obvious way so it shares no logic with the threaded
implementation under test. Returns scan_idx => Set of candidate precursor ids.
"""
function reference_expand(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, k)
    out = Dict{Int, Set{UInt32}}()
    nspec = length(spectra)
    for si in all_scan_idxs
        ci = getCycleIdx(spectra, si)
        s = Set{UInt32}()
        for j in -k:k
            sj = si + j
            (sj < 1 || sj > nspec) && continue
            getMsOrder(spectra, sj) == 2 || continue
            getCycleIdx(spectra, sj) == ci || continue
            rng = scan_to_prec_idx[sj]
            ismissing(rng) && continue
            for r in rng
                push!(s, precursors_passed[r])
            end
        end
        out[si] = s
    end
    return out
end

"""Read back per-scan candidate sets after an in-place expansion."""
function scan_sets(scan_to_prec_idx, precursors_passed, all_scan_idxs)
    out = Dict{Int, Set{UInt32}}()
    for si in all_scan_idxs
        rng = scan_to_prec_idx[si]
        out[si] = ismissing(rng) ? Set{UInt32}() : Set(precursors_passed[r] for r in rng)
    end
    return out
end

"""Random candidacy over the MS2 scans; ~15% of scans get no candidates."""
function random_candidacy(spectra, all_scan_idxs, rng_seed::Int)
    rng = MersenneTwister(rng_seed)
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    fill!(scan_to_prec_idx, missing)
    precursors_passed = UInt32[]
    for si in all_scan_idxs
        if rand(rng) < 0.15
            scan_to_prec_idx[si] = missing
            continue
        end
        m = rand(rng, 1:8)
        start = length(precursors_passed) + 1
        for _ in 1:m
            push!(precursors_passed, UInt32(rand(rng, 1:50)))
        end
        scan_to_prec_idx[si] = start:length(precursors_passed)
    end
    return scan_to_prec_idx, precursors_passed
end

@testset "ZT expand_to_metascans! matches a serial reference" begin
    dir = mktempdir()
    spectra = BasicMassSpecData(write_zt_ms_arrow(joinpath(dir, "zt.arrow");
                                                 n_cycles = 5, bins = 20))
    all_scan_idxs = [i for i in 1:length(spectra) if getMsOrder(spectra, i) == 2]
    @test length(all_scan_idxs) == 100

    @testset "k = $k" for k in (1, 2, 6)
        s2p, pp = random_candidacy(spectra, all_scan_idxs, 42 + k)
        want = reference_expand(copy(s2p), copy(pp), spectra, all_scan_idxs, k)

        new_pp = expand_to_metascans!(s2p, pp, spectra, all_scan_idxs, k)
        got = scan_sets(s2p, new_pp, all_scan_idxs)

        # 1. same candidate set on every scan
        @test got == want

        # 2. sorted and duplicate-free within each scan
        for si in all_scan_idxs
            rng = s2p[si]
            ismissing(rng) && continue
            v = new_pp[rng]
            @test issorted(v)
            @test length(unique(v)) == length(v)
        end

        # 3. ranges tile new_pp exactly, in scan order, with no gaps or overlap
        expected_len = sum(length(want[si]) for si in all_scan_idxs)
        @test length(new_pp) == expected_len
        off = 0
        for si in all_scan_idxs
            rng = s2p[si]
            if ismissing(rng)
                @test isempty(want[si])
            else
                @test first(rng) == off + 1
                off = last(rng)
            end
        end
        @test off == length(new_pp)

        # 4. expansion never crosses a cycle boundary
        for si in all_scan_idxs
            rng = s2p[si]
            ismissing(rng) && continue
            ci = getCycleIdx(spectra, si)
            for p in new_pp[rng]
                @test any(-k:k) do j
                    sj = si + j
                    1 <= sj <= length(spectra) &&
                        getMsOrder(spectra, sj) == 2 &&
                        getCycleIdx(spectra, sj) == ci
                end
            end
        end
    end

    @testset "determinism" begin
        a2p, app = random_candidacy(spectra, all_scan_idxs, 7)
        b2p, bpp = random_candidacy(spectra, all_scan_idxs, 7)
        ra = expand_to_metascans!(a2p, app, spectra, all_scan_idxs, 3)
        rb = expand_to_metascans!(b2p, bpp, spectra, all_scan_idxs, 3)
        @test ra == rb
        @test isequal(a2p, b2p)
    end

    @testset "k <= 0 is a no-op" begin
        s2p, pp = random_candidacy(spectra, all_scan_idxs, 11)
        before = copy(s2p)
        out = expand_to_metascans!(s2p, pp, spectra, all_scan_idxs, 0)
        @test out === pp
        @test isequal(s2p, before)
    end
end

@testset "_sort_dedup!" begin
    v = UInt32[5, 1, 5, 3, 1, 1]
    m = _sort_dedup!(v)
    @test m == 3
    @test v[1:m] == UInt32[1, 3, 5]

    empty_v = UInt32[]
    @test _sort_dedup!(empty_v) == 0

    single = UInt32[9]
    @test _sort_dedup!(single) == 1
    @test single[1] == UInt32(9)

    same = UInt32[4, 4, 4, 4]
    @test _sort_dedup!(same) == 1
    @test same[1] == UInt32(4)
end
