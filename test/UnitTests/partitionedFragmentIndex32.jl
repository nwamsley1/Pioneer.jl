using Pioneer: SoAFragBins, LocalFragment, LocalFragment32, Counter, LocalPartition, LocalPartition32,
    LocalPartitionedFragmentIndex, LocalPartitionedFragmentIndex32, FragIndexBin, MassErrorModel,
    AbstractLocalPartition, AbstractLocalPartitionedFragmentIndex, MAX_LOCAL_PRECS,
    local_id_type, max_local_precs, local_fragment_type, local_partition_type, local_index_type,
    getFragBins, getRTBins, getFragments, getSkipHints, getPartitions, getPartition, getNPartitions,
    get_partition_range, getPrecID, getScore, _build_local_partition, _score_partition_hinted!,
    SimpleFrag, EmitToBuffer, CountFilter, emit_candidates!, FragIndexScratch, prepare!, scratch_counters

# UInt32-ID variant of the partitioned fragment index. The UInt16 variant is covered by
# partitionedFragmentIndex.jl / buildPartitionedIndex.jl; here each check compares the two.

"Global precursor ID => score for every precursor the counter scored."
function global_scores(lc::Counter{I, UInt8}, l2g::Vector{UInt32}) where {I}
    out = Dict{UInt32, UInt8}()
    @inbounds for i in 1:(lc.size - 1)
        lid = lc.ids[i]
        s = lc.counts[lid]
        s > 0 && (out[l2g[lid]] = s)
    end
    out
end

"`n_precs` precursors, `n_frags` fragments each (rank bitmask 1 << (r-1)), fragment m/z spread over 150-1500."
function synthetic_frags(n_precs::Int, n_frags::Int; seed::Int = 7)
    rng = Random.MersenneTwister(seed)
    frags = SimpleFrag{Float32}[]
    for pid in 1:n_precs, r in 1:n_frags
        push!(frags, SimpleFrag{Float32}(Float32(150 + 1350 * rand(rng)), UInt32(pid), 500.0f0,
                                         Float32(10 * rand(rng)), UInt8(0), UInt8(1) << UInt8(r - 1)))
    end
    frags
end

@testset "PartitionedFragmentIndex UInt32 local IDs" begin

    @testset "types and ID-type mapping" begin
        @test sizeof(LocalFragment) == 4 && sizeof(LocalFragment32) == 8
        f = LocalFragment32(UInt32(70_000), UInt8(5))
        @test getPrecID(f) === UInt32(70_000) && getScore(f) === UInt8(5)
        @test local_id_type(LocalFragment) === UInt16 && local_id_type(LocalFragment32) === UInt32
        @test local_fragment_type(UInt16) === LocalFragment && local_fragment_type(UInt32) === LocalFragment32
        @test local_partition_type(UInt16) === LocalPartition && local_partition_type(UInt32) === LocalPartition32
        @test local_index_type(UInt16) === LocalPartitionedFragmentIndex && local_index_type(UInt32) === LocalPartitionedFragmentIndex32
        @test max_local_precs(UInt16) == MAX_LOCAL_PRECS == 65535
        @test max_local_precs(UInt32) > 4_000_000_000
        @test LocalPartition{Float32} <: AbstractLocalPartition{Float32} && LocalPartition32{Float32} <: AbstractLocalPartition{Float32}
        @test LocalPartitionedFragmentIndex32{Float32} <: AbstractLocalPartitionedFragmentIndex{Float32}
        @test local_id_type(LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32, Float32}[], 0)) === UInt16
        @test local_id_type(LocalPartitionedFragmentIndex32{Float32}(LocalPartition32{Float32}[], Tuple{Float32, Float32}[], 0)) === UInt32
    end

    @testset "UInt16 and UInt32 partitions score identically" begin
        n = 400
        frags = synthetic_frags(n, 6)
        l2g = UInt32.(1001:1000 + n)                 # non-trivial local -> global map
        p16 = _build_local_partition(copy(frags), l2g, UInt16(n), 0.0f0, 2.0f0, 3.0f0)
        p32 = _build_local_partition(copy(frags), l2g, UInt32(n), 0.0f0, 2.0f0, 3.0f0)
        @test p16 isa LocalPartition{Float32} && p32 isa LocalPartition32{Float32}
        @test eltype(getFragments(p32)) === LocalFragment32 && p32.n_local_precs === UInt32(n)
        # identical layout: bins, RT bins, hints, and fragment (id, score) sequence
        @test getFragBins(p16).lows == getFragBins(p32).lows && getFragBins(p16).highs == getFragBins(p32).highs
        @test getRTBins(p16) == getRTBins(p32) && getSkipHints(p16) == getSkipHints(p32)
        @test [Int(getPrecID(f)) for f in getFragments(p16)] == [Int(getPrecID(f)) for f in getFragments(p32)]
        @test getScore.(getFragments(p16)) == getScore.(getFragments(p32))

        rng = Random.MersenneTwister(11)
        mem = MassErrorModel(0.0f0, (10.0f0, 10.0f0))
        for trial in 1:20
            masses = Union{Missing, Float32}[sort(Float32.(150 .+ 1350 .* rand(rng, 200)))...]
            ints = Union{Missing, Float32}[ones(Float32, 200)...]
            irt_lo = Float32(10 * rand(rng) - 2); irt_hi = irt_lo + 4.0f0
            c16 = Counter(UInt16, UInt8, n + 1); c32 = Counter(UInt32, UInt8, n + 1)
            _score_partition_hinted!(c16, p16, irt_lo, irt_hi, masses, ints, mem)
            _score_partition_hinted!(c32, p32, irt_lo, irt_hi, masses, ints, mem)
            @test global_scores(c16, l2g) == global_scores(c32, l2g)
        end
    end

    @testset "local IDs above 65,535" begin
        n = 70_000
        frags = synthetic_frags(n, 1)                 # one fragment each: rank bit 0
        l2g = UInt32.(1:n) .+ UInt32(5_000_000)
        p32 = _build_local_partition(copy(frags), l2g, UInt32(n), 0.0f0, 2.0f0, 1.0f6)
        ids = [getPrecID(f) for f in getFragments(p32)]
        @test maximum(ids) == UInt32(n) && sort(ids) == UInt32.(1:n)       # every local ID kept, none wrapped
        # the UInt16 builder refuses instead of wrapping
        @test_throws InexactError _build_local_partition(copy(frags), l2g, UInt16(MAX_LOCAL_PRECS), 0.0f0, 2.0f0, 1.0f6)
        # query the fragment m/z of a few precursors above 65,535 and check they come back with the right global ID
        targets = [65_536, 66_000, 69_999, 70_000]
        tmz = Float32[frags[t].mz for t in targets]
        o = sortperm(tmz)
        masses = Union{Missing, Float32}[tmz[o]...]; ints = Union{Missing, Float32}[ones(Float32, length(o))...]
        c = Counter(UInt32, UInt8, n + 1)
        _score_partition_hinted!(c, p32, -1.0f0, 11.0f0, masses, ints, MassErrorModel(0.0f0, (0.1f0, 0.1f0)))
        g = global_scores(c, l2g)
        for t in targets
            @test get(g, UInt32(t) + UInt32(5_000_000), 0x00) == 0x01
        end
    end

    @testset "get_partition_range on a UInt32 index" begin
        parts = [LocalPartition32{Float32}(SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[]),
                     FragIndexBin{Float32}[], LocalFragment32[], UInt32[], UInt32(0), UInt16[]) for _ in 1:3]
        pfi = LocalPartitionedFragmentIndex32{Float32}(parts, [(400.0f0, 410.0f0), (410.0f0, 420.0f0), (420.0f0, 430.0f0)], 3)
        @test get_partition_range(pfi, 405.0f0, 415.0f0) == (1, 2)
        @test get_partition_range(pfi, 421.0f0, 425.0f0) == (3, 3)
        @test getNPartitions(pfi) == 3 && getPartition(pfi, 2) === parts[2]
    end

    @testset "emit_candidates! with a UInt32 counter" begin
        n = 70_000
        c = Counter(UInt32, UInt8, n + 1)
        Pioneer.or!(c, UInt32(3), UInt8(0x07))          # 3 fragments
        Pioneer.or!(c, UInt32(69_999), UInt8(0x0f))     # 4 fragments
        Pioneer.or!(c, UInt32(12), UInt8(0x01))         # 1 fragment (filtered)
        l2g = UInt32.(1:n) .* UInt32(2)
        si = Int32[0 for _ in 1:4]; pid = UInt32[0 for _ in 1:4]
        wp = emit_candidates!(EmitToBuffer(CountFilter(UInt8(3)), nothing), c, l2g, 5, 0.0f0, 0.0f0,
                              0.0f0, 0.0f0, Float32[], 1, si, pid, 0)
        @test wp == 2 && sort(pid[1:wp]) == UInt32[6, 139_998] && all(==(Int32(5)), si[1:wp])
    end

    @testset "FragIndexScratch UInt32 counters" begin
        s = FragIndexScratch(2)
        prepare!(s; n_threads = 2, est_per_thread = 10, counter_size = 70_001, int_buf_size = 0, id_type = UInt32)
        @test scratch_counters(s, UInt32) isa Vector{Counter{UInt32, UInt8}}
        @test all(c -> length(c.counts) >= 70_001, scratch_counters(s, UInt32))
        @test all(c -> length(c.counts) == 1, scratch_counters(s, UInt16))   # UInt16 counters untouched
        prepare!(s; n_threads = 2, est_per_thread = 10, counter_size = 100, int_buf_size = 0)
        @test all(c -> length(c.counts) >= 100, scratch_counters(s, UInt16))
    end
end
