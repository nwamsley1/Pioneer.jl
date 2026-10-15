using Test, Pioneer, Arrow, DataFrames
include(joinpath(@__DIR__, "..", "..", "reference", "library_storage", "CompactLibraryStorage.jl"))
const CLS = CompactLibraryStorage

function storage_test_index(::Type{I}, delta) where {I}
    low = 500f0
    high = reinterpret(Float32, reinterpret(UInt32, low) + UInt32(delta))
    bins = Pioneer.SoAFragBins(Float32[low], vcat(Float32[high], fill(Inf32, 7)), UInt32[1], UInt32[4])
    F = Pioneer.local_fragment_type(I)
    entries = F[F(I(id), mask) for (id, mask) in zip(1:4, UInt8[0x01, 0x80, 0xff, 0x00])]
    Q = Pioneer.local_partition_type(I){Float32}
    partition = Q(bins, [Pioneer.FragIndexBin(0f0, 10f0, UInt32(1), UInt32(1))],
        entries, UInt32[4, 3, 2, 1], I(4), UInt16[1])
    return Pioneer.local_index_type(I){Float32}([partition], [(500f0, 501f0)], 1)
end

function test_fragment_bits(actual, expected)
    @test actual.prec_id == expected.prec_id
    @test reinterpret(UInt32, actual.mz) == reinterpret(UInt32, expected.mz)
    if actual.intensity isa Float16
        @test reinterpret(UInt16, actual.intensity) == reinterpret(UInt16, expected.intensity)
    else
        @test map(x -> reinterpret(UInt32, x), actual.intensity) ==
            map(x -> reinterpret(UInt32, x), expected.intensity)
    end
    @test actual.packed_a == expected.packed_a
    @test actual.packed_b == expected.packed_b
end

@testset "Compact library storage reference" begin
    @testset "Index widths, exact deltas, and full masks" begin
        for I in (UInt16, UInt32), (delta, D) in ((0, UInt8), (255, UInt8),
            (256, UInt16), (65535, UInt16), (65536, UInt32))
            original = storage_test_index(I, delta)
            stored = CLS.encode_index(original)
            @test eltype(stored.partitions[1].high_deltas) == D
            restored = CLS.decode_index(stored)
            @test typeof(restored) == typeof(original)
            @test restored.partition_bounds == original.partition_bounds
            a, b = restored.partitions[1], original.partitions[1]
            @test a.fragments == b.fragments
            @test a.rt_bins == b.rt_bins
            @test a.local_to_global == b.local_to_global
            @test a.skip_hints == b.skip_hints
            for field in fieldnames(typeof(a.fragment_bins))
                @test isequal(getfield(a.fragment_bins, field), getfield(b.fragment_bins, field))
            end
            counter = Pioneer.Counter(I, UInt8, 5)
            Pioneer.searchFragmentBinUnconditional!(counter, a.fragments, UInt32(1):UInt32(4))
            @test counter.counts[1:4] == UInt8[1, 128, 255, 0]
            mktempdir() do directory
                file = joinpath(directory, "index.jls")
                CLS.write_index(file, original)
                @test CLS.read_index(file).partitions[1].fragments == b.fragments
                Pioneer.serialize_to_jls(file, original)
                @test_throws MethodError CLS.read_index(file)
            end
        end
        for I in (UInt16, UInt32)
            Q = Pioneer.local_partition_type(I){Float32}
            empty = Pioneer.local_index_type(I){Float32}(Q[], Tuple{Float32,Float32}[], 0)
            @test typeof(CLS.decode_index(CLS.encode_index(empty))) == typeof(empty)
        end
        bad = storage_test_index(UInt16, 1)
        bad.partitions[1].fragment_bins.first_bins[1] = 2
        @test_throws ArgumentError CLS.encode_index(bad)
    end

    @testset "Scalar and spline values, annotations, and empty precursor ranges" begin
        offsets = UInt64[1, 3, 3, 4]
        charges = UInt8[2, 3, 4]
        fragments = [Pioneer.CompactFrag{Float32}(UInt32(pid), mz, intensity, a, b)
            for (pid, mz, intensity, a, b) in [
                (1, -0f0, Float16(-0.0), UInt16(0x0fff & ~(0x7 << 6) | (2 << 6)), UInt16(0x07ff)),
                (1, 200f0, reinterpret(Float16, UInt16(0x7e05)), UInt16(2 << 6), UInt16(0)),
                (3, 400f0, Float16(500), UInt16(4 << 6 | 0x11), UInt16(63))]]
        for original in (fragments, [Pioneer.SplineCompactFrag{4,Float32}(f.prec_id, f.mz,
                (1f0, -0f0, reinterpret(Float32, UInt32(0x7fc00005)), 4f0), f.packed_a, f.packed_b) for f in fragments])
            encoded = CLS.encode_details(original, offsets, charges)
            @test eltype(encoded.counts) == UInt8
            decoded, ptr = CLS.decode_details(encoded, charges)
            @test ptr == offsets
            @test typeof(decoded) == typeof(original)
            for (actual, expected) in zip(decoded, original)
                test_fragment_bits(actual, expected)
            end
            mktempdir() do directory
                file = joinpath(directory, "details.jls")
                CLS.write_details(file, original, offsets, charges)
                loaded, ranges = CLS.read_details(file, charges)
                @test ranges == offsets
                for (actual, expected) in zip(loaded, original)
                    test_fragment_bits(actual, expected)
                end
            end
            @test_throws ArgumentError CLS.encode_details(original, UInt64[1, 3, 2, 4], charges)
            @test_throws ArgumentError CLS.encode_details(original, offsets, UInt8[1, 3, 4])
        end
        empty = Pioneer.CompactFrag{Float32}[]
        @test CLS.decode_details(CLS.encode_details(empty, UInt64[1], UInt8[]), UInt8[]) == (empty, UInt64[1])
        @test_throws ArgumentError CLS.encode_details(empty, UInt64[1, 1], UInt8[8])
        for (n, C) in ((255, UInt8), (256, UInt16), (65535, UInt16), (65536, UInt32))
            f = Pioneer.CompactFrag{Float32}(UInt32(1), 500f0, Float16(1), UInt16(2 << 6), UInt16(0))
            stored = CLS.encode_details(fill(f, n), UInt64[1, n + 1], UInt8[2])
            @test eltype(stored.counts) == C
        end
    end

    @testset "Mapped precursor dictionaries and derived prediction inputs" begin
        @test CLS.koina_sequence("ACMK", "(2,C,Unimod:4)(3,M,Unimod:35)") == "AC[UNIMOD:4]M[UNIMOD:35]K"
        @test CLS.koina_sequence("AK", "(0,N-term,Unimod:1)(2,K,Unimod:737)") == "[UNIMOD:1]AK[UNIMOD:737]"
        @test CLS.koina_sequence("AK", missing) == "AK"
        @test_throws ArgumentError CLS.koina_sequence("AK", "(3,X,Unimod:1)")
        @test_throws ArgumentError CLS.koina_sequence("AK", "malformed")
        mktempdir() do directory
            original = DataFrame(sequence = ["ASCK", "ASCK"],
                structural_mods = Union{Missing,String}["(3,C,Unimod:4)", "(2,S,Unimod:21)(3,C,Unimod:4)"],
                koina_sequence = ["ASC[UNIMOD:4]K", "AS[UNIMOD:21]C[UNIMOD:4]K"],
                accession_numbers = ["P1", "P2"], proteome_identifiers = ["HUMAN", "HUMAN"],
                isotopic_mods = Union{Missing,String}[missing, missing], prec_charge = UInt8[2, 2],
                mz = Float32[500, 501], is_decoy = [false, true], start_idx = [UInt32[1], UInt32[1]])
            path = joinpath(directory, "precursors.arrow")
            CLS.write_precursors(path, original)
            table = CLS.read_precursors(path)
            @test !hasproperty(table, :koina_sequence)
            @test table.sequence isa Arrow.DictEncoded
            for name in filter(!=(:koina_sequence), propertynames(original))
                @test isequal(collect(getproperty(table, name)), original[!, name])
            end
            native_path = joinpath(directory, "native.arrow")
            Arrow.write(native_path, original)
            native = Pioneer.SetPrecursors(Arrow.Table(native_path))
            compact = Pioneer.SetPrecursors(table)
            @test Pioneer.getSequence(compact) isa Arrow.DictEncoded
            @test compact.pid_to_cv_fold == native.pid_to_cv_fold
            a, b = Pioneer.getPrecursorTextIds(native), Pioneer.getPrecursorTextIds(compact)
            for field in fieldnames(typeof(a))
                @test isequal(getfield(a, field), getfield(b, field))
            end
            original.koina_sequence[1] = "not derivable"
            @test_throws ArgumentError CLS.write_precursors(path, original)
            Arrow.write(path, original)
            @test_throws ArgumentError CLS.read_precursors(path)
        end
    end
end
