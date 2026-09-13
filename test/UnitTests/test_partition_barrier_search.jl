using Test
using Pioneer

struct BarrierSearchSpectra <: Pioneer.MassSpecData
    masses::Vector{Union{Missing, Float32}}
    intensities::Vector{Union{Missing, Float32}}
    fail_scan::Int
end

Base.length(::BarrierSearchSpectra) = 6
Pioneer.getMzArray(s::BarrierSearchSpectra, ::Integer) = s.masses
function Pioneer.getIntensityArray(s::BarrierSearchSpectra, i::Integer)
    i == s.fail_scan && error("Injected index worker failure")
    return s.intensities
end
Pioneer.getRetentionTime(::BarrierSearchSpectra, i::Integer) = Float32(i / 2)
Pioneer.getCenterMz(::BarrierSearchSpectra, i::Integer) = Float32(490 + i * 5)
Pioneer.getIsolationWidthMz(::BarrierSearchSpectra, ::Integer) = 20.0f0

function barrier_search_index(; no_overlap=false)
    bins = Pioneer.SoAFragBins{Float32}(
        Float32[100, 200, 300], vcat(Float32[101, 201, 301], fill(Inf32, 7)),
        UInt32[1, 2, 3], UInt32[1, 2, 3])
    rt_bins = [Pioneer.FragIndexBin{Float32}(0f0, 10f0, UInt32(1), UInt32(3))]
    fragments = [Pioneer.LocalFragment(UInt16(id), UInt8(score))
                 for (id, score) in ((1, 1), (2, 2), (1, 4))]
    parts = [
        Pioneer.LocalPartition{Float32}(bins, rt_bins, fragments,
            UInt32[1, 2], UInt16(2), UInt16[1, 1, 1]),
        Pioneer.LocalPartition{Float32}(Pioneer.SoAFragBins{Float32}(
            Float32[], fill(Inf32, 7), UInt32[], UInt32[]),
            Pioneer.FragIndexBin{Float32}[], Pioneer.LocalFragment[],
            UInt32[], UInt16(0), UInt16[]),
        Pioneer.LocalPartition{Float32}(bins, rt_bins, fragments,
            UInt32[3, 4], UInt16(2), UInt16[1, 1, 1]),
    ]
    bounds = no_overlap ? [(900f0, 910f0), (910f0, 920f0), (920f0, 930f0)] :
                          [(490f0, 500f0), (500f0, 510f0), (510f0, 520f0)]
    return Pioneer.LocalPartitionedFragmentIndex{Float32}(parts, bounds, 3)
end

function run_barrier_search(params, n_workers; scratch=nothing, no_overlap=false,
        max_peaks=0, min_score=1, fail_scan=0)
    spectra = BarrierSearchSpectra(Union{Missing, Float32}[100.5, 200.5, 300.5],
        Union{Missing, Float32}[20, 20, 1], fail_scan)
    # Preserve original scan IDs despite a noncontiguous, reordered scan selection.
    scan_ids = [6, 2, 4, 1]
    ranges = fill!(Vector{Union{Missing, UnitRange{Int64}}}(undef, 6), missing)
    result = Pioneer.searchFragmentIndexPartitionMajorHinted(
        ranges, barrier_search_index(; no_overlap), spectra, scan_ids, n_workers,
        params, Pioneer.SquareQuadModel(0f0), Pioneer.MassErrorModel(0f0, (20f0, 20f0)),
        identity, 2f0, Float32[491, 499, 511, 519];
        score_filter=Pioneer.CountFilter(UInt8(min_score)), max_peaks, scratch)
    return result, ranges
end

@testset "Partition barrier preserves index results" begin
    mktempdir() do dir
        path = joinpath(dir, "params.json")
        write(path, """{"paths":{"ms_data":"unused","library":"unused.poin","results":"unused"}}""")
        parsed = Pioneer.parse_pioneer_parameters(path)
        main_params = Pioneer.get_parameters(Pioneer.MainSearch(), parsed)
        tuning_params = Pioneer.get_parameters(Pioneer.ParameterTuningSearch(), parsed)

        for params in (main_params, tuning_params), n_workers in (1, 2, 3, 8)
            scratch = Pioneer.FragIndexScratch(n_workers)
            for max_peaks in (0, 1), min_score in (1, 2), no_overlap in (false, true)
                result, ranges = run_barrier_search(params, n_workers;
                    max_peaks, min_score, no_overlap)
                ids, scores = result
                @test isempty(scores)
                if no_overlap || (max_peaks == 1 && min_score == 2)
                    @test isempty(ids)
                    @test all(ismissing, ranges)
                else
                    # Two workers encounter partition 3 before partition 1 for scan 4
                    # when the stable merge concatenates their private output buffers.
                    expected = n_workers == 2 ?
                        UInt32[3, 4, 1, 2, 3, 4, 3, 4, 1, 2, 1, 2] :
                        UInt32[3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2]
                    min_score == 2 && filter!(isodd, expected)
                    @test ids == expected
                    expected_ranges = min_score == 1 ?
                        Union{Missing, UnitRange{Int64}}[11:12, 3:6, missing, 7:10, missing, 1:2] :
                        Union{Missing, UnitRange{Int64}}[6:6, 2:3, missing, 4:5, missing, 1:1]
                    @test isequal(ranges, expected_ranges)
                end
                # Reuse scratch across populated and empty calls and require exact order.
                @test isequal(run_barrier_search(params, n_workers;
                    max_peaks, min_score, no_overlap, scratch), (result, ranges))
            end
        end

        # Failure occurs inside a worker, after the initial partition rendezvous.
        # The caller must receive the original error after all peers exit.
        err = try
            run_barrier_search(main_params, 3; fail_scan=4)
            nothing
        catch error
            error
        end
        @test err isa CompositeException
        @test occursin("Injected index worker failure", sprint(showerror, err))
        @test !isempty(first(first(run_barrier_search(main_params, 3))))
    end
end
