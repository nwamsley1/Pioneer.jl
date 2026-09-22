@testset "global Huber delta calibration" begin
    delta_grid = Float32[100, 200, 300]
    tuning_psms = DataFrame(
        ms_file_idx = UInt32[1, 1, 1, 2, 2, 2],
        precursor_idx = UInt32[42, 42, 42, 42, 42, 42],
        scan_idx = UInt32[7, 7, 7, 7, 7, 7],
        huber_delta = Float32[100, 200, 300, 100, 200, 300],
        weight = Float32[1, 3, 3, 1, 1, 3],
    )

    @test Pioneer.estimate_optimal_huber_delta(tuning_psms, delta_grid, 10.0f0) == 200.0f0
end

@testset "chromatogram integration uses calibrated Huber delta" begin
    base_solver = Pioneer.default_chromatogram_integration_huber_solver()
    calibrated = Pioneer.with_chromatogram_huber_delta(base_solver, 777.0f0)

    @test calibrated isa Pioneer.HuberSolver
    @test calibrated.delta == 777.0f0
    @test base_solver.delta == 300.0f0

    pmm_solver = Pioneer.PoissonMMSolver()
    @test Pioneer.with_chromatogram_huber_delta(pmm_solver, 777.0f0) === pmm_solver
end

function _write_huber_tuning_params(path::AbstractString; deconvolution_solver = "huber")
    write(path, """
    {
        "paths": {
            "ms_data": "/tmp/pioneer-ms-data",
            "library": "/tmp/pioneer-lib.poin",
            "results": "/tmp/pioneer-results"
        },
        "optimization": {
            "chromatogram_integration": {
                "deconvolution_solver": "$(deconvolution_solver)"
            }
        }
    }
    """)
end

@testset "Huber tuning follows chromatogram solver config" begin
    tmp = mktempdir()

    huber_path = joinpath(tmp, "huber.json")
    _write_huber_tuning_params(huber_path)
    huber_params = Pioneer.HuberTuningSearchParameters(Pioneer.parse_pioneer_parameters(huber_path))
    @test huber_params.enabled
    @test huber_params.base_solver.delta == 300.0f0
    @test 300.0f0 in huber_params.delta_grid

    pmm_path = joinpath(tmp, "pmm.json")
    _write_huber_tuning_params(pmm_path; deconvolution_solver = "pmm")
    pmm_params = Pioneer.HuberTuningSearchParameters(Pioneer.parse_pioneer_parameters(pmm_path))
    @test !pmm_params.enabled
end

# Reference reduction keeps the all-runs grouping used before streaming summaries.
function _reference_huber_histogram(psms, grid, threshold)
    delta_column = hasproperty(psms, :huber_delta) ? :huber_delta : Symbol("huber_δ")
    keys = hasproperty(psms, :ms_file_idx) ?
        [:ms_file_idx, :precursor_idx, :scan_idx] : [:precursor_idx, :scan_idx]
    curves = combine(groupby(psms, keys)) do group
        Pioneer.process_huber_curve(group.weight, group[!, delta_column])
    end
    filter!(row -> row.n == length(grid), curves)
    filter!(row -> row.wdiff > threshold / 100, curves)
    filter!(row -> !ismissing(row.huber50), curves)
    curves.huber50 = ceil.(Int, curves.huber50)
    hist = combine(groupby(curves, :huber50), nrow)
    sort!(hist, :huber50)
    return hist
end

@testset "streamed Huber summaries preserve global calibration" begin
    grid = Float32[100, 200, 300]
    runs = DataFrame[]
    for run in 1:12
        frame = Pioneer.empty_huber_tuning_results()
        curves = (
            Float32[1, 3, 3], Float32[1, 1, 3], Float32[3, 1, 1],
            Float32[1, 1, 1], Float32[0, 0, 0], Float32[NaN, 1, 3],
            Float32[1, 3], Float32[1, 2, 1],
        )
        for (pid, weights) in enumerate(curves)
            for i in reverse(eachindex(weights))
                push!(frame, (UInt32(run), UInt32(pid), UInt32(7), grid[i], weights[i]))
            end
        end
        # Unequal run contributions must retain curve weighting, not average run medians.
        if iseven(run)
            push!(frame, (UInt32(run), UInt32(99), UInt32(8), 100f0, 1f0))
            push!(frame, (UInt32(run), UInt32(99), UInt32(8), 200f0, 1f0))
            push!(frame, (UInt32(run), UInt32(99), UInt32(8), 300f0, 3f0))
        end
        push!(runs, frame)
    end
    combined = vcat(runs...)
    expected = _reference_huber_histogram(combined, grid, 10f0)
    expected_counts = Dict(zip(expected.huber50, expected.nrow))
    expected_delta = Pioneer.get_median_huber_delta(
        Float32.(cumsum(expected.nrow ./ sum(expected.nrow))), expected.huber50,
    )
    for ordered_runs in (runs, reverse(runs))
        histogram = Dict{Int, Int}()
        for frame in ordered_runs
            Pioneer.accumulate_huber_histogram!(histogram, frame, grid, 10f0)
        end
        @test histogram == expected_counts
        @test Pioneer.estimate_optimal_huber_delta(histogram) == expected_delta
    end
    @test Pioneer.estimate_optimal_huber_delta(combined, grid, 10f0) == expected_delta
    legacy = select(runs[1], Not(:ms_file_idx))
    rename!(legacy, :huber_delta => Symbol("huber_δ"))
    legacy_expected = _reference_huber_histogram(legacy, grid, 10f0)
    histogram = Dict{Int, Int}()
    Pioneer.accumulate_huber_histogram!(histogram, legacy, grid, 10f0)
    @test histogram == Dict(zip(legacy_expected.huber50, legacy_expected.nrow))
    @test Pioneer.accumulate_huber_histogram!(histogram,
        Pioneer.empty_huber_tuning_results(), grid, 10f0) == histogram
    rejected = combined[combined.precursor_idx .== UInt32(4), :]
    @test isempty(Pioneer.accumulate_huber_histogram!(Dict{Int, Int}(), rejected, grid, 10f0))
    @test_throws ArgumentError Pioneer.estimate_optimal_huber_delta(Dict{Int, Int}())
    @test_throws ArgumentError Pioneer.estimate_optimal_huber_delta(rejected, grid, 10f0)
end

@testset "Huber retained memory does not grow with run count" begin
    grid = Float32[100, 200, 300]
    frame = DataFrame(ms_file_idx=fill(UInt32(1), 3), precursor_idx=fill(UInt32(42), 3),
        scan_idx=fill(UInt32(7), 3), huber_delta=grid, weight=Float32[1, 3, 3])
    results = Pioneer.HuberTuningSearchResults(Ref(300f0), Dict{Int, Int}(), Ref(0))
    retained_bytes = 0
    for run in 1:6000
        frame.ms_file_idx .= UInt32(run)
        Pioneer.accumulate_huber_histogram!(results.huber_histogram, frame, grid, 10f0)
        results.n_observations[] += nrow(frame)
        Pioneer.reset_results!(results)
        run == 100 && (retained_bytes = Base.summarysize(results))
    end
    @test results.huber_histogram == Dict(150 => 6000)
    @test results.n_observations[] == 18_000
    @test Base.summarysize(results) == retained_bytes
    @test Pioneer.estimate_optimal_huber_delta(results.huber_histogram) == 150f0
end
