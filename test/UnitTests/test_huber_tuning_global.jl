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

@testset "Huber selects one global observation per precursor" begin
    empty_winner = Pioneer.HuberCalibrationWinner(-Inf32, 0, 0)
    winners = fill(empty_winner, 8)
    psms = DataFrame(precursor_idx=UInt32[1, 1, 1, 1, 2, 3, 4, 5, 6],
        prec_prob=Float32[0.8, 0.9, 0.9, 0.9, 0.95, 0.99, 0.99, NaN, Inf],
        ms_file_idx=UInt32[1, 3, 2, 2, 4, 1, 1, 1, 1],
        scan_idx=UInt32[1, 8, 7, 6, 9, 2, 3, 4, 5],
        mbr_recovered=Bool[0, 0, 0, 0, 0, 1, 0, 0, 0],
        MBR_transfer_candidate=Bool[0, 0, 0, 0, 0, 0, 1, 0, 0])
    @test Pioneer.collect_huber_winners!(winners, psms) === psms
    @test winners[1] == Pioneer.HuberCalibrationWinner(0.9f0, 2, 6)
    @test winners[2] == Pioneer.HuberCalibrationWinner(0.95f0, 4, 9)
    @test all(w -> w.file_idx == 0, winners[3:end])
    @test Pioneer.huber_winner_files(winners) == UInt32[2, 4]
    reversed_winners = fill(empty_winner, 8)
    Pioneer.collect_huber_winners!(reversed_winners, reverse(psms))
    @test reversed_winners == winners
    # Global selection precedes the scan-density cap: losing observations
    # in dense scans cannot displace a winner in a less crowded scan.
    file_psms = psms[psms.ms_file_idx .== 2, :]
    selected = Pioneer.select_huber_calibration_psms(
        Pioneer.global_huber_psms(file_psms, winners, 2), 1)
    @test selected.scan_idx == UInt32[6]
    @test isempty(Pioneer.global_huber_psms(psms[psms.ms_file_idx .== 1, :], winners, 1))
    @test sizeof(empty_winner) == 12
end

@testset "Duplicate runs do not multiply Huber calibration curves" begin
    winners = fill(Pioneer.HuberCalibrationWinner(-Inf32, 0, 0), 3)
    frame = DataFrame(precursor_idx=UInt32[1, 2, 3], prec_prob=Float32[0.9, 0.8, 0.7],
        ms_file_idx=fill(UInt32(1), 3), scan_idx=UInt32[5, 5, 6])
    retained_bytes = 0
    for run in 6000:-1:1
        frame.ms_file_idx .= UInt32(run)
        Pioneer.collect_huber_winners!(winners, frame)
        run == 6000 && (retained_bytes = Base.summarysize(winners))
    end
    @test Base.summarysize(winners) == retained_bytes
    @test Pioneer.huber_winner_files(winners) == UInt32[1]
    frame.ms_file_idx .= 1
    selected = Pioneer.global_huber_psms(frame, winners, 1)
    @test nrow(selected) == 3
    @test nrow(Pioneer.select_huber_calibration_psms(selected, 1)) == 2 # whole scan
    @test isempty(Pioneer.select_huber_calibration_psms(selected, 0))
    @test nrow(Pioneer.select_huber_calibration_psms(selected, 10)) == 3
end

@testset "Huber winner collection follows confidence filtering" begin
    mktempdir() do dir
        path = joinpath(dir, "psms.arrow")
        psms = DataFrame(precursor_idx=UInt32[1, 2, 3], prec_prob=Float32[0.9, 0.99, 0.99],
            ms_file_idx=UInt32[1, 1, 1], scan_idx=UInt32[1, 2, 3],
            global_qval=Float32[0.001, 0.001, 0.1], qval=Float32[0.001, 0.1, 0.001])
        Arrow.write(path, psms)
        winners = fill(Pioneer.HuberCalibrationWinner(-Inf32, 0, 0), 3)
        pipeline = Pioneer.TransformPipeline() |>
            Pioneer.filter_by_multiple_thresholds([(:global_qval, 0.01f0), (:qval, 0.01f0)]) |>
            ("collect_global_huber_winners" => (df -> Pioneer.collect_huber_winners!(winners, df)))
        refs = Pioneer.apply_pipeline_batch([Pioneer.PSMFileReference(path)], pipeline, joinpath(dir, "passing"))
        @test winners[1].file_idx == 1
        @test all(w -> w.file_idx == 0, winners[2:3])
        @test DataFrame(Arrow.Table(Pioneer.file_path(only(refs)))).precursor_idx == UInt32[1]
    end
end

struct HuberTestLibrary <: Pioneer.SpectralLibrary end
struct HuberTestDataReference <: Pioneer.MassSpecDataReference
    loaded_files::Vector{Int}
end
function Pioneer.getMSData(ref::HuberTestDataReference, idx::Int64)
    push!(ref.loaded_files, idx)
    error("test spectrum load")
end

@testset "Huber skips spectra without global winners and releases selection" begin
    mktempdir() do dir
        path = joinpath(dir, "params.json")
        _write_huber_tuning_params(path)
        params = Pioneer.parse_pioneer_parameters(path)
        ref = HuberTestDataReference(Int[])
        context = Pioneer.SearchContext(HuberTestLibrary(), Pioneer.SearchDataStructures[], ref, 1, 3, 1)
        context.huber_calibration_winners = fill(Pioneer.HuberCalibrationWinner(-Inf32, 0, 0), 3)
        Pioneer.execute_search(Pioneer.HuberTuningSearch(), context, params)
        @test isempty(ref.loaded_files)
        @test isempty(context.huber_calibration_winners)
        @test context.huber_delta[] == 300f0
        context.huber_calibration_winners = [Pioneer.HuberCalibrationWinner(0.9f0, 2, 7)]
        @test_throws ErrorException Pioneer.execute_search(Pioneer.HuberTuningSearch(), context, params)
        @test ref.loaded_files == [2]
        @test isempty(context.huber_calibration_winners)
        @test Pioneer.PrecursorScoringSearchParameters(params).calibrate_huber
        _write_huber_tuning_params(path; deconvolution_solver="pmm")
        @test !Pioneer.PrecursorScoringSearchParameters(Pioneer.parse_pioneer_parameters(path)).calibrate_huber
    end
end
