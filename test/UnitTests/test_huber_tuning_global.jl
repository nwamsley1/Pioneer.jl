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
