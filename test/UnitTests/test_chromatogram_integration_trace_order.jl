function _trace_order_test_chromatograms()
    return DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10, 11],
        isotopes_captured = [(Int8(0), Int8(1)), (Int8(0), Int8(1)),
                             (Int8(1), Int8(2)), (Int8(1), Int8(2)),
                             (Int8(0), Int8(1))],
        rt = Float32[2.0, 4.0, 1.0, 3.0, 1.5],
        scan_idx = UInt32[20, 40, 10, 30, 15],
        intensity = Float32[2.0, 4.0, 1.0, 3.0, 5.0],
        precursor_fraction_transmitted = Float32[0.7, 0.7, 0.35, 0.35, 0.8],
    )
end

@testset "combined chromatogram integration trace ordering" begin
    chromatograms = _trace_order_test_chromatograms()

    Pioneer.sort_chromatograms_for_integration!(chromatograms, Pioneer.CombineTraces(0.25f0))

    target_rows = findall(==(UInt32(10)), chromatograms.precursor_idx)
    @test chromatograms.rt[target_rows] == Float32[1.0, 2.0, 3.0, 4.0]
    @test chromatograms.scan_idx[target_rows] == UInt32[10, 20, 30, 40]
    @test chromatograms.isotopes_captured[target_rows] == [
        (Int8(1), Int8(2)),
        (Int8(0), Int8(1)),
        (Int8(1), Int8(2)),
        (Int8(0), Int8(1)),
    ]

    chrom_index, max_chrom_len = Pioneer.build_chrom_index(
        chromatograms,
        Pioneer.CombineTraces(0.25f0),
    )
    @test max_chrom_len == 4
    @test chrom_index[UInt32(10)] == 1:4
end

@testset "separate chromatogram integration trace ordering" begin
    chromatograms = _trace_order_test_chromatograms()

    Pioneer.sort_chromatograms_for_integration!(chromatograms, Pioneer.SeperateTraces())

    target_rows = findall(==(UInt32(10)), chromatograms.precursor_idx)
    @test chromatograms.rt[target_rows] == Float32[2.0, 4.0, 1.0, 3.0]
    @test chromatograms.isotopes_captured[target_rows] == [
        (Int8(0), Int8(1)),
        (Int8(0), Int8(1)),
        (Int8(1), Int8(2)),
        (Int8(1), Int8(2)),
    ]

    chrom_index, max_chrom_len = Pioneer.build_chrom_index(
        chromatograms,
        Pioneer.SeperateTraces(),
    )
    @test max_chrom_len == 2
    @test chrom_index[(UInt32(10), (Int8(0), Int8(1)))] == 1:2
    @test chrom_index[(UInt32(10), (Int8(1), Int8(2)))] == 3:4
end

@testset "separate trace quantification uses highest transmission trace" begin
    chromatograms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10, 20],
        isotopes_captured = [(Int8(1), Int8(2)), (Int8(0), Int8(1)),
                             (Int8(1), Int8(2)), (Int8(0), Int8(1)),
                             (Int8(0), Int8(1))],
        rt = Float32[1.0, 1.0, 2.0, 2.0, 1.0],
        scan_idx = UInt32[10, 10, 20, 20, 10],
        intensity = Float32[1.0, 10.0, 2.0, 20.0, 3.0],
        precursor_fraction_transmitted = Float32[0.35, 0.82, 0.40, 0.80, 0.55],
    )
    passing_psms = DataFrame(
        precursor_idx = UInt32[10, 20],
        isotopes_captured = [(Int8(1), Int8(2)), (Int8(1), Int8(2))],
        precursor_fraction_transmitted = Float32[0.35, 0.1],
    )

    selected = Pioneer.select_quant_trace_by_transmission(chromatograms)
    Pioneer.apply_quant_trace_selection!(passing_psms, selected)

    @test passing_psms.isotopes_captured == [(Int8(0), Int8(1)), (Int8(0), Int8(1))]
    @test passing_psms.precursor_fraction_transmitted == Float32[0.82, 0.55]
end

@testset "combined trace seed selection uses highest score then lower scan" begin
    passing_psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 20, 20],
        scan_idx = UInt32[40, 20, 30, 80, 60],
        lgbm_score = Float32[0.80, 0.95, 0.95, 0.70, 0.70],
        precursor_fraction_transmitted = Float32[0.95, 0.40, 0.80, 0.80, 0.70],
    )

    selected_rows = Pioneer.select_combined_trace_seed_rows_by_score(passing_psms)
    selected_psms = Pioneer.select_combined_trace_seed_psms_by_score(passing_psms)

    @test selected_rows == [2, 5]
    @test selected_psms.precursor_idx == UInt32[10, 20]
    @test selected_psms.scan_idx == UInt32[20, 60]
    @test selected_psms.lgbm_score == Float32[0.95, 0.70]
end

@testset "legacy chromatogram integration sort defaults to combined" begin
    chromatograms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10, 11],
        isotopes_captured = [(Int8(0), Int8(1)), (Int8(0), Int8(1)),
                             (Int8(1), Int8(2)), (Int8(1), Int8(2)),
                             (Int8(0), Int8(1))],
        rt = Float32[2.0, 4.0, 1.0, 3.0, 1.5],
        scan_idx = UInt32[20, 40, 10, 30, 15],
        intensity = Float32[2.0, 4.0, 1.0, 3.0, 5.0],
        precursor_fraction_transmitted = Float32[0.7, 0.7, 0.35, 0.35, 0.8],
    )

    Pioneer.sort_chromatograms_for_integration!(chromatograms)

    target_rows = findall(==(UInt32(10)), chromatograms.precursor_idx)
    @test chromatograms.rt[target_rows] == Float32[1.0, 2.0, 3.0, 4.0]
    @test chromatograms.scan_idx[target_rows] == UInt32[10, 20, 30, 40]
    @test chromatograms.isotopes_captured[target_rows] == [
        (Int8(1), Int8(2)),
        (Int8(0), Int8(1)),
        (Int8(1), Int8(2)),
        (Int8(0), Int8(1)),
    ]
end

function _write_trace_mode_params(
    path::AbstractString;
    trace_mode = nothing,
    deconvolution_solver = nothing,
)
    chrom_fields = String[]
    trace_mode !== nothing && push!(chrom_fields, "\"trace_mode\": \"$(trace_mode)\"")
    deconvolution_solver !== nothing && push!(chrom_fields, "\"deconvolution_solver\": \"$(deconvolution_solver)\"")

    trace_mode_json = isempty(chrom_fields) ? "" : """
        ,
        "optimization": {
            "chromatogram_integration": {
                $(join(chrom_fields, ",\n                "))
            }
        }
    """

    write(path, """
    {
        "paths": {
            "ms_data": "/tmp/pioneer-ms-data",
            "library": "/tmp/pioneer-lib.poin",
            "results": "/tmp/pioneer-results"
        }
        $(trace_mode_json)
    }
    """)
end

@testset "chromatogram integration trace mode config" begin
    tmp = mktempdir()

    default_path = joinpath(tmp, "default.json")
    _write_trace_mode_params(default_path)
    default_params = Pioneer.parse_pioneer_parameters(default_path)
    default_integration = Pioneer.IntegrateChromatogramSearchParameters(default_params)
    @test default_integration.isotope_tracetype isa Pioneer.CombineTraces
    @test isdefined(Pioneer, :HuberSolver)
    if isdefined(Pioneer, :HuberSolver)
        HuberSolver = getproperty(Pioneer, :HuberSolver)
        @test default_integration.deconvolution_solver isa HuberSolver
        @test default_integration.deconvolution_solver.delta == 300.0f0
        @test Pioneer.chromatogram_integration_solver_label(default_integration.deconvolution_solver) ==
            "HuberSolver(delta=300.0)"
    end

    pmm_path = joinpath(tmp, "pmm.json")
    _write_trace_mode_params(pmm_path; deconvolution_solver = "pmm")
    checked_pmm = Pioneer.checkParams(pmm_path)
    @test checked_pmm["optimization"]["chromatogram_integration"]["deconvolution_solver"] == "pmm"
    pmm_params = Pioneer.parse_pioneer_parameters(pmm_path)
    pmm_integration = Pioneer.IntegrateChromatogramSearchParameters(pmm_params)
    @test pmm_integration.deconvolution_solver isa Pioneer.PoissonMMSolver
    @test Pioneer.chromatogram_integration_solver_label(pmm_integration.deconvolution_solver) ==
        "PoissonMMSolver"

    separate_path = joinpath(tmp, "separate.json")
    _write_trace_mode_params(separate_path; trace_mode = "separate")
    checked = Pioneer.checkParams(separate_path)
    @test checked["optimization"]["chromatogram_integration"]["trace_mode"] == "separate"

    separate_params = Pioneer.parse_pioneer_parameters(separate_path)
    separate_integration = Pioneer.IntegrateChromatogramSearchParameters(separate_params)
    @test separate_integration.isotope_tracetype isa Pioneer.SeperateTraces
end

@testset "chromatogram transmission annotation keeps combined mode weighted" begin
    @test !Pioneer.compute_chromatogram_isotope_sets(Pioneer.CombineTraces(0.25f0))
    @test Pioneer.compute_chromatogram_isotope_sets(Pioneer.SeperateTraces())
end
