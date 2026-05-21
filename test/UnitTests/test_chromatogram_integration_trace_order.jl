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
    boundary_selection_method = nothing,
    use_boundary_candidate_scoring = nothing,
)
    chrom_fields = String[]
    trace_mode !== nothing && push!(chrom_fields, "\"trace_mode\": \"$(trace_mode)\"")
    deconvolution_solver !== nothing && push!(chrom_fields, "\"deconvolution_solver\": \"$(deconvolution_solver)\"")
    boundary_selection_method !== nothing &&
        push!(chrom_fields, "\"boundary_selection_method\": \"$(boundary_selection_method)\"")
    use_boundary_candidate_scoring !== nothing &&
        push!(chrom_fields, "\"use_boundary_candidate_scoring\": $(use_boundary_candidate_scoring)")

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
    @test default_integration.boundary_selection_method == "learned"
    @test !hasproperty(default_integration, :use_boundary_candidate_scoring)
    @test !hasproperty(default_integration, :learned_boundary_max_train_groups)
    @test default_integration.learned_boundary_max_isotope_trace_groups_per_file == 2000
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

    legacy_candidate_path = joinpath(tmp, "legacy_candidate.json")
    _write_trace_mode_params(legacy_candidate_path; use_boundary_candidate_scoring = true)
    legacy_candidate_params = Pioneer.parse_pioneer_parameters(legacy_candidate_path)
    legacy_candidate_integration = Pioneer.IntegrateChromatogramSearchParameters(legacy_candidate_params)
    @test legacy_candidate_integration.boundary_selection_method == "learned"
    @test !hasproperty(legacy_candidate_integration, :use_boundary_candidate_scoring)
    @test !Pioneer.write_intermediate_chromatogram_debug_plots(legacy_candidate_integration)

    second_derivative_path = joinpath(tmp, "second_derivative.json")
    _write_trace_mode_params(second_derivative_path; boundary_selection_method = "second_derivative")
    second_derivative_params = Pioneer.parse_pioneer_parameters(second_derivative_path)
    second_derivative_integration = Pioneer.IntegrateChromatogramSearchParameters(second_derivative_params)
    @test second_derivative_integration.boundary_selection_method == "second_derivative"
    @test Pioneer.write_intermediate_chromatogram_debug_plots(second_derivative_integration)

    explicit_candidate_path = joinpath(tmp, "explicit_candidate.json")
    _write_trace_mode_params(explicit_candidate_path; boundary_selection_method = "candidate_score")
    explicit_candidate_params = Pioneer.parse_pioneer_parameters(explicit_candidate_path)
    @test_throws ArgumentError Pioneer.IntegrateChromatogramSearchParameters(explicit_candidate_params)

    separate_path = joinpath(tmp, "separate.json")
    _write_trace_mode_params(separate_path; trace_mode = "separate")
    checked = Pioneer.checkParams(separate_path)
    @test checked["optimization"]["chromatogram_integration"]["trace_mode"] == "separate"

    separate_params = Pioneer.parse_pioneer_parameters(separate_path)
    separate_integration = Pioneer.IntegrateChromatogramSearchParameters(separate_params)
    @test separate_integration.isotope_tracetype isa Pioneer.SeperateTraces
end
