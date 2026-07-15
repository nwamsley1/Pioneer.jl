@testset "single-scan peak integrates at least three points" begin
    rt = Float32.(1:5)
    scan_idx = UInt32.(1:5)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 0, 10, 0, 0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)

    _, best_scan, _, start_scan, stop_scan = Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        3,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
    )

    @test best_scan == UInt32(3)
    @test start_scan == UInt32(2)
    @test stop_scan == UInt32(4)
    @test debug_plot_data[].scan_range == 2:4
end

@testset "baseline subtraction uses selected boundary endpoints" begin
    rt = Float32.(1:8)
    scan_idx = UInt32.(1:8)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[10, 12, 18, 40, 22, 18, 16, 14]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        4,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        forced_boundary_start_scan = UInt32(1),
        forced_boundary_stop_scan = UInt32(8),
        debug_plot_data = debug_plot_data,
    )

    baseline_subtracted = debug_plot_data[].baseline_subtracted
    @test baseline_subtracted[1] ≈ 0.0f0 atol = 1.0f-4
    @test baseline_subtracted[end] ≈ 0.0f0 atol = 1.0f-4
    @test debug_plot_data[].scan_range == 1:8
end

@testset "precursor integration persists selected boundary scans" begin
    chromatograms = DataFrame(
        precursor_idx = fill(UInt32(1), 5),
        rt = Float32.(1:5),
        scan_idx = UInt32.(10:10:50),
        intensity = Float32[0, 0, 10, 0, 0],
        precursor_fraction_transmitted = fill(1.0f0, 5),
    )
    peak_area = zeros(Float32, 1)
    new_best_scan = zeros(UInt32, 1)
    points_integrated = zeros(UInt32, 1)
    integration_start_scan = zeros(UInt32, 1)
    integration_stop_scan = zeros(UInt32, 1)

    Pioneer.integrate_precursors(
        chromatograms,
        0.25f0,
        UInt32[1],
        UInt32[30],
        peak_area,
        new_best_scan,
        points_integrated,
        integration_start_scan,
        integration_stop_scan;
        λ = 0.0f0,
    )

    @test new_best_scan == UInt32[30]
    @test integration_start_scan == UInt32[20]
    @test integration_stop_scan == UInt32[40]
end
