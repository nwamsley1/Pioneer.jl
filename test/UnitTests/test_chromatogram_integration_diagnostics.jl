@testset "window_trapezoid_area reproduces the quantified area" begin
    # The diagnostic area is only comparable to peak_area if the two use the same
    # quadrature. Integrating the *baseline-subtracted* trace with the diagnostic
    # helper must therefore reproduce what fillState! + integrateTrapezoidal
    # produced for the quantified result.
    rt = Float32.(1:9)
    scan_idx = UInt32.(1:9)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[5, 8, 20, 60, 140, 70, 25, 9, 6]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)

    peak_area, _, _, _, _, _, _, _, width = Pioneer.integrate_chrom(
        rt, scan_idx, intensity, fraction, 5, ws, state, 1.0f0, 0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
    )

    scan_range = debug_plot_data[].scan_range
    subtracted = debug_plot_data[].baseline_subtracted
    replayed = Pioneer.window_trapezoid_area(
        Vector{Float32}(subtracted), rt, first(scan_range), last(scan_range), 0, 1.0f0,
    )

    @test peak_area > 0.0f0
    @test replayed ≈ peak_area rtol = 1.0f-4
    @test width == UInt32(length(scan_range))
end

@testset "integration window is at least five scans wide" begin
    # getIntegrationBounds! seeds its search at apex +/- 2 and both searches only
    # move outward, so a lone spike still yields a five-scan window. A
    # points_integrated below five therefore means the trace collapsed inside a
    # full-width window, not that a narrower window was chosen.
    rt = Float32.(1:11)
    scan_idx = UInt32.(1:11)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)

    _, _, points_integrated, _, _, _, _, _, width = Pioneer.integrate_chrom(
        rt, scan_idx, intensity, fraction, 6, ws, state, 1.0f0, 0.0f0;
        min_fraction_transmitted = 0.25f0,
    )

    @test width >= UInt32(5)
    @test points_integrated < width
end

@testset "baseline subtraction on a slope leaves almost nothing" begin
    # The failure mode seen in the survey: the window sits on a monotone shoulder
    # rather than over a peak, so the endpoint-anchored baseline is nearly the
    # signal itself. peak_area collapses while peak_area_unsubtracted does not,
    # which is exactly the discrimination these diagnostics are meant to provide.
    rt = Float32.(1:9)
    scan_idx = UInt32.(1:9)
    fraction = fill(1.0f0, length(rt))
    slope = Float32[900, 800, 700, 600, 500, 400, 300, 200, 100]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)

    area, _, _, _, _, apex_smoothed, apex_subtracted, area_unsubtracted, _ =
        Pioneer.integrate_chrom(
            rt, scan_idx, slope, fraction, 5, ws, state, 1.0f0, 0.0f0;
            min_fraction_transmitted = 0.25f0,
            forced_boundary_start_scan = UInt32(2),
            forced_boundary_stop_scan = UInt32(8),
        )

    @test apex_smoothed > 0.0f0
    @test area_unsubtracted > 0.0f0
    # A straight line is entirely baseline: essentially none of it survives.
    @test apex_subtracted / apex_smoothed < 0.05f0
    @test area / area_unsubtracted < 0.05f0

    # A real peak over the same window keeps most of its signal, so the ratio
    # separates the two cases rather than being small everywhere.
    peak = Float32[100, 120, 300, 900, 1600, 900, 300, 120, 100]
    Pioneer.reset!(state)
    peak_area, _, _, _, _, peak_apex_smoothed, peak_apex_subtracted, peak_unsub, _ =
        Pioneer.integrate_chrom(
            rt, scan_idx, peak, fraction, 5, ws, state, 1.0f0, 0.0f0;
            min_fraction_transmitted = 0.25f0,
            forced_boundary_start_scan = UInt32(2),
            forced_boundary_stop_scan = UInt32(8),
        )

    @test peak_apex_subtracted / peak_apex_smoothed > 0.5f0
    @test peak_area / peak_unsub > 0.5f0
end
