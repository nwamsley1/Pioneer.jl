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

@testset "not-quantifiable rule withholds a window that is mostly baseline" begin
    # The rule compares the window area before baseline subtraction against what
    # survives it, and withholds the area when the unsubtracted area is at least
    # QUANT_MIN_AREA_SURVIVING_RATIO times larger. `quant_withheld` is
    # integrate_chrom's record of that decision, set where the rule fires, so it is
    # what these tests assert -- a zero peak_area on its own does not imply the rule
    # fired, which is the whole reason the flag replaced the unsubtracted area.
    rt = Float32.(1:9)
    scan_idx = UInt32.(1:9)
    fraction = fill(1.0f0, length(rt))
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)

    run_trace(trace) = begin
        Pioneer.reset!(state)
        Pioneer.integrate_chrom(
            rt, scan_idx, trace, fraction, 5, ws, state, 1.0f0, 0.0f0;
            min_fraction_transmitted = 0.25f0,
            forced_boundary_start_scan = UInt32(2),
            forced_boundary_stop_scan = UInt32(8),
        )
    end

    # A modest peak riding on a tall pedestal: the endpoint-anchored baseline takes
    # nearly the whole window, so the surviving area is a small fraction of it
    # (~8x here) and the area must be withheld rather than reported.
    on_pedestal = Float32[2000, 2020, 2200, 2800, 3500, 2800, 2200, 2020, 2000]
    area, _, _, _, _, _, _, withheld, _ = run_trace(on_pedestal)
    @test withheld
    @test area == 0.0f0                     # withheld, not merely small

    # A clean peak keeps most of its area and must be unaffected.
    peak = Float32[100, 120, 300, 900, 1600, 900, 300, 120, 100]
    peak_area, _, _, _, _, _, _, peak_withheld, _ = run_trace(peak)
    @test peak_area > 0.0f0
    @test !peak_withheld

    # The cut is on the ratio, so a window losing some but not most of its area
    # still quantifies. The same peak shape on a lower pedestal sits at ~3.5x,
    # below the cut, and brackets it against the withheld case above.
    pedestal = Float32[700, 720, 900, 1500, 2200, 1500, 900, 720, 700]
    ped_area, _, _, _, _, _, _, ped_withheld, _ = run_trace(pedestal)
    @test !ped_withheld
    @test ped_area > 0.0f0

    # A pure slope never reaches the rule at all: subtraction drives the apex to
    # zero and the `_apex_val <= 0` guard skips integration first. Its peak_area is
    # zero for a different reason, and quant_withheld says so -- the separation the
    # flag exists to provide.
    slope = Float32[900, 800, 700, 600, 500, 400, 300, 200, 100]
    slope_area, _, _, _, _, _, slope_apex_sub, slope_withheld, _ = run_trace(slope)
    @test slope_area == 0.0f0
    @test slope_apex_sub == 0.0f0
    @test !slope_withheld
end

@testset "baseline subtraction on a slope leaves almost nothing" begin
    # The failure mode seen in the survey: the window sits on a monotone shoulder
    # rather than over a peak, so the endpoint-anchored baseline is nearly the
    # signal itself and both the apex and the window area collapse relative to
    # their pre-subtraction values. That collapse is what the diagnostics are meant
    # to expose and what the not-quantifiable rule tests.
    #
    # integrate_chrom returns the `quant_withheld` decision rather than the two
    # areas the rule compares, so recover both here the way the first testset
    # recovers the quantified area: replay window_trapezoid_area over the captured
    # traces and the window that was actually integrated. The surviving area has to
    # come from the trace, not from the returned `area`, which is already zero
    # whenever the rule fires or the zero-apex guard trips.
    rt = Float32.(1:9)
    scan_idx = UInt32.(1:9)
    fraction = fill(1.0f0, length(rt))
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)

    window_area(trace, dbg) = Pioneer.window_trapezoid_area(
        Vector{Float32}(trace), rt,
        first(dbg.scan_range), last(dbg.scan_range), 0, 1.0f0,
    )

    run_trace(trace) = begin
        Pioneer.reset!(state)
        dbg = Ref{Any}(nothing)
        result = Pioneer.integrate_chrom(
            rt, scan_idx, trace, fraction, 5, ws, state, 1.0f0, 0.0f0;
            min_fraction_transmitted = 0.25f0,
            forced_boundary_start_scan = UInt32(2),
            forced_boundary_stop_scan = UInt32(8),
            debug_plot_data = dbg,
        )
        (result, dbg[])
    end

    slope = Float32[900, 800, 700, 600, 500, 400, 300, 200, 100]
    (area, _, _, _, _, apex_smoothed, apex_subtracted, withheld, _), slope_dbg =
        run_trace(slope)

    @test apex_smoothed > 0.0f0
    @test window_area(slope_dbg.wh_smoothed, slope_dbg) > 0.0f0   # the window held signal
    # A straight line is entirely baseline: none of it survives, so the apex is
    # driven to exactly zero and integration is skipped before the rule is reached.
    @test apex_subtracted == 0.0f0
    @test window_area(slope_dbg.baseline_subtracted, slope_dbg) == 0.0f0
    @test slope_dbg.status == "skipped_zero_apex"
    @test area == 0.0f0
    @test !withheld

    # A real peak over the same window keeps most of its signal, so the ratio
    # separates the two cases rather than being small everywhere.
    peak = Float32[100, 120, 300, 900, 1600, 900, 300, 120, 100]
    (peak_area, _, _, _, _, peak_apex_smoothed, peak_apex_subtracted, peak_withheld, _),
        peak_dbg = run_trace(peak)
    peak_unsubtracted = window_area(peak_dbg.wh_smoothed, peak_dbg)
    peak_surviving = window_area(peak_dbg.baseline_subtracted, peak_dbg)

    @test peak_apex_subtracted / peak_apex_smoothed > 0.5f0
    @test peak_surviving ≈ peak_area rtol = 1.0f-4
    @test peak_surviving / peak_unsubtracted > 0.5f0
    @test peak_unsubtracted / peak_surviving < Pioneer.QUANT_MIN_AREA_SURVIVING_RATIO
    @test !peak_withheld

    # And the rule's own case sits on the far side of that cut: the same peak on a
    # tall pedestal integrates, but the baseline takes enough of the window that
    # the ratio clears QUANT_MIN_AREA_SURVIVING_RATIO and the area is withheld.
    on_pedestal = Float32[2000, 2020, 2200, 2800, 3500, 2800, 2200, 2020, 2000]
    (ped_area, _, _, _, _, _, _, ped_withheld, _), ped_dbg = run_trace(on_pedestal)
    ped_unsubtracted = window_area(ped_dbg.wh_smoothed, ped_dbg)
    ped_surviving = window_area(ped_dbg.baseline_subtracted, ped_dbg)

    @test ped_surviving > 0.0f0
    @test ped_unsubtracted / ped_surviving >= Pioneer.QUANT_MIN_AREA_SURVIVING_RATIO
    @test ped_withheld
    @test ped_area == 0.0f0
end
