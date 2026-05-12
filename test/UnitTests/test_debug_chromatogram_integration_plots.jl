function _with_debug_chromatogram_plot_state(f::Function)
    old_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
    old_target = Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDX[]

    try
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDX[] = UInt32(909488)
        withenv("GKSwstype" => "100", "GKS_WSTYPE" => "100") do
            f()
        end
    finally
        Pioneer.DEBUG_CONSOLE_LEVEL[] = old_level
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDX[] = old_target
    end
end

@testset "targeted chromatogram integration debug plots" begin
    chromatograms = DataFrame(
        precursor_idx = UInt32[909488, 909488, 909488, 123],
        rt = Float32[1.0, 2.0, 3.0, 1.0],
        scan_idx = UInt32[10, 20, 30, 10],
        intensity = Float32[0.5, 3.0, 0.7, 1.0],
        precursor_fraction_transmitted = Float32[1.0, 1.0, 1.0, 1.0],
    )
    passing_psms = DataFrame(
        precursor_idx = UInt32[909488, 123],
        scan_idx = UInt32[20, 10],
    )

    _with_debug_chromatogram_plot_state() do
        tmp = mktempdir()
        Pioneer.debug_write_target_chromatogram_plots(
            chromatograms,
            passing_psms,
            0.25f0,
            1.0f-6,
            tmp,
            2,
            "run A.raw",
        )

        plot_dir = joinpath(tmp, "qc_plots", "chromatogram_integration_debug")
        files = filter(path -> endswith(path, ".png"), readdir(plot_dir; join = true))
        @test length(files) == 1
        @test occursin("file_2", basename(only(files)))
        @test occursin("run_A_raw", basename(only(files)))
        @test occursin("precursor_909488", basename(only(files)))
        @test filesize(only(files)) > 0
    end

    _with_debug_chromatogram_plot_state() do
        tmp = mktempdir()
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 0
        Pioneer.debug_write_target_chromatogram_plots(
            chromatograms,
            passing_psms,
            0.25f0,
            1.0f-6,
            tmp,
            2,
            "run A.raw",
        )

        plot_dir = joinpath(tmp, "qc_plots", "chromatogram_integration_debug")
        @test !isdir(plot_dir)
    end
end

@testset "chromatogram debug plot data labels actual algorithm outputs" begin
    rt = Float32[1.0, 2.0, 3.0]
    scan_idx = UInt32[10, 20, 30]
    intensity = Float32[10.0, 20.0, 10.0]
    fraction = Float32[1.0, 1.0, 1.0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        2,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.0f0,
        debug_plot_data = debug_plot_data,
    )

    @test debug_plot_data[] !== nothing
    @test debug_plot_data[].wh_smoothed == intensity
    @test debug_plot_data[].baseline_subtracted != debug_plot_data[].wh_smoothed
    @test debug_plot_data[].baseline_subtracted == Float32[0.0, 10.0, 0.0]
end
