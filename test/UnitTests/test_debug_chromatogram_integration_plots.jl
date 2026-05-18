function _with_debug_chromatogram_plot_state(f::Function)
    old_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
    old_targets = copy(Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[])

    try
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[] = Set(UInt32[909488])
        withenv("GKSwstype" => "100", "GKS_WSTYPE" => "100") do
            f()
        end
    finally
        Pioneer.DEBUG_CONSOLE_LEVEL[] = old_level
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[] = old_targets
    end
end

@testset "targeted chromatogram integration debug plots" begin
    chromatograms = DataFrame(
        precursor_idx = UInt32[909488, 909488, 909488, 591443, 591443, 591443, 123],
        rt = Float32[1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0],
        scan_idx = UInt32[10, 20, 30, 10, 20, 30, 10],
        intensity = Float32[0.5, 3.0, 0.7, 0.4, 2.5, 0.6, 1.0],
        precursor_fraction_transmitted = Float32[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    )
    passing_psms = DataFrame(
        precursor_idx = UInt32[909488, 591443, 123],
        scan_idx = UInt32[20, 20, 10],
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
