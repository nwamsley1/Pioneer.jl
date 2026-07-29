using DataFrames
using Test
using Pioneer

@testset "MBR chromatogram evidence does not inflate non-MBR rows" begin
    lean = Pioneer.MS2ChromObject(
        1.0f0,
        2.0f0,
        UInt32(3),
        UInt32(4),
    )
    rich = Pioneer.MS2MBRChromObject(
        1.0f0,
        2.0f0,
        UInt32(3),
        UInt32(4),
        ntuple(_ -> 0.0f0, 16)...,
    )
    @test sizeof(lean) == 16
    # spectrum_intensity is a per-SCAN quantity and now lives in a per-scan vector rather than
    # being copied onto every chromatogram row.
    @test sizeof(rich) == 80
    @test !hasproperty(DataFrame([rich]), :spectrum_intensity)
end

@testset "integrated MBR temporal evidence uses the selected peak bounds" begin
    psms = DataFrame(
        precursor_idx = UInt32[1],
        scan_idx = UInt32[11],
        new_best_scan = UInt32[11],
        peak_area = Float32[100],
        integration_start_scan = UInt32[10],
        integration_stop_scan = UInt32[12],
    )
    chromatograms = DataFrame(
        precursor_idx = fill(UInt32(1), 4),
        scan_idx = UInt32[9, 10, 11, 12],
        rt = Float32[1, 2, 3, 4],
        intensity = Float32[1, 2, 4, 2],
    )
    for rank in 1:8
        chromatograms[!, Symbol("shadow_frag$(rank)_int")] =
            rank == 1 ? Float32[100, 100, 100, 100] : zeros(Float32, 4)
        chromatograms[
            !,
            Pioneer.FITTED_FRAGMENT_INTENSITY_COLUMNS[rank],
        ] = rank == 1 ? Float32[100, 100, 100, 100] : zeros(Float32, 4)
    end

    scan_tic = zeros(Float32, 20)
    scan_tic[9:12] .= 100.0f0

    Pioneer.add_mbr_integrated_spectra_to_psms!(
        psms,
        chromatograms,
        identity;
        scan_tic = scan_tic,
    )

    trace = psms[1, Pioneer.MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN]
    @test length(trace) == 3 * Pioneer.MBR_TEMPORAL_TRACE_STRIDE
    @test trace[1:Pioneer.MBR_TEMPORAL_TRACE_STRIDE:end] ==
        Float32[2, 4, 2]
    @test psms[1, Pioneer.MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS[1]] ≈
        1.0f0
    @test psms[1, Pioneer.MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS[2]] ≈
        0.0f0
    @test psms[1, Pioneer.MBR_INTEGRATED_APEX_IRT_COLUMN] == 3.0f0
    @test psms[1, Pioneer.MBR_INTEGRATED_FITTED_HELLINGER_COLUMN] > 20.0f0
end

@testset "temporal Hellinger is zero for a matching donor spectrum" begin
    trace = Float32[
        2, 10, 0, 0, 0, 0, 0, 0, 0,
        1, 5, 0, 0, 0, 0, 0, 0, 0,
    ]
    donor = (
        1.0f0, 0.0f0, 0.0f0, 0.0f0,
        0.0f0, 0.0f0, 0.0f0, 0.0f0,
    )
    @test Pioneer._mbr_temporal_masked_hellinger(
        trace,
        donor,
        UInt8(0x03),
    ) ≈ 0.0f0
end
