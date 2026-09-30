# Acquisition metadata: the user-declared scan mode, its sanity check, and what is recorded about a ZT run.

@testset "ZT Scan DIA window table and naming" begin
    W = S.SwathWindow
    step = 1.0221
    zt = [W(392.5 + step * (k - 1), 392.5 + step * k) for k in 1:496]      # contiguous ~1 Da Q1 bins
    @test S.zt_candidate_windows(zt)
    @test !S.zt_candidate_windows(zt[1:99])                                  # too few bins
    swath = [W(400.0 + 24.0 * (k - 1), 400.0 + 24.0 * k + 1.0) for k in 1:40]   # tens of wide, overlapping windows
    @test !S.zt_candidate_windows(swath)
    gappy = [W(392.5 + 1.1 * (k - 1), 392.5 + 1.1 * (k - 1) + 1.0) for k in 1:496]   # 1 Da bins with gaps
    @test !S.zt_candidate_windows(gappy)
    # stepped windows can be contiguous too (nSWATH: 173 × 2.9 Da), so the table alone never decides
    nswath = [W(400.0 + 2.9 * (k - 1), 400.0 + 2.9 * k) for k in 1:173]
    @test S.zt_candidate_windows(nswath)

    @test S.output_stem("run1", false) == "run1"
    @test S.output_stem("run1", true) == "run1.zt"
end

if !isempty(ZT_DATA_DIR)
    @testset "real data: ZT Scan DIA metadata" begin
        for wiff in filter(endswith(".wiff"), readdir(ZT_DATA_DIR; join = true))
            run = WiffRun(wiff)
            @test S.acquisition_metadata(run)["acquisition_type"] == "swath"      # never guessed
            m = S.acquisition_metadata(run; zt_scan = true)
            @test m["acquisition_type"] == "zt_scan_dia"
            @test m["acquisition_type_source"] == "user"
            @test !isempty(m["acquisition_method"])
            @test parse(Int, m["q1_bins_per_cycle"]) == length(run.windows)
            @test 0.5 < parse(Float64, m["q1_bin_step_mz"]) < 5
            @test 0.5 < parse(Float64, m["q1_bin_dwell_ms"]) < 20
        end
    end
end

if !isempty(DATA_DIR)
    @testset "real data: SWATH runs" begin
        for wiff in filter(endswith(".wiff"), readdir(DATA_DIR; join = true))
            run = WiffRun(wiff)
            m = S.acquisition_metadata(run)
            @test m["acquisition_type"] == "swath"
            @test !haskey(m, "q1_bin_step_mz")
            # declaring ZT is refused when the window table cannot be ZT
            S.zt_candidate_windows(run.windows) ||
                @test_throws ArgumentError S.acquisition_metadata(run; zt_scan = true)
        end
    end
end
