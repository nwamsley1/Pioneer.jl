using Test
using Arrow
using DataFrames
using Pioneer

function make_two_peak_sparse_array(observed::Vector{Float32})
    H = Pioneer.SparseArray(UInt32(10))
    H.n_vals = 2
    H.m = 2
    H.n = 1
    H.rowval[1:2] .= UInt32[1, 2]
    H.colval[1:2] .= UInt16[1, 1]
    H.nzval[1:2] .= Float32[75, 25]
    H.x[1:2] .= observed
    H.matched[1:2] .= true
    H.isotope[1:2] .= UInt8[0, 0]
    H.colptr[1] = UInt32(1)
    H.colptr[2] = UInt32(3)
    return H
end

@testset "Dirichlet-multinomial tail feature" begin
    fill_time_ms = 1510.0f0
    weights = Float32[1]

    H_good = make_two_peak_sparse_array(Float32[75, 25])
    residuals_good = zeros(Float32, 2)
    good_tail = Pioneer.dirichlet_multinomial_tail_probability(
        weights,
        H_good,
        residuals_good,
        1,
        [1, 2],
        fill_time_ms
    )
    @test good_tail > 0.95f0

    residuals_bad = Float32[50, -50]
    bad_tail = Pioneer.dirichlet_multinomial_tail_probability(
        weights,
        H_good,
        residuals_bad,
        1,
        [1, 2],
        fill_time_ms
    )
    @test bad_tail < 0.05f0

    @test Pioneer.dirichlet_multinomial_tail_probability(
        weights,
        H_good,
        residuals_good,
        1,
        [1, 2],
        0.0f0
    ) == 1.0f0

    scores = Vector{Pioneer.SpectralScoresComplex{Float16}}(undef, 1)
    Pioneer.getDistanceMetrics(
        weights,
        residuals_good,
        H_good,
        scores;
        fill_time_ms = fill_time_ms
    )
    @test scores[1].dm_tail_prob > 0.95f0
end

@testset "DM tail probability is stored on complex PSMs" begin
    scored_psms = Vector{Pioneer.ComplexScoredPSM{Float32, Float16}}(undef, 1)
    unscored_psms = [
        Pioneer.ComplexUnscoredPSM(
            UInt8(1),
            UInt8(1),
            UInt8(2),
            UInt8(0),
            UInt8(2),
            UInt8(2),
            UInt8(0),
            UInt8(2),
            100.0f0,
            UInt8(2),
            200.0f0,
            UInt8(0),
            UInt8(0),
            1.0f0,
            UInt32(1),
            UInt32(1)
        )
    ]
    spectral_scores = [
        Pioneer.SpectralScoresComplex(
            Float16(0.9),
            Float16(0.8),
            Float16(10.0),
            Float16(10.0),
            Float16(10.0),
            Float16(10.0),
            Float16(1.0),
            Float16(2.0),
            Float16(0.0),
            0.25f0
        )
    ]
    id_to_col = Pioneer.ArrayDict(UInt32, UInt16, 4)
    Pioneer.update!(id_to_col, UInt32(1), UInt16(1))

    last_val = Pioneer.Score!(
        scored_psms,
        unscored_psms,
        spectral_scores,
        Float32[1.0],
        id_to_col,
        1,
        0.5,
        0,
        1,
        1000.0f0,
        1;
        min_spectral_contrast = 0.0f0,
        min_log2_matched_ratio = -10.0f0,
        min_y_count = 1,
        min_frag_count = 1,
        max_best_rank = 1,
        min_topn = 1
    )

    @test last_val == 1
    @test scored_psms[1].dm_tail_prob == 0.25f0
end

@testset "fillTimeMs accessors" begin
    mktempdir() do dir
        with_fill = joinpath(dir, "with_fill.arrow")
        Arrow.write(with_fill, DataFrame(
            mz_array = [Float32[100], Float32[200]],
            msOrder = UInt8[2, 2],
            fillTimeMs = Float32[12.5, 20.0],
        ))
        ms_data = Pioneer.BasicMassSpecData(with_fill)
        @test Pioneer.hasFillTimeMs(ms_data)
        @test Pioneer.getFillTimeMs(ms_data, 2) == 20.0f0

        without_fill = joinpath(dir, "without_fill.arrow")
        Arrow.write(without_fill, DataFrame(
            mz_array = [Float32[100], Float32[200]],
            msOrder = UInt8[2, 2],
        ))
        ms_data_no_fill = Pioneer.BasicMassSpecData(without_fill)
        @test !Pioneer.hasFillTimeMs(ms_data_no_fill)
        @test Pioneer.getFillTimeMs(ms_data_no_fill, 1) == 0.0f0
    end
end
