import Pioneer
using Arrow
using DataFrames
using Test

struct MBREmpiricalMockPrecursors <: Pioneer.LibraryPrecursors
    mz::Vector{Float32}
    irt::Vector{Float32}
end

Pioneer.getMz(p::MBREmpiricalMockPrecursors) = p.mz
Pioneer.getIrt(p::MBREmpiricalMockPrecursors) = p.irt
Pioneer.getCharge(p::MBREmpiricalMockPrecursors) = fill(UInt8(2), length(p.mz))
Pioneer.getLength(p::MBREmpiricalMockPrecursors) = fill(UInt8(9), length(p.mz))

function _mbr_empirical_fragment_lookup()
    frags = Pioneer.CompactFrag{Float32}[
        Pioneer.CompactFrag(
            UInt32(1), 500.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(3), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(1), 600.0f0, Float16(1),
            false, true, false, false,
            UInt8(1), UInt8(4), UInt8(2), UInt8(2), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(2), 600.0f0, Float16(1),
            false, true, false, false,
            UInt8(1), UInt8(4), UInt8(2), UInt8(1), UInt8(0),
        ),
        Pioneer.CompactFrag(
            UInt32(2), 500.0f0, Float16(1),
            true, false, false, false,
            UInt8(1), UInt8(3), UInt8(2), UInt8(2), UInt8(0),
        ),
    ]
    return Pioneer.StandardFragmentLookup{Float32}(frags, UInt64[1, 3, 5])
end

function _mbr_empirical_main_table(;
    precursor_idx,
    scan_idx,
    ms_file_idx,
    target,
    weight = ones(Float32, length(precursor_idx)),
    log2_intensity_explained = fill(Float16(1), length(precursor_idx)),
    n_scans = ones(Float32, length(precursor_idx)),
    main_pep = fill(0.5f0, length(precursor_idx)),
    library_hellinger = fill(0.25f0, length(precursor_idx)),
    frag_corr_bitvec = fill(UInt8(0x03), length(precursor_idx)),
    frag_corr_rank = fill(UInt16(11), length(precursor_idx)),
    frag1,
    frag2,
    frag3 = zeros(Float32, length(precursor_idx)),
)
    n = length(precursor_idx)
    return DataFrame(
        precursor_idx = UInt32.(precursor_idx),
        scan_idx = UInt32.(scan_idx),
        weight = Float32.(weight),
        log2_intensity_explained = Float16.(log2_intensity_explained),
        irt_pred = zeros(Float32, n),
        irt_obs = zeros(Float32, n),
        log_by_ratio_m0 = zeros(Float16, n),
        rt = zeros(Float32, n),
        n_scans = Float32.(n_scans),
        main_pep = Float32.(main_pep),
        smoothed_2d_shadow_hellinger = Float32.(library_hellinger),
        frag_corr_bitvec = UInt8.(frag_corr_bitvec),
        n_correlated_fragments_bitvec_rank = UInt16.(frag_corr_rank),
        ms_file_idx = fill(UInt32(ms_file_idx), n),
        target = Bool.(target),
        cv_fold = zeros(UInt8, n),
        frag1_smoothed_intensity = Float32.(frag1),
        frag2_smoothed_intensity = Float32.(frag2),
        frag3_smoothed_intensity = Float32.(frag3),
        frag4_smoothed_intensity = zeros(Float32, n),
        frag5_smoothed_intensity = zeros(Float32, n),
        frag6_smoothed_intensity = zeros(Float32, n),
        frag7_smoothed_intensity = zeros(Float32, n),
        frag8_smoothed_intensity = zeros(Float32, n),
    )
end

function _mbr_empirical_pass1_table(; precursor_idx, scan_idx, score)
    return DataFrame(
        precursor_idx = UInt32.(precursor_idx),
        scan_idx = UInt32.(scan_idx),
        trace_prob_prepass = Float32.(score),
        trace_prob_infold = Float32.(score),
    )
end

function _mbr_add_integrated_fixture_columns!(df, templates, temporal_means)
    n = nrow(df)
    df[!, Pioneer.MBR_INTEGRATED_WEIGHT_COLUMN] = Float32.(df.weight)
    df[!, Pioneer.MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN] =
        Float32.(df.log2_intensity_explained)
    df[!, Pioneer.MBR_INTEGRATED_APEX_IRT_COLUMN] = Float32.(df.irt_obs)
    df[!, Pioneer.MBR_INTEGRATED_N_SCANS_COLUMN] = UInt32.(df.n_scans)
    temporal_traces = Vector{Vector{Float32}}(undef, n)
    @inbounds for row in 1:n
        trace = Float32[1.0f0]
        for rank in 1:8
            push!(trace, Float32(temporal_means[row][rank]^2))
        end
        temporal_traces[row] = trace
    end
    df[!, Pioneer.MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN] = temporal_traces
    @inbounds for rank in 1:8
        df[!, Pioneer.MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS[rank]] =
            Float32[templates[row][rank] for row in 1:n]
        df[!, Pioneer.MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS[rank]] =
            Float32[temporal_means[row][rank] for row in 1:n]
    end
    return df
end

@testset "MBR temporal fragment summary uses integration bounds" begin
    passing_psms = DataFrame(
        precursor_idx = UInt32[1],
        scan_idx = UInt32[20],
        new_best_scan = UInt32[20],
        peak_area = Float32[4],
        integration_start_scan = UInt32[20],
        integration_stop_scan = UInt32[30],
    )
    chromatograms = DataFrame(
        precursor_idx = fill(UInt32(1), 4),
        scan_idx = UInt32[10, 20, 30, 40],
        rt = Float32[1, 2, 3, 4],
        intensity = Float32[100, 1, 3, 100],
        spectrum_intensity = fill(10.0f0, 4),
        shadow_frag1_int = Float32[0, 10, 0, 0],
        shadow_frag2_int = Float32[0, 0, 10, 0],
        shadow_frag3_int = Float32[1000, 0, 0, 1000],
        shadow_frag4_int = zeros(Float32, 4),
        shadow_frag5_int = zeros(Float32, 4),
        shadow_frag6_int = zeros(Float32, 4),
        shadow_frag7_int = zeros(Float32, 4),
        shadow_frag8_int = zeros(Float32, 4),
    )

    Pioneer.add_mbr_integrated_spectra_to_psms!(
        passing_psms,
        chromatograms,
        identity,
    )

    temporal_mean = ntuple(
        rank -> passing_psms[1, Pioneer.MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS[rank]],
        8,
    )
    @test temporal_mean[1] ≈ 0.25f0
    @test temporal_mean[2] ≈ 0.75f0
    @test temporal_mean[3] == 0.0f0
    temporal_trace = passing_psms[1, Pioneer.MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN]
    @test length(temporal_trace) == 2 * Pioneer.MBR_TEMPORAL_TRACE_STRIDE
    @test temporal_trace[[1, 2, 10, 12]] == Float32[1, 10, 3, 10]
    @test Pioneer._mbr_temporal_fragment_hellinger_from_summaries(
        temporal_mean,
        (1.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0),
    ) ≈ sqrt(0.75f0)
end

@testset "MBR masked temporal Hellinger normalizes each scan within the mask" begin
    trace = Float32[
        1, 1, 0, 0, 0, 0, 0, 0, 0,
        3, 0, 0, 9, 0, 0, 0, 0, 0,
    ]
    donor = (1.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0)

    @test Pioneer._mbr_temporal_corr_masked_hellinger_from_trace(
        trace,
        donor,
        UInt8(0x03),
    ) ≈ sqrt(0.75f0)
    @test Pioneer._mbr_temporal_corr_masked_hellinger_from_trace(
        trace,
        donor,
        UInt8(0x01),
    ) == 1.0f0
end

@testset "MBR temporal Hellinger compares every donor with the same recipient summary" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver_fold0.arrow")
        donor_path = joinpath(dir, "donor_fold0.arrow")
        frag1 = (1.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0)
        frag2 = (0.0f0, 1.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0)

        receiver = _mbr_empirical_main_table(
            precursor_idx = [1],
            scan_idx = [1001],
            ms_file_idx = 1,
            target = [true],
            frag1 = [100],
            frag2 = [0],
        )
        _mbr_add_integrated_fixture_columns!(receiver, [frag1], [frag1])
        receiver[!, Pioneer.MBR_INTEGRATED_APEX_IRT_COLUMN] = Float32[12.0]
        Arrow.write(receiver_path, receiver)
        Arrow.write(receiver_path * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [1],
            scan_idx = [1001],
            score = [0.40],
        ))

        donor = _mbr_empirical_main_table(
            precursor_idx = [1, 2],
            scan_idx = [2001, 2002],
            ms_file_idx = 2,
            target = [true, false],
            frag1 = [100, 0],
            frag2 = [0, 100],
        )
        _mbr_add_integrated_fixture_columns!(donor, [frag1, frag2], [frag1, frag2])
        donor[!, Pioneer.MBR_INTEGRATED_APEX_IRT_COLUMN] = Float32[8.0, 20.0]
        Arrow.write(donor_path, donor)
        Arrow.write(donor_path * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [1, 2],
            scan_idx = [2001, 2002],
            score = [0.90, 0.80],
        ))

        precursors = MBREmpiricalMockPrecursors(Float32[500, 501], Float32[10, 10.01])
        fragment_keys = Pioneer.build_mbr_fragment_annotation_keys(_mbr_empirical_fragment_lookup())
        partner_pools = Pioneer.build_counterfactual_partner_pools(
            [receiver_path, donor_path],
            precursors,
        )
        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1(
            [receiver_path, donor_path];
            passing_score_floor = 0.20f0,
            prefer_integrated_spectra = true,
            prefer_integrated_quant = true,
            require_integrated_irt = true,
        )
        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            fragment_keys;
            passing_score_floor = 0.20f0,
            use_integrated_mbr_features = true,
        )

        mbr = DataFrame(Arrow.Table(receiver_path * Pioneer.MBR_SIDECAR_SUFFIX))
        @test mbr.MBR_best_temporal_frag_hellinger_true[1] ≈ 0.0f0
        @test mbr.MBR_best_temporal_frag_hellinger_false[1] ≈ 1.0f0
        @test mbr.MBR_best_temporal_corr_frag_hellinger_true[1] ≈ 0.0f0
        @test mbr.MBR_best_temporal_corr_frag_hellinger_false[1] ≈ 1.0f0
        @test mbr.MBR_best_temporal_receiver_corr_frag_hellinger_true[1] ≈ 0.0f0
        @test mbr.MBR_best_temporal_receiver_corr_frag_hellinger_false[1] ≈ 1.0f0
        @test mbr.MBR_best_temporal_shared_corr_frag_hellinger_true[1] ≈ 0.0f0
        @test mbr.MBR_best_temporal_shared_corr_frag_hellinger_false[1] ≈ 1.0f0
    end
end

@testset "MBR correlated-fragment Hellinger uses the donor mask" begin
    recipient = (
        sqrt(0.5f0), sqrt(0.5f0), 0.0f0, 0.0f0,
        0.0f0, 0.0f0, 0.0f0, 0.0f0,
    )
    donor = (
        sqrt(0.5f0), 0.0f0, sqrt(0.5f0), 0.0f0,
        0.0f0, 0.0f0, 0.0f0, 0.0f0,
    )
    fragment_keys = Pioneer._MBRFragmentAnnotationKeys(UInt16[])

    all_rank_hellinger = Pioneer._mbr_smoothed_spectrum_hellinger_from_sqrt(
        recipient,
        donor,
        fragment_keys,
        UInt32(1),
        UInt32(2),
    )
    masked_hellinger = Pioneer._mbr_corr_masked_smoothed_spectrum_hellinger_from_sqrt(
        recipient,
        donor,
        UInt8(0x05),
    )

    expected = sqrt(0.5f0 * ((1.0f0 - sqrt(0.5f0))^2 + (0.0f0 - sqrt(0.5f0))^2))
    @test masked_hellinger ≈ expected
    @test masked_hellinger < all_rank_hellinger
    @test Pioneer._mbr_corr_masked_smoothed_spectrum_hellinger_from_sqrt(
        recipient,
        donor,
        UInt8(0x01),
    ) == 1.0f0
end

@testset "MBR receiver-correlated Hellinger uses the receiver mask" begin
    mktempdir() do dir
        receiver_path = joinpath(dir, "receiver_fold0.arrow")
        donor_path = joinpath(dir, "donor_fold0.arrow")

        Arrow.write(receiver_path, _mbr_empirical_main_table(
            precursor_idx = [1],
            scan_idx = [1001],
            ms_file_idx = 1,
            target = [true],
            frag_corr_bitvec = [0x03],
            frag_corr_rank = [42],
            frag1 = [100],
            frag2 = [100],
        ))
        Arrow.write(receiver_path * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [1],
            scan_idx = [1001],
            score = [0.40],
        ))

        Arrow.write(donor_path, _mbr_empirical_main_table(
            precursor_idx = [1],
            scan_idx = [2001],
            ms_file_idx = 2,
            target = [true],
            frag_corr_bitvec = [0x07],
            frag_corr_rank = [77],
            frag1 = [100],
            frag2 = [100],
            frag3 = [100],
        ))
        Arrow.write(donor_path * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [1],
            scan_idx = [2001],
            score = [0.90],
        ))

        precursors = MBREmpiricalMockPrecursors(Float32[500], Float32[10])
        fragment_keys = Pioneer.build_mbr_fragment_annotation_keys(_mbr_empirical_fragment_lookup())
        partner_pools = Pioneer.build_counterfactual_partner_pools([receiver_path, donor_path], precursors)
        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1(
            [receiver_path, donor_path];
            passing_score_floor = 0.20f0,
        )
        receiver_rank_table = zeros(UInt16, 256)
        receiver_rank_table[Int(0x03) + 1] = UInt16(66)

        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            receiver_path,
            donor_dict,
            partner_pools,
            fragment_keys,
            passing_score_floor = 0.20f0,
            bitvec_rank_tables_by_file = Dict(UInt32(1) => receiver_rank_table),
        )

        mbr = DataFrame(Arrow.Table(receiver_path * Pioneer.MBR_SIDECAR_SUFFIX))
        @test :MBR_best_temporal_frag_hellinger_true in propertynames(mbr)
        @test :MBR_best_receiver_corr_frag_hellinger_true in propertynames(mbr)
        @test :MBR_receiver_frag_corr_bitvec_rank in propertynames(mbr)
        @test :MBR_best_shared_corr_frag_hellinger_true in propertynames(mbr)
        @test :MBR_best_shared_corr_frag_bitvec_rank_true in propertynames(mbr)
        @test mbr.MBR_receiver_frag_corr_bitvec_rank[1] == 42.0f0
        @test isapprox(mbr.MBR_best_receiver_corr_frag_hellinger_true[1], 0.0f0; atol = 1.0f-6)
        @test mbr.MBR_best_corr_frag_hellinger_true[1] > 0.25f0
        @test isapprox(mbr.MBR_best_shared_corr_frag_hellinger_true[1], 0.0f0; atol = 1.0f-6)
        @test mbr.MBR_best_shared_corr_frag_bitvec_rank_true[1] == 66.0f0
    end
end

@testset "MBR empirical smoothed Hellinger aligns spectra by rank" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        f2 = joinpath(dir, "run2_fold0.arrow")
        f3 = joinpath(dir, "run3_fold0.arrow")

        Arrow.write(f1, _mbr_empirical_main_table(
            precursor_idx = [1],
            scan_idx = [1001],
            ms_file_idx = 1,
            target = [true],
            weight = [8],
            n_scans = [9],
            main_pep = [0.125],
            frag1 = [100],
            frag2 = [0],
        ))
        Arrow.write(f1 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [1],
            scan_idx = [1001],
            score = [0.40],
        ))

        Arrow.write(f2, _mbr_empirical_main_table(
            precursor_idx = [1, 2],
            scan_idx = [2000, 2001],
            ms_file_idx = 2,
            target = [true, false],
            weight = [8, 2],
            log2_intensity_explained = [1.0, 0.75],
            n_scans = [9, 5],
            library_hellinger = [0.10, 0.25],
            frag1 = [100, 0],
            frag2 = [0, 100],
        ))
        Arrow.write(f2 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [1, 2],
            scan_idx = [2000, 2001],
            score = [0.80, 0.91],
        ))

        Arrow.write(f3, _mbr_empirical_main_table(
            precursor_idx = [2],
            scan_idx = [3001],
            ms_file_idx = 3,
            target = [false],
            weight = [4],
            log2_intensity_explained = [0.25],
            n_scans = [2],
            library_hellinger = [0.75],
            frag1 = [100],
            frag2 = [0],
        ))
        Arrow.write(f3 * Pioneer.PASS1_SIDECAR_SUFFIX, _mbr_empirical_pass1_table(
            precursor_idx = [2],
            scan_idx = [3001],
            score = [0.30],
        ))

        precursors = MBREmpiricalMockPrecursors(
            Float32[500, 501],
            Float32[10, 10.01],
        )
        fragment_keys = Pioneer.build_mbr_fragment_annotation_keys(_mbr_empirical_fragment_lookup())
        partner_pools = Pioneer.build_counterfactual_partner_pools([f1, f2, f3], precursors)
        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1(
            [f1, f2, f3];
            passing_score_floor = 0.20f0,
        )

        Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            f1,
            donor_dict,
            partner_pools,
            fragment_keys,
            passing_score_floor = 0.20f0,
            lod_log2_weight_by_file = Dict(UInt32(1) => log2(4.0f0)),
            lod_log2_weight_global = log2(4.0f0),
        )
        mbr = DataFrame(Arrow.Table(f1 * Pioneer.MBR_SIDECAR_SUFFIX))

        @test isapprox(mbr.MBR_best_smoothed_frag_hellinger_false[1], 1.0f0; atol = 1.0f-6)
        @test :MBR_best_corr_frag_hellinger_true in propertynames(mbr)
        @test :MBR_best_corr_frag_hellinger_false in propertynames(mbr)
        @test :MBR_best_donor_frag_corr_bitvec_rank_true in propertynames(mbr)
        @test :MBR_best_donor_frag_corr_bitvec_rank_false in propertynames(mbr)
        @test :MBR_best_receiver_corr_frag_hellinger_true in propertynames(mbr)
        @test :MBR_best_receiver_corr_frag_hellinger_false in propertynames(mbr)
        @test :MBR_receiver_frag_corr_bitvec_rank in propertynames(mbr)
        @test :MBR_best_shared_corr_frag_hellinger_true in propertynames(mbr)
        @test :MBR_best_shared_corr_frag_hellinger_false in propertynames(mbr)
        @test :MBR_best_shared_corr_frag_bitvec_rank_true in propertynames(mbr)
        @test :MBR_best_shared_corr_frag_bitvec_rank_false in propertynames(mbr)
        @test mbr.MBR_best_donor_frag_corr_bitvec_rank_false[1] == 11.0f0
        @test mbr.MBR_receiver_frag_corr_bitvec_rank[1] == 11.0f0
        @test isapprox(mbr.MBR_worst_smoothed_frag_hellinger_false[1], 0.0f0; atol = 1.0f-6)
        @test mbr.MBR_best_is_missing_false[1] == false
        @test !(:MBR_max_pair_prob_false in propertynames(mbr))
        @test mbr.MBR_best_pair_prob_false[1] == 0.91f0
        @test mbr.MBR_worst_pair_prob_false[1] == 0.30f0
        @test mbr.MBR_log2_weight_lod_ratio[1] == 1.0f0
        @test mbr.MBR_best_log2_weight_ratio_false[1] == 2.0f0
        @test mbr.MBR_worst_log2_weight_ratio_false[1] == 1.0f0
        @test mbr.MBR_best_log2_explained_ratio_false[1] == 0.25f0
        @test mbr.MBR_worst_log2_explained_ratio_false[1] == 0.75f0
        @test mbr.MBR_best_abs_n_scans_diff_false[1] == 4.0f0
        @test mbr.MBR_worst_abs_n_scans_diff_false[1] == 7.0f0
        @test isapprox(mbr.MBR_best_log2_n_scans_ratio_false[1], log2(10.0f0 / 6.0f0); atol = 1.0f-6)
        @test mbr.MBR_best_donor_library_hellinger_false[1] == 0.25f0
        @test mbr.MBR_worst_donor_library_hellinger_false[1] == 0.75f0
        all_features = Set(vcat(Pioneer.FTR_FEATURES_F_TRUE, Pioneer.FTR_FEATURES_F_FALSE))
        @test !(:MBR_best_smoothed_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:MBR_best_smoothed_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE)
        @test :MBR_best_temporal_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_temporal_frag_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_best_temporal_frag_hellinger_rank_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_temporal_frag_hellinger_margin_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_temporal_corr_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_temporal_receiver_corr_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_temporal_shared_corr_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
        @test !(:MBR_worst_smoothed_frag_hellinger_true in all_features)
        @test !(:MBR_worst_smoothed_frag_hellinger_false in all_features)
        @test !(:main_search_prob in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:main_search_prob in Pioneer.FTR_FEATURES_F_FALSE)
        @test !(:trace_prob_infold in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:trace_prob_infold in Pioneer.FTR_FEATURES_F_FALSE)
        @test :MBR_best_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_worst_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_worst_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE
        @test !(:MBR_max_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:MBR_max_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE)
        @test :MBR_log2_weight_lod_ratio in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_log2_weight_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_log2_weight_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_worst_log2_weight_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_worst_log2_weight_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_worst_log2_explained_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_worst_log2_explained_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_worst_abs_n_scans_diff_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_worst_abs_n_scans_diff_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_worst_log2_n_scans_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_worst_log2_n_scans_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_best_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_worst_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_worst_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_worst_observed_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_worst_observed_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_single_donor_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_single_donor_false in Pioneer.FTR_FEATURES_F_FALSE
        @test !(:MBR_best_rt_diff_true in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:MBR_best_rt_diff_false in Pioneer.FTR_FEATURES_F_FALSE)
        @test !(:MBR_best_log_by_diff_true in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:MBR_best_log_by_diff_false in Pioneer.FTR_FEATURES_F_FALSE)
        @test :MBR_best_observed_irt_diff_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_observed_irt_diff_false in Pioneer.FTR_FEATURES_F_FALSE
        @test !(:MBR_best_donor_library_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:MBR_best_donor_library_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE)
        @test !(:MBR_worst_donor_library_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE)
        @test !(:MBR_worst_donor_library_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE)
    end
end
