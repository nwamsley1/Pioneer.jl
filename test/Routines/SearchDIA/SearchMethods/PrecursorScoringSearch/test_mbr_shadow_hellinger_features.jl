@testset "fragment annotation ids are added from library fragment metadata" begin
    frag_id(series::Integer, position::Integer, charge::Integer = 1) =
        UInt16(position | (series << 6) | (charge << 8))
    b5_id = frag_id(1, 5)
    y3_id = frag_id(2, 3)
    y3_z2_id = frag_id(2, 3, 2)

    fragments = Pioneer.CompactFrag{Float32}[
        Pioneer.CompactFrag(UInt32(1), 500.0f0, Float16(100), false, true, false, false,
                            UInt8(1), UInt8(5), UInt8(2), UInt8(1), UInt8(0)),
        Pioneer.CompactFrag(UInt32(1), 600.0f0, Float16(80), true, false, false, false,
                            UInt8(1), UInt8(3), UInt8(2), UInt8(2), UInt8(0)),
        Pioneer.CompactFrag(UInt32(2), 700.0f0, Float16(90), true, false, false, false,
                            UInt8(2), UInt8(3), UInt8(2), UInt8(1), UInt8(0)),
    ]
    lookup = Pioneer.StandardFragmentLookup{Float32}(fragments, UInt64[1, 3, 4])
    tbl = DataFrame(precursor_idx = UInt32[1, 2])

    Pioneer._add_fragment_annotation_id_columns!(tbl, lookup, fragments)

    @test tbl.frag1_id == UInt16[b5_id, y3_z2_id]
    @test tbl.frag2_id == UInt16[y3_id, 0]
    @test all(tbl[!, col] == UInt16[0, 0] for col in Pioneer.FRAGMENT_ANNOTATION_ID_COLUMNS[3:8])
end

@testset "MBR shadow Hellinger features are emitted for true and false donors" begin
    @test :empirical_frag_best_hellinger in Pioneer.ADVANCED_FEATURE_SET
    @test :MBR_shadow_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
    true_idx = findfirst(==(:MBR_shadow_hellinger_true), Pioneer.FTR_FEATURES_F_TRUE)
    @test true_idx !== nothing
    @test Pioneer.FTR_FEATURES_F_FALSE[true_idx] === :MBR_shadow_hellinger_false

    frag_id(series::Integer, position::Integer, charge::Integer = 1) =
        UInt16(position | (series << 6) | (charge << 8))
    b5_id = frag_id(1, 5)
    y3_id = frag_id(2, 3)

    mktempdir() do dir
        acceptor_path = joinpath(dir, "acceptor.arrow")
        donor_path = joinpath(dir, "donor.arrow")

        Arrow.write(acceptor_path, DataFrame(
            precursor_idx = UInt32[10],
            scan_idx = UInt32[100],
            ms_file_idx = UInt32[1],
            weight = Float32[10],
            log2_intensity_explained = Float16[1],
            irt_pred = Float32[50],
            irt_obs = Float32[48],
            log_by_ratio_m0 = Float16[0],
            rt = Float32[1],
            cv_fold = UInt8[0],
            target = Bool[true],
            shadow_frag1_int = Float32[100],
            shadow_frag2_int = Float32[0],
            shadow_frag3_int = Float32[0],
            shadow_frag4_int = Float32[0],
            shadow_frag5_int = Float32[0],
            shadow_frag6_int = Float32[0],
            shadow_frag7_int = Float32[0],
            shadow_frag8_int = Float32[0],
            frag1_id = UInt16[b5_id],
            frag2_id = UInt16[0],
            frag3_id = UInt16[0],
            frag4_id = UInt16[0],
            frag5_id = UInt16[0],
            frag6_id = UInt16[0],
            frag7_id = UInt16[0],
            frag8_id = UInt16[0],
        ))
        Arrow.write(acceptor_path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
            precursor_idx = UInt32[10],
            scan_idx = UInt32[100],
            trace_prob_prepass = Float32[0.4],
            trace_prob_infold = Float32[0.4],
        ))

        Arrow.write(donor_path, DataFrame(
            precursor_idx = UInt32[10, 20],
            scan_idx = UInt32[200, 201],
            ms_file_idx = UInt32[2, 2],
            weight = Float32[5, 5],
            log2_intensity_explained = Float16[1, 1],
            irt_pred = Float32[51, 51],
            irt_obs = Float32[49, 49],
            log_by_ratio_m0 = Float16[0, 0],
            rt = Float32[1, 1],
            cv_fold = UInt8[0, 0],
            target = Bool[true, true],
            shadow_frag1_int = Float32[100, 0],
            shadow_frag2_int = Float32[0, 100],
            shadow_frag3_int = Float32[0, 0],
            shadow_frag4_int = Float32[0, 0],
            shadow_frag5_int = Float32[0, 0],
            shadow_frag6_int = Float32[0, 0],
            shadow_frag7_int = Float32[0, 0],
            shadow_frag8_int = Float32[0, 0],
            frag1_id = UInt16[b5_id, y3_id],
            frag2_id = UInt16[0, b5_id],
            frag3_id = UInt16[0, 0],
            frag4_id = UInt16[0, 0],
            frag5_id = UInt16[0, 0],
            frag6_id = UInt16[0, 0],
            frag7_id = UInt16[0, 0],
            frag8_id = UInt16[0, 0],
        ))
        Arrow.write(donor_path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
            precursor_idx = UInt32[10, 20],
            scan_idx = UInt32[200, 201],
            trace_prob_prepass = Float32[0.95, 0.90],
            trace_prob_infold = Float32[0.95, 0.90],
        ))

        donor_dict = Pioneer.build_mbr_donor_dict_streaming_with_pass1([acceptor_path, donor_path])
        partner_col = zeros(UInt32, 20)
        partner_col[10] = UInt32(20)

        side_path = Pioneer.compute_mbr_features_per_file_to_sidecar_with_pass1!(
            acceptor_path,
            donor_dict,
            partner_col,
        )
        side = DataFrame(Arrow.Table(side_path))

        @test hasproperty(side, :MBR_shadow_hellinger_true)
        @test hasproperty(side, :MBR_shadow_hellinger_false)
        @test isapprox(side.MBR_shadow_hellinger_true[1], 0.0f0; atol = 1.0f-6)
        @test isapprox(side.MBR_shadow_hellinger_false[1], 0.0f0; atol = 1.0f-6)
    end
end
