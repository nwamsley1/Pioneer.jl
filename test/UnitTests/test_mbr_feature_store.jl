using Test, Pioneer, DataFrames, Arrow, Random

function mbr_store_fixture(n=1200)
    rng = MersenneTwister(83)
    frame = DataFrame(precursor_idx=UInt32.(1:n), scan_idx=UInt32.(1001:1000+n),
        ms_file_idx=ones(UInt32, n), cv_fold=UInt8.(mod.(1:n, 2)),
        target=mod.(1:n, 3) .!= 0, qval=fill(0.2f0, n), global_qval=fill(0.005f0, n),
        trace_prob_prepass=rand(rng, Float32, n), trace_prob_infold=rand(rng, Float32, n),
        MBR_best_is_missing_true=falses(n))
    for feature in Pioneer.MBR_RECEIVER_FEATURES
        frame[!, feature] = rand(rng, Float32, n)
    end
    for stem in Pioneer.MBR_PAIRED_FEATURE_STEMS
        frame[!, Pioneer._mbr_true_feature(stem)] = rand(rng, Float32, n)
        for k in 1:3
            frame[!, Pioneer._mbr_false_feature(stem, k)] = rand(rng, Float32, n)
        end
    end
    for k in 1:3
        frame[!, Pioneer._mbr_missing_feature(k)] = BitVector(mod.(1:n, k+2) .== 0)
    end
    return frame
end

@testset "Bounded MBR feature store" begin
    frame = mbr_store_fixture()
    Pioneer._mbr_add_hellinger_contrasts!(frame)
    tf, ff = Pioneer._mbr_available_feature_sets(frame)
    rows = vcat([1, 1200, 1201, 4800, 499, 500, 501, 1], rand(MersenneTwister(3), 1:4800, 200))
    expected = Pioneer._mbr_gather_feature_rows(frame, tf, ff, rows, 1200)
    mktempdir() do dir
        for budget in (64*1024^2, 1)
            store = Pioneer._MBRFeatureStore(joinpath(dir, "features.bin"); budget)
            try
                for first in 1:137:nrow(frame)
                    Pioneer._mbr_store_features!(store, frame[first:min(first+136, nrow(frame)), :])
                end
                @test (store.data.io === nothing) == (budget > 1)
                @test store.data.bytes <= budget
                @test readdir(dir) == (budget > 1 ? String[] : ["features.bin"])
                @test Pioneer._mbr_available_feature_sets(store.schema) == (tf, ff)
                @test Pioneer._mbr_gather_feature_rows(store, tf, ff, rows, 1200) == expected
                @test size(Pioneer._mbr_gather_feature_rows(store, tf, ff, Int[], 1200)) == (0, length(tf))
                @test_throws BoundsError Pioneer._mbr_gather_feature_rows(store, tf, ff, [4801], 1200)
                budget == 1 && @test isempty(store.data.blocks)
            finally
                close(store)
            end
        end
    end
end

@testset "Streamed MBR candidate loading and scoring" begin
    mktempdir() do dir
        frame = mbr_store_fixture()
        paths = String[]
        for (i, rows) in enumerate((1:503, 504:999, 1000:1200))
            part = frame[rows, :]
            part.ms_file_idx .= UInt32(i)
            part.qval[1:2] .= 0.005f0
            part.target[1:2] .= true
            path = joinpath(dir, "run_$i.arrow")
            push!(paths, path)
            Arrow.write(path, part)
            Arrow.write(path * Pioneer.PASS1_SIDECAR_SUFFIX,
                select(part, :precursor_idx, :scan_idx, :trace_prob_prepass, :trace_prob_infold))
            sidecar = select(part, [:precursor_idx, :scan_idx,
                filter(c -> startswith(string(c), "MBR_"), propertynames(part))...])
            # Exercise both dense and indexed sidecars.
            if i != 2
                sidecar = sidecar[3:end, :]
                sidecar[!, :row_idx] = collect(3:nrow(part))
            end
            Arrow.write(path * Pioneer.MBR_SIDECAR_SUFFIX, sidecar)
        end
        dense = Pioneer.load_postintegration_mbr_candidates(paths, 0.01f0)
        reference = copy(dense.candidates)
        counts = (dense.base_targets, dense.base_decoys)
        summary = Pioneer.apply_postintegration_mbr_rescoring!(reference;
            alpha=0.01f0, q_value_threshold=0.01f0, baseline_counts=counts, frame_is_candidates=true)
        for budget in (64*1024^2, 1)
            mktempdir(dir) do scratch
                store = Pioneer._MBRFeatureStore(joinpath(scratch, "features.bin"); budget)
                try
                    loaded = Pioneer.load_postintegration_mbr_candidates(paths, 0.01f0;
                        feature_store=store, feature_batch_size=97)
                    @test loaded.masks == dense.masks
                    @test loaded.n_rows == dense.n_rows
                    @test (loaded.base_targets, loaded.base_decoys) == counts
                    @test !hasproperty(loaded.candidates, :trace_prob_infold)
                    @test nrow(loaded.candidates) == nrow(reference)
                    result = Pioneer.apply_postintegration_mbr_rescoring!(loaded.candidates;
                        alpha=0.01f0, q_value_threshold=0.01f0, baseline_counts=counts,
                        frame_is_candidates=true, feature_source=store)
                    @test isequal(result, summary)
                    for col in propertynames(loaded.candidates)
                        @test isequal(loaded.candidates[!, col], reference[!, col])
                    end
                    @test Pioneer._write_mbr_recovery_sidecars_from_candidates!(
                        loaded.candidates, loaded.masks, loaded.n_rows, paths) == 3
                    @test readdir(scratch) == (budget > 1 ? String[] : ["features.bin"])
                finally
                    close(store)
                end
            end
        end
    end
end
