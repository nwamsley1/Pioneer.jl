using Test, Arrow, DataFrames, Dictionaries, Random
import Pioneer

struct ProteinPoolTestPrecursors <: Pioneer.LibraryPrecursors
    pid_to_cv_fold::Vector{UInt8}
end

function protein_pool_test_rows(n)
    df = DataFrame(
        protein_name = ["P$i" for i in 1:n],
        target = isodd.(1:n),
        n_non_mbr_peptides = fill(Int64(2), n),
        species = fill("TEST", n),
        file_idx = cld.(1:n, 100),
        peptide_list = fill("output only", n),
    )
    for feature in Pioneer.run_level_protein_feature_names()
        df[!, feature] = zeros(Float32, n)
    end
    df.pg_score = ifelse.(df.target, 2.0f0, 0.1f0) .+ Float32.(1:n) .* 1.0f-7
    return df
end

function protein_pool_test_refs(dir, rows, ranges)
    return map(enumerate(ranges)) do (i, indices)
        path = joinpath(dir, "proteins_$i.arrow")
        Arrow.write(path, rows[indices, :])
        Pioneer.ProteinGroupFileReference(path)
    end
end

function protein_pool_test_folds(rows)
    return Dictionary(
        rows.protein_name,
        [(best_score = 0.0f0, cv_fold = UInt8((i - 1) ÷ 2 % 2))
         for i in 1:nrow(rows)],
    )
end

@testset "Bounded run-level protein training pool" begin
    @testset "sample identities, shadow labels, and source order" begin
        mktempdir() do dir
            actual = protein_pool_test_rows(96)
            refs = protein_pool_test_refs(dir, actual, [1:12, 1:0, 13:44, 45:96])
            pushfirst!(refs, Pioneer.ProteinGroupFileReference(joinpath(dir, "missing.arrow")))
            shadows = actual[1:2:31, :]
            shadows.pg_score .= -1.0f0
            original_shadows = copy(shadows)
            full = Pioneer.prepare_run_level_protein_training_rows(
                Pioneer.load_run_level_protein_training_rows(refs), shadows,
            )
            for cap in (112, 150)
                pool = Pioneer.load_run_level_protein_training_pool(refs, shadows, cap)
                @test isequal(pool.rows, full.rows)
                @test pool.n_shadows == 16
            end

            # Compare to a reservoir over the fully loaded population.
            cap = 32
            selected = collect(1:cap)
            rng = MersenneTwister(1776)
            for row in (cap + 1):nrow(full.rows)
                slot = rand(rng, 1:row)
                slot <= cap && (selected[slot] = row)
            end
            sort!(selected)
            pool = Pioneer.load_run_level_protein_training_pool(refs, shadows, cap)
            @test isequal(pool.rows, full.rows[selected, :])
            @test nrow(pool.rows) == cap
            @test pool.n_shadows == count(>(96), selected) > 0
            @test any(i -> 44 < i <= 96, selected)
            @test isequal(pool.rows, Pioneer.load_run_level_protein_training_pool(refs, shadows, cap).rows)
            @test !isequal(pool.rows, Pioneer.load_run_level_protein_training_pool(
                refs, shadows, cap; rng = MersenneTwister(22),
            ).rows)
            @test isequal(shadows, original_shadows)
            @test :peptide_list ∉ propertynames(pool.rows)
            @test :species ∉ propertynames(pool.rows)
            @test :file_idx ∉ propertynames(pool.rows)

            qc_pool = Pioneer.load_run_level_protein_training_pool(
                refs, shadows, cap; include_qc_plot_columns = true,
            )
            qc_full = Pioneer.prepare_run_level_protein_training_rows(
                Pioneer.load_run_level_protein_training_rows(refs; include_qc_plot_columns = true),
                shadows,
            )
            @test isequal(qc_pool.rows, qc_full.rows[selected, :])
            folds = protein_pool_test_folds(actual)
            Pioneer.assign_protein_group_cv_folds!(pool.rows, folds)
            @test pool.rows.cv_fold == [folds[name].cv_fold for name in pool.rows.protein_name]
            @test all(pool.rows.target[pool.rows.protein_shadow_negative])
            @test !any(pool.rows.training_label[pool.rows.protein_shadow_negative])
            shadow_only = Pioneer.load_run_level_protein_training_pool(
                refs, shadows, 1; rng = MersenneTwister(17),
            )
            @test shadow_only.n_shadows == 1
            @test isequal(shadow_only.rows, full.rows[106:106, :])
        end
    end

    @testset "models trained on the pool score every source row" begin
        mktempdir() do dir
            actual = protein_pool_test_rows(1200)
            refs = protein_pool_test_refs(dir, actual, [1:300, 301:1200])
            pool = Pioneer.load_run_level_protein_training_pool(refs, DataFrame(), 320).rows
            pool_names = copy(pool.protein_name)
            folds = protein_pool_test_folds(actual)
            model = Pioneer.fit_protein_lightgbm_semisupervised(
                select(pool, :pg_score), pool.training_label, pool.pg_score,
                pool.precursor_consensus_prefix_shape, pool.n_non_mbr_peptides;
                n_iterations = 1,
            )
            @test model.booster !== nothing
            @test isempty(model.booster.booster.datasets)
            predictions = Pioneer.lightgbm_predict(model, select(pool, :pg_score))
            @test minimum(predictions[pool.training_label]) > maximum(predictions[.!pool.training_label])
            @test Pioneer.perform_protein_scoring_multifold(
                pool, dir, refs, ProteinPoolTestPrecursors(UInt8[0, 1]);
                protein_to_cv_fold = folds, write_qc_plots = false,
            )
            scored = vcat([DataFrame(Arrow.Table(ref.file_path)) for ref in refs]...)
            @test nrow(scored) == nrow(actual)
            @test Set(scored.protein_name) == Set(actual.protein_name)
            @test all(score -> 0.0f0 < score < 1.0f0, scored.pg_score)
            @test all(scored.peptide_list .== "output only")
            @test Set(pool_names) ⊊ Set(scored.protein_name)
            by_name = Dict(scored.protein_name .=> scored.pg_score)
            @test pool.pg_score ≈ [by_name[name] for name in pool_names]
            original_scores = Dict(actual.protein_name .=> actual.pg_score)
            @test scored.old_pg_score == [original_scores[name] for name in scored.protein_name]
            @test :training_label ∉ propertynames(scored)
            @test :protein_shadow_negative ∉ propertynames(scored)
            @test :cv_fold ∉ propertynames(scored)
        end
    end

    @testset "production path exceeds the previous limit and retains fallback" begin
        mktempdir() do dir
            actual = protein_pool_test_rows(100_001)
            actual.protein_name .= "P1"
            actual.target .= true
            actual.pg_score .= 2.0f0
            refs = protein_pool_test_refs(dir, actual, [1:50_000, 50_001:100_001])
            folds = protein_pool_test_folds(actual[1:1, :])
            @test !Pioneer.perform_run_level_protein_scoring(
                refs, 100_000, dir, ProteinPoolTestPrecursors(UInt8[0, 1]);
                protein_to_cv_fold = folds, write_qc_plots = false,
            )
            for ref in refs
                scored = DataFrame(Arrow.Table(ref.file_path))
                @test nrow(scored) == ref.row_count
                @test all(scored.pg_score .≈ 1.0f0 - exp(-2.0f0))
                @test all(scored.old_pg_score .== 2.0f0)
            end
            @test sum(ref.row_count for ref in refs) == 100_001
        end
    end
end
