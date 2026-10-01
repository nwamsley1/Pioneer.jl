using Test, Random, Arrow, DataFrames
using Pioneer

@testset "Grouped score statistics and bounded sorting" begin
    function reference_statistics(scores, labels, scale)
        keys = sort!(unique(Pioneer._score_key.(scores)); rev=true)
        ts = [count(i -> labels[i] && Pioneer._score_key(scores[i]) == key, eachindex(scores)) for key in keys]
        ds = [count(i -> !labels[i] && Pioneer._score_key(scores[i]) == key, eachindex(scores)) for key in keys]
        cumulative_t, cumulative_d = cumsum(ts), cumsum(ds)
        raw = [t == 0 ? Inf32 : Float32(d * Float64(scale) / t) for (t,d) in zip(cumulative_t,cumulative_d)]
        q = [minimum(raw[i:end]) for i in eachindex(raw)]
        weights = ts .+ scale .* ds
        probabilities = Pioneer._weighted_pava(vcat(0.0, scale .* ds ./ weights), vcat(1.0, weights))[2:end]
        pep = clamp.(probabilities ./ (1 .- probabilities), 0, 1)
        index = Dict(k => i for (i,k) in enumerate(keys))
        return Float32[q[index[Pioneer._score_key(s)]] for s in scores],
               Float32[pep[index[Pioneer._score_key(s)]] for s in scores]
    end

    rng = MersenneTwister(123)
    cases = [
        (Float32[3,3,2,1], Bool[1,0,1,1]),
        (fill(0.5f0, 301), [isodd(i) for i in 1:301]),
        (rand(rng, Float32[0.1,0.2,0.4,0.6,0.7,0.9], 2000), rand(rng, Bool, 2000)),
        (rand(rng, Float32, 1500), rand(rng, Bool, 1500)),
        (Float32[Inf,0.0,-0.0,NaN,-Inf,Inf], Bool[1,1,0,1,0,0]),
        (Float32.(1:201), trues(201)),
        (Float32.(1:201), falses(201)),
    ]
    for (scores, labels) in cases, scale in (0.5f0, 1.0f0, 2.0f0)
        expected_q, expected_pep = reference_statistics(scores, labels, scale)
        for budget in (4096, 8192, 1024^2)
            q, pep = similar(scores), similar(scores)
            Pioneer.get_score_statistics!(scores, labels, q, pep;
                memory_budget_bytes=budget, fdr_scale_factor=scale)
            @test q == expected_q
            @test pep ≈ expected_pep
            perm = randperm(rng, length(scores))
            q2, p2 = similar(q), similar(pep)
            Pioneer.get_score_statistics!(scores[perm], labels[perm], q2, p2;
                memory_budget_bytes=budget, fdr_scale_factor=scale)
            @test q2 == q[perm]
            @test p2 ≈ pep[perm]
            producer = emit -> Pioneer._emit_score_arrays(emit, scores, labels)
            for threshold in (0.0f0, 0.01f0, 0.5f0, 1.0f0)
                selected = [Pioneer._score_key(scores[i]) for i in eachindex(scores) if labels[i] && q[i] <= threshold]
                floor = Pioneer.qvalue_score_cutoff(producer; q_threshold=threshold,
                    memory_budget_bytes=budget, max_fanin=2, fdr_scale_factor=scale)
                @test floor === (isempty(selected) ? nothing : minimum(selected))
            end
        end
    end

    @testset "PAVA stack spills and merges across buffers" begin
        # Increasing decoy fractions build a deep stack; trailing targets force
        # it to merge back through multiple on-disk buffers.
        scores = Float32[]; labels = Bool[]
        for i in 1:400
            append!(scores, fill(Float32(401-i), 20))
            append!(labels, [j > (i ÷ 20) for j in 1:20])
        end
        append!(scores, fill(0.0f0, 8000)); append!(labels, trues(8000))
        q = similar(scores); pep = similar(scores)
        expected_q, expected_pep = reference_statistics(scores, labels, 1.0f0)
        Pioneer.get_score_statistics!(scores, labels, q, pep; memory_budget_bytes=4096)
        @test q == expected_q
        @test pep ≈ expected_pep
    end

    @testset "Compaction and temporary-file lifetime" begin
        fit = Pioneer.build_score_calibration(; memory_budget_bytes=4096, max_fanin=2) do emit
            for i in 1:10_000
                emit(Float32(i), true)
            end
        end
        @test fit.qval_spline.store.n == 2
        @test all(iszero, fit.qval_spline.([0,1,20,9999,10000,10001]))
        @test all(iszero, fit.pep_interp.([0,1,20,9999,10000,10001]))
        @test length(fit.qval_spline.store.cache) <= fit.qval_spline.store.max_pages
        path = fit.qval_spline.store.path
        Pioneer._close_score_calibration(fit)
        @test !ispath(dirname(path))
        @test isnothing(Pioneer.build_score_calibration(emit -> nothing; memory_budget_bytes=4096))
        @test isnothing(Pioneer.qvalue_score_cutoff(emit -> nothing))
        @test_throws ErrorException Pioneer.build_score_calibration(; memory_budget_bytes=4096) do emit
            for i in 1:1000
                emit(Float32(i), true)
            end
            error("Interrupted producer")
        end
    end

    @testset "Arrow batches, protein refs, and MBR initial-pass filter" begin
        mktempdir() do dir
            scores = repeat(Float32[0.9,0.9,0.8,0.7], 40)
            labels = repeat(Bool[1,0,1,1], 40)
            path = joinpath(dir, "scores.arrow")
            table = DataFrame(prec_prob=scores, trace_prob_prepass=scores, target=labels,
                qval=fill(0.01f0, length(scores)), global_qval=fill(0.01f0,length(scores)))
            open(Arrow.Writer, path; file=true) do writer
                Arrow.write(writer, table[1:70,:])
                Arrow.write(writer, table[71:end,:])
            end
            q, pep = reference_statistics(scores, labels, 1.0f0)
            for refs in ([Pioneer.PSMFileReference(path)], [Pioneer.ProteinGroupFileReference(path)])
                fit = Pioneer.build_qvalue_spline_from_refs(refs, :prec_prob, joinpath(dir,"unused.arrow");
                    compute_pep=true, memory_budget_bytes=4096)
                @test Float32.(fit.qval_spline.(scores)) == q
                @test Float32.(fit.pep_interp.(scores)) ≈ pep
                Pioneer._close_score_calibration(fit)
            end
            for initial_pass in (false, true)
                @test Pioneer._mbr_donor_score_floor([path]; donor_q_threshold=0.5f0,
                    require_initial_pass=initial_pass) == 0.7f0
                @test Pioneer._mbr_donor_score_floor([path]; donor_q_threshold=0.01f0,
                    require_initial_pass=initial_pass) == Inf32
            end
        end
    end
end

@testset "Calibration cache eviction and failure cleanup" begin
    scores = Float32.(5000:-1:1)
    labels = BitVector(i <= 1000 for i in 1:5000)
    q = similar(scores); pep = similar(scores)
    Pioneer.get_score_statistics!(scores, labels, q, pep)
    fit = Pioneer.build_score_calibration(emit -> Pioneer._emit_score_arrays(emit, scores, labels);
        memory_budget_bytes=4096, max_fanin=2)
    actual_q, actual_pep = similar(q), similar(pep)
    Threads.@threads for i in eachindex(scores)
        actual_q[i] = fit.qval_spline(scores[i])
        actual_pep[i] = fit.pep_interp(scores[i])
    end
    @test actual_q == q
    @test actual_pep == pep
    @test fit.qval_spline.store.n > fit.qval_spline.store.page_size * fit.qval_spline.store.max_pages
    @test length(fit.qval_spline.store.cache) <= fit.qval_spline.store.max_pages
    Pioneer._close_score_calibration(fit)
    mktempdir() do directory
        @test_throws ErrorException Pioneer.build_score_calibration(; memory_budget_bytes=4096, temp_parent=directory) do emit
            for i in 1:1000
                emit(Float32(i), true)
            end
            error("Interrupted after spilling")
        end
        @test isempty(readdir(directory))
    end
end

@testset "MBR stops with valid statistics when tied scores yield no positives" begin
    frame = DataFrame(target=trues(4), cv_fold=UInt8[0,1,0,1], constant=zeros(Float32,4))
    for i in 1:Pioneer.MBR_N_COUNTERFACTUALS
        frame[!, Pioneer._mbr_missing_feature(i)] = falses(4)
    end
    state, _ = Pioneer._mbr_semisupervised_oof(frame, [:constant],
        [[:constant] for _ in 1:Pioneer.MBR_N_COUNTERFACTUALS])
    @test state.iteration == 1
    @test state.metrics.n_positive == 0
    @test all(state.metrics.eval_mask[1:4])
    @test state.metrics.qvalues[1:4] == ones(Float32, 4)
    @test state.metrics.peps[1:4] == ones(Float32, 4)
end

@testset "Protein annotation pipelines accept bounded calibration mappings" begin
    mktempdir() do dir
        scores = Float32[0.95, 0.9, 0.8, 0.7]
        targets = Bool[1, 1, 0, 1]
        table = DataFrame(protein_name=["p1", "p2", "p3", "p4"],
            target=targets, entrap_id=zeros(UInt8, 4), pg_score=scores)
        path = joinpath(dir, "proteins.arrow")
        Arrow.write(path, table)
        refs = [Pioneer.ProteinGroupFileReference(path)]
        fit = Pioneer.build_qvalue_spline_from_refs(refs, :pg_score, joinpath(dir, "sorted.arrow");
            compute_pep=true, memory_budget_bytes=4096)
        try
            qvals = Float32.(fit.qval_spline.(scores))
            peps = Float32.(fit.pep_interp.(scores))
            global_scores = Dict((table.protein_name[i], targets[i], UInt8(0)) => scores[i] for i in 1:4)
            global_qvals = Pioneer.build_protein_global_qval_dict(global_scores)
            pipeline = Pioneer.TransformPipeline() |>
                Pioneer.add_dict_column_composite_key(:global_pg_score, [:protein_name, :target, :entrap_id], global_scores) |>
                Pioneer.add_dict_column_composite_key(:global_pg_qval, [:protein_name, :target, :entrap_id], global_qvals) |>
                Pioneer.add_interpolated_column(:pg_qval, :pg_score, fit.qval_spline) |>
                Pioneer.add_interpolated_column(:pg_pep, :pg_score, fit.pep_interp) |>
                Pioneer.filter_by_multiple_thresholds([(:global_pg_qval, 0.01f0), (:pg_qval, 0.01f0)])
            output = Pioneer.apply_pipeline_batch(refs, pipeline, joinpath(dir, "passing"))
            actual = DataFrame(Arrow.Table(Pioneer.file_path(only(output))))
            keep = qvals .<= 0.01f0
            @test actual.protein_name == table.protein_name[keep]
            @test actual.pg_qval == qvals[keep]
            @test actual.pg_pep == peps[keep]
            @test eltype(actual.pg_qval) == Float32
            @test eltype(actual.pg_pep) == Float32
            # The subsequent q-value-only recalibration uses the same helper.
            recalibrated = Pioneer.build_qvalue_spline_from_refs(output, :pg_score, joinpath(dir, "recalc.arrow");
                memory_budget_bytes=4096)
            try
                recalc_pipeline = Pioneer.TransformPipeline() |>
                    Pioneer.add_interpolated_column(:pg_qval, :pg_score, recalibrated.qval_spline)
                recalc_output = Pioneer.apply_pipeline_batch(output, recalc_pipeline, joinpath(dir, "recalibrated"))
                @test DataFrame(Arrow.Table(Pioneer.file_path(only(recalc_output)))).pg_qval ==
                    Float32.(recalibrated.qval_spline.(actual.pg_score))
            finally
                Pioneer._close_score_calibration(recalibrated)
            end
        finally
            Pioneer._close_score_calibration(fit)
        end
    end
end

@testset "Interpolation column helper retains legacy and empty-table behavior" begin
    interp = Pioneer.linear_interpolation(Float32[0, 1], Float32[1, 0];
        extrapolation_bc=Pioneer.Interpolations.Flat())
    operation = last(Pioneer.add_interpolated_column(:qval, :score, interp))
    frame = DataFrame(score=Float32[-1, 0.25, 2])
    @test operation(frame).qval == Float32[1, 0.75, 0]
    empty_frame = DataFrame(score=Float32[])
    @test isempty(operation(empty_frame).qval)
    @test eltype(empty_frame.qval) == Float32
end

@testset "Bulk score assignment across memory and disk orders" begin
    rng = MersenneTwister(877)
    scores = rand(rng, Float64[Inf, 0.9, 0.5, 0.0, -0.0, -1, -Inf, NaN], 1001)
    labels = rand(rng, Bool, length(scores))
    for scale in (0.5f0, 1f0, 2f0)
        expected_q = zeros(Float32, length(scores))
        expected_pep = similar(expected_q)
        Pioneer.get_score_statistics!(scores, labels, expected_q, expected_pep;
            fdr_scale_factor=scale, memory_budget_bytes=1_000_000)
        q, pep = similar(expected_q), similar(expected_pep)
        Pioneer.get_score_statistics!(scores, labels, q, pep;
            fdr_scale_factor=scale, memory_budget_bytes=4096)
        @test q == expected_q
        @test pep == expected_pep
        Pioneer.get_qvalues!(scores, labels, q; fdr_scale_factor=scale, memory_budget_bytes=4096)
        Pioneer.get_PEP!(scores, labels, pep; fdr_scale_factor=scale, memory_budget_bytes=4096)
        @test q == expected_q
        @test pep == expected_pep
        order = Pioneer._score_order(scores)
        Pioneer._get_PEP_from_order!(scores, labels, pep, Int32.(order), scale; memory_budget_bytes=4096)
        @test pep == expected_pep
        Pioneer.get_qvalues!(scores[order], labels[order], q; doSort=false,
            fdr_scale_factor=scale, memory_budget_bytes=4096)
        Pioneer.get_PEP!(scores[order], labels[order], pep; doSort=false,
            fdr_scale_factor=scale, memory_budget_bytes=4096)
        @test q == expected_q[order]
        @test pep == expected_pep[order]
    end
    mktempdir() do dir
        Pioneer.with_score_order(scores; memory_budget_bytes=4096, max_fanin=2, temp_parent=dir) do order
            @test order isa Pioneer.DiskScoreOrder
            indices = collect(order)
            @test sort(indices) == collect(eachindex(scores))
            @test issorted(Pioneer._score_key.(scores[indices]); rev=true)
            @test [order[i] for i in length(order):-1:1] == reverse(indices)
        end
        @test isempty(readdir(dir))
        @test_throws ErrorException Pioneer.with_score_order(_ -> error("consumer failure"), scores;
            memory_budget_bytes=4096, max_fanin=2, temp_parent=dir)
        @test isempty(readdir(dir))
    end
    # A monotonic sequence of group decoy fractions keeps enough PAVA blocks
    # to spill, including leading zero blocks and the unit pseudocount.
    scores = Float32.(2000:-1:1)
    for labels in (falses(2000), trues(2000), vcat(trues(1000), falses(1000)), vcat(falses(1000), trues(1000)))
        expected = zeros(Float32, 2000)
        actual = similar(expected)
        Pioneer._get_PEP_from_order!(scores, labels, expected, 1:2000, 1f0; memory_budget_bytes=1_000_000)
        Pioneer._get_PEP_from_order!(scores, labels, actual, 1:2000, 1f0; memory_budget_bytes=4096)
        @test actual == expected
    end
    # Exercise the cutoff on either side with the same grouped results.
    for n in (63, 64, 65)
        scores = fill(0.5f0, n)
        labels = [isodd(i) for i in 1:n]
        q, pep = zeros(Float32, n), zeros(Float32, n)
        Pioneer.get_score_statistics!(scores, labels, q, pep; memory_budget_bytes=4096)
        @test all(==(Float32(count(!, labels) / count(identity, labels))), q)
        @test all(==(Float32(count(!, labels) / count(identity, labels))), pep)
    end
end
