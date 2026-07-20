# Unit tests for the two-round precursor cross-run scoring machinery in
# src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/two_round_scoring.jl
#
# All functions under test are module-internal (not exported) and are referenced
# as `Pioneer._name`. Fixtures are tiny synthetic fold Arrow files + matching
# `.pass1_sidecar.arrow` sidecars, exercising the real inject/remove/group-feature
# code paths. Randomness in the clustering path uses fixed internal seeds, so the
# tests assert on partition structure and remain deterministic.

using Test
using Pioneer
using DataFrames
using Arrow

const _P = Pioneer

@testset "two-round group scoring" begin

    # ---------------------------------------------------------------
    # 1. _group_size schedule
    # ---------------------------------------------------------------
    @testset "_group_size schedule" begin
        @test _P._group_size(0) == 0
        @test _P._group_size(8) == 0        # R<9 -> cluster off
        @test _P._group_size(9) == 4        # min(9, 9÷2=4)
        @test _P._group_size(17) == 8       # 17÷2 = 8
        @test _P._group_size(18) == 9       # 18÷2 = 9, caps at 9
        @test _P._group_size(40) == 9       # capped
        @test _P._group_size(100) == 9      # capped
    end

    # ---------------------------------------------------------------
    # 2. _pos_logit
    # ---------------------------------------------------------------
    @testset "_pos_logit" begin
        @test _P._pos_logit(0.5f0) == 0.0f0
        @test _P._pos_logit(0.3f0) == 0.0f0          # x<=0.5 -> 0
        @test _P._pos_logit(0.9f0) > 0.0f0
        @test _P._pos_logit(0.9f0) ≈ log(0.9f0 / 0.1f0) atol=1e-4
        # x = 1.0 must be finite (clamped away from the asymptote).
        v = _P._pos_logit(1.0f0)
        @test isfinite(v)
        @test v > 0.0f0
        # monotone increasing above 0.5
        @test _P._pos_logit(0.8f0) < _P._pos_logit(0.95f0)
    end

    # ---------------------------------------------------------------
    # 3. TopK4 + leave-one-out helpers
    # ---------------------------------------------------------------
    @testset "TopK4 leave-one-out" begin
        # EMPTY sentinel
        @test _P.EMPTY_TOPK4.v == (0f0, 0f0, 0f0, 0f0)
        @test _P.EMPTY_TOPK4.r == (UInt32(0), UInt32(0), UInt32(0), UInt32(0))

        # Insert distinct (score, run) pairs; runs 1..5 with increasing scores.
        t = _P.EMPTY_TOPK4
        t = _P._tk_insert(t, 0.60f0, UInt32(1))
        t = _P._tk_insert(t, 0.99f0, UInt32(2))   # top
        t = _P._tk_insert(t, 0.90f0, UInt32(3))
        t = _P._tk_insert(t, 0.70f0, UInt32(4))
        t = _P._tk_insert(t, 0.80f0, UInt32(5))
        # Only the top-4 scores survive (0.60 from run 1 is dropped).
        @test t.v == (0.99f0, 0.90f0, 0.80f0, 0.70f0)
        @test t.r == (UInt32(2), UInt32(3), UInt32(5), UInt32(4))

        # loo_max excluding the top run (run 2) drops it -> 0.90
        @test _P._tk_loo_max(t, UInt32(2)) ≈ 0.90f0
        # loo_max excluding a non-top run keeps the top -> 0.99
        @test _P._tk_loo_max(t, UInt32(3)) ≈ 0.99f0

        # loo_top3: sum of positive logits of top-3 scores excluding `own`.
        # Excluding run 2 (0.99): top-3 of remaining = 0.90, 0.80, 0.70.
        expect = _P._pos_logit(0.90f0) + _P._pos_logit(0.80f0) + _P._pos_logit(0.70f0)
        @test _P._tk_loo_top3(t, UInt32(2)) ≈ expect atol=1e-5
        # Excluding a run not in the top-3 window still only sums 3 entries.
        # own = run 4 (0.70, the 4th): top-3 = 0.99, 0.90, 0.80.
        expect2 = _P._pos_logit(0.99f0) + _P._pos_logit(0.90f0) + _P._pos_logit(0.80f0)
        @test _P._tk_loo_top3(t, UInt32(4)) ≈ expect2 atol=1e-5

        # Scores <= 0.5 contribute 0 logit; a record with low scores sums to 0.
        tlow = _P.EMPTY_TOPK4
        tlow = _P._tk_insert(tlow, 0.40f0, UInt32(1))
        tlow = _P._tk_insert(tlow, 0.30f0, UInt32(2))
        @test _P._tk_loo_top3(tlow, UInt32(99)) == 0.0f0
    end

    # ---------------------------------------------------------------
    # 4. _cluster_runs
    # ---------------------------------------------------------------
    @testset "_cluster_runs" begin
        # (a) k<=0 -> single cluster of all ones, length R.
        rp_small = [UInt32[1,2,3] for _ in 1:6]
        lab0 = _P._cluster_runs(rp_small, 6, 0)
        @test lab0 == ones(Int, 6)
        @test length(lab0) == 6

        # (b) Two obvious disjoint blocks. Runs 1..10 share precursor set A
        #     (ids 1..50); runs 11..20 share disjoint set B (ids 51..100).
        R = 20
        A = collect(UInt32, 1:50)
        B = collect(UInt32, 51:100)
        run_present = Vector{Vector{UInt32}}(undef, R)
        for r in 1:10;  run_present[r] = copy(A); end
        for r in 11:20; run_present[r] = copy(B); end
        # k=10 so neither block (size 10) is split by _split_oversized!.
        lab = _P._cluster_runs(run_present, R, 10)
        @test length(lab) == R
        # Partition structure: every run in block A shares one label, every run
        # in block B shares another, and the two labels differ.
        blockA = lab[1:10]
        blockB = lab[11:20]
        @test all(==(blockA[1]), blockA)
        @test all(==(blockB[1]), blockB)
        @test blockA[1] != blockB[1]
    end

    # ---------------------------------------------------------------
    # 5. _nearest_sorted_idx binary search
    # ---------------------------------------------------------------
    @testset "_nearest_sorted_idx" begin
        empty_v = Float32[]
        @test _P._nearest_sorted_idx(empty_v, 1.0f0) == 0

        v = Float32[10.0, 20.0, 30.0, 40.0]
        # exact hits
        @test _P._nearest_sorted_idx(v, 10.0f0) == 1
        @test _P._nearest_sorted_idx(v, 30.0f0) == 3
        @test _P._nearest_sorted_idx(v, 40.0f0) == 4
        # below-min -> 1, above-max -> n
        @test _P._nearest_sorted_idx(v, 5.0f0) == 1
        @test _P._nearest_sorted_idx(v, 100.0f0) == 4
        # nearer to lower / upper neighbor
        @test _P._nearest_sorted_idx(v, 21.0f0) == 2
        @test _P._nearest_sorted_idx(v, 29.0f0) == 3
        # exact midpoint -> tie breaks to the LOWER index
        @test _P._nearest_sorted_idx(v, 25.0f0) == 2
    end

    # ---------------------------------------------------------------
    # Helpers for the Arrow-file fixtures (#6 and #7)
    # ---------------------------------------------------------------
    # Minimal fold-file with the columns inject_shadow_decoys! reads/grafts.
    function _write_fixture!(dir, name, df::DataFrame, s1::Vector{Float32})
        path = joinpath(dir, name)
        _P.writeArrow(path, df)
        side = DataFrame(; (_P.SHADOW_SIDECAR_SCORE_COL => s1,)...)
        _P.writeArrow(path * _P.PASS1_SIDECAR_SUFFIX, side)
        return path
    end

    # ---------------------------------------------------------------
    # 6. inject_shadow_decoys! / remove_shadow_decoys! round-trip
    # ---------------------------------------------------------------
    @testset "shadow decoy inject/remove round-trip" begin
        mktempdir() do dir
            # 4 rows: 2 targets (irt 10,20), 2 decoys (irt 10.5,20.5).
            # global_max is unique per row so we can verify the graft is the
            # ORIGINAL row's value (multiset preserved across shadows).
            n = 4
            df = DataFrame(
                precursor_idx = UInt32[1, 2, 3, 4],
                target        = Bool[true, true, false, false],
                irt_obs       = Float32[10.0, 20.0, 10.5, 20.5],
                cv_fold       = UInt8[0, 0, 0, 0],
            )
            for c in _P.GROUP_COLS
                df[!, c] = Float32.(1:n)   # distinct, identical per column
            end
            df[!, :global_max] = Float32[1.0, 2.0, 3.0, 4.0]
            s1 = Float32[0.8, 0.7, 0.6, 0.5]
            path = _write_fixture!(dir, "run1_fold0.arrow", df, s1)

            orig_gmax = copy(df.global_max)

            injected = _P.inject_shadow_decoys!([path])
            @test injected == 4

            out = DataFrame(Arrow.Table(path))
            # Row count doubles (both classes present -> a shadow per real row).
            @test nrow(out) == 8
            @test hasproperty(out, :is_shadow)
            real_rows   = out[.!out.is_shadow, :]
            shadow_rows = out[out.is_shadow, :]
            @test nrow(real_rows) == 4
            @test nrow(shadow_rows) == 4
            # Originals keep is_shadow == false.
            @test all(.!Bool.(real_rows.is_shadow))

            # 1:1 target/decoy marginal after injection (2T+2D real, opposite
            # shadows add 2D+2T).
            @test count(out.target) == 4
            @test count(.!out.target) == 4
            # Shadows carry the OPPOSITE class: 2 targets + 2 decoys among shadows.
            @test count(shadow_rows.target) == 2
            @test count(.!shadow_rows.target) == 2

            # The grafted GROUP feature (global_max) on shadows equals the set of
            # ORIGINAL rows' values (each real row grafts exactly one shadow).
            @test sort(shadow_rows.global_max) == sort(orig_gmax)

            # Targeted check: real target row (irt 10) -> nearest decoy is
            # precursor 3 (irt 10.5); its shadow is a decoy grafted with row1's
            # global_max (1.0).
            @test any(r -> !r.target && r.global_max == 1.0f0, eachrow(shadow_rows))

            # Sidecar stays aligned (row count matches the main file).
            side_after = DataFrame(Arrow.Table(path * _P.PASS1_SIDECAR_SUFFIX))
            @test nrow(side_after) == nrow(out)

            # remove restores the original layout.
            removed = _P.remove_shadow_decoys!([path])
            @test removed == 4
            restored = DataFrame(Arrow.Table(path))
            @test nrow(restored) == 4
            @test !hasproperty(restored, :is_shadow)   # flag column dropped
            @test sort(restored.precursor_idx) == UInt32[1, 2, 3, 4]
            side_restored = DataFrame(Arrow.Table(path * _P.PASS1_SIDECAR_SUFFIX))
            @test nrow(side_restored) == 4
        end
    end

    # ---------------------------------------------------------------
    # 7. _write_group_features! leave-one-out (cluster OFF, R<9)
    # ---------------------------------------------------------------
    @testset "_write_group_features! LOO" begin
        mktempdir() do dir
            # Fold 0, three runs (files), a single precursor (id=1) present in
            # all three with distinct s1 and irt. R=3 < 9 -> cluster OFF, so
            # cluster cols mirror global cols. Everything is hand-checkable.
            files = String[]
            specs = [(0.80f0, 100.0f0), (0.90f0, 110.0f0), (0.99f0, 120.0f0)]
            for (ri, (s1v, irtv)) in enumerate(specs)
                df = DataFrame(
                    precursor_idx = UInt32[1],
                    cv_fold       = UInt8[0],
                    irt_obs       = Float32[irtv],
                )
                push!(files, _write_fixture!(dir, "run$(ri)_fold0.arrow", df, Float32[s1v]))
            end

            _P._write_group_features!(files)

            # Read back per-file single-row values.
            gm = Float32[]; gt3 = Float32[]; np = Float32[]; nc = Float32[]
            cm = Float32[]; ct3 = Float32[]; cnp = Float32[]; cnc = Float32[]
            di = Float32[]
            for p in files
                t = DataFrame(Arrow.Table(p))
                @test nrow(t) == 1
                for c in _P.GROUP_COLS
                    @test hasproperty(t, c)
                end
                @test hasproperty(t, :delta_irt)
                push!(gm, t.global_max[1]);            push!(gt3, t.global_top3_logodds[1])
                push!(np, t.global_n_present[1]);      push!(nc, t.global_n_confident[1])
                push!(cm, t.cluster_max[1]);           push!(ct3, t.cluster_top3_logodds[1])
                push!(cnp, t.cluster_n_present[1]);    push!(cnc, t.cluster_n_confident[1])
                push!(di, t.delta_irt[1])
            end
            @test all(isfinite, gm) && all(isfinite, gt3) && all(isfinite, di)

            # global_max = max s1 of the OTHER runs (leave-one-out).
            @test gm[1] ≈ 0.99f0   # run1(0.80): max(0.90,0.99)
            @test gm[2] ≈ 0.99f0   # run2(0.90): max(0.80,0.99)
            @test gm[3] ≈ 0.90f0   # run3(0.99): max(0.80,0.90)

            # n_present LOO = (runs present) - 1 = 3 - 1 = 2 everywhere.
            @test np == Float32[2, 2, 2]

            # n_confident LOO: GROUP_CONF_S1 = 0.99. Only run3 (0.99) is confident.
            # run1/run2 exclude own (not confident) -> 1; run3 excludes own -> 0.
            @test nc == Float32[1, 1, 0]

            # top3_logodds LOO = sum of positive logits of other runs' scores.
            @test gt3[1] ≈ _P._pos_logit(0.99f0) + _P._pos_logit(0.90f0) atol=1e-5
            @test gt3[2] ≈ _P._pos_logit(0.99f0) + _P._pos_logit(0.80f0) atol=1e-5
            @test gt3[3] ≈ _P._pos_logit(0.90f0) + _P._pos_logit(0.80f0) atol=1e-5

            # delta_irt = |own irt - ref irt|, ref = irt of highest-s1 run (120).
            @test di[1] ≈ 20.0f0
            @test di[2] ≈ 10.0f0
            @test di[3] ≈ 0.0f0

            # Cluster OFF (R<9): cluster cols mirror global cols exactly.
            @test cm == gm
            @test ct3 == gt3
            @test cnp == np
            @test cnc == nc
        end
    end

end
