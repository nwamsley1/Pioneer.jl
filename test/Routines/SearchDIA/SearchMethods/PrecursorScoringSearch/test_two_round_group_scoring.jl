# Unit tests for the two-round precursor cross-run scoring machinery in
# src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/two_round_scoring.jl
#
# All functions under test are module-internal (not exported) and are referenced
# as `Pioneer._name`. Fixtures are tiny synthetic fold Arrow files + matching
# `.pass1_sidecar.arrow` sidecars, exercising the real group-feature code path.
# Randomness in the clustering path uses fixed internal seeds, so the tests
# assert on partition structure and remain deterministic.

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
    # Minimal fold-file plus the Pass-1 sidecar holding its round-1 OOF scores
    # (`trace_prob_prepass`) — what _write_group_features! reads as s1.
    function _write_fixture!(dir, name, df::DataFrame, s1::Vector{Float32})
        path = joinpath(dir, name)
        _P.writeArrow(path, df)
        side = DataFrame(trace_prob_prepass = s1)
        _P.writeArrow(path * _P.PASS1_SIDECAR_SUFFIX, side)
        return path
    end

    # ---------------------------------------------------------------
    # 6. _shadow_partners — nearest-iRT opposite-class pairing
    # ---------------------------------------------------------------
    # Shadow decoys are no longer round-tripped through the fold files: as of
    # 2093a9477 they are built in memory inside the training sample
    # (`_sample_both_folds(...; inject_shadows=true)`), so
    # inject_shadow_decoys! / remove_shadow_decoys! and their `is_shadow`
    # column no longer exist. What survives — and what the old round-trip test
    # was really pinning down — is the pairing rule that decides which
    # opposite-class row each shadow is drawn from.
    @testset "_shadow_partners" begin
        # 2 targets (irt 10, 20), 2 decoys (irt 10.5, 20.5): each row's nearest
        # opposite-class neighbour is the one it shares an iRT with.
        targets = Bool[true, true, false, false]
        irt     = Float32[10.0, 20.0, 10.5, 20.5]
        partner = _P._shadow_partners(targets, irt)
        @test partner == Int32[3, 4, 1, 2]
        # Every partner is genuinely the opposite class.
        @test all(targets[i] != targets[partner[i]] for i in eachindex(partner))

        # Exact midpoint ties break to the LOWER index, inherited from
        # _nearest_sorted_idx (see testset 5): target irt 10 sits halfway
        # between decoys at 5 and 15 -> picks the one at 5 (row 2).
        @test _P._shadow_partners(Bool[true, false, false],
                                  Float32[10.0, 5.0, 15.0]) == Int32[2, 1, 1]

        # Single-class input has no possible partner -> all zeros. The caller
        # skips such files rather than injecting (`pairable`).
        @test _P._shadow_partners(Bool[true, true], Float32[1.0, 2.0]) == Int32[0, 0]
        @test _P._shadow_partners(Bool[false, false], Float32[1.0, 2.0]) == Int32[0, 0]

        # Many-to-one is allowed: one decoy can be the nearest for several
        # targets. The lone decoy at 2.5 serves all three targets; it in turn
        # ties between targets 2 and 3 and takes the lower index.
        @test _P._shadow_partners(Bool[true, true, true, false],
                                  Float32[1.0, 2.0, 3.0, 2.5]) == Int32[4, 4, 4, 2]
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

            _P._write_group_features!(files, 1)   # max_pid = 1 (single precursor id in this fixture)

            # Read back per-file single-row values.
            gm = Float32[]; gt3 = Float32[]; np = Float32[]; nc = Float32[]
            cm = Float32[]; ct3 = Float32[]; cnp = Float32[]; cnc = Float32[]
            di = Float32[]
            for p in files
                t = DataFrame(Arrow.Table(p))
                @test nrow(t) == 1
                # The 9 GROUP columns go to a row-aligned `.group.sidecar.arrow`
                # rather than being written back into the ~80-column fold file.
                # Round-2 read sites resolve them via `_feature_columns`, and
                # Step 1b's PSMFileReference auto-discovers the sidecar.
                g = DataFrame(Arrow.Table(p * _P.GROUP_SIDECAR_SUFFIX))
                @test nrow(g) == nrow(t)
                for c in _P.GROUP_COLS
                    @test hasproperty(g, c)
                end
                @test hasproperty(g, :delta_irt)
                push!(gm, g.global_max[1]);            push!(gt3, g.global_top3_logodds[1])
                push!(np, g.global_n_present[1]);      push!(nc, g.global_n_confident[1])
                push!(cm, g.cluster_max[1]);           push!(ct3, g.cluster_top3_logodds[1])
                push!(cnp, g.cluster_n_present[1]);    push!(cnc, g.cluster_n_confident[1])
                push!(di, g.delta_irt[1])
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
