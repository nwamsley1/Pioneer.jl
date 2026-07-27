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
    # 3b. GTopK4: donor evidence stays aligned with the top-2 scores
    # ---------------------------------------------------------------
    @testset "GTopK4 donor evidence" begin
        @test _P.EMPTY_GTOPK4.t === _P.EMPTY_TOPK4
        @test _P.EMPTY_GTOPK4.e == (_P.EMPTY_EV, _P.EMPTY_EV)

        # Evidence is tagged with the run number so we can assert which run's
        # values ended up in each slot. Insert out of score order on purpose.
        ev(run) = _P.DonorEv(Float32(10 * run), Float32(run),
                             ntuple(i -> Float16(i == run ? 1 : 0), 8))
        g = _P.EMPTY_GTOPK4
        g = _P._gt_insert(g, 0.60f0, UInt32(1), ev(1))
        g = _P._gt_insert(g, 0.99f0, UInt32(2), ev(2))   # new top; 1 shifts to slot 2
        g = _P._gt_insert(g, 0.90f0, UInt32(3), ev(3))   # lands in slot 2; 1 falls out of evidence
        g = _P._gt_insert(g, 0.70f0, UInt32(4), ev(4))
        g = _P._gt_insert(g, 0.80f0, UInt32(5), ev(5))

        # Scores/runs behave exactly as the wrapped TopK4.
        @test g.t.v == (0.99f0, 0.90f0, 0.80f0, 0.70f0)
        @test g.t.r == (UInt32(2), UInt32(3), UInt32(5), UInt32(4))
        # Evidence slots track the top-2 runs (2 then 3).
        @test g.e[1] === ev(2)
        @test g.e[2] === ev(3)

        # A score below the current 4th never touches evidence.
        g2 = _P._gt_insert(g, 0.10f0, UInt32(9), ev(9))
        @test g2 === g

        # LOO donor: excluding a non-top run keeps slot 1 (run 2).
        s, d = _P._gt_loo_donor(g, UInt32(5))
        @test s ≈ 0.99f0
        @test d === ev(2)
        # Excluding the top run itself falls through to slot 2 (run 3).
        s, d = _P._gt_loo_donor(g, UInt32(2))
        @test s ≈ 0.90f0
        @test d === ev(3)
        # No other run at all -> score 0 signals "no donor".
        s, d = _P._gt_loo_donor(_P.EMPTY_GTOPK4, UInt32(1))
        @test s == 0.0f0
        @test d === _P.EMPTY_EV

        # Insertion order must not matter: descending order gives the same record.
        h = _P.EMPTY_GTOPK4
        for (sc, ru) in ((0.99f0, 2), (0.90f0, 3), (0.80f0, 5), (0.70f0, 4), (0.60f0, 1))
            h = _P._gt_insert(h, sc, UInt32(ru), ev(ru))
        end
        @test h.t.v == g.t.v && h.t.r == g.t.r && h.e == g.e
    end

    # ---------------------------------------------------------------
    # 3c. Fragment-spectrum normalization + donor Hellinger
    # ---------------------------------------------------------------
    @testset "donor frag Hellinger" begin
        f(v...) = _P._frag_sqrt_tuple(NTuple{8,Float32}(Float32.(v)))

        # L1-normalize then sqrt: a single non-zero rank -> unit mass on that rank,
        # regardless of its raw magnitude.
        @test f(4, 0, 0, 0, 0, 0, 0, 0) == f(1, 0, 0, 0, 0, 0, 0, 0)
        @test Float32(f(4, 0, 0, 0, 0, 0, 0, 0)[1]) ≈ 1.0f0
        # Two equal ranks -> sqrt(1/2) each.
        two = f(1, 1, 0, 0, 0, 0, 0, 0)
        @test Float32(two[1]) ≈ sqrt(0.5f0) atol=1e-3
        @test Float32(two[2]) ≈ sqrt(0.5f0) atol=1e-3
        # No intensity -> the empty tuple, which is NOT a valid spectrum.
        @test f(0, 0, 0, 0, 0, 0, 0, 0) == _P.EMPTY_FRAG
        @test !_P._frag_is_valid(_P.EMPTY_FRAG)
        @test _P._frag_is_valid(two)

        a = f(1, 0, 0, 0, 0, 0, 0, 0)
        b = f(0, 1, 0, 0, 0, 0, 0, 0)
        # Identical spectra -> 0; fully disjoint -> 1 (the Hellinger bounds).
        @test _P._donor_hellinger(a, a) ≈ 0.0f0 atol=1e-3
        @test _P._donor_hellinger(a, b) ≈ 1.0f0 atol=1e-3
        # Partial overlap sits strictly between, and is symmetric.
        h_ab = _P._donor_hellinger(two, b)
        @test 0.0f0 < h_ab < 1.0f0
        @test h_ab ≈ _P._donor_hellinger(b, two) atol=1e-4
        @test h_ab ≈ sqrt(0.5f0 * ((sqrt(0.5f0))^2 + (sqrt(0.5f0) - 1f0)^2)) atol=1e-3
        # Either side empty -> sentinel, never 0 (an absent spectrum must not read
        # as perfect agreement).
        @test _P._donor_hellinger(a, _P.EMPTY_FRAG) == _P.DONOR_SENTINEL
        @test _P._donor_hellinger(_P.EMPTY_FRAG, a) == _P.DONOR_SENTINEL
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
        # `SHADOW_SIDECAR_SCORE_COL` was removed in 2093a9477; the Pass-1 sidecar's
        # round-1 OOF score column is `trace_prob_prepass` (see pass1_oom.jl).
        side = DataFrame(trace_prob_prepass = s1)
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
            # (s1, irt_obs, weight, log2_explained, frag-spectrum) — the last three feed DONOR_COLS.
            # Spectra are chosen so the Hellingers are hand-checkable: run1 and run3 are disjoint
            # (-> 1.0), run2 half-overlaps each.
            specs = [(0.80f0, 100.0f0, 1.0f0, 2.0f0, Float32[1, 0, 0, 0, 0, 0, 0, 0]),
                     (0.90f0, 110.0f0, 2.0f0, 3.0f0, Float32[1, 1, 0, 0, 0, 0, 0, 0]),
                     (0.99f0, 120.0f0, 4.0f0, 5.0f0, Float32[0, 1, 0, 0, 0, 0, 0, 0])]
            for (ri, (s1v, irtv, wv, ev, fr)) in enumerate(specs)
                df = DataFrame(
                    precursor_idx = UInt32[1],
                    cv_fold       = UInt8[0],
                    irt_obs       = Float32[irtv],
                    weight        = Float32[wv],
                    log2_intensity_explained = Float32[ev],
                )
                for (r, c) in enumerate(_P.SMOOTHED_FRAGMENT_INTENSITY_COLUMNS)
                    df[!, c] = Float32[fr[r]]
                end
                push!(files, _write_fixture!(dir, "run$(ri)_fold0.arrow", df, Float32[s1v]))
            end

            _P._write_group_features!(files, 1)   # max_pid = 1 (single precursor id in this fixture)

            # Read back per-file single-row values. The columns are written to a
            # row-aligned `.group.sidecar.arrow`, not the main fold file.
            gm = Float32[]; gt3 = Float32[]; np = Float32[]; nc = Float32[]
            cm = Float32[]; ct3 = Float32[]; cnp = Float32[]; cnc = Float32[]
            di = Float32[]
            dh = Float32[]; dlw = Float32[]; dle = Float32[]
            for p in files
                @test nrow(DataFrame(Arrow.Table(p))) == 1     # main file untouched
                t = DataFrame(Arrow.Table(p * _P.GROUP_SIDECAR_SUFFIX))
                @test nrow(t) == 1
                for c in vcat(_P.GROUP_COLS, _P.DONOR_COLS)
                    @test hasproperty(t, c)
                end
                @test hasproperty(t, :delta_irt)
                push!(gm, t.global_max[1]);            push!(gt3, t.global_top3_logodds[1])
                push!(np, t.global_n_present[1]);      push!(nc, t.global_n_confident[1])
                push!(cm, t.cluster_max[1]);           push!(ct3, t.cluster_top3_logodds[1])
                push!(cnp, t.cluster_n_present[1]);    push!(cnc, t.cluster_n_confident[1])
                push!(di, t.delta_irt[1])
                push!(dh, t.donor_frag_hellinger[1])
                push!(dlw, t.donor_log2_weight_ratio[1])
                push!(dle, t.donor_log2_explained_ratio[1])
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

            # DONOR_COLS compare each row against its LEAVE-ONE-OUT top-scoring run:
            # run1 and run2 both donate from run3 (0.99); run3 donates from run2 (0.90).
            @test dlw[1] ≈ log2(1.0f0 / 4.0f0)   # -2
            @test dlw[2] ≈ log2(2.0f0 / 4.0f0)   # -1
            @test dlw[3] ≈ log2(4.0f0 / 2.0f0)   # +1

            @test dle[1] ≈ 2.0f0 - 5.0f0         # -3
            @test dle[2] ≈ 3.0f0 - 5.0f0         # -2
            @test dle[3] ≈ 5.0f0 - 3.0f0         # +2

            # run1's spectrum is disjoint from its donor run3's -> max distance.
            @test dh[1] ≈ 1.0f0 atol=1e-3
            # run2 half-overlaps run3, and run3 half-overlaps run2 -> equal, in (0,1).
            half = sqrt(0.5f0 * ((sqrt(0.5f0))^2 + (sqrt(0.5f0) - 1f0)^2))
            @test dh[2] ≈ half atol=1e-3
            @test dh[3] ≈ half atol=1e-3
            @test 0.0f0 < dh[2] < 1.0f0
        end
    end

    # ---------------------------------------------------------------
    # 8. DONOR_COLS sentinels: no other run, and a non-positive weight
    # ---------------------------------------------------------------
    @testset "DONOR_COLS sentinels" begin
        mktempdir() do dir
            # Fold 0, two runs. Precursor 1 is in both; precursor 2 only in run1
            # (no donor). Run1's precursor 1 has weight 0 -> ratio undefined.
            df1 = DataFrame(
                precursor_idx = UInt32[1, 2],
                cv_fold       = UInt8[0, 0],
                irt_obs       = Float32[100.0, 50.0],
                weight        = Float32[0.0, 5.0],
                log2_intensity_explained = Float32[4.0, 4.0],
            )
            df2 = DataFrame(
                precursor_idx = UInt32[1],
                cv_fold       = UInt8[0],
                irt_obs       = Float32[110.0],
                weight        = Float32[2.0],
                log2_intensity_explained = Float32[1.0],
            )
            # Row 1 of df1 carries a valid spectrum (so its Hellinger is real even though
            # its weight is 0); row 2 has none, which is a second, independent sentinel path.
            for (r, c) in enumerate(_P.SMOOTHED_FRAGMENT_INTENSITY_COLUMNS)
                df1[!, c] = Float32[r == 1 ? 1.0 : 0.0, 0.0]
                df2[!, c] = Float32[r == 1 ? 1.0 : 0.0]
            end
            f1 = _write_fixture!(dir, "run1_fold0.arrow", df1, Float32[0.80, 0.95])
            f2 = _write_fixture!(dir, "run2_fold0.arrow", df2, Float32[0.90])

            _P._write_group_features!([f1, f2], 2)

            s1 = DataFrame(Arrow.Table(f1 * _P.GROUP_SIDECAR_SUFFIX))
            # Row 1 (precursor 1): a donor exists, but own weight is 0 -> ONLY the weight
            # ratio takes the sentinel. The spectra are identical (both unit mass on rank 1)
            # so the Hellinger is 0, and the explained ratio is a real value.
            @test s1.donor_log2_weight_ratio[1] == _P.DONOR_SENTINEL
            @test s1.donor_frag_hellinger[1] ≈ 0.0f0 atol=1e-3
            @test s1.donor_log2_explained_ratio[1] ≈ 4.0f0 - 1.0f0
            # Row 2 (precursor 2): present in no other run -> all three sentinel, and
            # global_n_present == 0 is what disambiguates them.
            @test s1.global_n_present[2] == 0.0f0
            @test s1.donor_log2_weight_ratio[2] == _P.DONOR_SENTINEL
            @test s1.donor_frag_hellinger[2] == _P.DONOR_SENTINEL
            @test s1.donor_log2_explained_ratio[2] == _P.DONOR_SENTINEL
        end
    end

end
