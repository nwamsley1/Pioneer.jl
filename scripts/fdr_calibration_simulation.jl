#!/usr/bin/env julia
#=
FDR Calibration Simulation — CV vs No-CV
=========================================
Tests whether cross-validation is necessary for calibrated FDR.

Design:
  - Distribution A: MVN(+μ, I) — all labeled "target"
  - Distribution B: MVN(-μ, I) — half "target" (false targets), half "decoy"
  - With CV: train on fold 1, predict fold 0; train on fold 0, predict fold 1
  - Without CV: train on ALL data, predict on ALL data (overfitting possible)

If inflation > 1.0: false targets pass more than decoys → FDR is anti-conservative.
=#

using Random, Statistics, Printf
using LightGBM

const NT = Threads.nthreads()

# ─── Q-values ───────────────────────────────────────────────────────────────

function get_qvalues!(probs::Vector{Float32}, labels::Vector{Bool})::Vector{Float64}
    n = length(probs)
    qvals = Vector{Float64}(undef, n)
    order = sortperm(probs; rev=true)
    t = 0; d = 0
    @inbounds for i in 1:n
        idx = order[i]
        if labels[idx]; t += 1 else d += 1 end
        qvals[idx] = t > 0 ? Float64(d) / Float64(t) : 1.0
    end
    min_q = qvals[order[n]]
    @inbounds for i in n:-1:1
        idx = order[i]; min_q = min(min_q, qvals[idx]); qvals[idx] = min_q
    end
    return qvals
end

# ─── LightGBM scoring ──────────────────────────────────────────────────────

function make_estimator(; iters::Int=50, lr::Float64=0.2, depth::Int=3, leaves::Int=10, minleaf::Int=5)
    LGBMClassification(;
        num_iterations=iters, learning_rate=lr, max_depth=depth, num_leaves=leaves,
        min_data_in_leaf=minleaf, feature_fraction=0.5, bagging_fraction=0.5,
        bagging_freq=1, objective="binary", metric=["binary_logloss"], num_class=1,
        is_unbalance=false, verbosity=-1, num_threads=NT,
        deterministic=true, force_row_wise=true, seed=1776)
end

function lgbm_score(X::Matrix{Float32}, targets::Vector{Bool}, cv::Vector{UInt8}, use_cv::Bool;
                    iters::Int=50, lr::Float64=0.2, depth::Int=3, leaves::Int=10, minleaf::Int=5)
    n = size(X, 1)
    scores = Vector{Float32}(undef, n)

    if use_cv
        for (tr, te) in ((UInt8(1), UInt8(0)), (UInt8(0), UInt8(1)))
            tidx = findall(cv .== tr)
            eidx = findall(cv .== te)
            y = Int.(targets[tidx])
            np = count(==(1), y)
            if np < 10 || (length(y) - np) < 10
                scores[eidx] .= Float32(np / length(y)); continue
            end
            est = make_estimator(; iters, lr, depth, leaves, minleaf)
            LightGBM.fit!(est, X[tidx, :], y; verbosity=-1)
            raw = LightGBM.predict(est, X[eidx, :])
            scores[eidx] .= Float32.(ndims(raw) == 2 ? vec(raw) : raw)
        end
    else
        # No CV: train and predict on same data
        y = Int.(targets)
        est = make_estimator(; iters, lr, depth, leaves, minleaf)
        LightGBM.fit!(est, X, y; verbosity=-1)
        raw = LightGBM.predict(est, X)
        scores .= Float32.(ndims(raw) == 2 ? vec(raw) : raw)
    end
    return scores
end

# ─── Data generation ────────────────────────────────────────────────────────

function gen_data(rng::AbstractRNG, n_true::Int, n_false::Int, nf::Int, sep::Float64)
    n = n_true + n_false
    nft = n_false ÷ 2
    X = randn(rng, Float32, n, nf)
    @inbounds for i in 1:n_true; X[i, 1] += Float32(sep / 2); end
    @inbounds for i in n_true+1:n; X[i, 1] -= Float32(sep / 2); end

    targets = Vector{Bool}(undef, n)
    is_A = Vector{Bool}(undef, n)
    cv = Vector{UInt8}(undef, n)
    @inbounds for i in 1:n
        is_A[i] = i <= n_true
        targets[i] = i <= n_true + nft
        cv[i] = UInt8(i % 2)
    end

    perm = randperm(rng, n)
    return X[perm, :], targets[perm], is_A[perm], cv[perm]
end

# ─── Analysis ───────────────────────────────────────────────────────────────

function analyze(passing::BitVector, is_A::Vector{Bool}, targets::Vector{Bool})
    nft = count(passing .& .!is_A .& targets)
    nd = count(passing .& .!is_A .& .!targets)
    tot_ft = count(.!is_A .& targets)
    tot_d = count(.!is_A .& .!targets)
    ft_r = nft / max(1, tot_ft)
    d_r = nd / max(1, tot_d)
    n_tgt_pass = count(passing .& targets)
    true_fdr = n_tgt_pass > 0 ? Float64(nft) / Float64(n_tgt_pass) : 0.0
    tda_fdr = n_tgt_pass > 0 ? Float64(nd) / Float64(n_tgt_pass) : 0.0
    return (pass=count(passing), A=count(passing .& is_A), ft=nft, dec=nd,
            inflation=ft_r / max(1e-10, d_r), true_fdr=true_fdr, tda_fdr=tda_fdr)
end

function pr(label::String, r)
    @printf("  %-42s  pass=%6d  A=%6d  falseT=%4d  decoy=%4d  inflation=%6.2fx  true_FDR=%.4f  tda_FDR=%.4f\n", label, r.pass, r.A, r.ft, r.dec, r.inflation, r.true_fdr, r.tda_fdr)
end

# ─── Experiment runner ──────────────────────────────────────────────────────

function run_exp(; nt::Int=50000, nf_total::Int=100000, nfeat::Int=15, sep::Float64=2.0,
                  use_cv::Bool=true, q::Float64=0.01, seed::Int=42,
                  iters::Int=50, depth::Int=3, leaves::Int=10, minleaf::Int=5)
    X, tgt, isA, cv = gen_data(MersenneTwister(seed), nt, nf_total, nfeat, sep)
    sc = lgbm_score(X, tgt, cv, use_cv; iters, depth, leaves, minleaf)
    qv = get_qvalues!(sc, tgt)
    return analyze(BitVector(qv .<= q), isA, tgt)
end

# ─── PEP calibration (PAVA isotonic regression) ──────────────────────────────

function pep_calibrate(scores::Vector{Float32}, labels::Vector{Bool})::Vector{Float64}
    n = length(scores)
    order = sortperm(scores; rev=true)
    raw = Vector{Float64}(undef, n)
    @inbounds for i in 1:n; raw[i] = labels[order[i]] ? 0.0 : 1.0; end

    bval = Vector{Float64}(undef, n)
    bwt = Vector{Int}(undef, n)
    nb = 0
    @inbounds for i in 1:n
        nb += 1; bval[nb] = raw[i]; bwt[nb] = 1
        while nb > 1 && bval[nb] < bval[nb-1]
            w = bwt[nb-1] + bwt[nb]
            bval[nb-1] = (bval[nb-1] * bwt[nb-1] + bval[nb] * bwt[nb]) / w
            bwt[nb-1] = w; nb -= 1
        end
    end

    pep = Vector{Float64}(undef, n)
    bi = 1; pos = 0
    @inbounds for i in 1:n
        pos += 1
        if pos > bwt[bi]; bi += 1; pos = 1; end
        pep[order[i]] = clamp(bval[bi], 1e-10, 1.0 - 1e-10)
    end
    return pep
end

# ─── Three-way single-file comparison ────────────────────────────────────────

function run_threeway(; nt::Int=50000, nf_total::Int=100000, nfeat::Int=15, sep::Float64=2.5,
                      q::Float64=0.01, seed::Int=42,
                      iters::Int=50, depth::Int=3, leaves::Int=10, minleaf::Int=5)
    X, tgt, isA, cv = gen_data(MersenneTwister(seed), nt, nf_total, nfeat, sep)

    sc_nocv = lgbm_score(X, tgt, cv, false; iters, depth, leaves, minleaf)
    qv = get_qvalues!(sc_nocv, tgt)
    r_nocv = analyze(BitVector(qv .<= q), isA, tgt)

    sc_cv = lgbm_score(X, tgt, cv, true; iters, depth, leaves, minleaf)

    qv = get_qvalues!(sc_cv, tgt)
    r_pooled = analyze(BitVector(qv .<= q), isA, tgt)

    idx0 = findall(cv .== UInt8(0)); idx1 = findall(cv .== UInt8(1))
    qv0 = get_qvalues!(sc_cv[idx0], tgt[idx0])
    qv1 = get_qvalues!(sc_cv[idx1], tgt[idx1])
    passing = falses(length(tgt))
    passing[idx0] .= qv0 .<= q; passing[idx1] .= qv1 .<= q
    r_perfold = analyze(passing, isA, tgt)

    return r_nocv, r_pooled, r_perfold
end

# ─── Multi-file experiment with PEP calibration ─────────────────────────────

function run_multifile(; n_true::Int=5000, n_false::Int=10000, n_files::Int=5,
                        nfeat::Int=15, sep::Float64=2.5, q::Float64=0.01, seed::Int=42,
                        iters::Int=50, depth::Int=3, leaves::Int=10, minleaf::Int=5,
                        cv_mode::Symbol=:cv_perfold)
    rng = MersenneTwister(seed)
    n = n_true + n_false
    nft = n_false ÷ 2

    is_A_raw = Vector{Bool}(vcat(trues(n_true), falses(n_false)))
    tgt_raw = Vector{Bool}(vcat(trues(n_true), trues(nft), falses(nft)))
    cv_raw = UInt8[UInt8(i % 2) for i in 1:n]
    perm = randperm(rng, n)
    is_A = is_A_raw[perm]; tgt = tgt_raw[perm]; cv = cv_raw[perm]

    logodds = zeros(Float64, n)

    for f in 1:n_files
        X = randn(rng, Float32, n, nfeat)
        @inbounds for i in 1:n
            X[i, 1] += is_A[i] ? Float32(sep / 2) : Float32(-sep / 2)
        end

        use_cv = (cv_mode != :nocv)
        sc = lgbm_score(X, tgt, cv, use_cv; iters, depth, leaves, minleaf)

        if cv_mode == :cv_perfold
            for fold in UInt8[0, 1]
                idx = findall(cv .== fold)
                pep = pep_calibrate(sc[idx], tgt[idx])
                logodds[idx] .+= log.((1.0 .- pep) ./ pep)
            end
        else
            pep = pep_calibrate(sc, tgt)
            logodds .+= log.((1.0 .- pep) ./ pep)
        end
    end

    global_prob = Float32.(1.0 ./ (1.0 .+ exp.(-logodds ./ n_files)))

    if cv_mode == :cv_perfold
        passing = falses(n)
        for fold in UInt8[0, 1]
            idx = findall(cv .== fold)
            qv = get_qvalues!(global_prob[idx], tgt[idx])
            passing[idx] .= qv .<= q
        end
    else
        qv = get_qvalues!(global_prob, tgt)
        passing = BitVector(qv .<= q)
    end

    return analyze(passing, is_A, tgt)
end

# ─── Main ───────────────────────────────────────────────────────────────────

function main()
    print("Warming up JIT... ")
    run_exp(; nt=200, nf_total=400, nfeat=3, sep=3.0, seed=99)
    run_exp(; nt=200, nf_total=400, nfeat=3, sep=3.0, seed=99, use_cv=false)
    run_threeway(; nt=200, nf_total=400, nfeat=3, sep=3.0, seed=99)
    run_multifile(; n_true=100, n_false=200, n_files=2, nfeat=3, sep=3.0, seed=99, cv_mode=:cv_perfold)
    run_multifile(; n_true=100, n_false=200, n_files=2, nfeat=3, sep=3.0, seed=99, cv_mode=:cv_pooled)
    run_multifile(; n_true=100, n_false=200, n_files=2, nfeat=3, sep=3.0, seed=99, cv_mode=:nocv)
    println("done.\n")

    N = 50000; NF = 100000

    println("=" ^ 120)
    println("FDR CALIBRATION: CROSS-VALIDATION vs NO CROSS-VALIDATION")
    println("=" ^ 120)
    println("A=target, B=half target/half decoy. inflation=1.0 → calibrated. >1.0 → anti-conservative.\n")

    # ── Vary separation (difficulty) ──
    println("─" ^ 120)
    println("EXP 1: Vary separation (n_true=50k, n_false=100k, 15 features, 50 trees depth=3)")
    println("─" ^ 120)
    for sep in [1.5, 2.0, 2.5, 3.0, 4.0]
        rc = run_exp(; nt=N, nf_total=NF, sep, use_cv=true)
        rn = run_exp(; nt=N, nf_total=NF, sep, use_cv=false)
        pr(@sprintf("sep=%.1f  WITH CV", sep), rc)
        pr(@sprintf("sep=%.1f  NO CV", sep), rn)
        println()
    end

    # ── Vary number of features ──
    println("─" ^ 120)
    println("EXP 2: Vary features (separation=2.5, 50 trees depth=3)")
    println("─" ^ 120)
    for nfeat in [1, 3, 5, 10, 15, 30, 50]
        rc = run_exp(; nt=N, nf_total=NF, nfeat, sep=2.5, use_cv=true)
        rn = run_exp(; nt=N, nf_total=NF, nfeat, sep=2.5, use_cv=false)
        pr(@sprintf("features=%-3d  WITH CV", nfeat), rc)
        pr(@sprintf("features=%-3d  NO CV", nfeat), rn)
        println()
    end

    # ── Vary model complexity (deeper trees = more overfitting risk) ──
    println("─" ^ 120)
    println("EXP 3: Vary model complexity (separation=2.5, 15 features)")
    println("─" ^ 120)
    for (d, l, label) in [(2, 4, "depth=2/leaves=4"), (3, 10, "depth=3/leaves=10"),
                           (5, 31, "depth=5/leaves=31"), (10, 63, "depth=10/leaves=63")]
        rc = run_exp(; nt=N, nf_total=NF, sep=2.5, depth=d, leaves=l, use_cv=true)
        rn = run_exp(; nt=N, nf_total=NF, sep=2.5, depth=d, leaves=l, use_cv=false)
        pr(@sprintf("%-20s  WITH CV", label), rc)
        pr(@sprintf("%-20s  NO CV", label), rn)
        println()
    end

    # ── Vary iterations (more rounds = more overfitting risk) ──
    println("─" ^ 120)
    println("EXP 4: Vary boosting rounds (separation=2.5, 15 features, depth=3)")
    println("─" ^ 120)
    for it in [10, 50, 200, 500]
        rc = run_exp(; nt=N, nf_total=NF, sep=2.5, iters=it, use_cv=true)
        rn = run_exp(; nt=N, nf_total=NF, sep=2.5, iters=it, use_cv=false)
        pr(@sprintf("iters=%-4d  WITH CV", it), rc)
        pr(@sprintf("iters=%-4d  NO CV", it), rn)
        println()
    end

    # ── Vary sample size (small = more overfitting risk) ──
    println("─" ^ 120)
    println("EXP 5: Vary sample size (separation=2.5, 15 features, depth=3)")
    println("─" ^ 120)
    for (nt, nf) in [(1000, 2000), (5000, 10000), (10000, 20000), (50000, 100000)]
        rc = run_exp(; nt, nf_total=nf, sep=2.5, use_cv=true)
        rn = run_exp(; nt, nf_total=nf, sep=2.5, use_cv=false)
        pr(@sprintf("n=%dk  WITH CV", (nt+nf)÷1000), rc)
        pr(@sprintf("n=%dk  NO CV", (nt+nf)÷1000), rn)
        println()
    end

    # ── Stress test: many features + deep trees + small sample ──
    println("─" ^ 120)
    println("EXP 6: Stress test — conditions that maximize overfitting risk")
    println("─" ^ 120)
    for (label, nt, nf, nfeat, d, l, it) in [
        ("baseline",       50000, 100000, 15, 3,  10,  50),
        ("50 features",    50000, 100000, 50, 3,  10,  50),
        ("deep trees",     50000, 100000, 15, 10, 63,  50),
        ("many rounds",    50000, 100000, 15, 3,  10,  500),
        ("small sample",   2000,  4000,   15, 3,  10,  50),
        ("small+deep",     2000,  4000,   15, 10, 63,  50),
        ("small+50feat",   2000,  4000,   50, 3,  10,  50),
        ("worst case",     2000,  4000,   50, 10, 63,  500),
    ]
        rc = run_exp(; nt, nf_total=nf, nfeat, sep=2.5, depth=d, leaves=l, iters=it, use_cv=true)
        rn = run_exp(; nt, nf_total=nf, nfeat, sep=2.5, depth=d, leaves=l, iters=it, use_cv=false)
        pr(@sprintf("%-18s WITH CV", label), rc)
        pr(@sprintf("%-18s NO CV", label), rn)
        println()
    end

    # ── Single-file: pooled vs per-fold q-values ──
    println("─" ^ 120)
    println("EXP 7: Single file — CV+pooled q-values vs CV+per-fold q-values")
    println("─" ^ 120)
    for (label, nt, nf, nfeat, d, l, it) in [
        ("baseline",     50000, 100000, 15, 3,  10,  50),
        ("deep trees",   50000, 100000, 15, 10, 63,  50),
        ("small sample", 2000,  4000,   15, 3,  10,  50),
        ("small+deep",   2000,  4000,   15, 10, 63,  50),
        ("worst case",   2000,  4000,   50, 10, 63,  500),
    ]
        rn, rp, rpf = run_threeway(; nt, nf_total=nf, nfeat, sep=2.5, depth=d, leaves=l, iters=it)
        pr(@sprintf("%-18s NO CV", label), rn)
        pr(@sprintf("%-18s CV+POOLED qv", label), rp)
        pr(@sprintf("%-18s CV+PERFOLD qv", label), rpf)
        println()
    end

    # ── Multi-file with PEP calibration (Pioneer-like pipeline) ──
    println("─" ^ 120)
    println("EXP 8: Multi-file (5 files, PEP calibration + log-odds, Pioneer-like pipeline)")
    println("─" ^ 120)
    for (label, nt, nf, nfeat, d, l, it) in [
        ("baseline",     5000, 10000, 15, 3,  10,  50),
        ("deep trees",   5000, 10000, 15, 10, 63,  50),
        ("small sample", 1000, 2000,  15, 3,  10,  50),
        ("small+deep",   1000, 2000,  15, 10, 63,  50),
        ("worst case",   1000, 2000,  50, 10, 63,  500),
    ]
        rn = run_multifile(; n_true=nt, n_false=nf, nfeat, sep=2.5, depth=d, leaves=l, iters=it, cv_mode=:nocv)
        rp = run_multifile(; n_true=nt, n_false=nf, nfeat, sep=2.5, depth=d, leaves=l, iters=it, cv_mode=:cv_pooled)
        rpf = run_multifile(; n_true=nt, n_false=nf, nfeat, sep=2.5, depth=d, leaves=l, iters=it, cv_mode=:cv_perfold)
        pr(@sprintf("%-18s NO CV", label), rn)
        pr(@sprintf("%-18s CV+POOLED agg", label), rp)
        pr(@sprintf("%-18s CV+PERFOLD agg", label), rpf)
        println()
    end

    println("=" ^ 120)
end

main()
