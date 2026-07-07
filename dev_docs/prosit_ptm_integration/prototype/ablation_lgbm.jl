# Phase-2 ablation: does site-determining S add signal beyond frac + C(c,n)?
# Trains an LGBM target-vs-decoy classifier on the SDS feature table; 70/30 split; reports AUC and
# held-out target IDs at matched sample decoy-FLR (D/R). NOTE: sample is ~1:1 target:decoy, so the
# absolute D/R is not the calibrated FLR -- the ablation *comparison* is what's meaningful.
using Arrow, DataFrames, Statistics, Random, LightGBM

feat = DataFrame(Arrow.Table(ARGS[1]))
feat = feat[.!isnan.(feat.S_vs_targets), :]
feat.logC = log.(Float64.(feat.C))
y = Float64.(.!feat.is_loc_decoy)                    # 1 = target, 0 = decoy

sets = ["frac"           => [:frac],
        "frac+C"         => [:frac, :C, :logC],
        "frac+C+S"       => [:frac, :C, :logC, :S_vs_targets, :self_frac, :gsz]]  # n_comp dropped (leaks label via C)

rng = MersenneTwister(7); perm = shuffle(rng, 1:nrow(feat))
ntr = round(Int, 0.7*length(perm)); tr = perm[1:ntr]; te = perm[ntr+1:end]

# AUC via Mann-Whitney (P(score_target > score_decoy))
function auc(scores, labels)
    pos = scores[labels .== 1]; neg = scores[labels .== 0]
    isempty(pos) || isempty(neg) && return NaN
    c = 0.0; for p in pos, n in neg; c += p>n ? 1.0 : (p==n ? 0.5 : 0.0); end
    c / (length(pos)*length(neg))
end
# target IDs on test at decoy-FLR (D/R) <= tgt
function ids_at(scores, labels, tgt)
    o = sortperm(scores, rev=true); R=0; D=0; best=0
    for i in o
        labels[i]==1 ? (R+=1) : (D+=1)
        if R>0 && D/R <= tgt; best=max(best,R); end
    end
    best
end

println("train=",ntr,"  test=",length(te),"  (targets ",Int(sum(y)),", decoys ",Int(length(y)-sum(y)),")\n")
println(rpad("features",14), rpad("AUC",8), rpad("test IDs@10%",14), "test IDs@5%")
for (name, cols) in sets
    X = Matrix{Float64}(feat[:, cols])
    est = LGBMClassification(objective="binary", num_class=1, num_iterations=300, learning_rate=0.05,
                             num_leaves=31, min_data_in_leaf=30, feature_fraction=0.9, bagging_fraction=0.9,
                             bagging_freq=1, metric=["auc"])
    LightGBM.fit!(est, X[tr,:], y[tr]; verbosity=-1)
    pr = LightGBM.predict(est, X[te,:]); p = pr isa AbstractMatrix ? vec(pr) : pr
    a = auc(p, y[te]); i10 = ids_at(p, y[te], 0.10); i5 = ids_at(p, y[te], 0.05)
    println(rpad(name,14), rpad(round(a,digits=3),8), rpad(i10,14), i5)
end
