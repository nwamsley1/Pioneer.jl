
using Arrow, Tables, DataFrames, Dictionaries, Combinatorics, StatsBase, NMF, JLD2, LinearAlgebra, Random, DecisionTree, LoopVectorization, Splines2, ProgressBars, GLM, RobustModels, LoopVectorization
include("src/precursor.jl")
#include("src/buildFragmentIndex.jl")
include("src/Routines/LibrarySearch/buildFragmentIndex.jl")
include("src/Routines/LibrarySearch/matchpeaksLIB.jl")
include("src/Routines/LibrarySearch/buildDesignMatrix.jl")
include("src/Routines/LibrarySearch/spectralDistanceMetrics.jl")
include("src/Routines/LibrarySearch/spectralContrast.jl")
include("src/Routines/LibrarySearch/searchRAW.jl")
include("src/Routines/LibrarySearch/counter.jl")
#include("src/Routines/LibrarySearch/queryFragmentIndex.jl")
#include("src/Routines/LibrarySearch/queryFragmentInc.jl")
include("src/Routines/LibrarySearch/queryFragmentArr.jl")
include("src/PSM_TYPES/PSM.jl")
include("src/PSM_TYPES/LibraryXTandem.jl")
include("src/Routines/LibrarySearch/selectTransitions.jl")
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_detailed.jld2"  prosit_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_index_intensities.jld2"  prosit_index_intensities
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_precs.jld2"  prosit_precs

MS_TABLE = Arrow.Table("/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA.arrow")
prosit_totals = zeros(Float32, length(prosit_detailed))
i = 1
for prec in ProgressBar(prosit_detailed)
    prosit_totals[i] = sum([getIntensity(frag) for frag in prec]./norm([getIntensity(frag) for frag in prec]))
    i += 1
end
@time PSMs = SearchRAW(MS_TABLE, prosit_totals, prosit_index_intensities, prosit_detailed, UInt32(1), 
                        min_frag_count = 4, 
                        topN = 200, 
                        fragment_tolerance = 15.6, 
                        lambda = 1e5, 
                        max_peaks = 1000, 
                        scan_range = (0, 300000), 
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.65)
                        )
X, H, UNMATCHED, IDtoROW, fragmentMatches, fragmentMisses = SearchRAW(MS_TABLE, prosit_totals, prosit_index_intensities, prosit_detailed, UInt32(1), 
                        min_frag_count = 4, 
                        topN = 200, 
                        fragment_tolerance = 15.6, 
                        lambda = 1e5, 
                        max_peaks = 1000, 
                        scan_range = (101357, 101357), 
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.65)
                        )


function refinePSMs!(PSMs::DataFrame, precursors::Vector{LibraryPrecursor}; loss::AbstractEstimator = TauEstimator{TukeyLoss}(), maxiter = 200, min_spectral_contrast::AbstractFloat = 0.8)
    transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors[psm[:precursor_idx]])) => :decoy)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(precursors[psm[:precursor_idx]]))) => :iRT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:weight] < 10.0) => :nmf)
    function predictRTs!(PSMs::DataFrame; loss::AbstractEstimator = TauEstimator{TukeyLoss}(), maxiter = 200, min_spectral_contrast::AbstractFloat = 0.8)
        best_PSMs = combine(sdf -> sdf[argmax(sdf.spectral_contrast_all), :], groupby(PSMs, [:scan_idx]))
        targets = best_PSMs[:,:decoy].==false
        
        spectral_contrast = best_PSMs[:,:spectral_contrast_all].>=0.95#min_spectral_contrast
        
        best_matches = targets .& spectral_contrast
    
        #Predicted iRTs
        iRTs = best_PSMs[:,:iRT][best_matches]#,# ones(Float32, sum(best_matches)))
        #Emperical retention times
        RTs = best_PSMs[:,:RT][best_matches]
        #plot(iRTs, RTs, seriestype=:scatter)
        ns1 = Splines2.ns_(collect(range(-25.0, length=50, stop=160.0)),df=4,intercept=true);
        X = ns1(iRTs);
        fit1 = lm(X,RTs);
        PSMs[:,:RT_pred] =  GLM.predict(fit1, ns1(PSMs[:,:iRT]))
        PSMs[:,:RT_error] = abs.(PSMs[:,:RT_pred] .- PSMs[:,:RT])#abs.(PSMs[:,:RT] .- GLM.predict(fit1, ns1(PSMs[:,:iRT])))
        #slope, intercept = RobustModels.coef(rlm(iRTs, RTs, loss, initial_scale=:mad, maxiter = maxiter))
        #println("slope ", slope)
        #println("intercept ", intercept)
    
        #transform!(PSMs, AsTable(:) => ByRow(psm -> abs((psm[:iRT]*slope + intercept) - psm[:RT])) => :RT_error)
        return # RTs, iRTs
    end

    predictRTs!(PSMs, loss = loss, maxiter = maxiter, min_spectral_contrast = min_spectral_contrast)
    sort!(PSMs, [:scan_idx, :total_ions]);

    # Group DataFrame by "day" column
    grouped_df = groupby(PSMs, :scan_idx);


    PSMs[:,:next_best] = (combine(grouped_df) do sub_df
        pushfirst!(diff(sub_df.total_ions), zero(UInt32))
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    PSMs[:,:diff_hyper] = (combine(grouped_df) do sub_df
        sort!(sub_df, :hyperscore)
        pushfirst!(diff(sub_df.hyperscore), zero(Float64))
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    PSMs[:,:diff_scribe] = (combine(grouped_df) do sub_df
        sort!(sub_df, :scribe_score)
        pushfirst!(diff(sub_df.scribe_score), zero(Float64))
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    PSMs[:,:median_ions] = (combine(grouped_df) do sub_df
        #sort!(sub_df, :hyperscore)
        #pushfirst!(diff(sub_df.hyperscore), zero(Float64))
        repeat([median(sub_df.total_ions)], size(sub_df)[1])
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    grouped_df = groupby(PSMs, :precursor_idx);

    PSMs[:,:n_obs] = (combine(grouped_df) do sub_df
        #sort!(sub_df, :hyperscore)
        #pushfirst!(diff(sub_df.hyperscore), zero(Float64))
        repeat([size(sub_df)[1]], size(sub_df)[1])
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]
    transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(prosit_precs[psm[:precursor_idx]])) => :charge)

end
function rankPSMs!(PSMs::DataFrame; n_folds::Int = 3, n_trees::Int = 500, features::Int = 10, max_depth::Int = 10, fraction::AbstractFloat = 0.1)
   
    #X = Matrix(PSMs[:,[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:chebyshev,:city_block,:matched_ratio,:intensity,:count]])
    X = Matrix(PSMs[:,[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio]])
    X_labels = PSMs[:, :decoy]
    permutation = randperm(size(PSMs)[1])
    fold_size = length(permutation)÷n_folds

    folds = [((n-1)*fold_size + 1):(n*fold_size) for n in range(1, n_folds)]

    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    model = ""
    for test_fold_idx in range(1, n_folds)
        println(test_fold_idx)
        train_fold_idxs = vcat([folds[fold] for fold in range(1, length(folds)) if fold != test_fold_idx]...)
        #println("train_fold_idxs ", train_fold_idxs)
        train_features = X[train_fold_idxs,:]
        train_classes = X_labels[train_fold_idxs,1]
        #fraction = min(30000/length(X_labels), 0.15)
        model = build_forest(train_classes, train_features, features, n_trees, fraction, max_depth)
        probs = apply_forest_proba(model, X[folds[test_fold_idx],:],[true, false])
        PSMs[folds[test_fold_idx],:prob] = probs[:,2]
        println([:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio][sortperm(split_importance(model))])
    end
    #model = build_forest(X_labels, X', 4, 2000, 0.5, 3)
    #probs = apply_forest_proba(model, X',[true, false])
    #PSMs[:,:prob] = probs[:,2]
    return model
end
function getQvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool})
    #Could bootstratp to get more reliable values. 
    q_values = zeros(Float64, (length(probs),))
    order = reverse(sortperm(probs))
    targets = 0
    decoys = 0
    for i in order
        if labels[i] == true
            decoys += 1
        else
            targets += 1
        end
        q_values[i] = decoys/(targets + decoys)
    end
    PSMs[:,:q_values] = q_values;
end


refinePSMs!(PSMs, prosit_precs)


rt_index = buildRTIndex(best_psms)

@time rankPSMs!(PSMs, n_folds = 2, n_trees = 200, max_depth = 10, features = 10, fraction = 0.001)
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);

for i in [0.001, 0.002, 0.0005]
    @time rankPSMs!(PSMs, n_folds = 2, n_trees = 100, max_depth = 10, features = 10, fraction = i)
    @time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
    println("TEST ", length(unique(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :precursor_idx])))
end

for i in [5, 10, 15, 20]
    @time rankPSMs!(PSMs, n_folds = 2, n_trees = 100, max_depth = i, features = 10, fraction = 0.001)
    @time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
    println("TEST ", length(unique(PSMs[(PSMs[:,:q_values].<=0.1) .& (PSMs[:,:decoy].==false), :precursor_idx])))
end


for i in [3, 5, 7, 10]
    @time rankPSMs!(PSMs, n_folds = 2, n_trees = 100, max_depth = 10, features = i, fraction = 0.001)
    @time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
    println("TEST ", length(unique(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :precursor_idx])))
end

PSMs_RF = PSMs
length(unique(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :precursor_idx]))
@sk_import discriminant_analysis: (LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis)
@sk_import svm: SVC
PSMs = filter(row -> row.total_ions >=4, PSMs)
X = Matrix(PSMs[:,[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:RT,:diff_hyper,:diff_scribe,:city_block,:matched_ratio]])
X_labels = PSMs[:, :decoy]
replace!(X, NaN=>Inf)
replace!(X, -Inf=>0.0)
replace!(X, Inf=>0.0)
permutation = randperm(size(X)[1])
testnb = ScikitLearn.fit!(SVC(gamma=2, C=1), X[permutation[1:10000],:], X_labels[permutation[1:10000]])

@time model, coeffs = build_adaboost_stumps(X_labels, X, 10);
PSMs[:,:prob] = apply_adaboost_stumps_proba(model, coeffs, X, [true, false])[:,2]
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
length(unique(PSMs[(PSMs[:,:q_values].<=0.1) .& (PSMs[:,:decoy].==false), :precursor_idx]))

train_features = X[train_fold_idxs,:]
train_classes = X_labels[train_fold_idxs,1]
fraction = min(30000/length(X_labels), 0.15)
model = build_forest(train_classes, train_features, 10, 1000, fraction, 10)
probs = apply_forest_proba(model, X[folds[test_fold_idx],:],[true, false])

@time PSMs[:,:gnb] = ScikitLearn.predict_proba(testnb, X)[:,1]
@time PSMs[:,:gnb] = ScikitLearn.decision_function(testnb, X)
#@time getQvalues!(PSMs, PSMs[:,:gnb], PSMs[:,:decoy]);

stephist(PSMs[(PSMs[:,:decoy].==false), :prob], normalize = :pdf, alpha = 0.5)
stephist!(PSMs[(PSMs[:,:decoy].==true), :prob], normalize = :pdf, alpha = 0.5)

size(PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:prob].<0.2),:])
size(PSMs[(PSMs[:,:decoy].==true) .& (PSMs[:,:prob].<0.2),:])

size(PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:gnb].<0.05),:])
size(PSMs[(PSMs[:,:decoy].==true) .& (PSMs[:,:gnb].<0.05),:])

stephist(log2.(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :matched_ratio]), normalize = :pdf, alpha = 0.5)
stephist!(log2.(PSMs[(PSMs[:,:q_values].>0.1) .& (PSMs[:,:decoy].==true), :matched_ratio]),  normalize = :pdf, alpha = 0.5)


histogram(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :total_ions], normalize = :pdf, alpha = 0.5)
histogram!(PSMs[(PSMs[:,:q_values].>0.1) .& (PSMs[:,:decoy].==true), :total_ions],  normalize = :pdf, alpha = 0.5)


stephist(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :poisson], normalize = :pdf, alpha = 0.5)
stephist!(PSMs[(PSMs[:,:q_values].>0.1) .& (PSMs[:,:decoy].==true), :poisson],  normalize = :pdf, alpha = 0.5)


stephist(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :scribe_score], normalize = :pdf, alpha = 0.5)
stephist!(PSMs[(PSMs[:,:q_values].>0.1) .& (PSMs[:,:decoy].==true), :scribe_score],  normalize = :pdf, alpha = 0.5)


stephist(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :city_block], normalize = :pdf, alpha = 0.5)
stephist!(PSMs[(PSMs[:,:q_values].>0.1) .& (PSMs[:,:decoy].==true), :city_block],  normalize = :pdf, alpha = 0.5)

stephist(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false), :RT_error], normalize = :pdf, alpha = 0.5)
stephist!(PSMs[(PSMs[:,:q_values].>0.1) .& (PSMs[:,:decoy].==true), :RT_error],  normalize = :pdf, alpha = 0.5)

stephist(PSMs[(PSMs[:,:decoy].==false), :gnb], alpha = 0.5)
stephist!(PSMs[(PSMs[:,:decoy].==true), :gnb], alpha = 0.5)
#[:b_ladder, :charge, :y_ladder, :total_ions, :median_ions, :n_obs, :spectral_contrast_matched, :city_block, :diff_hyper, :spectral_contrast_all, :error, :diff_scribe, :hyperscore, :RT, :poisson, :intensity_explained, :matched_ratio, :scribe_score, :RT_error]
PSMs = 0
refinePSMs!(PSMs, prosit_precs)
X = Matrix(PSMs[:,[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio]])
replace!(X, NaN=>0.0)
replace!(X, -Inf=>0.0)
replace!(X, Inf=>0.0)


lda = fit(MulticlassLDA, X', X_labels; outdim=1)
X_labels = PSMs[:, :decoy]
Ylda = MultivariateStats.predict(lda, X')
transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(prosit_precs[psm[:precursor_idx]])) => :decoy)
  
test_ = ("target","decoy")
test_ = [test_[Int(x) + 1] for x in PSMs[:,:decoy]]
decoy_cv = glmnetcv(X, test_)

logit = round.(GLMNet.predict(decoy_cv, X, outtype = :prob), digits=3)

histogram(PSMs[PSMs[:,:decoy].==false,:logit], alpha = 0.5)
histogram!(PSMs[PSMs[:,:decoy].==true,:logit], alpha = 0.5)
sum(PSMs[PSMs[:,:logit].<0.4,:decoy])
sum(PSMs[PSMs[:,:logit].<0.4,:decoy].==false)
PSMs[:,:lda] = Ylda[1,:]
sum(PSMs[PSMs[:,:lda].>0.005,:decoy])

sum(PSMs[PSMs[:,:lda].>0.005,:decoy].==false)

@time PSMs = SearchRAW(MS_TABLE, prosit_index_intensities, prosit_detailed, UInt32(1), min_frag_count = 4, topN = i, fragment_tolerance = 15.6, lambda = 1e5, max_peaks = 1000, scan_range = (70000, 70000), precursor_tolerance = 20.0)
   
for i in [20]
    @time PSMs = SearchRAW(MS_TABLE, prosit_index_intensities, prosit_detailed, UInt32(1), min_frag_count = 4, topN = i, fragment_tolerance = 15.6, lambda = 1e5, max_peaks = 1000, scan_range = (0, 300000), precursor_tolerance = 20.0)
    refinePSMs!(PSMs, prosit_precs)
    rankPSMs!(PSMs, 2)
    @time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
    println("HITS ", length(unique(PSMs[(PSMs[:,:q_values].<0.01) .& (PSMs[:,:decoy].==false), :][:,:precursor_idx])))
end

combine(sdf -> sdf[argmax(sdf.prob), :], groupby(PSMs, [:precursor_idx]))
for i in [50, 100]
    @time PSMs = SearchRAW(MS_TABLE, prosit_index_intensities, prosit_detailed, UInt32(1), min_frag_count = 4, topN = i, fragment_tolerance = 15.6, lambda = 1e5, max_peaks = 1000, scan_range = (70001, 80000), precursor_tolerance = 30.0)
    refinePSMs!(PSMs, prosit_precs)
    rankPSMs!(PSMs, 2)
    @time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
    println("HITS ", length(unique(PSMs[(PSMs[:,:q_values].<0.01) .& (PSMs[:,:decoy].==false), :][:,:precursor_idx])))
end

for bin in prosit_index_intensities.precursor_bins
    for prec in bin.precs
        if getPrecID() == 

for i in [1000]
    @time PSMs = SearchRAW(MS_TABLE, prosit_index_intensities, prosit_detailed, UInt32(1), min_frag_count = 4, topN = i, fragment_tolerance = 15.6, lambda = 1e5, max_peaks = 1000, scan_range = (70001, 80000), precursor_tolerance = 30.0)
    refinePSMs!(PSMs, prosit_precs)
    rankPSMs!(PSMs, 2)
    @time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
    println("HITS ", length(unique(PSMs[(PSMs[:,:q_values].<0.01) .& (PSMs[:,:decoy].==false), :][:,:precursor_idx])))
end
ecdfplot(log2.(PSMs[(PSMs[:,:q_values].>0.01), :][:,:matched_ratio]))
ecdfplot!(log2.(PSMs[(PSMs[:,:q_values].<0.01) .& (PSMs[:,:decoy].==false), :][:,:matched_ratio]))
PSMs = 0 

length(unique(PSMs[(PSMs[:,:q_values].<0.1) .& (PSMs[:,:decoy].==false), :][:,:precursor_idx]))
hits_per_scan = combine(sdf->nrow(sdf), groupby(PSMs[PSMs[:,:q_values].<=0.05,:],:scan_idx))

prosit_totals = zeros(Float32, length(prosit_detailed))
i = 1
for prec in ProgressBar(prosit_detailed)
    prosit_totals[i] = sum([getIntensity(frag) for frag in prec]./norm([getIntensity(frag) for frag in prec]))
    i += 1
end

prosit_totals = zeros(Float32, length(prosit_detailed))
i = 1
for prec in ProgressBar(prosit_detailed)
    prosit_totals[i] = sum([getIntensity(frag) for frag in prec]./norm([getIntensity(frag) for frag in prec]))
    i += 1
end
2210

norm([getIntensity(x) for x in prosit_detailed[4211123]])

refinePSMs!(PSMs, prosit_precs)

@time PSMs_ = SearchRAW(MS_TABLE, prosit_index_intensities, prosit_detailed, UInt32(1), min_frag_count = 4, topN = 10, fragment_tolerance = 20.0, lambda = 1e4, max_peaks = 200, scan = 101358)
transform!( PSMs_, AsTable(:) => ByRow(psm -> Float64(getIRT(prosit_precs[psm[:precursor_idx]]))) => :iRT)
transform!( PSMs_, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
PSMs_[:,:RT_pred] = GLM.predict(fit1, ns1(PSMs_[:,:iRT]))
PSMs_[:,:RT_error] = abs.(PSMs_[:,:RT] - PSMs_[:,:RT_pred])
probs = apply_forest_proba(predict_targets_model, Matrix(PSMs_[:,[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast,:weight,:RT_error,:scribe_score,:y_ladder,:b_ladder]]), [true, false])
PSMs_[:,:prob] = probs[:,2]
PSMs_[:,:hit] = PSMs_[:,:prob].>=0.767
sort(PSMs_,:prob)
#RTs, iRTs = refinePSMs!(PSMs_40_tol10 , prosit_precs)

# Example DataFrame
df = DataFrame(day = [1, 1, 2, 2, 2, 3, 3, 3],
               score = [5, 10, 8, 12, 15, 4, 2, 6])

fit1, ns1 = refinePSMs!(PSMs, prosit_precs)
# Sort DataFrame by "day" and "score" columns


# Print the resulting DataFrame
println(df)
targets = PSMs[:,:decoy].==false
spectral_contrast = PSMs[:,:spectral_contrast].>=0.8#min_spectral_contrast

sum(spectral_contrast .& targets)

randperm(size(PS)[1])[1:60000]
test = combine(sdf -> sdf[argmin(sdf.q_value), :], groupby(PSMs, [:scan_idx]))
transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(prosit_precs[psm[:precursor_idx]])) => :charge)
########
#
########


best_PSMs = combine(sdf -> sdf[argmin(sdf.q_values), :], groupby(PSMs, [:precursor_idx]))


best = (best_PSMs[:,:q_values].<=0.01) .& (best_PSMs[:,:decoy].!=true)
plot(best_PSMs[best,:RT_pred], best_PSMs[best,:RT], seriestype=:scatter)
plot!([0, 150], [0, 150])

best_PSMs = combine(sdf -> sdf[argmin(sdf.q_values), :], groupby(PSMs, [:precursor_idx]))
best = (best_PSMs[:,:q_values].<=0.01) .& (best_PSMs[:,:decoy].!=false)

histogram(best_PSMs[best,:RT_error])

random = randperm(size(PSMs)[1])[1:sum(best)]
plot(PSMs[random,:iRT], PSMs[random,:RT], seriestype=:scatter)
@time rankPSMs!(PSMs)
train = randperm(size(PSMs)[1])[1:60000]
#train = randperm(size(PSMs)[1])[1:10000]
X = Matrix(PSMs[train,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error]])
X_labels = PSMs[train, :decoy]
model = build_forest(X_labels, X, 4, 1000, 0.5, 3)
X = Matrix(PSMs[:,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error]])
probs = apply_forest_proba(model, X,[true, false])
PSMs[:,:prob] = probs[:,2]
function getQvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool})
    #Could bootstratp to get more reliable values. 
    q_values = zeros(Float64, (length(probs),))
    order = reverse(sortperm(probs))
    targets = 0
    decoys = 0
    for i in order
        if labels[i] == true
            decoys += 1
        else
            targets += 1
        end
        q_values[i] = decoys/(targets + decoys)
    end
    PSMs[:,:q_values] = q_values;
end
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
length(unique(PSMs[PSMs[:,:q_values].<0.01,:][:,:precursor_idx]))
@time getQvalues!(PSMs_40, PSMs_40[:,:prob], PSMs_40[:,:decoy]);


@time getQvalues!(PSMs_40_tol10, PSMs_40_tol10[:,:prob], PSMs_40_tol10[:,:decoy]);



histogram(PSMs[PSMs[:,:decoy].==true,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[PSMs[:,:decoy].==false,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)


histogram(log2.(PSMs[PSMs[:,:decoy].==true,:int_ion]), alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(log2.(PSMs[PSMs[:,:decoy].==false,:int_ion]), alpha = 0.5)#, bins = -0.06:0.01:0.0

histogram(log2.(PSMs[PSMs[:,:decoy].==true,:total_ions]), alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(log2.(PSMs[PSMs[:,:decoy].==false,:total_ions]), alpha = 0.5)#, bins = -0.06:0.01:0.0


histogram(PSMs[PSMs[:,:decoy].==false,:q_values], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[PSMs[:,:decoy].==true,:q_values], alpha = 0.5)#, bins = -0.06:0.01:0.0)



histogram(PSMs[PSMs[:,:decoy].==true,:], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[PSMs[:,:decoy].==false,:q_values], alpha = 0.5)#, bins = -0.06:0.01:0.0)

histogram(PSMs[X_labels.==true,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[X_labels.==false,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)

value_counts(df, col) = combine(groupby(df, col), nrow)
counts = sort(value_counts(PSMs[PSMs[:,:q_values].<=0.01,:], :precursor_idx),:nrow)
unique(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==true),:][:,:precursor_idx])
unique(PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false),:][:,:precursor_idx])


PSMs[PSMs[:,:precursor_idx] .== 3768665,[:weight,:RT]]

a = sort(PSMs[PSMs[:,:precursor_idx] .==    3619976 ,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)


a = sort(PSMs[PSMs[:,:precursor_idx] .==  0x00382352,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)

a = sort(PSMs[PSMs[:,:precursor_idx] .==   4466503 ,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)

a = sort(PSMs[PSMs[:,:precursor_idx] .==   774613,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)

a = sort(PSMs[PSMs[:,:precursor_idx] .==     0x003d423a,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)

a = sort(PSMs[PSMs[:,:precursor_idx] .==  0x0059a963,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)

sort(counts,:nrow)

a = sort(PSMs[PSMs[:,:precursor_idx] .==  0x00245ae9 ,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)

a = sort(PSMs[PSMs[:,:precursor_idx] .==  0x0009d576,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)

a = sort(PSMs[PSMs[:,:precursor_idx] .==        0x000b848d,[:weight,:RT]], :RT)
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)


a = sort(PSMs[PSMs[:,:precursor_idx] .==    0x0030a4bd   ,[:weight,:RT]], :RT);
plot(a[:,:RT], a[:,:weight], seriestype=:scatter)


sum(PSMs[PSMs[:,:decoy].==true,:prob].>0.89)
sum(PSMs[PSMs[:,:decoy].==false,:prob].>0.89)


sum(PSMs[PSs[:,:decoy].==true,:q_values].<=0.01)
sum(PSMs[PSMs[:,:decoy].==false,:q_values].<=0.01)

sum(PSMs[PSMs[:,:decoy].==true,:int_ion].<=0.01)
sum(PSMs[PSMs[:,:decoy].==false,:int_ion].<=0.01)

sum(PSMs[PSMs[:,:decoy].==false,:int_ion].>0.02)
sum(PSMs[PSMs[:,:decoy].==true,:int_ion].>0.02)

PSMs[(PSMs[:,:q_values].<=0.01) .& (PSMs[:,:decoy].==false),:precursor_idx]

PSMs_40[(PSMs_40[:,:q_values].<=0.01) .& (PSMs_40[:,:decoy].==false),:precursor_idx]


X = Matrix(PSMs[1:1:end,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RTdiff]])
X_labels = PSMs[:, :decoy]
@time model = build_forest(X_labels, X, 4, 100, 0.5, 3)
probs = apply_forest_proba(model, X',[true, false])
PSMs[:,:prob] = probs[:,2]
PSMs[:,:q_value] = getQvalues(PSMs[:,:prob], PSMs[:,:decoy])
histogram(PSMs[PSMs[:,:decoy].==true,:q_values], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[PSMs[:,:decoy].==false,:q_values], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)

histogram(PSMs[PSMs[:,:decoy].==true,:q_values], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[PSMs[:,:decoy].==false,:q_values], alpha = 0.5)#, bins = 

histogram(PSMs[PSMs[:,:decoy].==true,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[PSMs[:,:decoy].==false,:prob], alpha = 0.5)#, bins = 

function getQvalues(probs::Vector{Float64}, labels::Vector{Bool}, bootstrap_n::Int = 1000)
    #Could bootstratp to get more reliable values. 
    order = reverse(sortperm(probs))
    q_values = zeros(Float64, (length(probs),))
    targets = 0
    decoys = 0
    for n in boostrap_n
        println(n)
        for i in order
            if labels[i] == true
                decoys += 1
            else
                targets += 1
            end
            q_values[i] += decoys/(targets + decoys)
        end
    end
    return q_values./bootstrap_n
end
PSMs[:,:q_values] = getQvalues(PSMs[:,:prob], PSMs[:,:decoy])

histogram(Ylda[X_labels.==true], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)
histogram!(Ylda[X_labels.==false], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)

plot(PSMs[X_labels.==true,:spectral_contrast], PSMs[X_labels.==true,:RTdiff], alpha = 0.5, seriestype=:scatter)#, bins = -0.06:0.01:0.0)

plot!(PSMs[X_labels.==false,:spectral_contrast], PSMs[X_labels.==false,:RTdiff], alpha = 0.5, seriestype=:scatter)#, bins = -0.06:0.01:0.0)


sum(Ylda[X_labels.==true].<-0.0085)
sum(Ylda[X_labels.==false].<-0.0085)

model = build_forest(X_labels, X', 4, 2000, 0.5, 3)
probs = apply_forest_proba(model, X',[true, false])
PSMs[:,:prob] = probs[:,2]


sortPSMs = sort(PSMs, :prob)

for 

histogram(PSMs[X_labels.==true,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[X_labels.==false,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)


sum(PSMs[X_labels.==true,:prob].>0.92)
sum(PSMs[X_labels.==false,:prob].>0.92)

unique(PSMs[PSMs[:,:prob].>0.92,:precursor_idx])




using Dictionaries
using CSV, Arrow, Tables, DataFrames, StatsBase
include("src/precursor.jl")
#include("src/buildFragmentIndex.jl")
ms2s = MS_TABLE[:msOrder].==2

mean([length(MS_TABLE[:masses][i]) for (i, order) in enumerate(MS_TABLE[:msOrder]) if (order == 2)])


sort([length(MS_TABLE[:masses][i]) for (i, order) in enumerate(MS_TABLE[:msOrder]) if (order == 1)])
[i for i in  MS_TABLE[:msOrder]]

test = ones(UInt32, 9387261)

function makezeros(list::Vector{UInt32})
    for i ∈ eachindex(list)
        if list[i] == 1
            list[i] = UInt32(0)
        end
    end
end
test = ones(UInt32, 9387261÷2)
filter(x->x==0, test)

pre_aloc_size = 1000
test_dict = Dict{UInt32, UInt8}(zeros(UInt8, pre_aloc_size), zeros(UInt32, pre_aloc_size), zeros(UInt8, pre_aloc_size), 0, 0, 0, pre_aloc_size, 0)

struct MyKey
    val::UInt
end

import Base.hash 
Base.hash(a::UInt32, h::UInt32) = xor(a, h)
#Base.==(a::UInt, b::UInt) = a.val == b.val
import Base.empty!
empty!(h::Accumulator{K,V}) where V where K = empty!(h.map)


test_dict = Dict{UInt32, UInt32}()
ints = [UInt32(x) for x in range(1, 200000)]

test_dict = Accumulator{UInt32, UInt32}()
for i in 1:1000
    for i in vals
        inc!(test_dict, i)#test_dict[i] = UInt32(1)
    end
    empty!(test_dict)
end


function makedicttest1(vals::Vector{UInt32})
    test_dict = Dict{UInt32, UInt32}()
    for i in 1:1
        #test_dict = Dict{UInt32, UInt32}()
        for i in vals
            if !haskey(test_dict, i)
                test_dict[i] = UInt32(1)
            else
                test_dict[i] += UInt32(1)
            end
        end
    end
    return test_dict
end

function makedicttest2(vals::Vector{UInt32})
    test_dict = Dict{UInt32, UInt32}()
    for i in 1:1000
        for i in vals
            if !haskey(test_dict, i)
                test_dict[i] = UInt32(1)
            else
                test_dict[i] += UInt32(1)
            end
        end
        empty!(test_dict)
    end
    return test_dict
end

function makedicttest3(vals::Vector{UInt32})
    test_dict = Accumulator{UInt32, UInt32}()
    for i in 1:1000
        for i in vals
            inc!(test_dict, i)#test_dict[i] = UInt32(1)
        end
        empty!(test_dict)
    end
end

function makedicttest4(vals::Vector{UInt32})
    test_dict = SortedDict{UInt32, UInt32}()
    for i in 1:1000
        for i in vals
            if !haskey(test_dict, i)
                test_dict[i] = UInt32(1)
            else
                test_dict[i] += UInt32(1)
            end
        end
        empty!(test_dict)
    end
end

test_acc = Accumulator{UInt32, UInt32}()
    for i in rand(UInt32, 1, 10000)
        for n in range(1, rand([1, 2, 3, 4, 5, 6]))
            inc!(test_acc, i)#test_dict[i] = UInt32(1)
        end
    end

sort!(test_acc)
test_iterator = (first(pair) for pair in sort(filter(pair -> val(pair) > 1, pairs(test_acc)), by = pair -> last(pair)))

@time begin
    makedicttest(ints)
end

(first(pair) for pair in sort(filter(kv -> last(kv) > 1, pairs(test_acc)), by = x -> last(x)))

import Base.hash
Base.hash(x::UInt32) = hash(x, zero(UInt))
Base.hash(x::UInt32, h::UInt) = x #hash(Int64(x), h)
#hash(x::Any) = hash(x, zero(UInt))
import Base.hashindex
function Base.hashindex(key::UInt32, sz)
    hsh = hash(key)
    idx = (((hsh % Int) & (sz-1)) + 1)::Int
    return idx, _shorthash7(hsh)
end
import Base._shorthash7
_shorthash7(hsh::UInt32) = (hsh >> UInt(25))%UInt8 | 0x80
import Base.ht_keyindex
function Base.ht_keyindex(h::Dict{K,V}, key) where V where K
    isempty(h) && return -1
    sz = length(h.keys)
    iter = 0
    maxprobe = h.maxprobe
    index, sh = hashindex(key, sz)
    keys = h.keys
    println("index ", index)
    println("sh ", sh)
    @inbounds while true
        isslotempty(h,index) && return -1
        if h.slots[index] == sh
            k = keys[index]
            if (key ===  k || isequal(key, k))
                return index
            end
        end

        index = (index & (sz-1)) + 1
        (iter += 1) > maxprobe && return -1
    end
    # This line is unreachable
end

precs.size 54622
frag_min929.4830146472168
frag_max 929.5016044934083
frag_bin 96106
frag_bin 96107
frag_bin 96108
prec_min 775.8543701171875
prec_max  784.3543701171875
precs.size 54703
frag_min 930.4836149926758
frag_max 930.5022248510743
frag_bin 96206
frag_bin 96207
frag_bin 96208
prec_min 775.8543701171875
prec_max  784.3543701171875
