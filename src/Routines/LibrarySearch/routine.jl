
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_detailed.jld2"  prosit_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_index_all.jld2"  prosit_index_all
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_precs.jld2"  prosit_precs

@time PSMs = SearchRAW(MS_TABLE, prosit_index, prosit_list_detailed, UInt32(1))
@time PSMs = SearchRAW(MS_TABLE, prosit_index, UInt32(1))

@time PSMs = SearchRAW(MS_TABLE, prosit_index_all, prosit_detailed, UInt32(1))

##########
#Get features
##########
function refinePSMs!(PSMs::DataFrame, precursors::Vector{LibraryPrecursor}; loss::AbstractEstimator = TauEstimator{TukeyLoss}(), maxiter = 200, min_spectral_contrast::AbstractFloat = 0.8)
    transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors[psm[:precursor_idx]])) => :decoy)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getIRT(precursors[psm[:precursor_idx]])) => :iRT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> MS_TABLE[:retentionTime][psm[:scan_idx]]) => :RT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:weight] == 0) => :nmf)
    function predictRTs!(PSMs::DataFrame; loss::AbstractEstimator = TauEstimator{TukeyLoss}(), maxiter = 200, min_spectral_contrast::AbstractFloat = 0.8)
    
        targets = PSMs[:,:decoy].==false
        spectral_contrast = PSMs[:,:decoy].>=min_spectral_contrast
    
        best_matches = targets .& spectral_contrast
    
        #Predicted iRTs
        iRTs = hcat(PSMs[:,:iRT][best_matches], ones(sum(best_matches)))
        #Emperical retention times
        RTs = PSMs[:,:RT][best_matches]
    
        slope, intercept = RobustModels.coef(rlm(iRTs, RTs, loss, initial_scale=:mad, maxiter = maxiter))
    
        transform!(PSMs, AsTable(:) => ByRow(psm -> abs((psm[:iRT]*slope + intercept) - psm[:RT])) => :RT_error)

    end
    predictRTs!(PSMs, loss = loss, maxiter = maxiter, min_spectral_contrast = min_spectral_contrast)
end
refinePSMs!(PSMs, prosit_precs)
########
#
########
function rankPSMs(PSMs::DataFrame)
    X = Matrix(PSMs[1:1:end,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RTdiff]])'
    X_labels = PSMs[:, :decoy]
    model = build_forest(X_labels, X', 4, 2000, 0.5, 3)
    probs = apply_forest_proba(model, X',[true, false])
    PSMs[:,:prob] = probs[:,2]
end

function getFDR(probs::Matrix{Float64}, threshold)
    sorted_probs = probs[sortperm(probs[:1]),:]
    PSMs[:m:decoy][sortperm(probs[:1])]
    targets = 0
    decoys = 0
    #use bisection algorithm
    N = length(sorted_probs)
    eachrow(sorted_probs[N:end,:])
    lo = Int64(N/2), hi = N
    while lo <= hi#eachrow(sorted_probs[N:end,:])
        sum(PSMs[:,:decoy][lo:end].>)
    end
end
histogram(Ylda[X_labels.==true], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)
histogram!(Ylda[X_labels.==false], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)

plot(PSMs[X_labels.==true,:spectral_contrast], PSMs[X_labels.==true,:RTdiff], alpha = 0.5, seriestype=:scatter)#, bins = -0.06:0.01:0.0)

plot!(PSMs[X_labels.==false,:spectral_contrast], PSMs[X_labels.==false,:RTdiff], alpha = 0.5, seriestype=:scatter)#, bins = -0.06:0.01:0.0)


sum(Ylda[X_labels.==true].<-0.0085)
sum(Ylda[X_labels.==false].<-0.0085)

model = build_forest(X_labels, X', 4, 2000, 0.5, 3)
probs = apply_forest_proba(model, X',[true, false])
PSMs[:,:prob] = probs[:,2]


histogram(PSMs[X_labels.==true,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[X_labels.==false,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)


sum(PSMs[X_labels.==true,:prob].>0.91)
sum(PSMs[X_labels.==false,:prob].>0.91)

unique(PSMs[PSMs[:,:prob].>0.91,:precursor_idx])




using Dictionaries
using CSV, Arrow, Tables, DataFrames, StatsBase
include("src/precursor.jl")
#include("src/buildFragmentIndex.jl")

