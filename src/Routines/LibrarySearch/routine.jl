
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_detailed.jld2"  prosit_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_index_all.jld2"  prosit_index_all
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_precs.jld2"  prosit_precs

@time PSMs = SearchRAW(MS_TABLE, prosit_index, prosit_list_detailed, UInt32(1))
@time PSMs = SearchRAW(MS_TABLE, prosit_index, UInt32(1))

@time test_table = SearchRAW(MS_TABLE, prosit_index_all, prosit_detailed, UInt32(1))
PSMs = test_table
transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(prosit_precs[psm[:precursor_idx]])) => :decoy)
transform!(PSMs, AsTable(:) => ByRow(psm -> getIRT(prosit_precs[psm[:precursor_idx]])) => :iRT)
transform!(PSMs, AsTable(:) => ByRow(psm -> MS_TABLE[:retentionTime][psm[:scan_idx]]) => :RT)
transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:weight] == 0) => :nmf)


best_PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, [:scan_idx])) 
#transform!(best_PSMs, AsTable(:) => ByRow(psm -> MS_TABLE[:retentionTime][psm[:scan_idx]]) => :RT)
plot(best_PSMs[:,:iRT][best_PSMs[:,:decoy].==false], best_PSMs[:,:RT][best_PSMs[:,:decoy].==false], seriestype = :scatter)
plot!(best_PSMs[:,:iRT][best_PSMs[:,:decoy].==true], best_PSMs[:,:RT][best_PSMs[:,:decoy].==true], seriestype = :scatter)


function predictRTs(PSMs::DataFrame, precursors::Vector{LibraryPrecursor}; loss::AbstractEstimator = TauEstimator{TukeyLoss}(), maxiter = 200, min_spectral_contrast::AbstractFloat = 0.8)
    targets = PSMs[:,:decoy].==false
    spectral_contrast = PSMs[:,:decoy].>=min_spectral_contrast

    best_matches = targets .& spectral_contrast

    iRTs = hcat(PSMs[:,:iRT][best_matches], ones(sum(best_matches)))
    RTs = PSMs[:,:RT][best_matches]

    return RobestModels.coef(rlm(iRTs, RTs, loss, initial_scale=:mad, maxiter = maxiter))
end

iRTs = best_PSMs[:,:iRT][(best_PSMs[:,:decoy].==false) .& (best_PSMs[:,:spectral_contrast].>=0.7)]
iRTs = reshape(iRTs, (length(iRTs), 1))
RTs = best_PSMs[:,:RT][(best_PSMs[:,:decoy].==false) .& (best_PSMs[:,:spectral_contrast].>=0.7)]
RTs  = [Float64(x) for x in RTs]
iRTs= hcat(iRTs, ones(length(iRTs)))
rlm(iRTs, RTs, TauEstimator{TukeyLoss}(), initial_scale=:mad, maxiter = 200)


plot(iRTs, RTs, seriestype = :scatter)
plot!([0, 150], [18.6381, 150*0.263446 + 18.6381])

transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:iRT]*0.263446 + 18.6381) => :predRT)

transform!(PSMs, AsTable(:) => ByRow(psm -> abs(psm[:RT] - psm[:predRT])) => :RTdiff)
X = Matrix(PSMs[1:1:end,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RTdiff]])'
X_labels = Vector(PSMs[1:1:end, :decoy])
lda = fit(MulticlassLDA, X, X_labels; outdim=1)
Ylda = predict(lda, X)
PSMs[:,:score] = Ylda[1,:]
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

