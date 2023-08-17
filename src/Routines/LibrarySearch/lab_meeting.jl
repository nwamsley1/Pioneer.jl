
##########
#Import Libraries
##########
using CSV, Arrow, Tables, DataFrames, Dictionaries, Combinatorics, StatsBase, NMF, JLD2, LinearAlgebra, Random, DecisionTree, LoopVectorization, Splines2, ProgressBars, GLM, RobustModels, LoopVectorization, SparseArrays, Interpolations, MultiKDE, XGBoost, SavitzkyGolay, NumericalIntegration
##########
#Import files
##########
include("src/precursor.jl")
#include("src/Routines/LibrarySearch/parsePrositLib.jl")
include("src/Routines/LibrarySearch/buildFragmentIndex.jl")
include("src/Routines/LibrarySearch/matchpeaksLib.jl")
include("src/Routines/LibrarySearch/buildDesignMatrix.jl")
include("src/Routines/LibrarySearch/spectralDistanceMetrics.jl")
include("src/Routines/LibrarySearch/searchRAW.jl")
include("src/Routines/LibrarySearch/counter.jl")
include("src/Routines/LibrarySearch/ML.jl")
include("src/Routines/LibrarySearch/refinePSMs.jl")
include("src/Routines/LibrarySearch/buildRTIndex.jl")
include("src/Routines/LibrarySearch/selectTransitions.jl")
include("src/Routines/LibrarySearch/NMF.jl")
include("src/Routines/LibrarySearch/integratePrecursors.jl")
##########
include("src/Routines/LibrarySearch/queryFragmentArr.jl")
include("src/PSM_TYPES/PSM.jl")
include("src/PSM_TYPES/LibraryXTandem.jl")


@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_detailed.jld2"  frag_detailed
@load  "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_15irt.jld2" prosit_index_5ppm_15irt 
#@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/precursor_list.jld2" precursor_list
#@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/prosit_index_5ppm_15irt.jld2" prosit_index_5ppm_15irt
#new_prosit_index =  prosit_index_5ppm_15irt
#@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/frags_detailed.jld2" frags_detailed
#@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/precursors.jld2" precursors
#precursor_list = precursors

#@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_simple_33NCEdynamic.jld2" frags_mouse_simple_33NCEdynamic
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_detailed_33NCEdynamic.jld2" frags_mouse_detailed_33NCEdynamic
#@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/precursors_mouse_detailed_33NCEdynamic.jld2" precursors_mouse_detailed_33NCEdynamic
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/prosit_mouse_33NCEdynamic_5ppm_15irt.jld2" prosit_mouse_33NCEdynamic_5ppm_15irt


@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_detailed_33NCEfixed.jld2" frags_mouse_detailed_33NCEfixed
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/precursors_mouse_detailed_33NCEfixed.jld2" precursors_mouse_detailed_33NCEfixed
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/prosit_mouse_33NCEfixed_5ppm_15irt.jld2" prosit_mouse_33NCEfixed_5ppm_15irt



@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_detailed_33NCEcorrected_start1.jld2" frags_mouse_detailed_33NCEcorrected_start1
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/precursors_mouse_detailed_33NCEcorrected_start1.jld2" precursors_mouse_detailed_33NCEcorrected_start1
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/prosit_mouse_33NCEcorrected_start1_5ppm_15irt.jld2" prosit_mouse_33NCEcorrected_start1_5ppm_15irt
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/linear_spline.jld2" linear_spline
MS_TABLE = Arrow.Table("/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA.arrow")

#@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_detailed_33NCEcorrected_chronologer.jld2" frags_mouse_detailed_33NCEcorrected_chronologer
#@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/precursors_mouse_detailed_33NCEcorrected_chronologer.jld2" precursors_mouse_detailed_33NCEcorrected_chronologer
#@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/prosit_mouse_33NCEcorrected_chronologer_5ppm_15irt.jld2" prosit_mouse_33NCEcorrected_chronologer_5ppm_15irt
#@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/rt_to_hi_spline.jld2" rt_to_hi_spline
#MS_TABLE = Arrow.Table("/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA.arrow")

#@load  "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_15irt.jld2" prosit_index_5ppm_15irt
#@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_detailed.jld2"  frag_detailed
###########
#Main Search
###########
   PSMs = SearchRAW(MS_TABLE, prosit_index_5ppm_15irt, frags_detailed, UInt32(1), linear_spline,
                        min_frag_count = 4, 
                        topN = 200, 
                        fragment_tolerance = 15.6, 
                        λ = Float32(1e3), 
                        γ =Float32(1),
                        max_peaks = 1000, 
                        scan_range = (0, 300000), #101357 #22894
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.5),
                        min_matched_ratio = Float32(0.45),
                        rt_tol = Float32(100.0)
                        )
newPSMs = SearchRAW(MS_TABLE, prosit_33NCEfixed_5ppm_15irt, frags_detailed_33NCEfixed, UInt32(1), linear_spline,
                        min_frag_count = 4, 
                        topN = 200, 
                        fragment_tolerance = 15.6, 
                        λ = Float32(1e3), 
                        γ =Float32(1),
                        max_peaks = 1000, 
                        scan_range = (0, 300000), #101357 #22894
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.5),
                        min_matched_ratio = Float32(0.45),
                        rt_tol = Float32(20.0)
                        )

newPSMs = SearchRAW(MS_TABLE, prosit_mouse_35NCEcorrected_start3_5ppm_15irt,  frags_mouse_detailed_35NCEcorrected_start3, UInt32(1), linear_spline,
                        min_frag_count = 4, 
                        topN = 200, 
                        fragment_tolerance = 15.6, 
                        λ = Float32(1e3), 
                        γ =Float32(1),
                        max_peaks = 1000, 
                        scan_range = (101357, 101357), #101357 #22894
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.5),
                        min_matched_ratio = Float32(0.45),
                        rt_tol = Float32(20.0)
                        )

newPSMs = SearchRAW(MS_TABLE, prosit_mouse_33NCEcorrected_chronologer_5ppm_15irt,  frags_mouse_detailed_33NCEcorrected_chronologer, UInt32(1), rt_to_hi_spline,
                        min_frag_count = 4, 
                        topN = 1000, 
                        fragment_tolerance = 15.6, 
                        λ = Float32(1e3), 
                        γ =Float32(1),
                        max_peaks = 10000, 
                        scan_range = (0, 300000), #101357 #22894
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.5),
                        min_matched_ratio = Float32(0.45),
                        rt_tol = Float32(20.0)
                        )

PSMs = SearchRAW(MS_TABLE, prosit_index_5ppm_15irt, frag_detailed, UInt32(1), linear_spline,
                        min_frag_count = 4, 
                        topN = 200, 
                        fragment_tolerance = 15.6, 
                        λ = Float32(1e3), 
                        γ =Float32(1),
                        max_peaks = 1000, 
                        scan_range = (101357, 101357), #101357 #22894
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.5),
                        min_matched_ratio = Float32(0.45),
                        rt_tol = Float32(20.0)
                        )


#=
current_peptide _GGPGSAVSPYPSFNVSSDVAALHK_
CSV.Row2{Any, PosLenString}:
 :RelativeIntensity  "0.09786233305931091"
 :FragmentMz         "147.11280822753906"
 :ModifiedPeptide    "_GGPGSAVSPYPSFNVSSDVAALHK_"
 :PrecursorCharge    "3"
 :PrecursorMz        "782.0498415566667"
 :iRT                "88.7706298828125"
 :FragmentNumber     "1"
 :FragmentType       "y"
 :FragmentCharge     "1"
=#
CSV.write("/Users/n.t.wamsley/Desktop/PSMs_080923.csv", PSMs)
PSMs_old = PSMs


tolerance = Float64[3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 100]
params_est = Vector{Tuple{Float64, Float64}}()
tolerance = Float64[30]
for tol in tolerance
    rtPSMs, all_matches = SearchRAW(MS_TABLE, 
                        prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
                        frags_mouse_detailed_33NCEcorrected_start1, 
                        UInt32(1), 
                        x->x, #Mapp RT to iRT
                        min_frag_count = 7, 
                        topN = 5, 
                        fragment_tolerance = tol,#20.0, 
                        λ = Float32(0), 
                        γ =Float32(0),
                        max_peaks = 10000, 
                        scan_range = (0, length(MS_TABLE[:scanNumber])), #All Scans
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.95),
                        min_matched_ratio = Float32(.6),
                        rt_tol = Float32(1e6), #Set arbitrarily high
                        sample_rate = 0.01, #Sampling rate
                        frag_ppm_err = 0.0,
                        collect_frag_errs = true
                        )

    frag_ppm_errs = [getPPM(x.theoretical_mz, x.match_mz) for x in all_matches]
    mix_guess = MixtureModel([Normal(0, 3), Uniform(-100, 100)], [0.5, 1 - 0.5])
    mix_mle = fit_mle(mix_guess, (frag_ppm_errs); display = :iter, atol = 1e-5, robust = true, infos = false)
    push!(params_est, (params(mix_mle)[1][1][1], params(mix_mle)[1][1][2]))
end

#Plots.plot(tolerance, [first(x) for x in params_est], seriestype=:scatter)
#Plots.plot(tolerance, [last(x) for x in params_est], seriestype=:scatter)


transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(precursors_mouse_detailed_33NCEcorrected_start1[psm[:precursor_idx]]))) => :iRT)
transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
transform!(rtPSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors_mouse_detailed_33NCEcorrected_start1[psm[:precursor_idx]])) => :decoy)
rtPSMs = rtPSMs[rtPSMs[:,:decoy].==false,:]
p = Plots.plot(rtPSMs[:,:RT], rtPSMs[:,:iRT], seriestype=:scatter,
            xlabel = "Retention Time RT (min)",
            ylabel ="Indexed Retention Time iRT (min)",
            label = nothing,
            size = 100*[13.3, 7.5]
,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm)

#Plots.plot(best_psms[best_psms[:,:q_value].<=0.01,:RT], best_psms[best_psms[:,:q_value].<=0.01,:iRT], seriestype=:scatter)
rt_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT], n = 50)
Plots.plot!(p, LinRange(minimum(rtPSMs[:,:RT]), maximum(rtPSMs[:,:RT]), 100), 
            rt_map.(LinRange(minimum(rtPSMs[:,:RT]), maximum(rtPSMs[:,:RT]), 100)),
            lw = 6.0,
            label = "RT Spline")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/rt_spline.pdf")
function getPPM(a::T, b::T) where {T<:AbstractFloat}
    (a-b)/(a/1e6)
end

frag_ppm_errs = [getPPM(x.theoretical_mz, x.match_mz) for x in all_matches]

mix_guess = MixtureModel([Normal(0, 3), Uniform(-100, 100)], [0.5, 1 - 0.5])
mix_mle = fit_mle(mix_guess, (frag_ppm_errs); display = :iter, atol = 1e-5, robust = true, infos = false)

b = Distributions.Uniform(params(mix_mle)[1][2][1], params(mix_mle)[1][2][2])
a = Distributions.Normal(params(mix_mle)[1][1][1], params(mix_mle)[1][1][2])
c = Distributions.Bernoulli(params(mix_mle)[2][1])
dist = zeros(Float64, 30000)
for i in eachindex(dist)
    w = rand(c)*1.0
    dist[i] = w*rand(a) + (1 - w)*rand(b)
end


p = Plots.histogram(dist, normalize = :probability, bins = 100, alpha = 0.5, label = "Samples from Fitted Dist.",
xlabel = "Mass Error (ppm)", ylabel = "Probability",
size = 100*[13.3, 7.5]
,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 18,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm, legend = :topleft)
Plots.histogram!(p, (frag_ppm_errs), normalize = :probability, bins = 100, alpha = 0.5, label = "Observed Mass Errors (ppm)")
#Plots.vline!([median(frag_ppm_errs)])
Plots.vline!(p, [params(mix_mle)[1][1][1]], lw = 6.0, color = :black, label = "Estimated Mass Error (ppm)")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/mass_errors.pdf")



mix_guess = MixtureModel([Normal(0, 3), Laplace(0, 3), Uniform(-100, 100)], [0.3, 0.3, 0.4])
mix_mle = fit_mle(mix_guess, (frag_ppm_errs); display = :iter, atol = 1e-5, robust = true, infos = false)

a = Distributions.Normal(params(mix_mle)[1][1][1], params(mix_mle)[1][1][2])
b = Distributions.Laplace(params(mix_mle)[1][1][1], params(mix_mle)[1][1][2])
c = Distributions.Uniform(params(mix_mle)[1][2][1], params(mix_mle)[1][2][2])
d = Distributions.Bernoulli(params(mix_mle)[2][1])
dist = zeros(Float64, 30000)
for i in eachindex(dist)
    w = rand(c)*1.0
    dist[i] = w*rand(a) + (1 - w)*rand(b)
end
Plots.histogram(dist, normalize = :pdf, bins = 100, alpha = 0.5)
Plots.histogram!((frag_ppm_errs), normalize = :pdf, bins = 100, alpha = 0.5)
Plots.vline!([median(frag_ppm_errs)])
Plots.vline!([params(mix_mle)[1][1][1]])



frag_ppm_errs = [getPPM(x.theoretical_mz, x.match_mz) for x in all_matches]

mix_guess = MixtureModel([Normal(0, 3), Laplace(0, 100)], [0.5, 1 - 0.5])
mix_mle = fit_mle(mix_guess, (frag_ppm_errs); display = :iter, atol = 1e-5, robust = true, infos = false)

b = Distributions.Normal(params(mix_mle)[1][2][1], params(mix_mle)[1][2][2])
a = Distributions.Laplace(params(mix_mle)[1][1][1], params(mix_mle)[1][1][2])
c = Distributions.Bernoulli(params(mix_mle)[2][1])
dist = zeros(Float64, 30000)
for i in eachindex(dist)
    w = rand(c)*1.0
    dist[i] = w*rand(a) + (1 - w)*rand(b)
end
Plots.histogram(dist, normalize = :pdf, bins = 50, alpha = 0.5)
Plots.histogram!((frag_ppm_errs), normalize = :pdf, bins = 50, alpha = 0.5)


test_dist = Normal(3.242241988271811, 2.923584485136043)




fragment_tolerance = params(mix_mle)[1][1][2]*3
#fragment_tolerance = 15.6
@time begin
    newPSMs = SearchRAW(MS_TABLE, prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  frags_mouse_detailed_33NCEcorrected_start1, UInt32(1), rt_map,
                            min_frag_count = 4, 
                            topN = 1000, 
                            fragment_tolerance = fragment_tolerance, 
                            λ = Float32(0), 
                            γ =Float32(0),
                            max_peaks = 10000, 
                            scan_range = (0, 300000), #101357 #22894
                            precursor_tolerance = 20.0,
                            min_spectral_contrast =  Float32(0.5),
                            min_matched_ratio = Float32(0.45),
                            rt_tol = Float32(20.0),
                            frag_ppm_err = 3.34930002879957
                            )

    PSMs = newPSMs
    PSMs = PSMs[PSMs[:,:weight].>100.0,:]
    @time refinePSMs!(PSMs, precursors_mouse_detailed_33NCEcorrected_start1)
    features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn]
    PSMs[isnan.(PSMs[:,:matched_ratio]),:matched_ratio] .= Inf
    PSMs[(PSMs[:,:matched_ratio]).==Inf,:matched_ratio] .= maximum(PSMs[(PSMs[:,:matched_ratio]).!=Inf,:matched_ratio])
    replace!(PSMs[:,:city_block], -Inf => minimum(PSMs[PSMs[:,:city_block].!=-Inf,:city_block]))
    replace!(PSMs[:,:scribe_score], Inf => minimum(PSMs[PSMs[:,:scribe_score].!=Inf,:scribe_score]))
    #PSMs = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/PSMs_080423.csv"))
    transform!(PSMs, AsTable(:) => ByRow(psm -> length(collect(eachmatch(r"ox", psm[:sequence])))) => [:Mox]);
    features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn]
end
features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all,:kendall,:spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn]


features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all,:kendall,:spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn]



@time rankPSMs!(PSMs, features, colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 0.25, n_folds = 2, num_round = 200, eta = 0.0375, max_depth = 5)
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:]

PSMs_second[(PSMs_second[:,:q_value].<=0.01).&(PSMs_second[:,:decoy].==false),:]


Plots.histogram(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)
Plots.histogram!(PSMs[(PSMs[:,:decoy].==true),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)

@time getQvalues!(PSMs, Float64.(PSMs[:,:kendall]), PSMs[:,:decoy]);
size(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:]) #205589


@time getQvalues!(PSMs, -1*Float64.(PSMs[:,:city_block]), PSMs[:,:decoy]);
size(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:]) #108781

@time getQvalues!(PSMs, Float64.(PSMs[:,:scribe_score]), PSMs[:,:decoy]);
size(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:]) #130952

sum((PSMs[:,:city_block].<-5).&(PSMs[:,:decoy]))/sum((PSMs[:,:city_block].<-5).&(PSMs[:,:decoy].==false))

Plots.histogram(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:city_block], normalize = :pdf, bins = 100, alpha = 0.5)
Plots.histogram!(PSMs[(PSMs[:,:decoy].==true),:city_block], normalize = :pdf, bins = 100, alpha = 0.5)


Plots.histogram(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:city_block], normalize = :pdf, bins = 100, alpha = 0.5)
Plots.histogram!(PSMs[(PSMs[:,:decoy].==true),:city_block], normalize = :pdf, bins = 100, alpha = 0.5)


targets = (PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false)
p = Plots.plot(PSMs[targets,:RT], PSMs[targets,:iRT], seriestype=:scatter,
            xlabel = "Retention Time RT (min)",
            ylabel ="Indexed Retention Time iRT (min)",
            label = nothing,
            size = 100*[13.3, 7.5]
,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm, 
alpha = 0.01)
#Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/bigplot.png")

#Plots.plot(best_psms[best_psms[:,:q_value].<=0.01,:RT], best_psms[best_psms[:,:q_value].<=0.01,:iRT], seriestype=:scatter)
rt_map_2 = KDEmapping(PSMs[targets,:RT], PSMs[targets,:iRT], n = 50)
Plots.plot!(p, LinRange(minimum(PSMs[targets,:RT]), maximum(PSMs[targets,:RT]), 100), 
            rt_map_2.(LinRange(minimum(PSMs[targets,:RT]), maximum(PSMs[targets,:RT]), 100)),
            lw = 6.0,
            label = "RT Spline Final")
 
Plots.plot!(p, LinRange(minimum(rtPSMs[:,:RT]), maximum(rtPSMs[:,:RT]), 100), 
            rt_map.(LinRange(minimum(rtPSMs[:,:RT]), maximum(rtPSMs[:,:RT]), 100)),
            lw = 6.0,
            label = "RT Spline Initial")

Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/bigplot.png")


t = 0.99
p = Plots.plot(title = "Target vs. Decoy PSMs", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)
r = ROCAnalysis.roc(PSMs[PSMs[:,:decoy],:prob], PSMs[PSMs[:,:decoy].==false,:prob]) 

Plots.plot!(p, r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], ylabel = "Precision", xlabel = "Recall", legend = :outertopright, label = "XGBoost", lw = 6)

r = ROCAnalysis.roc(PSMs[PSMs[:,:decoy],:kendall], PSMs[PSMs[:,:decoy].==false,:kendall]) 
Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], label = "Entropy", lw = 6)

r = ROCAnalysis.roc(PSMs[PSMs[:,:decoy],:scribe_score], PSMs[PSMs[:,:decoy].==false,:scribe_score]) 
Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], label = "Scribe", lw = 6)

r = ROCAnalysis.roc(PSMs[PSMs[:,:decoy],:spectral_contrast_all], PSMs[PSMs[:,:decoy].==false,:spectral_contrast_all]) 
Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], label = "Cosine", lw = 6)

r = ROCAnalysis.roc(PSMs[PSMs[:,:decoy],:hyperscore], PSMs[PSMs[:,:decoy].==false,:hyperscore]) 
Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], label = "XTandem!", lw = 6)


r = ROCAnalysis.roc(PSMs[PSMs[:,:decoy],:total_ions], PSMs[PSMs[:,:decoy].==false,:total_ions]) 
Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], label = "Total Ions", lw = 6)

r = ROCAnalysis.roc(PSMs[PSMs[:,:decoy],:matched_ratio], PSMs[PSMs[:,:decoy].==false,:matched_ratio]) 
Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], label = "Matched Ratio", lw = 6)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/PSMs_PR_Curve.pdf")


function signifChop(num, digits)
    if num == 0.0 then
        return num
    else
        e = ceil(log10(abs(num)))
        scale = 10^(digits - e)
        return trunc(num * scale) / scale
    end
end

#PSMs_second = PSMs
#PSMs_first = PSMs

#########
#save psms
#########
#CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv", PSMs)
#Get best scoring psm for each precursor at 10% FDR
###########
#precursors_mouse_detailed_33NCEcorrected_chronologe
@time begin
    best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[PSMs[:,:q_value].<=0.1,:], :precursor_idx))
    transform!(best_psms, AsTable(:) => ByRow(psm -> precursors_mouse_detailed_33NCEcorrected_start1[psm[:precursor_idx]].mz) => :prec_mz)
    sort!(best_psms,:RT, rev = false)
    rt_index = buildRTIndex(best_psms)
    size(best_psms)
    #λ, γ = optimizePenalty(λs, γs) #Get optinal penalty
    using DataStructures
    include("src/Routines/LibrarySearch/integratePrecursors.jl")
    @time chroms = integrateRAW(MS_TABLE, rt_index, frags_mouse_detailed_33NCEcorrected_start1, 
                    one(UInt32), 
                    fragment_tolerance=fragment_tolerance, 
                    frag_ppm_err = 3.34930002879957,
                    λ=zero(Float32), #λs[11], 
                    γ=zero(Float32), #γs[11], 
                    max_peak_width = 2.0, 
                    scan_range = (0, 300000)
                    );
    transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(chroms, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity, :count, :SN, :slope, :peak_error,:apex,:fwhm]);

    sum((best_psms[:,:intensity].>0).&(best_psms[:,:count].>=6))
    best_psms = best_psms[(best_psms[:,:intensity].>0).&(best_psms[:,:count].>=6),:];
    best_psms[:,:RT_error] = abs.(best_psms[:,:apex] .- best_psms[:,:RT_pred])

    function selectIsotopes(prec_list::Vector{Tuple{Float64, UInt32}}, isotope_dict::UnorderedDictionary{UInt32, Vector{Isotope{U}}}, rt::T, rt_tol::T) where {T,U<:AbstractFloat}
        isotopes = Vector{Isotope{U}}()
        i = 1
        rt_start = searchsortedfirst(prec_list, rt - rt_tol, lt=(r,x)->first(r)<x) #First RT bin to search
        rt_stop = searchsortedlast(prec_list, rt + rt_tol, lt=(x, r)->first(r)>x) #Last RT bin to search 
        #return rt_start, rt_stop
        for i in range(rt_start, rt_stop)
            append!(isotopes, isotope_dict[last(prec_list[i])])
        end
        return sort(isotopes, by = x->getMZ(x))
    end

    function integrateRAW(
                        spectra::Arrow.Table, 
                        #rt_index::retentionTimeIndex{T, U},
                        prec_list::Vector{Tuple{Float64, UInt32}},
                        isotopes::UnorderedDictionary{UInt32, Vector{Isotope{Float32}}},
                        ms_file_idx::UInt32;
                        precursor_tolerance::Float64 = 20.0,
                        quadrupole_isolation_width::Float64 = 8.5,
                        max_peak_width::Float64 = 2.0,
                        λ::Float32 = Float32(2e12),
                        γ::Float32 = Float32(1/2),
                        max_iter::Int = 1000,
                        nmf_tol::Float32 = Float32(100.0),
                        scan_range::Tuple{Int64, Int64} = (0, 0), 
                        ) where {T,U<:Real}
        
        ms1 = 0
        nmf = Dict(:precursor_idx => UInt32[], :weight => Float32[], :rt => Float32[])
        matches = ""
        misses = ""
        for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))

            if spectrum[:msOrder] == 2
                continue
            else
                ms1 += 1
            end
            if scan_range != (0, 0)
                i < first(scan_range) ? continue : nothing
                i > last(scan_range) ? continue : nothing
            end
            #Get peptides that could be in the spectra
            #transitions = selectTransitions(fragment_list, rt_index, Float64(spectrum[:retentionTime]), max_peak_width/2.0, spectrum[:precursorMZ], Float32(quadrupole_isolation_width/2.0))
            #isotopes[0x0006bbe9]
            #Match fragments to peaks
            iso = selectIsotopes(prec_rt_table, isotopes, Float64(spectrum[:retentionTime]), 1.0)
            matches, misses = matchPeaks(iso,
                                spectrum[:masses],
                                spectrum[:intensities],
                                PrecursorMatch{Float32},
                                count_unmatched=true,
                                δs = zeros(Float64, (1,)),
                                ppm = precursor_tolerance
                                )

            #=fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                        spectrum[:masses], 
                                        spectrum[:intensities], 
                                        count_unmatched =true,
                                        δs = zeros(T, (1,)),
                                        scan_idx = UInt32(i),
                                        ms_file_idx = ms_file_idx,
                                        min_intensity = zero(Float32),
                                        ppm = fragment_tolerance
                                        )=#


            if iszero(length(matches))
                continue
            end

            #Build templates for regrssion. 
            #Do we need to remove precursors with less than N matched fragments?
            #if spectrum[:retentionTime] < 49.8483 
            #    continue
            #elseif spectrum[:retentionTime] > 49.8992 
            #end
            X, Hs, Hst, IDtoROW = buildDesignMatrix(matches, misses)
            #return X, Hs, Hst, IDtoROW
            weights = sparseNMF(Hst, Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)

            for key in keys(IDtoROW)
                push!(nmf[:precursor_idx], key)
                push!(nmf[:weight], weights[IDtoROW[key]])
                push!(nmf[:rt], spectrum[:retentionTime])
            end
        end
        nmf = DataFrame(nmf)
        sort!(nmf, [:precursor_idx,:rt]);
        return groupby(nmf, :precursor_idx)
        #matches, misses
    end

    @time isotopes =  getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], best_psms[:,:charge], QRoots(4), 4)
    prec_rt_table = sort(collect(zip(best_psms[:,:RT], UInt32.(best_psms[:,:precursor_idx]))), by = x->first(x))
    test_df = integrateRAW(MS_TABLE, prec_rt_table, isotopes, one(UInt32), precursor_tolerance = 6.5, scan_range = (0, 300000), λ = Float32(0), γ = Float32(0))


    transform!(best_psms, AsTable(:) => ByRow(psm -> getCrossCorr(test_df, chroms, UInt32(psm[:precursor_idx]))) => [:offset,:cross_cor]);
    transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(test_df, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity_ms1, :count_ms1, :SN_ms1, :slope_ms1, :peak_error_ms1,:apex_ms1,:fwhm_ms1]);

    features = [:hyperscore, :total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:kendall, :missed_cleavage,:Mox,:best_rank,:topn]

        
    features = [:hyperscore, :total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all,:RT_error,:y_ladder,:b_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:kendall, :scribe_score, :missed_cleavage,:Mox,:best_rank,:topn]

    append!(features, [:intensity, :count, :SN, :peak_error,:fwhm,:offset, :cross_cor])
    #append!(features, [:intensity, :offset, :cross_cor])
    #append!(features, [:intensity, :count, :SN, :peak_error,:fwhm])
    @time bst = rankPSMs!(best_psms, features,colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 0.5, n_folds = 5, num_round = 200, max_depth = 10, eta = 0.0375)
    #@time bst = rankPSMs!(best_psms, features,colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 0.5, n_folds = 5, num_round = 300, eta = 0.025)
    @time getQvalues!(best_psms, best_psms[:,:prob], best_psms[:,:decoy]);
    length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:precursor_idx])) #42348, 43224
end

p = Plots.plot(title = "Target vs. Decoy Precursors", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)

t = 0.99
r = ROCAnalysis.roc(best_psms[best_psms[:,:decoy],:prob], best_psms[best_psms[:,:decoy].==false,:prob]) 
Plots.plot!(p, r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], ylabel = "Precision", xlabel = "Recall", legend = :outertopright, label = "XGBoost", lw = 6)

r = ROCAnalysis.roc(best_psms[best_psms[:,:decoy],:kendall], best_psms[best_psms[:,:decoy].==false,:kendall]) 
p = Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], legend = :outertopright, label = "Entropy", lw = 6)

r = ROCAnalysis.roc(best_psms[best_psms[:,:decoy],:spectral_contrast_all], best_psms[best_psms[:,:decoy].==false,:spectral_contrast_all]) 
p = Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], legend = :outertopright, label = "Cosine", lw = 6)

r = ROCAnalysis.roc(best_psms[best_psms[:,:decoy],:hyperscore], best_psms[best_psms[:,:decoy].==false,:hyperscore]) 
p = Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], legend = :outertopright, label = "XTandem!", lw = 6)

r = ROCAnalysis.roc(best_psms[best_psms[:,:decoy],:matched_ratio], best_psms[best_psms[:,:decoy].==false,:matched_ratio]) 
p = Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], legend = :outertopright, label = "Matched Ratio", lw = 6)

r = ROCAnalysis.roc(best_psms[best_psms[:,:decoy],:scribe_score], best_psms[best_psms[:,:decoy].==false,:scribe_score]) 
p = Plots.plot!(r.pfa[r.pmiss.>=t], r.pmiss[r.pmiss.>=t], legend = :outertopright, label = "Scribe", lw = 6)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/Precursors_PR_Curve.pdf")


Plots.histogram(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)
Plots.histogram!(PSMs[(PSMs[:,:decoy].==true),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)


Plots.histogram(PSMs[(PSMs[:,:decoy].==false),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)
Plots.histogram!(PSMs[(PSMs[:,:decoy].==true),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)


p = Plots.plot(title = "Target vs. Decoy PSMs", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)
Plots.histogram!(p, (PSMs[(PSMs[:,:decoy].==false),:prob]), normalize = :pdf, bins = 100, alpha = 0.5,
label = "Targets")
Plots.histogram!((PSMs[(PSMs[:,:decoy].==true),:prob]), normalize = :pdf, bins = 100, alpha = 0.5,
label = "Decoys")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_prob.pdf")


p = Plots.plot(title = "MS1 vs MS2 Aligment", 
xlabel = "Offset", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)

best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:offset]
Plots.histogram!(p, (best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:offset]), normalize = :pdf, bins = 100, alpha = 0.5,
label = "Targets")
Plots.histogram!((best_psms[(best_psms[:,:decoy].==true),:offset]), normalize = :pdf, bins = 100, alpha = 0.5,
label = "Decoys")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_offset.pdf")

p = Plots.plot(title = "MS1 vs MS2 Aligment", 
xlabel = "Offset", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)

best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:offset]
Plots.plot!(p, log10.(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:intensity]), log10.(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:intensity_ms1]), normalize = :pdf, bins = 100, alpha = 0.5,seriestype=:scatter,ylim = (0, 10),
label = "Targets")
Plots.plot!(log10.(best_psms[(best_psms[:,:decoy].==true),:intensity]), log10.(best_psms[(best_psms[:,:decoy].==true),:intensity_ms1]), normalize = :pdf, bins = 100, alpha = 0.5,seriestype=:scatter,
label = "Decoys", ylim = (0, 10))
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_ms1_ms2.png")

p = Plots.plot(title = "Log10 Intensities", 
xlabel = "Log10 Intensity", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)

best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:offset]
Plots.histogram!(p, log10.(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:intensity]), normalize = :pdf, bins = 100, alpha = 0.5,seriestype=:scatter,
label = "Targets")
Plots.histogram!(log10.(best_psms[(best_psms[:,:decoy].==true),:intensity]), normalize = :pdf, bins = 100, alpha = 0.5,seriestype=:scatter,
label = "Decoys")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_logint.png")





p = Plots.plot(title = "Target vs. Decoy PSMs", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)
Plots.histogram!(p, min.(log2.(PSMs[(PSMs[:,:decoy].==false),:matched_ratio]), 5), normalize = :pdf, bins = 100, alpha = 0.5,
label = "Targets", xlabel = "Matched Ratio", ylabel = "PDF")
Plots.histogram!(p, min.(log2.(PSMs[(PSMs[:,:decoy].==true),:matched_ratio]), 5), normalize = :pdf, bins = 100, alpha = 0.5,
label = "Decoys")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_matched_ratio.pdf")

p = Plots.plot(title = "Target vs. Decoy PSMs", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)
Plots.histogram!(p, min.(log2.(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:matched_ratio]), 5), 
normalize = :pdf, bins = 100, alpha = 0.5,
label = "Targets <= 1% FDR", xlabel = "Matched Ratio", ylabel = "PDF")
Plots.histogram!(p, min.(log2.(PSMs[(PSMs[:,:decoy].==true),:matched_ratio]), 5), normalize = :pdf, bins = 100, alpha = 0.5,
label = "Decoys")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_matched_ratio_1fdr.pdf")



p = Plots.plot(title = "Target vs. Decoy PSMs", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)
bins =  LinRange(0, 25, 100)
Plots.histogram!(p, PSMs[(PSMs[:,:decoy].==false),:scribe_score], normalize = :pdf, alpha = 0.5,
label = "Targets", xlabel = "Scribe Score", ylabel = "PDF",  xlim = (0, 25), bins = bins)
Plots.histogram!(p, PSMs[(PSMs[:,:decoy].==true),:scribe_score], normalize = :pdf, alpha = 0.5,
label = "Decoys", xlim = (0, 25), bins = bins)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_scribe.pdf")

p = Plots.plot(title = "Target vs. Decoy PSMs", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)
bins =  LinRange(0, 1, 100)
Plots.histogram!(p, PSMs[(PSMs[:,:decoy].==false),:kendall], normalize = :pdf, alpha = 0.5,
label = "Targets", xlabel = "Entropy Simmilarity", ylabel = "PDF",  xlim = (0, 1), bins = bins)
Plots.histogram!(p, PSMs[(PSMs[:,:decoy].==true),:kendall], normalize = :pdf, alpha = 0.5,
label = "Decoys", xlim = (0, 1), bins = bins)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_entropy.pdf")

p = Plots.plot(title = "Target vs. Decoy PSMs", 
xlabel = "XGBoost Prob", 
ylabel = "PDF", 
size = 100*[13.3, 7.5]
, thicknesss_scaling = 10,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm
)
bins =  LinRange(0, 1, 100)
Plots.histogram!(p, PSMs[(PSMs[:,:decoy].==false).&(PSMs[:,:q_value].<=0.01),:spectral_contrast_all], normalize = :pdf, alpha = 0.5,
label = "Targets <= 0.1% FDR", xlabel = "Spectral Contrast", ylabel = "PDF",  xlim = (0, 1), bins = bins)
Plots.histogram!(p, PSMs[(PSMs[:,:decoy].==true),:spectral_contrast_all], normalize = :pdf, alpha = 0.5,
label = "Decoys", xlim = (0, 1), bins = bins)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/target_decoy_hist_cosine.pdf")




Plots.histogram(best_psms[(best_psms[:,:decoy].==false),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)
Plots.histogram!(best_psms[(best_psms[:,:decoy].==true),:kendall], normalize = :pdf, bins = 100, alpha = 0.5)



#features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight]
#append!(features, [:intensity, :count, :SN, :slope, :peak_error,:apex,:fwhm])
#features = [:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all,:RT_error,:scribe_score,:RT,:diff_hyper,:Mox,:city_block,:matched_ratio,:weight,:intensity,:intensity_ms1,:peak_error,:fwhm,:count,:offset,:cross_cor,:missed_cleavage]
#@time bst = rankPSMs!(best_psms, features,colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 1.0, n_folds = 5, num_round = 200, eta = 0.0375)
#@time getQvalues!(best_psms, best_psms[:,:prob], best_psms[:,:decoy]);
#length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:precursor_idx])) #36782

transform!(best_psms, AsTable(:) => ByRow(psm -> test_table.id_to_pepGroup[test_table.id_to_pep[psm[:pep_id]].pepGroup_id].prot_ids) => :prot_ids)
#transform!(best_psms, AsTable(:) => ByRow(psm -> precursors_mouse_detailed_33NCEdynamic[psm[:precursor_idx]].p#ep_id) => :pep_id)

transform!(best_psms, AsTable(:) => ByRow(psm -> test_table.id_to_pepGroup[psm[:pep_id]]) => :prot_ids)
length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:sequence]))
prot_groups = best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:prot_ids]


best_psms[:,[:sequence,:q_value,:intensity,:count,:SN,:peak_error,:fwhm,:precursor_idx,:scan_idx,:prot_ids]]
best_psms[(best_psms[:,:intensity].>0.0) .& (best_psms[:,:q_value].<=0.01),[:precursor_idx,:q_value,:decoy,:intensity,:peak_error,:fwhm,:SN,:RT,:count,:slope,:apex]][30000:40000,:]
pid = UInt32(4294810)
integratePrecursor(chroms, pid, isplot = true)
chroms[(precursor_idx=pid,)]
pid = UInt32(  3470831 )
integratePrecursor(chroms, pid, isplot = true)
chroms[(precursor_idx=pid,)]
pid = UInt32(3531828) ####
integratePrecursor(chroms, pid, isplot = true)
chroms[(precursor_idx=pid,)]
pid = UInt32( 935412 )
integratePrecursor(chroms, pid, isplot = true)
pid = UInt32(1780792) #Check again
integratePrecursor(chroms, pid, isplot = true)
plot!(chroms[(precursor_idx=pid,)][:,:rt], chroms[(precursor_idx=pid,)][:,:weight], seriestype=:scatter)
chroms[(precursor_idx=pid,)]
pid = UInt32(445492) #Check again
integratePrecursor(chroms, pid, isplot = true)
pid = UInt32( 5385393) #Check again
integratePrecursor(chroms, pid, isplot = true)
pid = UInt32(4585201  ) #Check again
integratePrecursor(chroms, pid, isplot = true)

pid = UInt32(1370088) #Check again
integratePrecursor(chroms, pid, isplot = true)
plot(test_df[(precursor_idx=pid,)][:,:rt], test_df[(precursor_idx=pid,)][:,:weight], seriestype=:scatter)
plot!(chroms[(precursor_idx=pid,)][:,:rt], chroms[(precursor_idx=pid,)][:,:weight]*10, seriestype=:scatter)




chroms[(precursor_idx=pid,)]
integratePrecursor(test_df, pid, isplot = true)
#=
@time getQvalues!(non_zero, non_zero[:,:prob], non_zero[:,:decoy]);
integratePrecursor(chroms, UInt32(4171604), isplot = true)
savefig("/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_a.pdf")
integratePrecursor(chroms, UInt32(   3574807 ), isplot = true)
savefig("/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_b.pdf")
integratePrecursor(chroms, UInt32(  2642152), isplot = true)
savefig("/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_c.pdf")
integratePrecursor(chroms, UInt32(   508178 ), isplot = true)
savefig("/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_d.pdf")
merge_pdfs(["/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_a.pdf",
"/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_b.pdf",
"/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_c.pdf",
"/Users/n.t.wamsley/Projects/TEST_DATA/chromatogram_d.pdf"], "/Users/n.t.wamsley/Projects/TEST_DATA/chromatograms.pdf")
CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/best_psms_072223_03.csv", best_psms)

p = plot(title="Log10 of Peak Areas for Quantified Precursors");
histogram(p, log10.(non_zero[:,:intensity]))
savefig("/Users/n.t.wamsley/Projects/TEST_DATA/peak_area_hist.pdf")

frags = Int64[]
rt = Float32[]
for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(MS_TABLE)))
    if spectrum[:msOrder] == 2
        push!(frags, length(spectrum[:masses]))
        push!(rt, spectrum[:retentionTime])
    end
end

function getSeqs()
    rows = CSV.Rows("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv", reusebuffer=false, stringtype=String, select = [:ModifiedPeptide])
        target_seqs = Vector{String}()
        for (i, row) in ProgressBar(enumerate(rows))
            seq = row.ModifiedPeptide
            if i < 2
                push!(target_seqs, seq)
            end
            if target_seqs[end] != seq
                push!(target_seqs, seq)
            end 
        end
    return target_seqs
end

@time target_seqs = CSV.File("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv", header=1, stringtype=String, select=[:ModifiedPeptide], ntasks = 24)
target_seqs = [x[1] for x in target_seqs]
target_seqs = Set(target_seqs)

@time decoy_seqs = CSV.File("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_decoys.csv", header=1, stringtype=String, select=[:ModifiedPeptide], ntasks = 24)
@time decoy_seqs  = [x[1] for x in  decoy_seqs ]
@time decoy_seqs = Set( decoy_seqs )
=#

test = integrateRAW(MS_TABLE, rt_index, frags_mouse_detailed_33NCEcorrected_start1, 
one(UInt32), 
fragment_tolerance=15.6, 
λ=λs[11], 
γ=γs[11], 
max_peak_width = 2.0, 
scan_range = (0, 300000)
);

testdict = Accumulator{UInt32, Int64}()
a = counter(UInt32)
for match in fragmentMatches
    DataStructures.inc!(a, match.prec_id)
end
filter(x->(a[x.prec_id].>=4), fragmentMatches)
filter(x->(a[x.prec_id].>=4), fragmentMisses)
keys(filter(((k,v),) -> v >= 4, a))

countmap([x.prec_id for x in fragmentMatches])
filter(((k,v),) -> k == 1 || v == "c", dict)

times = DataFrame(Dict(
    :σ => [0, 1, 1.5, 2, 2.5, 3, 4, 4.28, 5, 6],
    :PSMs => [0, 219235, 289468, 330926,357066,371218,379823,380115,376461,367167],
    :Precursors => [0, 31529, 37837, 41877, 43939, 45580, 46498, 46776, 47047, 46934],
    :Seed => [0, 284489, 479391, 705481, 978078, 1311887, 2128633, 2404872, 3136221,4246771],
    :main_time => [0, 161, 198, 238, 279, 316, 413, 463, 528, 739]
))

p = Plots.plot()
p_twin = twinx(p)
Plots.plot!(p, times[:,:σ], times[:,:Seed], label="Seed Size", legend=:bottomright, seriestype=:scatter, ylabel="Seed Size")
Plots.plot!(p, times[:,:σ], times[:,:Seed], label=nothing, color = palette(:tab10)[1])
Plots.plot!(p_twin,times[:,:σ], times[:,:Precursors],xticks=:none,legend=:topleft, label=nothing, seriestype=:scatter, color = palette(:tab10)[2],
ylabel = "Precursors", bottom_margin = 10*Plots.mm, margin = 100cm)
Plots.plot!(p_twin,times[:,:σ], times[:,:Precursors],label = nothing, xticks=:none, color = palette(:tab10)[2])
annotate!(3.5, -0.5e6, "σ Mass Tolerance")
Plots.plot!([], [], label = "Precursors at 1% FDR", seriestype=:scatter, color= palette(:tab10)[2])
Plots.vline!([3], label = "3σ", lw = 3, color = :firebrick)
#lots.vline!([15.6/params(mix_mle)[1][1][2]], label = "DIA-NN")
Plots.vline!([4.287741187029869], label = "DIA-NN", lw = 3, color = :darkblue)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/seed_size_precursors.pdf")

p = Plots.plot()
p_twin = twinx(p)
Plots.plot!(p, times[:,:σ], times[:,:Seed], label="Seed Size", legend=:bottomright, seriestype=:scatter, ylabel="Seed Size")
Plots.plot!(p, times[:,:σ], times[:,:Seed], label=nothing, color = palette(:tab10)[1])
Plots.plot!(p_twin,times[:,:σ], times[:,:PSMs],xticks=:none,legend=:topleft, label=nothing, seriestype=:scatter, color = palette(:tab10)[2],
ylabel = "PSMs", bottom_margin = 10*Plots.mm, margin = 100cm)
Plots.plot!(p_twin,times[:,:σ], times[:,:PSMs],label = nothing, xticks=:none, color = palette(:tab10)[2])
annotate!(3.5, -0.5e6, "σ Mass Tolerance")
Plots.plot!([], [], label = "PSMs at 1% FDR", seriestype=:scatter, color= palette(:tab10)[2])
Plots.vline!([3], label = "3σ", lw = 3, color = :firebrick)
#lots.vline!([15.6/params(mix_mle)[1][1][2]], label = "DIA-NN")
Plots.vline!([4.287741187029869], label = "DIA-NN", lw = 3, color = :darkblue)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/seed_size_psms.pdf")


p = Plots.plot()
p_twin = twinx(p)
Plots.plot!(p, times[:,:σ], times[:,:main_time], label="Main Search Time (s)", legend=:bottomright, seriestype=:scatter, ylabel="Main Search Time (s)")
Plots.plot!(p, times[:,:σ], times[:,:main_time], label=nothing, color = palette(:tab10)[1])
#annotate!(3.5, -100.0, "σ Mass Tolerance")
Plots.plot!(p_twin,times[:,:σ], times[:,:PSMs],xticks=:none,legend=:topleft, label=nothing, seriestype=:scatter, color = palette(:tab10)[2],
ylabel = "PSMs", bottom_margin = 10*Plots.mm, margin = 100cm)

Plots.plot!(p_twin,times[:,:σ], times[:,:PSMs],label = nothing, xticks=:none, color = palette(:tab10)[2])


Plots.plot!([], [], label = "PSMs at 1% FDR", seriestype=:scatter, color= palette(:tab10)[2])
Plots.vline!([3], label = "3σ", lw = 3, color = :firebrick)
#lots.vline!([15.6/params(mix_mle)[1][1][2]], label = "DIA-NN")
Plots.vline!([4.287741187029869], label = "DIA-NN", lw = 3, color = :darkblue)
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/time_psms.pdf")