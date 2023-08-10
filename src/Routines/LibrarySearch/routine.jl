
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

newPSMs = SearchRAW(MS_TABLE, prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  frags_mouse_detailed_33NCEcorrected_start1, UInt32(1), linear_spline,
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

PSMs = newPSMs
@time refinePSMs!(PSMs, precursors_mouse_detailed_33NCEcorrected_start1)


#iRT_to_RT = refinePSMs!(PSMs, precursors_mouse_detailed_33NCEcorrected_start1)
#features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn]
features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn]
PSMs[isnan.(PSMs[:,:matched_ratio]),:matched_ratio] .= Inf
PSMs[(PSMs[:,:matched_ratio]).==Inf,:matched_ratio] .= maximum(PSMs[(PSMs[:,:matched_ratio]).!=Inf,:matched_ratio])
replace!(PSMs[:,:city_block], -Inf => minimum(PSMs[PSMs[:,:city_block].!=-Inf,:city_block]))
replace!(PSMs[:,:scribe_score], Inf => minimum(PSMs[PSMs[:,:scribe_score].!=Inf,:scribe_score]))
#PSMs = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/PSMs_080423.csv"))
transform!(PSMs, AsTable(:) => ByRow(psm -> length(collect(eachmatch(r"ox", psm[:sequence])))) => [:Mox]);

@time rankPSMs!(PSMs, features, colsample_bytree = 0.5, min_child_weight = 10, gamma = 10, subsample = 1.0, n_folds = 2, num_round = 200, eta = 0.0375)
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);


#best_psms[:,:old] = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs_old[PSMs_old[:,:q_value].<=0.01,:], :precursor_idx))
#best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[PSMs[:,:q_value].<=0.01,:], :precursor_idx)) #53917
best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[PSMs[:,:q_value].<=0.1,:], :precursor_idx)) #119107

#value_counts(df, col) = combine(groupby(df, col), nrow)
#histogram(value_counts(PSMs, :scan_idx)[:,:nrow])
describe(value_counts(PSMs, :scan_idx)[:,:nrow])
#########
#save psms
#########
#CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv", PSMs)

##########
#Integrate precursors
##########
#[5e2*(1.25^n) for n in range(1, 15)]
λs = Float32[1e3*(1.5^n) for n in range(1, 40)]
#Float32[1e3*(1.5^n) for n in range(1, 30)]
γs = Float32[1/2 for n in range(1, 40)]


λs = Float32[1e1*(1.5^n) for n in range(1, 40)]
#Float32[1e3*(1.5^n) for n in range(1, 30)]
γs = Float32[0 for n in range(1, 40)]


λs = Float32[1e5*(1.5^n) for n in range(1, 40)]
#Float32[1e3*(1.5^n) for n in range(1, 30)]
γs = Float32[1 for n in range(1, 40)]


λs = Float32[1e9*(1.5^n) for n in range(1, 40)]
#Float32[1e3*(1.5^n) for n in range(1, 30)]
γs = Float32[2 for n in range(1, 40)]

λs = Float32[1e12*(1.5^n) for n in range(1, 40)]
#Float32[1e3*(1.5^n) for n in range(1, 30)]
γs = Float32[3 for n in range(1, 40)]

#Float32[1e5*(1.7^n) for n in range(1, 50)]
λs = Float32[1e9*(1.5^n) for n in range(0, 35)]
#Float32[1e3*(1.5^n) for n in range(1, 30)]
γs = Float32[2 for n in range(0, 35)]
function optimizePenalty(λs, γs)
    targets = []
    decoys = []
    targets_SN = []
    decoys_SN = []
    for i in ProgressBar(eachindex(λs))
        chroms = integrateRAW(MS_TABLE, rt_index, frags_mouse_detailed_33NCEcorrected_start1, 
                            one(UInt32), 
                            fragment_tolerance=15.6, 
                            λ=λs[i], 
                            γ=γs[i], 
                            max_peak_width = 2.0, 
                            scan_range = (0, 300000),
                            mz_range = (706.0, 709.0));
        transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(chroms, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity, :count, :SN, :slope, :peak_error,:apex,:fwhm]);
        #push!(targets, sum((best_psms[(best_psms[:,:decoy].==false).&(best_psms[:,:q_value].<=0.01),:intensity].>0.0)))
        #push!(targets_SN, mean((best_psms[(best_psms[:,:decoy].==false).&(best_psms[:,:intensity].>0.0).&(best_psms[:,:q_value].<=0.01),:SN])))
        push!(targets, sum((best_psms[(best_psms[:,:decoy].==false),:intensity].>0.0)))
        push!(targets_SN, mean((best_psms[(best_psms[:,:decoy].==false).&(best_psms[:,:intensity].>0.0),:SN])))
      
        #push!(decoys, sum((best_psms[(best_psms[:,:decoy].==true),:intensity].>0.0)))
        #push!(decoys_SN, mean((best_psms[(best_psms[:,:decoy].==true).&(best_psms[:,:intensity].>0.0),:SN])))
    end
    tsn = targets.*targets_SN
    p = plot(log2.(λs) ,tsn, seriestype =:scatter)
    plot!(p, log2.(λs) ,tsn)
    plot!(p, log2.(λs) ,targets, seriestype =:scatter)
    plot!(p, log2.(λs) ,targets)
    vline!(p, [log2.(λs)[argmax(tsn)]], show = true)
    hline!([maximum(tsn)*0.98])
    println(maximum(tsn))
    λs[argmax(tsn)], γs[argmax(tsn)]
end

@time chroms = integrateRAW(MS_TABLE, rt_index, frags_mouse_detailed_33NCEcorrected_start1, 
                one(UInt32), 
                fragment_tolerance=15.6, 
                λ=zero(Float32), 
                γ=zero(Float32), 
                max_peak_width = 2.0, 
                scan_range = (0, 300000)
                );
#Remove peaks with zero intensity after no penalty
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(chroms, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity, :count, :SN, :slope, :peak_error,:apex,:fwhm]);
best_psms = best_psms[(best_psms[:,:intensity].>0).&(best_psms[:,:count].>=5),:];


###########
#Get best scoring psm for each precursor at 10% FDR
###########
#precursors_mouse_detailed_33NCEcorrected_chronologer
best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[PSMs[:,:q_value].<=0.1,:], :precursor_idx))
transform!(best_psms, AsTable(:) => ByRow(psm -> precursors_mouse_detailed_33NCEcorrected_start1[psm[:precursor_idx]].mz) => :prec_mz)
sort!(best_psms,:RT, rev = false)
rt_index = buildRTIndex(best_psms)
size(best_psms)
λ, γ = optimizePenalty(λs, γs) #Get optinal penalty
using DataStructures
@time chroms = integrateRAW(MS_TABLE, rt_index, frags_mouse_detailed_33NCEcorrected_start1, 
                one(UInt32), 
                fragment_tolerance=15.6, 
                λ=zero(Float32), #λs[11], 
                γ=zero(Float32), #γs[11], 
                max_peak_width = 2.0, 
                scan_range = (0, 300000)
                );
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(chroms, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity, :count, :SN, :slope, :peak_error,:apex,:fwhm]);

sum((best_psms[:,:intensity].>0).&(best_psms[:,:count].>=5))
best_psms = best_psms[(best_psms[:,:intensity].>0).&(best_psms[:,:count].>=5),:];
best_psms[:,:RT_error] = abs.(best_psms[:,:apex] .- best_psms[:,:RT_pred])
transform!(best_psms, AsTable(:) => ByRow(psm -> getCrossCorr(test_df, chroms, UInt32(psm[:precursor_idx]))) => [:offset,:cross_cor]);
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(test_df, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity_ms1, :count_ms1, :SN_ms1, :slope_ms1, :peak_error_ms1,:apex_ms1,:fwhm_ms1]);
features = [:hyperscore, :total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:median_ions,:n_obs,:charge,:city_block,:matched_ratio, :missed_cleavage,:Mox,:best_rank,:topn]
append!(features, [:intensity, :count, :SN, :peak_error,:fwhm,:offset, :cross_cor])
#append!(features, [:intensity, :offset, :cross_cor])
#append!(features, [:intensity, :count, :SN, :peak_error,:fwhm])
@time bst = rankPSMs!(best_psms, features,colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 0.5, n_folds = 5, num_round = 200, eta = 0.0375)
#@time bst = rankPSMs!(best_psms, features,colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 0.5, n_folds = 5, num_round = 300, eta = 0.025)
@time getQvalues!(best_psms, best_psms[:,:prob], best_psms[:,:decoy]);
length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:precursor_idx])) #42348, 43224

histogram(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:count])
histogram(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:peak_error])
histogram(log2.(abs.(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:fwhm])))
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
