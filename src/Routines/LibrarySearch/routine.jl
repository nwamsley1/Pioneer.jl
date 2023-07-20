
##########
#Import Libraries
##########
using CSV, Arrow, Tables, DataFrames, Dictionaries, Combinatorics, StatsBase, NMF, JLD2, LinearAlgebra, Random, DecisionTree, LoopVectorization, Splines2, ProgressBars, GLM, RobustModels, LoopVectorization
##########
#Import files
##########
include("src/precursor.jl")
include("src/Routines/LibrarySearch/buildFragmentIndex.jl")
include("src/Routines/LibrarySearch/matchpeaksLIB.jl")
include("src/Routines/LibrarySearch/buildDesignMatrix.jl")
include("src/Routines/LibrarySearch/spectralDistanceMetrics.jl")
include("src/Routines/LibrarySearch/spectralContrast.jl")
include("src/Routines/LibrarySearch/searchRAW.jl")
include("src/Routines/LibrarySearch/counter.jl")
include("src/Routines/LibrarySearch/ML.jl")
include("src/Routines/LibrarySearch/refinePSMs.jl")
##########
include("src/Routines/LibrarySearch/queryFragmentArr.jl")
include("src/PSM_TYPES/PSM.jl")
include("src/PSM_TYPES/LibraryXTandem.jl")
include("src/Routines/LibrarySearch/selectTransitions.jl")

@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_detailed.jld2"  prosit_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_index_intensities.jld2"  prosit_index_intensities
@load "/Users/n.t.wamsley/Projects/PROSIT/prosit_precs.jld2"  prosit_precs

MS_TABLE = Arrow.Table("/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA.arrow")


#This needs to be included in prosit_precs
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

refinePSMs!(PSMs, precursor_list)

PSMs[:,:iRT_pred] = [linear_spline(x) for x in PSMs[:,:RT]]
PSMs[:,:RT_error] = abs.(PSMs[:,:iRT] .- PSMs[:,:iRT_pred])
features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight]
@time rankPSMs!(PSMs, features, n_folds = 2, n_trees = 100, max_depth = 15, n_features = 10, fraction = 0.001)
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
#########PS
#save psms
#########
#CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv", PSMs)
###########
#Get best scoring psm for each precursor at 10% FDR
###########
best_psms = combine(sdf -> sdf[argmax(sdf.hyperscore),:], groupby(PSMs[PSMs[:,:q_values].<=0.1,:], :precursor_idx))
transform!(best_psms, AsTable(:) => ByRow(psm -> prec_mzs[psm[:precursor_idx]]) => :prec_mz)
rt_index = buildRTIndex(best_psms)
##########
#Integrate precursors
##########
chroms = integrateRAW(MS_TABLE, rt_index, prosit_detailed, 
                        one(UInt32), 
                        fragment_tolerance=15.6, 
                        λ=Float32(5e3), 
                        γ=Float32(0), 
                        max_peak_width = 2.0, 
                        scan_range = (0, 300000)
                        );

transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(chroms, psm[:precursor_idx], isplot = false)) => [:intensity, :count, :SN]);
non_zero = best_psms[(best_psms[:,:intensity].>0).&(best_psms[:,:count].>5),:];
@time getQvalues!(non_zero, non_zero[:,:prob], non_zero[:,:decoy]);
integratePrecursor(chroms, UInt32(4171604), isplot = true)
