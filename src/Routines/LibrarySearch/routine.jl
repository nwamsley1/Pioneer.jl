
##########
#Import Libraries
##########
using CSV, Arrow, Tables, DataFrames, Dictionaries, Combinatorics, StatsBase, NMF, JLD2, LinearAlgebra, Random, DecisionTree, LoopVectorization, Splines2, ProgressBars, GLM, RobustModels, LoopVectorization, SparseArrays, Interpolations, MultiKDE, XGBoost, SavitskyGolay, NumericalIntegration
##########
#Import files
##########
include("src/precursor.jl")
include("src/Routines/LibrarySearch/parsePrositLib.jl")
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
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/precursor_list.jld2" precursor_list
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/linear_spline.jld2" linear_spline

MS_TABLE = Arrow.Table("/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA.arrow")

###########
#Main Search
###########

PSMs = SearchRAW(MS_TABLE, prosit_index_5ppm_15irt, frag_detailed, UInt32(1), linear_spline,
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
@time refinePSMs!(PSMs, precursor_list)
features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight]

PSMs[isnan.(PSMs[:,:matched_ratio]),:matched_ratio] .= Inf
PSMs[(PSMs[:,:matched_ratio]).==Inf,:matched_ratio] .= maximum(PSMs[(PSMs[:,:matched_ratio]).!=Inf,:matched_ratio])
#replace!(PSMs[:,:city_block], -Inf => minimum(PSMs[PSMs[:,:city_block].!=-Inf,:city_block]))
#replace!(PSMs[:,:scribe_score], Inf => minimum(PSMs[PSMs[:,:scribe_score].!=Inf,:scribe_score]))
@time rankPSMs!(PSMs, features, colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 1.0, n_folds = 2, num_round = 200, eta = 0.0375)
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
#########
#save psms
#########
#CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv", PSMs)
###########
#Get best scoring psm for each precursor at 10% FDR
###########
best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[PSMs[:,:q_value].<=0.1,:], :precursor_idx))
transform!(best_psms, AsTable(:) => ByRow(psm -> precursor_list[psm[:precursor_idx]].mz) => :prec_mz)
sort!(best_psms,:RT, rev = false)
rt_index = buildRTIndex(best_psms)
##########
#Integrate precursors
##########
chroms = integrateRAW(MS_TABLE, rt_index, frag_detailed, 
                        one(UInt32), 
                        fragment_tolerance=15.6, 
                        λ=Float32(5e3), 
                        γ=Float32(0), 
                        max_peak_width = 2.0, 
                        scan_range = (0, 300000)
                        );

transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(chroms, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity, :count, :SN, :slope, :peak_error,:apex]);
best_psms = best_psms[(best_psms[:,:intensity].>0).&(best_psms[:,:count].>=5),:];
append!(features, [:intensity, :count, :SN, :slope, :peak_error,:apex])
@time bst = rankPSMs!(best_psms, features,colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 1.0, n_folds = 5, num_round = 200, eta = 0.0375)
@time getQvalues!(best_psms, best_psms[:,:prob], best_psms[:,:decoy]);
length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:precursor_idx]))



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