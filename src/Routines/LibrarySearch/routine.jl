
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
#Pre-Search
#Need to Estimate the following from a random sample of high-confidence targets
#1) Fragment Mass error/correction
#2) Fragment Mass tolerance
#3) iRT to RT conversion spline
###########
init_frag_tol = 30.0 #Initial tolerance should probably be pre-determined for each different instrument and resolution. 
rtPSMs, all_matches = SearchRAW(MS_TABLE, 
                    prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
                    frags_mouse_detailed_33NCEcorrected_start1, 
                    UInt32(1), 
                    x->x, #Mapp RT to iRT
                    min_frag_count = 7, 
                    topN = 5, 
                    fragment_tolerance = init_frag_tol,#20.0, 
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

function getPPM(a::T, b::T) where {T<:AbstractFloat}
    (a-b)/(a/1e6)
end
rtPSMs = rtPSMs[rtPSMs[:,:decoy].==false,:]
best_precursors = Set(rtPSMs[:,:precursor_idx])
best_matches = [x for x in all_matches if x.prec_id ∈ best_precursors]

#Model fragment errors with a mixture model of a uniform and laplace distribution 
@time frag_err_dist = estimateErrorDistribution(frag_ppm_errs, Laplace{Float64}, 0.0, 3.0, 30.0)

transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(precursors_mouse_detailed_33NCEcorrected_start1[psm[:precursor_idx]]))) => :iRT)
transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
transform!(rtPSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors_mouse_detailed_33NCEcorrected_start1[psm[:precursor_idx]])) => :decoy)

RT_to_iRT_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT], n = 50)

#fragment_tolerance = 15.6
@time begin
    newPSMs = SearchRAW(MS_TABLE, 
                            prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
                            frags_mouse_detailed_33NCEcorrected_start1, 
                            UInt32(1), #MS_FILE_IDX
                            RT_to_iRT_map, #RT to iRT map
                            min_frag_count = 4, 
                            topN = 1000, 
                            fragment_tolerance = quantile(frag_err_dist, 0.975), 
                            λ = Float32(0), 
                            γ =Float32(0),
                            max_peaks = 10000, 
                            scan_range = (0, 300000), #101357 #22894
                            precursor_tolerance = 20.0,
                            min_spectral_contrast =  Float32(0.5),
                            min_matched_ratio = Float32(0.45),
                            rt_tol = Float32(20.0),
                            frag_ppm_err = frag_err_dist.μ
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



@time rankPSMs!(PSMs, features, colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 0.25, n_folds = 2, num_round = 200, eta = 0.0375, max_depth = 5)
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:]

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