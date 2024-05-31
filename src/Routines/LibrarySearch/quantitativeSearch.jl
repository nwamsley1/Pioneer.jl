function quantitationSearch(
    #Mandatory Args
    spectra::Arrow.Table,
    precursors::Arrow.Table,
    ion_list::LibraryFragmentLookup{Float32},
    rt_index::retentionTimeIndex{Float32, Float32},
    ms_file_idx::UInt32,
    rt_to_irt::UniformSpline,
    err_dist::MassErrorModel,
    irt_tol::Float64,
    params::NamedTuple,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel,
    chromatograms::Vector{Vector{ChromObject}},
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}}
    ) where {S<:ScoredPSM{Float32, Float16},
                                            Q<:UnscoredPSM{Float32},
                                            R<:SpectralScores{Float16}}

    #fragment_tolerance = quantile(err_dist, params[:frag_tol_quantile])
    return quantitationSearch(
        spectra, 
        precursors,
        ion_list, 
        ms_file_idx,
        rt_to_irt,
        err_dist,     
        ionMatches,
        ionMisses,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        chromatograms,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        
        isotope_err_bounds = params[:isotope_err_bounds],

        δ = Float32(params[:deconvolution_params]["huber_delta"]),
        λ = Float32(params[:deconvolution_params]["lambda"]),
        max_iter_newton = Int64(params[:deconvolution_params]["max_iter_newton"]),
        max_iter_bisection = Int64(params[:deconvolution_params]["max_iter_bisection"]),
        max_iter_outer = Int64(params[:deconvolution_params]["max_iter_outer"]),
        accuracy_newton = Float32(params[:deconvolution_params]["accuracy_newton"]),
        accuracy_bisection = Float32(params[:deconvolution_params]["accuracy_bisection"]),
        max_diff = Float32(params[:deconvolution_params]["max_diff"]),

        min_topn_of_m = Tuple([Int64(x) for x in params[:quant_search_params]["min_topn_of_m"]]),
        min_frag_count = Int64(params[:quant_search_params]["min_frag_count"]),
        min_log2_matched_ratio = Float32(params[:quant_search_params]["min_log2_matched_ratio"]),
        min_spectral_contrast = Float32(params[:quant_search_params]["min_spectral_contrast"]),
        min_max_ppm = Tuple([Float32(x) for x in params[:frag_tol_params]["frag_tol_bounds"]]),
        n_frag_isotopes = Int64(params[:quant_search_params]["n_frag_isotopes"]),
        max_best_rank = Int64(params[:quant_search_params]["max_best_rank"]),

        rt_index = rt_index,
        irt_tol = irt_tol,
    )
end
function quantitationSearch(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    precursors::Union{Arrow.Table, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    ms_file_idx::UInt32,
                    rt_to_irt::UniformSpline,
                    mass_err_model::MassErrorModel,
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
                    ionTemplates::Vector{Vector{L}},
                    iso_splines::IsotopeSplineModel,
                    chromatograms::Vector{Vector{ChromObject}},
                    scored_PSMs::Vector{Vector{S}},
                    unscored_PSMs::Vector{Vector{Q}},
                    spectral_scores::Vector{Vector{R}},
                    precursor_weights::Vector{Vector{Float32}};
                    #keyword args
                    frag_ppm_err::Float32 = 0.0f0,

                    δ::Float32 = 10000f0,
                    λ::Float32 = 0f0,
                    max_iter_newton::Int64 = 100,
                    max_iter_bisection::Int64 = 100,
                    max_iter_outer::Int64 = 100,
                    accuracy_newton::Float32 = 100f0,
                    accuracy_bisection::Float32 = 100000f0,
                    max_diff::Float32 = 0.01f0,

                    isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
                    min_frag_count::Int64 = 1,
                    min_spectral_contrast::Float32  = 0f0,
                    min_log2_matched_ratio::Float32 = -Inf32,
                    min_topn_of_m::Tuple{Int64, Int64} = (2, 3),
                    min_max_ppm::Tuple{Float32, Float32} = (-Inf, Inf),
                    max_best_rank::Int64 = one(Int64),
                    n_frag_isotopes::Int64 = 1,
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    irt_tol::Float64 = Inf,
                    spec_order::Set{Int64} = Set(2)
                    ) where {T,U<:AbstractFloat, 
                                                            L<:LibraryIon{Float32}, 
                                                            S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}


    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
                                       
    
    thread_tasks, total_peaks = partitionScansToThreads2(spectra[:masses],
                                                        spectra[:retentionTime],
                                                        spectra[:centerMass],
                                                        spectra[:msOrder],

                                                        Threads.nthreads(),
                                                        1)
    
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return quantPSMs(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                precursors,
                                ion_list, 
                                ms_file_idx,
                                rt_to_irt,
                                mass_err_model,
                                frag_ppm_err,
                                δ,
                                λ,
                                max_iter_newton,
                                max_iter_bisection,
                                max_iter_outer,
                                accuracy_newton,
                                accuracy_bisection,
                                max_diff,
                                ionMatches[thread_id],
                                ionMisses[thread_id],
                                IDtoCOL[thread_id],
                                ionTemplates[thread_id],
                                iso_splines,
                                chromatograms[thread_id],
                                scored_PSMs[thread_id],
                                unscored_PSMs[thread_id],
                                spectral_scores[thread_id],
                                precursor_weights[thread_id],
                                isotope_err_bounds,
                                min_frag_count,
                                min_spectral_contrast,
                                min_log2_matched_ratio,
                                min_topn_of_m,
                                min_max_ppm,
                                max_best_rank,
                                n_frag_isotopes,
                                rt_index, 
                                irt_tol,
                                spec_order
                            )
        end
    end
    psms = fetch.(tasks)
    return psms
end
BPSMS = Dict{Int64, DataFrame}()
#PSMS_DIR = joinpath(MS_DATA_DIR,"Search","RESULTS")
#PSM_PATHS = [joinpath(PSMS_DIR, file) for file in filter(file -> isfile(joinpath(PSMS_DIR, file)) && match(r".jld2$", file) != nothing, readdir(PSMS_DIR))];



quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    params_[:deconvolution_params]["huber_delta"] = median(
        [quantile(x, 0.25) for x in MS_TABLE[:intensities]])*params_[:deconvolution_params]["huber_delta_prop"] 
        @time RESULT = quantitationSearch(MS_TABLE, 
            prosit_lib["precursors"],
            prosit_lib["f_det"],
            RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
            UInt32(ms_file_idx), 
            RT_iRT[file_id_to_parsed_name[ms_file_idx]],
            frag_err_dist_dict[ms_file_idx],
            irt_err,#irt_errs[ms_file_idx]/3,
            params_,  
            ionMatches,
            ionMisses,
            IDtoCOL,
            ionTemplates,
            iso_splines,
            chromatograms,
            complex_scored_PSMs,
            complex_unscored_PSMs,
            complex_spectral_scores,
            precursor_weights,
            );

        #Format Chromatograms 
        chroms = vcat([last(x) for x in RESULT]...);
        #sort!(chroms,:rt,alg=QuickSort)
        getIsotopesCaptured!(chroms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
        gchroms = groupby(chroms,[:precursor_idx,:isotopes_captured])

        #Format spectral scores 
        psms = vcat([first(x) for x in RESULT]...);
        addSecondSearchColumns!(psms, 
                                        MS_TABLE, 
                                        prosit_lib["precursors"][:mz],
                                        prosit_lib["precursors"][:prec_charge], 
                                        prosit_lib["precursors"][:is_decoy],
                                        precID_to_cv_fold);
        features = [:intercept, :charge, :total_ions, :err_norm, 
        :scribe, :city_block, :city_block_fitted, 
        :spectral_contrast, :entropy_score, :weight]
        if sum(psms[!,:p_count])>0
            push!(features, :p_count)
        end
        getIsotopesCaptured!(psms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
        psms[!,:prob] = zeros(Float32, size(psms, 1));
        scoreSecondSearchPSMs!(psms,features);
        spectral_scores = groupby(psms, [:precursor_idx,:isotopes_captured]);

        #Integrate Chromatograms and Combine with Spectral Simmilarity Scores
        precs_to_integrate = keys(spectral_scores) ∩ keys(gchroms)
        n_precs = length(precs_to_integrate)
        psms_integrated = initQuantScans(precs_to_integrate)


        integratePrecursors(spectral_scores,
                            gchroms,
                            precs_to_integrate,
                            psms_integrated,
                            λ=Float32(params_[:quant_search_params]["WH_smoothing_strength"]))
        #Remove examples where maximum point on chromatogram doesn't have a psm
        filter!(x->!ismissing(x.scan_idx), psms_integrated)
        #Add additional features 
        addPostIntegrationFeatures!(psms_integrated, 
                                    MS_TABLE, 
                                    prosit_lib["precursors"][:sequence],
                                    prosit_lib["precursors"][:structural_mods],
                                    prosit_lib["precursors"][:mz],
                                    prosit_lib["precursors"][:irt],
                                    prosit_lib["precursors"][:prec_charge],
                                    prosit_lib["precursors"][:missed_cleavages],
                                    ms_file_idx,
                                    file_id_to_parsed_name,
                                    RT_iRT,
                                    precID_to_iRT
                                    );
        psms_integrated[!,:file_name].=file_id_to_parsed_name[ms_file_idx]
        BPSMS[ms_file_idx] = psms_integrated;
end

best_psms = vcat(values(BPSMS)...)