function firstSearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Vector{LibraryPrecursorIon{Float32}},
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, UInt8}}
    ) where {S<:ScoredPSM{Float32, Float16},
    Q<:UnscoredPSM{Float32},
    R<:SpectralScores{Float16}}#where {S<:ScoredPSM{Float32, Float16}, LibraryIon{Float32}}
    err_dist = MassErrorModel((zero(Float32), zero(Float32), zero(Float32)), zero(Float32))
    return searchRAW(
        spectra, 
        frag_index,
        precursors, 
        ion_list,
        iRT_to_RT_spline,
        ms_file_idx,
        err_dist,
        searchScan!,

        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        precs,

        collect_fmatches = true,
        expected_matches = params[:expected_matches],
        frag_ppm_err = Float32(params[:frag_ppm_err]),
        min_frag_count = params[:min_frag_count],
        min_log2_matched_ratio = params[:min_log2_matched_ratio],
        min_index_search_score = params[:min_index_search_score],
        min_topn_of_m = params[:min_topn_of_m],
        min_max_ppm = (40.0f0, 40.0f0),
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_bounds = params[:rt_bounds],
        irt_tol = params[:irt_tol],
        sample_rate = params[:sample_rate],
    )
end

RT_to_iRT_map_dict = Dict{Int64, Any}()
frag_err_dist_dict = Dict{Int64,MassErrorModel}()
lk = ReentrantLock()
for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    #Randomly sample spectra to search and retain only the 
    #most probable psms as specified in "first_seach_params"
    @time RESULT =  firstSearch(
                                            MS_TABLE,
                                            prosit_lib["f_index"],
                                            prosit_lib["precursors"],
                                            prosit_lib["f_det"],
                                            x->x, #RT to iRT map'
                                            UInt32(ms_file_idx), #MS_FILE_IDX
                                            Laplace(zero(Float64), first_search_params[:frag_tol_presearch]),
                                            first_search_params,
                                            ionMatches,
                                            ionMisses,
                                            all_fmatches,
                                            IDtoCOL,
                                            ionTemplates,
                                            iso_splines,
                                            scored_PSMs,
                                            unscored_PSMs,
                                            spectral_scores,
                                            precursor_weights,
                                            precs,
                                            #scan_range = (100000, 110000)
                                            scan_range = (1, length(MS_TABLE[:masses]))
                                            );
    rtPSMs = vcat([first(result) for result in RESULT]...)
    all_matches = vcat([last(result) for result in RESULT]...)
    #@time begin
    @time addPreSearchColumns!(rtPSMs, MS_TABLE, precursors)
    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        (a-b)/(a/1e6)
    end
    #Get Retention Times and Target/Decoy Status 
    ####################
    #Use best_psms to estimate 
    #1) RT to iRT curve and 
    #2) mass error (ppm) distribution 
    @time begin
    best_precursors = Set(rtPSMs[:,:precursor_idx]);
    best_matches = [match for match in all_matches if match.prec_id ∈ best_precursors];
    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in best_matches];
    frag_ppm_intensities = [match.intensity for match in best_matches];
    end

    mass_err_model = ModelMassErrs(
        frag_ppm_intensities,
        frag_ppm_errs,
        40.0
    )
    
    #Model fragment errors with a mixture model of a uniform and laplace distribution 
    lock(lk) do 
        PLOT_PATH = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS", split(splitpath(MS_TABLE_PATH)[end],".")[1])
        @time RT_to_iRT_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT_predicted], n = 50, bandwidth = 4.0);
        @time plotRTAlign(rtPSMs[:,:RT], rtPSMs[:,:iRT_predicted], RT_to_iRT_map, 
                    f_out = PLOT_PATH);
        RT_to_iRT_map_dict[ms_file_idx] = RT_to_iRT_map
        frag_err_dist_dict[ms_file_idx] = mass_err_model
    end
end

jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "rt_map_dict_020824.jld2"); RT_to_iRT_map_dict)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "frag_err_dist_dict_020824.jld2"); frag_err_dist_dict)

