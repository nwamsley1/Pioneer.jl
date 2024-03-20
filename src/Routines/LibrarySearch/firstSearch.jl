
function mainLibrarySearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Arrow.Table,
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    err_dist::MassErrorModel{Float32},
    irt_tol::Float64,
    params::NamedTuple,
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
    precs::Vector{Counter{UInt32, UInt8}}) where {S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}


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

        filter_by_count = Bool(params[:first_search_params]["filter_on_frag_count"]),
        filter_by_rank = Bool(params[:first_search_params]["filter_by_rank"]),


        isotope_err_bounds = Tuple([Int64(x) for x in params[:isotope_err_bounds]]),
        expected_matches = params[:expected_matches],
        min_frag_count = Int64(params[:first_search_params]["min_frag_count"]),
        min_log2_matched_ratio = Float32(params[:first_search_params]["min_log2_matched_ratio"]),
        min_index_search_score = UInt8(params[:first_search_params]["min_index_search_score"]),
        min_topn_of_m = Tuple([Int64(x) for x in params[:first_search_params]["min_topn_of_m"]]),
        min_max_ppm = Tuple([Float32(x) for x in params[:frag_tol_params]["frag_tol_bounds"]]),#(10.0f0, 30.0f0),
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        irt_tol = irt_tol
       
    )
end

PSMs_Dict = Dictionary{String, DataFrame}()
main_search_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
    @time PSMs = vcat(mainLibrarySearch(
        MS_TABLE,
        prosit_lib["f_index"],
        prosit_lib["precursors"],
        library_fragment_lookup_table,
        RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
        UInt32(ms_file_idx), #MS_FILE_IDX
        frag_err_dist_dict[ms_file_idx],
        irt_errs[ms_file_idx],
        params_,
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
        precs
    #scan_range = (100000, 100010)
    )...);

    addMainSearchColumns!(PSMs, MS_TABLE, 
                        prosit_lib["precursors"][:sequence],
                        prosit_lib["precursors"][:missed_cleavages],
                        prosit_lib["precursors"][:is_decoy],
                        prosit_lib["precursors"][:irt],
                        prosit_lib["precursors"][:prec_charge]);
    #Observed iRT estimates based on pre-search
    PSMs[!,:iRT_observed] = RT_to_iRT_map_dict[ms_file_idx](PSMs[!,:RT])
    PSMs[!,:iRT_error] = Float16.(abs.(PSMs[!,:iRT_observed] .- PSMs[!,:iRT_predicted]))

    column_names = [:spectral_contrast,:scribe,:city_block,:entropy_score,
                    :iRT_error,:missed_cleavage,:Mox,:charge,:TIC,
                    :y_count,:err_norm,:spectrum_peak_count,:intercept]

    scoreMainSearchPSMs!(PSMs,
                                column_names,
                                n_train_rounds = params_[:first_search_params]["n_train_rounds_probit"],
                                max_iter_per_round = params_[:first_search_params]["max_iter_probit"],
                                max_q_value = params_[:first_search_params]["max_q_value_probit_rescore"]);

    getProbs!(PSMs);
    
    getBestPSMs!(PSMs,
                    prosit_lib["precursors"][:mz],
                    max_q_value = Float64(params_[:first_search_params]["max_q_value_filter"]),
                    max_psms = Int64(params_[:first_search_params]["max_precursors_passing"])
                )

    insert!(PSMs_Dict, 
        file_id_to_parsed_name[ms_file_idx], 
        PSMs
    );
end


println("Finished main search in ", main_search_time.time, "seconds")
println("Finished main search in ", main_search_time, "seconds")


