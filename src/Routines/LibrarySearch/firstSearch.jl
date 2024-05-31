function mainLibrarySearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Arrow.Table,
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    err_dist::MassErrorModel,
    irt_tol::Float64,
    params::NamedTuple,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel,
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, UInt8}}) where {S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}


    return mainLibrarySearch(
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
        min_frag_count = Int64(params[:first_search_params]["min_frag_count"]),
        n_frag_isotopes = Int64(params[:first_search_params]["n_frag_isotopes"]),
        min_log2_matched_ratio = Float32(params[:first_search_params]["min_log2_matched_ratio"]),
        min_spectral_contrast = Float32(params[:first_search_params]["min_spectral_contrast"]),
        min_index_search_score = UInt8(params[:first_search_params]["min_index_search_score"]),
        min_topn_of_m = Tuple([Int64(x) for x in params[:first_search_params]["min_topn_of_m"]]),
        min_max_ppm = Tuple([Float32(x) for x in params[:frag_tol_params]["frag_tol_bounds"]]),#(10.0f0, 30.0f0),
        irt_tol = irt_tol
       
    )
end
function mainLibrarySearch(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Arrow.Table, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt_spline::Any,
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel,
                    searchScan!::Union{Function, Missing},
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
                    ionTemplates::Vector{Vector{L}},
                    iso_splines::IsotopeSplineModel,
                    scored_PSMs::Vector{Vector{S}},
                    unscored_PSMs::Vector{Vector{Q}},
                    spectral_scores::Vector{Vector{R}},
                    precursor_weights::Vector{Vector{Float32}},
                    precs::Union{Missing, Vector{Counter{UInt32, UInt8}}};
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
                    min_index_search_score::UInt8 = zero(UInt8),
                    min_topn_of_m::Tuple{Int64, Int64} = (2, 3),
                    min_max_ppm::Tuple{Float32, Float32} = (-Inf, Inf),
                    filter_by_rank::Bool = false, 
                    filter_by_count::Bool = true,
                    max_best_rank::Int64 = one(Int64),
                    n_frag_isotopes::Int64 = 1,
                    quadrupole_isolation_width::Float64 = 8.5,
                    irt_tol::Float64 = Inf,
                    sample_rate::Float64 = Inf,
                    spec_order::Set{Int64} = Set(2)
                    ) where { 
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

    if ismissing(precs)
        precs = [missing for _ in range(1, Threads.nthreads())]
    end
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return searchFragmentIndex(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                frag_index,
                                scan_to_prec_idx,
                                rt_to_irt_spline,
                                mass_err_model,
                                searchScan!,
                                frag_ppm_err,
                                precs[thread_id],
                                isotope_err_bounds,
                                min_index_search_score,
                                min_max_ppm,
                                quadrupole_isolation_width,
                                irt_tol,
                                sample_rate,
                                spec_order
                            )
        end
    end
    precursors_passed_scoring = fetch.(tasks)
    #return scan_to_prec_idx, precursors_passed_scoring, thread_tasks
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getPSMS(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                precursors,
                                scan_to_prec_idx,
                                precursors_passed_scoring[thread_id],
                                ion_list, 
                                rt_to_irt_spline,
                                ms_file_idx,
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
                                scored_PSMs[thread_id],
                                unscored_PSMs[thread_id],
                                spectral_scores[thread_id],
                                precursor_weights[thread_id],
                                precs[thread_id],
                                isotope_err_bounds,
                                min_frag_count,
                                min_spectral_contrast,
                                min_log2_matched_ratio,
                                min_topn_of_m,
                                min_max_ppm,
                                filter_by_rank,
                                filter_by_count,
                                max_best_rank,
                                n_frag_isotopes,
                                irt_tol,
                                spec_order
                            )
        end
    end
    psms = fetch.(tasks)
    return psms
end
#frag_err_dist_dict[1] = MassErrorModel(frag_err_dist_dict[1].mass_offset, frag_err_dist_dict[1].mass_tolerance, 7.0f0, 15.0f0)
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
                        prosit_lib["precursors"][:structural_mods],
                        prosit_lib["precursors"][:missed_cleavages],
                        prosit_lib["precursors"][:is_decoy],
                        prosit_lib["precursors"][:irt],
                        prosit_lib["precursors"][:prec_charge]);
    
    #Observed iRT estimates based on pre-search
    PSMs[!,:iRT_observed] = RT_to_iRT_map_dict[ms_file_idx].(PSMs[!,:RT])
    PSMs[!,:iRT_error] = Float16.(abs.(PSMs[!,:iRT_observed] .- PSMs[!,:iRT_predicted]))
    
    column_names = [:spectral_contrast,:scribe,:city_block,:entropy_score,
                    :iRT_error,:missed_cleavage,:Mox,
                    :charge,:TIC,
                    :y_count,:err_norm,:spectrum_peak_count,:intercept]
    if sum(PSMs[!,:p_count])>0
        push!(column_names, :p_count)
    end
    scoreMainSearchPSMs!(PSMs,
                                column_names,
                                n_train_rounds = params_[:first_search_params]["n_train_rounds_probit"],
                                max_iter_per_round = params_[:first_search_params]["max_iter_probit"],
                                max_q_value = params_[:first_search_params]["max_q_value_probit_rescore"]);

    getProbs!(PSMs);
    #=
    bins = LinRange(0, 2, 100)
    histogram(PSMs[PSMs[!,:target].&(PSMs[!,:q_value].<=0.01), :entropy_score], alpha = 0.5, bins = bins, normalize = :pdf)
    histogram!(PSMs[PSMs[!,:target].==false, :entropy_score], alpha = 0.5, bins = bins, normalize=:pdf)
    

    bins = LinRange(0, 2, 100)
    histogram(PSMs[PSMs[!,:target].&(PSMs[!,:q_value].<=0.01), :spectral_contrast], alpha = 0.5, bins = bins, normalize = :pdf)
    histogram!(PSMs[PSMs[!,:target].==false, :spectral_contrast], alpha = 0.5, bins = bins, normalize=:pdf)

    bins = LinRange(0, 25, 100)
    histogram(PSMs[PSMs[!,:target].&(PSMs[!,:q_value].<=0.01), :total_ions], alpha = 0.5, bins = bins, normalize = :pdf)
    histogram!(PSMs[PSMs[!,:target].==false, :total_ions], alpha = 0.5, bins = bins, normalize=:pdf)

    bins = LinRange(-1, 5, 100)
    histogram(PSMs[PSMs[!,:target].&(PSMs[!,:q_value].<=0.01), :matched_ratio], alpha = 0.5, bins = bins, normalize = :pdf)
    histogram!(PSMs[PSMs[!,:target].==false, :matched_ratio], alpha = 0.5, bins = bins, normalize=:pdf)

    bins = LinRange(0, 25, 100)
    histogram(PSMs[PSMs[!,:target].&(PSMs[!,:q_value].<=0.01), :topn], alpha = 0.5, bins = bins, normalize = :pdf)
    histogram!(PSMs[PSMs[!,:target].==false, :topn], alpha = 0.5, bins = bins, normalize=:pdf)
    =#
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


#println("Finished main search in ", main_search_time.time, "seconds")
#println("Finished main search in ", main_search_time, "seconds")
#=
@time scan_to_prec_idx, precursors_passed_scoring, thread_tasks = mainLibrarySearch(
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
);

@time subset_lookup_table = buildSubsetLibraryFragmentLookupTable!(
    precursors_passed_scoring,
    library_fragment_lookup_table,
    length(precursors[:irt])
)
precursor_sets = [x for x in precursors_passed_scoring]
union(precursor_sets...)
[scan_to_prec_idx[thread_tasks[1][end][x]] for x in range(3012, 4000)]
scan_to_prec_idx[]
prec_diffs = zeros(Union{Missing, Int64}, length(thread_tasks[1][end]) - 1)
prec_totals = zeros(Union{Missing, Int64}, length(thread_tasks[1][end]) - 1)
for i in range(1, length(thread_tasks[1][end]) - 1)
    if iszero(thread_tasks[1][end][i])
        prec_diffs[i] = missing
        continue
    end
    prec_range_a = scan_to_prec_idx[thread_tasks[1][end][i]]
    prec_range_b = scan_to_prec_idx[thread_tasks[1][end][i+1]]
    if ismissing(prec_range_a) | ismissing(prec_range_b)
        prec_diffs[i] = missing
        continue
    end
    if (last(prec_range_a) - first(prec_range_a) > 0) & (last(prec_range_b) - first(prec_range_b) > 0)
        a = Set(precursors_passed_scoring[1][prec_range_a])
        b = Set(precursors_passed_scoring[1][prec_range_b])
        prec_totals[i] = length(a) + length(b)
        println(length(setdiff(b, a)))
        prec_diffs[i] = length(setdiff(b, a))
    end
end
to_keep = iszero.(prec_totals).==false

to_keep = prec_totals.>30
prec_totals = prec_totals[to_keep]
prec_diffs = prec_diffs[to_keep]
histogram(1.0 .- prec_diffs./prec_totals)

a = Set(precursors_passed_scoring[1][scan_to_prec_idx[thread_tasks[1][end][3013]]])
b = Set(precursors_passed_scoring[1][scan_to_prec_idx[thread_tasks[1][end][3014]]])
length(b)
length(a)
length(setdiff(a, b))
length(setdiff(b, a))

=#