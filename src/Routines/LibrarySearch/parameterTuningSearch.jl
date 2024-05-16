function rtAlignSearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Arrow.Table,
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    min_max_ppm::Tuple{Float32, Float32},
    params::NamedTuple,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    collect_fmatches::Bool,
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

    #constant_spline =  approximate(x->0.0f0, BSplineBasis(BSplineOrder(1), range(0.0, 1.0)), MinimiseL2Error())

    err_dist = MassErrorModel(
                0.0f0,
                (first(min_max_ppm), last(min_max_ppm))
    )                                
    return rtAlignSearch(
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
        collect_fmatches = collect_fmatches,
        expected_matches = params[:expected_matches],
        isotope_err_bounds = Tuple([Int64(x) for x in params[:isotope_err_bounds]]),
        min_index_search_score = UInt8(params[:presearch_params]["min_index_search_score"]),
        min_frag_count = Int64(params[:presearch_params]["min_frag_count"]),
        min_log2_matched_ratio = Float32(params[:presearch_params]["min_log2_matched_ratio"]),
        min_spectral_contrast = Float32(params[:presearch_params]["min_spectral_contrast"]),
        min_topn_of_m = Tuple([Int64(x) for x in params[:presearch_params]["min_topn_of_m"]]),
        min_max_ppm = min_max_ppm,
        #min_max_ppm = ( Float32(params[:presearch_params]["frag_tol_ppm"]),  
        #                Float32(params[:presearch_params]["frag_tol_ppm"])
        #             ),
        max_best_rank = Int64(params[:presearch_params]["max_best_rank"]),
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        sample_rate = Float64(params[:presearch_params]["sample_rate"]),
    )
end
function rtAlignSearch(
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
                    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
                    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
                    ionTemplates::Vector{Vector{L}},
                    iso_splines::IsotopeSplineModel,
                    scored_PSMs::Vector{Vector{S}},
                    unscored_PSMs::Vector{Vector{Q}},
                    spectral_scores::Vector{Vector{R}},
                    precursor_weights::Vector{Vector{Float32}},
                    precs::Union{Missing, Vector{Counter{UInt32, UInt8}}};
                    #keyword args
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float32 = 0.0f0,

                    δ::Float32 = 10000f0,
                    λ::Float32 = 0f0,
                    max_iter_newton::Int64 = 100,
                    max_iter_bisection::Int64 = 100,
                    max_iter_outer::Int64 = 100,
                    accuracy_newton::Float32 = 100f0,
                    accuracy_bisection::Float32 = 100000f0,
                    max_diff::Float32 = 0.01f0,

                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
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
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    irt_tol::Float64 = Inf,
                    sample_rate::Float64 = Inf,
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
    thread_tasks, total_peaks = partitionScansToThreads(spectra[:masses],
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
                                min_max_ppm,#(5.0f0, 5.0f0),#min_max_ppm,
                                quadrupole_isolation_width,
                                irt_tol,
                                sample_rate,
                                spec_order
                            )
        end
    end
    precursors_passed_scoring = fetch.(tasks)
    #return precursors_passed_scoring, scan_to_prec_idx
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
function massErrorSearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    scan_idxs::Vector{UInt32},
    precursors_passed_scoring::Vector{UInt32},
    library_fragment_lookup::LibraryFragmentLookup{Float32},
    ms_file_idx::UInt32,
    min_max_ppm::Tuple{Float32, Float32},
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}}
    )

    #Sort scans and precursors 
    sorted_scan_indices = sortperm(scan_idxs)
    scan_idxs = scan_idxs[sorted_scan_indices]
    precursors_passed_scoring = precursors_passed_scoring[sorted_scan_indices]
    #constant_spline =  approximate(x->0.0f0, BSplineBasis(BSplineOrder(1), range(0.0, 1.0)), MinimiseL2Error())
    function getScanToPrecIdx(scan_idxs::Vector{UInt32}, n_scans::Int64)
        scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, n_scans)
        start_idx, stop_idx = 1, 1
        for i in range(1, length(scan_idxs))
            stop_idx = i 
            if scan_idxs[start_idx] == scan_idxs[stop_idx]
                scan_to_prec_idx[scan_idxs[i]] = start_idx:stop_idx
            else
                scan_to_prec_idx[scan_idxs[i]] = i:i
                start_idx = i
            end
        end
        return scan_to_prec_idx
    end

    scan_to_prec_idx = getScanToPrecIdx(scan_idxs, length(spectra[:masses]))

    mass_err_model = MassErrorModel(
                0.0f0,
                (first(min_max_ppm), last(min_max_ppm))
    )               
    return massErrorSearch(
        spectra, 
        library_fragment_lookup,
        scan_idxs,
        scan_to_prec_idx,
        precursors_passed_scoring,
        ms_file_idx,
        mass_err_model,
        ionMatches,
        ionMisses,
        all_fmatches,
        ionTemplates,
    )
end

function massErrorSearch(
                    #Mandatory Args
                    spectra::Arrow.Table,
                    library_fragment_lookup::LibraryFragmentLookup{Float32},
                    scan_idxs::Vector{UInt32},
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    precursors_passed_scoring::Vector{UInt32},
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel,
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionTemplates::Vector{Vector{L}}
                    ) where {L<:LibraryIon{Float32}}


    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
    
    #Build thread tasks 
    thread_task_size = length(scan_idxs)÷Threads.nthreads()
    thread_tasks = []
    start_idx, stop_idx =1, min(thread_task_size - 1, length(scan_idxs))
    for i in range(1, Threads.nthreads())
        push!(thread_tasks, (i, start_idx:stop_idx))
        start_idx = stop_idx + 1
        stop_idx =  min(start_idx + thread_task_size - 1, length(scan_idxs))
    end

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getMassErrors(
                                spectra,
                                library_fragment_lookup,
                                last(thread_task),
                                scan_idxs,
                                scan_to_prec_idx,
                                precursors_passed_scoring,
                                ms_file_idx,
                                mass_err_model,
                                ionMatches[thread_id],
                                ionMisses[thread_id],
                                all_fmatches[thread_id],
                                ionTemplates[thread_id]
                            )
        end
    end
    psms = fetch.(tasks)
    return psms
end

#Clear directory of previous QC Plots
[rm(joinpath(rt_alignment_folder, x)) for x in readdir(rt_alignment_folder)]
[rm(joinpath(mass_err_estimation_folder, x)) for x in readdir(mass_err_estimation_folder)]

test_time = @time begin
RT_to_iRT_map_dict = Dict{Int64, Any}()
frag_err_dist_dict = Dict{Int64,MassErrorModel}()
irt_errs = Dict{Int64, Float64}()
#MS_TABLE_PATH = "/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E30H50Y20_30SPD_DIA_4.arrow"
for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    out_fname = String(first(split(splitpath(MS_TABLE_PATH)[end],".")));
    #Randomly sample spectra to search and retain only the 
    #most probable psms as specified in "first_seach_params"
    start_ppm = Float32(params_[:presearch_params]["frag_tol_ppm"])
    data_points = 0
    n = 0
    rtPSMs = nothing
    while n <= params_[:presearch_params]["max_presearch_iters"]
        @time RESULT =  rtAlignSearch(
                                                MS_TABLE,
                                                prosit_lib["presearch_f_index"],
                                                prosit_lib["precursors"],
                                                library_fragment_lookup_table,
                                                x->x, #RT to iRT map'
                                                UInt32(ms_file_idx), #MS_FILE_IDX
                                                (start_ppm, start_ppm),
                                                params_,
                                                ionMatches,
                                                ionMisses,
                                                false,
                                                all_fmatches,
                                                IDtoCOL,
                                                ionTemplates,
                                                iso_splines,
                                                scored_PSMs,
                                                unscored_PSMs,
                                                spectral_scores,
                                                precursor_weights,
                                                precs,
                                                );
        PSMs = vcat([result for result in RESULT]...)
        addPreSearchColumns!(PSMs, 
                                    MS_TABLE, 
                                    prosit_lib["precursors"][:is_decoy],
                                    prosit_lib["precursors"][:irt],
                                    prosit_lib["precursors"][:prec_charge]
                                )
        if rtPSMs === nothing
            rtPSMs = PSMs#vcat([result for result in RESULT]...)
        else
            rtPSMs = vcat(rtPSMs, PSMs)
        end            

        scorePresearch!(rtPSMs)
        getQvalues!(rtPSMs[!,:prob], rtPSMs[!,:target], rtPSMs[!,:q_value])
        #println("TEST ",  sum(rtPSMs[!,:q_value].<=params_[:presearch_params]["max_qval"]))
        if sum(rtPSMs[!,:q_value].<=params_[:presearch_params]["max_qval"]) >= params_[:presearch_params]["min_samples"]
            filter!(:q_value => x -> x<=params_[:presearch_params]["max_qval"], rtPSMs)
            rtPSMs[!,:best_psms] .= false
            grouped_psms = groupby(rtPSMs,:precursor_idx)
            for psms in grouped_psms
                best_idx = argmax(psms.prob)
                psms[best_idx,:best_psms] = true
            end
            filter!(x->x.best_psms, rtPSMs)
            break
        end
        n += 1
    end
    if n >= params_[:presearch_params]["max_presearch_iters"]
        min_samples = params_[:presearch_params]["min_samples"]
        max_iters = params_[:presearch_params]["max_presearch_iters"]
        @warn "Presearch did not find $min_samples precursors at the specified fdr in $max_iters iterations"
    end

    RT_to_iRT_map = UniformSpline( 
                                    rtPSMs[!,:iRT_predicted], 
                                    rtPSMs[!,:RT], 
                                    3, #Degree is always three
                                    5 #Spline knots fixed. currently not tunable parameter
                                    );

    plotRTAlign(rtPSMs[:,:RT], 
                rtPSMs[:,:iRT_predicted], 
                RT_to_iRT_map, 
                out_fdir = rt_alignment_folder,
                out_fname = out_fname
                );
    rtPSMs[!,:iRT_observed] = RT_to_iRT_map.(rtPSMs[!,:RT])
    irt_MAD = mad(rtPSMs[!,:iRT_observed] .- rtPSMs[!,:iRT_predicted])

    irt_σ = irt_MAD #Robust estimate of standard deviation

    irt_errs[ms_file_idx] = params_[:irt_err_sigma]*irt_σ
    RT_to_iRT_map_dict[ms_file_idx] = RT_to_iRT_map


    matched_fragments = vcat(massErrorSearch(
        MS_TABLE,
        rtPSMs[!,:scan_idx],
        rtPSMs[!,:precursor_idx],
        library_fragment_lookup_table,
        UInt32(ms_file_idx),
        (start_ppm, start_ppm),
        ionMatches,
        ionMisses,
        all_fmatches,
        ionTemplates
    )...);

    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        Float32((a-b)/(a/1e6))
    end
    #Get Retention Times and Target/Decoy Status 
    ####################
    #Use best_psms to estimate 
    #1) RT to iRT curve and 
    #2) mass error (ppm) distribution 

    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in matched_fragments];
    mass_err_model = ModelMassErrs(
        frag_ppm_errs,
        frag_err_quantile = Float32(params_[:presearch_params]["frag_err_quantile"]),
        out_fdir = mass_err_estimation_folder,
        out_fname = out_fname
    )

    #PLOT_PATH = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS", split(splitpath(MS_TABLE_PATH)[end],".")[1])
    #File name but remove file type
    frag_err_dist_dict[ms_file_idx] = mass_err_model
end
end

#Merge Quality Control PDFs 
merge_pdfs([joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) if endswith(x, ".pdf")], 
                joinpath(rt_alignment_folder, "rt_alignment_plots.pdf"), cleanup=true)
merge_pdfs([joinpath(mass_err_estimation_folder, x) for x in readdir(mass_err_estimation_folder) if endswith(x, ".pdf")], 
    joinpath(mass_err_estimation_folder, "mass_err_estimation_folder.pdf"), cleanup=true)
