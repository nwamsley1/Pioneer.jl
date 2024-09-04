function huberLossSearch(
    frag_err_dist_dict,
    rt_index_paths,
    prec_set,
    scan_idxs,
    δs,
    bin_rt_size,
    rt_irt,
    irt_errs,
    chromatograms,
    file_path_to_parsed_name,
    ms_file_idx,
    MS_TABLE_PATH,
    params_,
    spec_lib,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    complex_scored_psms,
    complex_unscored_psms,
    complex_spectral_scores,
    precursor_weights
    )

    function huberLossSearch(
                        #Mandatory Args
                        spectra::Arrow.Table, 
                        params::Any,
                        prec_set::Set{Tuple{UInt32, UInt32}},
                        scan_idxs::Set{UInt32},
                        δs::Vector{Float32};
                        kwargs...)
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
        
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin 
                thread_id = first(thread_task)
                return huberTuningSearch(
                                    spectra,
                                    last(thread_task), #getRange(thread_task),
                                    prec_set,
                                    scan_idxs,
                                    spec_lib["precursors"],
                                    kwargs[:fragment_lookup_table], 
                                    kwargs[:ms_file_idx],
                                    kwargs[:rt_to_irt_spline],
                                    kwargs[:mass_err_model],
                                    kwargs[:quad_transmission_func],
                                    δs,
                                    Float32(params[:deconvolution_params]["lambda"]),
                                    Int64(params[:deconvolution_params]["max_iter_newton"]),
                                    Int64(params[:deconvolution_params]["max_iter_bisection"]),
                                    Int64(params[:deconvolution_params]["max_iter_outer"]),
                                    Float32(params[:deconvolution_params]["accuracy_newton"]),
                                    Float32(params[:deconvolution_params]["accuracy_bisection"]),
                                    Float32(params[:deconvolution_params]["max_diff"]),
                                    kwargs[:ion_matches][thread_id],
                                    kwargs[:ion_misses][thread_id],
                                    kwargs[:id_to_col][thread_id],
                                    kwargs[:ion_templates][thread_id],
                                    kwargs[:iso_splines],
                                    kwargs[:unscored_psms][thread_id],
                                    kwargs[:spectral_scores][thread_id],
                                    kwargs[:precursor_weights][thread_id],
                                    (1, 0),
                                    1,#params[:quant_search_params]["n_frag_isotopes"],
                                    kwargs[:rt_index], 
                                    kwargs[:irt_err],
                                    Set(2),
                                )
            end
        end
        psms = fetch.(tasks)
        return psms
    end

    ms_table_path_to_psms_path = Dict{String, String}()
    parsed_fname = file_path_to_parsed_name[MS_TABLE_PATH]
    rt_df = DataFrame(Arrow.Table(rt_index_paths[parsed_fname]))
    rt_index = buildRtIndex(rt_df,
                            bin_rt_size = bin_rt_size)
    rt_irt = rt_irt[parsed_fname]
    MS_TABLE = Arrow.Table(MS_TABLE_PATH);
        params_[:quant_search_params]["min_y_count"] = 1
        precursors = spec_lib["precursors"]
        psms = vcat(huberLossSearch(
            MS_TABLE, 
            params_,
            prec_set,
            scan_idxs,
            δs;
            precursors = spec_lib["precursors"],
            fragment_lookup_table = spec_lib["f_det"],
            rt_index = rt_index,
            ms_file_idx = UInt32(ms_file_idx), 
            rt_to_irt_spline = rt_irt,
            mass_err_model = frag_err_dist_dict[ms_file_idx],
            irt_err = irt_errs[parsed_fname],#irt_errs[ms_file_idx]/3,
            ion_matches = ionMatches,
            ion_misses = ionMisses,
            id_to_col = IDtoCOL,
            ion_templates = ionTemplates,
            iso_splines = iso_splines,
            chromatograms = chromatograms,
            scored_psms = complex_scored_psms,
            unscored_psms = complex_unscored_psms,
            spectral_scores = complex_spectral_scores,
            precursor_weights = precursor_weights,
            quad_transmission_func = QuadTransmission(params_[:quad_transmission]["overhang"], params_[:quad_transmission]["smoothness"])
            )...);
    psms[!,:ms_file_idx] .= UInt32(ms_file_idx)
    return psms
end

function getHuberLossParam(
    huber_δs::Vector{Float32},
    gbpsms;
    minimum_percent_diff::Float32 = 10.0f0)
    psms = []
    @time for (key, sub_bpsms) in pairs(gbpsms)
        ms_file_idx = key[:ms_file_idx]
        MS_TABLE_PATH = MS_TABLE_PATHS[ms_file_idx]
        prec_set = Set(zip(sub_bpsms[!,:precursor_idx], sub_bpsms[!,:scan_idx]))
        scan_idxs = Set(sub_bpsms[!,:scan_idx])

        push!(psms, huberLossSearch(
            frag_err_dist_dict,
            rt_index_paths,
            prec_set,
            scan_idxs,
            huber_δs,
            bin_rt_size,
            rt_irt,
            irt_errs,
            chromatograms,
            file_path_to_parsed_name,
            ms_file_idx,
            MS_TABLE_PATH,
            params_,
            spec_lib,
            ionMatches,
            ionMisses,
            IDtoCOL,
            ionTemplates,
            iso_splines,
            complex_scored_PSMs,
            complex_unscored_PSMs,
            complex_spectral_scores,
            precursor_weights
            ));
    end
    psms = vcat(psms...);
    gpsms = groupby(psms, [:precursor_idx,:ms_file_idx,:scan_idx])


    function processHuberLossCurve(
        weights::AbstractVector{Float32},
        huber_δs::AbstractVector{Float32})
        min_w,max_w = minimum(weights),maximum(weights)
        huber50 = missing 
        w50 = min_w + (max_w - min_w)/2
        if length(weights)>1
            for i in range(1, length(weights)-1)
                if (w50 >= weights[i]) & (w50 <= weights[i + 1])
                    huber50 =  huber_δs[i] + (huber_δs[i + 1] - huber_δs[i])/2 
                end
            end
        end
        return (min = minimum(weights), max = maximum(weights), n = length(weights), huber50 = huber50, w50 = w50, wdiff = (max_w - min_w)/min_w)
    end

    function getMedianHuberDelta(
        cum_prob::Vector{Float32},
        δ::Vector{Int64};
        cum_prob_threshold = 0.5)
        N = length(cum_prob)
        if N == 1 return first(δ) end
        for i in range(1, N-1)
            if cum_prob[i+1]>=cum_prob_threshold
                x1, x2 = cum_prob[i], cum_prob[i+1]
                y1, y2 = δ[i], δ[i+1]
                slope = (y2 - y1)/(x2 - x1)
                midpoint = ((x2 + x1)/2)
                return y1 + slope*(midpoint - x1)
            end
        end
        @warn "Could not estimate huber delta"
        return first(δ)
    end

    combdf = combine(gpsms) do sdf
        processHuberLossCurve(sdf[!,:weight],sdf[!,:huber_δ])
    end
    #Filter examples with incomplete curves
    filter!(x->x.n==length(huber_δs), combdf)
    #Filter on minimum fold difference between high/low huber_δ value
    filter!(x->x.wdiff>(minimum_percent_diff/100), combdf)
    filter!(x->!ismissing(x.huber50), combdf)
    #Round to integer values 
    combdf[!,:huber50] = ceil.(Int, combdf[!,:huber50])
    #Count examples in each huber_δ bin 
    value_counts(df, col) = combine(groupby(df, col), nrow)
 
    huber_hist = sort(value_counts(combdf,:huber50),:huber50)
    huber_hist[!,:prob] = huber_hist[!,:nrow]./sum(huber_hist[!,:nrow])
    huber_hist[!,:cum_prob] = Float32.(cumsum(huber_hist[!,:prob]))
    return getMedianHuberDelta(huber_hist[!,:cum_prob],huber_hist[!,:huber50])
end

function getPsmsForHuberEstimation(
    prec_to_irt::Dictionary{UInt32, 
    @NamedTuple{best_prob::Float32, 
                best_ms_file_idx::UInt32, 
                best_scan_idx::UInt32, 
                best_irt::Float32, 
                mean_irt::Union{Missing, Float32}, 
                var_irt::Union{Missing, Float32}, 
                n::Union{Missing, UInt16},
                mz::Float32}},
    is_decoy::AbstractVector{Bool};
    q_value_threshold::Float32 = 0.01f0)
    N = length(prec_to_irt)
    precursor_idx = zeros(UInt32, N)
    ms_file_idx = zeros(UInt32, N)
    scan_idx = zeros(UInt32, N)
    best_prob = zeros(Float32, N)
    for (i, (pid, value)) in enumerate(pairs(prec_to_irt))
        precursor_idx[i] = pid
        ms_file_idx[i] = value[:best_ms_file_idx]
        scan_idx[i] = value[:best_scan_idx]
        best_prob[i] = value[:best_prob]
    end
    targets = [is_decoy[pid]==false for pid in precursor_idx]
    target_count=0
    qvalue = zeros(Float32, N)
    for i in range(1, N)
        target_count += targets[i]
        qvalue[i] = Float32((i-target_count)/i)
    end
    best_psms = DataFrame(Dict(
        :precursor_idx => precursor_idx,
        :ms_file_idx => ms_file_idx,
        :scan_idx => scan_idx,
        :best_prob => best_prob,
        :target => targets,
        :q_value => qvalue
        ));
    filter(x->x.q_value<=q_value_threshold,best_psms)
end

function getHuberLossParam(
    huber_δs::Vector{Float32},
    prec_to_irt::Dictionary{UInt32, 
    @NamedTuple{best_prob::Float32, 
                best_ms_file_idx::UInt32, 
                best_scan_idx::UInt32, 
                best_irt::Float32, 
                mean_irt::Union{Missing, Float32}, 
                var_irt::Union{Missing, Float32}, 
                n::Union{Missing, UInt16},
                mz::Float32}},
    is_decoy::AbstractVector{Bool})
    best_psms =  getPsmsForHuberEstimation(prec_to_irt,is_decoy)
    gbpsms = groupby(best_psms,:ms_file_idx)
    return getHuberLossParam(huber_δs,gbpsms)
end
