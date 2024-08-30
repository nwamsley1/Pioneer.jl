function huberLossSearch(
    frag_err_dist_dict,
    rt_index_paths,
    prec_set,
    scan_idxs,
    δs,
    bin_rt_size,
    RT_iRT,
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
                                    (3, 0),
                                    params[:quant_search_params]["n_frag_isotopes"],
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
    rt_index = buildRTIndex(rt_df,
                            bin_rt_size = bin_rt_size)
    rt_irt = RT_iRT[parsed_fname]
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
            scored_psms = complex_scored_PSMs,
            unscored_psms = complex_unscored_PSMs,
            spectral_scores = complex_spectral_scores,
            precursor_weights = precursor_weights,
            quad_transmission_func = QuadTransmission(1.0f0, 1000.0f0)
            )...);
    psms[!,:ms_file_idx] .= UInt32(ms_file_idx)
    return psms
end

function getHuberLossParam(
    huber_δs::Vector{Float32})
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
        huber_δs::AbstractVector{Float32}
        )
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
    combdf = combine(gpsms) do sdf
        processHuberLossCurve(sdf[!,:weight],sdf[!,:huber_δ])
    end
    filter!(x->x.n==length(huber_δs), combdf)
    filter!(x->x.wdiff>0.1, combdf)
    filter!(x->!ismissing(x.huber50), combdf)

    #histogram_plot = histogram(log2.(combdf[!,:huber50]))
    #edges = histogram_plot.series_list[1][:x]
    #counts = histogram_plot.series_list[1][:y]
    return combdf#edges, counts, exp2(mean(log2.(combdf[!,:huber50]))), histogram_plot

end

nan_bins = (isnan.(bin_counts).==false).*(bin_counts.>0.0f0)
bin_counts = bin_counts[nan_bins]
bin_edges = bin_edges[nan_bins]

gbpsms = groupby(best_psms,:ms_file_idx)
huber_δs = Float32[300*(1.5^i) for i in range(1, 15)];#Float32[200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400]
@time combdf_test3 = getHuberLossParam(huber_δs);
combdf_test3[!,:huber50] = ceil.(Int, combdf_test3[!,:huber50]);
huber_hist = sort(value_counts(combdf_test3,:huber50),:huber50);
huber_hist[!,:prob] = huber_hist[!,:nrow]./sum(huber_hist[!,:nrow]);
huber_hist[!,:cum_prob] = cumsum(huber_hist[!,:prob]);
hist_a = copy(huber_hist)


gbpsms = groupby(best_psms,:ms_file_idx)
huber_δs = Float32[300*(2^i) for i in range(1, 10)]#Float32[200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400]
combdf_test3 = getHuberLossParam(huber_δs)
combdf_test3[!,:huber50] = ceil.(Int, combdf_test3[!,:huber50])
huber_hist = sort(value_counts(combdf_test3,:huber50),:huber50)
huber_hist[!,:prob] = huber_hist[!,:nrow]./sum(huber_hist[!,:nrow])
huber_hist[!,:cum_prob] = cumsum(huber_hist[!,:prob])
hist_b = copy(huber_hist)


best_psms = DataFrame(Dict(
    :precursor_idx => collect(keys(prec_to_irt)),
    :ms_file_idx => [v[:best_ms_file_idx] for v in values(prec_to_irt)],
    :scan_idx => [v[:best_scan_idx] for v in values(prec_to_irt)],
    :best_prob => [v[:best_prob] for v in values(prec_to_irt)]
    ));
    best_psms[!,:target] = [precursors[:is_decoy][pid]==false for pid in best_psms[!,:precursor_idx]];
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
    targets = 0
    for i in range(1, size(best_psms, 1))
        targets += best_psms[i,:target]
        best_psms[i,:q_value] = Float32((i - targets)/targets)
    end
    filter!(x->x.q_value<=0.01, best_psms)