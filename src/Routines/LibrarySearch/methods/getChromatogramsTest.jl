precursors_passing

passing_rt_index = Dict{String, Any}()
for (file_path, rt_index) in RT_INDICES
    passing_rt_index[file_path] = retentionTimeIndex(Float32, UInt32)
    for rt_bin in rt_index.rt_bins
        push!(ppassing_rt_index[file_path].rt_bins,
        rtIndexBin(
            rt_bin.lb,
            rt_bin.ub,
            filter(x->first(x)∈precursors_passing)
        )
        )
end

sort!(best_psms, :RT)

best_psms_passing = copy(best_psms[(best_psms[!,:target]).&(best_psms[!,:q_value].<=0.01).&(best_psms[!,:best_trace]),:])
subset_rt_index = buildRTIndex(best_psms_passing[!,:irt_obs], 
             best_psms_passing[!,:prec_mz],
             best_psms_passing[!,:precursor_idx],
             0.25)

function getChromatograms(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    params::Any;
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
            return getChromatograms(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                precursors,
                                kwargs[:fragment_lookup_table], 
                                kwargs[:ms_file_idx],
                                kwargs[:rt_to_irt_spline],
                                kwargs[:mass_err_model],
                                Float32(params[:deconvolution_params]["huber_delta"]),
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
                                kwargs[:chromatograms][thread_id],
                                kwargs[:scored_psms][thread_id],
                                kwargs[:unscored_psms][thread_id],
                                kwargs[:spectral_scores][thread_id],
                                kwargs[:precursor_weights][thread_id],
                                (3, 0),
                                Int64(params[:quant_search_params]["min_frag_count"]),
                                Float32(params[:quant_search_params]["min_spectral_contrast"]),
                                Float32(params[:quant_search_params]["min_log2_matched_ratio"]),
                                Tuple([Int64(x) for x in params[:quant_search_params]["min_topn_of_m"]]),
                                Int64(params[:quant_search_params]["max_best_rank"]),
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


@time chroms = vcat(getChromatograms(
    MS_TABLE, 
    params_;
    precursors = prosit_lib["precursors"],
    fragment_lookup_table = library_fragment_lookup_table,
    rt_index = subset_rt_index,
    ms_file_idx = UInt32(ms_file_idx), 
    rt_to_irt_spline = RT_iRT[file_id_to_parsed_name[ms_file_idx]],
    mass_err_model = MassErrorModel{Float32}(3.195095f0, (15.593506f0*1.2f0, 7.698679f0*1.2f0)),#frag_err_dist_dict[ms_file_idx],
    irt_err = irt_err,#irt_errs[ms_file_idx]/3,
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
    )...);
getIsotopesCaptured!(chroms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
  

#=
uncorrected_chroms = copy(chroms)
uncorrected_gchroms = groupby(uncorrected_chroms, [:precursor_idx])
=#
subdf = uncorrected_gchroms[(precursor_idx = precs_passing[N],)]
gsubdf = groupby(subdf,:isotopes_captured)
p = plot(layout = (1, 2))
best_rt = subdf[argmax(subdf[!,:intensity]),:rt]
for (key, chrom) in pairs(gsubdf)
    plot_range = (chrom[!,:rt].>(best_rt - 0.1)).&(chrom[!,:rt].<(best_rt + 0.1))
    plot!(chrom[plot_range,:rt], 
            chrom[plot_range,:intensity], 
            subplot = 1,
            title = "Uncorrected \n"*precursors[:sequence][precs_passing[N]],
            seriestype=:scatter, 
            label = key[:isotopes_captured], 
            show = true)
end

gchroms = groupby(chroms, [:precursor_idx])
subdf = gchroms[(precursor_idx = precs_passing[N],)]
gsubdf = groupby(subdf,:isotopes_captured)
#plot()
for (key, chrom) in pairs(gsubdf)
    plot_range = (chrom[!,:rt].>(best_rt - 0.1)).&(chrom[!,:rt].<(best_rt + 0.1))
    plot!(chrom[plot_range,:rt], 
            chrom[plot_range,:intensity],             
            subplot = 2, title = "Corrected \n"*precursors[:sequence][precs_passing[N]], seriestype=:scatter, 
            label = key[:isotopes_captured], show = true)
end

N += 1


filter!(x->first(x.isotopes_captured)<2, chroms)
filter!(x->first(x.isotopes_captured)>-1, chroms)
correctPrecursorAbundances!(chroms[!,:intensity],
                            iso_splines,
                            chroms[!,:isotopes_captured],
                            chroms[!,:precursor_idx],
                            precursors[:mz],
                            precursors[:prec_charge],
                            precursors[:sulfur_count])
