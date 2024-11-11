
using Profile 
using PProf 
Profile.clear()
@profile begin
peak_fwhms, psms_paths = firstSearch(
    first_search_psms_folder,
    RT_to_iRT_map_dict,
    frag_err_dist_dict,
    irt_errs,
    quad_model_dict,
    file_id_to_parsed_name,
    MS_TABLE_PATHS,
    params_,
    spec_lib,
    ionMatches,
    ionMisses,
    all_fmatches,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    scored_PSMs,
    unscored_PSMs,
    spectral_scores,
    precs
);
end
pprof()

volumes = []
for ttask in test_tasks
    volume = 0.0
    for scan_id in ttask[2]
        if !iszero(scan_id)
            volume += ttable[:isolationWidthMz][scan_id]*length(ttable[:mz_array][scan_id])
        end
    end
    push!(volumes, volume)
end

lib_fragments = getFragments(spec_lib["f_det"])

frag_intensities = Float16[]
for prec in ProgressBar(first_psms[!,:precursor_idx])
    append!(frag_intensities, [x.intensity for x in lib_fragments[getPrecFragRange(spec_lib["f_det"], prec)] if (x.rank < 20)])
end


tscored_psms = DataFrame(Arrow.Table("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SCIEX_PXD050030/RESULTS/temp/first_search_psms/1fullscored.arrow"))