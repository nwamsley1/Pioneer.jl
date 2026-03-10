"""
    build_rt_indices!(search_context, valid_file_indices, passing_refs; min_prob=0.5f0)

Build per-file RT index Arrow files for IntegrateChromatogramsSearch.

Includes ALL library precursors. For precursors with passing PSMs,
uses observed iRT. For all others, uses library iRT.

Columns written: :irt (Float32), :prec_mz (Float32), :precursor_idx (UInt32), sorted by :irt.
"""
function build_rt_indices!(
    search_context::SearchContext,
    valid_file_indices::Vector{Int},
    passing_refs;
    min_prob::Float32 = 0.5f0
)
    t_total_start = time()

    precursors = getPrecursors(getSpecLib(search_context))
    n_precursors = length(precursors)
    lib_irts = getIrt(precursors)
    lib_mzs = getMz(precursors)

    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
    rt_indices_folder = joinpath(temp_folder, "rt_indices")
    mkpath(rt_indices_folder)

    for (file_idx, ref) in zip(valid_file_indices, passing_refs)
        t_file_start = time()

        # Load passing PSMs → observed iRT lookup
        tbl = Arrow.Table(file_path(ref))
        obs_irt = Dictionary{UInt32, Float32}()
        if length(tbl[:precursor_idx]) > 0
            pids = tbl[:precursor_idx]
            irts_col = tbl[:irt_obs]
            for i in eachindex(pids)
                insert!(obs_irt, UInt32(pids[i]), Float32(irts_col[i]))
            end
        end
        n_observed = length(obs_irt)

        # Build arrays with ALL library precursors
        out_irts = Vector{Float32}(undef, n_precursors)
        out_mzs = Vector{Float32}(undef, n_precursors)
        out_pids = Vector{UInt32}(undef, n_precursors)

        for i in 1:n_precursors
            pid = UInt32(i)
            out_pids[i] = pid
            out_mzs[i] = lib_mzs[i]
            out_irts[i] = haskey(obs_irt, pid) ? obs_irt[pid] : lib_irts[i]
        end

        rt_df = DataFrame(irt = out_irts, prec_mz = out_mzs, precursor_idx = out_pids)
        sort!(rt_df, :irt)

        file_name = getFileIdToName(getMSData(search_context), file_idx)
        rt_index_path = joinpath(rt_indices_folder, "$(file_name)_rt_index.arrow")
        Arrow.write(rt_index_path, rt_df)
        setRtIndex!(getMSData(search_context), file_idx, rt_index_path)

        elapsed = round(time() - t_file_start, digits=3)
        @info "  RT index $(file_name): $(n_precursors) precursors ($(n_observed) observed iRT, $(n_precursors - n_observed) library fallback) in $(elapsed)s"
    end

    total_elapsed = round(time() - t_total_start, digits=2)
    @info "Built RT indices for $(length(valid_file_indices)) files ($(n_precursors) precursors each) in $(total_elapsed)s"
end
