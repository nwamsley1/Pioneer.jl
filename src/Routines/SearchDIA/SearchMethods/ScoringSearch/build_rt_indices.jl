"""
    build_rt_indices!(search_context, valid_file_indices, passing_refs; min_prob=0.5f0)

Build per-file RT index Arrow files for IntegrateChromatogramsSearch.

Only includes precursors that passed final LightGBM 1% FDR in each specific file,
using their observed iRT values. No library iRT fallback needed.

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
    lib_mzs = getMz(precursors)

    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
    rt_indices_folder = joinpath(temp_folder, "rt_indices")
    mkpath(rt_indices_folder)

    for (file_idx, ref) in zip(valid_file_indices, passing_refs)
        t_file_start = time()

        # Read passing PSMs — these ARE the precursors for this file's RT index
        tbl = Arrow.Table(file_path(ref))
        pids = tbl[:precursor_idx]
        irts = tbl[:irt_obs]
        n_precs = length(pids)

        # Build RT index from passing precursors only
        out_irts = Vector{Float32}(undef, n_precs)
        out_mzs = Vector{Float32}(undef, n_precs)
        out_pids = Vector{UInt32}(undef, n_precs)

        for i in 1:n_precs
            pid = UInt32(pids[i])
            out_pids[i] = pid
            out_mzs[i] = lib_mzs[pid]
            out_irts[i] = Float32(irts[i])
        end

        rt_df = DataFrame(irt = out_irts, prec_mz = out_mzs, precursor_idx = out_pids)
        sort!(rt_df, :irt)

        file_name = getFileIdToName(getMSData(search_context), file_idx)
        rt_index_path = joinpath(rt_indices_folder, "$(file_name)_rt_index.arrow")
        Arrow.write(rt_index_path, rt_df)
        setRtIndex!(getMSData(search_context), file_idx, rt_index_path)

        elapsed = round(time() - t_file_start, digits=3)
        @info "  RT index $(file_name): $(n_precs) precursors (all observed iRT) in $(elapsed)s"
    end

    total_elapsed = round(time() - t_total_start, digits=2)
    @info "Built RT indices for $(length(valid_file_indices)) files in $(total_elapsed)s"
end
