"""
    build_rt_indices!(search_context, valid_file_indices, passing_refs; min_prob=0.5f0)

Build per-file RT index Arrow files for IntegrateChromatogramsSearch.

Only includes precursors that passed final LightGBM 1% FDR in each specific file,
using their observed **empirical RT** values (not iRT). This allows chromatogram
integration to use RT-space tolerances directly, avoiding non-linear iRT artifacts.

Columns written: :irt (Float32, actually empirical RT), :prec_mz (Float32),
:precursor_idx (UInt32), sorted by :irt.
"""
function build_rt_indices!(
    search_context::SearchContext,
    valid_file_indices::Vector{Int},
    passing_refs;
    min_prob::Float32 = 0.5f0
)
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
        rts = tbl[:rt]
        n_precs = length(pids)

        # Build RT index from passing precursors using empirical RT (not iRT)
        out_rts = Vector{Float32}(undef, n_precs)
        out_mzs = Vector{Float32}(undef, n_precs)
        out_pids = Vector{UInt32}(undef, n_precs)

        for i in 1:n_precs
            pid = UInt32(pids[i])
            out_pids[i] = pid
            out_mzs[i] = lib_mzs[pid]
            out_rts[i] = Float32(rts[i])
        end

        # Column named :irt for compatibility with buildRtIndex, but values are empirical RT
        rt_df = DataFrame(irt = out_rts, prec_mz = out_mzs, precursor_idx = out_pids)
        sort!(rt_df, :irt)

        file_name = getFileIdToName(getMSData(search_context), file_idx)
        rt_index_path = joinpath(rt_indices_folder, "$(file_name)_rt_index.arrow")
        Arrow.write(rt_index_path, rt_df)
        setRtIndex!(getMSData(search_context), file_idx, rt_index_path)

        elapsed = round(time() - t_file_start, digits=3)
        @debug_l1 "  RT index $(file_name): $(n_precs) precursors (empirical RT) in $(elapsed)s"
    end
end
