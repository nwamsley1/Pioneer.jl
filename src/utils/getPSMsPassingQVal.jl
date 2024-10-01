
function find_score_threshold(
    scores::Vector{@NamedTuple{score::Float32, target::Bool}},
    q_value_threshold::Float32
)
    # Second pass: Find the probability threshold
    targets_above = 0
    decoys_above = 0

    for (score, target) in reverse(scores)

        targets_above += target
        decoys_above += (1 - target)
        
        current_q_value = decoys_above / (targets_above + decoys_above)
        
        if current_q_value > q_value_threshold
           return score#return prob  # This is the probability threshold we're looking for
        end
    end

    return scores[end][:score]

end

function mergeSortedPSMScores(
                    input_dir::String, 
                    output_path::String,
                    ; N = 10000000)
    
    function fillColumn!(
        peptide_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    function addPrecursorToHeap!(
        psms_heap::BinaryMaxHeap{Tuple{Float32, Int64}},
        sort_key::AbstractVector{Float32},
        table_idx::Int64,
        row_idx::Int64)
        push!(
            psms_heap,
            (
            sort_key[row_idx],
            table_idx
            )
        )
    end

    input_paths = [path for path in readdir(input_dir,join=true) if endswith(path,".arrow")]
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))
    psms_batch = DataFrame()
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    psms_heap = BinaryMaxHeap{Tuple{Float32, Int64}}()
    psms_batch[!,:prob] = zeros(Float32, N)
    psms_batch[!,:target] = zeros(Bool, N)


    #Load psms_heap
    for (i, table) in enumerate(tables)
        addPrecursorToHeap!(
            psms_heap,
            table[:prob],
            i,
            1
        )
    end
    i = 1
    while length(psms_heap)>0
        _, table_idx = pop!(psms_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addPrecursorToHeap!(
            psms_heap,
            table[:prob],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in [:prob,:target]
                fillColumn!(
                    psms_batch[!,col],
                    col,
                    sorted_tuples,
                    tables,
                    i
                )
            end
            Arrow.append(
                output_path,
                psms_batch
            )
            i = 1
        end
    end
    for col in [:prob,:target]
        fillColumn!(
            psms_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        psms_batch[range(1, max(1, i - 1)),:]
    )
    return nothing
end

function getPEPSpline(
                            merged_psms_path::String,
                            score_col::Symbol;
                            min_pep_points_per_bin = 5000,
                            n_spline_bins = 20
                            )

    psms_scores = Arrow.Table(merged_psms_path)
    Q = length(psms_scores[1])
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_target_fraction, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets = 0.0f0, 0
    for i in range(1, Q)
        bin_size += 1
        targets += psms_scores[:target][i]
        mean_prob += psms_scores[score_col][i]
        if bin_size == min_pep_points_per_bin
            bin_idx += 1
            bin_target_fraction[bin_idx] = targets/bin_size
            bin_mean_prob[bin_idx] = mean_prob/bin_size
            bin_size, targets, mean_prob = zero(Int64), zero(Int64), zero(Float32)
        end
    end
    bin_target_fraction[end] = targets/bin_size
    bin_mean_prob[end] = mean_prob/bin_size
    #n_spline_bins = ceil(Int, length(bin_mean_prob)/15)
    #UniformSpline(bin_target_fraction, bin_mean_prob, 3, 20)
    return UniformSpline(bin_target_fraction, bin_mean_prob, 3, n_spline_bins)
end

function getQValueSpline(
                            merged_psms_path::String,
                            score_col::Symbol;
                            min_pep_points_per_bin = 1000
                            )

    psms_scores = Arrow.Table(merged_psms_path)
    Q = length(psms_scores[1])
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_qval, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets, decoys = 0.0f0, 0, 0
    targets, decoys = 0, 0
    for i in range(1, Q)
        targets += psms_scores[:target][i]
        decoys += (1 - psms_scores[:target][i])
    end

    min_q_val = typemax(Float32)
    for i in reverse(range(1, Q))
        bin_size += 1
        targets -= psms_scores[:target][i]
        decoys -= (1 - psms_scores[:target][i])
        mean_prob += psms_scores[score_col][i]
        if bin_size == min_pep_points_per_bin
            bin_idx += 1
            qval = decoys/(targets + decoys)
            if qval > min_q_val
                bin_qval[bin_idx] = min_q_val
            else
                min_q_val = qval
                bin_qval[bin_idx] = qval
            end
            bin_mean_prob[bin_idx] = mean_prob/bin_size
            bin_size, mean_prob = zero(Int64), zero(Float32)
        end
    end
    bin_qval[end] = targets/bin_size
    bin_mean_prob[end] = mean_prob/bin_size
    prepend!(bin_qval, 1.0f0)
    prepend!(bin_mean_prob, 0.0f0)
    bin_qval = bin_qval[isnan.(bin_mean_prob).==false]
    bin_mean_prob = bin_mean_prob[isnan.(bin_mean_prob).==false]
    #return bin_qval, bin_mean_prob
    return linear_interpolation(
        Interpolations.deduplicate_knots!(bin_mean_prob, move_knots=true),
        bin_qval, extrapolation_bc=Line())
end

function getPSMsPassingQVal(
                            quant_psms_folder::String, 
                            passing_psms_folder::String,
                            pep_spline::UniformSpline,
                            qval_interp::Interpolations.Extrapolation,
                            q_val_threshold::Float32,
                            )
    
    file_paths = [fpath for fpath in readdir(quant_psms_folder, join=true) if endswith(fpath,".arrow")]
    for file_path in file_paths
        # Read the Arrow table
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        passing_psms[!,:qval] = qval_interp.(passing_psms[!,:prob])
        passing_psms[!,:pep] = pep_spline.(passing_psms[!,:prob])
        # Sample the rows and convert to DataFrame
        select!(passing_psms,
        [
            :precursor_idx,
            :prob,
            :qval,
            :pep,
            :weight,
            :target,
            :irt_obs,
            :missed_cleavage,
            :isotopes_captured,
            :scan_idx,
            :ms_file_idx])
        filter!(x->x.qval<=q_val_threshold, passing_psms)
        # Append to the result DataFrame
        Arrow.write(
            joinpath(passing_psms_folder, basename(file_path)),
            passing_psms
        )
    end

    return
end