# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    build_protein_global_score_dicts(pg_refs, sqrt_n_runs, n_proteins)
    -> (global_pg_score_dict, pg_name_to_global_pg_score)

Stream per-file protein group files reading only `(protein_name, target, entrap_id, pg_score)`.
Returns composite-key dictionaries for protein global score computation while
retaining only the top `sqrt_n_runs` scores per protein.
"""
function _insert_top_protein_score!(
    scores::Vector{Float32},
    score::Float32,
    max_scores::Int
)
    max_scores <= 0 && return scores

    if (length(scores) == max_scores) && (score <= scores[end])
        return scores
    end

    insert_idx = length(scores) + 1
    @inbounds for i in eachindex(scores)
        if score > scores[i]
            insert_idx = i
            break
        end
    end

    insert!(scores, insert_idx, score)
    if length(scores) > max_scores
        pop!(scores)
    end

    return scores
end

function build_protein_global_score_dicts(
    pg_refs::Vector{ProteinGroupFileReference},
    sqrt_n_runs::Int,
    n_proteins::Int
)
    score_acc = Dict{Tuple{String,Bool,UInt8}, Vector{Float32}}()
    sizehint!(score_acc, n_proteins)

    for ref in pg_refs
        tbl = Arrow.Table(file_path(ref))
        n = length(tbl.protein_name)
        n == 0 && continue
        @inbounds for i in 1:n
            key = (tbl.protein_name[i], tbl.target[i], tbl.entrap_id[i])
            if !haskey(score_acc, key)
                score_acc[key] = sizehint!(Float32[], sqrt_n_runs)
            end
            _insert_top_protein_score!(score_acc[key], tbl.pg_score[i], sqrt_n_runs)
        end
    end

    n_observed = length(score_acc)
    global_pg_score_dict = Dict{Tuple{String,Bool,UInt8}, Float32}()
    sizehint!(global_pg_score_dict, n_observed)
    pg_name_to_global_pg_score = Dict{ProteinKey, Float32}()
    sizehint!(pg_name_to_global_pg_score, n_observed)

    for (key, scores) in score_acc
        gs = logodds(scores, sqrt_n_runs)
        global_pg_score_dict[key] = gs
        pg_name_to_global_pg_score[ProteinKey(key[1], key[2], key[3])] = gs
    end

    return global_pg_score_dict, pg_name_to_global_pg_score
end

"""
    build_protein_global_qval_dict(global_pg_score_dict)
    → Dict{Tuple{String,Bool,UInt8}, Float32}

Compute protein global q-values directly from score dictionary.
"""
function build_protein_global_qval_dict(
    global_pg_score_dict::Dict{Tuple{String,Bool,UInt8}, Float32}
)
    n = length(global_pg_score_dict)
    keys_vec = collect(keys(global_pg_score_dict))
    scores = Float32[global_pg_score_dict[k] for k in keys_vec]
    targets = Bool[k[2] for k in keys_vec]

    perm = sortperm(collect(zip(scores, targets)); by = x -> (-x[1], -x[2]))
    permute!(keys_vec, perm)
    permute!(scores, perm)
    permute!(targets, perm)

    qvals = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals)

    qval_dict = Dict{Tuple{String,Bool,UInt8}, Float32}()
    sizehint!(qval_dict, n)
    for i in 1:n
        qval_dict[keys_vec[i]] = qvals[i]
    end
    return qval_dict
end

"""
    update_psms_with_probit_scores_refs(paired_refs::Vector{PairedSearchFiles},
                                       pg_name_to_global_pg_score::Dict{ProteinKey,Float32},
                                       pg_score_to_qval::Interpolations.Extrapolation,
                                       global_pg_score_to_qval_dict::Dict{Tuple{String,Bool,UInt8}, Float32})

Update PSMs with probit-scored pg_score values and q-values using references.
"""
function update_psms_with_probit_scores_refs(
    paired_refs::Vector{PairedSearchFiles},
    pg_name_to_global_pg_score::Dict{ProteinKey,Float32},
    pg_score_to_qval::Interpolations.Extrapolation,
    global_pg_score_to_qval_dict::Dict{Tuple{String,Bool,UInt8}, Float32}
)
    total_psms_updated = 0
    files_processed = 0

    for paired_ref in paired_refs
        psm_ref = paired_ref.psm_ref
        pg_ref = paired_ref.protein_ref

        if !exists(psm_ref) || !exists(pg_ref)
            @user_warn "Missing file" psm_path = file_path(paired_ref.psm_ref) pg_path = file_path(paired_ref.protein_ref)
            continue
        end

        pg_table = Arrow.Table(file_path(pg_ref))
        pg_score_lookup = Dict{ProteinKey, Tuple{Float32, Float32}}()
        n_pg_rows = length(pg_table[:protein_name])

        for i in 1:n_pg_rows
            key = ProteinKey(
                pg_table[:protein_name][i],
                pg_table[:target][i],
                pg_table[:entrap_id][i]
            )
            pep_val = pg_table[:pg_pep][i]
            pg_score_lookup[key] = (pg_table[:pg_score][i], pep_val)
        end

        # Compute the 5 new score columns without rewriting the whole PSM file.
        # We read only the (small) columns we need to derive them, then write
        # the 5 new columns to a sidecar registered on psm_ref. The full file
        # rewrite happens later at the MaxLFQ-boundary sort, which is the
        # natural consolidation point for main + sidecars.
        needed = materialize_columns(psm_ref,
            [:inferred_protein_group, :use_for_protein_quant,
             :target, :entrapment_group_id])
        n_psms = nrow(needed)
        probit_pg_scores = Vector{Union{Missing, Float32}}(missing, n_psms)
        global_pg_scores = Vector{Union{Missing, Float32}}(missing, n_psms)
        pg_qvals = Vector{Union{Missing, Float32}}(missing, n_psms)
        global_pg_qvals = Vector{Union{Missing, Float32}}(missing, n_psms)
        pg_peps = Vector{Union{Missing, Float32}}(missing, n_psms)

        ipg_col = needed[!, :inferred_protein_group]
        uq_col  = needed[!, :use_for_protein_quant]
        tgt_col = needed[!, :target]
        ent_col = needed[!, :entrapment_group_id]

        for i in 1:n_psms
            ismissing(ipg_col[i]) && continue
            uq_col[i] == false && continue
            key = ProteinKey(ipg_col[i], tgt_col[i], ent_col[i])
            haskey(pg_score_lookup, key) || continue
            scores_tuple = pg_score_lookup[key]
            probit_pg_scores[i] = scores_tuple[1]
            pg_peps[i] = scores_tuple[2]
            haskey(pg_name_to_global_pg_score, key) ||
                throw("Missing global pg score lookup key!!!")
            global_pg_scores[i] = pg_name_to_global_pg_score[key]
            pg_qvals[i] = pg_score_to_qval(probit_pg_scores[i])
            dict_key = (key.name, key.is_target, key.entrap_id)
            global_pg_qvals[i] = get(global_pg_score_to_qval_dict, dict_key, missing)
        end

        add_columns_via_sidecar!(psm_ref,
            :pg_score        => probit_pg_scores,
            :global_pg_score => global_pg_scores,
            :pg_qval         => pg_qvals,
            :qlobal_pg_qval  => global_pg_qvals,
            :pg_pep          => pg_peps;
            tag = "pg_scores")
        total_psms_updated += n_psms

        files_processed += 1
    end
end
"""
    add_pep_column(new_col::Symbol, score_col::Symbol, target_col::Symbol;
                   doSort::Bool = false, fdr_scale_factor::Float32 = 1.0f0)

Add a column of posterior error probabilities using `get_PEP!`.
The input DataFrame must already be sorted by `score_col` if `doSort` is false.
"""
function add_pep_column(new_col::Symbol, score_col::Symbol, target_col::Symbol;
                        doSort::Bool=false, fdr_scale_factor::Float32=1.0f0)
    desc = "add_pep_column($new_col from $score_col)"
    op = function(df)
        scores = df[!, score_col]::AbstractVector{Float32}
        targets = df[!, target_col]::AbstractVector{Bool}
        pep_vals = Vector{Float32}(undef, length(scores))
        get_PEP!(scores, targets, pep_vals; doSort=doSort, fdr_scale_factor=fdr_scale_factor)
        df[!, new_col] = pep_vals
        return df
    end
    return desc => op
end

