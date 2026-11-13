# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
    build_protein_global_qval_dict(global_pg_score_dict)
    → Dict{Tuple{K, Bool, UInt8}, Float32}   (K: the pg_id, or a name)

Compute protein global q-values directly from score dictionary.
"""
function build_protein_global_qval_dict(
    global_pg_score_dict::Dict{Tuple{K, Bool, UInt8}, Float32}
) where {K}
    n = length(global_pg_score_dict)
    keys_vec = collect(keys(global_pg_score_dict))
    scores = Float32[global_pg_score_dict[k] for k in keys_vec]
    targets = Bool[k[2] for k in keys_vec]

    qvals = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals)

    qval_dict = Dict{Tuple{K, Bool, UInt8}, Float32}()
    sizehint!(qval_dict, n)
    for i in 1:n
        qval_dict[keys_vec[i]] = qvals[i]
    end
    return qval_dict
end

"""
    update_psms_with_protein_scores_refs(paired_refs::Vector{PairedSearchFiles},
                                         pg_name_to_global_pg_score::Dict{PGKey,Float32},
                                         pg_score_to_qval,
                                         global_pg_score_to_qval_dict::Dict{Tuple{UInt32, Bool, UInt8}, Float32})

Update PSMs with model-scored `pg_score` values and q-values using references.
"""
function update_psms_with_protein_scores_refs(
    paired_refs::Vector{PairedSearchFiles},
    pg_name_to_global_pg_score::Dict{PGKey,Float32},
    pg_score_to_qval,
    global_pg_score_to_qval_dict::Dict{Tuple{UInt32, Bool, UInt8}, Float32}
)
    for paired_ref in paired_refs
        psm_ref = paired_ref.psm_ref
        pg_ref = paired_ref.protein_ref

        if !exists(psm_ref) || !exists(pg_ref)
            @user_warn "Missing file" psm_path = file_path(paired_ref.psm_ref) pg_path = file_path(paired_ref.protein_ref)
            continue
        end

        pg_score_lookup = Dict{PGKey, Tuple{Float32, Float32}}()
        # Only scalars are kept, so the protein-group file is unmapped afterwards.
        with_arrow_table(file_path(pg_ref)) do pg_table
            _fill_pg_score_lookup!(pg_score_lookup, pg_table[:protein_name], pg_table[:target],
                                   pg_table[:entrap_id], pg_table[:pg_score], pg_table[:pg_pep])
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
        protein_scores = Vector{Union{Missing, Float32}}(missing, n_psms)
        global_pg_scores = Vector{Union{Missing, Float32}}(missing, n_psms)
        pg_qvals = Vector{Union{Missing, Float32}}(missing, n_psms)
        global_pg_qvals = Vector{Union{Missing, Float32}}(missing, n_psms)
        pg_peps = Vector{Union{Missing, Float32}}(missing, n_psms)

        _fill_psm_pg_scores!(protein_scores, global_pg_scores, pg_qvals, global_pg_qvals, pg_peps,
            needed[!, :inferred_protein_group], needed[!, :use_for_protein_quant],
            needed[!, :target], needed[!, :entrapment_group_id],
            pg_score_lookup, pg_name_to_global_pg_score, pg_score_to_qval, global_pg_score_to_qval_dict)

        add_columns_via_sidecar!(psm_ref,
            :pg_score        => protein_scores,
            :global_pg_score => global_pg_scores,
            :pg_qval         => pg_qvals,
            :global_pg_qval  => global_pg_qvals,
            :pg_pep          => pg_peps;
            tag = "pg_scores")
    end
end

# Column kernels for the protein-score annotation (function barriers: the columns are untyped).
function _fill_pg_score_lookup!(pg_score_lookup, names, targets, entrap_ids, scores, peps)
    for i in eachindex(names)
        pg_score_lookup[ProteinKey(names[i], targets[i], entrap_ids[i])] = (scores[i], peps[i])
    end
    return pg_score_lookup
end

function _fill_psm_pg_scores!(protein_scores, global_pg_scores, pg_qvals, global_pg_qvals, pg_peps,
                              ipg_col, uq_col, tgt_col, ent_col,
                              pg_score_lookup, pg_name_to_global_pg_score, pg_score_to_qval,
                              global_pg_score_to_qval_dict)
    for i in eachindex(ipg_col)
        ismissing(ipg_col[i]) && continue
        uq_col[i] == false && continue
        key = ProteinKey(ipg_col[i], tgt_col[i], ent_col[i])
        haskey(pg_score_lookup, key) || continue
        scores_tuple = pg_score_lookup[key]
        protein_scores[i] = scores_tuple[1]
        pg_peps[i] = scores_tuple[2]
        haskey(pg_name_to_global_pg_score, key) ||
            throw("Missing global pg score lookup key!!!")
        global_pg_scores[i] = pg_name_to_global_pg_score[key]
        pg_qvals[i] = pg_score_to_qval(protein_scores[i])
        dict_key = (key.name, key.is_target, key.entrap_id)
        global_pg_qvals[i] = get(global_pg_score_to_qval_dict, dict_key, missing)
    end
    return nothing
end
