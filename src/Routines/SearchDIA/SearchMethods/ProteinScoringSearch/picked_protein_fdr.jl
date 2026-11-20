# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

#==========================================================
Experimental picked protein-group FDR for the experiment-wide protein q-value.

The default experiment-wide protein q-value counts every target and decoy group. Picked
competition keeps only the better of each target/decoy pair before counting, so the decoy
of a protein that is present (which its target beats) no longer uses up the FDR budget.
Target and decoy groups are inferred separately and need not match one to one, so the
pairing rule matters. ENV["PIONEER_PROTEIN_PICK"] selects it (default off):

  off        all groups counted (current behaviour)
  name       a group competes with the opposite-label group of the same name (pg_id);
             a group with no such partner competes alone
  accession  a group survives only if its score beats every opposite-label group that
             shares an accession with it (ties go to the decoy)
  mirror     Osprey-style: units are target groups (plus accessions in no target group);
             each side's unit score is the best experiment-wide precursor score over the
             unit's accessions, target and decoy precursors paired through the library's
             shared accessions; the better side wins and q-values are computed on that score

Competition is always within one entrapment group, so A_trap;B_trap competes only with its
own decoys. q = decoys / targets among winners (get_qvalues!); a loser gets q = 1.
==========================================================#

const PROTEIN_PICK_MODES = (:off, :name, :accession, :mirror)

function protein_pick_mode()
    mode = Symbol(lowercase(strip(get(ENV, "PIONEER_PROTEIN_PICK", "off"))))
    mode in PROTEIN_PICK_MODES || throw(ArgumentError(
        "PIONEER_PROTEIN_PICK must be one of $(join(PROTEIN_PICK_MODES, ", ")); got $mode"))
    return mode
end

# q-values for the winners; every other key gets 1.
function _picked_qvals(keys_vec::Vector{K}, scores::Vector{Float32}, targets::Vector{Bool},
                       winner::BitVector) where {K}
    idx = findall(winner)
    q = ones(Float32, length(keys_vec))
    if !isempty(idx)
        qw = Vector{Float32}(undef, length(idx))
        get_qvalues!(scores[idx], targets[idx], qw)
        q[idx] .= qw
    end
    return Dict{K, Float32}(keys_vec[i] => q[i] for i in eachindex(keys_vec))
end

_accessions(name::AbstractString) = split(name, ';'; keepempty = false)

function _log_pick(mode, keys_vec, targets, winner, qd)
    nt = count(targets); nd = length(targets) - nt
    wt = count(winner .& targets); wd = count(winner .& .!targets)
    pass = count(i -> targets[i] && qd[keys_vec[i]] <= 0.01f0, eachindex(keys_vec))
    @user_info "Picked protein FDR ($mode): groups target=$nt decoy=$nd; winners target=$wt decoy=$wd; " *
               "targets at experiment-wide q<=1%: $pass"
end

"""
    build_protein_global_qval_dict_picked(global_pg_score_dict, mode; group_names, precursors,
                                          precursor_global_prob)

Experiment-wide protein q-values with picked target/decoy competition (see file header).
"""
function build_protein_global_qval_dict_picked(
    global_pg_score_dict::Dict{Tuple{K, Bool, UInt8}, Float32},
    mode::Symbol;
    group_names::Vector{String},
    precursors = nothing,
    precursor_global_prob::Dict{UInt32, Float32} = Dict{UInt32, Float32}(),
) where {K}
    keys_vec = collect(keys(global_pg_score_dict))
    scores = Float32[global_pg_score_dict[k] for k in keys_vec]
    targets = Bool[k[2] for k in keys_vec]
    winner = falses(length(keys_vec))

    if mode == :name
        n_unpaired = 0
        for (i, k) in enumerate(keys_vec)
            partner = get(global_pg_score_dict, (k[1], !k[2], k[3]), nothing)
            partner === nothing && (n_unpaired += 1)
            # Ties go to the decoy.
            winner[i] = partner === nothing || (k[2] ? scores[i] > partner : scores[i] >= partner)
        end
        @user_info "Picked protein FDR (name): $n_unpaired of $(length(keys_vec)) groups have no same-name partner"
        qd = _picked_qvals(keys_vec, scores, targets, winner)
        _log_pick(mode, keys_vec, targets, winner, qd)
        return qd

    elseif mode == :accession
        # Best score per (accession, target, entrap_id).
        best = Dict{Tuple{SubString{String}, Bool, UInt8}, Float32}()
        accs = Vector{Vector{SubString{String}}}(undef, length(keys_vec))
        for (i, k) in enumerate(keys_vec)
            accs[i] = _accessions(group_names[k[1]])
            for a in accs[i]
                key = (a, k[2], k[3])
                best[key] = max(get(best, key, -Inf32), scores[i])
            end
        end
        n_contested = 0
        for (i, k) in enumerate(keys_vec)
            opp = -Inf32
            for a in accs[i]
                opp = max(opp, get(best, (a, !k[2], k[3]), -Inf32))
            end
            isfinite(opp) && (n_contested += 1)
            winner[i] = k[2] ? scores[i] > opp : scores[i] >= opp
        end
        @user_info "Picked protein FDR (accession): $n_contested of $(length(keys_vec)) groups share an accession with an opposite-label group"
        qd = _picked_qvals(keys_vec, scores, targets, winner)
        _log_pick(mode, keys_vec, targets, winner, qd)
        return qd

    elseif mode == :mirror
        precursors === nothing && error("mirror picked protein FDR needs the library precursors")
        isempty(precursor_global_prob) && error("mirror picked protein FDR needs experiment-wide precursor scores")
        ids = getPrecursorTextIds(precursors)
        is_decoy = getIsDecoy(precursors)
        entrap = getEntrapmentGroupId(precursors)
        # Best experiment-wide precursor score per (accession index, target, entrap_id).
        acc_best = Dict{Tuple{UInt32, Bool, UInt8}, Float32}()
        for (pid, s) in precursor_global_prob
            t = !is_decoy[pid]
            e = entrap[pid]
            for a in ids.accession_set_members[ids.accession_set_id[pid]]
                key = (a, t, e)
                acc_best[key] = max(get(acc_best, key, -Inf32), s)
            end
        end
        acc_index = Dict(name => UInt32(i) for (i, name) in enumerate(ids.accession_names))
        # Units: target groups, then accessions that are in no target group.
        covered = Set{Tuple{UInt32, UInt8}}()
        unit_t = Float32[]; unit_d = Float32[]; unit_group = Int[]
        for (i, k) in enumerate(keys_vec)
            k[2] || continue
            bt = -Inf32; bd = -Inf32
            for a in _accessions(group_names[k[1]])
                ai = get(acc_index, String(a), nothing)
                ai === nothing && continue
                push!(covered, (ai, k[3]))
                bt = max(bt, get(acc_best, (ai, true, k[3]), -Inf32))
                bd = max(bd, get(acc_best, (ai, false, k[3]), -Inf32))
            end
            push!(unit_t, bt); push!(unit_d, bd); push!(unit_group, i)
        end
        for ((a, t, e), _) in acc_best
            (a, e) in covered && continue
            push!(covered, (a, e))
            push!(unit_t, get(acc_best, (a, true, e), -Inf32))
            push!(unit_d, get(acc_best, (a, false, e), -Inf32))
            push!(unit_group, 0)
        end
        # Pick: target wins ties, as in Osprey.
        n = length(unit_t)
        w_score = Vector{Float32}(undef, n); w_target = Vector{Bool}(undef, n)
        for u in 1:n
            tw = unit_t[u] >= unit_d[u]
            w_target[u] = tw
            w_score[u] = tw ? unit_t[u] : unit_d[u]
        end
        ok = findall(isfinite, w_score)
        qu = ones(Float32, n)
        if !isempty(ok)
            qw = Vector{Float32}(undef, length(ok))
            get_qvalues!(w_score[ok], w_target[ok], qw)
            qu[ok] .= qw
        end
        qd = Dict{Tuple{K, Bool, UInt8}, Float32}(k => 1.0f0 for k in keys_vec)
        for u in 1:n
            g = unit_group[u]
            (g == 0 || !w_target[u]) && continue
            qd[keys_vec[g]] = qu[u]
        end
        nt_units = count(!=(0), unit_group)
        pass = count(u -> unit_group[u] != 0 && w_target[u] && qu[u] <= 0.01f0, 1:n)
        @user_info "Picked protein FDR (mirror): units=$n (target groups $nt_units, ungrouped accessions " *
                   "$(n - nt_units)); target winners=$(count(w_target)) decoy winners=$(n - count(w_target)); " *
                   "target groups at experiment-wide q<=1%: $pass"
        return qd
    end
    error("unknown picked protein FDR mode $mode")
end
