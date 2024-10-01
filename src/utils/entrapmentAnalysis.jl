function entrapmentAnalysis(
    quant_psms_folder::String,
    base_pep_ids::AbstractVector{Int64},
    entrapment_group_ids::AbstractVector{Int64}
    )
    #Table is sorted by score in descending order
    merged_psms = DataFrame(Tables.columntable(Arrow.Table(readdir(quant_psms_folder,join=true))))
    merged_psms = merged_psms[!,[:precursor_idx,:prob,:target]]
    merged_psms[!,:base_pep_id] = [base_pep_ids[pid] for pid in merged_psms[!,:precursor_idx]]
    merged_psms[!,:entrapment_group_id] = [entrapment_group_ids[pid] for pid in merged_psms[!,:precursor_idx]]
    sort!(merged_psms, :prob, rev = true)
    entrapment_size = 0
    target_size = 0
    for groupid in entrapment_group_ids
        entrapment_size += (groupid > 0)
        target_size += iszero(groupid)
    end
    r = entrapment_size/target_size
    println("r $r")

    target_scores = Vector{Union{Missing, Float32}}(undef, maximum(base_pep_ids))
    for i in range(1, length(target_scores))
        target_scores[i] = missing
    end
    
    combined_fdp = zeros(Float32, size(merged_psms, 1))
    decoy_fdr = zeros(Float32, size(merged_psms, 1))
    Nϵ = 0 #number of entrapment targets at the current score
    Nτ = 0 #number of non-entrapment targets at the current score
    Nd = 0 #number of decoys at the current score
    for i in range(1, size(merged_psms, 1))
        is_target = merged_psms[i,:target]
        if is_target
            if merged_psms[i,:entrapment_group_id] == 0
                Nτ += 1
                #update target score
                #pid = merged_psms[i,:precursor_idx]
                #target_scores[base_pep_ids[pid]] = merged_psms[i,:prob]
            else
                Nϵ += 1
            end
        else
            Nd += 1
        end
        combined_fdp[i] = (Nϵ*(1 + 1/r))/(Nτ + Nϵ)
        decoy_fdr[i] = Nd/(Nτ + Nϵ)
    end
    #=
    paired_fdp = zeros(Float32, size(merged_psms, 1))
    Nϵst = 0
    Nϵts = 0
    for i in reverse(range(1, size(merged_psms, 1)))
        paired_fdp[i] = (Nϵ + Nϵst + 2*Nϵts)/(Nτ + Nϵ)
        is_target = merged_psms[i,:target]
        if is_target
            if iszero(merged_psms[i,:entrapment_group_ids])
                Nϵ += is_entrapment
                pid = merged_psms[i,:precursor_idx]
                target_score = target_scores[base_pep_ids[pid]]
                if ismissing(target_score) #Target 
                    Nϵst += 1
                elseif target_score < merged_psms[i,:prob]
                    Nϵts += 1
                end
            else
                Nτ += (1 - is_entrapment)
            end
        end
    end
    =#
    return DataFrame((decoy_fdr = decoy_fdr, combined_fdp = combined_fdp))
end

fdr_estimates = entrapmentAnalysis(
    quant_psms_folder,
    precursors[:base_pep_id],
    precursors[:entrapment_group_id]
)
plot(fdr_estimates[!,:decoy_fdr], fdr_estimates[!,:combined_fdp], xlim = (0, 0.02), ylim = (0, 0.02))
plot!([0, 0.5], [0, 0.5], color = :black, xlim = (0, 0.02), ylim = (0, 0.02))

test_precs = DataFrame(Tables.columntable(precursors))[!,[:proteome_identifiers,:accession_numbers,:sequence,:structural_mods,:is_decoy,:base_pep_id,:entrapment_group_id]]
sort!(test_precs,[:base_pep_id,:entrapment_group_id,:is_decoy])
