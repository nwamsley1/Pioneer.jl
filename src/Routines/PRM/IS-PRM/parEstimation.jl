function getPARs(ptable::ISPRMPrecursorTable, 
                scan_adresses::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}, 
                precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram};
                minimum_scans::Int = 3,
                ms_file_idx::UInt32 = UInt32(0))
    i = 0
    for (key, value) in pairs(getIDToLightHeavyPair(ptable))
        #if i%100 == 0
        #    print(i,",")
        #end
        #i+=1
        if !isassigned(precursor_chromatograms, value.light_prec_id)
            continue
        end
        if !isassigned(precursor_chromatograms, value.heavy_prec_id)
            continue
        end
        integration_bounds = Set(
                                getIntegrationBounds(
                                    getScanCycleUnion(
                                        getScanAdressesForPrecursor(scan_adresses, precursor_chromatograms, value.light_prec_id),
                                        getScanAdressesForPrecursor(scan_adresses, precursor_chromatograms, value.heavy_prec_id)
                                        )
                                    )
                                )
                            
       #= if value.light_sequence == "AANLLLGYK"
            println("AANLLLGYK")
            println("light scan adrewses ")
            [println(x) for x in getScanAdressesForPrecursor(scan_adresses, precursor_chromatograms, value.light_prec_id)]
            println("heavy scan adrewses ")
            [println(x) for x in getScanAdressesForPrecursor(scan_adresses, precursor_chromatograms, value.heavy_prec_id)]
            println("integration_bounds ", integration_bounds)
        end=#
        if length(integration_bounds) < minimum_scans + 1
            continue
        end

        light_scans = getScansInBounds(scan_adresses, precursor_chromatograms, value.light_prec_id, integration_bounds)
        heavy_scans = getScansInBounds(scan_adresses, precursor_chromatograms, value.heavy_prec_id, integration_bounds)

        betas, dev_ratio = fitPAR(light_scans, 
                heavy_scans, 
                integration_bounds, 
                precursor_chromatograms[value.light_prec_id].transitions, 
                precursor_chromatograms[value.heavy_prec_id].transitions
                )

        setParModel(value, coef = betas, dev_ratio = dev_ratio, ms_file_idx = ms_file_idx)
    end
end

function getScanAdressesForPrecursor(scan_adresses::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}},
                                        precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram},
                                        prec_id::UInt32)
    scan_adresses[precursor_chromatograms[prec_id].scan_idxs[2:end]]
end

function getScansInBounds(scan_adresses::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}},
                    precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram},
                    prec_id::UInt32, 
                    bounds::Set{Int64})
    scans = [(index, scan_address) for (index, scan_address) in enumerate(getScanAdressesForPrecursor(scan_adresses, precursor_chromatograms, prec_id)) if scan_address.ms1 in bounds]
    #println("scans ", scans)
    function keep_first_occurrences_sorted(tuples::Vector{<:Tuple})
        last_second = nothing
        new_tuples = Int64[]
        for tup in tuples
            if last_second != tup[2].ms1
                push!(new_tuples, tup[1])
                last_second = tup[2].ms1
            end
        end
        return new_tuples
    end
    keep_first_occurrences_sorted(scans)
end

function initDesignMatrix(transitions::Set{String}, bounds::AbstractArray)
    zeros(
                #length(transitions)*length(bounds), 
                length(transitions)*length(bounds),
                length(transitions) + 1
            )
end

function fillDesignMatrix(X::Matrix{Float64}, transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64}, bounds::AbstractArray)

    for (i, transition) in enumerate(transition_union)
        X[((i-1)*length(bounds) + 1):(i*length(bounds)),1] = transitions[transition][scans]
        X[((i-1)*length(bounds) + 1):(i*length(bounds)),i+1] = transitions[transition][scans]
    end
    X
end

function getResponseMatrix(transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64})
    y = Float64[]
    for (i, transition) in enumerate(transition_union)
        #println("transition ", transition)
        #println("transitions[transition] ", transitions[transition])
        #println("scans ", scans)
        append!(y, transitions[transition][scans])
    end
    y
end

function fitPAR(light_scans::Vector{Int64}, heavy_scans::Vector{Int64}, 
                bounds::Set{Int}, light_transitions::UnorderedDictionary{String, Vector{Float32}}, heavy_transitions::UnorderedDictionary{String, Vector{Float32}})

    transition_union = Set(keys(light_transitions))∩Set(keys(heavy_transitions))
    if length(transition_union) == 0
        return zeros((1, 3)), 0.0
    end

    function find_keys_with_greatest_sums(dict::UnorderedDictionary{String, Vector{Float32}}, key_set::Set{String})
        sorted_keys = sort(collect(key_set), by=x->sum(dict[x]), rev=true)

        return Set(sorted_keys[1:min(4, length(key_set))])
    end

    function find_keys_with_nonzero(dict::UnorderedDictionary{String, Vector{Float32}}, set::Set{String}, min_non_zero::Int)
        non_zero_keys = Set(String[])
        for key in set
            if sum(dict[key].!=Float32(0.0)) >= min_non_zero
                    push!(non_zero_keys, key)
            end
        end
        return non_zero_keys
    end
    transition_union = find_keys_with_greatest_sums(heavy_transitions, transition_union)
    transition_union = find_keys_with_nonzero(light_transitions, transition_union, 5)

    if length(transition_union) <= 1
        return zeros((1, 3)), 0.0
    end

    #println("NEW")
    X = fillDesignMatrix(initDesignMatrix(transition_union, light_scans), light_transitions, transition_union, light_scans, light_scans)
    y = getResponseMatrix(heavy_transitions, transition_union, heavy_scans)

    if sum(X)==0 #This is hacky. Need to fix. 
        X[1] = 1.0
        println("sum(X) == 0")
    end
    m =  glmnet(X, y, 
                    penalty_factor = append!([0.0], ones(length(transition_union))), 
                    intercept = false, 
                    lambda=Float64[mean(y)*mean(y)])
    #println(out.betas)
    #println(out.dev_ratio)
    return m.betas.ca, m.dev_ratio[1]
end

function getPAR(ptable::ISPRMPrecursorTable, prec_id::UInt32, ms_file_idx::UInt32)
    lh_pair = ptable.lh_pair_id_to_light_heavy_pair[ptable.prec_id_to_lh_pair_id[prec_id]]
    isotope = "light"
    if lh_pair.heavy_prec_id == prec_id
        isotope = "heavy"
    end
    if length(lh_pair.par_model) == 0
        return (missing, missing, isotope)
    elseif !isassigned(lh_pair.par_model, ms_file_idx)
        return (missing, missing, isotope)
    elseif isinf(1/lh_pair.par_model[ms_file_idx].par_model_coef[1])
        return (missing, missing, isotope)
    else
        #println("1/lh_pair.par_model[ms_file_idx].par_model_coef[1]", 1/lh_pair.par_model[ms_file_idx].par_model_coef[1])
        #println(lh_pair.par_model[ms_file_idx])
        return (1/lh_pair.par_model[ms_file_idx].par_model_coef[1],
                lh_pair.par_model[ms_file_idx].dev_ratio,
                isotope
                )
    end
end