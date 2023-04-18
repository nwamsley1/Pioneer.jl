function getPARs(ptable::ISPRMPrecursorTable, 
                scan_adresses::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}, 
                precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram};
                minimum_scans::Int = 3,
                ms_file_idx::UInt32 = UInt32(0))
    for (key, value) in pairs(getIDToLightHeavyPair(ptable))
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
        if length(integration_bounds) < minimum_scans + 1
            continue
        end
        light_scans = getScansInBounds(scan_adresses, precursor_chromatograms, value.light_prec_id, integration_bounds)
        println("light scans ", light_scans)
        heavy_scans = getScansInBounds(scan_adresses, precursor_chromatograms, value.light_prec_id, integration_bounds)
        println("heavy scans ", heavy_scans)
        println("integration bounds ", integration_bounds)
        m = fitPAR(light_scans, 
                heavy_scans, 
                integration_bounds, 
                precursor_chromatograms[value.light_prec_id].transitions, 
                precursor_chromatograms[value.heavy_prec_id].transitions
                )
        setParModel(value, coef = m.betas.ca, dev_ratio = m.dev_ratio[1], ms_file_idx = ms_file_idx)
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
    [index for (index, scan_address) in enumerate(getScanAdressesForPrecursor(scan_adresses, precursor_chromatograms, prec_id)) if scan_address.ms1 in bounds]
end

function initDesignMatrix(transitions::Set{String}, bounds::AbstractArray)
    zeros(
                #length(transitions)*length(bounds), 
                length(transitions)*length(bounds),
                length(transitions) + 1
            )
end

function fillDesignMatrix(X::Matrix{Float64}, transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64}, bounds::AbstractArray)
    println("FILLX")
    println(transitions)
    println(transition_union)
    println(scans)
    println(bounds)
    for (i, transition) in enumerate(transition_union)
        X[((i-1)*length(bounds) + 1):(i*length(bounds)),1] = transitions[transition][scans]
        X[((i-1)*length(bounds) + 1):(i*length(bounds)),i+1] = transitions[transition][scans]
    end
    X
end

function getResponseMatrix(transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64})
    y = Float64[]
    for (i, transition) in enumerate(transition_union)
        append!(y, transitions[transition][scans])
    end
    y
end

function fitPAR(light_scans::Vector{Int64}, heavy_scans::Vector{Int64}, 
                bounds::Set{Int}, light_transitions::UnorderedDictionary{String, Vector{Float32}}, heavy_transitions::UnorderedDictionary{String, Vector{Float32}})

    transition_union = Set(keys(light_transitions))∩Set(keys(heavy_transitions))
    if length(transition_union) == 0
        return 
    end
    
    X = fillDesignMatrix(initDesignMatrix(transition_union, light_scans), light_transitions, transition_union, light_scans, light_scans)
    y = getResponseMatrix(heavy_transitions, transition_union, heavy_scans)
    return glmnet(X, y, 
                    penalty_factor = append!([0.0], ones(length(transition_union))), 
                    intercept = false, 
                    lambda=Float64[mean(y)])
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
    else
        return (1/lh_pair.par_model[ms_file_idx].par_model_coef[1],
                lh_pair.par_model[ms_file_idx].dev_ratio,
                isotope
                )
    end
end