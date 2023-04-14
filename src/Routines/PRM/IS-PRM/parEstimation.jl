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
        heavy_scans = getScansInBounds(scan_adresses, precursor_chromatograms, value.light_prec_id, integration_bounds)
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

function initDesignMatrix(transitions::Set{String}, bounds::Set{Int})
    zeros(
                length(transitions)*length(bounds), 
                length(transitions) + 1
            )
end

function fillDesignMatrix(X::Matrix{Float64}, transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64}, bounds::Set{Int})
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
    
    X = fillDesignMatrix(initDesignMatrix(transition_union, bounds), light_transitions, transition_union, light_scans, bounds)
    y = getResponseMatrix(heavy_transitions, transition_union, heavy_scans)
    return glmnet(X, y, 
                    penalty_factor = append!([0.0], ones(length(transition_union))), 
                    intercept = false, 
                    lambda=Float64[median(y)])
end

#=
    glmnet(X, y, penalty_factor = Float64[0, 1, 1, 1, 1, 1], intercept = true, lambda=Float64[median(y)])

    
    I = 1e9*Diagonal(ones(6))
    I[1,1] = 0
    B = (transpose(X)*X + I)\(transpose(X)*y)
    1 - sum((X*B - y).^2)/sum((y .- mean(y)).^2)
    B = (transpose(X[:,1])*X[:,1])\(transpose(X[:,1])*y)
    1 - sum((X[:,1]*B - y).^2)/sum((y .- mean(y)).^2)
=#