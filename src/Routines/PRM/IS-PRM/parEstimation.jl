using RobustModels
using Suppressor
using Test
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
                            
        if length(integration_bounds) < minimum_scans + 1
            continue
        end

        light_scans = getScansInBounds(scan_adresses, precursor_chromatograms, value.light_prec_id, integration_bounds)
        heavy_scans = getScansInBounds(scan_adresses, precursor_chromatograms, value.heavy_prec_id, integration_bounds)

        par, goodness_of_fit = fitPAR(light_scans, 
                heavy_scans, 
                integration_bounds, 
                precursor_chromatograms[value.light_prec_id].transitions, 
                precursor_chromatograms[value.heavy_prec_id].transitions
                )

        setParModel(value, coef = par, goodness_of_fit = goodness_of_fit, ms_file_idx = ms_file_idx)
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
                #length(transitions) + 1
                1
            )
end

function fillDesignMatrix(X::Matrix{Float64}, transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64}, bounds::AbstractArray)
    i = 1
    for transition in transition_union
        for scan in transitions[transition][scan]
            X[i] = transitions[transition][scan]
            i += 1
        end
    end
    return X
end

function getResponseMatrix(transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64})
    y = zeros(scans*length(transition_union))
    i = 1
    for transition in transition_union
        for scan in transitions[transition][scans]
        y[i] = scan
        i += 1
        end
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
            if sum(dict[key].>=Float32(1e-12)) >= min_non_zero
                    push!(non_zero_keys, key)
            end
        end
        return non_zero_keys
    end

    transition_union = find_keys_with_greatest_sums(heavy_transitions, transition_union)
    transition_union = find_keys_with_nonzero(light_transitions, transition_union, 5)

    if length(transition_union) <= 2
        return 0.0, 0.0
    end

    #println("NEW")
    X = fillDesignMatrix(initDesignMatrix(transition_union, light_scans), light_transitions, transition_union, light_scans, light_scans)
    y = getResponseMatrix(heavy_transitions, transition_union, heavy_scans)

    if sum(X)==0 #This is hacky. Need to fix. 
        X[1] = 1.0
        println("sum(X) == 0")
    end

    function getModel(X::Matrix{T}, y::Vector{T}, I::BitVector, loss::AbstractEstimator) where T <: AbstractFloat
        rlm(X[I,:]./mean(y),y[I]./mean(y), loss, initial_scale=:mad, maxiter = 200)
    end

    par = 0.0
    GOF = 0.0

    try #If model fails to converge give default values. 
        @suppress begin #Suppress warnings about convergence 
            model =  getModel(X, y, ((X[:,1].!=0.0) .& (y.!=0.0)), TauEstimator{TukeyLoss}())
        end
        par = RobustModels.coef(model)[1]
        GOF = RobustModels.stderror(model)[1]/RobustModels.coef(model)[1]
    catch
        par = 0.0
        GOF = 0.0
    end

    return par, GOF
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
        return (1/getCoef(lh_pair.par_model[ms_file_idx]),
                getGoodnessOfFit(lh_pair.par_model[ms_file_idx]),
                isotope
                )
    end
end