using RobustModels
using Suppressor
using Test
function getPARs(ptable::ISPRMPrecursorTable,
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

        heavy_scans, light_scans = getScanPairs(getRTs(precursor_chromatograms[value.heavy_prec_id]),
                                                getRTs(precursor_chromatograms[value.light_prec_id]),
                                                max_diff = Float32(0.1))

        if length(heavy_scans) < minimum_scans + 1
            continue
        end

        par, goodness_of_fit = fitPAR(light_scans, 
                                      heavy_scans, 
                                      getTransitions(precursor_chromatograms[value.light_prec_id]), 
                                      getTransitions(precursor_chromatograms[value.heavy_prec_id])
                                      )

        setParModel(value, coef = par[1], goodness_of_fit = goodness_of_fit, ms_file_idx = ms_file_idx)
    end
end


function fillDesignMatrix(transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64})
    X = zeros(length(transition_union)*length(scans),1)
    i = 1
    for transition in transition_union
        for scan in scans #transitions[transition][scans]
            X[i] = transitions[transition][scan]
            i += 1
        end
    end
    return X
end

function getResponseMatrix(transitions::UnorderedDictionary{String, Vector{Float32}}, transition_union::Set{String}, scans::Vector{Int64})
    y = zeros(length(scans)*length(transition_union))
    i = 1
    for transition in transition_union
        for scan in scans
        y[i] = transitions[transition][scan]
        i += 1
        end
    end
    y
end

function fitPAR(light_scans::Vector{Int64}, 
                heavy_scans::Vector{Int64},
                light_transitions::UnorderedDictionary{String, Vector{Float32}}, 
                heavy_transitions::UnorderedDictionary{String, Vector{Float32}})

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
    X = fillDesignMatrix(light_transitions, transition_union, light_scans)
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
            par = RobustModels.coef(model)[1]
            GOF = RobustModels.stderror(model)[1]/RobustModels.coef(model)[1]
        end
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