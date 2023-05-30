using RobustModels
using Suppressor
using Test

"""
    getPARs(ptable::ISPRMPrecursorTable,precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram};minimum_scans::Int = 3,ms_file_idx::UInt32 = UInt32(0))

Estimates peak area ratios and goodness-of-fit for light/heavy precursor pairs. 

### Input

- `ptable::ISPRMPrecursorTable` -- See src/Routines/PRM/IS-PRM/buildPrecursorTable.
- `precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}` -- Dictionary of `PrecursorChromatograms` for each precursor. 
- `minimum_scans::Int` -- Minimum number of scans required to estimate a PAR. 
- `ms_file_idx::UInt32` -- Identifier for the MS raw data 

### Output
- Estimates a peak area ratio between the light (endogenous) and heavy (internal standard) peptides represented in 
`precursor_chromatograms`. Each `LightHeavyPrecursorPair` in `ptable::ISPRMPrecursorTable` has an attribute `par_model::Dictionary{UInt32, ParModel}` 
which maps ms_file_idx::UInt32 to `ParModel` objects. `ParModel`s have attributes par_model_coef and goodness_of_fit. `getPARs` estimates 
these attributes using `fitPAR` and assigns themr to the `par_model` attribute. 

### Notes

### Examples 

"""
function getPARs(ptable::ISPRMPrecursorTable,
                precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram};
                minimum_scans::Int = 3,
                ms_file_idx::UInt32 = UInt32(0))
    #Both light and heavy versions of the precursor neeed to have `PrecursorChromatograms`
    for (key, value) in pairs(getIDToLightHeavyPair(ptable))
        if !isassigned(precursor_chromatograms, value.light_prec_id)
            continue
        end
        if !isassigned(precursor_chromatograms, value.heavy_prec_id)
            continue
        end

        #See src/Routines/PRM/IS-PRM/getScanPairs.jl. Identifies 1-1 correspondene between scans 
        #targeting the light and heavy version of the precursor. Returns these as Vector{Int}. 
        #Int is the scan ID. 
        heavy_scans, light_scans = getScanPairs(getRTs(precursor_chromatograms[value.heavy_prec_id]),
                                                getRTs(precursor_chromatograms[value.light_prec_id]),
                                                max_diff = Float32(0.1))

        #Require at least `minimum_scans` targeting the heavy version of the precursor
        if length(heavy_scans) < minimum_scans + 1
            continue
        end

        #Estimate the peak area ratio by MM-Estimation without and intercept. Standard error of the slope 
        #estimate divided by the absolute value of the slope is the `goodness_of_fit`. See `fitPAR`
        par, goodness_of_fit = fitPAR(light_scans, 
                                      heavy_scans, 
                                      getTransitions(precursor_chromatograms[value.light_prec_id]), 
                                      getTransitions(precursor_chromatograms[value.heavy_prec_id])
                                      )

        #Method defind on `LightHeavyPrecursorPair` type. Assignes the `par_model` attribute. See 
        #src/Routines/PRM/IS-PRM/buildPrecursorTable.jl
        setParModel(value, coef = par[1], goodness_of_fit = goodness_of_fit, ms_file_idx = ms_file_idx)
    end
end

"""
    fillDesignMatrix(transitions::UnorderedDictionary{String, Vector{T}}, transition_union::Set{String}, scans::Vector{Int64}) where {T<:Real}

Sets design matrix for par estimation. 

### Input

- `transitions::UnorderedDictionary{String, Vector{T}}` -- Dicionary mapping transitions names to vectors of transitoin intensities
- `transition_union::Set{String}` -- Should be keys in `transitions`. These are the transitions to consider
- `scans::Int` -- Indices of selected scans . 

### Output
- Nx1 matrix of transition intensities

### Notes

### Examples 
julia> transitions = UnorderedDictionary(
               ["y3+1","y4+1","y5+1"],
               [
                   Float32[1.1273232f7, 8.631513f6, 5.4651805f6, 3.0412648f6, 1.6092901f6, 1.024025f6],
                   Float32[4.8260845f6, 3.7873248f6, 2.3309442f6, 1.3229029f6, 704660.44, 436613.62],
                   Float32[5.183467f6, 3.9857915f6, 2.6672412f6, 1.4236808f6, 772209.56, 519827.88]
               ]
           )
3-element UnorderedDictionary{String, Vector{Float32}}
 "y3+1" │ Float32[1.1273232f7, 8.631513f6, 5.4651805f6, 3.0412648f6, 1.6092901f6, 1.024025f6]
 "y4+1" │ Float32[4.8260845f6, 3.7873248f6, 2.3309442f6, 1.3229029f6, 704660.44, 436613.62]
 "y5+1" │ Float32[5.183467f6, 3.9857915f6, 2.6672412f6, 1.4236808f6, 772209.56, 519827.88]

julia> transition_union = Set(["y3+1","y4+1","y5+1"])
Set{String} with 3 elements:
  "y3+1"
  "y4+1"
  "y5+1"

julia> scans = [1, 2, 3]
3-element Vector{Int64}:
 1
 2
 3

julia> fillDesignMatrix(transitions, transition_union, scans)
9x1 Matrix{Float64}:
 1.1273232e7
 8.631513e6
 5.4651805e6
 4.8260845e6
 3.78732475e6
 2.33094425e6
 5.183467e6
 3.9857915e6
 2.66724125e6
"""
function fillDesignMatrix(transitions::UnorderedDictionary{String, Vector{T}}, transition_union::Set{String}, scans::Vector{Int64}) where {T<:Real}
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

"""
    getResponseMatrix(transitions::UnorderedDictionary{String, Vector{T}}, transition_union::Set{String}, scans::Vector{Int64}) where {T<:Real}

Sets response matrix for par estimation. 

### Input

- `transitions::UnorderedDictionary{String, Vector{T}}` -- Dicionary mapping transitions names to vectors of transitoin intensities
- `transition_union::Set{String}` -- Should be keys in `transitions`. These are the transitions to consider
- `scans::Int` -- Indices of selected scans . 

### Output
- Vector of transition intensiteis for selected scans/transitions

### Notes

### Examples 
transitions = UnorderedDictionary(
                      ["y3+1","y4+1","y5+1"],
                      [
                          Float32[1.1273232f7, 8.631513f6, 5.4651805f6, 3.0412648f6, 1.6092901f6, 1.024025f6],
                          Float32[4.8260845f6, 3.7873248f6, 2.3309442f6, 1.3229029f6, 704660.44, 436613.62],
                          Float32[5.183467f6, 3.9857915f6, 2.6672412f6, 1.4236808f6, 772209.56, 519827.88]
                      ]
                  )
3-element UnorderedDictionary{String, Vector{Float32}}
 "y3+1" │ Float32[1.1273232f7, 8.631513f6, 5.4651805f6, 3.0412648f6, 1.6092901f6, 1.024025f6]
 "y4+1" │ Float32[4.8260845f6, 3.7873248f6, 2.3309442f6, 1.3229029f6, 704660.44, 436613.62]
 "y5+1" │ Float32[5.183467f6, 3.9857915f6, 2.6672412f6, 1.4236808f6, 772209.56, 519827.88]

julia> transition_union = Set(["y3+1","y4+1","y5+1"])
Set{String} with 3 elements:
  "y3+1"
  "y4+1"
  "y5+1"

julia> scans = [1, 2, 3]
3-element Vector{Int64}:
 1
 2
 3

julia> getResponseMatrix(transitions, transition_union, scans)
9-element Vector{Float64}:
 1.1273232e7
 8.631513e6
 5.4651805e6
 4.8260845e6
 3.78732475e6
 2.33094425e6
 5.183467e6
 3.9857915e6
 2.66724125e6
"""
function getResponseMatrix(transitions::UnorderedDictionary{String, Vector{T}}, transition_union::Set{String}, scans::Vector{Int64}) where {T<:Real}
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

"""
    fitPAR(light_scans::Vector{Int64}, heavy_scans::Vector{Int64}, light_transitions::UnorderedDictionary{String, Vector{Float32}}, heavy_transitions::UnorderedDictionary{String, Vector{Float32}})

    Estimates peak the peak area ratio for a light heavy pair by MM-estimation. 

### Input

- `light_scans::Vector{Int64}` -- Scan indices to select for the light (endogenous) precursor. Must have same length as `heavy_scans`
- `heavy_scans::Vector{Int64}` -- Scan indices to select for the heavy (internal standard) precursor. Must have same length as `light_scans`
- `light_transitions::UnorderedDictionary{String, Vector{Float32}}` -- Dictionary mapping transitions to intensities. Must have same keys as `heavy_transitions`
- `heavy_transitions::UnorderedDictionary{String, Vector{Float32}}` -- Dictionary mapping transitions to intensities. Must have same keys as `light_transitions`

### Output
- Output is Tuple{Float64, Float64}. First entry is the peak area ratio, the second entry is the goodness_of_fit. 

### Notes

### Examples 

"""
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

    #Design matrix is a single column matrix with intensities of the light precursor fragment ions
    X = fillDesignMatrix(light_transitions, transition_union, light_scans)
    #A vector of intensities of the heavy precursor fragment ions. Has 1-1 correspondence with X. 
    y = getResponseMatrix(heavy_transitions, transition_union, heavy_scans)

    if sum(X)==0 #This is hacky. Need to fix. 
        X[1] = 1.0
        println("sum(X) == 0")
    end

    # Estimator from RobustModels package 
    function getModel(X::Matrix{T}, y::Vector{T}, I::BitVector, loss::AbstractEstimator) where T <: AbstractFloat
        rlm(X[I,:]./mean(y),y[I]./mean(y), loss, initial_scale=:mad, maxiter = 200)
    end

    par = 0.0
    GOF = 0.0

    try #If model fails to converge give default values. 
        @suppress begin #Suppress warnings about convergence 
            #Ignore entires where light or heavy peptide had zero intensity 
            model =  getModel(X, y, ((X[:,1].!=0.0) .& (y.!=0.0)), TauEstimator{TukeyLoss}())   
            #Slope of the linear model (no intercept) is the peak area ratio estimate
            par = RobustModels.coef(model)[1]
            #Std error of the slope coefficient divided by the slope
            GOF = RobustModels.stderror(model)[1]/RobustModels.coef(model)[1]
        end
    catch
        par = 0.0
        GOF = 0.0
    end

    return par, GOF
end

"""
    getPAR(ptable::ISPRMPrecursorTable, prec_id::UInt32, ms_file_idx::UInt32)

    Retrieves the PAR and goodness_of_fit for a given precursor and experiment

### Input

- `ptable::ISPRMPrecursorTable` -- See src/Routines/PRM/IS-PRM/buildPrecursorTable.jl
- `prec_id::UInt32` -- Precursor id
- `ms_file_idx::UInt32` -- experiment id 

### Output
- Output is Tuple{Float64, Float64}. First entry is the peak area ratio, the second entry is the goodness_of_fit. Can be missing. 

### Notes

### Examples 

"""
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
        return (1/getCoef(lh_pair.par_model[ms_file_idx]),
                getGoodnessOfFit(lh_pair.par_model[ms_file_idx]),
                isotope
                )
    end
end