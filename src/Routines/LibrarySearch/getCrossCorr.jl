function getCrossCorr(MS1::GroupedDataFrame{DataFrame}, MS2::GroupedDataFrame{DataFrame}, precursor_idx::I; lag::Int = 10) where {I<:Integer}
    function missingID(gdf::GroupedDataFrame{DataFrame}, precursor_idx::I)
        if !((precursor_idx=precursor_idx,) in keys(gdf)) #If the precursor is not found
            return true
        end
        return false
    end
    if missingID(MS1, precursor_idx) | missingID(MS2, precursor_idx)
        return (missing, missing)
    end
    return getCrossCorr(MS1[(precursor_idx =precursor_idx,)], MS2[(precursor_idx = precursor_idx, )], lag)
end

function getCrossCorr(MS1::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, MS2::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, lag::Int64 = 10)
    scan_pairs = getScanPairs(MS1[:,:rt], MS2[:,:rt])
    lag = min(lag, length(scan_pairs[1]) - 1)
    if lag < 3
        return (missing, missing)
    end
    lag = collect(range(-lag, lag))
    cross_cor = crosscor(MS1[scan_pairs[1],:weight], MS2[scan_pairs[2],:weight], lag)
    if all(isnan.(cross_cor))
        return (missing, missing)
    end
    return lag[argmax(cross_cor)], maximum(cross_cor)
end

