# Array management functions for MainSearch deconvolution.
# Used by dispatched helpers in process_scans.jl.

"""
    zero_negligible_weights!(weights, n_cols; eps=1e-6f0)

Zero out weights below `eps` before expensive spectral scoring.
"""
function zero_negligible_weights!(
    weights::Vector{Float32},
    n_cols::Int;
    eps::Float32 = Float32(1e-6)
)
    @inbounds for col in 1:n_cols
        if weights[col] < eps
            weights[col] = zero(Float32)
        end
    end
    return nothing
end

"""
    resize_arrays!(search_data::SearchDataStructures, weights::Vector{Float32})

Resize working arrays when needed for larger deconvolution problems.
"""
function resize_arrays!(search_data::SearchDataStructures, weights::Vector{Float32})
    new_entries = n_active(getIdToCol(search_data)) - length(weights) + 1000
    resize!(weights, length(weights) + new_entries)
    resize!(getColNorm2(search_data), length(getColNorm2(search_data)) + new_entries)
    resize!(getMainSearchSpectralScores(search_data), length(getMainSearchSpectralScores(search_data)) + new_entries)
    # Grow both PSM buffer families so any active search method (Main or
    # Tuning) has enough slots; cheap and only triggers on >5000 active prec.
    append!(getMainUnscoredPsms(search_data), [eltype(getMainUnscoredPsms(search_data))() for _ in 1:new_entries])
    append!(getTuningUnscoredPsms(search_data),  [eltype(getTuningUnscoredPsms(search_data))()  for _ in 1:new_entries])
end
