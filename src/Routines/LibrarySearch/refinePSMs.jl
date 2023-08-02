
function refinePSMs!(PSMs::DataFrame, precursors::Vector{LibraryPrecursor{T}}; min_spectral_contrast::AbstractFloat = 0.9,  n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors[psm[:precursor_idx]])) => :decoy)
    transform!(PSMs, AsTable(:) => ByRow(psm -> precursors[psm[:precursor_idx]].missed_cleavages) => :missed_cleavage)
    transform!(PSMs, AsTable(:) => ByRow(psm -> precursors[psm[:precursor_idx]].sequence) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(precursors[psm[:precursor_idx]]))) => :iRT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:weight] < 10.0) => :nmf)

    best_psms = combine(sdf -> sdf[argmax(sdf.matched_ratio), :], groupby(PSMs[(PSMs[:,:spectral_contrast_all].>min_spectral_contrast) .& (PSMs[:,:decoy].==false),:], [:scan_idx]))
    linear_spline = rtSpline(best_psms[:,:iRT], best_psms[:,:RT], n_bins = n_bins, granularity = granularity)
    PSMs[:,:RT_pred] = linear_spline(PSMs[:,:iRT])
    PSMs[:,:RT_error] = abs.(PSMs[:,:RT_pred] .- PSMs[:,:RT])

    sort!(PSMs, [:scan_idx, :total_ions]);
    # Group DataFrame by "day" column
    grouped_df = groupby(PSMs, :scan_idx);

    #PSMs[:,:next_best] = Vector{Union{Missing, UInt32}}(undef, size(PSMs)[1])
    PSMs[:,:next_best] = (combine(grouped_df) do sub_df
        pushfirst!(diff(sub_df.total_ions), zero(UInt32))
    end)[:,:x1]

    PSMs[:,:diff_hyper] = (combine(grouped_df) do sub_df
        sort!(sub_df, :hyperscore)
        pushfirst!(diff(sub_df.hyperscore), zero(Float64))
    end)[:,:x1]

    PSMs[:,:rank_hyper] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.hyperscore)
    end)[:,:x1]

    PSMs[:,:rank_scribe] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.scribe_score)
    end)[:,:x1]

    PSMs[:,:rank_poisson] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.poisson)
    end)[:,:x1]

    PSMs[:,:rank_total] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.total_ions)
    end)[:,:x1]

    PSMs[:,:diff_scribe] = (combine(grouped_df) do sub_df
        sort!(sub_df, :scribe_score)
        pushfirst!(diff(sub_df.scribe_score), zero(Float64))
    end)[:,:x1]

    PSMs[:,:median_ions] = (combine(grouped_df) do sub_df
        repeat([median(sub_df.total_ions)], size(sub_df)[1])
    end)[:,:x1]

    grouped_df = groupby(PSMs, :precursor_idx);

    PSMs[:,:n_obs] = (combine(grouped_df) do sub_df
        repeat([size(sub_df)[1]], size(sub_df)[1])
    end)[:,:x1]

    transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(precursors[psm[:precursor_idx]])) => :charge)

end

function rtSpline(X::Vector{T}, Y::Vector{T}; n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    sort_order = sortperm(X)

    #Divide RT space into estimation bins
    est_bins = [Int(bin÷1) for bin in range(1, length = n_bins, stop = length(sort_order))]

    #x and y values for each RT estimation bin
    xs = Vector{T}(undef, length(est_bins) - 1)
    ys = Vector{T}(undef, length(est_bins) - 1)
    for i in 1:(length(est_bins) - 1)

        #RTs for the i'th estimation bin
        obs = X[sort_order[est_bins[i]:est_bins[i + 1]]]
        x = Vector(LinRange(minimum(obs), maximum(obs), granularity))
        kde = KDEUniv(ContinuousDim(), 3.0, obs, MultiKDE.gaussian)
        y = [MultiKDE.pdf(kde, _x, keep_all=false) for _x in x]
        xs[i] = x[argmax(y)]
        ys[i] = mean(Y[sort_order[est_bins[i]]])#x[sort_order[bins[i]:bins[i + 1]]])
    end
    return LinearInterpolation(xs, ys, extrapolation_bc = Line() )
end