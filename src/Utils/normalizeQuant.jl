function groupPSMS(
                    psms::DataFrame;
                    max_q_value::AbstractFloat = 0.01,
                    min_points_above_FWHM::Int = 2)
    return groupby(psms[
        (psms[!,:q_value].<=max_q_value).&(psms[!,:target]),:],
        :file_name)
end

function getQuantSplines(
    gpsms::GroupedDataFrame{DataFrame},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7)
    quant_splines = Dictionary{String, UniformSpline}()
    for key in keys(gpsms)
        psms = gpsms[key]
        #psms = psms[psms[!,:species].=="HUMAN",:]
        sort!(psms, :RT, alg = QuickSort)
        nprecs = size(psms, 1)
        bin_size = nprecs÷N
        median_quant = zeros(Float64,N);
        median_rts = zeros(Float64, N);

        for bin_idx in range(1, N)
            bin_start = (bin_idx - 1)*bin_size + 1
            bin_stop = bin_idx*bin_size
            median_rts[bin_idx] = median(@view(psms[bin_start:bin_stop,:RT]));
            median_quant[bin_idx] =log2(median(@view(psms[bin_start:bin_stop,quant_col_name])));
        end
        
        splinefit = UniformSpline(
            median_quant,
            median_rts,
            3, spline_n_knots
        )

        insert!(quant_splines, key[:file_name], splinefit)
    end
    return quant_splines
end

function getQuantCorrections(
    quant_splines::Dictionary{String, UniformSpline},
    rt_range::Tuple{AbstractFloat, AbstractFloat};
    N = 100)
    median_quant = zeros(Float32, (length(quant_splines), N))
    rt_range = collect(LinRange(first(rt_range), last(rt_range), N))
    for (i, rt) in enumerate(rt_range)
        j = 1
        for (key, spline) in pairs(quant_splines)
            median_quant[j, i] = spline(rt_range[i])
            j += 1
        end
    end
    median_quant = reshape(median(median_quant, dims = 1), (N,))
    median_quant = LinearInterpolation(rt_range, median_quant)
    corrections = Dictionary{String, Any}()
    for (key, spline) in pairs(quant_splines)
        offset = zeros(Float32, N)
        for (i, rt) in enumerate(rt_range)
            offset[i] = spline(rt) - median_quant(rt)
        end
        insert!(
            corrections, 
            key,
            LinearInterpolation(rt_range, offset)
        )
    end
    return corrections 
end

function applyNormalization!(
    psms::DataFrame,
    quant_col::Symbol,
    corrections::Dictionary{String, Any}
)
    norm_quant_col = Symbol(String(quant_col)*"_normalized")
    psms[!,norm_quant_col] = zeros(Float32, size(psms, 1))
    for i in range(1, size(psms, 1))
        hc = corrections[psms[i,:file_name]](psms[i,:RT])
        psms[i,norm_quant_col] = 2^(log2(max(psms[i,quant_col], 0.0)) - hc)
    end
end

function normalizeQuant(
    psms::DataFrame,
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    max_q_value::AbstractFloat = 0.01,
    min_points_above_FWHM::Int = 2)

    rt_range = (
        minimum(psms[!,:RT]),
        maximum(psms[!,:RT]),
    )
    gpsms = groupPSMS(psms, 
                        max_q_value = max_q_value,
                        min_points_above_FWHM = min_points_above_FWHM)

    quant_splines_dict = getQuantSplines(
        gpsms,
        quant_col_name, 
        N = N,
        spline_n_knots = spline_n_knots)

    quant_corrections_dict = getQuantCorrections(
        quant_splines_dict,
        rt_range,
        N = N
    )

    applyNormalization!(
        psms, 
        quant_col_name,
        quant_corrections_dict
    )

end
