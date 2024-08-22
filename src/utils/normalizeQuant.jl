function groupPSMS(
                    psms::DataFrame;
                    max_q_value::AbstractFloat = 0.01,
                    min_points_above_FWHM::Int = 2)
    return groupby(psms[
        (psms[!,:q_value].<=max_q_value).&(psms[!,:target]),:],
        :file_name)
end

function getQuantSplines(psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7)
    quant_splines = Dictionary{String, UniformSpline}()
    min_rt, max_rt = typemax(Float32), typemin(Float32)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        sort!(psms, :irt_obs, alg = QuickSort)
        nprecs = size(psms, 1)
        bin_size = nprecs÷N
        median_quant = zeros(Float64,N);
        median_rts = zeros(Float64, N);

        for bin_idx in range(1, N)
            bin_start = (bin_idx - 1)*bin_size + 1
            bin_stop = bin_idx*bin_size
            median_rts[bin_idx] = median(@view(psms[bin_start:bin_stop,:irt_obs]));
            median_quant[bin_idx] =log2(median(@view(psms[bin_start:bin_stop,quant_col_name])));
        end
        
        splinefit = UniformSpline(
            median_quant,
            median_rts,
            3, spline_n_knots
        )
        if minimum(psms[!,:irt_obs])<min_rt
            min_rt = minimum(psms[!,:irt_obs])
        end
        if maximum(psms[!,:irt_obs])>max_rt
            max_rt = maximum(psms[!,:irt_obs])
        end
        insert!(quant_splines, fpath, splinefit)
    end
    return quant_splines, (min_rt, max_rt)
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
    psms_paths::Vector{String},
    quant_col::Symbol,
    corrections::Dictionary{String, Any}
)

    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        norm_quant_col = :peak_area_normalized
        correction_spline = corrections[fpath]
        for i in range(1, size(psms, 1))
            hc = correction_spline(psms[i,:irt_obs])
            psms[i,norm_quant_col] = 2^(log2(max(psms[i,quant_col], 0.0)) - hc)
        end
        Arrow.write(
            fpath,
            psms
        )
    end
end

function normalizeQuant(
    second_quant_folder::String,
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    max_q_value::AbstractFloat = 0.01,
    min_points_above_FWHM::Int = 2)


    psms_paths = [fpath for fpath in readdir(second_quant_folder, join=true) if endswith(fpath, ".arrow")]

    quant_splines_dict, rt_range = getQuantSplines(
        psms_paths,
        quant_col_name, 
        N = N,
        spline_n_knots = spline_n_knots)
    println("rt_range $rt_range")
    quant_corrections_dict = getQuantCorrections(
        quant_splines_dict,
        rt_range,
        N = N
    )

    applyNormalization!(
        psms_paths, 
        quant_col_name,
        quant_corrections_dict
    )
    return
end
