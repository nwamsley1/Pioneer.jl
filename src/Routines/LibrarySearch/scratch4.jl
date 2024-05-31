best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:file_name].=="01"),:]

best_precursors = unique(best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:file_name].=="01"),[:precursor_idx,:iso_rank]])
jldsave("/Users/n.t.wamsley/Desktop/best_precursors.jld2"; best_precursors)
best_precursors = load("/Users/n.t.wamsley/Desktop/best_precursors.jld2")["best_precursors"]
N = 10000

PSMS[!,:iRT] = RT_iRT["01"].(PSMS[!,:RT])

PSMS[!,:center_mass] = MS_TABLE[:centerMass][PSMS[!,:scan_idx]]

dtype = Float32;
gx, gw = gausslegendre(100);
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    );

    y = subchrom[1][!,:intensity]
b = zeros(Float32, 200)
A = getWittakerHendersonDesignMat(length(b), 0.5f0)
prob = LinearProblem(A, b)
linsolve = init(prob)

integratePrecursorMS2(MS2_CHROMS[(precursor_idx = best_precursors[N,1],iso_rank = best_precursors[N,2])],
state,
gx::Vector{Float64},
gw::Vector{Float64},
intensity_filter_fraction =  Float32(params_[:integration_params]["intensity_filter_threshold"]),
α = 0.001f0,
half_width_at_α = 0.05f0,
isplot = true
);

MS2_CHROMS[(precursor_idx = best_precursors[N,1],iso_rank = best_precursors[N,2])][!,
[:best_rank,:topn,:b_count,:y_count,:isotope_count,:scribe,:scribe_fitted,
:city_block_fitted,:spectral_contrast,:matched_ratio,:weight,:RT,:center_mass,:scan_idx,:iso_rank,:prob]]



PSMS[PSMS[!,:precursor_idx].==0x002cc660,#,best_precursors[N, 1],
[:best_rank,:topn,:b_count,:y_count,:isotope_count,:scribe,:scribe_fitted,
:city_block_fitted,:spectral_contrast,:matched_ratio,:weight,:RT,:scan_idx,:iso_rank,:prob]]
N += 1

rt = [16.287916
16.315779
16.343622
16.371744
16.399593]

wt = [
    804.01715
 3346.7717
 1581.042
 2357.8306
 2029.3892
]

pad!(wt, )
N += 1

diff(MS2_CHROMS[(precursor_idx = best_precursors[N,1],iso_rank = best_precursors[N,2])][!,:scan_idx])


getIsoRanks!(chroms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
gchroms = groupby(chroms,:precursor_idx)
#tet = 0x00368bc8
#=

julia> N
10153
=#
#wholechrom = gchroms[(precursor_idx = best_precursors[N,1],)]

b = zeros(Float32, 260);
A = getWittakerHendersonDesignMat(length(b), 1.0f0);
prob = LinearProblem(A, b);
linsolve = init(prob);
u2 = zeros(Float32, length(linsolve.b));

N = 0
for i in range(1, length(gchroms))
    if size(gchroms[i], 1) > N
        N =  size(gchroms[i], 1)
    end
end

subchrom = gchroms[N]#groupby(gchroms[(precursor_idx = best_precursors[N,1],)],:iso_rank);

dtype = Float32;
gx, gw = gausslegendre(100);
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    );
    integrateChrom(subchrom,
    linsolve,
    u2,
    state,
    gw,
    gx,
    α = 0.01f0,
    height_at_integration_width = 0.001f0,
    n_pad = 0,
    isplot = true
    )
    reset!(state)
integrateChrom(subchrom,
linsolve,
u2,
state,
gw,
gx,
α = 0.01f0,
height_at_integration_width = 0.001f0,
n_pad = 10,
isplot = true
)
reset!(state)
integrateChromFirst(subchrom,
       linsolve,
       u2,
       state,
       gw,
       gx,
       α = 0.01f0,
       height_at_integration_width = 0.001f0,
       #n_pad = 10,
       isplot = true
       )

N += 1

t = subchrom[1][!,:rt]
y = subchrom[1][!,:intensity]
plot(
t, y,
alpha = 0.5, seriestype=:scatter, color = :blue,
#xlim = (16.2, 16.5)
)
plot!(
t, y,
alpha = 0.5, color = :blue
)
sg = savitzky_golay(y, 11, 4)
plot!(
t, sg.y,
alpha = 0.5, seriestype=:scatter, color = :orange
)
PSMS[PSMS[!,:precursor_idx].==best_precursors[N, 1],
[:best_rank,:topn,:b_count,:y_count,:isotope_count,:scribe,:scribe_fitted,
:city_block_fitted,:spectral_contrast,:matched_ratio,:weight,:RT,:scan_idx,:iso_rank,:prob]]
N += 1

sg = savitzky_golay(y, 21, 9)
plot!(
t, sg.y,
alpha = 0.5, seriestype=:scatter, color = :orange
)
plot!(
t, sg.y,
alpha = 0.5, color = :orange
)



sg = savitzky_golay(y, 7, 3)
plot!(
t, sg.y,
alpha = 0.5, seriestype=:scatter, color = :green
)
plot!(
t, sg.y,
alpha = 0.5, color = :green
)



function initQuantScans(
    precs_to_integrate::Vector{DataFrames.GroupKey{GroupedDataFrame{DataFrame}}},
    )

    new_cols = [
        (:best_rank,                UInt8)
        (:topn,                     UInt8)
        (:longest_y,                UInt8)
        (:b_count,                  UInt8)
        (:y_count,                  UInt8)
        (:isotope_count,            UInt8)
        (:total_ions,               UInt8)
        (:poisson,                  Float16)
        (:hyperscore,               Float16)
        (:log2_intensity_explained,  Float16)
        (:error,                    Float32)
        (:error_norm,               Float32)

        (:scribe,                   Float16)
        (:scribe_fitted,            Float16)
        (:city_block,               Float16)
        (:city_block_fitted,        Float16)
        (:spectral_contrast,        Float16)
        (:matched_ratio,            Float16)
        (:entropy_score,            Float16)
        

        (:weight,                   Float32)
        (:peak_area,                Float32)
        (:trapezoid_area,           Float32)
        (:FWHM,                     Float16)
        (:FWHM_01,                  Float16)
        (:points_above_FWHM,        UInt16)
        (:points_above_FWHM_01,     UInt16)
        (:assymetry,                Float16)
        (:base_width_min,           Float16)

        (:max_score,                Float16)
        (:mean_score,               Float16)
        (:max_entropy,              Float16)
        (:max_scribe_score,         Float16)
        (:max_city_fitted,          Float16)
        (:mean_city_fitted,         Float16)
        (:y_ions_sum,               UInt16)
        (:max_y_ions,               UInt16)
        (:max_matched_ratio,        Float16)

        (:cv_fold,                  UInt8)
        (:target,                   Bool)
        (:RT,                       Float32)
        (:precursor_idx,            UInt32)
        #(:isotopes_captured,        Tuple{UInt8, UInt8})
        (:ms_file_idx,              UInt32)
        (:scan_idx,                 Union{Missing, UInt32})
        ];
    
        psms = DataFrame()
        N = length(precs_to_integrate)
        for column in new_cols
            col_type = last(column);
            col_name = first(column)
            psms[!,col_name] = zeros(col_type, N)
        end
        psms[!,:isotopes_captured] = Vector{Tuple{UInt8, UInt8}}(undef, N)
        return psms
end

sd = Vector{Union{Missing, Float64}}(undef, length(gpsms))
for i in range(1, length(gpsms))
    psms = gpsms[i]
    if size(psms, 1) < 3
        sd[i] = missing
        continue
    end
    sd[i] = std(psms[!,:iRT_observed])
end

test_psms = first(PSMs_Dict)

test_psms_pass = test_psms[test_psms[!,:q_value].<=0.01,:]
plot(test_psms_pass[!,:RT], 
        abs.(test_psms_pass[!,:iRT_predicted] .- test_psms_pass[!,:iRT_observed]),
        seriestype=:scatter,
        alpha = 0.01)


plot(test_psms_pass[!,:iRT_predicted], 
    test_psms_pass[!,:RT],
    seriestype=:scatter,
    alpha = 0.01)


plot!(LinRange(0, 30, 100), 
    first(iRT_RT).(LinRange(0, 30, 100)))

plot(test_psms_pass[!,:RT],
    test_psms_pass[!,:RT] .- first(iRT_RT).(test_psms_pass[!,:iRT_predicted]),
    seriestype=:scatter,
    alpha = 0.01)

test_psms = vcat(values(PSMs_Dict)...);
filter!(x->x.q_value<=0.01, test_psms)
psms_by_prec = groupby(test_psms, :precursor_idx)
rt_mads = Vector{Union{Missing, Float32}}(undef, length(psms_by_prec))

i = 1
for (prec, psms) in pairs(psms_by_prec)
    if size(psms, 1) > 2
        rt_mads[i] = maximum(psms[!,:iRT_observed]) -  minimum(psms[!,:iRT_observed])
    else
        rt_mads[i] = missing
    end
    i += 1
end

plot(skipmissing(test_rts), 
      skipmissing(rt_mads), 
      seriestype=:scatter,
      alpha = 0.01,
      ylim = (0, 0.5))


scans = coalesce.(abs.((MS_TABLE[:centerMass] .-  888.45f0)) .< 1e-6, false)
coalesce.(scans, false)

function integrateChromFirst(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, 
                                linsolve::LinearSolve.LinearCache,
                                u2::Vector{Float32},
                                state::GD_state{HuberParams{U}, V, I, J},
                                gw::Vector{Float64},
                                gx::Vector{Float64}; 
                                α::Float32 = 0.01f0, 
                                height_at_integration_width = 0.001f0,
                                isplot::Bool = false) where {U,V<:AbstractFloat, I,J<:Integer}
    
    #########
    #Helper Functions  
    #########
    function WHSmooth!( linsolve::LinearSolve.LinearCache, 
                        intensities::AbstractVector{Float32})
        #Reset linsolve and second derivative 
        @inbounds for i in range(1, length(linsolve.b))
            linsolve.b[i] = zero(Float32)
            linsolve.u[i] = zero(Float32)
            u2[i] = zero(Float32)
        end
        #Copy data to linsolve
        best_scan, max_intensity = 0, typemin(Float32)
        @inbounds for i in range(1, size(chrom, 1))
            if intensities[i]>max_intensity
                max_intensity = intensities[i]
                best_scan = i
            end
            linsolve.b[i] = intensities[i]
        end

        #WH smoothing 
        solve!(linsolve)

        #Best scan is the most intense 
        #best_scan = argmax(linsolve.b)
        
        return best_scan, linsolve.b[best_scan]
    end

    function fillU2!(
        u2::Vector{Float32},
        u::Vector{Float32})
        #Get second-order descrete derivative 
        u2[1], u2[end] = zero(Float32), zero(Float32)
        @inbounds @fastmath for i in range(2, length(linsolve.b) - 1)
            u2[i] = u[i + 1] - 2*u[i] + u[i - 1]
        end
    end

    function getIntegrationBounds!(u2::Vector{Float32},
                                   u::Vector{Float32},
                                   N::Int64,
                                   best_scan::Int64)
        start, stop = 1, 1
        start_search, stop_search = best_scan + 1, best_scan - 1

        #get RH boundary
        @inbounds @fastmath begin 
            for i in range(start_search, N-1)
                if (u2[i-1] < u2[i]) & (u2[i+1]<u2[i])
                    stop = min(i, N)
                    break
                end
            end
            for i in range(stop, N-1)
                if u[i + 1] > u[i]
                    break
                else
                    stop = i
                end
            end

            #get LH boundary 
            for i in reverse(range(2, stop_search))
                if (u2[i] > u2[i - 1]) & (u2[i+1]<u2[i])
                    start = max(i, 1)
                    break
                end
            end
            for i in reverse(range(2, start))
                if u[i - 1] > u[i]
                    break
                else
                    start = i
                end
            end
        end
        return range(start, stop)#range( min(best_scan-3, start), max(stop,best_scan+3))
    end

    function fillState!(state::GD_state{HuberParams{Float32}, Float32, Int64, Int64},
                        u::Vector{Float32},
                        rt::AbstractVector{Float16},
                        start::Int64, 
                        stop::Int64,
                        best_scan::Int64
                        )

        start_rt = rt[start]
        best_rt = rt[best_scan]
        #start_rt, best_rt = rt[start], rt[best_scan]
        rt_width = rt[stop] - start_rt
        norm_factor = u[best_scan]

        #Write data to state
        #Normalize so that maximum intensity is 1 
        #And time difference from start to finish is 1. 
        @inbounds @fastmath for i in range(1, stop - start + 1)
            n = start + i - 1
            state.t[i] = (chrom[n,:rt] - start_rt)/rt_width
            state.data[i] = u[n]/norm_factor
        end

        state.max_index = stop - start + 1
        best_rt = Float32((best_rt - start_rt)/rt_width)
        return norm_factor, start_rt, rt_width, best_rt
    end

    function fitEGH(state::GD_state{HuberParams{T}, U, I, J}, 
                    lower_bounds::HuberParams{T}, 
                    upper_bounds::HuberParams{T},
                    α::T,
                    half_width_at_α::T,
                    best_rt::T) where {T,U<:AbstractFloat, I,J<:Integer}

        #half_width_at_α = 0.15
        #Initial Parameter Guesses
        state.params = getP0(T(α), 
                            T(half_width_at_α), 
                            T(half_width_at_α),
                            T(best_rt),
                            T(0.6),
                            lower_bounds, upper_bounds)

        GD(state,
                lower_bounds,
                upper_bounds,
                tol = 1e-4, 
                max_iter = 300, 
                δ = 1e-5,#1e-3, #Huber loss parameter. 
                α=Float64(α),
                β1 = 0.9,
                β2 = 0.999,
                ϵ = 1e-8)
        
    end

    function getPeakProperties(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
        points_above_FWHM = zero(Int32)
        points_above_FWHM_01 = zero(Int32)
        for i in range(1, state.max_index)
            intensity = state.data[i]
            if intensity > (state.params.H*0.5)
                points_above_FWHM += one(Int32)
            end
            if intensity > (state.params.H*0.01)
                points_above_FWHM_01 += one(Int32)
            end 
        end

        FWHM = getFWHM(state, 0.5)
        FWHM_01 = getFWHM(state, 0.01)

        return FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01
    end
    #Baseline subtraction?
    function subtractBaseline!(
        u::Vector{Float32}, #smoothed data
        best_scan::Int64, #peak apex
        scan_range::UnitRange{Int64}) #start and stop of integration bounds 
        #Fine LH baseline 
        lmin,li = typemax(Float32),first(scan_range)
        @inbounds @fastmath for i in range(first(scan_range), best_scan)
            if u[i] < lmin
                lmin = u[i]
                li = i
            end
        end

        #Find RH baseline 
        rmin,ri = typemax(Float32),last(scan_range)
        @inbounds @fastmath for i in range(best_scan, last(scan_range))
            if u[i] < rmin
                rmin = u[i]
                ri = i
            end
        end


        #= Another option is to just use the extrema of the integration boundary 
        lmin = linsolve.u[first(scan_range)]
        li = first(scan_range)
        ri = last(scan_range)
        rmin = linsolve.u[last(scan_range)]
        =#
        #Subtract the baseline 
        h = (rmin - lmin)/(ri - li)
        @inbounds @fastmath for i in scan_range
            u[i] = u[i]-(lmin + (i - li)*h)
        end

    end

    function combinedShared(sa::SparseArray{Ti, T}) where {Ti<:Integer, T<:AbstractFloat}
        #Assumes sparse array is sorted in column major order. 
    end

    function integrateTrapezoidal(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
        retval = state.data[2]
        #Assumption that state.max_index is at least 4. 
        for i in range(3, state.max_index - 1)
            @inbounds retval = (state.t[2] - state.t[1])*(state.data[1] + state.data[2])
            @inbounds @fastmath for i in 2:(state.max_index - 1)
                retval += (state.t[i + 1] - state.t[i])*(state.data[i] + state.data[i + 1])
            end
        end
        return (1//2)*retval
    end
    #Whittaker Henderson Smoothing
    best_scan, max_intensity = WHSmooth!(
        linsolve,
        chrom[!,:intensity]
    )
    #Second discrete derivative of smoothed data
    fillU2!(
        u2,
        linsolve.u
    )
    
    #Integration boundaries based on smoothed second derivative 
    scan_range = getIntegrationBounds!(
        u2,
        linsolve.u,
        size(chrom, 1),
        best_scan,
    )

    subtractBaseline!(
        linsolve.u,
        best_scan,
        scan_range
    )
    
    #File `state` to fit EGH function. Get the inensity, and rt normalization factors 
    norm_factor, start_rt, rt_norm, best_rt = fillState!(
        state,
        linsolve.u,
        chrom[!,:rt],
        first(scan_range),
        last(scan_range),
        best_scan
    )
    
    #Initial estimate for FWHM
    a = max(best_scan - 1, 1)
    b = min(best_scan + 1, size(chrom, 1))

    half_width_at_α = Float32(((chrom[b,:rt]-start_rt)/rt_norm - (chrom[a,:rt] - start_rt)/rt_norm)/2)
    T = eltype(chrom.intensity)
    ##########
    #Fit EGH to data. 
    fitEGH(state, 
            HuberParams(T(0.001), T(0), T(-1), T(0.95)),
            HuberParams(T(1),  T(Inf), T(1), T(1.05)),
            α,
            half_width_at_α,
            best_rt
            )

    if isplot
        mi = state.max_index
        start = max(best_scan - 18, 1)
        stop = min(best_scan + 18, length(chrom.rt))
        plot(chrom.rt[start:stop], chrom.intensity[start:stop], seriestype=:scatter, alpha = 0.5, show = true)
        vline!([chrom.rt[first(scan_range)], chrom.rt[last(scan_range)]])
        plot!(state.t[1:mi].*rt_norm .+ start_rt, norm_factor.*state.data[1:mi], seriestype=:scatter, alpha = 0.5, show = true)
        xbins = LinRange(state.t[1]-0.5, state.t[state.max_index]+0.5, 100)
        plot!(xbins.*rt_norm .+ start_rt, [norm_factor*F(state, x) for x in xbins])
        plot!(chrom.rt[start:stop], u2[start:stop])
        hline!([norm_factor*0.95])
    end

    peak_area = rt_norm*norm_factor*Integrate(state, 
                                        gx,
                                        gw, 
                                        α = height_at_integration_width)


    trapezoid_area = rt_norm*norm_factor*integrateTrapezoidal(state)

    #trapezoid_area = 0.0f0
    FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01 = getPeakProperties(state)
    best_scan_idx = chrom[best_scan,:scan_idx]
    return best_scan_idx, peak_area, trapezoid_area, max_intensity,  FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01 
end