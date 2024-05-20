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

b = zeros(Float32, 200);
A = getWittakerHendersonDesignMat(length(b), 1.0f0);
prob = LinearProblem(A, b);
linsolve = init(prob);
u2 = zeros(Float32, length(linsolve.b));

subchrom = gchroms[(precursor_idx = 0x000c4cc1, iso_rank = (0x00, 0x03))]#groupby(gchroms[(precursor_idx = best_precursors[N,1],)],:iso_rank);

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
α = 0.01f0,
height_at_integration_width = 0.001f0,
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