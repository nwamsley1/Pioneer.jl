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

struct ChromObject
    rt::Float16
    intensity::Float32
    scan_idx::UInt32
    precursor_idx::UInt32
end

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

subchrom = groupby(gchroms[(precursor_idx = best_precursors[N,1],)],:iso_rank);

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

integrateChrom(subchrom[1],
linsolve,
u2,
state,
gx::Vector{Float64},
gw::Vector{Float64},
α = 0.01f0,
height_at_integration_width = 0.001f0,
isplot = true
);
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



function getIsoRanks!(chroms::DataFrame, 
                        prec_charge::AbstractArray{UInt8},
                        prec_mz::AbstractArray{Float32},
                        MS_TABLE::Arrow.Table)
    #sum(MS2_CHROMS.weight.!=0.0)
    chroms[!,:iso_rank] = Vector{Tuple{UInt8, UInt8}}(undef, size(chroms, 1))
    
    tasks_per_thread = 5
    chunk_size = max(1, size(chroms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(chroms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk

                prec_id = chroms[i,:precursor_idx]
                mz = prec_mz[prec_id]
                charge = prec_charge[prec_id]

                scan_id = chroms[i,:scan_idx]
                scan_mz = MS_TABLE[:centerMass][scan_id]
                window_width = MS_TABLE[:isolationWidth][scan_id]

                isotopes = getPrecursorIsotopeSet(mz, 
                                                    charge, 
                                                    Float32(scan_mz-window_width/2),
                                                    Float32(scan_mz+window_width/2) 
                                                    )                
                chroms[i,:iso_rank] = isotopes
            end
        end
    end
    fetch.(tasks)
    return nothing
end

getIsoRanks!(chroms, precursors[:prec_charge],precursors[:mz], MS_TABLE)



