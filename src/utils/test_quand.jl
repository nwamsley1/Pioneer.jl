function UniformSplineDesignMat(
                        u::Vector{T},
                        t0::Vector{T}, 
                        t1::Vector{T},
                        degree::I, #Degree of the piecewise polynomials
                        n_knots::I, #Number of control points
                        ) where {I<:Integer, T<:AbstractFloat}
    if degree != 3
        error("Non-cubic splines not yet implemented. Use a degree of 3")
    end
    if n_knots < 3
        error("need at least 3 knots")
    end
    #if length(u) != length(t)
    #    error("length(u) is not equal to length(t)")
    #end

    #Uniform B Spline basis for the given degree
    #only implemented for d=3 but coule expand in the future 
    function getSplineBasis(degree::I)
        return NTuple{4, Polynomial}([
            Polynomial([0, 0, 0, 1])/6, #b0
            Polynomial([1, 3, 3, -3])/6, #b1
            Polynomial([4, 0, -6, 3])/6, #b2
            Polynomial([1, -3, 3, -1])/6, #b3
        ])
    end

    function buildDesignMat(t::Vector{T}, #location of data points
                            knots::Vector{T},
                            bin_width::T,
                            spline_basis::NTuple{4, Polynomial}
                    ) where {T<:AbstractFloat}

        function fillDesignMatRow!(X::Matrix{T}, 
                                    row::Int,
                                    knot_idx::Int,
                                    u::T,
                                    spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}
            i = length(spline_basis)
            #println(" t - knot_val: ", t - knot_val)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i](u)
                i -= 1
            end
        end

        X = zeros(T, (length(t), length(knots) + 3))
        for (i, t) in enumerate(t)
            #knot_idx = min(Int64((t - first(knots))÷bin_width)+1, length(knots))
            knot_idx = min(
                            floor(Int32, (t - first(knots))/bin_width)+one(Int32), 
                            length(knots)
                            )
            fillDesignMatRow!(
                X,
                i,
                knot_idx,
                (t-knots[knot_idx])/bin_width,
                spline_basis
            )
        end

        return X
    end

    function buildPieceWise(
                    knots::Vector{T},
                    bin_width::T,
                    spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}

        function fillDesignMatRow!(X::Matrix{Polynomial}, 
                                    row::Int,
                                    knot_idx::Int,
                                    spline_basis::NTuple{4, Polynomial})
            i = length(spline_basis)
            #println(" t - knot_val: ", t - knot_val)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i]
                i -= 1
            end
        end

        X = zeros(Polynomial, (length(knots), length(knots) + length(spline_basis) - 1))
        for (i, t) in enumerate(knots)
            t = t + bin_width/2
            knot_idx = min(
                floor(Int32, (t - first(knots))/bin_width)+one(Int32), 
                length(knots)
                )
            fillDesignMatRow!(
                X,
                i,
                knot_idx,
                spline_basis
            )
        end

        return X
    end
    spline_basis = getSplineBasis(degree)
    _first = min(minimum(t0), minimum(t1))
    _last = max(maximum(t0), maximum(t1)) 
    bin_width = (_last - _first)/(n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X0 = buildDesignMat(t0, collect(knots), bin_width, spline_basis)
    X1 = buildDesignMat(t1, collect(knots), bin_width, spline_basis)
    X = X0 .- X1
    c = X\u 
    plot(X*c, u, seriestype=:scatter, show=true, alpha = 0.1)
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly*c
    n_coeffs = n_knots*(degree + 1)
    coeffs = SVector{n_coeffs}(vcat([polynomial.coeffs for polynomial in piecewise_polynomials]...))


    return UniformSpline{n_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
end
function simmulateQuad(
    qtf::QuadTransmission,
    window_half_width::Float32,
    mz_offset_range::Tuple{Float32, Float32},
    center_mz_range::Tuple{Float32, Float32},
    n::Int64)

    x0 = Float32[]
    x1 = Float32[]
    yt = Float32[]
    for i in range(1, n)
        #Random window center mz 
        center_mz = Float32(rand(
            Uniform(first(center_mz_range), last(center_mz_range))
        ))
        #Random offset of mz from the center 
        offset = Float32(rand(
            Uniform(first(mz_offset_range), last(mz_offset_range))
            ))

        #Precursor charge state, m0 and m0 m/z, and m0 mass 
        prec_charge_state = rand([2])
        m0_mz = center_mz + offset
        m1_mz =  m0_mz + NEUTRON/prec_charge_state
        mono_mass = Float32(m0_mz*prec_charge_state) #Approximate not accounting for proton 


        s_count = 0 #Assume zero sulfurs.
        #Simmulate true abundances of mono and m+1 precursor ions 
        ζ0 = 100000
        δ = iso_splines(s_count, 0, mono_mass)/iso_splines(s_count, 1, mono_mass)
        ζ1 =  ζ0/δ
        #Observed m0 and m+1 abundances are modified by the quad transmission funciton 
        y0_obs = ζ0*qtf(center_mz, window_half_width, Float32(m0_mz))
        y1_obs = ζ1*qtf(center_mz, window_half_width, Float32(m1_mz))
        #Get offsets 
        push!(x1,m1_mz - center_mz) 
        push!(x0, m0_mz - center_mz)
        #Get observed data
        push!(yt, log2.(y0_obs/(y1_obs*δ)))
    end
    return DataFrame((x0 = x0, x1 = x1, yt = yt))
end
using Distributions 
qtf = QuadTransmission(0.0f0, 5.0f0)
quad_results = simmulateQuad(
    qtf,
    1.0f0,
    (-2.0f0, 2.0f0),
    (400.0f0, 1000.0f0),
    1000
)

gx = UniformSplineDesignMat(
    quad_results[!,:yt], 
    quad_results[!,:x0], 
    quad_results[!,:x1],
    3, 
    30
)
#plot(LinRange(-4, 4, 100), gx.(LinRange(-4, 4, 100)))
log2y = gx.(LinRange(-4, 4, 100))
y = exp2.(log2y .- maximum(log2y))
plot(LinRange(-4, 4, 100), y)
plot!(LinRange(-4, 4, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(-4, 4, 100))))


#=

testquad = QuadTransmission(params_[:quad_transmission]["overhang"], params_[:quad_transmission]["smoothness"])
plot(LinRange(98, 102, 100), [testquad(100.0, 1.0, x) for x in LinRange(98, 102, 100)], seriestype=:scatter)
plot(LinRange(98, 102, 100), log2.([testquad(100.0, 1.0, x) for x in LinRange(98, 102, 100)]), seriestype=:scatter)


testquad = QuadTransmission(10.0f0, 10000.0f0)
plot(LinRange(98, 102, 100), [testquad(100.0, 1.0, x) for x in LinRange(98, 102, 100)], seriestype=:scatter)
plot(LinRange(98, 102, 100), log2.([testquad(100.0, 1.0, x) for x in LinRange(98, 102, 100)]), seriestype=:scatter)

X = UniformSplineDesignMat(
    Float32.(coalesce.(outf[!,:x0], 0.0)), 
    3, 
    10
)
IDtoCOL = [ArrayDict(UInt32, UInt16, n_precursors*2 +1 ) for _ in range(1, N)];
function solve_constrained_least_squares_regression(A::Matrix, X0::Matrix, X1::Matrix, b::Vector)
    m, n = size(A)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:n])
    @variable(model, residuals[1:m])
    @constraint(model, residuals == A * x - b)
    @constraint(model, X0*x .<= 0)
    @constraint(model, X1*x .<= 0)
    @objective(model, Min, sum(residuals.^2))
    optimize!(model)
    return JuMP.value.(x)
end
function solve_constrained_thing(yt::Vector{Float32}, x0::Vector{Float32}, x1::Vector{Float32}, b::Float32)
    m = length(yt)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, b)
    @variable(model, residuals[1:m])
    @constraint(model, residuals == ((1 .+ abs.(x1 ./1 ).^(2*b))./(1 .+ abs.(x0./1).^(2*b))) .- yt) 
    @constraint(model, b >= 0)
    @objective(model, Min, sum(residuals.^2))
    optimize!(model)
    return JuMP.value.(b)
end
@time begin 
@test abs(solve_constrained_thing(exp2.(yt), x0, x1, 2.0f0) - 5.0) < (1e-3)
@test abs(solve_constrained_thing(exp2.(yt), x0, x1, 1.0f0) - 5.0) < (1e-3)
@test abs(solve_constrained_thing(exp2.(yt), x0, x1, 0.5f0) - 5.0) < (1e-3)
@test abs(solve_constrained_thing(exp2.(yt), x0, x1, 100.0f0) - 5.0) < (1e-3)
end

function solve_constrained_thing_log(yt::Vector{Float32}, x0::Vector{Float32}, x1::Vector{Float32}, b::Float32)
    m = length(yt)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, b)
    @variable(model, residuals[1:m])
    @constraint(model, residuals == (log2.(1 .+ abs.(x1 ./1 ).^(2*b)) - log2.(1 .+ abs.(x0./1).^(2*b))) .- yt) 
    @constraint(model, b >= 0)
    @objective(model, Min, sum(residuals.^2))
    optimize!(model)
    return JuMP.value.(b)
end


@time begin 
@test abs(solve_constrained_thing_log(yt, x0, x1, 2.0f0) - 5.0) < (1e-3)
@test abs(solve_constrained_thing_log(yt, x0, x1, 1.0f0) - 5.0) < (1e-3)
@test abs(solve_constrained_thing_log(yt, x0, x1, 0.5f0) - 5.0) < (1e-3)
@test abs(solve_constrained_thing_log(yt, x0, x1, 100.0f0) - 5.0) < (1e-3)
end

=#

#X = X[:,1:]
#@time c = solve_constrained_least_squares_regression(X, u)


A, b = rand(1000, 8), rand(1000);
@time x = solve_constrained_least_squares_regression(A, b)



#=
IDtoCOL = [ArrayDict(UInt32, UInt16, n_precursors*2 +1 ) for _ in range(1, N)];
 precursor_weights = [zeros(Float32, n_precursors*2 + 1 ) for _ in range(1, N)];
test_qpsms = getQuadT(
        1500.0f0,
        0,
        precursor_dict,
        MS_TABLE_PATHS,
        precursors[:is_decoy],
        frag_err_dist_dict,
        rt_index_paths,
        bin_rt_size,
        rt_irt,
        irt_errs,
        chromatograms,
        file_path_to_parsed_name,
        params_,
        spec_lib,
        ionMatches,
        ionMisses,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        complex_scored_PSMs,
        complex_unscored_PSMs,
        complex_spectral_scores,
        precursor_weights);

sort!(test_qpsms, [:ms_file_idx, :scan_idx,:precursor_idx,:iso_idx])
test_qpsms[!,:scan_idx]
filter!(x->!iszero(x.weight), test_qpsms)
test_qpsms[!,:mono_mz] = [precursors[:mz][pid] for pid in test_qpsms[!,:precursor_idx]]
test_qpsms[!,:prec_charge] = [precursors[:prec_charge][pid] for pid in test_qpsms[!,:precursor_idx]]
test_qpsms[!,:sulfur_count] = [precursors[:sulfur_count][pid] for pid in test_qpsms[!,:precursor_idx]]

test_qpsms[!,:iso_mz] = test_qpsms[!,:mono_mz] .+ NEUTRON.*test_qpsms[!,:iso_idx]./test_qpsms[!,:prec_charge]
test_qpsms[!,:x] = test_qpsms[!,:iso_mz] .- test_qpsms[!,:center_mz]
test_qpsms[!,:δ] .= zero(Float32)
for i in range(1, size(test_qpsms, 1))
    s_count = min(Int64(test_qpsms[i,:sulfur_count]), 5)
    mono_mass =  test_qpsms[i,:mono_mz]*test_qpsms[i,:prec_charge]
    test_qpsms[i,:δ] = iso_splines(s_count, 0, mono_mass)/iso_splines(s_count, 1, mono_mass)
end
histogram(test_qpsms[!,:δ])

gqpsms = groupby(test_qpsms, [:ms_file_idx, :scan_idx,:precursor_idx])
outf = combine(gqpsms) do psms
    if size(psms, 1) == 2
       return (center_mz = psms[1,:center_mz],
               min_weight = minimum(psms[!,:weight]),
               max_weight = maximum(psms[!,:weight]),
       δ = psms[1,:δ], yt = log2((1/psms[1,:δ])*psms[1,:weight]/psms[2,:weight]), x0 = psms[1,:x], x1 = psms[2,:x])
    else
        return (center_mz = missing, 
        min_weight = missing,
        max_weight = missing,
        δ = missing, yt = missing, x0 = missing, x1 = missing)
    end
end
filter!(x->!ismissing(x.yt), outf)
filter!(x->!isinf(x.yt), outf)
filter!(x->!isnan(x.yt), outf)
filter!(x->x.max_weight>10000.0f0, outf)


yt = Float32.(coalesce.(outf[!,:yt], 0.0))
x0 = Float32.(coalesce.(outf[!,:x0], 0.0))
x1 = Float32.(coalesce.(outf[!,:x1], 0.0))
yt = append!(yt, ones(Float32, 100).*10000.0f0)
x0 = append!(x0, zeros(Float32, 100))
x1 = append!(x1, ones(Float32, 100).*2.0f0)
yt = append!(yt,  ones(Float32, 100).*0.00001f0)
x0 = append!(x0,  ones(Float32, 100).*(-2.0f0))
x1 = append!(x1, zeros(Float32, 100))

gx = UniformSplineDesignMat(
    yt, x0, x1,
    3, 
    5
)

#plot(LinRange(-4, 4, 100), gx.(LinRange(-4, 4, 100)))
norm_fac = maximum(gx.(LinRange(-4, 4, 100)))
y = gx.(LinRange(-4, 4, 100))./norm_fac .- 1
plot(LinRange(-4, 4, 100), exp2.(y))
hline!([0.0])
vline!([-1.5, 1.5])

plot(LinRange(-4, 4, 100), 
    gx.(LinRange(-4, 4, 100) .+ gx.(LinRange(-4, 4, 100).+))
    )


plot(LinRange(-2, 2, 100), exp2.(gx.(LinRange(-2, 2, 100))))
vline!([-1, 1])

X1 = UniformSplineDesignMat(
    Float32.(coalesce.(outf[!,:x1], 0.0)), 
    3, 
    10
)

A = X0 .+ X1
c = X\outf[!,:yt]


plot(x0, exp2.(yt), seriestype=:scatter, alpha = 0.1)

gx, X, u = UniformSplineDesignMat(
    yt, x0, x1,
    3, 
    5
)
plot(LinRange(-4, 4, 100), gx.(LinRange(-4, 4, 100)))
y = gx.(LinRange(-4, 4, 100))
norm_fac = maximum(y)
plot(LinRange(-4, 4, 100), exp2.(y./norm_fac .- 1))
plot!(LinRange(-4, 4, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(-4, 4, 100))))
vline!([-1.5, 1.5, -1.0, 1.0])


plot(LinRange(-4, 4, 100), exp2.(gx.(LinRange(-4, 4, 100))))

#plot(LinRange(-4, 4, 100), gx.(LinRange(-4, 4, 100)))
#norm_fac = maximum(gx.(LinRange(-4, 4, 100)))
#y = gx.(LinRange(-4, 4, 100))./norm_fac .- 1
#plot(LinRange(-4, 4, 100), exp2.(y))
#hline!([0.0])
#vline!([-1.5, 1.5])



=#

function getQuadT(
    δ::Float32,
    isotope_idx::Int64,
    gbpsms,
    MS_TABLE_PATHS,
    frag_err_dist_dict,
    rt_index_paths,
    bin_rt_size,
    rt_irt,
    irt_errs,
    chromatograms,
    file_path_to_parsed_name,
    params_,
    spec_lib,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    complex_scored_PSMs,
    complex_unscored_PSMs,
    complex_spectral_scores,
    precursor_weights;
    minimum_percent_diff::Float32 = 10.0f0)
    psms = []

    for (key, sub_bpsms) in ProgressBar(pairs(gbpsms))
        ms_file_idx = key[:ms_file_idx]
        MS_TABLE_PATH = MS_TABLE_PATHS[ms_file_idx]
        prec_set = Set(zip(sub_bpsms[!,:precursor_idx], sub_bpsms[!,:scan_idx]))
        scan_idxs = Set(sub_bpsms[!,:scan_idx])

        push!(psms, estimateQuadTransmissionSearch(
            isotope_idx,
            frag_err_dist_dict,
            rt_index_paths,
            prec_set,
            scan_idxs,
            δ,
            bin_rt_size,
            rt_irt,
            irt_errs,
            chromatograms,
            file_path_to_parsed_name,
            ms_file_idx,
            MS_TABLE_PATH,
            params_,
            spec_lib,
            ionMatches,
            ionMisses,
            IDtoCOL,
            ionTemplates,
            iso_splines,
            complex_scored_PSMs,
            complex_unscored_PSMs,
            complex_spectral_scores,
            precursor_weights
            ));
    end
    psms = vcat(psms...);

    return psms
end

function getQuadT(
    δ::Float32,
    isotope_idx::Int64,
    prec_to_irt::Dictionary{UInt32, 
    @NamedTuple{best_prob::Float32, 
                best_ms_file_idx::UInt32, 
                best_scan_idx::UInt32, 
                best_irt::Float32, 
                mean_irt::Union{Missing, Float32}, 
                var_irt::Union{Missing, Float32}, 
                n::Union{Missing, UInt16},
                mz::Float32}},
    MS_TABLE_PATHS::Vector{String},
    is_decoy::AbstractVector{Bool},
    frag_err_dist_dict,
    rt_index_paths,
    bin_rt_size,
    rt_irt,
    irt_errs,
    chromatograms,
    file_path_to_parsed_name,
    params_,
    spec_lib,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    complex_scored_PSMs,
    complex_unscored_PSMs,
    complex_spectral_scores,
    precursor_weights)

    best_psms = getPsmsForHuberEstimation(prec_to_irt,is_decoy)
    gbpsms = groupby(best_psms,:ms_file_idx)
    return getQuadT(
        δ,
        isotope_idx,
        gbpsms,
        MS_TABLE_PATHS,
        frag_err_dist_dict,
        rt_index_paths,
        bin_rt_size,
        rt_irt,
        irt_errs,
        chromatograms,
        file_path_to_parsed_name,
        params_,
        spec_lib,
        ionMatches,
        ionMisses,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        complex_scored_PSMs,
        complex_unscored_PSMs,
        complex_spectral_scores,
        precursor_weights
        )
end



function estimateQuadTransmissionSearch(
    isotope_idx,
    frag_err_dist_dict,
    rt_index_paths,
    prec_set,
    scan_idxs,
    δ,
    bin_rt_size,
    rt_irt,
    irt_errs,
    chromatograms,
    file_path_to_parsed_name,
    ms_file_idx,
    MS_TABLE_PATH,
    params_,
    spec_lib,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    complex_scored_psms,
    complex_unscored_psms,
    complex_spectral_scores,
    precursor_weights
    )

    function estimateQuadTransmissionSearch(
                        #Mandatory Args
                        spectra::Arrow.Table, 
                        params::Any,
                        prec_set::Set{Tuple{UInt32, UInt32}},
                        scan_idxs::Set{UInt32},
                        δ::Float32,
                        isotope_idx::Int64;
                        kwargs...)
        ########
        #Each thread needs to handle a similair number of peaks. 
        #For example if there are 10,000 scans and two threads, choose n so that
        #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
        #of scans have an equal number of fragment peaks in the spectra
        thread_tasks, total_peaks = partitionScansToThreads(spectra[:mz_array],
                                                            spectra[:retentionTime],
                                                            spectra[:centerMz],
                                                            spectra[:msOrder],

                                                            Threads.nthreads(),
                                                            1)
        
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin 
                thread_id = first(thread_task)
                return QuadTransmissionSearch(
                                    spectra,
                                    last(thread_task), #getRange(thread_task),
                                    prec_set,
                                    scan_idxs,
                                    spec_lib["precursors"],
                                    kwargs[:fragment_lookup_table], 
                                    kwargs[:ms_file_idx],
                                    kwargs[:rt_to_irt_spline],
                                    kwargs[:mass_err_model],
                                    kwargs[:quad_transmission_func],
                                    δ,
                                    Float32(params[:deconvolution_params]["lambda"]),
                                    Int64(params[:deconvolution_params]["max_iter_newton"]),
                                    Int64(params[:deconvolution_params]["max_iter_bisection"]),
                                    Int64(params[:deconvolution_params]["max_iter_outer"]),
                                    Float32(params[:deconvolution_params]["accuracy_newton"]),
                                    Float32(params[:deconvolution_params]["accuracy_bisection"]),
                                    Float32(params[:deconvolution_params]["max_diff"]),
                                    kwargs[:ion_matches][thread_id],
                                    kwargs[:ion_misses][thread_id],
                                    kwargs[:id_to_col][thread_id],
                                    kwargs[:ion_templates][thread_id],
                                    kwargs[:iso_splines],
                                    kwargs[:unscored_psms][thread_id],
                                    kwargs[:spectral_scores][thread_id],
                                    kwargs[:precursor_weights][thread_id],
                                    (3, 1),
                                    isotope_idx,#params[:quant_search_params]["n_frag_isotopes"],
                                    kwargs[:rt_index], 
                                    kwargs[:irt_err],
                                    Set(2),
                                )
            end
        end
        psms = fetch.(tasks)
        return psms
    end

    ms_table_path_to_psms_path = Dict{String, String}()
    parsed_fname = file_path_to_parsed_name[MS_TABLE_PATH]
    rt_df = DataFrame(Arrow.Table(rt_index_paths[parsed_fname]))
    rt_index = buildRtIndex(rt_df,
                            bin_rt_size = bin_rt_size)
    rt_irt = rt_irt[parsed_fname]
    MS_TABLE = Arrow.Table(MS_TABLE_PATH);
        params_[:quant_search_params]["min_y_count"] = 1
        precursors = spec_lib["precursors"]
        psms = vcat(estimateQuadTransmissionSearch(
            MS_TABLE, 
            params_,
            prec_set,
            scan_idxs,
            δ,
            isotope_idx;
            precursors = spec_lib["precursors"],
            fragment_lookup_table = spec_lib["f_det"],
            rt_index = rt_index,
            ms_file_idx = UInt32(ms_file_idx), 
            rt_to_irt_spline = rt_irt,
            mass_err_model = frag_err_dist_dict[ms_file_idx],
            irt_err = irt_errs[parsed_fname],#irt_errs[ms_file_idx]/3,
            ion_matches = ionMatches,
            ion_misses = ionMisses,
            id_to_col = IDtoCOL,
            ion_templates = ionTemplates,
            iso_splines = iso_splines,
            chromatograms = chromatograms,
            scored_psms = complex_scored_psms,
            unscored_psms = complex_unscored_psms,
            spectral_scores = complex_spectral_scores,
            precursor_weights = precursor_weights,
            quad_transmission_func = QuadTransmission(10.0f0, 10000.0f0)
            )...);
    psms[!,:ms_file_idx] .= UInt32(ms_file_idx)
    return psms
end

