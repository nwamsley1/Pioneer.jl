mzvals = copy(precursors[:mz])
filter!(x->x<604, mzvals)
filter!(x->x>598, mzvals)
histogram(mzvals, bins = LinRange(598, 604, 1000))


outf[!,:prec_charge] = [precursors[:prec_charge][pid] for pid in outf[!,:precursor_idx]]
#filter!(x->coalesce(x.x0, -Inf)>(0.0), outf)
outf2 = copy(outf[outf[!,:prec_charge].==2,:])
outf3 = copy(outf[outf[!,:prec_charge].==3,:])

add_equal_width_bins!(outf2, :x0, 25)
add_equal_width_bins!(outf3, :x0, 25)
outf = vcat([outf2, outf3]...)
comb_df = combine(groupby(outf,[:prec_charge,:x0_bin])) do outbin 
    return (n = size(outbin, 1), 
            median_yt = median(outbin[!,:yt]),
            mean_yt = mean(outbin[!,:yt]),
            median_x0 = median(outbin[!,:x0]),
            prec_charge = outbin.prec_charge[1]
            )
end
filter!(x->x.n>50,comb_df)

plot(comb_df[!,:median_x0], comb_df[!,:median_yt], seriestype=:scatter, alpha = 1.0)
yt = comb_df[!,:median_yt]
x0 = Float32.(comb_df[!,:median_x0])
x1 = Float32.(comb_df[!,:median_x0] .+ NEUTRON./comb_df[!,:prec_charge])
#=
bin_to_val = Dict(zip(comb_df[!,:x0_bin], comb_df[!,:median_yt]))
outf[!,:yt_corrected] .= -1.0f0
for i in range(1, size(outf, 1))
    bin_id = outf[i,:x0_bin]
    if bin_id  ∈ keys(bin_to_val)
        outf[i,:yt_corrected] = bin_to_val[bin_id]
    end
end
filter!(x->x.yt_corrected>-1.0, outf)
filter!(x->x.yt_corrected>-1.0, outf)

yt, x0, x1 = Float32.(coalesce.(outf[!,:yt_corrected], 0.0f0)), Float32.(coalesce.(outf[!,:x0], 0.0f0)),Float32.(coalesce.(outf[!,:x1], 0.0f0))

yt, x0, x1 = Float32.(coalesce.(outf[!,:yt], 0.0f0)), Float32.(coalesce.(outf[!,:x0], 0.0f0)),Float32.(coalesce.(outf[!,:x1], 0.0f0))

plot(x0, yt, seriestype=:scatter, alpha = 0.1)
=#


tspline = UniformSpline(yt,x0,3,5)
plot(x0, yt, seriestype=:scatter)
plot!(LinRange(-2, 2, 100), tspline.(LinRange(-2, 2, 100)))


tspline = UniformSpline(max.(outf[!,:yt], 0.0f0),outf[!,:x0],3,10)
plot(outf[!,:x0], max.(outf[!,:yt], 0.0f0), seriestype=:scatter, xlim = (0, 2), ylim = (-1, 5), alpha = 0.02)
plot!(LinRange(-2, 2, 100), tspline.(LinRange(-2, 2, 100)))


[x0.>0]
gx = UniformSplineDesignMat(
    yt, x0, x1,
    3, 
    5
)
#plot(LinRange(-4, 4, 100), gx.(LinRange(-4, 4, 100)))
log2y = gx.(LinRange(0, 4, 100))
y = exp2.(log2y .- maximum(log2y))
plot(LinRange(0, 4, 100), y)
qtf = QuadTransmission(0.4f0, 5.0f0)
plot!(LinRange(0, 4, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(0, 4, 100))))


qtf = QuadTransmission(0.5f0, 2.0f0)
plot!(LinRange(-4, 4, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(-4, 4, 100))))



ind_to_keep = x0 .> -0.5f0
gx = UniformSplineDesignMat(
    yt[ind_to_keep], x0[ind_to_keep], x1[ind_to_keep],
    3, 
    5
)
gx2 = UniformSplineDesignMat(
    yt[ind_to_keep], x0[ind_to_keep], x1[ind_to_keep],
    3, 
    30
)
#plot(LinRange(-4, 4, 100), gx.(LinRange(-4, 4, 100)))
log2y = gx.(LinRange(0, 4, 100))
y = exp2.(log2y .- maximum(log2y))
plot(LinRange(0, 4, 100), y)
log2y2= gx2.(LinRange(0, 4, 100))
y2 = exp2.(log2y2 .- maximum(log2y2))
plot!(LinRange(0, 4, 100), y2)
plot!(LinRange(0, 4, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(0, 4, 100))))



plot(LinRange(0, 4, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(0, 4, 100))))
test_tq = TQuad(1.0, 2.0, 2.0, 2.0, 0.0)
plot_bins = LinRange(-3, 3, 100)
plot!(
    plot_bins,
    test_tq.(plot_bins)
)

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

include("utils/isotopeSplines.jl")
test_qpsms = getQuadT(
        100000.0f0,
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
    iso_probs = [iso_splines(s_count, iso_idx, mono_mass) for iso_idx in range(0, 5)]
    #iso_probs = iso_probs./sum(iso_probs)
    test_qpsms[i,:δ] = iso_probs[1]/iso_probs[2]
end
histogram(test_qpsms[!,:δ])

gqpsms = groupby(test_qpsms, [:ms_file_idx, :scan_idx,:precursor_idx])
outf = combine(gqpsms) do psms
    if size(psms, 1) == 2
        sort!(psms, :iso_idx)
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
#filter!(x->x.x0>0.0, outf)
#plot(outf[!,:x0], outf[!,:yt], seriestype=:scatter, ylim = (-1, 5), alpha = 0.02)

outf[!,:prec_charge] = [precursors[:prec_charge][pid] for pid in outf[!,:precursor_idx]];
#plot(outf[outf[!,:prec_charge].==2,:x0], outf[outf[!,:prec_charge].==2,:yt], seriestype=:scatter, ylim = (-1, 5), alpha = 0.02)
#plot!(outf[outf[!,:prec_charge].==3,:x0], outf[outf[!,:prec_charge].==3,:yt], seriestype=:scatter, ylim = (-1, 5), alpha = 0.02)


outf2 = copy(outf[outf[!,:prec_charge].==2,:])
outf3 = copy(outf[outf[!,:prec_charge].==3,:])

add_equal_width_bins!(outf2, :x0, 25)
add_equal_width_bins!(outf3, :x0, 25)
outf = vcat([outf2, outf3]...)
comb_df = combine(groupby(outf,[:prec_charge,:x0_bin])) do outbin 
    return (n = size(outbin, 1), 
            median_yt = median(outbin[!,:yt]),
            mean_yt = mean(outbin[!,:yt]),
            median_x0 = median(outbin[!,:x0]),
            prec_charge = outbin.prec_charge[1]
            )
end
filter!(x->x.n>50,comb_df)
p = plot(
 ylabel = "log2(f(z0)/f(z1))",
 xlabel = "z0")
plot!(p, comb_df[comb_df[!,:prec_charge].==2,:median_x0], comb_df[comb_df[!,:prec_charge].==2,:median_yt], seriestype=:scatter,
label = "charge +2")
plot!(p, comb_df[comb_df[!,:prec_charge].==3,:median_x0], comb_df[comb_df[!,:prec_charge].==3,:median_yt], seriestype=:scatter,
label = "charge +3")


qtf = QuadTransmission(0.5f0, 7.0f0)
quad_results = simmulateQuad(
    qtf,
    1.0f0,
    (-2.0f0, 2.0f0),
    (400.0f0, 1000.0f0),
    1000
)
sort!(quad_results,:x0)
plot!(p,
quad_results[!,:x0],
quad_results[!,:yt],
label = "simmulated_data"
)
hline!([-1.0, 0.0, 1.0])

filter!(x->x.max_weight>1000.0f0, outf)


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
