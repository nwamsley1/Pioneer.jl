IDtoCOL = [ArrayDict(UInt32, UInt16, n_precursors*3 +1 ) for _ in range(1, N)];
precursor_weights = [zeros(Float32, n_precursors*3 + 1 ) for _ in range(1, N)];

include("utils/isotopeSplines.jl")

global global_varable = nothing
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

filter!(x->x.iso_idx<3, test_qpsms)
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
    test_qpsms[i,:δ] = (iso_splines(s_count, 0, mono_mass)/iso_splines(s_count, 1, mono_mass))
end
histogram(test_qpsms[!,:δ])

gqpsms = groupby(test_qpsms, [:ms_file_idx,:scan_idx,:precursor_idx])
outf = combine(gqpsms) do psms
    if size(psms, 1) == 2
       sort!(psms, :iso_idx)
       s_count = min(Int64(psms[1,:sulfur_count]), 5)
       mono_mass = psms[1,:mono_mz]*psms[1,:prec_charge]
            δ =  (iso_splines(s_count, 0, mono_mass)/iso_splines(s_count, 1, mono_mass))
            return (center_mz = psms[1,:center_mz],
                    min_weight = minimum(psms[!,:weight]),
                    max_weight = maximum(psms[!,:weight]),
                    δ = δ, yt = log2((1/δ)*psms[1,:weight]/psms[2,:weight]), x0 = psms[1,:x], x1 = psms[2,:x])
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
filter!(x->x.max_weight > 1000.0f0, outf)
#filter!(x->x.x0>0.0, outf)
#plot(outf[!,:x0], outf[!,:yt], seriestype=:scatter, ylim = (-1, 5), alpha = 0.02)

outf[!,:prec_charge] = [precursors[:prec_charge][pid] for pid in outf[!,:precursor_idx]];
plot(outf[outf[!,:prec_charge].==2,:x0], outf[outf[!,:prec_charge].==2,:yt], seriestype=:scatter, ylim = (-5, 5), alpha = 0.02)
plot!(outf[outf[!,:prec_charge].==3,:x0], outf[outf[!,:prec_charge].==3,:yt], seriestype=:scatter, ylim = (-5, 5), alpha = 0.02)


outf2 = copy(outf[(outf[!,:prec_charge].==2).&((outf[!,:center_mz].>800.0)),:])
outf3 = copy(outf[(outf[!,:prec_charge].==3).&((outf[!,:center_mz].>800.0)),:])

add_equal_width_bins!(outf2, :x0, 50)
add_equal_width_bins!(outf3, :x0, 50)
outf = vcat([outf2, outf3]...)
comb_df = combine(groupby(outf,[:prec_charge,:x0_bin])) do outbin 
    return (n = size(outbin, 1), 
            median_yt = median(outbin[!,:yt]),
            mean_yt = mean(outbin[!,:yt]),
            median_x0 = median(outbin[!,:x0]),
            prec_charge = outbin.prec_charge[1]
            )
end
filter!(x->x.n>25,comb_df)
comb_df[!,:median_x1] = comb_df[!,:median_x0] .+ NEUTRON./comb_df[!,:prec_charge]
#comb_df = comb_df[comb_df[!,:median_x0].>=0.0f0, :]
#comb_df2 = copy(comb_df)
#comb_df2[!,:median_x0] = comb_df2[!,:median_x0].*(-1.0) .- 0.5f0
#comb_df2[!,:median_x1] = comb_df2[!,:median_x1].*(-1.0) .- 0.5f0
#comb_df2[!,:median_yt] = comb_df2[!,:median_yt].*(-1.0)
#comb_df = vcat([comb_df, comb_df2]...)

p = plot(
 ylabel = "log2(f(z0)/f(z1))",
 xlabel = "z0")
plot!(p, comb_df[comb_df[!,:prec_charge].==2,:median_x0], comb_df[comb_df[!,:prec_charge].==2,:median_yt], seriestype=:scatter,
label = "charge +2")
plot!(p, comb_df[comb_df[!,:prec_charge].==3,:median_x0], comb_df[comb_df[!,:prec_charge].==3,:median_yt], seriestype=:scatter,
label = "charge +3")

x03, yt3 =  Float32.(comb_df[comb_df[!,:prec_charge].==3,:median_x0]), Float32.(comb_df[comb_df[!,:prec_charge].==3,:median_yt]);
gl = UniformSpline(yt3, x03, 3, 15)#linear_interpolation(x03,yt3,extrapolation_bc=Flat());
x03 =Float32.(LinRange(minimum(x03), maximum(x03), 500));
x13 = x03 .+ NEUTRON/3;
y03 = gl.(x03);

x02, yt2 =  Float32.(comb_df[comb_df[!,:prec_charge].==2,:median_x0]), Float32.(comb_df[comb_df[!,:prec_charge].==2,:median_yt]);
gl = UniformSpline(yt2, x02, 3, 15)#linear_interpolation(x02,yt2,extrapolation_bc=Flat());
x02 =Float32.(LinRange(minimum(x02), maximum(x02), 500));
x12 = x02 .+ NEUTRON/3;
y02 = gl.(x02);


x0 = vcat([x02, x03]...);
x1 = vcat([x12, x13]...);
yt = vcat([y02, y03]...);
idx = (x0.>-2.0).&(x0.<2.0)#typemin(Float32)
gx = UniformSplineDesignMat(
    Float32.(yt[idx]), Float32.(x0[idx]), Float32.(x1[idx]),
    3, 
    15#,λ = 10.0
)
log2y = gx.(LinRange(-3, 3, 100))
y = exp2.(log2y .- maximum(log2y))
plot(LinRange(-3, 3, 100), y)





comb_df[comb_df[!,:prec_charge].==3,:median_x0], comb_df[comb_df[!,:prec_charge].==3,:median_yt]

x0, x1, yt = copy(comb_df[!,:median_x0]), copy(comb_df[!,:median_x1]), copy(comb_df[!,:median_yt])
idx = (x0.>-2.0).&(x0.<2.0)#typemin(Float32)
gx = UniformSplineDesignMat(
    Float32.(yt[idx]), Float32.(x0[idx]), Float32.(x1[idx]),
    3, 
    15
)
log2y = gx.(LinRange(-3, 3, 100))
y = exp2.(log2y .- maximum(log2y))
plot!(LinRange(-3, 3, 100), y)



gx = UniformSplineDesignMat(
    Float32.(yt[idx]), Float32.(x0[idx]), Float32.(x1[idx]),
    3, 
    10
)
log2y = gx.(LinRange(-3, 3, 100))
y = exp2.(log2y .- maximum(log2y))
plot!(LinRange(-3, 3, 100), y)





x0, x1, yt = Float32.(outf[!,:x0]), Float32.(outf[!,:x1]), Float32.(outf[!,:yt])
idx = (x0.>-2.0).&(x0.<2.0)#typemin(Float32)
gx = UniformSplineDesignMat(
    Float32.(yt[idx]), Float32.(x0[idx]), Float32.(x1[idx]),
    3, 
    5
)
log2y = gx.(LinRange(-3, 3, 100))
y = exp2.(log2y .- maximum(log2y))
plot!(LinRange(-3, 3, 100), y)




df2 = comb_df[comb_df[!,:prec_charge].==2,:]
filter!(x->x.median_x0>-1.0, df2)
#gx = UniformSpline(
#    df2[!,:median_yt],Float32.(df2[!,:median_x0]),
#    3, 10
#)
#y = gx.(LinRange(-3, 3, 100))
gl = linear_interpolation(
    Float32.(df2[!,:median_x0]),df2[!,:median_yt],extrapolation_bc=Line())
y = gl.(LinRange(-3, 3, 100))
plot(df2[!,:median_x0], df2[!,:median_yt], seriestype=:scatter)
plot!(LinRange(-3, 3, 100), y)

x0 = Float32.(LinRange(-2, 2, 1000))
x1 = x0 .+ 0.5f0
yt = Float32.(gl.(x0))
gx = UniformSplineDesignMat(
    yt, x0, x1,
    3, 
    10
)
log2y = gx.(LinRange(-3, 3, 100))
y = exp2.(log2y .- maximum(log2y))
plot(LinRange(-3, 3, 100), y)
qtf = QuadTransmission(0.75f0, 5.0f0)
plot!(LinRange(-3, 3, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(-3, 3, 100))))


p = plot(
 ylabel = "log2(f(z0)/f(z1))",
 xlabel = "z0")
plot!(p, comb_df[comb_df[!,:prec_charge].==2,:median_x0], comb_df[comb_df[!,:prec_charge].==2,:median_yt], seriestype=:scatter,
label = "charge +2")
plot!(p, comb_df[comb_df[!,:prec_charge].==3,:median_x0], comb_df[comb_df[!,:prec_charge].==3,:median_yt], seriestype=:scatter,
label = "charge +3")


yt = Float32.(coalesce.(outf[!,:yt], 0.0f0))
x0 = Float32.(coalesce.(outf[!,:x0], 0.0f0))
x1 = Float32.(coalesce.(outf[!,:x1], 0.0f0))
gx = UniformSplineDesignMat(
    Float32.(yt), Float32.(x0), Float32.(x1),
    3, 
    20
)

log2y = gx.(LinRange(-3, 3, 100))
y = exp2.(log2y .- maximum(log2y))
plot(LinRange(-3, 3, 100), y)
qtf = QuadTransmission(0.75f0, 5.0f0)
plot!(LinRange(-3, 3, 100), qtf.(0.0f0, 1.0f0, Float32.(LinRange(-3, 3, 100))))





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

Hs = SparseArray{UInt32, Float32}(21, 14, 2, UInt32[0x00000001, 0x00000003, 0x00000005, 0x00000007, 0x00000009, 0x0000000b, 0x0000000d, 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e], UInt16[0x0001, 0x0001, 0x0001, 0x0001, 0x0001, 0x0001, 0x0001, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002, 0x0002], Float32[15.546875, 8.828125, 27.5625, 30.546875, 23.828125, 19.8125, 15.0234375, 15.765625, 1.7646484, 7.984375, 1.9824219, 23.625, 7.4921875, 24.03125, 10.4453125, 17.828125, 9.0625, 13.59375, 8.7578125, 8.9140625, 8.0234375], Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], UInt8[0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01], Float32[19070.264, 7259.259, 17350.783, 21396.404, 17471.098, 17425.406, 18178.537, 19070.264, 1839.0159, 7259.259, 1562.0409, 17350.783, 4562.439, 21396.404, 8440.88, 17471.098, 5573.7974, 17425.406, 5284.5195, 18178.537, 6727.6846], UInt32[0x00000001, 0x00000008, 0x00000015])
N = Hs.n_vals

#Get Data
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end
#Design Matrix 
H = Matrix(sparse(Hs.rowval[1:N],
Hs.colval[1:N],
Hs.nzval[1:N]))
#OLS Regression 
rowvals = copy(Hs.rowval)
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end

hcat(H, y, H*(H\y))
plot(y, H*(H\y), seriestype=:scatter)

hcat(H[[2, 4, 6, 8, 10, 12, 14],:], y[[2, 4, 6, 8, 10, 12, 14]])


H[[2, 4, 6, 8, 10, 12, 14],:]\y[[2, 4, 6, 8, 10, 12, 14]]



H[[1, 3, 5, 7, 9, 11, 13],:]\y[[1, 3, 5, 7, 9, 11, 13]]