function getAreaRatios(
    second_quant_path::String
    )
    #To get peak areas for seperate traces go to "mergePSMTables.jl"
    #Go to sortAndFilterQuandTables
    #comment ount `filter!(x->x.best_trace, psms_table)`
    quant_table = DataFrame(Tables.columntable(Arrow.Table(second_quant_path)))
    quant_table[!,:first_iso] = [first(x) for x in quant_table[!,:isotopes_captured]]
    sort!(quant_table, :first_iso)
    gquant = groupby(quant_table,:precursor_idx)
    area_ratios = Float32[]
    for (precursor_idx, psms) in pairs(gquant)
        if size(psms, 1) != 2
            continue
        end
        push!(area_ratios, log2(psms[1,:peak_area]/psms[2,:peak_area]))
    end
    return area_ratios
end

ar_fitted = getAreaRatios("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/ALTERNATING_WINDOW_TEST/NOV20_TEST_ECLIPSE/arrow_out/SeperateTraces_FittedQuad_Astral/RESULTS/temp/second_quant/03.arrow");
ar_default = getAreaRatios("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/ALTERNATING_WINDOW_TEST/NOV20_TEST_ECLIPSE/arrow_out/SeperateTraces_DefaultQuad_Astral/RESULTS/temp/second_quant/03.arrow");
ar_none = getAreaRatios("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/ALTERNATING_WINDOW_TEST/NOV20_TEST_ECLIPSE/arrow_out/SeperateTraces_NoQuad_Astral/RESULTS/temp/second_quant/03.arrow");
pbins = LinRange(-3, 3, 100)
histogram(ar_fitted, alpha = 0.5, bins = pbins, normalize = :probability)
histogram!(ar_default, alpha = 0.5, bins = pbins, normalize = :probability)
histogram!(ar_none, alpha = 0.5, bins = pbins, normalize = :probability)

p = plot(xlabel = "log2 Ratio of Isotope Trace Peak Areas", legend = :topleft)
histogram!(p, ar_fitted, alpha = 0.5, bins = pbins, color = 1, label = "Quad Fitted")
vline!([median(ar_fitted)], color = :black, lw = 4, label = nothing)
histogram!(p, ar_default, alpha = 0.5, bins = pbins, label = "Default Quad Model")
vline!([median(ar_default)], color = :black, lw = 4, label = nothing)
histogram!(p, ar_none, alpha = 0.5, bins = pbins, label = "No Fragment Isotope \n Correction")
vline!([median(ar_none)], color = :black, lw = 4, label = nothing)
vline!([0.0], label = "0", color = :black, lw = 4, linestyle = :dot)
savefig(p, "/Users/n.t.wamsley/Documents/ThermoMeeting_112524/IsotopeTraceRatios.pdf")





chroms_default = sort(load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/ALTERNATING_WINDOW_TEST/NOV20_TEST_ECLIPSE/arrow_out/CombinedTraces_DefaultQuad_Eclipse/RESULTS/03testchroms.jld2")["chroms"],:scan_idx);
chroms_no_quad = sort(load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/ALTERNATING_WINDOW_TEST/NOV20_TEST_ECLIPSE/arrow_out/CombinedTraces_NoQuad_Eclipse/RESULTS/03testchroms.jld2")["chroms"],:scan_idx);

gchroms_default = groupby(chroms_default, :precursor_idx);
gchroms_noquad = groupby(chroms_no_quad, :precursor_idx);


precs_to_keep =  [
3135834,
 2309582,
  3178163,
   2743662,
   3188585,
     3188636,
     3190697,
     3194879,
      2757932,
       2758057,
       2759938,
       2761811,
       3209938,
       3676643,
        2771567,
          3219743,
           3221227  ,
             3223309
]


key = (precursor_idx = precs_to_keep[n], )
default = gchroms_default[key]
no_quad = gchroms_noquad[key]
plot(default[!,:rt], default[!,:intensity], seriestype=:scatter, alpha = 0.5, label = "Uncorrected")
plot!(no_quad[!,:rt], no_quad[!,:intensity], seriestype=:scatter, alpha = 0.5, label = "Corrected")
n += 1


λ::Float32 = 1.0f0
N = 20 + 26
b = zeros(Float32, N);
A = getWittakerHendersonDesignMat(length(b), λ);
prob = LinearProblem(A, b);
linsolve = init(prob);
u2 = zeros(Float32, length(linsolve.b));
dtype = Float32
state = Chromatogram(
    zeros(dtype, N), #t
    zeros(dtype, N), #data
    N #max index
    )


apex_scan = findfirst(
        x->x==v2[!,:scan_idx][argmax(v2[!,:intensity])], v2[!,:scan_idx])
peak_area, _ = integrateChrom(
v2,
apex_scan,
linsolve,u2,state, n_pad = 20, max_apex_offset = 1,
isplot = true
)

plot(state.t[1:state.max_index], state.data[1:state.max_index])
reset!(state)
