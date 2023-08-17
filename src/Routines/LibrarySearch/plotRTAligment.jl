p = Plots.plot(rtPSMs[:,:RT], rtPSMs[:,:iRT], seriestype=:scatter,
            xlabel = "Retention Time RT (min)",
            ylabel ="Indexed Retention Time iRT (min)",
            label = nothing,
            size = 100*[13.3, 7.5]
,
fontsize = 24,
titlefontsize = 24,
legendfontsize = 24,
tickfontsize = 24,
guidefontsize = 24,
margin = 10mm)

Plots.vline!(p, [params(mix_mle)[1][1][1]], lw = 6.0, color = :black, label = "Estimated Mass Error (ppm)")
rt_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT], n = 50)
Plots.plot!(p, LinRange(minimum(rtPSMs[:,:RT]), maximum(rtPSMs[:,:RT]), 100), 
            rt_map.(LinRange(minimum(rtPSMs[:,:RT]), maximum(rtPSMs[:,:RT]), 100)),
            lw = 6.0,
            label = "RT Spline")
Plots.savefig(p, "/Users/n.t.wamsley/Documents/presentations/august_labmeeting/rt_spline.pdf")

function plotRTAlign()
end