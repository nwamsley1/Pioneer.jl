using Arrow, CairoMakie

function analyze_ms2_dynamic_range(filepath::String; label::String="", plot_suffix::String="")
    tbl = Arrow.Table(filepath)
    ms_order = tbl[:msOrder]
    intensity_arrays = tbl[:intensity_array]

    log10_ratios = Float64[]
    total_peaks = 0
    below_max_1e4 = 0
    below_max_1e3 = 0
    below_max_1e2 = 0
    below_sum_1e4 = 0
    below_sum_1e3 = 0
    below_sum_1e2 = 0
    n_ms2 = 0
    n_empty = 0

    for i in 1:length(ms_order)
        ms_order[i] == 2 || continue
        intensities = intensity_arrays[i]
        n = length(intensities)
        if n == 0
            n_empty += 1
            continue
        end
        n_ms2 += 1

        mx = Float64(-Inf)
        mn = Float64(Inf)
        sm = Float64(0)
        for j in 1:n
            v = Float64(intensities[j])
            v > mx && (mx = v)
            v < mn && (mn = v)
            sm += v
        end

        if mx > 0 && mn > 0
            push!(log10_ratios, log10(mx / mn))
        end

        total_peaks += n

        # Thresholds relative to scan max
        t_max_4 = mx * 1e-4
        t_max_3 = mx * 1e-3
        t_max_2 = mx * 1e-2
        # Thresholds relative to summed intensity (TIC)
        t_sum_4 = sm * 1e-4
        t_sum_3 = sm * 1e-3
        t_sum_2 = sm * 1e-2

        for j in 1:n
            v = Float64(intensities[j])
            v < t_max_4 && (below_max_1e4 += 1)
            v < t_max_3 && (below_max_1e3 += 1)
            v < t_max_2 && (below_max_1e2 += 1)
            v < t_sum_4 && (below_sum_1e4 += 1)
            v < t_sum_3 && (below_sum_1e3 += 1)
            v < t_sum_2 && (below_sum_1e2 += 1)
        end
    end

    println("=" ^ 70)
    println(label == "" ? filepath : label)
    println("=" ^ 70)
    println("Total MS2 scans: $n_ms2 (empty: $n_empty)")
    println("Total MS2 peaks: $total_peaks")
    println()

    println("--- Relative to SCAN MAX ---")
    println("  < 1e-4 of max: $below_max_1e4 / $total_peaks = $(round(100*below_max_1e4/total_peaks, digits=2))%")
    println("  < 1e-3 of max: $below_max_1e3 / $total_peaks = $(round(100*below_max_1e3/total_peaks, digits=2))%")
    println("  < 1e-2 of max: $below_max_1e2 / $total_peaks = $(round(100*below_max_1e2/total_peaks, digits=2))%")
    println()

    println("--- Relative to SUMMED INTENSITY (TIC) ---")
    println("  < 1e-4 of sum: $below_sum_1e4 / $total_peaks = $(round(100*below_sum_1e4/total_peaks, digits=2))%")
    println("  < 1e-3 of sum: $below_sum_1e3 / $total_peaks = $(round(100*below_sum_1e3/total_peaks, digits=2))%")
    println("  < 1e-2 of sum: $below_sum_1e2 / $total_peaks = $(round(100*below_sum_1e2/total_peaks, digits=2))%")
    println()

    sort!(log10_ratios)
    println("Log10(max/min) ratio quantiles:")
    for q in [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]
        idx = max(1, min(length(log10_ratios), round(Int, q * (length(log10_ratios)-1) + 1)))
        println("  $(Int(q*100))%: $(round(log10_ratios[idx], digits=2))")
    end

    # Plot histogram
    plotname = plot_suffix == "" ? "p_ms2_dynamic_range.png" : "p_ms2_dynamic_range_$(plot_suffix).png"
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1,1],
        xlabel = "log₁₀(max / min intensity) per MS2 scan",
        ylabel = "Count",
        title = label == "" ? "Dynamic Range (n=$(length(log10_ratios)))" : "$label (n=$(length(log10_ratios)))")
    hist!(ax, log10_ratios, bins=100, color=:steelblue)
    save(plotname, fig, px_per_unit=2)
    println("\nHistogram saved to $plotname")
    println()
end

# Astral file
analyze_ms2_dynamic_range(
    "/Users/nathanwamsley/Data/For_Figures/OlsenAstralThreeProteome200ng/20230324_OLEP08_200ng_30min_E5H50Y45_180K_2Th3p5ms_01.arrow";
    label="Astral ThreeProteome 200ng",
    plot_suffix="astral"
)

# Eclipse file
analyze_ms2_dynamic_range(
    "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow";
    label="Eclipse MixedSpecies 500ng",
    plot_suffix="eclipse"
)
