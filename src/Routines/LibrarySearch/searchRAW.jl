function SearchRAW(
                    spectra::Arrow.Table, 
                    #ptable::PrecursorDatabase,
                    frag_index::FragmentIndex{Float32},
                    fragment_list::Vector{Vector{LibraryFragment{Float32}}},
                    ms_file_idx::UInt32,
                    interpolation::Interpolations.Extrapolation{Float64, 1, Interpolations.GriddedInterpolation{Float64, 1, Float64, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Line{Nothing}};
                    isolation_width::Float64 = 4.25,
                    precursor_tolerance::Float64 = 5.0,
                    fragment_tolerance::Float64 = 20.0,
                    topN::Int64 = 20,
                    min_frag_count::Int64 = 4,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    λ::Float32 = Float32(1e3),
                    γ::Float32 = zero(Float32),
                    scan_range::Tuple{Int64, Int64} = (0, 0), 
                    max_peaks::Int = 200, 
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0),
                    rt_tol::Float32 = Float32(30.0),
                    #fragment_match_ppm::U,
                    data_type::Type{T} = Float64
                    ) where {T<:Real}
    
    scored_PSMs = makePSMsDict(XTandem(data_type))
    ms2, MS1, MS1_i = 0, 0, 0
    precs = Counter(UInt32, UInt8, Float32, length(fragment_list)) #Prec counter
    times = Dict(:counter => 0.0, :reset => 0.0, :nmf => 0.0, :metrics => 0.0, :match_peaks => 0.0, :build => 0.0, :score => 0.0)
    n = 0
    kt = KendallTau(Float32)
    for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    #for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))

        if spectrum[:msOrder] == 1
            MS1 = spectrum[:masses]
            MS1_i = spectrum[:intensities]
            continue
        else
            ms2 += 1
        end

        if scan_range != (0, 0)
            i < first(scan_range) ? continue : nothing
            i > last(scan_range) ? continue : nothing
        end


        min_intensity = spectrum[:intensities][sortperm(spectrum[:intensities], rev = true)[min(max_peaks, length(spectrum[:intensities]))]]

        iRT = interpolation(spectrum[:retentionTime])
        iRT_low = Float32(iRT - rt_tol)
        iRT_high = Float32(iRT + rt_tol)
        #reset!(precs)

        times[:counter] += @elapsed  prec_count, match_count = searchScan!(precs,
                    frag_index, 
                    min_intensity, spectrum[:masses], spectrum[:intensities], MS1, spectrum[:precursorMZ], iRT_low, iRT_high,
                    Float32(fragment_tolerance), 
                    Float32(precursor_tolerance),
                    Float32(isolation_width),
                    min_frag_count = min_frag_count, 
                    min_ratio = Float32(min_matched_ratio),
                    topN = topN
                    )
        
        #println("prec_count $prec_count")
        #println("match_count $match_count")
        #return precs
        #return precs
        #if getSize(precs) <= 1
            #times[:reset] += @elapsed reset!(precs)
        #    reset!(precs)
        #    continue
        #end

        transitions = selectTransitions(fragment_list, precs, topN)

        #println("length(transitions) ", length(transitions))
        #return transitions
        #times[:reset] += @elapsed reset!(precs)
        reset!(precs)
        
        if length(transitions) == 0
            continue
        end
        
       

        #times[:match_peaks] += @elapsed fragmentMatches, fragmentMisses = matchPeaks(transitions, 
        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    FragmentMatch{Float32},
                                    count_unmatched=true,
                                    δs = zeros(T, (1,)),
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity,
                                    ppm = fragment_tolerance
                                    )

        if iszero(length(fragmentMatches))
            continue
        end

        #times[:build] += @elapsed X, H, UNMATCHED, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses, topN)
        X, Hs, Hst, IDtoROW, matched_cols = buildDesignMatrix(fragmentMatches, fragmentMisses)
        
        times[:nmf] += @elapsed weights = sparseNMF(Hst, Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)[:]

        scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all, kt_pval = getDistanceMetrics(Hst, X, matched_cols, kt)
        
        #For progress and debugging. 

        unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)

        #times[:score] += @elapsed Score!(scored_PSMs, unscored_PSMs, 
        Score!(scored_PSMs, unscored_PSMs, 
                length(spectrum[:intensities]), 
                Float64(sum(spectrum[:intensities])), 
                match_count/prec_count, scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all, kt_pval, weights, IDtoROW,
                scan_idx = Int64(i),
                min_spectral_contrast = min_spectral_contrast,
                min_frag_count = min_frag_count
                )
        n += 1
    end

    println("processed $ms2 scans!")
    println("counter: ", times[:counter]/n)
    println("reset: ", times[:reset]/n)
    println("nmf : ", times[:nmf]/n)
    println("metrics : ", times[:metrics]/n)
    println("build : ", times[:build]/n)
    println("score : ", times[:score]/n)
    println("match_peaks : ", times[:match_peaks]/n)
    return DataFrame(scored_PSMs)
end

#=
CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/psms_071923.csv",PSMs)

test_psms = SearchRAW(MS_TABLE, prosit_index_5ppm_15irt, frag_detailed, UInt32(1), linear_spline,
       min_frag_count = 4, 
       topN = 200, 
       fragment_tolerance = 15.6, 
       λ = Float32(1e4), 
       γ =Float32(1),
       max_peaks = 1000, 
       scan_range = (0, 300000), #101357 
       precursor_tolerance = 20.0,
       min_spectral_contrast =  Float32(0.65),
       rt_tol = Float32(15.0)
       )
2×18 DataFrame
 Row │ b_ladder  city_block  entropy     error       hyperscore  intensity_explained  matched_ratio  ms_file_idx  poisson   precursor_idx  scan_idx  scribe_score  spectral_contrast_all  spectral_contrast_matched  spectrum_peaks  total_ions  ⋯
     │ Int8      Float32     Float64     Float64     Float64     Float64              Float32        UInt32       Float64   UInt32         Int64     Float32       Float32                Float32                    UInt32          UInt32      ⋯
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │        6    -3.64194  -5.08692e8  0.00273705     94.9202            0.233232        5.26499             1  -83.0898        4211123    101357      15.9464                0.939872                   0.979832             642          32  ⋯
   2 │        1    -2.55285  -6.99821e7  0.00474243     38.5233            0.0327072       0.751578            1  -16.2589        6522764    101357       9.14359               0.817562                   0.918068             642          10
=#
   #=

    prosit_index_5ppm_5irt.fragment_bins[50893]
FragBin{Float32}(397.25577f0, 397.2571f0, 0x0000c6cd)
PSMs = DataFrame(CSV.File("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv"))

 LibraryFragment{Float32}(397.25577f0, 0x01, true, 0x03, 0x01, 0.21178919f0, 0x03, 0x004041b3)
 LibraryFragment{Float32}(133.09012f0, 0x03, true, 0x03, 0x02, 0.0024077685f0, 0x03, 0x004041b3)
 LibraryFragment{Float32}(212.10297f0, 0x01, false, 0x03, 0x03, 0.19929862f0, 0x03, 0x004041b3)
 LibraryFragment{Float32}(468.2929f0, 0x01, true, 0x04, 0x04, 0.23923871f0, 0x03, 0x004041b3)

best_psms = combine(sdf -> sdf[argmin(sdf.RT_error),:], groupby(PSMs[(PSMs[:,:q_values].<=0.005).&(PSMs[:,:decoy].==false).&(PSMs[:,:RT_error].<1000.0),:], :precursor_idx))
best_psms = combine(sdf -> sdf[argmin(sdf.RT_error),:], groupby(PSMs[(PSMs[:,:matched_ratio].>10.0).&(PSMs[:,:decoy].==false).&(PSMs[:,:spectral_contrast_all].>0.97).&(PSMs[:,:q_values].<0.01),:], :precursor_idx))
plot(best_psms[:,:iRT], best_psms[:,:RT], seriestype=:scatter, xlabel = "iRT", ylabel = "RT")

ns1 = Splines2.ns_(collect(range(0.0, length=100, stop=130.0)),df=4,intercept=true);
X = ns1(best_psms[:,:RT]);
fit1 = lm(X,best_psms[:,:iRT]);
best_psms[:,:iRT_pred] =  GLM.predict(fit1, ns1(best_psms[:,:RT]))
plot(best_psms[:,:iRT], best_psms[:,:RT], seriestype=:scatter, xlabel = "iRT", ylabel = "RT")

#GLM.predict(fit1, ns1(collect(range(10.0, length = 100, stop = 130.0))))
plot!(GLM.predict(fit1, ns1(collect(range(0.0, length = 100, stop = 130.0)))),collect(range(10.0, length = 100, stop = 130.0)), xlabel = "iRT", ylabel = "RT")

sort_order = sortperm(best_psms[:,:iRT])
histogram(best_psms[sort_order[1000:1200],:RT])
#X, H, UNMATCHED, IDtoROW, fragmentMatches, fragmentMisses = 
test_psms = SearchRAW(MS_TABLE, prosit_index_5ppm_5irt, frag_detailed, UInt32(1), linear_spline,
       min_frag_count = 4, 
       topN = 200, 
       fragment_tolerance = 15.6, 
       λ = Float32(1e4), 
       γ =Float32(1),
       max_peaks = 1000, 
       scan_range = (101357, 101357), 
       precursor_tolerance = 20.0,
       min_spectral_contrast =  Float32(0.65)
       )


       using MultiKDE
using Distributions, Random, Plots
using MultiKDE
sort_order = sortperm(best_psms[:,:RT])
bins = [Int(x÷1) for x in range(1, length = 200, stop = length(sort_order))]
ys = []
xs = Float64[]
rts = Float64[]
for i in 1:(length(bins) - 1)
    observations = best_psms[sort_order[bins[i]:(bins[i + 1])],:iRT]
    x = Vector(LinRange(minimum(observations), maximum(observations), 50))
    kde = KDEUniv(ContinuousDim(), 3.0, observations, MultiKDE.gaussian)
    y = [MultiKDE.pdf(kde, _x, keep_all=false) for _x in x]
    push!(ys, y)
    push!(xs, x[argmax(y)])
    push!(rts, mean(best_psms[sort_order[bins[i]:(bins[i + 1])],:RT]))
    #push!(rts, best_psms[sort_order[bins[i]],:RT])
end
plot(best_psms[:,:RT], best_psms[:,:iRT], seriestype=:scatter, xlabel = "RT", ylabel = "iRT")
plot!(rts, xs)
linear_spline = LinearInterpolation(rts, xs, extrapolation_bc = Line() )
ns1 = Splines2.ns_(collect(range(0.0, length=30, stop=130.0)),df=10,intercept=false);
X = ns1(rts);
fit1 = lm(X,xs);
best_psms[:,:iRT_pred] =  linear_spline(best_psms[:,:RT])
plot!(best_psms[:,:RT], best_psms[:,:iRT_pred], xlabel = "RT", ylabel = "iRT")



highest = maximum([maximum(y) for y in ys])
plot(xs[1], ys[1], label=bws, fmt=:svg)
plot!(observations, [highest+0.05 for _ in 1:length(ys)], seriestype=:scatter, label="observations", size=(900, 450), legend=:outertopright)

# Simulation
# Simulation
bws = [1 3.0 5.0]
d = Normal(0, 1)
observations = best_psms[sort_order[1:200],:RT]
granularity_1d = 100
x = Vector(LinRange(0.0, 130.0, granularity_1d))
ys = []
for bw in bws
    kde = KDEUniv(ContinuousDim(), bw, observations, MultiKDE.gaussian)
    y = [MultiKDE.pdf(kde, _x, keep_all=false) for _x in x]
    push!(ys, y)
end

# Plot
highest = maximum([maximum(y) for y in ys])
plot(x, ys, label=bws, fmt=:svg)
plot!(observations, [highest+0.05 for _ in 1:length(ys)], seriestype=:scatter, label="observations", size=(900, 450), legend=:outertopright)
=#