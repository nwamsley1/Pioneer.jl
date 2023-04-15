function plotBestSpectra(matched_ions::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Float32, Int64, Vector{String}, Vector{Float32}, Vector{Float32}}}, RAW::Arrow.Table, title::String, out_path::String)

    function plotSpectra!(p::Plots.Plot{Plots.GRBackend}, masses, intensities)
        for (peak_idx, mass) in enumerate(masses)
            plot!(p, [mass, mass], [0, intensities[peak_idx]], legend = false, color = "black")
        end
    end
    #p = plot()
    #plotSpectra(p, NRF2_Survey.masses[2703], NRF2_Survey.intensities[2703])
    #display(p)
    function addFragmentIons!(p::Plots.Plot{Plots.GRBackend}, matched_ions::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Float32, Int64, Vector{String}, Vector{Float32}, Vector{Float32}}})
        i = 1
        for (mz, intensity) in zip(matched_ions[:mz], matched_ions[:intensity])
            plot!(p, [mz, mz], [0, -1*intensity], legend = false, color = "red")
            annotate!(mz + 1, -1*intensity, fontfamily = "helvetica", text(matched_ions[:name][i], :black, :top, :left, 8))
            i += 1
        end
    end

    p = plot(title = title, fontfamily="helvetica")
    plotSpectra!(p, RAW[:masses][matched_ions[:scan_idx]], RAW[:intensities][matched_ions[:scan_idx]])
    addFragmentIons!(p, matched_ions)
    savefig(joinpath(out_path,title*".pdf"))

end

function plotAllBestSpectra(matched_precursors::UnorderedDictionary{UInt32, PrecursorChromatogram}, ptable::PrecursorDatabase, MS_TABLE::Arrow.Table, out_path::String, fname::String)
    if !isdir(out_path)
        mkpath(out_path)
    end

    @time for key in keys(matched_precursors)
        peptide_sequence = ptable.id_to_pep[getPepIDFromPrecID(ptable, key)].sequence
        protein_name = join(getProtNamesFromPepSeq(ptable, peptide_sequence), "|")
        plotBestSpectra(matched_precursors[key].best_psm, MS_TABLE,
                        protein_name*"-"*peptide_sequence,
                        out_path)
    end
    files = filter(x -> isfile(joinpath(out_path, x)) && match(r"\.pdf$", x) != nothing, readdir(out_path))
    merge_pdfs(map(file -> joinpath(out_path,file), files), joinpath(out_path, fname), cleanup=true)

end

#plotAllBestSpectra(matched_precursors, testPtable, NRF2_Survey, "./data/figures/spectra/", "MERGED_NRF2_SURVEY_SPECTRA.pdf")

function plotFragmentIonChromatogram(transitions::UnorderedDictionary{String, Vector{Float32}}, rts::Vector{Float32}, title::String, out_path::String)
    p = plot(title = title, fontfamily="helvetica")
    for (color, t) in enumerate(keys(transitions))
        plot!(p, rts, transitions[t], color = color, legend = true, label = t)
        plot!(p, rts, transitions[t], seriestype=:scatter, color = color, label = nothing)
    end
    savefig(joinpath(out_path, title*".pdf"))
end

#testPtable.id_to_pep[1].sequence
function plotAllFragmentIonChromatograms(matched_precursors::UnorderedDictionary{UInt32, PrecursorChromatogram}, ptable::PrecursorDatabase, out_path::String, fname::String)

    if !isdir(out_path)
        mkpath(out_path)
    end

    @time for key in keys(matched_precursors)
        peptide_sequence = ptable.id_to_pep[getPepIDFromPrecID(ptable, key)].sequence
        protein_name = join(getProtNamesFromPepSeq(ptable, peptide_sequence), "|")

        plotFragmentIonChromatogram(matched_precursors[key].transitions, matched_precursors[key].rts, 
                                    protein_name*"-"*peptide_sequence, 
                                    out_path)
        # Get all files in the directory that match the pattern
    end
    files = filter(x -> isfile(joinpath(out_path, x)) && match(r"\.pdf$", x) != nothing, readdir(out_path))
    merge_pdfs(map(file -> joinpath(out_path,file), files), joinpath(out_path, fname), cleanup=true)
end

function plotPairedFragmentIonChromatogram(light_transitions::UnorderedDictionary{String, Vector{Float32}}, 
    heavy_transitions::UnorderedDictionary{String, Vector{Float32}}, light_rts::Vector{Float32}, heavy_rts::Vector{Float32}, title::String, out_path::String)
    p = plot(title = title, fontfamily="helvetica")
    twinx()
    max_light = 0.0
    max_heavy = 0.0

    for (color, t) in enumerate(keys(heavy_transitions))
        if isassigned(light_transitions, t)
            if maximum(light_transitions[t])>max_light
                max_light = maximum(light_transitions[t])
            end
        end
        if maximum(heavy_transitions[t])>max_heavy
            max_heavy = maximum(heavy_transitions[t])
        end
    end

    for (color, t) in enumerate(keys(heavy_transitions))

        if isassigned(light_transitions, t)
            plot!(p, light_rts, -1*(max_heavy/max_light)*light_transitions[t], color = color, legend = false, label = t, axis = 2)
            plot!(p, light_rts, -1*(max_heavy/max_light)*light_transitions[t], seriestype=:scatter, legend = false, color = color, label = nothing, axis = 2)
        end

        plot!(p, heavy_rts, heavy_transitions[t], color = color, legend = true, label = t, axis = 1)
        plot!(p, heavy_rts, heavy_transitions[t], seriestype=:scatter, color = color, label = nothing, axis = 1)
    end
    #ylims!((0, max_heavy), axis = 1)
    #ylims!((0, max_light), axis = 2)
    plot!(yticks = (max_heavy, string(max_heavy)), axis = 1)
    plot!(yticks = (-1*max_light, string(max_light)), axis = 2)
    yticks!(2, [5, 10, 15], ["Low", "Medium", "High"])
    savefig(joinpath(out_path, title*".pdf"))
end

ptable.lh_pair_id_to_light_heavy_pair[0x000005c]

plotPairedFragmentIonChromatogram(
precursor_chromatograms[1][0x000000b7].transitions,
precursor_chromatograms[1][0x000000b8].transitions,
precursor_chromatograms[1][0x000000b7].rts,
precursor_chromatograms[1][0x000000b8].rts,
"TEST",
"./")