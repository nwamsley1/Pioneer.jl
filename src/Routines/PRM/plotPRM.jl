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
using LaTeXStrings
function plotPairedFragmentIonChromatogram(light_transitions::UnorderedDictionary{String, Vector{Float32}}, 
    heavy_transitions::UnorderedDictionary{String, Vector{Float32}}, light_rts::Vector{Float32}, heavy_rts::Vector{Float32}, title::String, 
    par::Real, dev_ratio::Real, out_path::String)
    p = plot(title = title, fontfamily="helvetica", legend=:outertopright)

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
            plot!(p, light_rts, -1*(max_heavy/max_light)*light_transitions[t], color = color, legend = false, label = nothing)
            plot!(p, light_rts, -1*(max_heavy/max_light)*light_transitions[t], seriestype=:scatter, legend = false, color = color, label = nothing)
        end

        plot!(p, heavy_rts, heavy_transitions[t], color = color, legend = :outertopright, label = t)
        plot!(p, heavy_rts, heavy_transitions[t], seriestype=:scatter, color = color, label = nothing)
    end

    ticks = [1*max_heavy, 0, -1*(max_heavy/max_light)*max_light]
    ticklabels = [ @sprintf "%.2E" x for x in [max_heavy, 0, max_light]]
    par = @sprintf "%.2E" par
    dev_ratio = @sprintf "%.2E" dev_ratio
    annotate!(minimum(heavy_rts), max_heavy, fontfamily = "helvetica", text("PAR $par \nDev Ratio $dev_ratio", :black, :top, :left, 8))
    yticks!((ticks, 
             ticklabels), 
             axis = 1)
    #plot!(legend=:outertopright)
    savefig(joinpath(out_path, title*".pdf"))
end

ptable.lh_pair_id_to_light_heavy_pair[0x000005c]

function plotAllPairedFragmentIonChromatograms(ptable::ISPRMPrecursorTable, 
                                                chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, 
                                                ms_file_idx::UInt32, 
                                                out_path::String, 
                                                fname::String)
    if !isdir(out_path)
        mkpath(out_path)
    end
    for lh_pair in ptable.lh_pair_id_to_light_heavy_pair
        
        function getPrecursorChromatogram(prec_id::UInt32, chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram})
            if !isassigned(chromatograms, prec_id)
                return UnorderedDictionary{String, Vector{Float32}}(), Float32[]
            else
                return chromatograms[prec_id].transitions, chromatograms[prec_id].rts
            end
        end

        function getParModel(par_models::Dictionary{UInt32, ParModel}, ms_file_idx::UInt32)
            if !isassigned(par_models,  ms_file_idx)
                return 0.0, 0.0
            else
                par_models[ms_file_idx].par_model_coef[ms_file_idx], par_models[ms_file_idx].dev_ratio
            end
        end

        light_transitions, light_rts = getPrecursorChromatogram(lh_pair.light_prec_id, chromatograms)
        heavy_transitions, heavy_rts = getPrecursorChromatogram(lh_pair.heavy_prec_id, chromatograms)
        par, dev_ratio = getParModel(lh_pair.par_model, ms_file_idx)

        if (length(light_transitions) == 0) & (length(heavy_transitions) == 0)
            continue
        end
        plotPairedFragmentIonChromatogram(
            light_transitions, 
            heavy_transitions,
            light_rts, 
            heavy_rts,
            lh_pair.light_sequence,
            par, 
            dev_ratio,
            out_path)
    end
        files = filter(x -> isfile(joinpath(out_path, x)) && match(r"\.pdf$", x) != nothing, readdir(out_path))
        merge_pdfs(map(file -> joinpath(out_path,file), files), joinpath(out_path, fname), cleanup=true)

end

