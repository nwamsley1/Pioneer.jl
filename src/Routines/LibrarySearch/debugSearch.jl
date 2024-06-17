function plotPSM(
                    spectra::Arrow.Table,
                    scan_idx::Int64,#UnitRange{Int64},
                    precursor_idx::UInt32,
                    precursors::Union{Arrow.Table, Missing},
                    library_fragment_lookup::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt::UniformSpline,
                    mass_err_model::MassErrorModel,
                    δ::Float32,
                    λ::Float32,
                    max_iter_newton::Int64,
                    max_iter_bisection::Int64,
                    max_iter_outer::Int64,
                    accuracy_newton::Float32,
                    accuracy_bisection::Float32,
                    max_diff::Float32,
                    ionMatches::Vector{FragmentMatch{Float32}},
                    ionMisses::Vector{FragmentMatch{Float32}},
                    IDtoCOL::ArrayDict{UInt32, UInt16},
                    ionTemplates::Vector{L},
                    iso_splines::IsotopeSplineModel,
                    isotope_err_bounds::Tuple{Int64, Int64},
                    n_frag_isotopes::Int64,
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing},
                    irt_tol::Float64,
                    ) where {T,U<:AbstractFloat, L<:LibraryIon{Float32}}

    ##########
    #Initialize 
    prec_idx, ion_idx = 0, 0
    Hs = SparseArray(UInt32(5000));
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);
    isotopes = zeros(Float32, 5)
    irt_start, irt_stop = 1, 1
    prec_temp_size = 0
    precs_temp = Vector{UInt32}(undef, 50000)

    irt = rt_to_irt(spectra[:retentionTime][scan_idx])
    irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
    irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 
    #Candidate precursors and their retention time estimates have already been determined from
    #A previous serach and are incoded in the `rt_index`. Add candidate precursors that fall within
    #the retention time and m/z tolerance constraints
    precs_temp_size = 0
    ion_idx, prec_idx, prec_temp_size = selectRTIndexedTransitions!(
        ionTemplates,
        precs_temp,
        precs_temp_size,
        library_fragment_lookup,
        precursors[:mz],
        precursors[:prec_charge],
        precursors[:sulfur_count],
        iso_splines,
        isotopes,
        n_frag_isotopes,
        rt_index,
        irt_start,
        irt_stop,
        spectra[:centerMass][scan_idx] - spectra[:isolationWidth][scan_idx]/2.0f0,
        spectra[:centerMass][scan_idx] + spectra[:isolationWidth][scan_idx]/2.0f0,
        (
            spectra[:lowMass][scan_idx], spectra[:highMass][scan_idx]
        ),
        isotope_err_bounds,
        10000)

    ##########
    #Match sorted list of plausible ions to the observed spectra
    nmatches, nmisses = matchPeaks!(ionMatches, 
                                    ionMisses, 
                                    ionTemplates, 
                                    ion_idx, 
                                    spectra[:masses][scan_idx], 
                                    spectra[:intensities][scan_idx], 
                                    mass_err_model,
                                    spectra[:highMass][scan_idx],
                                    UInt32(scan_idx), 
                                    zero(UInt32))
    sort!(@view(ionMatches[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
    ##########
    #Spectral Deconvolution
        #Spectral deconvolution. Build sparse design/template matrix for regression 
        #Sparse matrix representation of templates written to Hs. 
        #IDtoCOL maps precursor ids to their corresponding columns. 
        buildDesignMatrix!(Hs, ionMatches, ionMisses, nmatches, nmisses, IDtoCOL)
        #Adjuste size of pre-allocated arrays if needed 
        if IDtoCOL.size > length(_weights_)
            new_entries = IDtoCOL.size - length(_weights_) + 1000 
            append!(_weights_, zeros(eltype(_weights_), new_entries))
        end
        #Get most recently determined weights for each precursors
        #Get initial residuals
        initResiduals!(_residuals_, Hs, _weights_);
        #Spectral deconvolution. Hybrid bisection/newtowns method
        solveHuber!(Hs, _residuals_, _weights_, 
                        δ, λ, 
                        max_iter_newton, 
                        max_iter_bisection,
                        max_iter_outer,
                        accuracy_newton,
                        accuracy_bisection,
                        10.0,#Hs.n/10.0,
                        max_diff
                        );

        p = plot(title = "scan id $scan_idx")
        reconstructed_spectra = zeros(Float32, length(spectra[:masses][scan_idx]))
        for i in range(1, nmatches)
            peak_id = getPeakInd(ionMatches[i])
            prec_id = getPrecID(ionMatches[i])
            intensity = ionMatches[i].predicted_intensity
            w = _weights_[IDtoCOL[prec_id]]
            reconstructed_spectra[peak_id] += w*intensity
        end
        w = _weights_[IDtoCOL[precursor_idx]]
        for i in range(1, nmatches)
            if getPrecID(ionMatches[i])==precursor_idx
                match = ionMatches[i]
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.intensity], color = 1, alpha = 1.0, lw = 3, label = nothing)
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, reconstructed_spectra[getPeakInd(match)]], color = 2, alpha = 0.5, lw = 3, label = nothing)
                    plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 3, alpha = 0.5, lw = 3, label = nothing)
            end
        end
        for i in range(1, nmisses)
            if getPrecID(ionMisses[i])==precursor_idx
                match = ionMisses[i]
                plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 3, alpha = 0.5, lw = 3, label = nothing)
                annotate!(p, match.theoretical_mz, match.predicted_intensity*w, "*")
            end
        end
        plot!(p, [ionMatches[1].theoretical_mz], [0], color = 1, label = "Experimental Spectrum")
        plot!(p, [ionMatches[1].theoretical_mz], [0], color = 2, label = "Reconstructed Spectrum")
        plot!(p, [ionMatches[1].theoretical_mz], [0], color = 3, label = "Target Spectrum")
        savefig(p, "/Users/n.t.wamsley/Desktop/prec_id_"*string(precursor_idx)*"_scan_id_"*string(scan_idx)*".pdf")
        reset!(IDtoCOL);
        reset!(Hs);
    return
end

function plotPSM(
    spectra::Arrow.Table,
    scan_idx::Int64,
    precursor_idx::UInt32,
    params::Any;
    kwargs...)

    plotPSM(
        spectra,
        scan_idx,
        precursor_idx,
        kwargs[:precursors],
        kwargs[:fragment_lookup_table],
        kwargs[:rt_to_irt_spline],
        kwargs[:mass_err_model],
        Float32(params[:deconvolution_params]["huber_delta"]),
        Float32(params[:deconvolution_params]["lambda"]),
        Int64(params[:deconvolution_params]["max_iter_newton"]),
        Int64(params[:deconvolution_params]["max_iter_bisection"]),
        Int64(params[:deconvolution_params]["max_iter_outer"]),
        Float32(params[:deconvolution_params]["accuracy_newton"]),
        Float32(params[:deconvolution_params]["accuracy_bisection"]),
        Float32(params[:deconvolution_params]["max_diff"]),
        kwargs[:ion_matches][1],
        kwargs[:ion_misses][1],
        kwargs[:id_to_col][1],
        kwargs[:ion_templates][1],
        kwargs[:iso_splines],
        (3, 0),
        params[:quant_search_params]["n_frag_isotopes"],
        kwargs[:rt_index], 
        kwargs[:irt_err]
    )
end

@time RESULT = plotPSM(
    MS_TABLE, 
    216885,
    UInt32(2664945),
    params_;
    precursors = prosit_lib["precursors"],
    fragment_lookup_table = library_fragment_lookup_table,
    rt_index = RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
    ms_file_idx = UInt32(ms_file_idx), 
    rt_to_irt_spline = RT_iRT[file_id_to_parsed_name[ms_file_idx]],
    mass_err_model =  MassErrorModel{Float32}(3.195095f0, (15.593506f0*1.5f0, 7.698679f0*1.5f0)),
    irt_err = irt_err,#irt_errs[ms_file_idx]/3,
    ion_matches = ionMatches,
    ion_misses = ionMisses,
    id_to_col = IDtoCOL,
    ion_templates = ionTemplates,
    iso_splines = iso_splines,
    chromatograms = chromatograms,
    scored_psms = complex_scored_PSMs,
    unscored_psms = complex_unscored_PSMs,
    spectral_scores = complex_spectral_scores,
    precursor_weights = precursor_weights,
    );


function plotChrom!(
                    p::Plots.Plot{Plots.GRBackend},
                    subplot_idx::Int64,
                    spectra::Arrow.Table,
                    scan_idxs::AbstractArray{UInt32},#UnitRange{Int64},
                    weight::AbstractArray{Float32},
                    fragment_ions::Vector{DetailedFrag{Float32}},
                    isotopes_captured::Tuple{Int8, Int8},
                    iso_splines::IsotopeSplineModel{40, Float32},
                    prec_mz::Float32,
                    prec_charge::UInt8,
                    prec_sulfur_count::UInt8,
                    ionMatches::Vector{FragmentMatch{Float32}},
                    ionMisses::Vector{FragmentMatch{Float32}};
                    precursor_name = "TEST",
                    max_frag_rank = 6)

    ##########
    #Initialize 
    isotopes = zeros(Float32, 5)
    fragment_ions = copy(fragment_ions)
    sort!(fragment_ions, by = x->x.intensity, rev = true)
    ###########
    #Get Isotopes and Isotope Correction
    fragment_ions_isotoped = Vector{DetailedFrag{Float32}}()
    for frag in fragment_ions
        getFragIsotopes!(
            isotopes,
            iso_splines,
            prec_mz,
            prec_charge,
            prec_sulfur_count,
            frag,
            Tuple([Int64(x) for x in isotopes_captured])
        )
        for iso_idx in range(0, 1)
            frag_mz = Float32(frag.mz + iso_idx*NEUTRON/frag.frag_charge)
            push!(fragment_ions_isotoped,
            DetailedFrag(
                frag.prec_id,

                frag_mz, #Estimated isotopic m/z
                Float16(isotopes[iso_idx + 1]), #Estimated relative abundance 

                frag.ion_type,
                iso_idx>0, #Is the fragment an isotope?

                frag.frag_charge,
                frag.ion_position,
                frag.prec_charge,
                frag.rank,
                frag.sulfur_count
            )
            )
        end
    end
    fragment_ions = fragment_ions_isotoped[1:min(length(fragment_ions_isotoped), max_frag_rank*2)]

    function getIsotopeFragments(fragment_ions::Vector{DetailedFrag{Float32}},
                                 is_isotope::Bool)
        fragment_ions = [frag for frag in fragment_ions if frag.is_isotope == is_isotope]
        ion_types = [x.ion_type for x in fragment_ions]
        peak_ind = [x.ion_position for x in fragment_ions]
        frag_charge = [x.frag_charge for x in fragment_ions]
        ion_names = String[]
        types = ('b','y','p')
        for i in range(1, length(fragment_ions))
            if fragment_ions[i].is_isotope == is_isotope
                push!(ion_names, types[ion_types[i]]*string(peak_ind[i])*"^"*string(UInt8(is_isotope))*"+"*string(frag_charge[i]))
            end
        end
        #Need to sort fragments by mz. 
        sort!(fragment_ions, by = x->x.mz)
        return ion_names, fragment_ions
    end

    function buildChromatogram(
        spectra::Arrow.Table,
        scan_idxs::AbstractArray{UInt32},
        fragment_ions::Vector{DetailedFrag{Float32}},
        ionMatches::Vector{FragmentMatch{Float32}},
        ionMisses::Vector{FragmentMatch{Float32}})
        rts = zeros(Float32, length(scan_idxs))
        empirical_intensities = zeros(Float32, (length(scan_idxs), length(fragment_ions)))
        fitted_intensities = zeros(Float32, (length(scan_idxs), length(fragment_ions)))
        for (i, scan_idx) in enumerate(scan_idxs)
            ##########
            #Match sorted list of plausible ions to the observed spectra
            nmatches, nmisses = matchPeaks!(ionMatches, 
                                            ionMisses, 
                                            fragment_ions, 
                                            length(fragment_ions), 
                                            spectra[:masses][scan_idx], 
                                            spectra[:intensities][scan_idx], 
                                            mass_err_model,
                                            spectra[:highMass][scan_idx],
                                            UInt32(scan_idx), 
                                            zero(UInt32))
            rts[i] = spectra[:retentionTime][scan_idx]
            for j in range(1, nmatches)
                m = ionMatches[j]
                empirical_intensities[i,m.predicted_rank] = m.intensity
                fitted_intensities[i,m.predicted_rank] = m.predicted_intensity*weight[i]
            end
            for j in range(1, nmisses)
                m = ionMisses[j]
                empirical_intensities[i,m.predicted_rank] = 0.0
                fitted_intensities[i,m.predicted_rank] = m.predicted_intensity*weight[i]
            end
        end
        return empirical_intensities, fitted_intensities, rts
    end

    function plotChromatogram!(
        p::Plots.Plot{Plots.GRBackend},
        ion_names::Vector{String},
        empirical_intensities::Matrix{<:AbstractFloat},
        fitted_intensities::Matrix{<:AbstractFloat},
        weight::AbstractVector{<:AbstractFloat},
        rts::AbstractVector{<:AbstractFloat},
        subplot_idx::Int64)

        norm = maximum(fitted_intensities)
        empirical_intensities ./= norm;
        fitted_intensities ./= norm;
        peak_apex = argmax(weight)
        start = max(peak_apex-5, 1)
        stop = min(peak_apex + 5, length(weight))
        for i in range(1, size(empirical_intensities, 2))
            plot!(p, rts[start:stop], empirical_intensities[start:stop,i], show = true, color = i, seriestype=:scatter, label = nothing, subplot = subplot_idx)
            plot!(p, rts[start:stop], empirical_intensities[start:stop,i], show = true, color = i, lw = 2, label = ion_names[i], subplot = subplot_idx)
        end
        for i in range(1, size(empirical_intensities, 2))
            plot!(p, rts[start:stop], (-1)*fitted_intensities[start:stop,i], show = true, color = i, seriestype=:scatter, label = nothing, subplot = subplot_idx)
            plot!(p, rts[start:stop],(-1)*fitted_intensities[start:stop,i], show = true, color = i, lw = 2, label = nothing, subplot = subplot_idx)
        end
        #println(weight.*((-1)/max(weight)))
        plot!(p, rts[start:stop], (weight.*((-1)/maximum(weight)))[start:stop], show = true, seriestype=:scatter, label = "weight", color = :grey,subplot = subplot_idx)
        plot!(p, rts[start:stop], (weight.*((-1)/maximum(weight)))[start:stop], show = true, label = nothing, color = :grey,subplot = subplot_idx)
        hline!(p, [1.0, -1.0], color = :black, lw = 2, label = nothing, subplot = subplot_idx)
    end

    ion_names_m0, fragment_ions_m0 = getIsotopeFragments(fragment_ions,
                                                    false)

    empirical_intensities, fitted_intensities, rts = buildChromatogram(
        spectra,
        scan_idxs,
        fragment_ions_m0,
        ionMatches,
        ionMisses
    )

    plotChromatogram!(
        p, 
        ion_names_m0,
        empirical_intensities,
        fitted_intensities,
        weight,
        rts,
        subplot_idx,
    )

    ion_names_m1, fragment_ions_m1 = getIsotopeFragments(fragment_ions,
                                                    true)

    empirical_intensities, fitted_intensities, rts = buildChromatogram(
        spectra,
        scan_idxs,
        fragment_ions_m1,
        ionMatches,
        ionMisses
    )

    plotChromatogram!(
        p, 
        ion_names_m1,
        empirical_intensities,
        fitted_intensities,
        weight,
        rts,
        subplot_idx + 1
    )
  
    return
end

subdf = uncorrected_gchroms[(precursor_idx = precs_passing[N],)]
sort!(subdf, :rt)
gsubdf = groupby(subdf,:isotopes_captured)
p=plot(ylim = (-1.2, 1.2), layout = (size(gsubdf, 1), 2), title =  precursors[:sequence][precs_passing[N]])
plotChrom!(
    p, 
    1,
    MS_TABLE,
    gsubdf[1][!,:scan_idx],
    gsubdf[1][!,:intensity],
    library_fragment_lookup_table.frags[library_fragment_lookup_table.prec_frag_ranges[gsubdf[1][1,:precursor_idx]]],
    gsubdf[1][1,:isotopes_captured],
    iso_splines,
    precursors[:mz][precs_passing[N]],
    precursors[:prec_charge][precs_passing[N]],
    precursors[:sulfur_count][precs_passing[N]],
    ionMatches[1],
    ionMisses[1],
    precursor_name = string(gsubdf[1][1,:isotopes_captured])
);
plotChrom!(
    p, 
    3,
    MS_TABLE,
    gsubdf[2][!,:scan_idx],
    gsubdf[2][!,:intensity],
    library_fragment_lookup_table.frags[library_fragment_lookup_table.prec_frag_ranges[gsubdf[1][1,:precursor_idx]]],
    gsubdf[2][1,:isotopes_captured],
    iso_splines,
    precursors[:mz][precs_passing[N]],
    precursors[:prec_charge][precs_passing[N]],
    precursors[:sulfur_count][precs_passing[N]],
    ionMatches[1],
    ionMisses[1],
    precursor_name = precursors[:sequence][precs_passing[N]]
);
N += 1
#=
plotChrom(
    MS_TABLE,
    gsubdf[2][!,:scan_idx],
    gsubdf[2][!,:intensity],
    library_fragment_lookup_table.frags[library_fragment_lookup_table.prec_frag_ranges[gsubdf[1][1,:precursor_idx]]],
    gsubdf[2][1,:isotopes_captured],
    iso_splines,
    precursors[:mz][precs_passing[N]],
    precursors[:prec_charge][precs_passing[N]],
    precursors[:sulfur_count][precs_passing[N]],
    ionMatches[2],
    ionMisses[2],
    precursor_name = precursors[:sequence][precs_passing[N]]
);
=#
N += 1


    histogram(frag_ppm_errs, xlim = (-20, 20))
    frag_abs_errs = [(match.theoretical_mz-match.match_mz)/match.charge for match in matched_fragments];
    histogram(frag_abs_errs, xlim = (-0.02, 0.02))
    vline!([quantile(frag_abs_errs, 0.99)], color = :black, lw = 3)
    vline!([quantile(frag_abs_errs, 0.01)], color = :black, lw = 3)
    frag_abs_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in matched_fragments if match.frag_charge == 1];
    histogram(frag_abs_errs, xlim = (-20, 20), alpha = 0.5, normalize = :pdf)
    frag_abs_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in matched_fragments if match.frag_charge == 2];
    histogram!(frag_abs_errs, xlim = (-20, 20), alpha = 0.5, normalize = :pdf)
    vline!([quantile(frag_abs_errs, 0.99)], color = :black, lw = 3)
    vline!([quantile(frag_abs_errs, 0.01)], color = :black, lw = 3)


    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz)for match in matched_fragments];
    histogram2d(frag_ppm_errs, [match.theoretical_mz for match in matched_fragments], xlim = (-20, 20), alpha = 1.0, normalize = :pdf)