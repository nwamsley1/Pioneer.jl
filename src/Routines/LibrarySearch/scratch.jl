pbar = ProgressBar(total=100)
       update(pbar)     # Advances pbar count by 1
       update(pbar, 2) # Can giv


function distributeScans(N::Int64, m::Int64)
scan_idx = 0
scans = Vector{NamedTuple{(:first_scan, :last_scan, :thread_id), Tuple{Int64, Int64, Int64}}}()
proc = 0
i = 1
first_scan = 1
last_scan = 1
while scan_idx <= N
if (i > m) | (scan_idx == N)
i = 0
proc += 1
push!(scans, (first_scan = first_scan, 
              last_scan = scan_idx,
              thread_id = proc%Threads.nthreads() + 1)
      )
scan_idx += 1
first_scan = scan_idx
continue
end
scan_idx += 1
i += 1
end
return scans
end



precs = [Counter(UInt32, Float32, length(ion_list)) for _ in range(1, Threads.nthreads())]
println(length(precs))
##########
#Initialize 

println("min_spectral_contrast ", min_spectral_contrast) #Remove precursors with spectral contrast lower than this ammount
println("min_matched_ratio ", min_matched_ratio) 
println("min_frag_count ", min_frag_count) 
println("min_weight ", min_weight) 
println("min_topn ", min_topn) 
###########
#Pre-allocate Arrays to save (lots) of time in garbage collection. 
all_fmatches = Vector{IonMatchType}()
collect_fmatches ? all_fmatches = [[IonMatchType() for x in range(1, expected_matches)] for _ in range(1, Threads.nthreads())] : nothing

#These are overwritten for every searched spectrum. "expected_matches"
#is a guess for the largest array size that would be needed for any spectrum. 
#If the guess is too small, the arrays will simply be increased in size as needed
#by a pre-determined block-size. 
println("START")
ionMatches = [[IonMatchType() for _ in range(1, expected_matches)]  for _ in range(1, Threads.nthreads())]#IonMatchType is something that inherits from the "Match" class. 
ionMisses = [[IonMatchType() for _ in range(1, expected_matches)] for _ in range(1, Threads.nthreads())]
ionTemplates = [[IonTemplateType() for _ in range(1, expected_matches)]  for _ in range(1, Threads.nthreads())]
prec_ids = [[zero(UInt32) for _ in range(1, expected_matches)] for _ in range(1, Threads.nthreads())]

IDtoCOL = nothing
if ismissing(precs)
IDtoCOL = [ArrayDict(UInt32, UInt16, length(precursors)) for _ in range(1, Threads.nthreads())]
else
IDtoCOL = [ArrayDict(UInt32, UInt16, length(first(precs).ids)) for _ in range(1, Threads.nthreads())]
end
#H_COLS, H_ROWS, H_VALS, H_MASK = zeros(Int64, expected_matches), zeros(Int64, expected_matches), zeros(Float32, expected_matches), zeros(Float32, expected_matches)
scored_PSMs = [Vector{LibPSM{Float32, Float16}}(undef, 5000) for _ in range(1, Threads.nthreads())];
unscored_PSMs = [[LXTandem(Float32) for _ in 1:5000] for _ in range(1, Threads.nthreads())];
spectral_scores = [Vector{SpectralScores{Float16}}(undef, 5000) for _ in range(1, Threads.nthreads())];
Hs = [SparseArray(50000) for _ in range(1, Threads.nthreads())];
println("STOP")

#weights
precursor_weights = ""
if ismissing(ion_list)
precursor_weights = [zeros(Float32, maximum(keys(isotope_dict))) for _ in range(1, Threads.nthreads())]
else
precursor_weights = [zeros(Float32, length(ion_list)) for _ in range(1, Threads.nthreads())]
end
_weights_ = [zeros(Float32, 5000) for _ in range(1, Threads.nthreads())];
_residuals_ = [zeros(Float32, 5000) for _ in range(1, Threads.nthreads())];
last_vals = [0 for _ in range(1, Threads.nthreads())];
#fragment_intensities = Dictionary{String, Vector{Tuple{Float32, Float32}}}()
##########
#Iterate through spectra
scan_batches = distributeScans(length(spectra[:masses]), 10000)
Threads.@threads :static for i in ProgressBar(1:length(scan_batches))
ThreadID = scan_batches[i][:thread_id] 
scan_range = (scan_batches[i][:first_scan], scan_batches[i][:last_scan])
msms_counts = Dict{Int64, Int64}()
frag_err_idx = 1
prec_idx = 0
ion_idx = 0
cycle_idx = 0 
minimum_rt, maximum_rt = first(rt_bounds), last(rt_bounds)
println("scan range ", scan_range, " assigned to thread ", ThreadID)
for i in range(scan_batches[i][:first_scan], scan_batches[i][:last_scan])
###########
#Scan Filtering
msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
if msn == 1
cycle_idx += 1
end
msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type

(i >= first(scan_range)) & (i <= last(scan_range)) ? nothing : continue #Skip if outside the scan range
first(rand(1)) <= sample_rate ? nothing : continue #dice-roll. Usefull for random sampling of scans. 

min_intensity = getMinIntensity(spectra[:intensities][i], max_peaks) #Ignore peaks in the spectrum below this minimum intensity

iRT_low, iRT_high = getRTWindow(iRT_to_RT_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, maximum_rt, minimum_rt, rt_tol) #Convert RT to expected iRT window

##########
#Ion Template Selection
#SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  
if !ismissing(searchScan!) | !ismissing(frag_index)
prec_count, match_count = searchScan!(
            precs[ThreadID], #counter which keeps track of plausible matches 
            frag_index, 
            min_intensity, spectra[:masses][i], spectra[:intensities][i], spectra[:precursorMZ][i], 
            iRT_low, iRT_high,
            Float32(fragment_tolerance), 
            Float32(precursor_tolerance),
            Float32(quadrupole_isolation_width/2.0),
            min_frag_count = min_frag_count_index_search, 
            min_ratio = Float32(min_index_search_score),
            topN = topN_index_search,#topN
            )

end
#selectIons! 
#Get a sorted list by m/z of ion templates (fills ionTemplates). The spectrum will be searched for matches to these ions only.
if !ismissing(precs) 
ion_idx, prec_idx = selectIons!(ionTemplates[ThreadID], 
                                                            precursors,
                                                            ion_list,
                                                            precs[ThreadID],
                                                            topN,
                                                            Float32(iRT_to_RT_spline(spectra[:retentionTime][i])),
                                                            Float32(rt_tol), #rt_tol
                                                            (
                                                            spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                                            spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                                            )
                                                            )::Tuple{Int64, Bool}
else
ion_idx, prec_idx = selectIons!(
                                ionTemplates[ThreadID],
                                precursors,
                                ion_list,
                                prec_ids[ThreadID],
                                rt_index,
                                spectra[:retentionTime][i],
                                Float32(max_peak_width/2),
                                spectra[:precursorMZ][i],
                                Float32(quadrupole_isolation_width/2.0),
                                (
                                    spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                    spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                    )
                                )
end
ion_idx < 2 ? continue : nothing 

##########
#Match sorted list of plausible ions to the observed spectra
nmatches, nmisses = matchPeaks(ionTemplates[ThreadID], #Search the spectra for these ions 
                                                ion_idx, #search ionTemplates[1:ion_idx]
                                                ionMatches[ThreadID], #Fill with matched ions 
                                                ionMisses[ThreadID], #Fill with unmatched ions 
                                                spectra[:masses][i], 
                                                spectra[:intensities][i], 
                                                count_unmatched=true, #Should we fill "ionMisses"?
                                                δs = [frag_ppm_err], #Mass offsets 
                                                scan_idx = UInt32(i),
                                                ms_file_idx = ms_file_idx,
                                                min_intensity = min_intensity, #Ignore peaks below this intensity
                                                ppm = fragment_tolerance, #Fragment match tolerance in ppm
                                                most_intense = most_intense
                                                )

##########
#Spectral Deconvolution and Distance Metrics 
if nmatches > 2 #Few matches to do not perform de-convolution 

#Spectral deconvolution. Build sparse design/template matrix for regression 
#Sparse matrix representation of templates written to Hs. 
#IDtoCOL maps precursor ids to their corresponding columns. 
buildDesignMatrix!(Hs[ThreadID], ionMatches[ThreadID], ionMisses[ThreadID], nmatches, nmisses, IDtoCOL[ThreadID])

if IDtoCOL[ThreadID].size > length(_weights_[ThreadID])
    append!(_weights_[ThreadID], zeros(eltype(_weights_[ThreadID]), IDtoCOL[ThreadID].size - length(_weights_[ThreadID]) + 1000 ))
    append!(spectral_scores[ThreadID], Vector{SpectralScores{Float16}}(undef, IDtoCOL[ThreadID].size - length(spectral_scores[ThreadID]) + 1000 ))
    append!(unscored_PSMs[ThreadID], [LXTandem(Float32) for _ in 1:(IDtoCOL[ThreadID].size - length(unscored_PSMs[ThreadID]) + 1000)]);
end
#Get most recently determined weights for eahc precursors
#"Hot" start
for i in range(1, IDtoCOL[ThreadID].size)#pairs(IDtoCOL)
    _weights_[ThreadID][IDtoCOL[ThreadID][IDtoCOL[ThreadID].keys[i]]] = precursor_weights[ThreadID][IDtoCOL[ThreadID].keys[i]]
end

#Get initial residuals
initResiduals!(_residuals_[ThreadID], Hs[ThreadID], _weights_[ThreadID]);
if ismissing(precs) 
    #Spectral deconvolution.
    solveHuber!(Hs[ThreadID], _residuals_[ThreadID], _weights_[ThreadID], 
                                        huber_δ, max_iter_outer = 100, max_iter_inner = 20, tol = Hs[ThreadID].n);
end
#Remember determined weights for eahc precursors
for i in range(1, IDtoCOL[ThreadID].size)
    precursor_weights[ThreadID][IDtoCOL[ThreadID].keys[i]] = _weights_[ThreadID][IDtoCOL[ThreadID][IDtoCOL[ThreadID].keys[i]]]# = precursor_weights[id]
end

if ismissing(isotope_dict) 
    getDistanceMetrics(_weights_[ThreadID], Hs[ThreadID], spectral_scores[ThreadID])
end

##########
#Scoring and recording data
if !ismissing(scored_PSMs)

    ScoreFragmentMatches!(unscored_PSMs[ThreadID],
                        IDtoCOL[ThreadID],
                        ionMatches[ThreadID], 
                        nmatches, 
                        err_dist)

    last_vals[ThreadID] = Score!(scored_PSMs[ThreadID], 
        unscored_PSMs[ThreadID],
        spectral_scores[ThreadID],
        _weights_[ThreadID],
        #match_count/prec_count,
        nmatches/(nmatches + nmisses),
        last_vals[ThreadID],
        Hs[ThreadID].n,
        Float32(sum(spectra[:intensities][i])), 
        i,
        min_spectral_contrast = min_spectral_contrast, #Remove precursors with spectral contrast lower than this ammount
        min_matched_ratio = min_matched_ratio,
        min_frag_count = min_frag_count, #Remove precursors with fewer fragments 
        min_weight = min_weight,
        min_topn = min_topn,
        block_size = 500000,
        )
end
end
#Add fragment matches to all_fmatches 
frag_err_idx = collectFragErrs(all_fmatches[ThreadID], ionMatches[ThreadID], nmatches, frag_err_idx, collect_fmatches)

##########
#Update Chromatograms 
#=
if !ismissing(chromatograms)
frag_counts = counter(UInt32) #Basically a Dictionary that counts the number of matched ions (values) for each precursor (keys)
for match_idx in range(1, nmatches) #fragmentMatches
    DataStructures.inc!(frag_counts, ionMatches[match_idx].prec_id)
end
#Add precursor templates with their weights and retention times to the chromatogram table 
chrom_idx = fillChroms!(chromatograms, IDtoCOL, chrom_idx, i, cycle_idx, prec_ids, prec_idx, frag_counts, weights, spectra[:retentionTime][i])
end
=#
##########
#Reset pre-allocated arrays 
reset!(ionTemplates[ThreadID], ion_idx)
reset!(ionMatches[ThreadID], nmatches)
reset!(ionMisses[ThreadID], nmisses)
fill!(prec_ids[ThreadID], zero(UInt32))
for i in range(1, Hs[ThreadID].n)
unscored_PSMs[ThreadID][i] = LXTandem(Float32)
end
reset!(IDtoCOL[ThreadID]);
reset!(Hs[ThreadID]);
end
end
#return fragment_intensities
############
#Return Chromatograms and Score/Feature Table
#return all_fmatches
if collect_fmatches
#return DataFrame(scored_PSMs), all_fmatches
return  DataFrame(@view(scored_PSMs[1][1:last_vals[1]])), all_fmatches[1]
else
return DataFrame(@view(scored_PSMs[1][1:last_vals[1]]))
end
#=
if ismissing(chromatograms)
#return all_fmatches
#return DataFrame(scored_PSMs)
DataFrame(@view(scored_PSMs[1:last_val]))
elseif ismissing(scored_PSMs)
chromatograms = DataFrame(chromatograms)
sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
return groupby(DataFrame(chromatograms), :precursor_idx)
else
chromatograms = DataFrame(chromatograms)
sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
return DataFrame(scored_PSMs), groupby(DataFrame(chromatograms), :precursor_idx)
end
=#
end

@everywhere begin 
    function testParallel(spectra::Arrow.Table, scan_idxs::Vector{Int64}, f_det::Vector{Vector{LibraryFragment{Float32}}})
        scan_sum = 0.0
        for scan_idx in scan_idxs
            for i in 1:1
                scan_sum += log2(sum(log2.(spectra[:intensities][scan_idx])))
            end
        end
        return scan_sum, f_det[1][1].frag_mz
    end
end

function distributeScans(N::Int64, m::Int64)
    scan_idx = 1
    scans = [ones(Int64, 1) for _ in range(1, Threads.nthreads())]
    proc = 0
    i = 1
    while scan_idx <= N
        push!(scans[proc%Threads.nthreads() + 1], scan_idx)
        if i > m
            i = 0
            proc += 1
        end
        i += 1
        scan_idx += 1
    end
    return scans
end

MS_TABLE = Arrow.Table(MS_TABLE_PATHS[1])   
distributed_scans = distributeScans(length(MS_TABLE[:masses]), 10000)
proc_sums = Vector{Float64}(undef, nprocs())

Threads.@threads for i in 1:length(distributed_scans)
    println(testParallel(MS_TABLE, distributed_scans[i], f_det), " + ", Threads.threadid())
end


pmap((distributed_scans)->testParallel(MS_TABLE, distributed_scans, f_det_flat_shared), [distributed_scans[i] for i in 1:4])


@distributed (append!) for i in 1:4
    testParallel(MS_TABLE, distributed_scans[i], f_det)
end

status = pmap(1:nprocs()) do i
    testParallel(MS_TABLE, distributed_scans[i])
end

f_det_sub = f_det[1:1000000]

pmap((distributed_scans)->testParallel(MS_TABLE, distributed_scans, f_det_flat_shared), [distributed_scans[i] for i in 1:4])

test_f_det_shared = [SharedArray(_) for _ in f_det]

test_f_det_shared = SharedArray(test_f_det_shared)


f_det_shared = SharedArray([SharedArray(sub_array) for sub_array in f_det])

status = pmap(1:nprocs()) do i
    testParallel(MS_TABLE, distributed_scans[i])
end