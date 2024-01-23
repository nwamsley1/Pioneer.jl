function makeRTIndices(PSMS::DataFrame, 
                       precID_to_iRT::Dictionary{UInt32, Tuple{Float32, Float32}},
                       iRT_RT::Any;
                       min_prob::AbstractFloat = 0.75)
    #Maps filepath to a retentionTimeIndex (see buildRTIndex.jl)
    rt_indices = Dictionary{String, retentionTimeIndex{Float32, Float32}}()
    filter!(x->x.prob>=min_prob, PSMS);
    grouped_psms = groupby(PSMS,:file_path) 
    #iRT dict
    #iRT_dict = Dictionary{String, UnorderedDictionary{UInt32, Float32}}()
    for (groupkey, psms) in ProgressBar(pairs(grouped_psms)) #For each file in the experiment
        #insert!(iRT_dict, file_path, UnorderedDictionary{UInt32, Float32}())
        #Impute empirical iRT value for psms with probability lower than the threshold
        file_path = groupkey.file_path
        RTs, mzs, prec_ids = zeros(Float32, length(precID_to_iRT)), zeros(Float32, length(precID_to_iRT)), zeros(UInt32, length(precID_to_iRT))

        prec_set = Dictionary(psms[!,:precursor_idx],
                            zip(psms[!,:RT], 
                                #psms[!,:iRT]::Vector{Float32},
                                psms[!,:prec_mz])
                            )::Dictionary{UInt32, Tuple{Float32, Float32}}
        #Set of precursors where empirical iRTs can be used and do not need to be imputed 
        pset = Set(psms[!,:precursor_idx])
        #Loop through all precursors in the seed 
        Threads.@threads for (i, (prec_id, irt_mz_imputed)) in collect(enumerate(pairs(precID_to_iRT)))
            prec_ids[i] = prec_id
            if prec_id ∈ pset::Set{UInt32} #Use empirical iRT/RT
                rt, mz = prec_set[prec_id]::Tuple{Float32, Float32}
                RTs[i] = rt
                mzs[i] = mz
            else #Impute empirical iRT from the best observed psm for the precursor accross the experiment 
                irt, mz = irt_mz_imputed::Tuple{Float32, Float32}
                rt = iRT_RT[file_path](irt)::Float64 #Convert iRT to RT 
                RTs[i] = rt
                mzs[i] = mz
                #insert!(iRT_dict[file_path], prec_id, abs(irt - prec_set[prec_id][2]))
            end
        end
        
        #Build RT index 
        rt_df = DataFrame(Dict(:RT => RTs,
                                :prec_mz => mzs,
                                :precursor_idx => prec_ids))
        sort!(rt_df, :RT)
        rt_index = buildRTIndex(rt_df);
        insert!(rt_indices, file_path, rt_index)
    end
    return rt_indices
end

function mapRTandiRT(psms_dict::Dictionary{String, DataFrame}; min_prob::AbstractFloat = 0.9)

    #Dictionaries mapping fild_id names to data
    prec_ids = Dictionary{String, Set{UInt32}}() #File name => precursor id
    iRT_RT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    RT_iRT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    for (key, psms) in ProgressBar(pairs(psms_dict)) #For each data frame 
        psms[!,:file_path] .= key #Add column for file name

        insert!(prec_ids,key,Set(psms[!,:precursor_idx])) #Set of precursors ids in the dataframe

        best_hits = psms[!,:prob].>min_prob#Map RTs using only the best psms
        
        #Build RT=>iRT and iRT=> RT mappings for the file and add to the dictionaries 
        insert!(iRT_RT, key, KDEmapping(psms[best_hits,:iRT_predicted],
        psms[best_hits,:RT]
        ))
        insert!(RT_iRT, key, KDEmapping(psms[best_hits,:RT],
        psms[best_hits,:iRT_predicted]
        ))
    end
    return iRT_RT, RT_iRT
end

function getPrecIDtoiRT(psms_dict::Dictionary{String, DataFrame}, RT_iRT::Any; max_q_value::AbstractFloat = 0.1)
    psms = vcat(values(psms_dict)...); #Combine data from all files in the experiment
    #Only consider precursors with a q_value passing the threshold
    # at least once accross the entire experiment
    #filter!(x->x.q_value<=max_q_value,psms); 

    #For each precursor, what was the best scan, maximum q-value, and observed iRT for the best scan
    psms[!,:best_scan] .= false #Best scan for a given precursor_idx accross the experiment
    psms[!,:iRT_observed] .= zero(Float64) #Observed iRT for each psm
    psms[!,:iRT_diff_from_median] .= zero(Float64) #Observed iRT for each psm
    #psms[!,:iRT_observed_best] .= zero(Float64) #Best iRT accross the experiment for a given precursor
    grouped_psms = groupby(psms, :precursor_idx); #Groupby precursor idx

    precIDtoiRT = Dictionary{UInt32, Tuple{Float32, Float32}}()
    for i in ProgressBar(range(1, length(grouped_psms))) #For each precursor 

        #Best scan for the precursor accross the experiment
        #best_scan = argmax(grouped_psms[i].prob) 
        #grouped_psms[i].best_scan[best_scan] = true #Makr columsn containing the best psm for a precursor
        #irt = RT_iRT[grouped_psms[i].file_path[best_scan]](grouped_psms[i].RT[best_scan])
        #grouped_psms[i].iRT_observed[best_scan] = irt
        #Could alternatively estimate the iRT using the best-n scans 
        for j in range(1, size(grouped_psms[i])[1])
            #Get empirical iRT for the psm given the empirical RT
            irt = RT_iRT[grouped_psms[i].file_path[j]](grouped_psms[i].RT[j])
            grouped_psms[i].iRT_observed[j] = irt
        end

        median_irt = median(grouped_psms[i].iRT_observed)
        
        for j in range(1, size(grouped_psms[i])[1])
            grouped_psms[i].iRT_diff_from_median[j] = abs(median_irt - grouped_psms[i].iRT_observed[j])
        end

        insert!(precIDtoiRT, #dict
                first(grouped_psms[i].precursor_idx), #key
                (Float32(median_irt),  first(grouped_psms[i].prec_mz)) #value
                )
    end

    #filter!(x->x.best_scan, psms); #Filter out non-best scans 
    #Map the precursor idx to the best empirical iRT extimate and the precursor m/z
    return precIDtoiRT, psms #Dictionary(psms[!,:precursor_idx], zip(psms[!,:iRT_observed], psms[!,:prec_mz]))
end


function KZfilter!(X::Vector{R}, m::Int64, k::Int64) where {R<:Real}
    #Kolmogorov-Zurbenko filter
    #zero padding on both sides   
    #implement as k applications of a moving average filter

    #k - number of MA applications
    #m - width of filter
    T = length(X)
    for n in range(1, k)
        for t in range(1, T)
            #w = (m - 1)÷2
            S = zero(T)
            start, stop = max(1, t - m), min(T, t + m)
            for s in range(start, stop)
                S += X[s]
            end
            X[t] = S/(2*m + 1)
        end
    end
end

N = 10000
N = 1

prec_id = precs[N]
MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)][!,[:precursor_idx,:decoy,:RT,:weight,:total_ions,:topn,:entropy_score]]

CHROM = MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)]
plot(CHROM[!,:RT],
CHROM[!,:weight],
seriestype=:scatter)
weights = copy(CHROM[!,:weight]);
KZfilter!(weights, 2, 3);
plot!(CHROM[!,:RT],
weights,
seriestype=:scatter)
weights = copy(CHROM[!,:weight]);
KZfilter!(weights, 1,3);
plot!(CHROM[!,:RT],
weights,
seriestype=:scatter)

N += 1

precs = unique(best_psms[(best_psms[!,:q_value].<0.01).&(best_psms[!,:decoy].==false),:precursor_idx])
prec_id = precs[N]

dtype = Float32
gx, gw = gausslegendre(100)
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    )
reset!(state)
#prec_id = 10429654
#integratePrecursorMS2(MS2_CHROMS_GROUPED[N],
prec_id = precs[N]
reset!(state)
integratePrecursorMS2(MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)],
                        state,
                        gx,
                        gw,
                        intensity_filter_fraction = 0.01f0,
                        α = 0.001f0,
                        half_width_at_α = 0.15f0,
                        LsqFit_tol = 1e-3,
                        Lsq_max_iter = 100,
                        tail_distance = 0.25f0,
                        isplot = true
)
MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)][!,[:precursor_idx,:total_ions,:isotope_count,:entropy_score,
:scribe,:spectral_contrast,:weight,:RT,:scan_idx,:isotopes,:matched_ratio,:minMZ,:prec_mz,:maxMZ,:scanMZ,:iso_rank,:decoy]]

#good_scans = (abs.(MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)][!,[:maxMZ]] .- 616.53).<0.1)[!,:maxMZ];
#good_scans = (MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)][!,:isotopes] .== (0.0f0, 5.0f0))[!,:isotopes];

#good_scans = MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)][!,:iso_rank] .== 1
#good_scans = [true, true, true, false, true, false, false, true, false, true, false]
#MS2_CHROMS_GROUPED[(precursor_idx = prec_id,)][good_scans,[:precursor_idx,:total_ions,:isotope_count,:entropy_score,
#:scribe,:spectral_contrast,:weight,:RT,:scan_idx,:matched_ratio,:isotopes,:minMZ,:prec_mz,:maxMZ,:scanMZ,:decoy]]
N += 1

#MS2_CHROMS[!,:isotopes] .= [(0.0, 0.0) for x in range(1, size(MS2_CHROMS, 1))]
MS2_CHROMS[!,:iso_rank] .= zero(UInt8)
for i in ProgressBar(range(1, size(MS2_CHROMS,1)))
    charge = MS2_CHROMS[i,:charge]
    mz = MS2_CHROMS[i,:prec_mz]
    window = (Float32(MS2_CHROMS[i,:minMZ]), Float32(MS2_CHROMS[i,:maxMZ]))
    isotopes = getPrecursorIsotopeSet(mz, charge, window)
    rank = zero(UInt8)
    if iszero(first(isotopes))
        if last(isotopes) > 1
            rank = UInt8(1)
        elseif last(isotopes) == 1
            rank = UInt8(2)
        else
            rank = UInt8(3)
        end
    else
        rank = UInt8(4)
    end
    MS2_CHROMS[i,:iso_rank] = rank
end

MS2_CHROMS[!,:scanMZ] = [MS_TABLE[:precursorMZ][scan_id] for scan_id in MS2_CHROMS[!,:scan_idx]]
MS2_CHROMS[!,:maxMZ] = MS2_CHROMS[!,:scanMZ] .+ 8.0036/2
MS2_CHROMS[!,:minMZ] = MS2_CHROMS[!,:scanMZ] .-8.0036/2

N +=1
test_psms = vcat(values(PSMs_Dict)...);
test = groupby(test_psms,:precursor_idx)

Threads.@threads for i in ProgressBar(range(1, length(test_psms[!,:file_path])))
    words = split(split(test_psms[i,:file_path], "\\")[end], "_")
    test_psms[i,:name] = words[5]*"_"*words[7]*"_"*words[8]
end
test[N][!,[:precursor_idx,:q_value,:topn,:y_count,:b_count,:entropy_score,:RT,:iRT_observed,:name]]
describe(test[N][!,[:iRT_observed]])
N += 1

grouped_psms = groupby(best_psms, :precursor_idx)
for i in ProgressBar(range(1, length(grouped_psms)))
    iso_ranks = unique(grouped_psms[i][!,:iso_rank])
    if length(iso_ranks)<=1
        continue
    end
    best_rank, best_score = 0, 0
    score = 0
    best_rank = minimum(iso_ranks)
    
    for rank in iso_ranks
        for j in range(1, size(grouped_psms[i], 1))
            if grouped_psms[i][j,:iso_rank] == rank
                #if grouped_psms[i][j,:q_value] <=0.01
                    score += grouped_psms[i][j,:weight]
                #end
            end
        end
        if score == best_score
            if rank < best_rank
                best_score = score
                best_rank = rank
            end
            score = 0
        elseif score > best_score
            best_score = score
            best_rank = rank
            score = 0
        else
            score = 0
        end
    end
    
    for j in range(1, size(grouped_psms[i], 1))
        if grouped_psms[i][j,:iso_rank] != best_rank
            grouped_psms[i][j,:best_scan] = false
        end
    end
end

N = 10000
sort(grouped_psms[N], :iso_rank)[!,[:sequence,:iso_rank,:decoy,:q_value,:H,:best_scan]]
N += 1
#=
PSMs_Dict = load("C:\\Users\\n.t.wamsley\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\PSMs_Dict_isotopes_010224_isotopes_M0M1.jld2")["PSMs_Dict"];
iRT_RT, RT_iRT = mapRTandiRT(PSMs_Dict)
precID_to_iRT = getPrecIDtoiRT(PSMs_Dict, RT_iRT)
RT_INDICES = makeRTIndices(PSMs_Dict,precID_to_iRT,iRT_RT)
#N += 1
#precID_to_iRT[N][!,[:sequence,:decoy,:prob,:q_value,:iRT,:iRT_observed,:RT]]
#println("mean ", mean(precID_to_iRT[N][!,:iRT_observed]))
#println("median ", median(precID_to_iRT[N][!,:iRT_observed]))
#println("std ", std(precID_to_iRT[N][!,:iRT_observed]))
#

test_sub = last(test,10000)[!,[:sequence,:prob,:RT,:iRT]]
plot(test_sub[!,:iRT], test_sub[!,:RT], seriestype=:scatter, alpha = 0.1)
test_map = KDEmapping(test_sub[!,:iRT], test_sub[!,:RT])
plot!(collect(LinRange(-50, 150, 100)), [test_map(_) for _ in LinRange(-50, 150, 100)])
PSMs_Dict = load("C:\\Users\\n.t.wamsley\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\PSMs_Dict_isotopes_010224_isotopes_M0M1.jld2")["PSMs_Dict"];
iRT_RT, RT_iRT = mapRTandiRT(PSMs_Dict)
prec_ID_to_iRT = getPrecIDtoiRT(PSMs_Dict, RT_iRT)

=#
#RT_INDICES = makeRTIndices(PSMs_Dict,precID_to_iRT,iRT_RT)

PSMs_Dict = load("C:\\Users\\n.t.wamsley\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\PSMs_Dict_isotopes_010224_isotopes_M0M1.jld2")["PSMs_Dict"];

psms = vcat(values(psms_dict)...); #Combine data from all files in the experiment
grouped_psms = groupby(psms, :precursor_idx); #Groupby precursor idx
psms[!,:best_scan] .= false
for i in ProgressBar(range(1, length(grouped_psms))) #For each precursor 

    #Best scan for the precursor accross the experiment
    best_scan = argmax(grouped_psms[i].prob) 
    grouped_psms[i].best_scan[best_scan] = true #Makr columsn containing the best psm for a precursor
end
filter!(x->x.best_scan, psms); #Fil


psms[1:250000,:]
"C:\\Users\\n.t.wamsley\\Desktop\\good_precursors_a.jld2"
spectronaut_precs = load("C:\\Users\\n.t.wamsley\\Desktop\\good_precursors_a.jld2")["combined_set"];
SEQS = "_".*psms[!,:sequence].*"_.".*string.(psms[!,:charge]);
SEQS = replace.(SEQS, "C" => "C[Carbamidomethyl (C)]");
psms[!,:modified_sequence] = SEQS;


for i in (50000, 100000, 200000, 250000, 500000)
    println(length(setdiff(spectronaut_precs, Set(psms[1:i,:modified_sequence]))))
end