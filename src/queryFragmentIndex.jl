mutable struct IonIndexMatch{T<:AbstractFloat}
    summed_intensity::T
    count::Int64
end

getCount(i::IonIndexMatch{T}) where {T<:AbstractFloat} = i.count
getIntensity(i::IonIndexMatch{T}) where {T<:AbstractFloat} = i.summed_intensity

function getScore(i::IonIndexMatch{T}) where {T<:AbstractFloat}
    function logfac(N)
        N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
    end
    return getCount(i)#logfac(getCount(i)) + log(getIntensity(i))
end    

function addMatch!(im::IonIndexMatch{T}, int::T) where {T<:AbstractFloat}
    im.summed_intensity += int
    im.count += 1
end

function findFirstFragmentBin(frag_index::Vector{FragBin{T}}, frag_min::T, frag_max::T) where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    potential_match = nothing
    while lo <= hi

        mid = (lo + hi) ÷ 2

        if (frag_min) < getHighMZ(frag_index[mid])
            if (frag_max) > getHighMZ(frag_index[mid])
                potential_match = mid
            end
            hi = mid - 1
        elseif (frag_max) > getLowMZ(frag_index[mid])
            if (frag_min) < getLowMZ(frag_index[mid])
                potential_match = mid
            end
            lo = mid + 1
        end
    end

    return potential_match#, Int64(getPrecBinID(frag_index[potential_match]))
end

function searchPrecursorBin!(precs::Dictionary{UInt32, IonIndexMatch{Float32}}, precursor_bin::PrecursorBin{T}, intensity::Float32, window_min::U, window_max::U) where {T,U<:AbstractFloat}
   
    N = getLength(precursor_bin)

    if N>10000
        return nothing, nothing
    end

    lo, hi = 1, N

    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) < window_min
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    window_start = (lo <= N ? lo : return nothing, nothing)

    if getPrecMZ(getPrecursor(precursor_bin, window_start)) > window_max
        return nothing, nothing
    end

    lo, hi = window_start, N

    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) > window_max
            hi = mid - 1
        else
            lo = mid + 1
        end
    end

    window_stop = hi

    function addFragmentMatches!(precs::Dictionary{UInt32, IonIndexMatch{Float32}}, precursor_bin::PrecursorBin{T}, start::Int, stop::Int) where {T<:AbstractFloat}
        for precursor_idx in start:stop

            prec_id = getPrecID(getPrecursor(precursor_bin, precursor_idx))
    
            if haskey(precs, prec_id)
                addMatch!(precs[prec_id], intensity)
            else
                insert!(precs, prec_id, IonIndexMatch(intensity, 1))
            end
        end
    end

    addFragmentMatches!(precs, precursor_bin, window_start, window_stop)

    return

end

function queryFragment!(precs::Dictionary{UInt32, IonIndexMatch{Float32}}, frag_index::FragmentIndex{T}, min_frag_bin::Int64, intensity::Float32, frag_min::U, frag_max::U, prec_min::U, prec_max::U) where {T,U<:AbstractFloat}
    frag_bin = findFirstFragmentBin(getFragBins(frag_index), frag_min, frag_max)

    #No fragment bins contain the fragment m/z
    if (frag_bin === nothing)
        return min_frag_bin
    end

    while getLowMZ(getFragmentBin(frag_index, frag_bin)) < frag_max
        if frag_bin > min_frag_bin
            searchPrecursorBin!(precs, getPrecursorBin(frag_index, frag_bin), intensity, prec_min, prec_max)
        end
        frag_bin += 1

        if frag_bin > length(getFragBins(frag_index))
            break
        end
    end

    return frag_bin - 1
end

function searchScan!(precs::Dictionary{UInt32, IonIndexMatch{U}}, f_index::FragmentIndex{T}, massess::Vector{Union{Missing, U}}, intensities::Vector{Union{Missing, U}}, precursor_window::U, ppm::T, width::T; topN::Int = 100, min_frag_count::Int = 4) where {T,U<:AbstractFloat}
    
    getFragTol(mass::U, ppm::T) = mass*(1 - ppm/1e6), mass*(1 + ppm/1e6)

    function filterPrecursorMatches!(precs::Dictionary{UInt32, IonIndexMatch{T}}, topN::Int = 30, min_frag_count::Int = 4) where {T<:AbstractFloat}
        #Do not consider peptides wither fewer than 
        filter!(prec->getCount(prec)>min_frag_count, precs)
        sort!(precs, by = prec->getScore(prec), rev = true)
        #Iterator of Peptide ID's for the `topN` scoring peptides
        return Iterators.take(keys(precs), min(topN, length(keys(precs))))
        #=if length(precs) >0
            return first(pairs(precs))
        end=#
    end

    min_frag_bin = 0

    for (mass, intensity) in zip(massess, intensities)

        mass, intensity = coalesce(mass, 0.0),  coalesce(intensity, 0.0)

        FRAGMIN, FRAGMAX = getFragTol(mass, ppm) 

        min_frag_bin = queryFragment!(precs, f_index, min_frag_bin, intensity, FRAGMIN, FRAGMAX, precursor_window - width, precursor_window + width)
    end 

    return filterPrecursorMatches!(precs, topN, min_frag_count)
end

function selectTransitions(ptable::PrecursorTable, pep_ids::Base.Iterators.Take{Indices{UInt32}}, charges::Vector{UInt8}, isotopes::Vector{UInt8}; y_start::Int = 3, b_start::Int = 3, ppm::T = 20.0) where {T<:AbstractFloat}
    transitions = Vector{Transition}()
    for pep_id in pep_ids
        append!(transitions, getTransitions(getPep(ptable,pep_id), pep_id, charges, isotopes, y_start = y_start, b_start = b_start, ppm = ppm))
    end
    return sort!(transitions, by = x->getMZ(x))
end

function getTransitions(peptide::Peptide, pep_id::UInt32, charges::Vector{UInt8}, isotopes::Vector{UInt8}; y_start::Int = 3, b_start::Int = 3, ppm::T = 20.0) where {T<:AbstractFloat}
    getTransitions(Precursor(getSeq(peptide), prec_id = pep_id), charges, isotopes, y_start = y_start, b_start = b_start, ppm = ppm)
end

function SearchRAW(
                   spectra::Arrow.Table, 
                   ptable::PrecursorDatabase,
                   frag_index::FragmentIndex{T},
                   #selectTransitions, 
                   #right_precursor_tolerance::U,
                   #left_precursor_tolerance::U,
                   #transition_charges::Vector{UInt8},
                   #transition_isotopes::Vector{UInt8},
                   #b_start::Int64,
                   #y_start::Int64,
                   #fragment_match_ppm::U,
                   ms_file_idx::UInt32;
                   data_type::Type{T} = Float64
                   ) where {T,U<:Real}
    
    scored_PSMs = makePSMsDict(FastXTandem(data_type))
    #precursorList needs to be sorted by precursor MZ. 
    #Iterate through rows (spectra) of the .raw data. 
    i = 0
    ms2 = 0
    for spectrum in Tables.namedtupleiterator(spectra)
        i += 1
        if spectrum[:msOrder] != 2
            continue
        end
        ms2 += 1
        precs = Dictionary{UInt32, IonIndexMatch{Float32}}()
        pep_id_iterator = searchScan!(precs, frag_index, spectrum[:masses], spectrum[:intensities], spectrum[:precursorMZ], 20.0, 0.5)

        #transitions = selectTransitions(ptable, pep_id_iterator, UInt8[1, 2], UInt8[0, 1])
        transitions = selectTransitions(ptable, pep_id_iterator, UInt8[1], UInt8[0])

        fragmentMatches = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    #δs = params[:δs],
                                    δs = zeros(T, (1,)),#[Float64(0)],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx)

        unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)
        
        #Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(spectrum[:scanNumber]))
        Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(i))
    end
    println("processed $ms2 scans!")
    return scored_PSMs
end

@time test = SearchRAW(MS_TABLE, test_table, f_index, UInt32(1))
#=
PSMs = DataFrame(test)

transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(getPep(test_table, psm[:precursor_idx]))) => :decoy)
transform!(PSMs, AsTable(:) => ByRow(psm -> length(getSeq(getPep(test_table, psm[:precursor_idx])))) => :length)
PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, [:scan_idx])) 
decoy_scores = filter(row -> row.decoy, PSMs)[!, :hyperscore]
target_scores = filter(row -> !row.decoy, PSMs)[!, :hyperscore]

high_score = filter(row -> row.hyperscore >= 120, PSMs)#[!,!]
sum(high_score[!,:decoy])/size(high_score)[1]
PSMs = combine(sdf -> sdf[argmax(sdf.y_ladder), :], groupby(PSMs, [:scan_idx])) 
decoy_scores = filter(row -> row.decoy, PSMs)[!, :y_ladder]
target_scores = filter(row -> !row.decoy, PSMs)[!, :y_ladder]


println(mean(target_scores) - mean(decoy_scores))
MannWhitneyUTest(decoy_scores, target_scores)
histogram(decoy_scores, alpha = 0.5)#, normalize = :pdf)
histogram!(target_scores, alpha = 0.5)#, normalize = :pdf)

P

p = plot(layout=(1,2), size=(800,300))
targets = X[X_labels.==false,:]
decoys = X[X_labels.==true,:]
X = Matrix(PSMs[1:1:end,[:hyperscore,:total_ions,:y_ladder,:length,:error]])'
#X = transpose(Matrix(PSMs[!,[:hyperscore,:total_ions,:y_ladder,:length,:error]]))
X_labels = Vector(PSMs[1:1:end, :decoy])
#X_labels = [Dict(false=>0.0, true=>1.0)[x] for x in X_labels]
pca = fit(PCA, X; maxoutdim=2)
Ypca = predict(pca, X)
lda = fit(MulticlassLDA, X, X_labels; outdim=1)
Ylda = predict(lda, X)
p = plot(layout=(1,2), size=(800,300))

for s in [true, false]

    points = Ypca[:,X_labels.==s]
    scatter!(p[1], points[1,:],points[2,:], label=s, legend=:bottomleft)
    points = Ylda[:,X_labels.==s]
    scatter!(p[2], points[1,:],points[2,:], label=s, legend=:bottomleft)

end

plot!(p[1], title="PCA")
plot!(p[2], title="LDA")

histogram(Ylda[X_labels.==true], alpha = 0.5)
histogram!(Ylda[X_labels.==false], alpha = 0.5)

sum(Ylda[X_labels.==true].>0.05)
sum(Ylda[X_labels.==false].>0.05)

sum(Ylda[X_labels.==true].>0.08)/sum(Ylda[X_labels.==false].<-0.08)
#sum(Ylda[X_labels.==false].>0.08)

i = 0
count = 0
good_scans = 0
#precs = Dictionary{UInt32, IonIndexMatch{Float32}}()
best = []
len = 0
@time begin
    for scan in 1:length(MS_TABLE[:msOrder])#[4473]
        precs = Dictionary{UInt32, IonIndexMatch{Float32}}()
        if (MS_TABLE[:msOrder][scan] == 1) #& !ismissing(MS_TABLE[:precursorMZ])
            continue
        end
        i += 1
        #println(scan)
        len = searchScan!(precs, f_index, MS_TABLE[:masses][scan], MS_TABLE[:intensities][scan], MS_TABLE[:precursorMZ][scan], 20.0, 0.5)
        #len = length(filterPrecursorMatches!(precs, topN = 20, min_frag_count = 4))
        
        #if len != nothing
        #    push!(best, len)
        #end
        #if scan == 4473
        #    [println(x) for x in len]
        #end

    end
end

scores = []
for pair in best
    push!(scores, (isDecoy(test_table.id_to_pep[first(pair)]), getScore(last(pair))))
end

targets = filter(x->!first(x), scores)
decoys = filter(x->first(x), scores)
histogram([last(x) for x in decoys], alpha = 0.5)
histogram!([last(x) for x in targets], alpha = 0.5)
precs = Dictionary{UInt32, IonIndexMatch{Float32}}()
test_masses = MS_TABLE[:masses][scan]
test_intensities = MS_TABLE[:intensities][scan]
i = 1
for (mass, intensity) in zip(test_masses, test_intensities)
    mass = coalesce(mass, 0.0)
    FRAGMIN = mass - 10*(mass/1e6)
    FRAGMAX = mass + 10*(mass/1e6)
    min_frag_bin = queryFragment!(precs, f_index, min_frag_bin, intensity, FRAGMIN, FRAGMAX, 500.0, 508.0)
    i += 1

end

precs = Dictionary{UInt32, Int64}()
@btime searchPrecursorBin(precs, getPrecursorBin(f_index, BIN), PRECMIN, PRECMAX)
precs = Dictionary{UInt32, Int64}()
@btime queryFragment!(precs, f_index, FRAGMIN, FRAGMAX, PRECMIN, PRECMAX)

MS_TABLE = Arrow.Table("/Users/n.t.wamsley/RIS_temp/ZOLKIND_MOC1_MAY23/parquet_out/MA5171_MOC1_DMSO_R01_PZ.arrow")


test_masses = MS_TABLE[:masses][10002]
test_intensities = MS_TABLE[:intensities][10002]




test = UInt32[]
a = Dictionary{UInt32, IonIndexMatch{Float32}}()
@time begin
    for scan in [10001]#1:length(MS_TABLE[:msOrder])#[10001]
        if MS_TABLE[:msOrder] == 1
            continue
        end
        precs = Dictionary{UInt32, IonIndexMatch{Float32}}()
        test_masses = MS_TABLE[:masses][scan]
        test_intensities = MS_TABLE[:intensities][scan]
        i = 1
        min_frag_bin = 0
        for (mass, intensity) in zip(test_masses, test_intensities)
            mass = coalesce(mass, 0.0)
            FRAGMIN = mass - 10*(mass/1e6)
            FRAGMAX = mass + 10*(mass/1e6)
            min_frag_bin = queryFragment!(precs, f_index, min_frag_bin, intensity, FRAGMIN, FRAGMAX, 500.0, 508.0)
            println(length(filterPrecursorMatches!(precs, topN = 20, min_frag_count = 4)))
            i += 1

        end
        if length(precs) > 0
            push!(test, keys(sort(precs, by = x->x.count))[end])
        end
        a = precs
    end
end

searchPrecursorBin!(precs::Dictionary{UInt32, IonIndexMatch{Float32}}, precursor_bin::PrecursorBin{T}, intensity::Float32, window_min::U, window_max::U)

precs = Dictionary{UInt32, IonIndexMatch{Float32}}()
test_masses = MS_TABLE[:masses][43]
test_intensities = MS_TABLE[:intensities][43]
i = 1
for (mass, intensity) in zip(test_masses, test_intensities)
    mass = coalesce(mass, 0.0)
    FRAGMIN = mass - 10*(mass/1e6)
    FRAGMAX = mass + 10*(mass/1e6)
    #println("FRAGMIN ", FRAGMIN)
    #println("FRAGMAX ", FRAGMAX)
    #println(i)
    min_frag_bin = queryFragment!(precs, f_index, min_frag_bin, intensity, FRAGMIN, FRAGMAX, 500.0, 508.0)
    i += 1
end
#searchPrecursorBin!(precs, f_index.precursor_bins[44], Float32(10.0), 500.0, 501.0)
=#