#=function findFirstFragmentBin(frag_index::Vector{FragBin{T}}, frag_min::T, frag_max::T) where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    potential_match = nothing 
    i = 1
    while lo <= hi
        i += 1
        mid = (lo + hi) ÷ 2
        #If either frag_min or frag_max are in the FragBin
        if ((frag_min > getLowMZ(frag_index[mid])) & (frag_min < getHighMZ(frag_index[mid]))) | (frag_max > getLowMZ(frag_index[mid])) & (frag_max < getHighMZ(frag_index[mid]))
            potential_match = mid
            hi = mid - 1
        #frag_max is below the FragBin
        elseif (frag_max) < getLowMZ(frag_index[mid])
            hi = mid - 1
        #frag_min is above the frag_bin 
        else
            lo = mid + 1
        end
    end
    println(i)
    return potential_match
end=#

function findFirstFragmentBin(frag_index::Vector{FragBin{T}}, frag_min::T, frag_max::T) where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    potential_match = nothing
    i = 1
    while lo <= hi
        i += 1
        mid = (lo + hi) ÷ 2
        #If either frag_min or frag_max are in the FragBin
        if ((frag_min > getLowMZ(frag_index[mid])) & (frag_min < getHighMZ(frag_index[mid])))
            potential_match = mid
            hi = mid - 1
        #frag_max is below the FragBin
        elseif (frag_max) < getLowMZ(frag_index[mid])
            hi = mid - 1
        #frag_min is above the frag_bin 
        else
            lo = mid + 1
        end
    end
    return potential_match#, Int64(getPrecBinID(frag_index[potential_match]))
end

function searchPrecursorBin!(precs::Dictionary{UInt32, IonIndexMatch{Float32}}, precursor_bin::PrecursorBin{T}, intensity::Float32, window_min::U, window_max::U) where {T,U<:AbstractFloat}
   
    N = getLength(precursor_bin)

    if N>1000
        return nothing, nothing
    end

    lo, hi = 1, N + 1

    while lo < hi
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

    window_stop = window_start

    prec_id = getPrecID(getPrecursor(precursor_bin, window_stop))

    if haskey(precs, prec_id)
        addMatch!(precs[prec_id], intensity)
    else
        insert!(precs, prec_id, IonIndexMatch(intensity, 1))
    end

    if (window_stop + 1) > N
        return
    end

    while (getPrecMZ(getPrecursor(precursor_bin, window_stop + 1)) < window_max)

        window_stop += 1
        prec_id = getPrecID(getPrecursor(precursor_bin, window_stop))
        if haskey(precs, prec_id)
            addMatch!(precs[prec_id], intensity)
        else
            insert!(precs, prec_id, IonIndexMatch(intensity, 1))
        end
        if (window_stop + 1) > N
            return #window_start, window_stop
        end

    end

    return window_stop, window_start#window_start, window_stop#window_stop
end

function searchPrecursorBin!(precs::Dictionary{UInt32, IonIndexMatch{Float32}}, precursor_bin::PrecursorBin{T}, intensity::Float32, window_min::U, window_max::U) where {T,U<:AbstractFloat}
   
    N = getLength(precursor_bin)

    if N>1000
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

    while lo < hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) > window_max
            hi = mid - 1
        else
            lo = mid + 1
        end
    end

    window_stop = hi

    for precursor_idx in window_start:window_stop

        prec_id = getPrecID(getPrecursor(precursor_bin, precursor_idx))

        if haskey(precs, prec_id)
            addMatch!(precs[prec_id], intensity)
        else
            insert!(precs, prec_id, IonIndexMatch(intensity, 1))
        end
    end

    return window_stop, window_start#window_start, window_stop#window_stop
end

function queryFragment!(precs::Dictionary{UInt32, IonIndexMatch{Float32}}, frag_index::FragmentIndex{T}, min_frag_bin::Int64, intensity::Float32, frag_min::U, frag_max::U, prec_min::U, prec_max::U) where {T,U<:AbstractFloat}
    first_frag_bin = findFirstFragmentBin(getFragBins(frag_index), frag_min, frag_max)
    if (first_frag_bin === nothing)
        return min_frag_bin
    end
    if (min_frag_bin == first_frag_bin)
        return min_frag_bin
    end
    while getLowMZ(getFragmentBin(frag_index, first_frag_bin)) < frag_max
        searchPrecursorBin!(precs, getPrecursorBin(frag_index, first_frag_bin), intensity, prec_min, prec_max)
        first_frag_bin += 1
    end
    return first_frag_bin - 1
end

FRAGMIN = 400.0
FRAGMAX = 400.05
PRECMIN = 500.0
PRECMAX = 508.0
BIN = 17957

precs = Dictionary{UInt32, Int64}()
@btime searchPrecursorBin(precs, getPrecursorBin(f_index, BIN), PRECMIN, PRECMAX)
precs = Dictionary{UInt32, Int64}()
@btime queryFragment!(precs, f_index, FRAGMIN, FRAGMAX, PRECMIN, PRECMAX)

MS_TABLE = Arrow.Table("/Users/n.t.wamsley/RIS_temp/ZOLKIND_MOC1_MAY23/parquet_out/MA5171_MOC1_DMSO_R01_PZ.arrow")


test_masses = MS_TABLE[:masses][10002]
test_intensities = MS_TABLE[:intensities][10002]

mutable struct IonIndexMatch{T<:AbstractFloat}
    summed_intensity::T
    count::Int64
end

function addMatch!(im::IonIndexMatch{T}, int::T) where {T<:AbstractFloat}
    im.summed_intensity += int
    im.count += 1
end

#min_frag_bin = 0
test = UInt32[]
@time begin
    for scan in 1:length(MS_TABLE[:msOrder])#[10001]
        if MS_TABLE[:msOrder] == 1
            continue
        end
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
        if length(precs) > 0
            push!(test, keys(sort(precs, by = x->x.count))[end])
        end
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