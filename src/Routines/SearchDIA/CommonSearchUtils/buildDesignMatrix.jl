# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    MzGroupingMap

Data structure for managing m/z-based grouping of precursors in MS1 design matrix construction.
Groups precursors with identical (rounded) m/z values into shared design matrix columns.
"""
mutable struct MzGroupingMap
    # Maps rounded m/z → design matrix column
    mz_to_col::Dict{UInt32, UInt16}

    # Maps precursor ID → rounded m/z (for coefficient retrieval)
    precid_to_mz_group::Dict{UInt32, UInt32}

    # Maps rounded m/z → list of precursor IDs in that group
    mz_group_to_precids::Dict{UInt32, Vector{UInt32}}

    # Precision for m/z rounding (default: 100000 for 5 decimal places)
    mz_precision::UInt32

    # Current column counter
    current_col::UInt16
end

"""
    MzGroupingMap(mz_precision::UInt32 = 100000)

Constructor for MzGroupingMap with configurable m/z precision.
Default precision of 100000 groups masses within ~0.01 ppm at 1000 Da.
"""
function MzGroupingMap(mz_precision::UInt32 = UInt32(100000))
    return MzGroupingMap(
        Dict{UInt32, UInt16}(),
        Dict{UInt32, UInt32}(),
        Dict{UInt32, Vector{UInt32}}(),
        mz_precision,
        UInt16(0)
    )
end

"""
    round_mz_for_grouping(mz::Float32, precision::UInt32) -> UInt32

Round m/z to specified precision for grouping identical masses.
Example: round_mz_for_grouping(1000.12345f0, 100000) = 100012345
"""
function round_mz_for_grouping(mz::Float32, precision::UInt32)::UInt32
    return round(UInt32, mz * precision)
end

"""
    reset!(grouping::MzGroupingMap)

Reset the grouping map for reuse across scans.
"""
function reset!(grouping::MzGroupingMap)
    empty!(grouping.mz_to_col)
    empty!(grouping.precid_to_mz_group)
    empty!(grouping.mz_group_to_precids)
    grouping.current_col = UInt16(0)
    return nothing
end

function buildDesignMatrix!(H::SparseArray{UInt32,Float32}, 
                            matches::Vector{m},  
                            misses::Vector{m}, 
                            nmatches::Int64, 
                            nmisses::Int64, 
                            precID_to_col::ArrayDict{UInt32, UInt16}; 
                            block_size = 10000) where {m<:MatchIon{Float32}}
    T = Float32
    #Number of rows equals the number of unique matched peaks
    #Remember "getPeakInd(x)" is hte index of the matched peak in the MS2 spectrum.
    M = 1
    for i in range(2, nmatches)
        if getPeakInd(matches[i])!= getPeakInd(matches[i - 1])
            M += 1
        end
    end
    M += nmisses

    #If M exceeds pre-allocated size
    if (nmatches + nmisses) >= length(H.colval) - 1
        block_size = max(block_size,nmatches + nmisses - length(H.colval))
        append!(H.colval, zeros(eltype(H.colval), block_size))
        append!(H.rowval, zeros(eltype(H.rowval), block_size))
        append!(H.nzval, zeros(eltype(H.nzval), block_size))
        append!(H.x, zeros(eltype(H.x), block_size))
        append!(H.matched, zeros(eltype(H.matched), block_size))
        append!(H.isotope, zeros(eltype(H.isotope), block_size))
        append!(H.colptr, Vector{UInt32}(undef, block_size))
    end
    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    #X = zeros(T, M)
    #Maps a precursor id to a row of H. 

    #Current highest row encountered
    prec_col = zero(UInt32)
    row = 0
    #Number of unique peaks encountered. 
    last_peak_ind = -1
    last_row, last_col = -1, -1
    j = 0
    H.n_vals = 0
    for i in range(1, nmatches)#matches)
        match = matches[i]
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if iszero(precID_to_col[getPrecID(match)])
            prec_col += one(UInt32)
            update!(precID_to_col, getPrecID(match), UInt16(prec_col))
        end

        #If this peak has not been encountered yet, then start filling a new 
        if getPeakInd(match) != last_peak_ind
            row += 1
            last_peak_ind = getPeakInd(match)
        end

        col =  precID_to_col[getPrecID(match)]
        if (col != last_col) | (row != last_row)
            j += 1
            H.n_vals += 1
        end
        H.colval[j] = col
        H.rowval[j] = row
        H.nzval[j] += getPredictedIntensity(match)
        H.x[j] = getIntensity(match)
        H.matched[j] = true
        H.isotope[j] = getIsoIdx(match)

        last_col = col
        last_row = row

    end

    i = j + 1#nmatches + 1
    for j in range(1, nmisses)#range(nmatches + 1, nmatches + nmisses)
        miss = misses[j]#j - nmatches]
        #If a match for this precursor hasn't been encountered yet, then there were no matched ions for this
        #precursor and it should not be included in the design matrix. 
        if !iszero(precID_to_col[getPrecID(miss)])
            
            row += 1
            H.colval[i] = precID_to_col[getPrecID(miss)]
            H.rowval[i] = row
            H.nzval[i] += getPredictedIntensity(miss)
            H.x[i] = zero(Float32)
            H.matched[i] = false
            H.isotope[i] = getIsoIdx(miss)
            i += 1
            H.n_vals += 1
        end
    end
    sortSparse!(H)
end

"""
    buildDesignMatrixMS1!(H::SparseArray{UInt32,Float32},
                          matches::Vector{m}, misses::Vector{m},
                          nmatches::Int64, nmisses::Int64,
                          mz_grouping::MzGroupingMap,
                          precursors::LibraryPrecursors;
                          block_size = 10000) where {m<:MatchIon{Float32}}

MS1-specific design matrix construction that groups precursors by m/z values.
Identical to buildDesignMatrix! except precursors with identical m/z are
assigned to the same design matrix column to avoid multicollinearity.
"""
function buildDesignMatrixMS1!(
    H::SparseArray{UInt32,Float32},
    matches::Vector{m},
    misses::Vector{m},
    nmatches::Int64,
    nmisses::Int64,
    mz_grouping::MzGroupingMap,
    precursors::LibraryPrecursors;
    block_size = 10000
) where {m<:MatchIon{Float32}}

    T = Float32

    # Get precursor m/z values for grouping
    precursor_mzs = getMz(precursors)

    # Calculate number of rows (same logic as original buildDesignMatrix!)
    M = 1
    for i in range(2, nmatches)
        if getPeakInd(matches[i]) != getPeakInd(matches[i - 1])
            M += 1
        end
    end
    M += nmisses

    # Handle array resizing (same logic as original)
    if (nmatches + nmisses) >= length(H.colval) - 1
        block_size = max(block_size, nmatches + nmisses - length(H.colval))
        append!(H.colval, zeros(eltype(H.colval), block_size))
        append!(H.rowval, zeros(eltype(H.rowval), block_size))
        append!(H.nzval, zeros(eltype(H.nzval), block_size))
        append!(H.x, zeros(eltype(H.x), block_size))
        append!(H.matched, zeros(eltype(H.matched), block_size))
        append!(H.isotope, zeros(eltype(H.isotope), block_size))
        append!(H.colptr, Vector{UInt32}(undef, block_size))
    end

    # Initialize tracking variables
    row = 0
    last_peak_ind = -1
    last_row, last_col = -1, -1
    j = 0
    H.n_vals = 0

    # Process matches with m/z grouping
    for i in range(1, nmatches)
        match = matches[i]
        prec_id = getPrecID(match)
        prec_mz = precursor_mzs[prec_id]

        # Get or create m/z group column
        mz_group = round_mz_for_grouping(prec_mz, mz_grouping.mz_precision)

        if !haskey(mz_grouping.mz_to_col, mz_group)
            # Create new column for this m/z group
            mz_grouping.current_col += 1
            mz_grouping.mz_to_col[mz_group] = mz_grouping.current_col
            mz_grouping.mz_group_to_precids[mz_group] = UInt32[]
        end

        # Track precursor → m/z group mapping
        mz_grouping.precid_to_mz_group[prec_id] = mz_group

        # Add precursor to group if not already present
        if prec_id ∉ mz_grouping.mz_group_to_precids[mz_group]
            push!(mz_grouping.mz_group_to_precids[mz_group], prec_id)
        end

        # Use m/z group column instead of individual precursor column
        col = mz_grouping.mz_to_col[mz_group]

        # Handle peak indexing (same as original)
        if getPeakInd(match) != last_peak_ind
            row += 1
            last_peak_ind = getPeakInd(match)
        end

        # Add entry to design matrix (same as original)
        if (col != last_col) | (row != last_row)
            j += 1
            H.n_vals += 1
        end
        H.colval[j] = col
        H.rowval[j] = row
        H.nzval[j] += getPredictedIntensity(match)
        H.x[j] = getIntensity(match)
        H.matched[j] = true
        H.isotope[j] = getIsoIdx(match)

        last_col = col
        last_row = row
    end

    # Process misses with m/z grouping
    i = j + 1
    for j_miss in range(1, nmisses)
        miss = misses[j_miss]
        prec_id = getPrecID(miss)

        # Only add miss if precursor was seen in matches (same logic as original)
        if haskey(mz_grouping.precid_to_mz_group, prec_id)
            mz_group = mz_grouping.precid_to_mz_group[prec_id]
            col = mz_grouping.mz_to_col[mz_group]

            row += 1
            H.colval[i] = col
            H.rowval[i] = row
            H.nzval[i] += getPredictedIntensity(miss)
            H.x[i] = zero(Float32)
            H.matched[i] = false
            H.isotope[i] = getIsoIdx(miss)
            i += 1
            H.n_vals += 1
        end
    end

    # Sort sparse matrix (same as original)
    sortSparse!(H)
end

"""
    distribute_ms1_coefficients!(precursor_weights::Vector{Float32},
                                grouped_weights::Vector{Float32},
                                mz_grouping::MzGroupingMap)

Distribute fitted coefficients from m/z groups back to individual precursors.
After solving the grouped design matrix, this assigns the fitted coefficient
for each m/z group to all precursors in that group.

Args:
- precursor_weights: Array indexed by precursor ID to store individual weights
- grouped_weights: Array indexed by column number containing fitted group coefficients
- mz_grouping: Mapping between precursors and m/z groups
"""
function distribute_ms1_coefficients!(
    precursor_weights::Vector{Float32},
    grouped_weights::Vector{Float32},
    mz_grouping::MzGroupingMap
)
    # Iterate through all m/z groups and distribute their coefficients
    for (mz_group, col_idx) in mz_grouping.mz_to_col
        # Get the fitted coefficient for this m/z group
        group_coefficient = col_idx <= length(grouped_weights) ? grouped_weights[col_idx] : 0.0f0

        # Assign this coefficient to all precursors in the group
        if haskey(mz_grouping.mz_group_to_precids, mz_group)
            for prec_id in mz_grouping.mz_group_to_precids[mz_group]
                # Make sure the array is large enough
                if prec_id <= length(precursor_weights)
                    precursor_weights[prec_id] = group_coefficient
                end
            end
        end
    end
end
