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
Sort state management for file references.

Provides utilities for tracking, validating, and managing
the sort state of files in the FileOperations system.
"""

# Forward declaration - FileReference is defined in FileReferences.jl
# This file will be included after FileReferences.jl is loaded

#==========================================================
Sort State Management
==========================================================#

"""
    mark_sorted!(ref::FileReference, keys::Symbol...)

Mark a file reference as sorted by the specified keys.
"""
function mark_sorted!(ref::FileReference, keys::Symbol...)
    ref.sorted_by = keys
    return ref
end

"""
    is_sorted_by(ref::FileReference, keys::Symbol...) -> Bool

Check if a file is sorted by the specified keys in the exact order.
"""
function is_sorted_by(ref::FileReference, keys::Symbol...)
    return sorted_by(ref) == keys
end

"""
    ensure_sorted!(ref::FileReference, keys::Symbol...)

Ensure a file is sorted by the specified keys, sorting if necessary.
"""
function ensure_sorted!(ref::FileReference, keys::Symbol...)
    if !is_sorted_by(ref, keys...)
        sort_file_by_keys!(ref, keys...)
    end
    return ref
end

"""
    Base.sort!(ref::FileReference, cols::Vector{Symbol}; kwargs...)

Override Base.sort! to automatically track sort state for FileReferences.
This ensures that any sorting operation updates the FileReference metadata.
"""
function Base.sort!(ref::FileReference, cols::Vector{Symbol}; 
                   rev::Union{Bool, Vector{Bool}}=false, kwargs...)
    # Note: This is a placeholder that delegates to sort_file_by_keys!
    # The actual implementation is in io/ArrowOperations.jl
    # We just ensure the interface is consistent with Base.sort!
    
    # Convert single bool to appropriate format for sort_file_by_keys!
    reverse = if rev isa Bool
        rev  # sort_file_by_keys! expects a single bool
    else
        # If mixed directions, we can't use sort_file_by_keys! directly
        error("Mixed sort directions not yet supported. Use all ascending or all descending.")
    end
    
    # Delegate to the existing function
    sort_file_by_keys!(ref, cols...; reverse=reverse)
    
    return ref
end

# Export sort state functions
export mark_sorted!, is_sorted_by, ensure_sorted!