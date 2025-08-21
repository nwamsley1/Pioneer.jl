# WarningTracker.jl
# Manages warning collection and categorization

export WarningTracker, WarningEntry, WarningCategory
export add_warning!, get_warnings, get_warnings_by_category, clear_warnings!
export categorize_warning

# Warning categories
@enum WarningCategory begin
    MASS_ERROR
    RT_CALIBRATION
    MISSING_DATA
    PARAMETER_TUNING
    PROTEIN_INFERENCE
    QUANTIFICATION
    FILE_IO
    PERFORMANCE
    CONFIGURATION
    OTHER
end

struct WarningEntry
    timestamp::DateTime
    message::String
    category::WarningCategory
    location::Union{Nothing, String}
    context::Dict{Symbol, Any}
end

mutable struct WarningTracker
    warnings::Vector{WarningEntry}
    category_counts::Dict{WarningCategory, Int}
    lock::ReentrantLock
    max_warnings::Int
    
    function WarningTracker(; max_warnings::Int = 10000)
        new(
            WarningEntry[],
            Dict{WarningCategory, Int}(),
            ReentrantLock(),
            max_warnings
        )
    end
end

"""
    add_warning!(tracker::WarningTracker, message::String; 
                 category::WarningCategory = OTHER,
                 location::Union{Nothing, String} = nothing,
                 context::Dict{Symbol, Any} = Dict{Symbol, Any}())

Add a warning to the tracker with automatic categorization if not specified.
"""
function add_warning!(tracker::WarningTracker, message::String; 
                     category::Union{Nothing, WarningCategory} = nothing,
                     location::Union{Nothing, String} = nothing,
                     context::Dict{Symbol, Any} = Dict{Symbol, Any}())
    
    # Auto-categorize if not specified
    if category === nothing
        category = categorize_warning(message)
    end
    
    entry = WarningEntry(
        now(),
        message,
        category,
        location,
        context
    )
    
    lock(tracker.lock) do
        # Check if we've hit the limit
        if length(tracker.warnings) >= tracker.max_warnings
            # Keep only the most recent warnings
            deleteat!(tracker.warnings, 1:100)  # Remove oldest 100
        end
        
        push!(tracker.warnings, entry)
        tracker.category_counts[category] = get(tracker.category_counts, category, 0) + 1
    end
    
    return entry
end

"""
    categorize_warning(message::String) -> WarningCategory

Automatically categorize a warning based on its message content.
"""
function categorize_warning(message::String)
    msg_lower = lowercase(message)
    
    # Check for specific patterns
    if occursin(r"mass\s*(error|tolerance|accuracy)", msg_lower) || 
       occursin(r"ppm|dalton|m/z", msg_lower)
        return MASS_ERROR
    elseif occursin(r"rt|retention\s*time|calibration|alignment", msg_lower)
        return RT_CALIBRATION
    elseif occursin(r"missing|not\s*found|empty|no\s*data", msg_lower)
        return MISSING_DATA
    elseif occursin(r"parameter|tuning|optimization|convergence", msg_lower)
        return PARAMETER_TUNING
    elseif occursin(r"protein|inference|grouping|parsimony", msg_lower)
        return PROTEIN_INFERENCE
    elseif occursin(r"quant|intensity|abundance|maxlfq", msg_lower)
        return QUANTIFICATION
    elseif occursin(r"file|read|write|i/o|permission", msg_lower)
        return FILE_IO
    elseif occursin(r"performance|memory|slow|timeout", msg_lower)
        return PERFORMANCE
    elseif occursin(r"config|setting|parameter|option", msg_lower)
        return CONFIGURATION
    else
        return OTHER
    end
end

"""
    get_warnings(tracker::WarningTracker) -> Vector{WarningEntry}

Get all warnings from the tracker.
"""
function get_warnings(tracker::WarningTracker)
    lock(tracker.lock) do
        return copy(tracker.warnings)
    end
end

"""
    get_warnings_by_category(tracker::WarningTracker, category::WarningCategory) -> Vector{WarningEntry}

Get warnings filtered by category.
"""
function get_warnings_by_category(tracker::WarningTracker, category::WarningCategory)
    lock(tracker.lock) do
        return filter(w -> w.category == category, tracker.warnings)
    end
end

"""
    clear_warnings!(tracker::WarningTracker)

Clear all warnings from the tracker.
"""
function clear_warnings!(tracker::WarningTracker)
    lock(tracker.lock) do
        empty!(tracker.warnings)
        empty!(tracker.category_counts)
    end
end

"""
    get_warning_summary(tracker::WarningTracker) -> Dict

Get a summary of warnings by category.
"""
function get_warning_summary(tracker::WarningTracker)
    lock(tracker.lock) do
        return Dict(
            :total => length(tracker.warnings),
            :by_category => copy(tracker.category_counts),
            :recent => tracker.warnings[max(1, end-9):end]  # Last 10 warnings
        )
    end
end