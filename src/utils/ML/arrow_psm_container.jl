#=
ArrowFilePSMContainer: file-backed PSM container for out-of-memory scoring.

Each (ms_file_idx, cv_fold) pair is stored as a separate Arrow file.
Mutable scoring columns (trace_prob, q_value, MBR_*) live in a separate
"sidecar" file so scoring iterations only rewrite ~32 bytes/PSM instead
of the full ~300-500 bytes/PSM of immutable feature data.

File layout after prepare_training_data!:
  {base}_fold0.arrow          ← sorted, with pair_id + irt_bin_idx (immutable)
  {base}_fold0_scores.arrow   ← trace_prob, q_value, MBR_* columns (mutable)
=#

#############################################################################
# ArrowFileGroup — one (ms_file_idx, cv_fold) pair
#############################################################################

"""
    ArrowFileGroup

Represents a single (ms_file_idx, cv_fold) Arrow file pair:
- `data_path`: Main data file (immutable after `prepare_training_data!`)
- `scores_path`: Sidecar with mutable scoring columns
- `n_rows`: Cached row count
"""
struct ArrowFileGroup
    data_path::String
    scores_path::String
    n_rows::Int
end

"""
    ArrowFileGroup(data_path::String)

Construct from a data file path. Derives scores_path by replacing
`_foldN.arrow` with `_foldN_scores.arrow`. Row count read from file.
"""
function ArrowFileGroup(data_path::String)
    scores_path = _derive_scores_path(data_path)
    tbl = Arrow.Table(data_path)
    n = length(tbl[1])
    return ArrowFileGroup(data_path, scores_path, n)
end

function _derive_scores_path(data_path::String)
    m = match(r"(_fold\d+)(\.arrow)$", data_path)
    if m !== nothing
        return data_path[1:m.offset-1] * m.captures[1] * "_scores" * m.captures[2]
    end
    # Fallback: insert _scores before .arrow
    return replace(data_path, r"\.arrow$" => "_scores.arrow")
end

"""
    get_fold_number(group::ArrowFileGroup) -> UInt8

Parse the CV fold number from the file path suffix (e.g. `_fold0.arrow` → 0).
"""
function get_fold_number(group::ArrowFileGroup)
    m = match(r"_fold(\d+)\.arrow$", group.data_path)
    m === nothing && error("Cannot parse fold number from path: $(group.data_path)")
    return parse(UInt8, m.captures[1])
end

#############################################################################
# ArrowFilePSMContainer
#############################################################################

"""
    ArrowFilePSMContainer <: AbstractPSMContainer

File-backed PSM container for out-of-memory percolator scoring.
Data lives on disk as Arrow files grouped by (ms_file_idx, cv_fold).

Only a minimal subset of the AbstractPSMContainer interface is implemented.
Operations that require all data in memory will error with guidance messages.
"""
struct ArrowFilePSMContainer <: AbstractPSMContainer
    file_groups::Vector{ArrowFileGroup}
    total_rows::Int
    max_training_psms::Int
end

"""
    ArrowFilePSMContainer(data_paths::Vector{String}; max_training_psms::Int=typemax(Int))

Construct from a vector of data file paths (e.g. `*_fold0.arrow`, `*_fold1.arrow`).
`max_training_psms` sets the maximum number of PSMs to sample for training.
"""
function ArrowFilePSMContainer(data_paths::Vector{String}; max_training_psms::Int=typemax(Int))
    groups = [ArrowFileGroup(p) for p in data_paths]
    total = sum(g.n_rows for g in groups)
    return ArrowFilePSMContainer(groups, total, max_training_psms)
end

#############################################################################
# Implemented Interface Methods
#############################################################################

nrows(c::ArrowFilePSMContainer) = c.total_rows

function get_cv_folds(c::ArrowFilePSMContainer)
    folds = Set{UInt8}()
    for g in c.file_groups
        push!(folds, get_fold_number(g))
    end
    return sort!(collect(folds))
end

"""
    get_file_groups_for_fold(c::ArrowFilePSMContainer, fold::UInt8) -> Vector{ArrowFileGroup}

Return all file groups belonging to the given CV fold.
"""
function get_file_groups_for_fold(c::ArrowFilePSMContainer, fold::UInt8)
    return [g for g in c.file_groups if get_fold_number(g) == fold]
end

function has_column(c::ArrowFilePSMContainer, col::Symbol)
    isempty(c.file_groups) && return false
    tbl = Arrow.Table(c.file_groups[1].data_path)
    return col in Arrow.names(tbl)
end

#############################################################################
# Error Stubs (data is on disk, not in memory)
#############################################################################

const _ARROW_STUB_MSG = "ArrowFilePSMContainer stores data on disk. " *
    "Use file-group iterators (get_file_groups_for_fold) for OOM processing."

function get_column(::ArrowFilePSMContainer, col::Symbol)
    error("get_column(:$col): $_ARROW_STUB_MSG")
end

function set_column!(::ArrowFilePSMContainer, col::Symbol, data)
    error("set_column!(:$col): $_ARROW_STUB_MSG")
end

function to_dataframe(::ArrowFilePSMContainer)
    error("to_dataframe: $_ARROW_STUB_MSG")
end

function sort_container!(::ArrowFilePSMContainer, cols::Vector{Symbol})
    error("sort_container!($cols): $_ARROW_STUB_MSG")
end

function get_view(::ArrowFilePSMContainer, indices::AbstractVector{<:Integer})
    error("get_view: $_ARROW_STUB_MSG")
end

function copy_container(::ArrowFilePSMContainer)
    error("copy_container: $_ARROW_STUB_MSG")
end

function get_fold_indices(::ArrowFilePSMContainer, fold)
    error("get_fold_indices($fold): $_ARROW_STUB_MSG")
end

function get_train_indices(::ArrowFilePSMContainer, test_fold)
    error("get_train_indices($test_fold): $_ARROW_STUB_MSG")
end
