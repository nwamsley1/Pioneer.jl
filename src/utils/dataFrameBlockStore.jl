# Copyright (C) 2024 Nathan Wamsley
# Licensed under AGPL v3+; see LICENSE.

"""
    DataFrameBlockStore(path, budget)

Retain owned DataFrame blocks up to `budget` bytes, then spill to one seekable
file at `path`. Callers must bound individual blocks and close the store before
removing its scratch directory. Reads may cache one block outside the budget.
This store is single-threaded; appending transfers ownership of the block.
"""
mutable struct DataFrameBlockStore
    path::String
    budget::Int
    bytes::Int
    blocks::Vector{DataFrame}
    offsets::Vector{Int64}
    row_ends::Vector{Int}
    io::Union{Nothing, IOStream}
    cached_index::Int
    cached_block::DataFrame
end

function DataFrameBlockStore(path, budget)
    budget > 0 || throw(ArgumentError("DataFrame block cache budget must be positive"))
    return DataFrameBlockStore(path, budget, 0, DataFrame[], Int64[], Int[],
                              nothing, 0, DataFrame())
end
Base.close(store::DataFrameBlockStore) = store.io === nothing ? nothing : close(store.io)
Base.flush(store::DataFrameBlockStore) = store.io === nothing ? nothing : flush(store.io)

function _store_dataframe_block!(store::DataFrameBlockStore, block::DataFrame)
    bytes = store.io === nothing ? Base.summarysize(block) : 0
    if store.io === nothing && store.bytes + bytes > store.budget
        store.io = open(store.path, "w+")
        for saved in store.blocks
            push!(store.offsets, position(store.io))
            Serialization.serialize(store.io, saved)
        end
        empty!(store.blocks)
        store.bytes = 0
    end
    if store.io === nothing
        push!(store.blocks, block)
        store.bytes += bytes
    else
        seekend(store.io)
        push!(store.offsets, position(store.io))
        Serialization.serialize(store.io, block)
    end
    push!(store.row_ends, (isempty(store.row_ends) ? 0 : last(store.row_ends)) + nrow(block))
    return nothing
end

function _dataframe_block(store::DataFrameBlockStore, index; cache::Bool=false)
    store.io === nothing && return store.blocks[index]
    if cache && store.cached_index == index
        return store.cached_block
    end
    store.cached_index = 0
    store.cached_block = DataFrame()
    seek(store.io, store.offsets[index])
    block = Serialization.deserialize(store.io)::DataFrame
    if cache
        store.cached_block = block
        store.cached_index = index
    end
    return block
end
