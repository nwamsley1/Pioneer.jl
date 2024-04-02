abstract type Node end
getVal(nd::Node) = nd.val
getBinID(nd::Node) = nd.bin_id
getTopID(nd::Node) = nd.top_id
getLastID(nd::Node) = nd.last_id

struct TTreeNode{I<:Unsigned,T<:Real} <: Node
    val::T
    bin_id::I
    top_id::I
    last_id::I
end

mutable struct TournamentTree{I<:Unsigned,T<:Real}
    nodes::Vector{TTreeNode{I,T}}
end
getNode(ttree::TournamentTree, node_id::Int64) = ttree.nodes[node_id]
setNode!(ttree::TournamentTree, node::Node, node_id::Int64) = ttree.nodes[node_id] = node


function fillNode!(ttree::TournamentTree{I,T}, node_id::Int64, bin_id::I, top_id::I, bin_len::I) where {I<:Unsigned, T<:Real}
    setNode!(ttree, 
                TTreeNode(val, bin_id, top_id, bin_len),
                node_id
    )
end

function growTTree!(ttree::TournamentTree{I,T}, N::Integer) where {I<:Unsigned, T<:Real}
    ttree.nodes = Vector{TTreeNode{I, T}}(undef, N + N - 1)
end

function getNearestPowerOf2(n::UInt32)
    n -= 1
    n |= n >> 1
    n |= n >> 2
    n |= n >> 4
    n |= n >> 8
    n |= n >> 16
    n += 1
    return n
end


function buildTTree!(
    ttree::TournamentTree{I, T},
    sub_bin_ranges::Vector{UnitRange{I}},
    n_sub_bins::UInt32,
    values::Vector{T}
    ) where {I<:Unsigned,T<:Real}

    #Number of leaves needs to be the next highest power of 2. 
    N = UInt32(getNearestPowerOf2(n_sub_bins))
    if N > length(getLeaves(ttree))
        growTTree!(ttree, N)
    end

    #Fill leaves 
    i = one(UInt32)
    for sub_bin_range in sub_bin_ranges
        top_id = first(sub_bin_range)
        last_id = last(sub_bin_range)
        fillNode!(ttree, 
                    values[first(sub_bin_range)], #val
                    i, #bin_id
                    top_id,
                    last_id
                    )
        i += one(UInt32)
    end
    for i in range(n_sub_bins, N) #Fill placeholder leaves (if last_sub_bin_idx is not a power of 2)
        fillNode!(ttree, 
                typemax(T),
                i, #leaf_id
                one(I),
                one(I)
        )
    end
    #Fill remaining nodes 
    n = N >> 1 #Number of nodes in base layer of the tree 
    start = node_id
    lower_node_id = 1
    upper_node_id = n + 1
    n = n >> 1
    while n > 1#node_id <= (N - 1) #While there are nodes left to fill
        start = upper_node_id
        while upper_node_id <= start + n #Fill current layer of nodes 

            node_a, node_b = ttree.nodes[lower_node_id], ttree.nodes[lower_node_id+1]

            if getVal(node_a) <= getVal(node_b)
                ttree.nodes[upper_node_id] = node_a
            else
                ttree.nodes[upper_node_id] = node_b
            end
            lower_node_id += 2
            upper_node_id += 1
        end
        n = n >> 1
    end

end

function removeSmallestElement(
    ttree::TournamentTree{I, T},
    values_in::Vector{T},
    values_out::Vector{T},
    out_idx::Int64,
    leading_node_id::Int64
) where {I<:Unsigned,T<:Real}
    node = getNode(ttree, leading_node_id)
    values_out[out_idx] = getVal(node)
    out_idx += 2
    leaf_id = getLeafID(node)
    leaf = getLeaf(ttree, leaf_id)
    if leaf.top_id === leaf.last_id
        setLeaf!(ttree, leaf, leaf_id) = Leaf(
            typemax(T),
            leaf_id,
            leaf_id
        )
    else
        setLeaf!(ttree, leaf, leaf_id) = Leaf(
            values_in[leaf.top_id],
            leaf.top_id + one(I),
            leaf.last_id
        )
    end
end
#=
    #Fill first layer of nodes 
    leaf_id = 1
    node_id = 1
    while leaf_id <= N
        #Compare leaves and assign winner to the node 
        val_a, val_b = getVal(getLeaf(ttree,leaf_id)), getVal(getLeaf(ttree,leaf_id+1))
        if val_a <= val_b
            ttree.nodes[node_id] = TTreeNode(val_a, UInt32(leaf_id))
        else
            ttree.nodes[node_id] = TTreeNode(val_b, UInt32(leaf_id + 1))
        end
        #Increment leaf and node counters 
        leaf_id += 2
        node_id += 1
    end
=#