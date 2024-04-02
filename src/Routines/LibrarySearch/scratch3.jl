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
    n_bins::I
end
getNode(ttree::TournamentTree, node_id::UInt32) = ttree.nodes[node_id]
getNodes(ttree::TournamentTree) = ttree.nodes
setNode!(ttree::TournamentTree, node::TTreeNode, node_id::UInt32) = ttree.nodes[node_id] = node


function fillNode!(ttree::TournamentTree{I,T}, node_id::UInt32, val::T, bin_id::I, top_id::I, bin_len::I) where {I<:Unsigned, T<:Real}
    setNode!(ttree, 
                TTreeNode(val, bin_id, top_id, bin_len),
                node_id
    )
end

function growTTree!(ttree::TournamentTree{I,T}, N::I) where {I<:Unsigned, T<:Real}
    ttree.nodes = Vector{TTreeNode{I, T}}(undef, N + N - 1)
    ttree.n_bins = N
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
    if N+N-1 > length(getNodes(ttree))
        growTTree!(ttree, N)
    end
    println("N $N")
    #Fill leaves 
    i = one(UInt32)
    for sub_bin_range in sub_bin_ranges
        top_id = first(sub_bin_range)
        last_id = last(sub_bin_range)
        fillNode!(ttree, 
                    i, #node_id
                    values[first(sub_bin_range)], #val
                    i, #bin_id
                    top_id,
                    last_id
                    )
        i += one(UInt32)
    end
    for i in range(n_sub_bins+one(UInt32), N) #Fill placeholder leaves (if last_sub_bin_idx is not a power of 2)
        fillNode!(ttree, 
                i,
                typemax(T),
                i, #leaf_id
                one(I),
                one(I)
        )
    end
    #Fill remaining nodes 
    lower_node_id = 1
    upper_node_id = N + 1
    n = N >> 1
    while n > 1#node_id <= (N - 1) #While there are nodes left to fill
        start = upper_node_id
        n = n >> 1
        while upper_node_id <= start + n + 1 #Fill current layer of nodes 

            node_a, node_b = ttree.nodes[lower_node_id], ttree.nodes[lower_node_id+1]

            if getVal(node_a) <= getVal(node_b)
                ttree.nodes[upper_node_id] = node_a
            else
                ttree.nodes[upper_node_id] = node_b
            end
            lower_node_id += 2
            upper_node_id += 1
        end
    end

end

function removeSmallestElement!(
    ttree::TournamentTree{I, T},
    values_in::Vector{T},
    values_out::Vector{T},
    out_idx::Int64
) where {I<:Unsigned,T<:Real}
    #Number of bins in the merge 
    n = ttree.n_bins
    #Get node containing minimum value 
    node_id = n + n - one(UInt32)
    node = getNode(ttree, node_id)
    #Write minimum value to output array and increase index 
    values_out[out_idx] = getVal(node)
    out_idx += 1
    #Repair tree 
    top_id = getTopID(node)
    last_id = getLastID(node)
    bin_id = getBinID(node)
    if top_id > last_id
        top_id += one(I)
        setNode!(ttree, 
                TTreeNode(
                typemax(T),
                bin_id,
                top_id,
                last_id),
                bin_id)
    else
        top_id += one(I)
        setNode!(ttree, 
                TTreeNode(
                values_in[top_id],
                bin_id,
                top_id,
                last_id),
                bin_id)
    end
    println("node_id $node_id")
    lower_node_id = bin_id - (bin_id%I(2) === zero(I))
    println("lower_node_id $lower_node_id")
    N = ttree.n_bins
    layer_size = n
    n = (lower_node_id >> 2)
    while N < ttree.n_bins + ttree.n_bins - one(I)
        println("n $n ", N + n)
        node_a, node_b = getNode(ttree, lower_node_id), getNode(ttree, lower_node_id + one(I))
        if getVal(node_a) <= getVal(node_b)
            setNode!(ttree, node_a,  N + n + one(I))
        else
            setNode!(ttree, node_b,  N + n + one(I))
        end
        layer_size = layer_size >> 1
        N += layer_size
        n = n >> 1
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