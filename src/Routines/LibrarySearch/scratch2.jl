test_values = rand(Float32, 100)
sorted_out = zeros(Float32, 100)
sorted_sub_ranges = Vector{UnitRange{UInt32}}([1:10, 11:15, 12:12, 95:100, 85:87])
for sub_range in sorted_sub_ranges
    sort!(@view(test_values[sub_range]))
end
test_ttree = TournamentTree(
    Vector{TTreeNode{UInt32, Float32}}(undef, 1),
    Vector{Leaf{UInt32, Float32}}(undef, 2)
)

buildTTree!(
    test_ttree,
    sorted_sub_ranges,
    UInt32(length( sorted_sub_ranges)),
    test_values
)
