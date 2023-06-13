lowbound = [10.0, 10.005, 10.015, 11.0, 16.0]
highbound = [10.0, 10.01, 10.015, 15.0, 16.1]
test_fragbins = Vector{FragBin{Float64}}()

for i in 1:length(lowbound)
    push!(test_fragbins, FragBin(lowbound[i], highbound[i], UInt32(i)))
end

@test findFirstFragmentBin(test_fragbins, 10.0, 10.0) == 1
@test findFirstFragmentBin(test_fragbins, 10.001, 10.0) === nothing
@test findFirstFragmentBin(test_fragbins, 10.001, 100.0) == 2

test_precs = Dictionary{UInt32, IonIndexMatch{Float32}}()