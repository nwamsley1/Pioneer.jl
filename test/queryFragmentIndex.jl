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

i = 0
while i < 10
    println(i)
    i += 1
end

transitions = getTransitions(Precursor("IEKLFSQQAQIEAQVLDK"), UInt8[1], UInt8[0], y_start = 3, b_start = 3, ppm = 20.0)
i = 13213
_intensities = MS_TABLE[:intensities][i]
_masses = MS_TABLE[:masses][i]
_masses =  Vector{Union{Missing, Float32}}([x for x in [175.11836f0 181.62372f0 185.1641f0 197.1272f0 213.1566f0 215.14012f0 242.14981f0 243.13153f0 249.1609f0 314.15317f0 356.2205f0 371.20816f0 406.19803f0 415.21606f0 431.2224f0 471.22815f0 488.25815f0 489.25125f0 559.296f0 560.2905f0 647.34467f0 657.3512f0 672.37476f0 673.37476f0 783.4229f0 801.4133f0 802.43634f0 914.503f0 915.4875f0 1027.5999f0 1061.5077f0 1062.5331f0 1156.6205f0 1328.6967f0]])
_intensities = _masses
matchPeaks(transitions, 
_intensities, 
_masses, 
#δs = params[:δs],
δs = zeros(Float32, (1,)),
scan_idx = UInt32(i),
ms_file_idx = UInt32(1),
min_intensity = Float32(0.0))