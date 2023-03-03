using Tables, Arrow
@time table = Arrow.Table("data/parquet/GAPDH_VGVNGFGR.arrow")

#dict = Dict("customer age" => [15, 20, 25],
#                   "first name" => ["Rohit", "Rahul", "Akshat"])
#DataFrame(dict)
#table.precursorMZ
#def getSub(x, low)
#    Set(findall(x->x<478, skipmissing(table.precursorMZ)));
#b = Set(findall(x->x>477, skipmissing(table.precursorMZ)));

function getSub(mean::Float64, ppm::Float64, array::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}})
    findall(x->coalesce(abs(mean-x)<((mean/1000000.0)*ppm), false), array)
end

#Function that given a sorted list of transitions gets all the hits
#of_eltype(Float32, table.intensities[1]))
#function getHits(mass_list::Vector{Float32}, ppm::Float64, masses::MappedArray{Float32, 1, Vector{Union{Missing, Float32}}, MappedArrays.var"#7#9"{Float32}, MappedArrays.var"#8#10"{Union{Missing, Float32}}}, 
#    intensities::MappedArray{Float32, 1, Vector{Union{Missing, Float32}}, MappedArrays.var"#7#9"{Float32}, MappedArrays.var"#8#10"{Union{Missing, Float32}}})
function getHits(mass_list::Vector{Float32}, ppm::Float32, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
    function getMzFeatures(mass_list::Vector{Float32}, ppm::Float32)
       # [(low=x-(x/1000000.0)*ppm, mass = x, high=x+(x/1000000.0)*ppm) for x in mass_list]
        #[(low=150, mass = 151, high=152), (low=893, mass = 894, high=895)]
        (map(x->(low=x-(x/1000000.0)*ppm, mass = x, high=x+(x/1000000.0)*ppm, test = ppm*10.0, testb = ppm*20.0), mass_list))
    end
    MzFeatures = getMzFeatures(mass_list, ppm)
    feature = 1
    peak = 1
    #println("test");
    while (peak <= length(masses)) & (feature <= length(mass_list))
        if masses[peak] > MzFeatures[feature].low
            if masses[peak] < MzFeatures[feature].high
                #"masses[peak] is in the range(low, high), "
                #There could be multiple peaks in the tolerance and we want to 
                #choose the one that is closest in mass
                smallest_diff = abs(masses[peak] - MzFeatures[feature].mass)
                i = 0
                @inbounds while (masses[peak+1+i] < MzFeatures[feature].high)
                    new_diff = abs(masses[peak+1+i] - MzFeatures[feature].mass)
                    if new_diff < smallest_diff
                        smallest_diff = new_diff
                        peak = peak+1+i
                        i = 0
                    end
                    i+=1
                end
                #println(MzFeatures[feature].mass, " has mass ", masses[peak], " and intensity ", intensities[peak])
                feature += 1
                continue
            end
            feature += 1
            continue
        end
        peak+=1
        #println(peak)
        #println(feature)
    end
end


test_masses =  Vector{Union{Missing, Float32}}(
    [151.67221f0, missing, 894.0937f0, 894.0938f0, 894.0939f0]
    )
test_intensities =  Vector{Union{Missing, Float32}}([missing for x in test_masses])
test_ppm = Float32(20.0)
test_mz = Vector{Float32}([151.67221f0, 700.0, 894.0938f0])

getHits(test_mz, test_ppm, test_masses, test_intensities)
#import Base: <
#function <(a::Any, b::Missing)
#    false
#end
#<(a::Missing, b::Any) = <(b, a)

#import Base: -
#function -(a::Float64, b::Missing)
#    a
#end
#-(a::Missing, b::Float64) = -(b, a)