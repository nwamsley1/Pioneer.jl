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

function getHits(mass_list::Vector{Float32}, ppm::Float64, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
    MzFeatures = map(x->(low=x-(x/1000000.0)*ppm, mass = x, high=x+(x/1000000.0)*ppm), mass_list);
    feature = 1
    peak = 1
    #println("test");
    while (peak<=length(masses)) & (feature <=length(mass_list))
        #println(masses[peak])
        if masses[peak] > MzFeatures[feature].low
            if masses[peak] < MzFeatures[feature].high
                println(masses[peak])
                println(MzFeatures[feature])
                #"masses[peak] is in the range(low, high), "
                #There could be multiple peaks in the tolerance and we want to 
                #choose the one that is closest in mass
                diff = abs(masses[peak] - MzFeatures[feature].mass)
                i = 0
                @inbounds while (masses[peak+1+i] < MzFeatures[feature].high)
                    new_diff = abs(masses[peak+1+i] - MzFeatures[feature].mass) 
                    if new_diff < diff
                        diff = new_diff
                        peak = peak+1+i
                        i = 0
                    end
                    i+=1
                end
                println(MzFeatures[feature].mass, " has mass ", masses[peak], " and intensity ", intensities[peak])
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
    println(feature)
    println(peak)
end
getHits([151.67221f0, 894.0938f0], 20.0, [151.67221f0, missing, 894.0937f0, 894.0938f0, 894.0939f0], 
                                        [Float32(100.0), missing, Float32(100.0), missing ,missing])
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