using Tables, Arrow
@time table = Arrow.Table("data/parquet/GAPDH_VGVNGFGR.arrow")

#dict = Dict("customer age" => [15, 20, 25],
#                   "first name" => ["Rohit", "Rahul", "Akshat"])
#DataFrame(dict)
#table.precursorMZ
def getSub(x, low)
    Set(findall(x->x<478, skipmissing(table.precursorMZ)));
b = Set(findall(x->x>477, skipmissing(table.precursorMZ)));

function getSub(mean::Float64, diff::Float64, array::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}})
    findall(x->coalesce(abs(mean-x)<diff, false), array)
end


findall(x->coalesce(abs(477.0-x)<1, false), table.precursorMZ)

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