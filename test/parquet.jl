
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