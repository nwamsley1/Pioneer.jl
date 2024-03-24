function branchless_binary(t::Vector{Float32},
                           x::Float32,
                           y::Float32,
                           lo::UInt32,
                           hi::UInt32)
    #hi_f = hi
    base = lo
    @fastmath begin
        len = hi - lo + UInt32(1)

        while len > 1
            mid = len>>>0x01
            base += (t[base + mid - UInt32(1)] < x)*mid
            len -= mid# - UInt32(1)
        end
        window_start = base

        #Make initial guess closer to lhs
        init_guess = min(window_start + UInt32(2),hi) #Could exceed bounds. 
        tf = t[init_guess]>y  #Is the initial guess greater than the query
        hi = tf*init_guess + (one(UInt32)-tf)*hi
        len = hi - base  + one(UInt32)
        base = hi
        #If guess wass correct, adjust starting point of search and search length. 
        #base = tf*(init_guess - hi) + hi#tf*(hi - lo) + lo#(#What would mid be if the initial guess was true?)
        #len = (UInt32(4))*(tf) + (1-tf)*((hi - window_start + UInt32(1)))
        #println("tf $tf base $base len $len window_start $window_start hi $hi")
        while len > 1
            mid = len>>>0x01
            base -= (t[base - mid + UInt32(1)] > y)*mid
            len -= mid
        end
        window_stop = base
        #println("tf $tf base $base len $len window_start $window_start window_stop $window_stop hi $hi")
        #println("test ", window_start === window_stop)
        if window_start === window_stop
            if t[window_start]>y
                return one(UInt32), zero(UInt32)
            end
            if t[window_stop]<x
                return one(UInt32), zero(UInt32)
            end
        end
    end
    return window_start, window_stop
end


x = 4.1f0
y = 4.5f0
test_t = 100.0f0.*sort(rand(Float32, 100000)); 
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]
test_t[first(answer)-1:last(answer)+1]


x = 4.1f0
y = 4.5f0
test_t = Float32[0.0, 4.1, 5.0, 5.2]
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(3))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]
test_t[first(answer)-1:last(answer)+1]



x = 4.1f0
y = 4.5f0
test_t = Float32[0.0, 4.1, 5.0]
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(1))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]
test_t[first(answer)-1:last(answer)+1]


x = 4.1f0
y = 4.5f0
test_t = 100.0f0.*sort(rand(Float32, 100000))   
test_t = test_t[(test_t .< 4.0) .| (test_t .> 6.0)]
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]


