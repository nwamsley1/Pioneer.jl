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
        len = hi - base + UInt32(1)
        base = hi
        mid = len>>>0x01
        init_guess = window_start + UInt32(1)

        tf = t[init_guess]>y
        base -= tf*mid
        len = hi - tf*init_guess + one(UInt32)#hi - base + UInt32(1)
        #println("lo $lo hi $hi init_guess $init_guess mid $mid")
        while len > 1
            mid = len>>>0x01
            base -= (t[base - mid + UInt32(1)] > y)*mid
            len -= mid
            println("base $base len $len")
        end
        window_stop = base

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
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(4))
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


