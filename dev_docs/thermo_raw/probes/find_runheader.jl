# Locate the RunHeader by its self-pointer (UInt64 at RunHeader+7472 == RunHeader address), then find where else
# that address is stored (the FileHeader / RawFileInfo pointer).
using Mmap
function find_rh(path)
    b = Mmap.mmap(path); n = length(b)
    rd(o) = unsafe_load(Ptr{UInt64}(pointer(b, o + 1)))
    hits = Int[]
    GC.@preserve b for o in 7472:n-8
        v = rd(o)
        v == o - 7472 && push!(hits, o - 7472)
    end
    ptrs = Dict{Int, Vector{Int}}()
    for h in hits
        pat = collect(reinterpret(UInt8, [UInt64(h)])); i = 1; locs = Int[]
        while length(locs) < 6
            r = findnext(pat, b, i); r === nothing && break
            push!(locs, first(r) - 1); i = first(r) + 1
        end
        ptrs[h] = locs
    end
    hits, ptrs
end
for path in ARGS
    @time hits, ptrs = find_rh(path)
    println(basename(path), ": RunHeader candidates ", hits, "; stored at ", ptrs)
end
