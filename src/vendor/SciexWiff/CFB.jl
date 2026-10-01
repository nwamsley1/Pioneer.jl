"""
    CFB

Minimal read-only reader for Microsoft Compound File Binary (CFB/OLE2) containers,
written from the public MS-CFB specification. `.wiff` files are CFB containers.

    cf = CFB.CompoundFile(path)
    CFB.stream_paths(cf)                       # every stream, as "A/B/C" paths
    CFB.read_stream(cf, "SampleSubtree/Sample1/Idx")
"""
module CFB

export CompoundFile, stream_paths, read_stream, has_stream

const MAGIC = (0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1)

# Special sector numbers (MS-CFB §2.1)
const MAXREGSECT = 0xfffffffa
const DIFSECT    = 0xfffffffc
const FATSECT    = 0xfffffffd
const ENDOFCHAIN = 0xfffffffe
const FREESECT   = 0xffffffff
const NOSTREAM   = 0xffffffff

# Directory entry object types (MS-CFB §2.6.1)
const TYPE_UNALLOCATED = 0x00
const TYPE_STORAGE     = 0x01
const TYPE_STREAM      = 0x02
const TYPE_ROOT        = 0x05

struct CFBError <: Exception
    msg::String
end
Base.showerror(io::IO, e::CFBError) = print(io, "CFBError: ", e.msg)

struct DirEntry
    name::String
    type::UInt8
    left::UInt32
    right::UInt32
    child::UInt32
    start_sector::UInt32
    size::UInt64
end

struct CompoundFile
    data::Vector{UInt8}
    sector_size::Int
    mini_sector_size::Int
    mini_cutoff::Int
    fat::Vector{UInt32}
    minifat::Vector{UInt32}
    entries::Vector{DirEntry}
    ministream::Vector{UInt8}
    paths::Dict{String, Int}   # "A/B/C" => entry index (1-based)
end

@inline u16(d, o) = ltoh(reinterpret(UInt16, view(d, o+1:o+2))[1])
@inline u32(d, o) = ltoh(reinterpret(UInt32, view(d, o+1:o+4))[1])
@inline u64(d, o) = ltoh(reinterpret(UInt64, view(d, o+1:o+8))[1])

CompoundFile(path::AbstractString) = CompoundFile(read(path))

function CompoundFile(data::Vector{UInt8})
    length(data) >= 512 || throw(CFBError("file shorter than a CFB header"))
    Tuple(data[1:8]) == MAGIC || throw(CFBError("bad CFB signature"))

    major = u16(data, 0x1a)
    u16(data, 0x1c) == 0xfffe || throw(CFBError("bad byte-order mark"))
    sector_shift = u16(data, 0x1e)
    (major == 3 && sector_shift == 9) || (major == 4 && sector_shift == 12) ||
        throw(CFBError("unsupported version $major / sector shift $sector_shift"))
    sector_size = 1 << sector_shift
    mini_sector_size = 1 << u16(data, 0x20)
    n_fat = Int(u32(data, 0x2c))
    first_dir = u32(data, 0x30)
    mini_cutoff = Int(u32(data, 0x38))
    first_minifat = u32(data, 0x3c)
    first_difat = u32(data, 0x44)
    n_difat = Int(u32(data, 0x48))

    n_sectors = (length(data) - sector_size) ÷ sector_size
    sector_offset(s) = (Int(s) + 1) * sector_size   # 0-based byte offset of sector s
    function sector(s)
        s < n_sectors || throw(CFBError("sector $s out of range ($n_sectors sectors)"))
        o = sector_offset(s)
        view(data, o+1:o+sector_size)
    end

    # FAT sector list: 109 entries in the header, then the DIFAT chain.
    fat_sectors = UInt32[]
    for k in 0:108
        s = u32(data, 0x4c + 4k)
        s <= MAXREGSECT && push!(fat_sectors, s)
    end
    per_sector = sector_size ÷ 4
    s = first_difat
    for _ in 1:n_difat
        s <= MAXREGSECT || break
        blk = sector(s)
        for k in 0:per_sector-2
            f = u32(blk, 4k)
            f <= MAXREGSECT && push!(fat_sectors, f)
        end
        s = u32(blk, 4 * (per_sector - 1))
    end
    length(fat_sectors) >= n_fat ||
        throw(CFBError("header says $n_fat FAT sectors, found $(length(fat_sectors))"))

    fat = UInt32[]
    sizehint!(fat, n_fat * per_sector)
    for fs in fat_sectors[1:n_fat]
        blk = sector(fs)
        for k in 0:per_sector-1
            push!(fat, u32(blk, 4k))
        end
    end

    chain(start) = follow_chain(fat, start, length(fat))
    function read_chain(start)
        out = UInt8[]
        for c in chain(start)
            append!(out, sector(c))
        end
        out
    end

    # Directory
    dirbytes = read_chain(first_dir)
    entries = DirEntry[]
    for o in 0:128:length(dirbytes)-128
        namelen = Int(u16(dirbytes, o + 64))
        nchars = max(namelen ÷ 2 - 1, 0)          # length includes the UTF-16 NUL
        name = transcode(String, [u16(dirbytes, o + 2k) for k in 0:nchars-1])
        size = major == 3 ? UInt64(u32(dirbytes, o + 120)) : u64(dirbytes, o + 120)
        push!(entries, DirEntry(name, dirbytes[o+67], u32(dirbytes, o + 68),
            u32(dirbytes, o + 72), u32(dirbytes, o + 76), u32(dirbytes, o + 116), size))
    end
    isempty(entries) && throw(CFBError("empty directory"))
    entries[1].type == TYPE_ROOT || throw(CFBError("first directory entry is not the root"))

    # Mini FAT and mini stream (the root entry's chain holds the mini stream).
    minifat = UInt32[]
    if first_minifat <= MAXREGSECT
        mf = read_chain(first_minifat)
        minifat = [u32(mf, 4k) for k in 0:length(mf)÷4-1]
    end
    root = entries[1]
    ministream = root.start_sector <= MAXREGSECT ? read_chain(root.start_sector) : UInt8[]
    length(ministream) >= root.size ||
        throw(CFBError("mini stream shorter than root entry size"))

    paths = Dict{String, Int}()
    collect_paths!(paths, entries, root.child, "")

    CompoundFile(data, sector_size, mini_sector_size, mini_cutoff, fat, minifat,
        entries, ministream, paths)
end

"""Follow a FAT (or mini-FAT) chain from `start`, with a cycle guard."""
function follow_chain(table::Vector{UInt32}, start::UInt32, limit::Int)
    out = UInt32[]
    s = start
    while s != ENDOFCHAIN
        s <= MAXREGSECT || throw(CFBError("bad sector $(repr(s)) in chain"))
        s < length(table) || throw(CFBError("chain sector $s beyond table"))
        push!(out, s)
        length(out) <= limit || throw(CFBError("cycle in sector chain"))
        s = table[s+1]
    end
    out
end

"""Walk the red-black sibling tree under `id`, recording streams and storages by path."""
function collect_paths!(paths, entries, id::UInt32, prefix::String, depth::Int = 0)
    depth > length(entries) && throw(CFBError("cycle in directory tree"))
    id == NOSTREAM && return
    id < length(entries) || throw(CFBError("directory id $id out of range"))
    e = entries[id+1]
    collect_paths!(paths, entries, e.left, prefix, depth + 1)
    path = isempty(prefix) ? e.name : prefix * "/" * e.name
    if e.type in (TYPE_STREAM, TYPE_STORAGE)
        paths[path] = id + 1
        e.type == TYPE_STORAGE && collect_paths!(paths, entries, e.child, path, depth + 1)
    end
    collect_paths!(paths, entries, e.right, prefix, depth + 1)
end

normpath_cfb(p::AbstractString) = strip(replace(p, '\\' => '/'), '/')

"""All stream paths (storages excluded), sorted."""
stream_paths(cf::CompoundFile) =
    sort!([p for (p, i) in cf.paths if cf.entries[i].type == TYPE_STREAM])

has_stream(cf::CompoundFile, path::AbstractString) =
    (i = get(cf.paths, normpath_cfb(path), 0); i > 0 && cf.entries[i].type == TYPE_STREAM)

"""Size in bytes of the stream at `path`."""
function stream_size(cf::CompoundFile, path::AbstractString)
    i = get(cf.paths, normpath_cfb(path), 0)
    i > 0 || throw(KeyError(path))
    Int(cf.entries[i].size)
end

"""Read the full contents of the stream at `path` (e.g. `"SampleSubtree/Sample1/Idx"`)."""
function read_stream(cf::CompoundFile, path::AbstractString)
    i = get(cf.paths, normpath_cfb(path), 0)
    i > 0 || throw(KeyError(path))
    e = cf.entries[i]
    e.type == TYPE_STREAM || throw(CFBError("$path is not a stream"))
    size = Int(e.size)
    size == 0 && return UInt8[]
    out = Vector{UInt8}(undef, size)
    if size < cf.mini_cutoff
        ss = cf.mini_sector_size
        pos = 0
        for s in follow_chain(cf.minifat, e.start_sector, length(cf.minifat))
            n = min(ss, size - pos)
            n <= 0 && break
            o = Int(s) * ss
            o + n <= length(cf.ministream) || throw(CFBError("mini sector $s beyond mini stream"))
            copyto!(out, pos + 1, cf.ministream, o + 1, n)
            pos += n
        end
    else
        ss = cf.sector_size
        pos = 0
        for s in follow_chain(cf.fat, e.start_sector, length(cf.fat))
            n = min(ss, size - pos)
            n <= 0 && break
            o = (Int(s) + 1) * ss
            o + n <= length(cf.data) || throw(CFBError("sector $s beyond end of file"))
            copyto!(out, pos + 1, cf.data, o + 1, n)
            pos += n
        end
    end
    pos == size || throw(CFBError("stream $path truncated: $pos of $size bytes"))
    out
end

end # module
