module PDFGenerator
import ..Plots

# Parse a PNG produced by Plots/GR (color type 2 / RGB, 8-bit, non-interlaced)
# and return (width, height, idat_bytes). The IDAT zlib stream is exactly what
# PDF's /FlateDecode + /Predictor 15 (PNG predictor) expects, so we can embed
# it directly with no inflate/deflate round-trip.
function _parse_png_rgb(path::String)
    bytes = read(path)
    @assert bytes[1:8] == UInt8[0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A] "not a PNG"
    width = 0
    height = 0
    color_type = UInt8(0)
    bit_depth = UInt8(0)
    interlace = UInt8(0)
    idat = IOBuffer()
    i = 9
    while i <= length(bytes)
        len = (UInt32(bytes[i])   << 24) |
              (UInt32(bytes[i+1]) << 16) |
              (UInt32(bytes[i+2]) << 8)  |
              (UInt32(bytes[i+3]))
        ctype = String(bytes[i+4:i+7])
        data_start = i + 8
        data_end   = i + 7 + Int(len)
        if ctype == "IHDR"
            width      = (UInt32(bytes[data_start])   << 24) |
                         (UInt32(bytes[data_start+1]) << 16) |
                         (UInt32(bytes[data_start+2]) << 8)  |
                         (UInt32(bytes[data_start+3]))
            height     = (UInt32(bytes[data_start+4]) << 24) |
                         (UInt32(bytes[data_start+5]) << 16) |
                         (UInt32(bytes[data_start+6]) << 8)  |
                         (UInt32(bytes[data_start+7]))
            bit_depth  = bytes[data_start+8]
            color_type = bytes[data_start+9]
            interlace  = bytes[data_start+12]
        elseif ctype == "IDAT"
            write(idat, view(bytes, data_start:data_end))
        elseif ctype == "IEND"
            break
        end
        i = data_end + 5   # length + type(4) + data + CRC(4)
    end
    @assert color_type == 2 "PNG color type $color_type unsupported; expected 2 (RGB)"
    @assert bit_depth == 8 "PNG bit depth $bit_depth unsupported; expected 8"
    @assert interlace == 0 "PNG interlace $interlace unsupported; expected 0"
    return Int(width), Int(height), take!(idat)
end

# Minimal multipage PDF builder. Embeds one PNG per page using /FlateDecode
# with /Predictor 15 (PNG row predictor) — PDF readers decode this natively.
# Output page size = image size in PDF points (1 pt per pixel, ~96 dpi).
function write_pngs_to_pdf(png_paths::AbstractVector{String}, dest::String)
    io = IOBuffer()
    write(io, "%PDF-1.4\n%\xE2\xE3\xCF\xD3\n")  # binary marker
    offsets = Int[]
    function emit(obj_str::AbstractString)
        push!(offsets, position(io))
        write(io, obj_str)
    end
    function emit_stream(header::AbstractString, payload::Vector{UInt8})
        push!(offsets, position(io))
        write(io, header)
        write(io, "stream\n")
        write(io, payload)
        write(io, "\nendstream\nendobj\n")
    end

    # Object 1: Catalog
    # Object 2: Pages
    # Then for each page: Page (3, 6, 9,...), Image XObject (4, 7, 10,...), Contents (5, 8, 11,...)
    n_pages = length(png_paths)
    page_obj_nums    = [3 + 3*(i-1) for i in 1:n_pages]
    image_obj_nums   = [4 + 3*(i-1) for i in 1:n_pages]
    content_obj_nums = [5 + 3*(i-1) for i in 1:n_pages]

    emit("1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n")
    kids = join(["$(n) 0 R" for n in page_obj_nums], " ")
    emit("2 0 obj\n<< /Type /Pages /Kids [$kids] /Count $n_pages >>\nendobj\n")

    for (i, png_path) in enumerate(png_paths)
        w, h, idat = _parse_png_rgb(png_path)
        pn = page_obj_nums[i]
        im = image_obj_nums[i]
        cn = content_obj_nums[i]
        # Page
        emit("$pn 0 obj\n<< /Type /Page /Parent 2 0 R " *
             "/MediaBox [0 0 $w $h] " *
             "/Resources << /XObject << /Im0 $im 0 R >> >> " *
             "/Contents $cn 0 R >>\nendobj\n")
        # Image XObject
        img_header = "$im 0 obj\n<< /Type /XObject /Subtype /Image " *
                     "/Width $w /Height $h /ColorSpace /DeviceRGB " *
                     "/BitsPerComponent 8 /Filter /FlateDecode " *
                     "/DecodeParms << /Predictor 15 /Columns $w /Colors 3 /BitsPerComponent 8 >> " *
                     "/Length $(length(idat)) >>\n"
        emit_stream(img_header, idat)
        # Page Contents: paint image scaled to W×H pt
        content = "q\n$w 0 0 $h 0 0 cm\n/Im0 Do\nQ\n"
        content_bytes = Vector{UInt8}(content)
        content_header = "$cn 0 obj\n<< /Length $(length(content_bytes)) >>\n"
        emit_stream(content_header, content_bytes)
    end

    # xref
    xref_offset = position(io)
    n_objs = 2 + 3 * n_pages
    write(io, "xref\n0 $(n_objs+1)\n0000000000 65535 f \n")
    for off in offsets
        write(io, lpad(string(off), 10, '0'), " 00000 n \n")
    end
    write(io, "trailer\n<< /Size $(n_objs+1) /Root 1 0 R >>\nstartxref\n$xref_offset\n%%EOF\n")

    open(dest, "w") do f
        write(f, take!(io))
    end
    return dest
end
end  # module PDFGenerator

# Public helper: render each Plots.Plot to a temp PNG, then assemble into a
# multipage PDF. The temp PNGs are deleted on success. Output is one PDF per
# call, mirroring the historical save_multipage_pdf UX, but the pages are
# embedded raster (~30-150 KB each) instead of vector (~hundreds of KB-MB).
function save_multipage_pdf(plots::AbstractVector{<:Plots.Plot}, dest::String)
    ensure_directory_exists(dest)
    isempty(plots) && return dest
    mktempdir() do tmp
        png_paths = String[]
        for (i, p) in enumerate(plots)
            png = joinpath(tmp, "page_$i.png")
            Plots.savefig(p, png)
            push!(png_paths, png)
        end
        PDFGenerator.write_pngs_to_pdf(png_paths, dest)
    end
    return dest
end
