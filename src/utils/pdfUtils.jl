# Utility functions for working with PDFs

"""
    merge_pdfs_safe(files::Vector{String}, dest::String; cleanup::Bool=false)

Merge `files` into a single PDF at `dest` using the `pdfunite` tool. All
temporary files are created inside `dest`'s directory so write
permissions are respected on read-only installations.

If `cleanup` is true, the source files will be removed after merging.
"""
function merge_pdfs_safe(files::Vector{String}, dest::String; cleanup::Bool=false)
    ensure_directory_exists(dest)
    pdfunite = something(Sys.which("pdfunite"), "pdfunite")
    tmp_path, io = mktemp(dir=dirname(dest))
    close(io)
    try
        run(`$pdfunite $(files...) $tmp_path`)
        mv(tmp_path, dest; force=true)
    finally
        isfile(tmp_path) && rm(tmp_path, force=true)
    end
    if cleanup
        for f in files
            safeRm(f, nothing)
        end
    end
    return dest
end

