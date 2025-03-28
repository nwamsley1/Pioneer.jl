abstract type ChromObject end
struct MS2ChromObject <: ChromObject
    rt::Float32
    intensity::Float32
    scan_idx::UInt32
    precursor_idx::UInt32
end
struct MS1ChromObject <: ChromObject
    rt::Float32
    intensity::Float32
    m0::Bool
    n_iso::UInt8
    scan_idx::UInt32
    precursor_idx::UInt32
end

function growChromObjects!(chromatograms::Vector{ChromObject}, block_size::Int64)
    chromatograms = append!(chromatograms, Vector{ChromObject}(undef, block_size))
end

