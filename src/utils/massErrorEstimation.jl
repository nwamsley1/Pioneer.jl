
struct MassErrorModel{T<:AbstractFloat}
    mass_offset::T#UniformSpline{N, T}
    mass_tolerance::Tuple{T, T}#UniformSpline{N, T}
end

getRightTol(mem::MassErrorModel) = last(mem.mass_tolerance)
getLeftTol(mem::MassErrorModel) = first(mem.mass_tolerance)
getMassOffset(mem::MassErrorModel) = last(mem.mass_offset)

function getMassCorrection(mem::MassErrorModel)
    return mem.mass_offset
end

function getLocation(mem::MassErrorModel)
    return mem.location
end

#Adjust a theoretical mass to be
function (mem::MassErrorModel)(mass::Float32)
    ppm_norm = Float32(1e6)
    ppm = mass/ppm_norm
    mass += getMassOffset(mem)*ppm
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - l_tol), Float32(mass + r_tol)
end

#Correct an empeirical mass. 
function getCorrectedMz(mem::MassErrorModel, mz::Float32)
    return  Float32(mz - getMassOffset(mem)*(mz/1e6))
end

"""
   getMzBounds(mem::MassErrorModel, mass::Float32)

Given a theoretical mass and a `MassErrorModel`, gets the minimum and maximum bounds for an expected empirical mass.
Critically, assumes the empirical mass has already been corrected using mass offset. See `getCorrectedMz`. 

### Input

- `mem::MassErrorModel`: -- Model for the mass error of an ion 
- `mass::Float32` -- A theoretical mass/mz to get boundaries for 
### Output
Tuple{Float32, Float32}
A tuple with the lower and upper boundary respectively. 
### Notes

- Suppose the mass error was 3 ppm. And the tolerance was 10 ppm and 5ppm on the left and right hand sides respectively. 
A theoretical mass of 1000000.0f0 m/z, would have a tolerance of (999990.0f0, 1000005.0f0). A theoretical mass falling 
### Algorithm 

### Examples 

"""
function getMzBounds(mem::MassErrorModel, mass::Float32)
    ppm = mass/(1e6)
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - l_tol), Float32(mass + r_tol)
end

function ModelMassErrs(ppm_errs::Vector{Float32};
                       frag_err_quantile::Float32 = 0.01f0,
                       out_fdir::String = "./",
                       out_fname::String = "mass_err_estimate")

    bins = LinRange(minimum(ppm_errs), maximum(ppm_errs), 100)
    mass_err = median(ppm_errs)
    ppm_errs = ppm_errs .- mass_err
    l_bound = quantile(ppm_errs, frag_err_quantile)
    r_bound = quantile(ppm_errs, 1 - frag_err_quantile)
    errs = ppm_errs#ppm_errs .+ mass_err
    plot_title = ""
    n = 0
    for i in range(1, length(out_fname))
        n += 1
        if n > 24
            n = 1
            plot_title *= "\n"
        end
        plot_title *= out_fname[i]
    end
    n = length(errs)
    p = Plots.histogram(errs,
                    orientation = :h, 
                    yflip = true,
                    #seriestype=:scatter,
                    title = plot_title*"\n n = $n",
                    xlabel = "Count",
                    ylabel = "Mass Error (ppm)",
                    label = nothing,
                    bins = bins,
                    ylim = (minimum(ppm_errs)-2, maximum(ppm_errs)+2),
                    topmargin =15mm,
                    #bottommargin = 10mm,
                    )

    Plots.hline!([l_bound + mass_err, mass_err, r_bound + mass_err], label = nothing, color = :black, lw = 2)
    l_err = l_bound + mass_err
    Plots.annotate!(last(xlims(p)), l_err, text("$l_err", :black, :right, :bottom, 12))
    r_err = r_bound + mass_err
    Plots.annotate!(last(xlims(p)), r_err, text("$r_err", :black, :right, :bottom, 12))
    Plots.annotate!(last(xlims(p)), mass_err, text("$mass_err", :black, :right, :bottom, 12))
    savefig(p, joinpath(out_fdir, out_fname)*".pdf")

    MassErrorModel(
                    Float32(mass_err),
                    (Float32(abs(l_bound)), Float32(abs(r_bound)))
                    )
end

