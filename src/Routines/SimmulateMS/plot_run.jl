"""
    plot_ms_run_overview(run::SimulatedMSRun{T}) where {T<:AbstractFloat}

Create an overview plot of a simulated MS run showing MS1 and MS2 scans.

Returns:
- A Plots.jl plot object
"""
function plot_ms_run_overview(run::SimulatedMSRun{T}) where {T<:AbstractFloat}
    # Extract MS1 and MS2 data
    ms1_indices = findall(i -> run.spectra[i].ms_order == 0x01, 1:length(run.spectra))
    ms2_indices = findall(i -> run.spectra[i].ms_order == 0x02, 1:length(run.spectra))
    
    ms1_rt = run.rt_points[ms1_indices]
    ms2_rt = run.rt_points[ms2_indices]
    ms2_mz = [run.spectra[i].centerMz for i in ms2_indices]
    
    # Create the plot
    p = plot(
        xlabel = "Retention Time (min)",
        ylabel = "m/z",
        title = "MS Run Overview",
        legend = :topleft
    )
    
    # Plot MS1 scans
    scatter!(
        p,
        ms1_rt,
        ones(length(ms1_rt)) * 0,  # Place at bottom of plot
        marker = :square,
        color = :blue,
        markersize = 6,
        label = "MS1 Scans"
    )
    
    # Plot MS2 scans
    scatter!(
        p,
        ms2_rt,
        ms2_mz,
        marker = :circle,
        color = :red,
        markersize = 4,
        label = "MS2 Scans"
    )
    
    return p
end