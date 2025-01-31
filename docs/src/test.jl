function create_mass_spectrum_svg(
    masses1::Vector{Float64}, intensities1::Vector{Float64},  
    masses2::Vector{Float64}, intensities2::Vector{Float64},  
    output_file::String;
    x_min::Float64=400.0, x_max::Float64=1400.0,  # explicit x boundaries
    y_max::Float64=71000.0,
    orange_peaks=true)                       # explicit y boundary
    
    svg_header = """<?xml version="1.0" encoding="utf-8"?>
    <svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px"
        viewBox="0 0 800 400" style="enable-background:new 0 0 800 400;" xml:space="preserve">
    <style type="text/css">
        .st0{fill:none;stroke:#DDDDDD;}
        .st1{fill:#444444;}
        .st2{fill:#EE8866;}
        .st3{font-family:'Arial-BoldMT';}
        .st4{font-size:16px;}
        .st5{fill:none;stroke:#000000;stroke-width:2;}
        .st6{fill:none;stroke:#000000;stroke-width:4;}
        .st7{font-family:'ArialMT';}
        .st8{font-size:24px;}
    </style>"""

    grid_lines = """
    <!-- Grid lines -->
    <line class="st0" x1="50" y1="300" x2="750" y2="300"/>
    <line class="st0" x1="50" y1="250" x2="750" y2="250"/>
    <line class="st0" x1="50" y1="200" x2="750" y2="200"/>
    <line class="st0" x1="50" y1="150" x2="750" y2="150"/>
    <line class="st0" x1="50" y1="100" x2="750" y2="100"/>
    <line class="st0" x1="50" y1="50" x2="750" y2="50"/>"""

    axes_and_title = """
    <!-- Title -->
    <text transform="matrix(1 0 0 1 287.9766 30)" class="st3 st4">Mass Spectrum</text>

    <!-- Axes -->
    <line class="st5" x1="50" y1="350" x2="750" y2="350"/>
    <line class="st6" x1="50" y1="50" x2="50" y2="350"/>

    <!-- X-axis label -->
    <text transform="matrix(1 0 0 1 388.7241 394.0034)" class="st7 st8">m/z</text>

    <!-- X-axis labels -->
    <text transform="matrix(1 0 0 1 50 374.0034)" class="st7 st8">$(round(Int, x_min))</text>
    <text transform="matrix(1 0 0 1 400 374.0034)" class="st7 st8">$(round(Int, (x_min + x_max)/2))</text>
    <text transform="matrix(1 0 0 1 750 374.0034)" class="st7 st8">$(round(Int, x_max))</text>"""

        y_axis_labels = """
    <!-- Y-axis ticks and labels -->
    <line class="st5" x1="45" y1="350" x2="55" y2="350"/>
    <text transform="matrix(1 0 0 1 8.0833 356.0009)" class="st7 st8">0</text>

    <line class="st5" x1="45" y1="300" x2="55" y2="300"/>
    <text transform="matrix(1 0 0 1 8.0833 306.0009)" class="st7 st8">$(round(Int, y_max/5/1000))k</text>

    <line class="st5" x1="45" y1="250" x2="55" y2="250"/>
    <text transform="matrix(1 0 0 1 8.0833 256.0009)" class="st7 st8">$(round(Int, 2*y_max/5/1000))k</text>

    <line class="st5" x1="45" y1="200" x2="55" y2="200"/>
    <text transform="matrix(1 0 0 1 8.0833 206.0009)" class="st7 st8">$(round(Int, 3*y_max/5/1000))k</text>

    <line class="st5" x1="45" y1="150" x2="55" y2="150"/>
    <text transform="matrix(1 0 0 1 8.0833 156.0009)" class="st7 st8">$(round(Int, 4*y_max/5/1000))k</text>

    <line class="st5" x1="45" y1="100" x2="55" y2="100"/>
    <text transform="matrix(1 0 0 1 8.0833 106.0009)" class="st7 st8">$(round(Int, y_max/1000))k</text>"""

    # Scaling functions
    function scale_x(mz)
        return 50 + ((mz - x_min) * (700 / (x_max - x_min)))
    end

    function scale_y(intensity)
        # Clamp intensity to y_max
        clamped_intensity = min(intensity, y_max)
        return 350 - ((clamped_intensity / y_max) * 250)
    end

    # Generate first set of peaks (thin, grey)
    peaks1 = ["<!-- First set of peaks (grey) -->"]
    for i in 1:length(masses1)
        # Only plot if mass is within x boundaries
        if x_min ≤ masses1[i] ≤ x_max
            x = scale_x(masses1[i])
            y = scale_y(intensities1[i])
            peak = """<path class="st1" d="M$x,350 L$x,$y L$(x+2),$y L$(x+2),350"/>"""
            push!(peaks1, peak)
        end
    end
    # Generate second set of peaks (thick, orange)
    if orange_peaks
        peaks2 = ["<!-- Second set of peaks (orange) -->"]
        for i in 1:length(masses2)
            # Only plot if mass is within x boundaries
            if x_min ≤ masses2[i] ≤ x_max
                x = scale_x(masses2[i])
                y = scale_y(intensities2[i])
                peak = """<path class="st2" d="M$x,350 L$x,$y L$(x+4),$y L$(x+4),350"/>"""
                push!(peaks2, peak)
            end
        end
    else
        peaks2 = ["<!-- Second set of peaks (orange) -->"]
        for i in 1:length(masses2)
            # Only plot if mass is within x boundaries
            if x_min ≤ masses2[i] ≤ x_max
                x = scale_x(masses2[i])
                y = scale_y(intensities2[i])
                peak = """<path class="st1" d="M$x,350 L$x,$y L$(x+2),$y L$(x+2),350"/>"""
                push!(peaks2, peak)
            end
        end
    end

    # Combine all elements
    svg_content = join([
        svg_header,
        grid_lines,
        axes_and_title,
        y_axis_labels,
        join(peaks1, "\n"),
        join(peaks2, "\n"),
        "</svg>"
    ], "\n")

    # Write to file
    open(output_file, "w") do io
        write(io, svg_content)
    end
end

# Example usage:
# create_mass_spectrum_svg(
#     masses1, intensities1, 
#     masses2, intensities2, 
#     "mass_spectrum.svg",
#     x_min=400.0, x_max=1400.0, y_max=71000.0
# )
# Example usage:
# masses1 = [822.4246, 823.4255, 823.9402]  # first set of masses
# intensities1 = [71043.35, 53232.38, 21332.94]  # first set of intensities
# masses2 = [822.4246, 823.4255]  # second set of masses
# intensities2 = [65000.0, 48000.0]  # second set of intensities
# create_mass_spectrum_svg(masses1, intensities1, masses2, intensities2, "mass_spectrum.svg")
function find_masses_within_ppm(mass_list::Vector{Float64}, query_mass::Float64, ppm_tolerance::Float64=10.0)
    # Calculate absolute mass tolerance
    mass_tolerance = (ppm_tolerance * query_mass) / 1_000_000
    
    # Find indices where masses are within tolerance
    matching_indices = findall(m -> abs(m - query_mass) ≤ mass_tolerance, mass_list)
    
    return matching_indices
end

scan_idx = 50138
test_spectra = DataFrame(
    Dict(
        "mz"=>Float64.(coalesce.(ms_table[:mz_array][scan_idx], 0.0)), 
        "intensity"=>Float64.(coalesce.(ms_table[:intensity_array][scan_idx], 0.0))
        )
)
sort!(test_spectra, :intensity, rev=true)
prosit_masses = [
    1049.493164
    865.3720093
    1162.577271
    978.4560547
    1275.661255
    484.3129578
    636.2657471
    764.3243408
    668.4341431
    505.2252808
    371.2288818
    597.3970337
    258.1448364
    781.5181885
    262.13974
    390.1983337
    171.1128082
    147.1128082
    1475.777466
]
test_spectra_indices = [first(x) for x in [find_masses_within_ppm(
    Float64.(coalesce.(test_spectra[!,"mz"],0)),
    x
) for x in prosit_masses] if length(x)>0]

#=
create_mass_spectrum_svg(
    test_spectra[!,:mz][8:15:end], 
    test_spectra[!,:intensity][8:15:end], 
    test_spectra[test_spectra_indices ,:mz],
    test_spectra[test_spectra_indices ,:intensity],
    "/Users/n.t.wamsley/Documents/thesis_defense/figures/frag_plot/draft_plots/DDA_spectrum_unnanotated.svg",
    x_min=400.0, x_max=1400.0, y_max=6000.0, orange_peaks = false
)

create_mass_spectrum_svg(
    test_spectra[!,:mz][8:15:end], 
    test_spectra[!,:intensity][8:15:end], 
    test_spectra[test_spectra_indices ,:mz],
    test_spectra[test_spectra_indices ,:intensity],
    "/Users/n.t.wamsley/Documents/thesis_defense/figures/frag_plot/draft_plots/DDA_spectrum_anotated.svg",
    x_min=400.0, x_max=1400.0, y_max=6000.0, orange_peaks = true
)

mz = [100.0756912
147.1128082
171.1128082
258.1448364
262.13974
371.2288818
390.1983337
484.3129578
505.2252808
597.3970337
636.2657471
638.3342896
668.4341431
764.3243408
781.5181885
865.3720093
882.5658569
978.4560547
1010.624451
1049.493164
1141.664917
1162.577271
1256.691895
1275.661255
1384.750488
1388.745361
1475.777466
1499.777466]
intensities = [
0.5
1
0.5
0.5
1
0.5
1
0.5
1
0.5
1
1
0.5
1
0.5
1
0.5
1
0.5
1
0.5
1
0.5
1
0.5
1
1
0.5
]
create_mass_spectrum_svg(
    Float64[400], 
    Float64[0], 
    mz,
    intensities,
    "/Users/n.t.wamsley/Documents/thesis_defense/figures/frag_plot/draft_plots/dumb_spectra_julia.svg",
    x_min=400.0, x_max=1400.0, y_max=1.0, orange_peaks = true
)


=#

ms1_scan_idx = 50138 + 48
test_ms1 = DataFrame(
    Dict(
        "mz"=>Float64.(coalesce.(ms_table[:mz_array][ms1_scan_idx], 0.0)), 
        "intensity"=>Float64.(coalesce.(ms_table[:intensity_array][ms1_scan_idx], 0.0))
        )
)
sort!(test_ms1, :intensity, rev=true)
create_mass_spectrum_svg(
    test_spectra[!,:mz], 
    test_spectra[!,:intensity], 
    Float64[400],
    Float64[0],
    "/Users/n.t.wamsley/Documents/thesis_defense/figures/frag_plot/draft_plots/MS1_spectrum_unnanotated.svg",
    x_min=300.0, x_max=1000.0, orange_peaks = false
)


function tic_plot_svg(
    retention_times::Vector{Float64},  
    intensities::Vector{Float64},      
    output_file::String;
    x_max::Float64=ceil(maximum(retention_times)),  # automatically set x_max from data
    y_max::Float64=maximum(intensities))
    
    # Round x_max up to nearest 5 for nice axis labels
    x_max = ceil(x_max / 1) * 1
    
    # Calculate number of x-axis divisions (aim for roughly 5-minute intervals)
    n_x_divisions = Int(floor(x_max / 5))
    
    svg_header = """<?xml version="1.0" encoding="utf-8"?>
    <svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px"
        viewBox="0 0 800 400" style="enable-background:new 0 0 800 400;" xml:space="preserve">
    <style type="text/css">
        .grid{fill:none;stroke:#CCCCCC;stroke-width:0.5;}
        .axis{fill:none;stroke:#000000;stroke-width:1;}
        .peak{fill:#000000;}
        .text{font-family:Arial;}
        .label{font-size:12px;}
        .title{font-size:16px;font-weight:bold;}
    </style>

    <!-- Background -->
    <rect width="800" height="400" fill="white"/>"""

    # Rest of the header code remains the same...
    labels = """
        <!-- Title and Labels -->
        <text x="400" y="30" text-anchor="middle" class="text title">Total Ion Current</text>
        <text x="400" y="390" text-anchor="middle" class="text label">Time (min)</text>
        <text x="20" y="200" transform="rotate(-90,20,200)" text-anchor="middle" class="text label">Relative Abundance</text>"""

    # Grid lines
    grid_lines = ["""<!-- Grid lines -->"""]
    for i in 0:5
        y_pos = 350 - (i * 50)
        push!(grid_lines, """<line class="grid" x1="50" y1="$(y_pos)" x2="750" y2="$(y_pos)"/>""")
    end
    
    # Y-axis labels
    y_labels = ["""<!-- Y-axis labels -->"""]
    for i in 0:5
        y_pos = 350 - (i * 50)
        y_value = round(Int, (i * y_max) / 5)
        push!(y_labels, """<text x="45" y="$(y_pos+5)" text-anchor="end" class="text label">$(y_value)</text>""")
    end

    # X-axis labels with dynamic range
    x_labels = ["""<!-- X-axis labels -->"""]
    for i in 0:n_x_divisions
        x_pos = 50 + (i * (700/n_x_divisions))
        x_value = round(Int, i * (x_max/n_x_divisions))
        push!(x_labels, """<text x="$(x_pos)" y="370" text-anchor="middle" class="text label">$(x_value)</text>""")
    end

    # Rest of the code remains the same...
    axes = """
        <!-- Axes -->
        <line class="axis" x1="50" y1="350" x2="750" y2="350"/>
        <line class="axis" x1="50" y1="50" x2="50" y2="350"/>"""

    function scale_x(time)
        return 50 + ((time / x_max) * 700)
    end

    function scale_y(intensity)
        return 350 - ((intensity / y_max) * 250)
    end

    peaks = ["""<!-- TIC peaks -->"""]
    for i in 1:length(retention_times)
        x = scale_x(retention_times[i])
        y = scale_y(intensities[i])
        
        # Draw a vertical line from the intensity down to the x-axis (y=350)
        # Make the bars 1px wide
        peak = """<path class="st1" d="M$x,350 L$x,$y L$(x+1),$y L$(x+1),350 Z"/>"""
        push!(peaks, peak)
    end


    svg_content = join([
        svg_header,
        labels,
        join(grid_lines, "\n"),
        join(y_labels, "\n"),
        join(x_labels, "\n"),
        axes,
        join(peaks, "\n"),
        "</svg>"
    ], "\n")

    open(output_file, "w") do io
        write(io, svg_content)
    end
end
# Example usage:
# rt = collect(0.0:0.1:35.0)
# intensities = [exp(-(x-10)^2/5) + 0.8*exp(-(x-15)^2/2) + 1.2*exp(-(x-25)^2/8) + 0.2*rand() for x in rt]
# tic_plot_svg(rt, intensities, "tic_plot.svg")

# Example usage with toy data:
rt = collect(0.0:0.1:35.0)  # retention times from 0 to 35 minutes
# Generate some example peaks
intensities = [
    exp(-(x-10)^2/5) +     # peak at 10 minutes
    0.8*exp(-(x-15)^2/2) + # peak at 15 minutes
    1.2*exp(-(x-25)^2/8) + # peak at 25 minutes
    0.2*rand()             # random noise
    for x in rt
]
# Scale to percentage
intensities = 100 * intensities / maximum(intensities)

# Create the plot

rt = [Float64(ms_table[:retentionTime][i]) for i in range(1, length(ms_table[:retentionTime])) if (ms_table[:retentionTime][i] < 7.0)&(ms_table[:retentionTime][i] > 1.0)&(ms_table[:msOrder][i] == 1)]
TIC = [Float64(ms_table[:TIC][i]) for i in range(1, length(ms_table[:TIC])) if (ms_table[:retentionTime][i] < 7.0)&(ms_table[:retentionTime][i] > 1.0)&(ms_table[:msOrder][i] == 1)]
plot(rt, TIC)

tic_plot_svg(rt, TIC, "/Users/n.t.wamsley/Desktop/tic_plot2.svg", x_max = 7.0)

plot(rt, TIC)