abstract type AbstractPioneerPlotSpec end

struct PlotAnnotationSpec
    x::Float32
    y::Float32
    text::String
    halign::Symbol
    size::Int
end

struct PlotSeriesSpec
    x::Vector{Float32}
    y::Vector{Float32}
    label::Union{Nothing, String}
    seriestype::Symbol
    color::Union{Nothing, Symbol}
    alpha::Float32
    linewidth::Float32
    linestyle::Symbol
    markersize::Float32
end

struct MultiSeriesPlotSpec <: AbstractPioneerPlotSpec
    title::String
    xlabel::String
    ylabel::String
    series::Vector{PlotSeriesSpec}
    legend::Union{Nothing, Symbol}
    grid::Bool
    x_min::Union{Nothing, Float32}
    x_max::Union{Nothing, Float32}
    y_min::Union{Nothing, Float32}
    y_max::Union{Nothing, Float32}
    right_margin_px::Int
    annotations::Vector{PlotAnnotationSpec}
end

struct RtAlignmentPlotSpec <: AbstractPioneerPlotSpec
    title::String
    rt::Vector{Float32}
    irt::Vector{Float32}
    curve_rt::Vector{Float32}
    curve_irt::Vector{Float32}
    curve_label::Union{Nothing, String}
    curve_color::Symbol
    curve_style::Symbol
    curve_width::Float32
    annotation_text::Union{Nothing, String}
    annotation_x::Union{Nothing, Float32}
    annotation_y::Union{Nothing, Float32}
end

struct MassErrorPlotSpec <: AbstractPioneerPlotSpec
    title::String
    errs::Vector{Float32}
    center::Float32
    left_spread::Float32
    right_spread::Float32
    x_label::String
    y_label::String
    orientation::Symbol
    center_label::Union{Nothing, String}
    bound_label::Union{Nothing, String}
    annotation_text::Union{Nothing, String}
    annotation_x::Union{Nothing, Float32}
    annotation_y::Union{Nothing, Float32}
end
