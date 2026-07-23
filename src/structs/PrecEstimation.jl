# Precursor estimation strategy traits.
#
# Used by the fused path (run_fused!) to dispatch getFragIsotopes! between
# partial vs full precursor capture variants. Originally defined alongside
# the deleted selectTransitions! dispatcher; only the trait types survived.

"""
    PrecEstimation

Abstract type for precursor isotope estimation strategy. Drives how
fragment ion isotope patterns are calculated by `getFragIsotopes!`.
"""
abstract type PrecEstimation end

"""
    PartialPrecCapture <: PrecEstimation

Use when only part of the precursor isotope envelope is captured in the
isolation window. Detailed Goldfarb-style abundance calculation that
accounts for partial isotopologue capture.
"""
struct PartialPrecCapture <: PrecEstimation end

"""
    FullPrecCapture <: PrecEstimation

Use when the full precursor isotope envelope is captured. Simplified
direct spline-evaluation, assumes complete capture.
"""
struct FullPrecCapture <: PrecEstimation end

"""
    PartialPrecCaptureNorm <: PrecEstimation

Sciex ZT scanning-DIA variant of `PartialPrecCapture`: the quad transmission still
warps the fragment isotope *ratios* (via `getFragAbundance!`), but the fragment's
total predicted intensity is renormalized to be transmission-INDEPENDENT. This keeps
the empirical transmission out of the design-column magnitude, so the fitted per-bin
weight tracks the actual ion current (the meta-scan triangle the ZT shape features
key on) instead of being flattened by the transmission.
"""
struct PartialPrecCaptureNorm <: PrecEstimation end
