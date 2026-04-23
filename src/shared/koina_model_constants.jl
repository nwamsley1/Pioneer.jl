const MODEL_CONFIGS = Dict(
    "unispec" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = InstrumentSpecificModel("unispec"),
        instruments = Set(["QE", "QEHFX", "LUMOS", "ELITE", "VELOS", "NONE"]),
    ),
    "altimeter" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = SplineCoefficientModel("altimeter"),
        instruments = Set([]),
    ),
    "prosit_2020_hcd" => (
        annotation_type = GenericFragAnnotation("y1+1"),
        model_type = InstrumentAgnosticModel("prosit_2020_hcd"),
        instruments = Set([]),
    ),
    "AlphaPeptDeep" => (
        annotation_type = GenericFragAnnotation("y1+1"),
        model_type = InstrumentSpecificModel("AlphaPeptDeep"),
        instruments = Set(["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]),
    ),
)

const KOINA_URLS = Dict(
    "unispec" => "https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer",
    "prosit_2020_hcd" => "https://koina.wilhelmlab.org:443/v2/models/Prosit_2020_intensity_HCD/infer",
    "AlphaPeptDeep" => "https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_ms2_generic/infer",
    "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer",
    "altimeter" => "https://koina.wilhelmlab.org:443/v2/models/Altimeter_2024_splines_index/infer",
)
