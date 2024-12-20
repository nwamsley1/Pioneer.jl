# src/constants.jl

#=
Core physical constants
All masses in unified atomic mass units (Da)
=#
const PROTON_MASS = Float32(1.0072764)
const NEUTRON_MASS = Float32(1.00335)
const ELECTRON_MASS = Float32(0.0005486)
const WATER_MASS = Float32(18.010565)
const H2O_MASS = WATER_MASS
const NH3_MASS = Float32(17.026549)
const CO_MASS = Float32(27.994915)
const CO2_MASS = Float32(43.989830)
const H3PO4_MASS = Float32(97.976896)

# Amino acid monoisotopic masses
const AA_MASSES = Dict{Char, Float32}(
    'A' => Float32(71.03711),
    'R' => Float32(156.10111),
    'N' => Float32(114.04293),
    'D' => Float32(115.02694),
    'C' => Float32(103.00919),
    'E' => Float32(129.04259),
    'Q' => Float32(128.05858),
    'G' => Float32(57.02146),
    'H' => Float32(137.05891),
    'I' => Float32(113.08406),
    'L' => Float32(113.08406),
    'K' => Float32(128.09496),
    'M' => Float32(131.04049),
    'F' => Float32(147.06841),
    'P' => Float32(97.05276),
    'S' => Float32(87.03203),
    'T' => Float32(101.04768),
    'W' => Float32(186.07931),
    'Y' => Float32(163.06333),
    'V' => Float32(99.06841),
    'U' => Float32(150.95363),  # Selenocysteine
    'O' => Float32(237.14773)   # Pyrrolysine
)

# Element masses
const ELEMENT_MASSES = Dict{String, Float32}(
    "H" => Float32(1.007825035),
    "C" => Float32(12.0),
    "N" => Float32(14.003074),
    "O" => Float32(15.994915),
    "P" => Float32(30.973762),
    "S" => Float32(31.972071)
)

# Isotope abundances for mass spec calculations
const ISOTOPE_ABUNDANCES = Dict{String, Float32}(
    "C13" => Float32(0.0107),   # Carbon-13
    "N15" => Float32(0.00364),  # Nitrogen-15
    "O18" => Float32(0.00205),  # Oxygen-18
    "S34" => Float32(0.0442)    # Sulfur-34
)

# Prediction model configurations
const PREDICTION_MODELS = Set([
    "unispec",
    "prosit_2020_hcd", 
    "AlphaPeptDeep"
])

const MODEL_CONFIGS = Dict(
    "unispec" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = InstrumentSpecificModel("unispec"),
        instruments = Set(["QE", "QEHFX", "LUMOS", "ELITE", "VELOS", "NONE"])
    ),
    "prosit_2020_hcd" => (
        annotation_type = GenericFragAnnotation("y1+1"),
        model_type = InstrumentAgnosticModel("prosit_2020_hcd"),
        instruments = Set([])
    ),
    "AlphaPeptDeep" => (
        annotation_type = GenericFragAnnotation("y1+1"),
        model_type = InstrumentSpecificModel("AlphaPeptDeep"),
        instruments = Set(["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"])
    )
)

# Koina API endpoints
const KOINA_URLS = Dict(
    "unispec" => "https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer",
    "prosit_2020_hcd" => "https://koina.wilhelmlab.org:443/v2/models/Prosit_2020_intensity_HCD/infer",
    "AlphaPeptDeep" => "https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_ms2_generic/infer",
    "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"
)

# Charge state adjustment factors
const CHARGE_ADJUSTMENT_FACTORS = Float32[1.0, 0.9, 0.85, 0.8, 0.75]

# Common regular expressions
const REGEX_PATTERNS = Dict(
    "trypsin" => r"[KR][^P|$]",
    "lysC" => r"[K][^P|$]",
    "mod_pattern" => r"\((\d+),(\w+),([^)]+)\)",
    "internal_ion" => r"Int/([A-Z]+)([^/]+)/(\d+)",
    "immonium_ion" => r"^I[A-Z]{2}",
    "isotope_pattern" => r"\+([0-9]+)i"
)

# Default modification masses
const DEFAULT_MODS = Dict{String, Float32}(
    "Oxidation" => Float32(15.994915),
    "Carbamidomethyl" => Float32(57.021464),
    "Phospho" => Float32(79.966331),
    "Acetyl" => Float32(42.010565),
    "GlyGly" => Float32(114.042927),
    "TMT6plex" => Float32(229.162932),
    "TMT10plex" => Float32(229.162932),
    "iTRAQ4plex" => Float32(144.102063),
    "iTRAQ8plex" => Float32(304.205360)
)

# Fragment ion types
const FRAGMENT_TYPES = Dict(
    "a" => -CO_MASS,          # a-ion mass shift
    "b" => Float32(0.0),      # b-ion (no mass shift)
    "c" => Float32(17.026549),# c-ion mass shift
    "x" => CO_MASS + Float32(17.026549), # x-ion mass shift
    "y" => Float32(18.010565),# y-ion mass shift
    "z" => -NH3_MASS         # z-ion mass shift
)

# Common neutral losses
const NEUTRAL_LOSSES = Dict(
    "H2O" => -H2O_MASS,
    "NH3" => -NH3_MASS,
    "CO" => -CO_MASS,
    "H3PO4" => -H3PO4_MASS,
    "CH3SOH" => Float32(-63.998285),  # Common loss from oxidized methionine
    "HPO3" => Float32(-79.966331)     # Common loss from phosphorylated residues
)