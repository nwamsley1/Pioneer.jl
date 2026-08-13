# Constants for parameter tuning search methods.
# These values were previously read from JSON configuration and are now hardcoded
# since they rarely need user adjustment.

# Shared tuning fragment settings (used by ParameterTuning, NceTuning, QuadTuning)
const TUNING_MIN_SCORE = [UInt8(8), UInt8(7), UInt8(6), UInt8(5)]  # top4 required + count threshold tiers
const TUNING_MIN_SPECTRAL_CONTRAST = Float32(0.5)
const TUNING_MIN_LOG2_RATIO = Float32(1.5)
const TUNING_MIN_TOPN_OF_M = (Int64(3), Int64(3))
const TUNING_MAX_FRAG_RANK = UInt8(25)
const PARAM_TUNING_MAX_FRAG_RANK = UInt8(7)
const TUNING_RELATIVE_IMPROVEMENT_THRESHOLD = Float32(1.25)
const TUNING_N_FRAG_ISOTOPES = Int64(1)
const TUNING_INTENSITY_FILTER_QUANTILE = Float32(0.5)

# ParameterTuningSearch-specific
const TUNING_FRAG_ERR_QUANTILE = Float32(0.01)
const TUNING_MIN_SAMPLES = Int64(1200)          # PSM target for collection phase
const TUNING_MAX_PRESEARCH_ITERS = Int64(10)
const TUNING_MAX_Q_VALUE = Float32(0.01)
const TUNING_TOPN_PEAKS = Int64(200)            # Top-N intensity peak filter for wide scout
const TUNING_MAX_FRAGS_FOR_MASS_ERR = UInt8(3)  # Fragments per PSM for mass error estimation
const TUNING_MIN_COLLECT_SCANS = Int64(5000)    # Minimum initial scan estimate for collection
const TUNING_CALIBRATION_BIN_SIZE = Int64(200)  # Equal-count bins for mass-error calibration summaries
const TUNING_MZ_BIAS_BIN_SIZE = Int64(100)      # Denser bins for m/z-dependent bias medians
const TUNING_MZ_BIAS_KNOTS = Int64(24)          # Knots for binned m/z-dependent bias spline
const TUNING_MZ_BIAS_LAMBDA = Float64(0.1)      # Smoothness penalty for m/z-dependent bias spline
const TUNING_RT_BIAS_BIN_SIZE = Int64(100)      # Denser bins for RT-dependent bias medians
const TUNING_RT_BIAS_KNOTS = Int64(48)          # Knots for raw RT bias spline
const TUNING_RT_BIAS_LAMBDA = Float64(0.03)     # Smoothness penalty for RT-dependent bias spline
const TUNING_BIAS_EXTRAP_EDGE_POINTS = Int64(5) # Binned medians per edge for robust linear extrapolation

# Wide scout constants
const WIDE_SCOUT_TOL_PPM = Float32(100.0)       # ±100 ppm wide scout tolerance
const WIDE_SCOUT_TARGET_PSMS = Int64(1000)      # PSM target at 1% FDR for wide scout
const WIDE_SCOUT_INITIAL_SCANS = Int64(2000)    # Starting scan count for wide scout
const WIDE_SCOUT_FALLBACK_TOL_PPM = Float32(30.0) # Fallback tolerance when <50 fragments

# Scout model fitting
const SCOUT_MIN_FRAGS = Int64(50)               # Minimum fragments to fit scout model
const SCOUT_INTENSITY_THRESHOLD = Int64(1000)   # PSMs at 1% FDR needed for intensity correction

# Model fitting thresholds
const MIN_PSMS_FOR_RT = Int64(300)              # Minimum PSMs to fit RT alignment model
const MIN_FRAGS_FOR_INTENSITY_MODEL = Int64(2000) # Minimum fragments for IntensityMassErrorModel

# Gaussian coverage for IntensityMassErrorModel tolerance (95% → k=1.96)
const TUNING_GAUSSIAN_COVERAGE = Float64(0.95)

# iRT tolerance: σ_iRT × this multiplier
const TUNING_IRT_TOL_SIGMA = Int64(3)

# Score tier backoff and top-N fragment requirement
const TUNING_SCORE_TIERS = (UInt8(8), UInt8(7), UInt8(6), UInt8(5))
const TUNING_N_REQUIRED_TOP = Int64(4)

# RT alignment (consumed only by ParameterTuningSearch)
const TUNING_LAMBDA_PENALTY = Float32(0.1)
const TUNING_RANSAC_THRESHOLD_PSMS = Int64(500)
const TUNING_MIN_PSMS_FOR_SPLINE = Int64(10)

# Quad tuning
const QUAD_TUNING_MIN_PSMS_PER_THOMPSON = Int64(250)
const QUAD_TUNING_MIN_FRAGMENTS = Int64(3)
const QUAD_TUNING_INITIAL_PERCENT = Float32(2.5)
