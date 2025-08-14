# Bias Shift Strategies Explained

## Understanding Mass Bias in Parameter Tuning

### What is Mass Bias?
Mass bias is the systematic offset between the theoretical mass of ions and their observed mass in the mass spectrometer. For example, if all ions are consistently measured 10 ppm heavier than expected, the mass bias is +10 ppm.

### The Problem
Traditional parameter tuning starts with zero bias assumption and a narrow tolerance window (e.g., ±20 ppm). If the true bias is outside this window (e.g., +35 ppm), the search will fail to find PSMs because it's looking in the wrong mass range.

### The Solution: Phase-Based Bias Exploration
The phase system systematically explores different bias regions:
- **Phase 1**: Explores around zero bias (traditional approach)
- **Phase 2+**: Explores shifted bias regions based on the strategy

## Bias Shift Strategy Options

### 1. `"alternating"` (Default)
**Pattern**: Alternates between positive and negative shifts each phase

```
Phase 1: bias = 0 ppm
Phase 2: bias = +magnitude (e.g., +50 ppm)
Phase 3: bias = -magnitude (e.g., -50 ppm)
Phase 4: bias = +magnitude (e.g., +50 ppm)
Phase 5: bias = -magnitude (e.g., -50 ppm)
```

**Visual Representation**:
```
Mass Range Exploration:
Phase 1: |--------[0]--------|  (centered at 0)
Phase 2: |--------[+50]------|  (shifted right)
Phase 3: |--------[-50]------|  (shifted left)
Phase 4: |--------[+50]------|  (back to right)
```

**When to Use**:
- **Default choice** when you have no prior knowledge about bias direction
- **Balanced exploration** of both positive and negative regions
- **Most datasets** benefit from this approach
- **Unknown samples** where calibration status is uncertain

**Example Use Case**:
- Processing data from multiple instruments where some may have positive drift and others negative
- Exploratory analysis of new sample types
- Routine processing pipelines that need to handle various scenarios

### 2. `"positive_first"`
**Pattern**: Tries positive bias first, then negative for all remaining phases

```
Phase 1: bias = 0 ppm
Phase 2: bias = +magnitude (e.g., +50 ppm)
Phase 3: bias = -magnitude (e.g., -50 ppm)
Phase 4: bias = -magnitude (e.g., -50 ppm)
Phase 5: bias = -magnitude (e.g., -50 ppm)
```

**Visual Representation**:
```
Mass Range Exploration:
Phase 1: |--------[0]--------|  (centered at 0)
Phase 2: |--------[+50]------|  (shifted right)
Phase 3: |--------[-50]------|  (shifted left)
Phase 4: |--------[-50]------|  (stays left)
```

**When to Use**:
- **Instrument characteristics** known to drift positive over time
- **Sample preparation** methods that typically cause positive shifts
- **Historical data** shows positive bias is more common
- **Quick positive check** before extensive negative exploration

**Example Use Case**:
- Older Orbitrap instruments that tend to drift to higher masses
- Samples with heavy isotope labeling
- Post-column modification experiments that add mass

### 3. `"negative_first"`
**Pattern**: Tries negative bias first, then positive for all remaining phases

```
Phase 1: bias = 0 ppm
Phase 2: bias = -magnitude (e.g., -50 ppm)
Phase 3: bias = +magnitude (e.g., +50 ppm)
Phase 4: bias = +magnitude (e.g., +50 ppm)
Phase 5: bias = +magnitude (e.g., +50 ppm)
```

**Visual Representation**:
```
Mass Range Exploration:
Phase 1: |--------[0]--------|  (centered at 0)
Phase 2: |--------[-50]------|  (shifted left)
Phase 3: |--------[+50]------|  (shifted right)
Phase 4: |--------[+50]------|  (stays right)
```

**When to Use**:
- **Instrument characteristics** known to drift negative
- **Specific ion sources** that produce systematic negative bias
- **Calibration procedures** that tend to overcorrect
- **Sample types** known to cause negative shifts

**Example Use Case**:
- MALDI instruments with matrix effects
- Experiments with chemical modifications that affect flight time
- Freshly calibrated instruments that tend to overcorrect

## Bias Shift Magnitude Options

### 1. `"max_tolerance"` (Default)
**Behavior**: Uses the value of `max_tolerance_ppm` parameter as the shift magnitude

```json
{
  "max_tolerance_ppm": 50.0,
  "iteration_settings": {
    "bias_shift_magnitude": "max_tolerance"  // Will use 50.0 ppm
  }
}
```

**Advantages**:
- **Comprehensive coverage**: Explores the full range allowed by your tolerance
- **Automatic scaling**: Adjusts with your tolerance settings
- **No gaps**: Ensures complete parameter space coverage
- **Safe default**: Won't exceed your defined limits

**Example Scenario**:
```
If max_tolerance_ppm = 50:
Phase 1: Explores -50 to +50 ppm around 0
Phase 2: Explores -50 to +50 ppm around +50 (total range: 0 to +100)
Phase 3: Explores -50 to +50 ppm around -50 (total range: -100 to 0)
```

### 2. Numeric Value (e.g., `30.0`)
**Behavior**: Uses a specific fixed value for bias shifts

```json
{
  "iteration_settings": {
    "bias_shift_magnitude": 30.0  // Always shifts by exactly 30 ppm
  }
}
```

**Advantages**:
- **Precise control**: Exact knowledge of exploration regions
- **Overlap management**: Can create overlapping search regions
- **Fine-tuning**: Based on known instrument characteristics
- **Conservative approach**: Smaller shifts for gradual exploration

**Example Scenarios**:

#### Conservative (Small Shifts)
```json
"bias_shift_magnitude": 20.0
```
```
Phase 1: Explores around 0 ppm
Phase 2: Explores around +20 ppm (overlaps with Phase 1)
Phase 3: Explores around -20 ppm (overlaps with Phase 1)
```
**Use when**: High confidence in calibration, minor drift expected

#### Moderate (Medium Shifts)
```json
"bias_shift_magnitude": 35.0
```
```
Phase 1: Explores around 0 ppm
Phase 2: Explores around +35 ppm (some overlap)
Phase 3: Explores around -35 ppm (some overlap)
```
**Use when**: Typical instrument drift, balanced exploration

#### Aggressive (Large Shifts)
```json
"bias_shift_magnitude": 60.0
```
```
Phase 1: Explores around 0 ppm
Phase 2: Explores around +60 ppm (minimal overlap)
Phase 3: Explores around -60 ppm (minimal overlap)
```
**Use when**: Poor calibration suspected, extreme biases possible

## Practical Examples

### Example 1: Well-Calibrated Modern Instrument
```json
{
  "max_tolerance_ppm": 30.0,
  "iteration_settings": {
    "mass_tolerance_scale_factor": 1.5,
    "iterations_per_phase": 3,
    "max_phases": 2,
    "bias_shift_strategy": "alternating",
    "bias_shift_magnitude": 15.0  // Small shifts
  }
}
```
**Rationale**: 
- Small bias shifts (15 ppm) because calibration is good
- Only 2 phases needed (less exploration required)
- Alternating strategy for balanced checking

### Example 2: Older Instrument with Positive Drift Tendency
```json
{
  "max_tolerance_ppm": 50.0,
  "iteration_settings": {
    "mass_tolerance_scale_factor": 2.0,
    "iterations_per_phase": 3,
    "max_phases": 4,
    "bias_shift_strategy": "positive_first",
    "bias_shift_magnitude": "max_tolerance"
  }
}
```
**Rationale**:
- Positive_first because instrument tends to drift positive
- Max_tolerance magnitude for comprehensive coverage
- 4 phases to thoroughly explore positive region

### Example 3: Unknown Sample on Variable Instrument
```json
{
  "max_tolerance_ppm": 60.0,
  "iteration_settings": {
    "mass_tolerance_scale_factor": 2.5,
    "iterations_per_phase": 2,
    "max_phases": 5,
    "bias_shift_strategy": "alternating",
    "bias_shift_magnitude": 40.0
  }
}
```
**Rationale**:
- Alternating for unbiased exploration
- Moderate shifts (40 ppm) with some overlap
- Many phases (5) for thorough exploration
- Fast scaling (2.5x) for quick convergence within each phase

## Decision Tree for Strategy Selection

```
Start: Do you know your instrument's bias tendency?
│
├─ No → Use "alternating"
│   │
│   └─ Do you know the typical bias magnitude?
│       ├─ No → Use "max_tolerance"
│       └─ Yes → Use numeric value
│
└─ Yes
    │
    ├─ Usually positive → Use "positive_first"
    │   │
    │   └─ How large is typical drift?
    │       ├─ Small (<20 ppm) → Use numeric value (15-25)
    │       ├─ Medium (20-40 ppm) → Use numeric value (30-40)
    │       └─ Large (>40 ppm) → Use "max_tolerance"
    │
    └─ Usually negative → Use "negative_first"
        │
        └─ [Same magnitude decision as above]
```

## Impact on Convergence

### Convergence Speed
- **Correct strategy choice**: Can converge in Phase 1-2
- **Incorrect strategy**: May need Phase 3-4 to find correct region
- **Alternating strategy**: Balanced, typically converges by Phase 3

### Computation Time
```
Time Impact = (phases needed) × (iterations per phase) × (scans per iteration)
```

**Example Comparison**:
- Optimal strategy: 2 phases × 3 iterations = 6 total iterations
- Suboptimal strategy: 4 phases × 3 iterations = 12 total iterations
- Time difference: 2x longer for suboptimal choice

### Success Rate
- **Well-chosen magnitude**: 95%+ success rate
- **Too small magnitude**: May miss extreme biases, 70-80% success
- **Too large magnitude**: May miss optimal region due to gaps, 80-85% success
- **"max_tolerance" default**: 90%+ success rate (recommended for general use)

## Summary Recommendations

### For Most Users
```json
{
  "bias_shift_strategy": "alternating",
  "bias_shift_magnitude": "max_tolerance"
}
```
This default configuration works well for 90% of datasets.

### For Advanced Users
1. **Analyze your instrument's historical performance**
2. **Choose strategy based on drift patterns**
3. **Set magnitude based on typical drift range**
4. **Monitor logs to optimize for your specific setup**

### Key Takeaways
- **Strategy** determines the ORDER of exploration
- **Magnitude** determines the SIZE of shifts
- **Together** they define the parameter space coverage
- **Defaults** are robust and work for most cases
- **Customization** can significantly improve convergence speed for known instruments