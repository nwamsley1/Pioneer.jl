```@meta
CurrentModule = Pioneer
```

# Pioneer.jl

Pioneer.jl is a Julia package for analyzing Data-Independent Acquisition (DIA) mass spectrometry data.

## Features

- Fast and accurate peptide identification
- Sophisticated retention time alignment
- Advanced quantification methods
- Flexible parameter configuration
- Comprehensive result reporting

## Quick Links

- [Installation Guide](@ref installation-guide)
- [Quick Start Tutorial](@ref quick-start-tutorial)
- [API Reference](@ref api-reference)
- [Algorithm Documentation](@ref algorithm-documentation)

## Package Overview

Pioneer.jl implements a multi-stage search strategy for DIA data analysis:

1. Parameter Tuning
2. NCE/Quadrupole Optimization
3. First Pass Search
4. Quantification
5. Statistical Validation

For detailed information, see the [User Guide](@ref user-guide).

```@index
```

```@autodocs
Modules = [Pioneer]
Private = false
Order = [:module, :type, :function, :macro]
```
