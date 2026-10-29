# Getting Started

## Installation

XRayLabTool requires Julia 1.10 or later.

```julia
using Pkg
Pkg.add(url="https://github.com/imewei/XRayLabTool.git")
```

## Single material calculation

Use [`calculate_single_material_properties`](@ref) for one material:

```julia
using XRayLabTool

# SiO2 at energies from 5 to 20 keV, density 2.2 g/cm³
result = calculate_single_material_properties("SiO2", collect(5.0:1.0:20.0), 2.2)

# Access scalar properties
result.formula              # "SiO2"
result.molecular_weight     # 60.08 g/mol
result.number_of_electrons  # 30.0
result.electron_density     # electrons/Å³

# Access energy-dependent vector properties
result.energy               # [5.0, 6.0, ..., 20.0] keV
result.wavelength           # corresponding wavelengths in Å
result.dispersion           # δ values
result.absorption           # β values
result.f1                   # real scattering factor (summed over the formula unit)
result.f2                   # imaginary scattering factor (summed over the formula unit)
result.critical_angle       # critical angles in degrees
result.attenuation_length   # 1/e attenuation lengths in cm
result.real_sld             # real part of SLD in Å⁻²
result.imag_sld             # imaginary part of SLD in Å⁻²
```

## Batch calculation

Use [`calculate_xray_properties`](@ref) for multiple materials at once. This uses multi-threading when Julia is started with multiple threads (`julia --threads=4`):

```julia
formulas  = ["SiO2", "Al2O3", "Fe2O3", "Au"]
energies  = collect(1.0:0.1:20.0)
densities = [2.2, 3.95, 5.24, 19.3]

results = calculate_xray_properties(formulas, energies, densities)

# Access by index (same order as input formulas)
sio2 = results[1]   # SiO2
au   = results[4]   # Au

println("Au critical angle at 8 keV: ", au.critical_angle[71])
```

## Input constraints

| Parameter | Constraint |
|-----------|------------|
| Formula | Case-sensitive: `"CO"` = Carbon Monoxide, `"Co"` = Cobalt |
| Energy | Finite and in range **0.03 – 30.0 keV**; `result.energy` is returned sorted ascending |
| Density | Positive, finite `Float64` in g/cm³ |
| Formula count | Must match density count in batch mode |
| Element counts | Positive and finite; zero counts (`H0`) and counts that overflow or underflow the molecular weight are rejected |
| Nesting | At most 64 levels of parentheses |

All violations raise `ArgumentError`.

## Chemical formula syntax

The parser supports:

| Formula | Interpretation |
|---------|---------------|
| `"SiO2"` | 1 Si + 2 O |
| `"Al2O3"` | 2 Al + 3 O |
| `"CaCO3"` | 1 Ca + 1 C + 3 O |
| `"Ca(OH)2"` | Parenthesized groups: 1 Ca + 2 O + 2 H |
| `"Ca3(PO4)2"` | Nested groups: 3 Ca + 2 P + 8 O |
| `"Fe((OH)2)3"` | Deeply nested: 1 Fe + 6 O + 6 H |
| `"H0.5He0.5"` | Fractional stoichiometry |
| `"C100H200"` | Large subscripts |

## Multi-threading

For best batch performance, start Julia with multiple threads:

```bash
julia --threads=auto   # Uses all available cores
julia --threads=4      # Fixed 4 threads
```

The batch function pre-populates element caches on the main thread, then distributes formula calculations across threads using dynamic scheduling.

## Memory management

For long-running applications processing many materials, clear the element cache periodically:

```julia
XRayLabTool.clear_caches!()
```
