# XRayLabTool.jl

*Material property calculations for X-ray interactions*

XRayLabTool is a Julia package that calculates X-ray optical properties of materials from their chemical composition and density. It is designed for synchrotron scientists, materials researchers, and X-ray optics developers.

## What it computes

Given a chemical formula (e.g., `"SiO2"`), an energy range (0.03–30 keV), and a mass density, XRayLabTool returns:

| Property | Symbol | Units |
|----------|--------|-------|
| Dispersion coefficient | δ | dimensionless |
| Absorption coefficient | β | dimensionless |
| Scattering factors (summed over the formula unit) | f₁, f₂ | electrons/formula unit |
| Critical angle | θc | degrees |
| Attenuation length | L | cm |
| Scattering length density | SLD | Å⁻² |
| Electron density | ρe | Å⁻³ |
| Molecular weight | MW | g/mol |

## Quick example

```julia
using XRayLabTool

# Single material: SiO2 at Cu Kα energy
result = calculate_single_material_properties("SiO2", [8.04], 2.2)

println("δ = ", result.dispersion[1])         # 7.22e-6
println("θc = ", result.critical_angle[1])     # 0.218°
println("L = ", result.attenuation_length[1])  # attenuation length in cm

# Multiple materials in one call
results = calculate_xray_properties(
    ["SiO2", "Al2O3", "Au"],
    collect(5.0:0.5:20.0),
    [2.2, 3.95, 19.3],
)

sio2 = results[1]  # same order as the input formulas
au   = results[3]
```

## Data source

Atomic scattering factors are based on the [CXRO](http://www.cxro.lbl.gov) / Henke tables, covering all elements from H to U in the 30 eV – 30 keV range. The package uses PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) interpolation to evaluate f₁ and f₂ at arbitrary energies within this range.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/imewei/XRayLabTool.git")
```

## Contents

```@contents
Pages = ["getting_started.md", "physics.md", "api.md", "migration.md", "changelog.md"]
Depth = 2
```
