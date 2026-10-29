# XRayLabTool.jl

**Material property calculations for X-ray interactions**

XRayLabTool is a Julia package that calculates X-ray optical properties of materials from their chemical composition and density. It is designed for synchrotron scientists, materials researchers, and X-ray optics developers.

## Features

- Compute optical constants (delta, beta), scattering factors (f1, f2), critical angles, attenuation lengths, and scattering length densities
- Single and batch multi-material calculations with multi-threaded execution
- Based on CXRO/Henke atomic scattering factor tables (H through U, 30 eV -- 30 keV)
- PCHIP interpolation for smooth, monotonicity-preserving evaluation at arbitrary energies
- Parenthesized chemical formula support: `Ca(OH)2`, `Ca3(PO4)2`, `Fe(OH)2`, `Fe(OH)3`
- Debug logging with topic categories (enable with `ENV["JULIA_DEBUG"] = "XRayLabTool"`)
- Lightweight: only 3 third-party dependencies (PCHIPInterpolation, Mendeleev, Unitful) plus the `Logging` standard library

## Requirements

- Julia 1.10 or later (CI runs on Julia 1.10, the latest release, and the upcoming prerelease)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/imewei/XRayLabTool.git")
```

## Quick start

### Single material

```julia
using XRayLabTool

result = calculate_single_material_properties("SiO2", [8.04, 10.0, 12.0], 2.2)

result.formula              # "SiO2"
result.molecular_weight     # 60.08 g/mol
result.dispersion           # delta values at each energy
result.critical_angle       # critical angles in degrees
result.attenuation_length   # 1/e attenuation lengths in cm
result.real_sld             # real part of SLD in 1/A^2
```

### Multiple materials

```julia
results = calculate_xray_properties(
    ["SiO2", "Al2O3", "Au"],
    collect(5.0:0.5:20.0),
    [2.2, 3.95, 19.3],
)

sio2 = results[1]   # SiO2 (same order as input)
au   = results[3]   # Au
```

## Input parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `formula` | `String` or `Vector{String}` | Case-sensitive chemical formula(s). `"CO"` = Carbon Monoxide, `"Co"` = Cobalt |
| `energy` | `Vector{Float64}` | X-ray energies in keV (range: 0.03 -- 30.0) |
| `density` | `Float64` or `Vector{Float64}` | Positive, finite mass density in g/cm3 (one per formula in batch mode) |

Invalid input raises `ArgumentError`: energies outside 0.03 -- 30 keV or non-finite, non-positive densities, unknown elements, and formulas with zero, non-finite or unrepresentable element counts or more than 64 levels of nested parentheses. Results are returned with `energy` sorted in ascending order.

## Output: `XRayResult`

| Field | Type | Description |
|-------|------|-------------|
| `formula` | `String` | Chemical formula |
| `molecular_weight` | `Float64` | Molecular weight (g/mol) |
| `number_of_electrons` | `Float64` | Electrons per formula unit |
| `mass_density` | `Float64` | Mass density (g/cm3) |
| `electron_density` | `Float64` | Electron density (1/A3) |
| `energy` | `Vector{Float64}` | X-ray energy (keV) |
| `wavelength` | `Vector{Float64}` | Wavelength (A) |
| `dispersion` | `Vector{Float64}` | Dispersion coefficient delta |
| `absorption` | `Vector{Float64}` | Absorption coefficient beta |
| `f1` | `Vector{Float64}` | Real scattering factor, summed over the formula unit |
| `f2` | `Vector{Float64}` | Imaginary scattering factor, summed over the formula unit |
| `critical_angle` | `Vector{Float64}` | Critical angle (degrees) |
| `attenuation_length` | `Vector{Float64}` | Attenuation length (cm) |
| `real_sld` | `Vector{Float64}` | Real part of SLD (1/A2) |
| `imag_sld` | `Vector{Float64}` | Imaginary part of SLD (1/A2) |

## Multi-threading

The batch function distributes calculations across threads automatically. Start Julia with multiple threads for best performance:

```bash
julia --threads=auto
```

## Development

```bash
# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run tests with threading
julia --threads=4 --project=. -e 'using Pkg; Pkg.test()'

# Build documentation
julia --project=docs docs/make.jl

# Format source code
julia -e 'using JuliaFormatter; format("src")'
```

## Documentation

Full documentation including physics background, API reference, and migration guide is available at the [documentation site](https://imewei.github.io/XRayLabTool/).

## References

- B.L. Henke, E.M. Gullikson, and J.C. Davis, *X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92*, Atomic Data and Nuclear Data Tables **54**, 181-342 (1993).
- [CXRO - Center for X-ray Optics](http://www.cxro.lbl.gov)
- [NIST X-Ray Form Factor Tables](https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables)

> This package is translated from a MATLAB script originally developed by **Zhang Jiang** at the **Advanced Photon Source**, Argonne National Laboratory.

## License

MIT
