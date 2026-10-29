# Migration Guide

## v0.6 → v0.7

### Breaking changes

**`calculate_xray_properties` returns `Vector{XRayResult}` instead of `Dict{String, XRayResult}`**:
Previously, duplicate formulas with different densities would silently overwrite each other in the Dict. The new Vector return type preserves all entries in input order.

```julia
# Old (v0.6)
results = calculate_xray_properties(formulas, energies, densities)
sio2 = results["SiO2"]

# New (v0.7)
results = calculate_xray_properties(formulas, energies, densities)
sio2 = results[1]  # index-based access
```

**Batch validation is now strict**: Invalid formulas in a batch call now throw a single `ArgumentError` listing all invalid formulas before any computation begins. Previously, the first invalid formula crashed the call and subsequent ones were never checked.

**Single-material path now sorts energies**: `calculate_single_material_properties` now returns results with energies sorted in ascending order, matching the batch API behavior.

### Non-breaking improvements

- Parenthesized chemical formulas are now supported: `Ca(OH)2`, `Ca3(PO4)2`, `Fe((OH)2)3`
- Physical constants updated to CODATA 2018 / 2019 SI exact values (~2e-7 relative change)
- Thread-safe cache reads (locked Dict access prevents segfault under concurrent use)
- Critical angle returns 0.0 instead of throwing `DomainError` when delta is negative near absorption edges

## v0.5 → v0.6

### Breaking changes

**Removed dependencies**: `CSV.jl` and `DataFrames.jl` are no longer dependencies. If you accessed internal caches that returned `DataFrame` objects, they now return interpolator tuples instead.

**Internal cache API changed**:
- `F1F2_TABLE_CACHE` (stored `DataFrame`) → `INTERPOLATOR_CACHE` (stores `Tuple{PCHIPInterp, PCHIPInterp}`)
- `load_scattering_factor_table` removed → use `load_element_interpolators` instead

### Non-breaking improvements

- 4.5× faster precompilation (removed ~38 transitive dependencies)
- 3.6× faster single-material calculations (cached interpolators, fused loops)
- Thread-safe caching (the lock-free read path introduced here was later removed because it raced with `clear_caches!`; see the [Changelog](changelog.md))
- Proper `@threads :dynamic` scheduling for batch calculations
- Energy range validation on single-material path (was only on batch path before)
- Fixed Fluorine data file (`f.nff`) that had a typo preventing BaF₂, CaF₂, etc.

### Compat bounds tightened

`Project.toml` compat bounds changed from `">= X"` (allows any future breaking release) to standard SemVer caret bounds.

## v0.4 → v0.5

### Function renames

| Old (deprecated) | New (recommended) |
|-------------------|-------------------|
| `Refrac(formulas, energies, densities)` | `calculate_xray_properties(formulas, energies, densities)` |
| `SubRefrac(formula, energies, density)` | `calculate_single_material_properties(formula, energies, density)` |

Old names still work but emit deprecation warnings under `--depwarn=yes` (which `Pkg.test()` enables).

### Field renames

| Old (deprecated) | New (recommended) |
|-------------------|-------------------|
| `result.Formula` | `result.formula` |
| `result.MW` | `result.molecular_weight` |
| `result.Number_Of_Electrons` | `result.number_of_electrons` |
| `result.Density` | `result.mass_density` |
| `result.Electron_Density` | `result.electron_density` |
| `result.Energy` | `result.energy` |
| `result.Wavelength` | `result.wavelength` |
| `result.Dispersion` | `result.dispersion` |
| `result.Absorption` | `result.absorption` |
| `result.Critical_Angle` | `result.critical_angle` |
| `result.Attenuation_Length` | `result.attenuation_length` |
| `result.reSLD` | `result.real_sld` |
| `result.imSLD` | `result.imag_sld` |

Old field names still work via a `getproperty` override and will continue to work for the foreseeable future.
