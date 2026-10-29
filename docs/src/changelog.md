# Changelog

All notable changes to this project are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), adhering to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Internal refactor, no behavior change: the molecular-weight and electron-count accumulation shared by `_validate_formulas` and the per-material calculation is one helper (`_molecular_weight_and_electrons`), the 0.03 / 30 keV limits are the named constants `MIN_ENERGY_KEV` / `MAX_ENERGY_KEV`, and an unused intermediate array was removed

## [0.7.2] - 2026-10-29

### Fixed
- Data race on the element caches: the lock-free `CACHES_FROZEN` read path raced with `clear_caches!` and with concurrent cache population, causing crashes or spurious "Invalid formulas" errors under multithreaded use. All cache reads now take the lock; `CACHES_FROZEN` is removed
- `_validate_formulas` no longer reports unexpected internal errors as invalid formulas; only `ArgumentError`s are collected
- Formulas with a zero count (`H0`, `(SiO2)0`, `SiO0`) now raise `ArgumentError` instead of silently dropping the element or returning NaN results
- `XRayResult.energy` is now a private copy; it no longer aliases the caller's vector or other results from the same batch
- `@debug` output for derived quantities now shows values instead of `#NN (generic function ...)`
- Non-finite or unusable element counts now raise `ArgumentError` instead of returning NaN/Inf properties: a group multiplier that overflows to `Inf` (e.g. a 200-digit count inside a 200-digit group multiplier), and counts so small or large that the molecular weight or the density prefactor leaves the representable range (e.g. `H0.` followed by 319 zeros and a `1`, which gave `dispersion = [Inf]`)
- `atomic_number_and_mass` converts only a missing key (`KeyError`) into the "not found in periodic table" `ArgumentError`, and `load_element_interpolators` reports an unknown element by checking for the data file instead of catching every exception. An existing but unreadable or corrupt `.nff` table now raises an `ErrorException` naming the element and file (it was reported as "NOT in the table list", and an `ArgumentError` would be reported by `_validate_formulas` as an invalid formula)
- Deeply nested parentheses (more than 64 levels) now raise `ArgumentError`; a ~100k-deep formula previously overflowed the stack, and since `_validate_formulas` stopped swallowing non-`ArgumentError` exceptions it would have escaped as a `StackOverflowError`
- `result.<field>` access is type-stable again (the legacy-name lookup no longer sits on the hot path); legacy names such as `result.Energy` still work
- Docs: δ/β equations now show the 10⁶ g/cm³→g/m³ factor, the constants table matches the code, f₁/f₂ units are per formula unit, the quick-start no longer indexes batch results by formula string, and `ų` typos are corrected to `Å³`

### Performance
- Evaluate f₁ and f₂ of each element together: the energy interval is located once (a short forward walk over the sorted energies, with a binary-search relocation for the first point, sparse energies and out-of-order input) and the Hermite basis terms are shared. Results are bit-identical to the previous per-point `itp(x)` evaluation (verified for all 92 elements by `test/test_edge_cases.jl`). Single material, 1k energies: 70 → 24 µs; 8-formula batch, 1k energies: 706 → 231 µs (1 thread)
- `precompile` statements for the public entry points cut the first batch call from ~30-50 ms to ~1 ms
- `accumulate_optical_coefficients!` returns immediately for empty input

### Docs
- README, Getting Started and the public docstrings now document the input constraints (finite positive element counts, 64-level nesting limit, finite energies), the `ArgumentError` cases, and that `result.energy` is sorted ascending; f₁/f₂ and electron counts are described per formula unit
- Migration guide no longer claims the cache has lock-free reads

### Tests
- Added a closed-form δ/β check (rtol 1e-9); the previous absolute tolerance (1e-6) let a 4% error in the optical coefficients pass

## [0.7.1] - 2026-03-28

### Fixed
- `deploydocs` repo URL and canonical link now match actual GitHub repo name (`XRayLabTool` not `XRayLabTool.jl`)
- Race condition in `clear_caches!`: frozen flag was cleared before dicts were emptied, allowing lock-free reads on a mutating Dict
- `@debug` keyword args for derived quantities now use lazy closures to avoid eager evaluation when logging is disabled

### Performance
- Lock-free frozen cache reads after batch validation (CACHES_FROZEN atomic flag eliminates ~5μs lock overhead per warm call)
- `issorted` guard avoids redundant sort allocation when energies already sorted (batch path)
- Eliminated redundant `Vector{XRayResult}` copy in batch return path
- Named unit-conversion constants replace inline magic numbers
- Single-material allocations reduced 19% (23.5KB → 19.1KB per call)

### Changed
- Extracted `_compute_derived_quantities` function for readability
- Removed 3 untested deprecated internal wrappers (`get_atomic_data`, `load_f1f2_table`, `calculate_scattering_factors!`)
- Removed redundant `format_project.jl` (pre-commit hook handles formatting)
- Renamed `NS_TO_MS` → `NS_PER_MS` for clarity
- 419 tests (up from 408 in v0.7.0)

## [0.7.0] - 2026-03-28

### Breaking
- `calculate_xray_properties` now returns `Vector{XRayResult}` instead of `Dict{String, XRayResult}`. Use index-based access (`results[1]`) instead of key-based (`results["SiO2"]`). This preserves duplicate formulas with different densities.
- Batch validation is strict: all formulas are checked upfront and invalid ones are listed in a single `ArgumentError` before any computation begins.
- `calculate_single_material_properties` now sorts energies ascending, matching batch API behavior.

### Added
- Recursive descent formula parser supporting parenthesized groups: `Ca(OH)2`, `Ca3(PO4)2`, `Fe((OH)2)3`
- Comprehensive debug logging: 27 `@debug` log points across 6 categories (`:parser`, `:cache`, `:io`, `:validation`, `:computation`, `:batch`). Enable with `ENV["JULIA_DEBUG"] = "XRayLabTool"`.
- ASCII validation in formula parser and element interpolator loader
- Thread-safety test for concurrent cache access
- Negative delta guard: `critical_angle` returns 0.0 instead of throwing `DomainError` when delta is negative near absorption edges

### Fixed
- Thread-safe cache reads: replaced lock-free `Dict.get()` with locked reads (Julia's `Dict` is not safe for concurrent read + write)
- Physical constants updated to CODATA 2018 / 2019 SI exact values (Planck, elementary charge, Avogadro — ~2e-7 relative change)
- Formula parser no longer silently produces wrong results for parenthesized formulas

### Changed
- 408 tests (up from 357): formula parsing, batch validation, thread safety, negative delta, logging, performance
- `Logging` stdlib added to dependencies

## [0.6.0] - 2026-03-28

### Breaking
- Removed `CSV.jl` and `DataFrames.jl` dependencies (replaced with zero-dependency file parsing)
- Internal `F1F2_TABLE_CACHE` replaced with `INTERPOLATOR_CACHE` (caches interpolator pairs directly)
- Removed internal functions: `load_scattering_factor_table`, `pchip_interpolators`, `create_interpolators`

### Performance
- 4.5x faster precompilation (62 to 24 transitive packages)
- 3.6x faster single-material calculations (cached interpolators, fused computation loops)
- 4x fewer allocations per call (19 KB vs 80 KB for 191 energy points)
- Lock-free cache reads on warm-cache path (zero allocation on cache hits)
- `@threads :dynamic` scheduling for better load balancing in batch mode
- Non-allocating energy range validation
- Type-stable hot loop via concrete `ElementInterpolatorData` type alias

### Fixed
- Data race: element caches now pre-populated before `@threads` loop entry
- Fluorine data file (`f.nff`) typo: first energy was 129.3 instead of 29.3
- Energy range validation added to single-material path (previously only batch path validated)
- Circular deprecation architecture: `@deprecate` macros replaced with explicit `depwarn` wrappers

### Added
- `Base.propertynames` override for tab-completion on `XRayResult`
- Documenter.jl documentation site with physics background, API reference, and migration guide
- Consolidated CI pipeline with multi-threaded testing, formatting checks, and documentation deployment
- 357 tests (up from ~80): physical constants, numerical sanity, scaling properties, physics cross-validation, property-based formula parsing, performance regression guards

### Changed
- `Dict`-based `getproperty` with O(1) lookup replaces O(n) if-elseif chain
- Narrowed imports: `using Unitful: ustrip`, `using Base.Threads: @threads`
- SemVer-compatible compat bounds in Project.toml
- Batch results collected via lock-free `Vector` indexed writes instead of `Dict` + `ReentrantLock`

## [0.5.0] - 2025-03-28

### Added
- New descriptive API:
  - `calculate_xray_properties` replaces `Refrac` (multiple materials)
  - `calculate_single_material_properties` replaces `SubRefrac` (single material)
- New snake_case field names on `XRayResult` (see table below)
- Backward-compatible `getproperty` override: old field names still work
- Multi-threaded batch processing via `@threads`
- Comprehensive docstrings on all public and internal functions

### Deprecated
- `Refrac()` — use `calculate_xray_properties()` instead
- `SubRefrac()` — use `calculate_single_material_properties()` instead
- Old field names (still accessible, will continue to work):

| Old (deprecated) | New (recommended) |
|---|---|
| `Formula` | `formula` |
| `MW` | `molecular_weight` |
| `Number_Of_Electrons` | `number_of_electrons` |
| `Density` | `mass_density` |
| `Electron_Density` | `electron_density` |
| `Energy` | `energy` |
| `Wavelength` | `wavelength` |
| `Dispersion` | `dispersion` |
| `Absorption` | `absorption` |
| `Critical_Angle` | `critical_angle` |
| `Attenuation_Length` | `attenuation_length` |
| `reSLD` | `real_sld` |
| `imSLD` | `imag_sld` |

Note: `f1` and `f2` are unchanged.

## [0.4.0]

### Added
- Initial public release
- Core X-ray property calculations based on CXRO/Henke tables
- Support for all elements H through U, 30 eV to 30 keV
- Caching system for atomic data and scattering factor tables
- `Refrac` (batch) and `SubRefrac` (single material) API
