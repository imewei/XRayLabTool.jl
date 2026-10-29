"""
XRayLabTool - Material property calculations for X-ray interactions

This module provides functions to calculate X-ray optical properties of materials
based on their chemical composition and density.

**Main Functions (New API - Recommended):**
- `calculate_xray_properties`: Calculate properties for multiple chemical formulas
- `calculate_single_material_properties`: Calculate properties for a single chemical formula

**Usage Examples:**
```julia
# Single formula (new API)
result = calculate_single_material_properties("SiO2", [8.0, 10.0, 12.0], 2.2)

# Multiple formulas (new API) — returns Vector{XRayResult}
results = calculate_xray_properties(["SiO2", "Al2O3"], [8.0, 10.0, 12.0], [2.2, 3.95])
sio2 = results[1]  # access by index (same order as input)

# Access properties
println("Molecular weight: ", result.molecular_weight)
println("Critical angle: ", result.critical_angle[1])
```

!!! warning "Deprecated Functions"
    The old function names `Refrac` and `SubRefrac` are still available for backward
    compatibility but are deprecated. Please use the new API names:
    - `Refrac` → `calculate_xray_properties`
    - `SubRefrac` → `calculate_single_material_properties`

**Input Parameters:**
1. Chemical formula(s): Case-sensitive strings (e.g., "CO" for Carbon Monoxide vs "Co" for Cobalt)
2. Energy: Vector of X-ray energies in KeV (range: 0.03 - 30 KeV)
3. Mass density: Density in g/cm³ (single value or vector matching formula count)

**Output Properties:**
The functions return XRayResult struct(s) containing:
1. Chemical formula (`formula`)
2. Molecular weight in g/mol (`molecular_weight`)
3. Number of electrons per molecule (`number_of_electrons`)
4. Mass density in g/cm³ (`mass_density`)
5. Electron density in 1/Å³ (`electron_density`)
6. X-ray energy in KeV (`energy`)
7. X-ray wavelength in Å (`wavelength`)
8. Dispersion coefficient (`dispersion`)
9. Absorption coefficient (`absorption`)
10. Real part of the scattering factor, summed over the formula unit (`f1`)
11. Imaginary part of the scattering factor, summed over the formula unit (`f2`)
12. Critical angle in degrees (`critical_angle`)
13. Attenuation length in cm (`attenuation_length`)
14. Real part of scattering length density in Å⁻² (`real_sld`)
15. Imaginary part of SLD in Å⁻² (`imag_sld`)
"""
module XRayLabTool

using Logging
using PCHIPInterpolation
using Mendeleev: chem_elements
using Unitful: ustrip
using Base.Threads: @threads

# Export statements - both new and deprecated names for backward compatibility
export Refrac, SubRefrac, XRayResult
export calculate_xray_properties, calculate_single_material_properties

# Concrete type alias for the PCHIP interpolator used throughout
const PCHIPInterp = Interpolator{Vector{Float64}, Vector{Float64}, Vector{Float64}}

# Concrete element data tuple: (count, f1_interpolator, f2_interpolator)
const ElementInterpolatorData = Tuple{Float64, PCHIPInterp, PCHIPInterp}

# Deepest parenthesis nesting accepted by the formula parser; real formulas nest a few levels.
# Bounds the recursion so a hostile string raises ArgumentError instead of overflowing the stack.
const MAX_FORMULA_NESTING = 64

# =====================================================================================
# DATA STRUCTURES
# =====================================================================================

"""
    XRayResult

Struct to store complete X-ray optical property calculations for a material.

Contains all computed properties including scattering factors, optical constants,
and derived quantities like critical angles and attenuation lengths.
"""
struct XRayResult
    formula::String                      # Chemical formula
    molecular_weight::Float64            # Molecular weight (g/mol)
    number_of_electrons::Float64         # Electrons per molecule
    mass_density::Float64                # Mass density (g/cm³)
    electron_density::Float64            # Electron density (1/Å³)
    energy::Vector{Float64}              # X-ray energy (KeV)
    wavelength::Vector{Float64}          # X-ray wavelength (Å)
    dispersion::Vector{Float64}          # Dispersion coefficient
    absorption::Vector{Float64}          # Absorption coefficient
    f1::Vector{Float64}                  # Real part of scattering factor (formula-unit sum)
    f2::Vector{Float64}                  # Imaginary part of scattering factor (formula-unit sum)
    critical_angle::Vector{Float64}      # Critical angle (degrees)
    attenuation_length::Vector{Float64}  # Attenuation length (cm)
    real_sld::Vector{Float64}            # Real part of SLD (Å⁻²)
    imag_sld::Vector{Float64}            # Imaginary part of SLD (Å⁻²)
end

# Deprecated field name mappings (old name → new name) for backward compatibility
const _LEGACY_FIELD_MAP = Dict{Symbol, Symbol}(
    :Formula => :formula,
    :MW => :molecular_weight,
    :Number_Of_Electrons => :number_of_electrons,
    :Density => :mass_density,
    :Electron_Density => :electron_density,
    :Energy => :energy,
    :Wavelength => :wavelength,
    :Dispersion => :dispersion,
    :Absorption => :absorption,
    :Critical_Angle => :critical_angle,
    :Attenuation_Length => :attenuation_length,
    :reSLD => :real_sld,
    :imSLD => :imag_sld,
)

# Real fields take the `hasfield` branch, which constant-folds for literal names and keeps
# `result.field` type-stable; only legacy names consult the map.
function Base.getproperty(result::XRayResult, name::Symbol)
    hasfield(XRayResult, name) && return getfield(result, name)
    return getfield(result, get(_LEGACY_FIELD_MAP, name, name))
end

function Base.propertynames(::XRayResult, private::Bool = false)
    return (fieldnames(XRayResult)..., keys(_LEGACY_FIELD_MAP)...)
end

# =====================================================================================
# PHYSICAL CONSTANTS
# =====================================================================================

"""
Physical constants used in X-ray calculations.
All values are in SI units unless otherwise specified.
"""
const THOMSON_SCATTERING_LENGTH = 2.8179403205e-15     # Thomson scattering length (m) — CODATA 2018
const SPEED_OF_LIGHT = 299792458.0    # Speed of light (m/s) — exact by definition
const PLANCK = 6.62607015e-34         # Planck constant (J·s) — exact since 2019 SI redefinition
const ELEMENT_CHARGE = 1.602176634e-19 # Elementary charge (C) — exact since 2019 SI redefinition
const AVOGADRO = 6.02214076e23        # Avogadro's number (mol⁻¹) — exact since 2019 SI redefinition

# Pre-computed constants for efficiency
const HC_OVER_ELECTRON_CHARGE_keV = (SPEED_OF_LIGHT * PLANCK / ELEMENT_CHARGE) / 1000.0
const SCATTERING_PREFACTOR = THOMSON_SCATTERING_LENGTH * AVOGADRO * 1e6 / (2 * π)

# Unit conversion factors
const KEV_TO_EV = 1000.0
const M_TO_ANGSTROM = 1e10
const M_TO_CM = 1e2
const NS_PER_MS = 1e6
const ELECTRON_DENSITY_FACTOR = AVOGADRO * 1e6 / 1e30  # Nₐ × (cm³→m³) / (m³→Å³): converts (ρ/M)·Z to electrons/Å³
const SLD_FACTOR = 2 * π / 1e20                          # Å⁻² conversion

# Deprecation bridges for backward compatibility
const THOMPSON = THOMSON_SCATTERING_LENGTH  # Deprecated: use THOMSON_SCATTERING_LENGTH
const ENERGY_TO_WAVELENGTH_FACTOR = HC_OVER_ELECTRON_CHARGE_keV  # Deprecated: use HC_OVER_ELECTRON_CHARGE_keV
const SCATTERING_FACTOR = SCATTERING_PREFACTOR  # Deprecated: use SCATTERING_PREFACTOR

# =====================================================================================
# CACHING SYSTEM
# =====================================================================================

"""
Caches for performance optimization:
- ATOMIC_DATA_CACHE: Stores (atomic_number, atomic_mass) for elements
- INTERPOLATOR_CACHE: Stores pre-built PCHIP interpolator pairs (f1, f2) per element

Both caches are protected by locks for thread safety during parallel processing.
Every read and write goes through the lock: a lock-free read would race with
`clear_caches!` and with concurrent cache population from other tasks.
"""
const ATOMIC_DATA_CACHE = Dict{String, Tuple{Int, Float64}}()
const ATOMIC_DATA_LOCK = ReentrantLock()
const INTERPOLATOR_CACHE = Dict{String, Tuple{PCHIPInterp, PCHIPInterp}}()
const INTERPOLATOR_LOCK = ReentrantLock()

# =====================================================================================
# VALIDATION
# =====================================================================================

"""Non-allocating energy range validation (avoids temporary BitVector from broadcasting)."""
@inline function _validate_energy_range(energies_keV::Vector{Float64})
    if isempty(energies_keV)
        throw(ArgumentError("Energy vector must not be empty"))
    end
    if any(e -> !isfinite(e) || e < 0.03 || e > 30.0, energies_keV)
        throw(ArgumentError("Energy is out of range 0.03KeV ~ 30KeV (NaN/Inf not allowed)"))
    end
    @debug "Energy range validated" _group = :validation count = length(energies_keV) min_keV =
        minimum(energies_keV) max_keV = maximum(energies_keV)
end

"""Validate that mass density is a positive finite number."""
@inline function _validate_density(density::Float64, label::String = "")
    if !isfinite(density) || density <= 0
        ctx = isempty(label) ? "" : " for $label"
        throw(ArgumentError("Mass density must be positive and finite$ctx (got $density)"))
    end
end

"""Validate that the molecular weight implied by the element counts is usable.

Counts that are subnormal or huge push the molecular weight (or the density prefactor
derived from it) out of the representable range and would otherwise yield Inf results.
"""
@inline function _validate_molecular_weight(molecular_weight::Float64, formula::String)
    if !(
        isfinite(molecular_weight) &&
        molecular_weight > 0 &&
        isfinite(SCATTERING_PREFACTOR / molecular_weight)
    )
        throw(
            ArgumentError(
                "Invalid chemical formula: element counts in \"$formula\" give a molecular weight ($molecular_weight g/mol) outside the representable range",
            ),
        )
    end
end

"""Validate all formulas upfront. Throws ArgumentError listing all invalid formulas."""
function _validate_formulas(formulaList::Vector{String})
    invalid = String[]
    for formula in formulaList
        try
            symbols, counts = parse_formula(formula)
            molecular_weight = 0.0
            for (sym, count) in zip(symbols, counts)
                molecular_weight += count * atomic_number_and_mass(sym)[2]
                load_element_interpolators(sym)
            end
            _validate_molecular_weight(molecular_weight, formula)
        catch e
            # Only expected failures mean "invalid formula"; surface real internal errors
            e isa ArgumentError || rethrow()
            push!(invalid, formula)
        end
    end
    if !isempty(invalid)
        @debug "Formula validation failed" _group = :validation invalid_formulas =
            join(invalid, ", ") invalid_count = length(invalid)
        throw(
            ArgumentError(
                "Invalid formulas: $(join(invalid, ", ")). " *
                "Check element symbols and spelling.",
            ),
        )
    end
    @debug "All formulas validated" _group = :validation formula_count = length(formulaList)
end

# =====================================================================================
# HELPER FUNCTIONS
# =====================================================================================

"""
    atomic_number_and_mass(element_symbol::String) -> (Int, Float64)

Retrieve atomic number and atomic mass for a given element symbol.
Uses caching to avoid repeated lookups of the same element.

# Arguments
- 'element_symbol::String': Chemical symbol (e.g., "Si", "O")

# Returns
- 'Tuple{Int, Float64}': (atomic_number, atomic_mass_in_amu)

# Throws
- 'ArgumentError': If element symbol is not found in periodic table
"""
function atomic_number_and_mass(element_symbol::String)
    cached = lock(ATOMIC_DATA_LOCK) do
        get(ATOMIC_DATA_CACHE, element_symbol, nothing)
    end
    if cached !== nothing
        @debug "Atomic data cache hit" _group = :cache element = element_symbol Z =
            cached[1] mass = cached[2]
        return cached
    end

    # Search periodic table for element; only a missing key means "unknown element"
    elem = try
        chem_elements[Symbol(element_symbol)]
    catch e
        e isa KeyError || rethrow()
        throw(ArgumentError("Element $element_symbol not found in periodic table"))
    end
    atomic_number = Int(elem.atomic_number)
    atomic_mass = Float64(ustrip(elem.atomic_weight))

    result = (atomic_number, atomic_mass)

    # Locked write — rare path (first lookup per element)
    lock(ATOMIC_DATA_LOCK) do
        ATOMIC_DATA_CACHE[element_symbol] = result
    end
    @debug "Atomic data cache miss → stored" _group = :cache element = element_symbol Z =
        atomic_number mass = atomic_mass

    return result
end

"""
    parse_formula(formulaStr::String) -> (Vector{String}, Vector{Float64})

Parse a chemical formula string into element symbols and their counts.

# Arguments
- 'formulaStr::String': Chemical formula (e.g., "SiO2", "Al2O3", "Ca(OH)2", "Ca3(PO4)2")

# Returns
- 'Tuple{Vector{String}, Vector{Float64}}': (element_symbols, element_counts)

# Examples
```julia
symbols, counts = parse_formula("SiO2")
# symbols = ["Si", "O"], counts = [1.0, 2.0]

symbols, counts = parse_formula("Al2O3")
# symbols = ["Al", "O"], counts = [2.0, 3.0]

symbols, counts = parse_formula("Ca(OH)2")
# symbols = ["Ca", "O", "H"], counts = [1.0, 2.0, 2.0]

symbols, counts = parse_formula("Ca3(PO4)2")
# symbols = ["Ca", "P", "O"], counts = [3.0, 2.0, 8.0]
```

# Throws
- 'ArgumentError': If the formula string is invalid or empty, an element count is not positive and finite, or parentheses are nested deeper than `MAX_FORMULA_NESTING`
"""
function parse_formula(formulaStr::String)
    @debug "Parsing formula" _group = :parser input = formulaStr
    if isempty(formulaStr)
        @debug "Formula rejected: empty string" _group = :parser
        throw(ArgumentError("Invalid chemical formula: empty string"))
    end

    elements, counts, pos = _parse_group(formulaStr, 1)

    if pos <= lastindex(formulaStr)
        @debug "Formula rejected: unexpected character" _group = :parser formula =
            formulaStr position = pos character = formulaStr[pos]
        throw(
            ArgumentError(
                "Invalid chemical formula: unexpected character '$(formulaStr[pos])' at position $pos in \"$formulaStr\"",
            ),
        )
    end

    if isempty(elements)
        @debug "Formula rejected: no elements found" _group = :parser formula = formulaStr
        throw(ArgumentError("Invalid chemical formula: $formulaStr"))
    end

    if any(c -> !(isfinite(c) && c > 0), counts)
        throw(
            ArgumentError(
                "Invalid chemical formula: element counts must be positive and finite in \"$formulaStr\"",
            ),
        )
    end

    @debug "Formula parsed" _group = :parser formula = formulaStr elements =
        join(elements, ",") counts = join(counts, ",")
    return elements, counts
end

"""
    _parse_group(s, pos) -> (Vector{String}, Vector{Float64}, Int)

Recursive descent parser for chemical formula groups. Handles nested parentheses.
Returns (element_symbols, element_counts, next_position).
"""
function _parse_group(s::String, pos::Int, depth::Int = 0)
    if depth > MAX_FORMULA_NESTING
        throw(
            ArgumentError(
                "Invalid chemical formula: parentheses nested deeper than $MAX_FORMULA_NESTING levels",
            ),
        )
    end
    elements = String[]
    counts = Float64[]
    sizehint!(elements, 4)
    sizehint!(counts, 4)

    while pos <= lastindex(s)
        c = s[pos]

        if c == '('
            # Recurse into sub-group
            sub_elems, sub_counts, pos = _parse_group(s, pos + 1, depth + 1)

            if pos > lastindex(s) || s[pos] != ')'
                throw(
                    ArgumentError(
                        "Invalid chemical formula: unbalanced parenthesis in \"$s\"",
                    ),
                )
            end
            pos += 1  # skip ')'

            if isempty(sub_elems)
                throw(
                    ArgumentError("Invalid chemical formula: empty parentheses in \"$s\""),
                )
            end

            # Read optional multiplier after ')'
            multiplier, pos = _parse_number(s, pos)

            # Apply multiplier to all elements in the sub-group
            for i in eachindex(sub_counts)
                sub_counts[i] *= multiplier
            end
            @debug "Parenthesized group expanded" _group = :parser sub_elements =
                join(sub_elems, ",") multiplier = multiplier
            append!(elements, sub_elems)
            append!(counts, sub_counts)

        elseif c == ')'
            # End of current group — return to caller
            return elements, counts, pos

        elseif isuppercase(c) && isascii(c)
            # Parse element symbol: uppercase ASCII + optional lowercase ASCII
            sym_start = pos
            pos += 1
            while pos <= lastindex(s) && islowercase(s[pos]) && isascii(s[pos])
                pos += 1
            end
            symbol = s[sym_start:prevind(s, pos)]

            # Parse optional count after element
            count, pos = _parse_number(s, pos)

            push!(elements, symbol)
            push!(counts, count)

        else
            # Unexpected character (digits without element, special chars, etc.)
            throw(
                ArgumentError(
                    "Invalid chemical formula: unexpected character '$(c)' at position $pos in \"$s\"",
                ),
            )
        end
    end

    return elements, counts, pos
end

"""
    _parse_number(s, pos) -> (Float64, Int)

Parse an optional number (integer or decimal) at position `pos`.
Returns 1.0 if no number is found.
"""
function _parse_number(s::String, pos::Int)
    start = pos
    while pos <= lastindex(s) && (isdigit(s[pos]) || s[pos] == '.')
        pos += 1
    end
    if pos == start
        return 1.0, pos
    end
    return parse(Float64, s[start:prevind(s, pos)]), pos
end

"""
    _read_nff_interpolators(file, element_symbol) -> (f1_interpolator, f2_interpolator, n_points)

Read an `.nff` table (header line, then `energy_eV,f1,f2` rows) and build the PCHIP
interpolators. A failure here means the bundled data file is unreadable or corrupt, not that
the user supplied a bad formula, so it raises `ErrorException` (never `ArgumentError`, which
`_validate_formulas` would report as an invalid formula).
"""
function _read_nff_interpolators(file::String, element_symbol::String)
    try
        lines = readlines(file)
        n = length(lines) - 1  # Skip header
        energy = Vector{Float64}(undef, n)
        f1_data = Vector{Float64}(undef, n)
        f2_data = Vector{Float64}(undef, n)
        for j in 1:n
            parts = split(lines[j + 1], ',')
            energy[j] = parse(Float64, parts[1])
            f1_data[j] = parse(Float64, parts[2])
            f2_data[j] = parse(Float64, parts[3])
        end
        return Interpolator(energy, f1_data), Interpolator(energy, f2_data), n
    catch e
        error(
            "Scattering-factor data for element $element_symbol is unreadable or corrupt ($file): " *
            sprint(showerror, e),
        )
    end
end

"""
    load_element_interpolators(element_symbol::String) -> (PCHIPInterp, PCHIPInterp)

Load and cache PCHIP interpolators for an element's scattering factors.

Returns cached interpolator pair if available, otherwise loads the element's
.nff data file, builds PCHIP interpolators for f1 and f2, and caches them.

# Arguments
- `element_symbol::String`: Chemical symbol (e.g., "Si", "O")

# Returns
- `Tuple{PCHIPInterp, PCHIPInterp}`: (f1_interpolator, f2_interpolator)

# Throws
- `ArgumentError`: If the element has no data file (unknown element)
- `ErrorException`: If the data file exists but cannot be read or parsed
"""
function load_element_interpolators(element_symbol::String)
    cached = lock(INTERPOLATOR_LOCK) do
        get(INTERPOLATOR_CACHE, element_symbol, nothing)
    end
    if cached !== nothing
        @debug "Interpolator cache hit" _group = :cache element = element_symbol
        return cached
    end

    # Reject non-ASCII symbols before constructing file paths
    if !all(isascii, element_symbol)
        throw(ArgumentError("Element $element_symbol contains non-ASCII characters"))
    end

    # Construct filename (lowercase element symbol + .nff extension)
    fname = lowercase(element_symbol) * ".nff"
    file = normpath(joinpath(@__DIR__, "AtomicScatteringFactor", fname))
    @debug "Resolving .nff file" _group = :io element = element_symbol path = file

    # A missing data file is the expected "unknown element" case
    if !isfile(file)
        throw(
            ArgumentError(
                "Element $element_symbol is NOT in the table list (no file $fname)",
            ),
        )
    end

    # Load data file and build interpolators (outside lock)
    t_io = time_ns()
    itp_f1, itp_f2, n = _read_nff_interpolators(file, element_symbol)
    elapsed_io_ms = (time_ns() - t_io) / NS_PER_MS
    @debug "Loaded .nff file" _group = :io element = element_symbol data_points = n elapsed_ms =
        round(elapsed_io_ms; digits = 3)

    result = (itp_f1, itp_f2)

    # Cache under lock
    lock(INTERPOLATOR_LOCK) do
        INTERPOLATOR_CACHE[element_symbol] = result
    end
    @debug "Interpolator cache miss → stored" _group = :cache element = element_symbol

    return result
end

"""
    _hermite_pair(itp1, itp2, x, j) -> (f1, f2)

Evaluate two PCHIP interpolators that share one grid at `x`, which lies in interval `j`
(`xs[j] <= x <= xs[j+1]`). Same expressions, in the same order, as PCHIPInterpolation's own
evaluation, so the result is bit-identical to `(itp1(x), itp2(x))`; the basis terms are
computed once instead of twice.
"""
@inline function _hermite_pair(itp1::PCHIPInterp, itp2::PCHIPInterp, x::Float64, j::Int)
    xs = itp1.xs
    x1 = xs[j]
    x2 = xs[j + 1]
    h = x2 - x1
    t1 = (x2 - x) / h
    t2 = (x - x1) / h
    ϕ1 = 3t1^2 - 2t1^3
    ϕ2 = 3t2^2 - 2t2^3
    ψ1 = t1^3 - t1^2
    ψ2 = t2^3 - t2^2

    ys1 = itp1.ys
    ds1 = itp1.ds
    f1 = ys1[j] * ϕ1 + ys1[j + 1] * ϕ2 - ds1[j] * h * ψ1 + ds1[j + 1] * h * ψ2

    ys2 = itp2.ys
    ds2 = itp2.ds
    f2 = ys2[j] * ϕ1 + ys2[j + 1] * ϕ2 - ds2[j] * h * ψ1 + ds2[j + 1] * h * ψ2
    return f1, f2
end

"""
    accumulate_optical_coefficients!(...)

Vectorized calculation of X-ray scattering factors and optical properties.

This function performs the core calculation of delta (δ), beta (β), and total
scattering factors for a material based on its elemental composition.

# Arguments
- 'energies_eV::Vector{Float64}': X-ray energies in eV
- 'wavelengths_m::Vector{Float64}': Corresponding wavelengths in meters
- 'mass_density::Float64': Material density in g/cm³
- 'molecular_weight::Float64': Molecular weight in g/mol
- 'element_data::Vector{ElementInterpolatorData}': Element data (count, f1_interp, f2_interp)
- 'delta::Vector{Float64}': Output array for δ (real part of refractive index)
- 'beta::Vector{Float64}': Output array for β (imaginary part of refractive index)
- 'f1_total::Vector{Float64}': Output array for total f1 values
- 'f2_total::Vector{Float64}': Output array for total f2 values

# Mathematical Background
The optical coefficients delta and beta are calculated using:
- δ = (λ²/2π) × rₑ × (10⁶·ρ) × Nₐ × (Σᵢ nᵢ × f1ᵢ) / M
- β = (λ²/2π) × rₑ × (10⁶·ρ) × Nₐ × (Σᵢ nᵢ × f2ᵢ) / M

Where:
- λ: X-ray wavelength in meters
- rₑ: Thomson scattering length
- ρ: Mass density in g/cm³ (the factor 10⁶ converts to g/m³ for λ in metres)
- Nₐ: Avogadro's number
- nᵢ: Number of atoms of element i
- f1ᵢ, f2ᵢ: Atomic scattering factors for element i
- M: Molecular weight
"""
function accumulate_optical_coefficients!(
    energies_eV::Vector{Float64},
    wavelengths_m::Vector{Float64},
    mass_density::Float64,
    molecular_weight::Float64,
    element_data::Vector{ElementInterpolatorData},
    delta::Vector{Float64},
    beta::Vector{Float64},
    f1_total::Vector{Float64},
    f2_total::Vector{Float64},
)
    n_energies = length(energies_eV)
    n_energies == 0 && return nothing

    # Pre-compute density-dependent factor for efficiency
    # Factor includes: (λ²/2π) × rₑ × ρ × Nₐ / M
    common_factor = SCATTERING_PREFACTOR * mass_density / molecular_weight

    # Process each element in the formula
    for (count, itp1, itp2) in element_data
        # Element-specific contribution factor
        element_contribution_factor = common_factor * count

        # f1 and f2 of an element share one energy grid, so the interval is located once per
        # energy; energies arrive sorted, so a monotone walk replaces the binary search.
        xs = itp1.xs
        shared_grid = xs === itp2.xs
        x_lo = first(xs)
        x_hi = last(xs)
        last_interval = lastindex(xs) - 1
        j = firstindex(xs)

        @inbounds for i in 1:n_energies
            x = energies_eV[i]
            if shared_grid && x_lo <= x <= x_hi
                # Binary-search when moving backwards (unsorted input) or more than 8 nodes
                # ahead (sparse energies, or the first point); otherwise the walk below takes
                # at most 8 steps. Both routes select the same interval.
                if x < xs[j] || x >= xs[min(j + 8, last_interval + 1)]
                    j = clamp(searchsortedlast(xs, x), firstindex(xs), last_interval)
                end
                while j < last_interval && x >= xs[j + 1]
                    j += 1
                end
                f1_val, f2_val = _hermite_pair(itp1, itp2, x, j)
            else
                # Outside the tabulated range (or distinct grids): defer to the interpolators
                f1_val = itp1(x)
                f2_val = itp2(x)
            end

            # Calculate wavelength-dependent factors
            wave_sq = wavelengths_m[i]^2

            # Accumulate contributions to optical properties
            delta[i] += wave_sq * element_contribution_factor * f1_val
            beta[i] += wave_sq * element_contribution_factor * f2_val

            # Accumulate total scattering factors
            f1_total[i] += count * f1_val
            f2_total[i] += count * f2_val
        end
    end
end

# =====================================================================================
# DERIVED QUANTITIES
# =====================================================================================

"""
    _compute_derived_quantities(wavelengths_m, delta, beta)

Compute wavelength (Å), critical angle (deg), attenuation length (cm),
and SLD (Å⁻²) from optical coefficients delta and beta.
"""
function _compute_derived_quantities(
    wavelengths_m::Vector{Float64},
    delta::Vector{Float64},
    beta::Vector{Float64},
)
    n = length(wavelengths_m)
    wavelength_angstrom = Vector{Float64}(undef, n)
    critical_angle = Vector{Float64}(undef, n)
    attenuation_length = Vector{Float64}(undef, n)
    re_sld = Vector{Float64}(undef, n)
    im_sld = Vector{Float64}(undef, n)

    inv_4pi = 1.0 / (4 * π)
    rad_to_deg = 180.0 / π

    @inbounds for i in 1:n
        λ = wavelengths_m[i]
        λ_sq = λ * λ

        wavelength_angstrom[i] = λ * M_TO_ANGSTROM
        critical_angle[i] = delta[i] > 0 ? sqrt(2.0 * delta[i]) * rad_to_deg : 0.0
        attenuation_length[i] = λ * inv_4pi / beta[i] * M_TO_CM
        re_sld[i] = delta[i] * SLD_FACTOR / λ_sq
        im_sld[i] = beta[i] * SLD_FACTOR / λ_sq
    end

    return wavelength_angstrom, critical_angle, attenuation_length, re_sld, im_sld
end

# =====================================================================================
# MAIN CALCULATION FUNCTIONS
# =====================================================================================

"""
    _calculate_xray_properties_impl(formulaList, energies_keV, massDensityList) -> Vector{XRayResult}

Internal implementation for batch X-ray property calculation.
"""
function _calculate_xray_properties_impl(
    formulaList::Vector{String},
    energies_keV::Vector{Float64},
    massDensityList::Vector{Float64},
)
    t_batch_start = time_ns()
    @debug "Batch calculation started" _group = :batch formula_count = length(formulaList) energy_count =
        length(energies_keV) threads = Threads.nthreads()
    @debug "Batch input summary" _group = :validation formula_count = length(formulaList) energy_count =
        length(energies_keV) density_count = length(massDensityList)

    # ==================================================================================
    # INPUT VALIDATION
    # ==================================================================================

    # Check for empty inputs
    if isempty(formulaList) || isempty(energies_keV)
        throw(ArgumentError("Formula list and energy vector must not be empty"))
    end

    # Non-allocating energy range validation (functional form avoids temporary BitVector)
    _validate_energy_range(energies_keV)

    # Check length consistency
    if length(formulaList) != length(massDensityList)
        throw(ArgumentError("Formula list and mass density list must have the same length"))
    end

    # Validate all densities
    for (i, d) in enumerate(massDensityList)
        _validate_density(d, formulaList[i])
    end

    # ==================================================================================
    # STRICT VALIDATION (single-threaded — validates ALL formulas before any computation,
    # reports all errors at once, and warms caches as a side effect)
    # ==================================================================================

    energies_keV_sorted = sort(energies_keV)

    t_validate = time_ns()
    _validate_formulas(formulaList)
    elapsed_validate_ms = (time_ns() - t_validate) / NS_PER_MS
    @debug "Validation + cache warm-up complete" _group = :batch elapsed_ms =
        round(elapsed_validate_ms; digits = 3)

    # ==================================================================================
    # PARALLEL CALCULATION (each task writes its own results_vec slot; cache reads take a short lock)
    # ==================================================================================

    n_formulas = length(formulaList)
    results_vec = Vector{XRayResult}(undef, n_formulas)

    t_parallel = time_ns()
    @threads :dynamic for i in 1:n_formulas
        results_vec[i] = _calculate_single_material_impl(
            formulaList[i],
            energies_keV_sorted,
            massDensityList[i];
            _validated = true,
        )
    end
    elapsed_parallel_ms = (time_ns() - t_parallel) / NS_PER_MS
    @debug "Parallel computation complete" _group = :batch elapsed_ms =
        round(elapsed_parallel_ms; digits = 3)

    elapsed_total_ms = (time_ns() - t_batch_start) / NS_PER_MS
    @debug "Batch calculation complete" _group = :batch result_count = n_formulas elapsed_ms =
        round(elapsed_total_ms; digits = 3)

    return results_vec
end

"""
    _calculate_single_material_impl(formulaStr, energies_keV, massDensity) -> XRayResult

Internal implementation for single-material X-ray property calculation.
"""
function _calculate_single_material_impl(
    formulaStr::String,
    energies_keV::Vector{Float64},
    massDensity::Float64;
    _validated::Bool = false,
)
    # ==================================================================================
    # INPUT VALIDATION (skipped when called from batch path which validates upfront)
    # ==================================================================================

    if !_validated
        _validate_energy_range(energies_keV)
        _validate_density(massDensity, formulaStr)
    end

    t_start = time_ns()

    # Sorted private copy: the result must not alias the caller's vector (or sibling results)
    energies_keV = issorted(energies_keV) ? copy(energies_keV) : sort(energies_keV)

    @debug "Single-material calculation started" _group = :computation formula = formulaStr density =
        massDensity energy_count = length(energies_keV)

    # ==================================================================================
    # FORMULA PARSING AND ATOMIC DATA LOOKUP
    # ==================================================================================

    # Parse the chemical formula into elements and their counts
    element_symbols, element_counts = parse_formula(formulaStr)
    n_elements = length(element_symbols)
    n_energies = length(energies_keV)

    # Pre-allocate arrays for efficiency
    molecular_weight = 0.0
    number_of_electrons = 0.0
    atomic_data = Vector{Tuple{Int, Float64}}(undef, n_elements)

    # Look up atomic data for each element in the formula
    for i in 1:n_elements
        atomic_number, atomic_mass = atomic_number_and_mass(element_symbols[i])
        atomic_data[i] = (atomic_number, atomic_mass)

        # Accumulate molecular weight and total electrons
        molecular_weight += element_counts[i] * atomic_mass
        number_of_electrons += atomic_number * element_counts[i]
    end

    _validate_molecular_weight(molecular_weight, formulaStr)

    @debug "Molecular properties computed" _group = :computation formula = formulaStr MW =
        round(molecular_weight; digits = 3) electrons = number_of_electrons element_count =
        n_elements

    # ==================================================================================
    # ENERGY-WAVELENGTH CONVERSION + ARRAY INITIALIZATION (fused to avoid intermediates)
    # ==================================================================================

    wavelengths_m = Vector{Float64}(undef, n_energies)
    energies_eV = Vector{Float64}(undef, n_energies)
    @inbounds for i in 1:n_energies
        wavelengths_m[i] = HC_OVER_ELECTRON_CHARGE_keV / energies_keV[i]
        energies_eV[i] = energies_keV[i] * KEV_TO_EV
    end

    # Pre-allocate output arrays with zeros
    delta = zeros(Float64, n_energies)
    beta = zeros(Float64, n_energies)
    f1_total = zeros(Float64, n_energies)
    f2_total = zeros(Float64, n_energies)

    # ==================================================================================
    # SCATTERING FACTOR INTERPOLATOR LOADING (cached per element)
    # ==================================================================================

    element_data = Vector{ElementInterpolatorData}(undef, n_elements)

    for i in 1:n_elements
        itp_f1, itp_f2 = load_element_interpolators(element_symbols[i])
        element_data[i] = (element_counts[i], itp_f1, itp_f2)
        @debug "Element scattering data loaded" _group = :computation element =
            element_symbols[i] count = element_counts[i]
    end

    # ==================================================================================
    # MAIN SCATTERING CALCULATION
    # ==================================================================================

    # Calculate dispersion, absorption, and total scattering factors
    # This is the computationally intensive part
    accumulate_optical_coefficients!(
        energies_eV,
        wavelengths_m,
        massDensity,
        molecular_weight,
        element_data,
        delta,
        beta,
        f1_total,
        f2_total,
    )

    # ==================================================================================
    # DERIVED QUANTITY CALCULATIONS
    # ==================================================================================

    electron_density =
        massDensity / molecular_weight * ELECTRON_DENSITY_FACTOR * number_of_electrons

    wavelength_angstrom, critical_angle, attenuation_length, re_sld, im_sld =
        _compute_derived_quantities(wavelengths_m, delta, beta)

    @debug "Derived quantities computed" _group = :computation formula = formulaStr electron_density =
        round(electron_density; sigdigits = 6) critical_angle_range = (
        round(minimum(critical_angle); sigdigits = 4),
        round(maximum(critical_angle); sigdigits = 4),
    ) attenuation_length_range = (
        round(minimum(attenuation_length); sigdigits = 4),
        round(maximum(attenuation_length); sigdigits = 4),
    )

    @debug "Result summary" _group = :computation formula = formulaStr delta_range =
        (minimum(delta), maximum(delta)) beta_range = (minimum(beta), maximum(beta)) sld_range =
        (minimum(re_sld), maximum(re_sld))

    elapsed_ms = (time_ns() - t_start) / NS_PER_MS
    @debug "Material calculation complete" _group = :computation formula = formulaStr MW =
        round(molecular_weight; digits = 3) energy_count = n_energies elapsed_ms =
        round(elapsed_ms; digits = 3)

    # ==================================================================================
    # RESULT ASSEMBLY
    # ==================================================================================

    return XRayResult(
        formulaStr,
        molecular_weight,
        number_of_electrons,
        massDensity,
        electron_density,
        energies_keV,
        wavelength_angstrom,
        delta,
        beta,
        f1_total,
        f2_total,
        critical_angle,
        attenuation_length,
        re_sld,
        im_sld,
    )
end

# =====================================================================================
# NEW API FUNCTIONS WITH DESCRIPTIVE NAMES
# =====================================================================================

"""
    calculate_xray_properties(formulaList, energy, massDensityList) -> Vector{XRayResult}

Calculate X-ray optical properties for multiple chemical formulas (vector interface).

# Arguments
- `formulaList::Vector{String}`: List of chemical formulas
- `energy::Vector{Float64}`: X-ray energies in KeV (0.03 - 30 KeV)
- `massDensityList::Vector{Float64}`: Mass densities in g/cm³

# Returns
- `Vector{XRayResult}`: Results in same order as input formulas

# Examples
```julia
formulas = ["SiO2", "Al2O3", "Fe2O3"]
energies = [8.0, 10.0, 12.0, 15.0]
densities = [2.2, 3.95, 5.24]

results = calculate_xray_properties(formulas, energies, densities)
sio2_result = results[1]
println("SiO2 molecular weight: ", sio2_result.molecular_weight)
```
"""
function calculate_xray_properties(
    formulaList::Vector{String},
    energy::Vector{Float64},
    massDensityList::Vector{Float64},
)
    return _calculate_xray_properties_impl(formulaList, energy, massDensityList)
end

"""
    calculate_single_material_properties(formulaStr, energy, massDensity) -> XRayResult

Calculate X-ray optical properties for a single chemical formula.

# Arguments
- `formulaStr::String`: Chemical formula (e.g., "SiO2", "Al2O3")
- `energy::Vector{Float64}`: X-ray energies in KeV
- `massDensity::Float64`: Mass density in g/cm³

# Returns
- `XRayResult`: Complete set of calculated optical properties

# Examples
```julia
result = calculate_single_material_properties("SiO2", [8.0, 10.0, 12.0], 2.2)
println("Molecular weight: ", result.molecular_weight)
println("Critical angles: ", result.critical_angle)
```
"""
function calculate_single_material_properties(
    formulaStr::String,
    energy::Vector{Float64},
    massDensity::Float64,
)
    return _calculate_single_material_impl(formulaStr, energy, massDensity)
end

# =====================================================================================
# UTILITY FUNCTIONS
# =====================================================================================

"""
    clear_caches!()

Clear all cached data to free memory.

This function clears both the atomic data cache and the scattering-factor interpolator cache.
Useful for memory management in long-running applications or when processing
many different materials.

# Examples
```julia
# After processing many materials, clear caches to free memory
XRayLabTool.clear_caches!()  # not exported
```
"""
function clear_caches!()
    lock(ATOMIC_DATA_LOCK) do
        empty!(ATOMIC_DATA_CACHE)
    end
    lock(INTERPOLATOR_LOCK) do
        empty!(INTERPOLATOR_CACHE)
    end
    @debug "Caches cleared" _group = :cache
end

# =====================================================================================
# DEPRECATED ALIASES
# =====================================================================================

# Deprecated aliases for backward compatibility
# These will issue deprecation warnings when used

# Removed in v0.7: get_atomic_data, load_f1f2_table, calculate_scattering_factors!
# Removed in v0.6: load_scattering_factor_table, create_interpolators, pchip_interpolators

# =====================================================================================
# DEPRECATED PUBLIC API
# =====================================================================================

"""
    Refrac(formulaList, energy, massDensityList) -> Vector{XRayResult}

!!! warning "Deprecated"
    Use [`calculate_xray_properties`](@ref) instead.
"""
function Refrac(
    formulaList::Vector{String},
    energies_keV::Vector{Float64},
    massDensityList::Vector{Float64},
)
    Base.depwarn(
        "`Refrac` is deprecated, use `calculate_xray_properties` instead.",
        :Refrac,
    )
    return _calculate_xray_properties_impl(formulaList, energies_keV, massDensityList)
end

"""
    SubRefrac(formulaStr, energy, massDensity) -> XRayResult

!!! warning "Deprecated"
    Use [`calculate_single_material_properties`](@ref) instead.
"""
function SubRefrac(formulaStr::String, energies_keV::Vector{Float64}, massDensity::Float64)
    Base.depwarn(
        "`SubRefrac` is deprecated, use `calculate_single_material_properties` instead.",
        :SubRefrac,
    )
    return _calculate_single_material_impl(formulaStr, energies_keV, massDensity)
end

# Cache compiled specializations of the public entry points in the package image, cutting
# first-call latency. Plain `precompile` (no workload) so no cache state is serialized.
precompile(calculate_single_material_properties, (String, Vector{Float64}, Float64))
precompile(calculate_xray_properties, (Vector{String}, Vector{Float64}, Vector{Float64}))
precompile(Refrac, (Vector{String}, Vector{Float64}, Vector{Float64}))
precompile(SubRefrac, (String, Vector{Float64}, Float64))

end  # module
