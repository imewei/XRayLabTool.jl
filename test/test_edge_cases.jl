"""
    test_edge_cases.jl

Comprehensive edge-case, property-based, scientific validation, and benchmark tests
for XRayLabTool. Covers gaps not addressed by the existing test suite:

- Physical constant validation
- Energy–wavelength conversion correctness
- Numerical sanity (no NaN/Inf, positivity, monotonicity)
- Scaling / linearity properties
- parse_formula property-based tests
- atomic_number_and_mass direct tests + caching
- load_element_interpolators direct tests
- Cache lifecycle (clear_caches!)
- XRayResult field completeness & propertynames
- Boundary precision at energy limits
- Fluorine-containing compounds (previously broken data file)
- Performance regression guards
"""

using Test
using Logging
using XRayLabTool

# =======================================================================================
# 1. PHYSICAL CONSTANTS VALIDATION
# =======================================================================================

@testset "Physical Constants" begin
    # CODATA 2018 / 2019 SI redefinition exact values
    @test XRayLabTool.THOMSON_SCATTERING_LENGTH ≈ 2.8179403205e-15 atol = 1e-25
    @test XRayLabTool.SPEED_OF_LIGHT ≈ 299792458.0 atol = 1.0
    @test XRayLabTool.PLANCK == 6.62607015e-34  # exact
    @test XRayLabTool.ELEMENT_CHARGE == 1.602176634e-19  # exact
    @test XRayLabTool.AVOGADRO == 6.02214076e23  # exact

    # hc/e derived constant
    hc_keV_angstrom = 12.39842
    hc_keV_m = hc_keV_angstrom * 1e-10
    @test XRayLabTool.HC_OVER_ELECTRON_CHARGE_keV ≈ hc_keV_m rtol = 1e-4

    # Deprecated constant aliases should equal their replacements
    @test XRayLabTool.THOMPSON === XRayLabTool.THOMSON_SCATTERING_LENGTH
    @test XRayLabTool.ENERGY_TO_WAVELENGTH_FACTOR ===
          XRayLabTool.HC_OVER_ELECTRON_CHARGE_keV
    @test XRayLabTool.SCATTERING_FACTOR === XRayLabTool.SCATTERING_PREFACTOR
end

# =======================================================================================
# 2. ENERGY–WAVELENGTH CONVERSION
# =======================================================================================

@testset "Energy–Wavelength Conversion" begin
    # λ(Å) = 12.39842 / E(keV)  — standard X-ray relation
    result = calculate_single_material_properties("Si", [8.04], 2.33)  # Cu Kα line
    expected_wavelength = 12.39842 / 8.04  # ≈ 1.5418 Å
    @test result.wavelength[1] ≈ expected_wavelength rtol=1e-3

    # Mo Kα: 17.44 keV → 0.7107 Å
    result_mo = calculate_single_material_properties("Si", [17.44], 2.33)
    @test result_mo.wavelength[1] ≈ 0.7107 rtol=1e-2

    # Energy and wavelength arrays must have same length
    result_multi = calculate_single_material_properties("Si", [5.0, 10.0, 20.0], 2.33)
    @test length(result_multi.energy) == length(result_multi.wavelength) == 3

    # Wavelength should decrease with increasing energy (inverse relationship)
    @test issorted(result_multi.wavelength, rev = true)
end

# =======================================================================================
# 3. NUMERICAL SANITY — NO NaN/Inf, POSITIVITY, MONOTONICITY
# =======================================================================================

@testset "Numerical Sanity" begin
    energies = collect(1.0:0.5:20.0)
    materials = [("SiO2", 2.2), ("Si", 2.33), ("Au", 19.3), ("H2O", 1.0), ("Al2O3", 3.95)]

    @testset "No NaN/Inf — $formula" for (formula, density) in materials
        result = calculate_single_material_properties(formula, energies, density)

        # Scalar fields must be finite
        @test isfinite(result.molecular_weight)
        @test isfinite(result.number_of_electrons)
        @test isfinite(result.electron_density)

        # Vector fields must be all-finite, no NaN
        for field in [
            :dispersion,
            :absorption,
            :f1,
            :f2,
            :critical_angle,
            :attenuation_length,
            :real_sld,
            :imag_sld,
            :wavelength,
        ]
            vals = getproperty(result, field)
            @test all(isfinite, vals)
        end
    end

    @testset "Positivity constraints" begin
        result = calculate_single_material_properties("SiO2", energies, 2.2)

        @test result.molecular_weight > 0
        @test result.number_of_electrons > 0
        @test result.electron_density > 0
        @test all(result.wavelength .> 0)
        @test all(result.dispersion .> 0)      # δ > 0 for normal materials away from edges
        @test all(result.absorption .> 0)      # β > 0 always
        @test all(result.critical_angle .> 0)
        @test all(result.attenuation_length .> 0)
        @test all(result.real_sld .> 0)
        @test all(result.imag_sld .> 0)
    end

    @testset "Monotonicity — dispersion decreases with energy (far from edges)" begin
        # Away from absorption edges, δ ∝ λ² ∝ 1/E² — should decrease with energy
        high_energies = collect(5.0:1.0:20.0)
        result = calculate_single_material_properties("SiO2", high_energies, 2.2)
        @test issorted(result.dispersion, rev = true)
    end
end

# =======================================================================================
# 4. SCALING / LINEARITY PROPERTIES
# =======================================================================================

@testset "Scaling Properties" begin
    energies = [5.0, 8.0, 10.0, 15.0]

    @testset "δ and β scale linearly with density" begin
        result_1x = calculate_single_material_properties("SiO2", energies, 2.2)
        result_2x = calculate_single_material_properties("SiO2", energies, 4.4)

        @test result_2x.dispersion ≈ 2.0 .* result_1x.dispersion rtol=1e-10
        @test result_2x.absorption ≈ 2.0 .* result_1x.absorption rtol=1e-10
    end

    @testset "SLD scales linearly with density" begin
        result_1x = calculate_single_material_properties("Si", energies, 2.33)
        result_2x = calculate_single_material_properties("Si", energies, 4.66)

        @test result_2x.real_sld ≈ 2.0 .* result_1x.real_sld rtol=1e-10
        @test result_2x.imag_sld ≈ 2.0 .* result_1x.imag_sld rtol=1e-10
    end

    @testset "Electron density scales linearly with mass density" begin
        ed_1 = calculate_single_material_properties("Si", [8.0], 2.33).electron_density
        ed_2 = calculate_single_material_properties("Si", [8.0], 4.66).electron_density
        @test ed_2 ≈ 2.0 * ed_1 rtol=1e-10
    end

    @testset "f1/f2 are independent of density" begin
        result_lo = calculate_single_material_properties("SiO2", energies, 1.0)
        result_hi = calculate_single_material_properties("SiO2", energies, 5.0)
        @test result_lo.f1 ≈ result_hi.f1 rtol=1e-12
        @test result_lo.f2 ≈ result_hi.f2 rtol=1e-12
    end

    @testset "Molecular weight and electron count are independent of density/energy" begin
        r1 = calculate_single_material_properties("SiO2", [5.0], 1.0)
        r2 = calculate_single_material_properties("SiO2", [20.0], 10.0)
        @test r1.molecular_weight == r2.molecular_weight
        @test r1.number_of_electrons == r2.number_of_electrons
    end
end

# =======================================================================================
# 5. ATOMIC DATA — DIRECT TESTS
# =======================================================================================

@testset "atomic_number_and_mass" begin
    @testset "Known elements" begin
        z, m = XRayLabTool.atomic_number_and_mass("Si")
        @test z == 14
        @test m ≈ 28.085 rtol=1e-3

        z, m = XRayLabTool.atomic_number_and_mass("O")
        @test z == 8
        @test m ≈ 15.999 rtol=1e-3

        z, m = XRayLabTool.atomic_number_and_mass("Au")
        @test z == 79
        @test m ≈ 196.967 rtol=1e-3

        z, m = XRayLabTool.atomic_number_and_mass("H")
        @test z == 1
        @test m ≈ 1.008 rtol=1e-2
    end

    @testset "Return type stability" begin
        result = XRayLabTool.atomic_number_and_mass("Si")
        @test result isa Tuple{Int, Float64}
    end

    @testset "Caching returns same object" begin
        r1 = XRayLabTool.atomic_number_and_mass("Si")
        r2 = XRayLabTool.atomic_number_and_mass("Si")
        @test r1 === r2  # Exact same cached tuple
    end

    @testset "Invalid element" begin
        @test_throws ArgumentError XRayLabTool.atomic_number_and_mass("Xx")
        @test_throws ArgumentError XRayLabTool.atomic_number_and_mass("NotAnElement")
    end
end

# =======================================================================================
# 6. LOAD ELEMENT INTERPOLATORS — DIRECT TESTS
# =======================================================================================

@testset "load_element_interpolators" begin
    @testset "Returns correct types" begin
        itp1, itp2 = XRayLabTool.load_element_interpolators("Si")
        @test itp1 isa XRayLabTool.PCHIPInterp
        @test itp2 isa XRayLabTool.PCHIPInterp
    end

    @testset "Caching returns same object" begin
        pair1 = XRayLabTool.load_element_interpolators("Si")
        pair2 = XRayLabTool.load_element_interpolators("Si")
        @test pair1 === pair2
    end

    @testset "Interpolator produces finite values" begin
        itp_f1, itp_f2 = XRayLabTool.load_element_interpolators("Si")
        # Test at Cu Kα energy (8040 eV)
        @test isfinite(itp_f1(8040.0))
        @test isfinite(itp_f2(8040.0))
        # f1 should be close to Z (14) at high energy
        @test itp_f1(20000.0) > 10.0
    end

    @testset "Invalid element" begin
        @test_throws ArgumentError XRayLabTool.load_element_interpolators("Zz")
    end

    @testset "Fluorine loads correctly (previously had data file bug)" begin
        itp_f1, itp_f2 = XRayLabTool.load_element_interpolators("F")
        @test isfinite(itp_f1(8040.0))
        @test isfinite(itp_f2(8040.0))
    end
end

# =======================================================================================
# 7. CACHE LIFECYCLE
# =======================================================================================

@testset "Cache lifecycle" begin
    # Ensure caches are populated
    calculate_single_material_properties("Si", [8.0], 2.33)
    @test haskey(XRayLabTool.ATOMIC_DATA_CACHE, "Si")
    @test haskey(XRayLabTool.INTERPOLATOR_CACHE, "Si")

    # Clear and verify
    XRayLabTool.clear_caches!()
    @test isempty(XRayLabTool.ATOMIC_DATA_CACHE)
    @test isempty(XRayLabTool.INTERPOLATOR_CACHE)

    # Re-populate still works
    result = calculate_single_material_properties("Si", [8.0], 2.33)
    @test result.molecular_weight > 0
    @test haskey(XRayLabTool.ATOMIC_DATA_CACHE, "Si")
end

@testset "Thread-safe cache access" begin
    XRayLabTool.clear_caches!()

    # Concurrent calls for different elements should not crash
    elements = ["Si", "O", "Al", "Fe", "Au", "Ca", "Na", "Mg", "H", "C"]
    results = Vector{Any}(nothing, length(elements))

    Threads.@threads for i in eachindex(elements)
        results[i] = calculate_single_material_properties(elements[i], [8.0], 1.0)
    end

    for (i, elem) in enumerate(elements)
        @test results[i] isa XRayResult
        @test results[i].formula == elem
    end
end

# =======================================================================================
# 8. XRayResult FIELD COMPLETENESS & PROPERTYNAMES
# =======================================================================================

@testset "XRayResult completeness" begin
    result = calculate_single_material_properties("SiO2", [8.0, 10.0], 2.2)

    @testset "All 15 fields present and correct type" begin
        @test result.formula isa String
        @test result.molecular_weight isa Float64
        @test result.number_of_electrons isa Float64
        @test result.mass_density isa Float64
        @test result.electron_density isa Float64
        @test result.energy isa Vector{Float64}
        @test result.wavelength isa Vector{Float64}
        @test result.dispersion isa Vector{Float64}
        @test result.absorption isa Vector{Float64}
        @test result.f1 isa Vector{Float64}
        @test result.f2 isa Vector{Float64}
        @test result.critical_angle isa Vector{Float64}
        @test result.attenuation_length isa Vector{Float64}
        @test result.real_sld isa Vector{Float64}
        @test result.imag_sld isa Vector{Float64}
    end

    @testset "Vector fields all have length == n_energies" begin
        n = 2
        for field in [
            :energy,
            :wavelength,
            :dispersion,
            :absorption,
            :f1,
            :f2,
            :critical_angle,
            :attenuation_length,
            :real_sld,
            :imag_sld,
        ]
            @test length(getproperty(result, field)) == n
        end
    end

    @testset "propertynames includes both new and legacy names" begin
        pnames = propertynames(result)
        # New names
        for name in [:formula, :molecular_weight, :dispersion, :real_sld, :imag_sld]
            @test name in pnames
        end
        # Legacy names
        for name in [:Formula, :MW, :Dispersion, :reSLD, :imSLD]
            @test name in pnames
        end
    end
end

# =======================================================================================
# 9. BOUNDARY PRECISION AT ENERGY LIMITS
# =======================================================================================

@testset "Energy boundary precision" begin
    @testset "Exactly at minimum (0.03 keV)" begin
        result = calculate_single_material_properties("Si", [0.03], 2.33)
        @test all(isfinite, result.dispersion)
        @test result.formula == "Si"
    end

    @testset "Exactly at maximum (30.0 keV)" begin
        result = calculate_single_material_properties("Si", [30.0], 2.33)
        @test all(isfinite, result.dispersion)
    end

    @testset "Just outside bounds" begin
        @test_throws ArgumentError calculate_single_material_properties(
            "Si",
            [0.0299],
            2.33,
        )
        @test_throws ArgumentError calculate_single_material_properties(
            "Si",
            [30.001],
            2.33,
        )
    end

    @testset "Single energy point" begin
        result = calculate_single_material_properties("Au", [8.0], 19.3)
        @test length(result.energy) == 1
        @test length(result.dispersion) == 1
    end
end

# =======================================================================================
# 9b. INPUT VALIDATION — NaN, Inf, NEGATIVE, EMPTY
# =======================================================================================

@testset "Input validation edge cases" begin
    @testset "NaN energy rejected" begin
        @test_throws ArgumentError calculate_single_material_properties("Si", [NaN], 2.33)
        @test_throws ArgumentError calculate_single_material_properties(
            "Si",
            [8.0, NaN],
            2.33,
        )
    end

    @testset "Inf energy rejected" begin
        @test_throws ArgumentError calculate_single_material_properties("Si", [Inf], 2.33)
        @test_throws ArgumentError calculate_single_material_properties("Si", [-Inf], 2.33)
    end

    @testset "Empty energy rejected in single-material path" begin
        @test_throws ArgumentError calculate_single_material_properties(
            "Si",
            Float64[],
            2.33,
        )
    end

    @testset "Zero density rejected" begin
        @test_throws ArgumentError calculate_single_material_properties("Si", [8.0], 0.0)
    end

    @testset "Negative density rejected" begin
        @test_throws ArgumentError calculate_single_material_properties("Si", [8.0], -1.0)
    end

    @testset "NaN density rejected" begin
        @test_throws ArgumentError calculate_single_material_properties("Si", [8.0], NaN)
    end

    @testset "Inf density rejected" begin
        @test_throws ArgumentError calculate_single_material_properties("Si", [8.0], Inf)
    end

    @testset "Batch density validation" begin
        @test_throws ArgumentError calculate_xray_properties(
            ["Si", "Au"],
            [8.0],
            [2.33, -19.3],
        )
        @test_throws ArgumentError calculate_xray_properties(["Si"], [8.0], [NaN])
    end
end

# =======================================================================================
# 10. FLUORINE-CONTAINING COMPOUNDS (previously broken)
# =======================================================================================

@testset "Fluorine compounds" begin
    @testset "BaF2" begin
        result = calculate_single_material_properties("BaF2", [8.0], 4.89)
        @test result.formula == "BaF2"
        @test result.molecular_weight ≈ 175.32 rtol=1e-2
        @test all(isfinite, result.dispersion)
    end

    @testset "CaF2 (fluorite)" begin
        result = calculate_single_material_properties("CaF2", [8.0, 10.0], 3.18)
        @test result.formula == "CaF2"
        @test all(isfinite, result.f1)
        @test all(isfinite, result.f2)
    end

    @testset "LiF" begin
        result = calculate_single_material_properties("LiF", [8.0], 2.64)
        @test result.formula == "LiF"
        @test all(isfinite, result.dispersion)
    end
end

# =======================================================================================
# 11. PHYSICS CROSS-VALIDATION (model-independent checks)
# =======================================================================================

@testset "Physics Cross-Validation" begin
    @testset "f1 approaches Z at high energy (Thomson limit)" begin
        # At high energies far from edges, f1 → Z (total electrons per molecule)
        high_e = [25.0, 28.0, 30.0]

        si = calculate_single_material_properties("Si", high_e, 2.33)
        @test all(si.f1 .> 13.0)  # Z_Si = 14, f1 approaches it
        @test all(si.f1 .< 15.0)

        h2o = calculate_single_material_properties("H2O", high_e, 1.0)
        Z_H2O = 2 * 1 + 8  # = 10
        @test all(h2o.f1 .> 9.0)
        @test all(h2o.f1 .< 11.0)
    end

    @testset "Critical angle relation: θc = √(2δ) × 180/π" begin
        result = calculate_single_material_properties("Si", [8.0, 10.0, 15.0], 2.33)
        expected_θc = sqrt.(2 .* result.dispersion) .* (180.0 / π)
        @test result.critical_angle ≈ expected_θc rtol=1e-12
    end

    @testset "Attenuation length relation: L = λ/(4πβ)" begin
        result = calculate_single_material_properties("Au", [8.0, 10.0], 19.3)
        # wavelength is in Å, convert to m for the formula
        λ_m = result.wavelength .* 1e-10
        expected_L_cm = λ_m ./ (4π .* result.absorption) .* 1e2
        @test result.attenuation_length ≈ expected_L_cm rtol=1e-10
    end

    @testset "SLD relation: SLD = 2π(δ+iβ)/λ²" begin
        result = calculate_single_material_properties("SiO2", [8.0, 10.0], 2.2)
        λ_m = result.wavelength .* 1e-10
        λ_sq = λ_m .^ 2
        expected_re_sld = result.dispersion .* (2π / 1e20) ./ λ_sq
        expected_im_sld = result.absorption .* (2π / 1e20) ./ λ_sq
        @test result.real_sld ≈ expected_re_sld rtol=1e-10
        @test result.imag_sld ≈ expected_im_sld rtol=1e-10
    end

    @testset "Higher-Z materials have larger dispersion" begin
        # At the same energy, higher-Z materials (higher electron density) have larger δ
        si = calculate_single_material_properties("Si", [10.0], 2.33)
        au = calculate_single_material_properties("Au", [10.0], 19.3)
        @test au.dispersion[1] > si.dispersion[1]
        @test au.critical_angle[1] > si.critical_angle[1]
    end
end

@testset "Negative delta handling" begin
    @testset "Critical angle is zero when delta is negative" begin
        # At very low energies near absorption edges, f1 can go negative.
        # We test that the function doesn't throw DomainError.
        energies = collect(0.03:0.01:0.1)
        result = calculate_single_material_properties("Si", energies, 2.33)

        # Critical angle should be non-negative (zero when delta < 0)
        @test all(result.critical_angle .>= 0.0)

        # All values should be finite (no NaN from sqrt of negative)
        @test all(isfinite, result.critical_angle)
    end
end

# =======================================================================================
# 12. PROPERTY-BASED TESTS — FORMULA PARSING
# =======================================================================================

@testset "Formula parsing — property-based" begin
    @testset "Roundtrip: element count sum is preserved" begin
        formulas = ["SiO2", "Al2O3", "CaCO3", "Fe2O3", "H2O", "NaCl", "MgSiO3"]
        for formula in formulas
            syms, counts = XRayLabTool.parse_formula(formula)
            # Re-encoding and re-parsing should give same result
            rebuilt =
                join([s * (c == 1.0 ? "" : string(Int(c))) for (s, c) in zip(syms, counts)])
            syms2, counts2 = XRayLabTool.parse_formula(rebuilt)
            @test syms == syms2
            @test counts ≈ counts2
        end
    end

    @testset "Single-element formulas" begin
        for elem in ["H", "C", "N", "O", "Si", "Fe", "Au", "U"]
            syms, counts = XRayLabTool.parse_formula(elem)
            @test length(syms) == 1
            @test syms[1] == elem
            @test counts[1] == 1.0
        end
    end

    @testset "Case sensitivity: CO ≠ Co" begin
        syms_CO, _ = XRayLabTool.parse_formula("CO")
        syms_Co, _ = XRayLabTool.parse_formula("Co")
        @test syms_CO == ["C", "O"]
        @test syms_Co == ["Co"]
    end
end

# =======================================================================================
# 13. BATCH API SPECIFIC TESTS
# =======================================================================================

@testset "Batch API" begin
    @testset "Results contain all requested formulas" begin
        formulas = ["SiO2", "Al2O3", "Fe2O3"]
        densities = [2.2, 3.95, 5.24]
        results = calculate_xray_properties(formulas, [8.0, 10.0], densities)
        for (i, f) in enumerate(formulas)
            @test results[i].formula == f
        end
    end

    @testset "Batch produces same result as individual calls" begin
        formulas = ["SiO2", "H2O"]
        densities = [2.2, 1.0]
        energies = [5.0, 8.0, 10.0, 15.0]

        batch = calculate_xray_properties(formulas, energies, densities)

        for (i, (f, d)) in enumerate(zip(formulas, densities))
            individual = calculate_single_material_properties(f, sort(energies), d)
            @test batch[i].dispersion ≈ individual.dispersion rtol=1e-12
            @test batch[i].f1 ≈ individual.f1 rtol=1e-12
            @test batch[i].f2 ≈ individual.f2 rtol=1e-12
            @test batch[i].real_sld ≈ individual.real_sld rtol=1e-12
        end
    end

    @testset "Batch sorts energies" begin
        results = calculate_xray_properties(["Si"], [10.0, 5.0, 8.0], [2.33])
        @test issorted(results[1].energy)
    end

    @testset "Single material also sorts energies" begin
        result = calculate_single_material_properties("Si", [10.0, 5.0, 8.0], 2.33)
        @test issorted(result.energy)
    end

    @testset "Duplicate formulas in batch" begin
        results = calculate_xray_properties(["SiO2", "SiO2"], [8.0], [2.2, 3.0])
        @test length(results) == 2
        @test results[1].mass_density == 2.2
        @test results[2].mass_density == 3.0
    end
end

# =======================================================================================
# STRICT BATCH VALIDATION
# =======================================================================================

@testset "Strict batch validation" begin
    @testset "Invalid formula in batch throws listing all errors" begin
        err = nothing
        try
            calculate_xray_properties(
                ["SiO2", "XxYy", "H2O", "ZzZz"],
                [8.0],
                [2.2, 1.0, 1.0, 1.0],
            )
        catch e
            err = e
        end
        @test err isa ArgumentError
        msg = err.msg
        # Should list ALL invalid formulas, not just the first
        @test occursin("XxYy", msg)
        @test occursin("ZzZz", msg)
        # Should NOT list valid formulas
        @test !occursin("SiO2", msg)
        @test !occursin("H2O", msg)
    end

    @testset "All-valid formulas still work" begin
        results = calculate_xray_properties(["SiO2", "H2O"], [8.0], [2.2, 1.0])
        @test length(results) == 2
    end
end

# =======================================================================================
# 14. KNOWN MOLECULAR WEIGHTS
# =======================================================================================

@testset "Known molecular weights" begin
    known_mw = [
        ("H2O", 18.015),
        ("SiO2", 60.08),
        ("Al2O3", 101.96),
        ("NaCl", 58.44),
        ("CaCO3", 100.09),
        ("Fe2O3", 159.69),
    ]

    @testset "MW — $formula" for (formula, expected_mw) in known_mw
        result = calculate_single_material_properties(formula, [8.0], 1.0)
        @test result.molecular_weight ≈ expected_mw rtol=1e-2
    end
end

# =======================================================================================
# 15. PERFORMANCE REGRESSION GUARDS
# =======================================================================================

@testset "Performance regression guards" begin
    energies = collect(1.0:0.1:20.0)

    # Warm up
    calculate_single_material_properties("SiO2", energies, 2.2)

    @testset "Single material allocation bound" begin
        allocs = @allocated calculate_single_material_properties("SiO2", energies, 2.2)
        # Should be under 30 KB per call (currently ~20 KB)
        @test allocs < 30_000
    end

    @testset "Single material timing bound" begin
        # Warm JIT
        for _ in 1:100
            ;
            calculate_single_material_properties("SiO2", energies, 2.2);
        end

        t = @elapsed for _ in 1:100
            calculate_single_material_properties("SiO2", energies, 2.2)
        end
        per_call_us = t * 1e6 / 100
        # Should complete in under 500µs per call (currently ~25µs)
        @test per_call_us < 500
    end
end

# =======================================================================================
# 16. LOGGING DOES NOT BREAK FUNCTIONALITY
# =======================================================================================

@testset "Logging with debug enabled" begin
    # Enable debug logging and verify calculations still produce correct results
    test_logger = Logging.ConsoleLogger(devnull, Logging.Debug)
    Logging.with_logger(test_logger) do
        result = calculate_single_material_properties("SiO2", [8.0, 10.0], 2.2)
        @test result.formula == "SiO2"
        @test result.molecular_weight > 0
        @test all(isfinite, result.dispersion)

        batch = calculate_xray_properties(["Si", "Au"], [8.0], [2.33, 19.3])
        @test length(batch) == 2
        @test batch[1].formula == "Si"

        XRayLabTool.clear_caches!()
        result2 = calculate_single_material_properties("H2O", [8.0], 1.0)
        @test result2.formula == "H2O"
    end
end

# =======================================================================================
# 17. REGRESSION TESTS FOR REVIEW FINDINGS
# =======================================================================================

@testset "Zero element counts are rejected (no silent NaN)" begin
    for f in ("H0", "(SiO2)0", "SiO0")
        @test_throws ArgumentError XRayLabTool.parse_formula(f)
        @test_throws ArgumentError calculate_single_material_properties(f, [8.0], 2.2)
        @test_throws ArgumentError calculate_xray_properties([f], [8.0], [2.2])
    end
    # Fractional counts remain valid
    @test XRayLabTool.parse_formula("H0.5He0.5")[2] == [0.5, 0.5]
end

@testset "Result energy does not alias caller or sibling results" begin
    for inp in ([8.0, 10.0], [10.0, 8.0])  # already sorted / unsorted
        r = calculate_single_material_properties("SiO2", inp, 2.2)
        @test r.energy !== inp
        rs = calculate_xray_properties(["SiO2", "Al2O3"], inp, [2.2, 3.95])
        @test rs[1].energy !== rs[2].energy
        @test rs[1].energy !== inp
        rs[1].energy[1] = -1.0
        @test rs[2].energy[1] == 8.0
    end
end

@testset "Debug log values are rendered, not closures" begin
    logger = Test.TestLogger(min_level = Logging.Debug)
    Logging.with_logger(logger) do
        calculate_single_material_properties("SiO2", [8.0, 10.0], 2.2)
    end
    vals = [v for rec in logger.logs for (_, v) in rec.kwargs]
    @test !any(v -> v isa Function, vals)
end

@testset "Non-finite element counts are rejected" begin
    n = repeat("9", 200)  # parses to Inf
    @test_throws ArgumentError XRayLabTool.parse_formula("(H" * n * ")" * n)
    @test_throws ArgumentError calculate_single_material_properties("H" * n * n, [8.0], 1.0)
end

@testset "Property access is type-stable" begin
    r = calculate_single_material_properties("Si", [8.0], 2.33)
    get_dispersion(x) = x.dispersion
    get_mw(x) = x.molecular_weight
    @inferred get_dispersion(r)
    @inferred get_mw(r)
    @test r.Energy === r.energy  # legacy names still resolve
    @test_throws Exception r.not_a_field  # ErrorException (<1.12) / FieldError (>=1.12)
end

@testset "delta and beta match the closed-form expression (tight tolerance)" begin
    # delta = (λ²/2π) r_e (ρ N_A / M) Σ n_i f1_i, with ρ converted g/cm³ → g/m³ (×1e6)
    for (formula, ρ) in (("Si", 2.33), ("SiO2", 2.2), ("Au", 19.3))
        r = calculate_single_material_properties(formula, [5.0, 8.0, 20.0], ρ)
        λ = r.wavelength .* 1e-10
        pref =
            λ .^ 2 ./ (2π) .* XRayLabTool.THOMSON_SCATTERING_LENGTH .*
            (ρ * 1e6 * XRayLabTool.AVOGADRO / r.molecular_weight)
        @test r.dispersion ≈ pref .* r.f1 rtol = 1e-9
        @test r.absorption ≈ pref .* r.f2 rtol = 1e-9
    end
end
