using XRayLabTool
using Test
include("utils.jl")
using .Main: @expect_error
include("formula_parsing.jl")
include("test_edge_cases.jl")

# Test constants
const DEFAULT_TOL = 1e-6
const ENERGIES = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
const DENSITIES = [2.2, 1.0]
const MATERIALS = ["SiO2", "H2O"]

# Helper function for approximate equality testing
function test_approx_equal(actual, expected; tol = DEFAULT_TOL, description = "")
    if !isempty(description)
        @testset "$description" begin
            @test actual ≈ expected atol=tol
        end
    else
        @test actual ≈ expected atol=tol
    end
end

# Helper function to test array approximate equality
function test_array_approx_equal(actual, expected; tol = DEFAULT_TOL, description = "")
    if !isempty(description)
        @testset "$description" begin
            @test length(actual) == length(expected)
            @test all(abs.(actual .- expected) .< tol)
        end
    else
        @test length(actual) == length(expected)
        @test all(abs.(actual .- expected) .< tol)
    end
end

@testset "XRayLabTool.jl Tests" begin
    @testset "Basic Setup and Initialization - New API" begin
        # Create test data using new API
        data = calculate_xray_properties(MATERIALS, ENERGIES, DENSITIES)

        # Test that materials are properly initialized
        @test data[1].formula == "SiO2"
        @test data[2].formula == "H2O"

        # Test data structure integrity
        @test length(data) == 2
    end

    @testset "SiO2 Properties - New API" begin
        data = calculate_xray_properties(MATERIALS, ENERGIES, DENSITIES)
        sio2 = data[1]

        # Expected values for SiO2
        expected_dispersion = [(3, 9.451484792575434e-6), (5, 5.69919201789506e-06)]

        expected_f1 =
            [(1, 30.641090313037314), (3, 30.46419063207884), (5, 30.366953850108544)]

        expected_real_sld = [(3, 1.8929689855615698e-5), (5, 1.886926933936152e-5)]

        # Test dispersion values using new field name
        @testset "SiO2 Dispersion" begin
            for (idx, expected_val) in expected_dispersion
                test_approx_equal(
                    sio2.dispersion[idx],
                    expected_val,
                    description = "dispersion[$idx]",
                )
            end
        end

        # Test f1 values
        @testset "SiO2 f1 values" begin
            for (idx, expected_val) in expected_f1
                test_approx_equal(sio2.f1[idx], expected_val, description = "f1[$idx]")
            end
        end

        # Test real_sld values using new field name
        @testset "SiO2 real_sld values" begin
            for (idx, expected_val) in expected_real_sld
                test_approx_equal(
                    sio2.real_sld[idx],
                    expected_val,
                    description = "real_sld[$idx]",
                )
            end
        end
    end

    @testset "H2O Properties - New API" begin
        data = calculate_xray_properties(MATERIALS, ENERGIES, DENSITIES)
        h2o = data[2]

        # Expected values for H2O
        expected_dispersion = [(3, 4.734311949237782e-6), (5, 2.8574954896752405e-6)]

        expected_f1 =
            [(1, 10.110775776847062), (3, 10.065881494924541), (5, 10.04313810715771)]

        expected_real_sld = [(3, 9.482008260671003e-6), (5, 9.460584107129207e-6)]

        # Test dispersion values using new field name
        @testset "H2O Dispersion" begin
            for (idx, expected_val) in expected_dispersion
                test_approx_equal(
                    h2o.dispersion[idx],
                    expected_val,
                    description = "dispersion[$idx]",
                )
            end
        end

        # Test f1 values
        @testset "H2O f1 values" begin
            for (idx, expected_val) in expected_f1
                test_approx_equal(h2o.f1[idx], expected_val, description = "f1[$idx]")
            end
        end

        # Test real_sld values using new field name
        @testset "H2O real_sld values" begin
            for (idx, expected_val) in expected_real_sld
                test_approx_equal(
                    h2o.real_sld[idx],
                    expected_val,
                    description = "real_sld[$idx]",
                )
            end
        end
    end

    @testset "Single Material Silicon Properties - New API" begin
        si = calculate_single_material_properties("Si", [20.0], 2.33)

        # Expected values for Silicon
        expected_values = [
            (:dispersion, 1, 1.20966554922812e-06),
            (:f1, 1, 14.048053047106292),
            (:f2, 1, 0.053331074920700626),
            (:real_sld, 1, 1.9777910804587255e-5),
            (:imag_sld, 1, 7.508351793358633e-8),
        ]

        @testset "Silicon property: $(property)" for (property, idx, expected_val) in
                                                     expected_values

            actual_val = getproperty(si, property)[idx]
            test_approx_equal(actual_val, expected_val, description = "$(property)[$idx]")
        end
    end

    @testset "Edge Cases and Error Handling - New API" begin
        # Test with empty materials (if applicable)
        @test_throws Exception calculate_xray_properties(String[], ENERGIES, DENSITIES)

        # Test with mismatched array lengths
        @test_throws Exception calculate_xray_properties(MATERIALS, ENERGIES, [1.0])  # Wrong density count

        # Test that results have correct length
        data = calculate_xray_properties(MATERIALS, ENERGIES, DENSITIES)
        @test length(data) == 2
    end

    @testset "Property Consistency - New API" begin
        data = calculate_xray_properties(MATERIALS, ENERGIES, DENSITIES)

        # Test that all materials have the same energy array length
        for (i, material) in enumerate(MATERIALS)
            material_data = data[i]
            @test length(material_data.f1) == length(ENERGIES)
            @test length(material_data.dispersion) == length(ENERGIES)
            @test length(material_data.real_sld) == length(ENERGIES)
        end
    end

    @testset "Input Validation and Error Handling Tests - New API" begin
        @testset "Multi-Material Function Error Handling" begin
            # Test energy outside 0.03–30 keV range (below minimum)
            @testset "Energy below 0.03 keV" begin
                @expect_error calculate_xray_properties(["SiO2"], [0.02, 5.0], [2.2]) r"Energy is out of range.*0\.03KeV.*30\.0KeV"
            end

            # Test energy outside 0.03–30 keV range (above maximum)
            @testset "Energy above 30 keV" begin
                @expect_error calculate_xray_properties(["SiO2"], [5.0, 35.0], [2.2]) r"Energy is out of range.*0\.03KeV.*30\.0KeV"
            end

            # Test mismatched formulaList & massDensityList lengths
            @testset "Mismatched list lengths" begin
                @expect_error calculate_xray_properties(
                    ["SiO2", "Al2O3"],
                    [8.0, 10.0],
                    [2.2],
                ) r"Formula list and mass density list must have the same length"
            end

            # Test non-vector inputs (pass scalar formula)
            @testset "Non-vector formula input" begin
                @expect_error calculate_xray_properties("SiO2", [8.0, 10.0], [2.2]) r"MethodError.*calculate_xray_properties"
            end

            # Test non-vector inputs (pass scalar energy)
            @testset "Non-vector energy input" begin
                @expect_error calculate_xray_properties(["SiO2"], 8.0, [2.2]) r"MethodError.*calculate_xray_properties"
            end

            # Test non-vector inputs (pass scalar density)
            @testset "Non-vector density input" begin
                @expect_error calculate_xray_properties(["SiO2"], [8.0, 10.0], 2.2) r"MethodError.*calculate_xray_properties"
            end

            # Test empty formula list
            @testset "Empty formula list" begin
                @expect_error calculate_xray_properties(String[], [8.0, 10.0], Float64[]) r"Formula list and energy vector must not be empty"
            end

            # Test empty energy vector
            @testset "Empty energy vector" begin
                @expect_error calculate_xray_properties(["SiO2"], Float64[], [2.2]) r"Formula list and energy vector must not be empty"
            end
        end

        @testset "Single Material Function Error Handling" begin
            # Test energy range validation on single-material path
            @testset "Single material energy below 0.03 keV" begin
                @expect_error calculate_single_material_properties("SiO2", [0.01], 2.2) r"Energy is out of range.*0\.03KeV.*30\.0KeV"
            end

            @testset "Single material energy above 30 keV" begin
                @expect_error calculate_single_material_properties("SiO2", [50.0], 2.2) r"Energy is out of range.*0\.03KeV.*30\.0KeV"
            end

            # Test invalid chemical formula
            @testset "Invalid chemical formula" begin
                @expect_error calculate_single_material_properties(
                    "InvalidElement123",
                    [8.0],
                    2.2,
                ) r"Element.*not found"
            end

            # Test empty formula string
            @testset "Empty formula string" begin
                @expect_error calculate_single_material_properties("", [8.0], 2.2) r"Invalid chemical formula"
            end
        end

        @testset "Duplicated Formulas - Independence Test" begin
            # Test that duplicated formulas produce independent results
            @testset "Duplicate formulas yield consistent results" begin
                formulas = ["SiO2", "SiO2", "Al2O3"]
                energies = [8.0, 10.0, 12.0]
                densities = [2.2, 2.2, 3.95]

                results = calculate_xray_properties(formulas, energies, densities)

                # All entries should be present in order
                @test results[1].formula == "SiO2"
                @test results[3].formula == "Al2O3"

                # Results should be consistent for the same material
                sio2_result = results[1]
                @test sio2_result.formula == "SiO2"
                @test sio2_result.mass_density == 2.2
                @test length(sio2_result.energy) == length(energies)
            end
        end

        @testset "Edge Case Input Values" begin
            # Test with boundary energy values
            @testset "Boundary energy values" begin
                # Test exactly at boundaries (should work)
                result_min = calculate_xray_properties(["SiO2"], [0.03], [2.2])
                @test result_min[1].formula == "SiO2"

                result_max = calculate_xray_properties(["SiO2"], [30.0], [2.2])
                @test result_max[1].formula == "SiO2"
            end

            # Test with very small density
            @testset "Very small density" begin
                result = calculate_single_material_properties("H2O", [8.0], 0.001)
                @test result.mass_density == 0.001
                @test result.formula == "H2O"
            end

            # Test with very large density
            @testset "Very large density" begin
                result = calculate_single_material_properties("Au", [8.0], 19.3)  # Gold density
                @test result.mass_density == 19.3
                @test result.formula == "Au"
            end
        end
    end

    # =====================================================================================
    # LEGACY API COMPATIBILITY TESTS
    # =====================================================================================

    @testset "Legacy API Compatibility Tests" begin
        @testset "Legacy Function Names Issue Deprecation Warnings" begin
            # Test that old function names issue deprecation warnings
            @test_deprecated data_legacy = Refrac(MATERIALS, ENERGIES, DENSITIES)
            @test_deprecated si_legacy = SubRefrac("Si", [20.0], 2.33)
        end

        @testset "Legacy Results Match New API" begin
            # Compare results between old and new APIs
            data_new = calculate_xray_properties(MATERIALS, ENERGIES, DENSITIES)
            data_legacy = @test_deprecated Refrac(MATERIALS, ENERGIES, DENSITIES)

            # Results should be identical
            @test data_new[1].dispersion ≈ data_legacy[1].dispersion
            @test data_new[1].f1 ≈ data_legacy[1].f1

            # Test single material functions
            si_new = calculate_single_material_properties("Si", [20.0], 2.33)
            si_legacy = @test_deprecated SubRefrac("Si", [20.0], 2.33)

            @test si_new.dispersion ≈ si_legacy.dispersion
            @test si_new.f1 ≈ si_legacy.f1
        end

        @testset "Legacy Field Names Still Work" begin
            # Test that old field names still work via property accessor
            data = calculate_xray_properties(["SiO2"], [8.0, 10.0], [2.2])
            sio2 = data[1]

            # Test deprecated field names work
            @test_nowarn sio2.Formula
            @test_nowarn sio2.MW
            @test_nowarn sio2.Number_Of_Electrons
            @test_nowarn sio2.Density
            @test_nowarn sio2.Electron_Density
            @test_nowarn sio2.Energy
            @test_nowarn sio2.Wavelength
            @test_nowarn sio2.Dispersion
            @test_nowarn sio2.Absorption
            @test_nowarn sio2.Critical_Angle
            @test_nowarn sio2.Attenuation_Length
            @test_nowarn sio2.reSLD
            @test_nowarn sio2.imSLD

            # Test that they return the same values as new field names
            @test sio2.Formula == sio2.formula
            @test sio2.MW == sio2.molecular_weight
            @test sio2.Number_Of_Electrons == sio2.number_of_electrons
            @test sio2.Density == sio2.mass_density
            @test sio2.Electron_Density == sio2.electron_density
            @test sio2.Energy == sio2.energy
            @test sio2.Wavelength == sio2.wavelength
            @test sio2.Dispersion == sio2.dispersion
            @test sio2.Absorption == sio2.absorption
            @test sio2.Critical_Angle == sio2.critical_angle
            @test sio2.Attenuation_Length == sio2.attenuation_length
            @test sio2.reSLD == sio2.real_sld
            @test sio2.imSLD == sio2.imag_sld
        end

        @testset "Legacy Tests Using Old Field Names" begin
            # Run some original tests using old field names to ensure backward compatibility
            data = @test_deprecated Refrac(MATERIALS, ENERGIES, DENSITIES)
            sio2 = data[1]

            # Expected values for SiO2 (same as before)
            expected_dispersion = [(3, 9.451484792575434e-6), (5, 5.69919201789506e-06)]
            expected_f1 =
                [(1, 30.641090313037314), (3, 30.46419063207884), (5, 30.366953850108544)]
            expected_reSLD = [(3, 1.8929689855615698e-5), (5, 1.886926933936152e-5)]

            # Test using legacy field names
            @testset "SiO2 Legacy Dispersion" begin
                for (idx, expected_val) in expected_dispersion
                    test_approx_equal(
                        sio2.Dispersion[idx],
                        expected_val,
                        description = "Legacy Dispersion[$idx]",
                    )
                end
            end

            @testset "SiO2 Legacy reSLD" begin
                for (idx, expected_val) in expected_reSLD
                    test_approx_equal(
                        sio2.reSLD[idx],
                        expected_val,
                        description = "Legacy reSLD[$idx]",
                    )
                end
            end
        end
    end
end

# =====================================================================================
# BENCHMARK TESTS
# =====================================================================================

@testset "Performance benchmarks" begin
    println("Running performance tests...")

    # Time single calculation using new API
    @time single_result =
        calculate_single_material_properties("SiO2", collect(1.0:0.1:20.0), 2.2)

    # Time multi calculation using new API
    formulas = ["SiO2", "Al2O3", "Fe2O3", "CaCO3", "MgO"]
    densities = [2.2, 3.95, 5.24, 2.71, 3.58]
    @time multi_result =
        calculate_xray_properties(formulas, collect(1.0:0.1:20.0), densities)

    println("Performance tests completed!")
end
