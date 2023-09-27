using XRayLabTool
using Test

# Define a custom approximation testing function
function assert_approx_equal(actual, expected; tol = 1e-2)
    @test length(actual) == length(expected)
    @test all(abs.(actual .- expected) .< tol)
end

@testset "XRayLabTool.jl" begin
    data = Refrac(["SiO2", "H2O"], [1.0, 2.0, 3.5, 5.0, 6.0], [2.2, 1.0])

    @test haskey(data, "SiO2")
    @test haskey(data, "H2O")

    si = SubRefrac("Si", [20.0], 2.33)

    assert_approx_equal(data["SiO2"].Dispersion[3], 3.814915113516378e-5)
    assert_approx_equal(data["SiO2"].f1[1], 29.448492570232553)
    assert_approx_equal(data["SiO2"].reSLD[5], 1.8979332462750346e-5)

    assert_approx_equal(si.Dispersion[1], 1.2096870811637957e-6)
    assert_approx_equal(si.f1[1], 14.048053002299712)
    assert_approx_equal(si.reSLD[1], 1.977791074150516e-5)
end