using XRayLabTool
using Test

# Define a custom approximation testing function
function assert_approx_equal(actual, expected; tol = eps())
    @test length(actual) == length(expected)
    @test all(abs.(actual .- expected) .< tol)
end

@testset "XRayLabTool.jl" begin
    data = Refrac(["SiO2", "H2O"], [5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [2.2, 1.0])

    @test haskey(data, "SiO2")
    @test haskey(data, "H2O")

    si = SubRefrac("Si", [20.0], 2.33)

    assert_approx_equal(data["SiO2"].Dispersion[3], 9.4512802977202e-06)
    assert_approx_equal(data["SiO2"].Dispersion[5], 5.69919201789506e-06)
    assert_approx_equal(data["SiO2"].f1[1], 30.6410903130373)
    assert_approx_equal(data["SiO2"].f1[3], 30.4641906320788)
    assert_approx_equal(data["SiO2"].f1[5], 30.3669538501085)
    assert_approx_equal(data["SiO2"].reSLD[3], 1.89292802877783e-05)
    assert_approx_equal(data["SiO2"].reSLD[5], 1.88560424934754e-05)

    assert_approx_equal(si.Dispersion[1], 1.20966554922812e-06)
    assert_approx_equal(si.f1[1], 14.0480530471063)
    assert_approx_equal(si.f2[1], 0.0533310749207006)
    assert_approx_equal(si.reSLD[1], 1.97775587027766e-05)
    assert_approx_equal(si.imSLD[1], 7.50821812381753e-08)
end