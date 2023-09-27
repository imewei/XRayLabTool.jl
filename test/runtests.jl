using XRayLabTool
using Test

@testset "XRayLabTool.jl" begin
    data = Refrac(["SiO2","H2O"],[1,2,3.5,5,6],[2.2,1.0])
    si = SubRefrac("Si", [20.0], 2.33)
    @test data["SiO2"]["Dispersion"][3] - 3.814915113516378e-5 < eps()
    @test data["SiO2"]["f1"][1] - 29.448492570232553 < eps()
    @test data["SiO2"]["reSLD"][5] - 1.8979332462750346e-5 < eps()
    @test si["Dispersion"][1] - 1.2096870811637957e-6 < eps()
    @test si["f1"][1] - 14.048053002299712 < eps()
    @test si["reSLD"][1] - 1.977791074150516e-5 < eps()
end