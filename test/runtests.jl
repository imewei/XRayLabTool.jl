using XRayLabTool
using Test

@testset "Refrac Function Tests" begin
    data = Refrac(["SiO2", "H2O"], [1, 2, 3.5, 5, 6], [2.2, 1.0])

    @test begin
        disp_value = data["SiO2"]["Dispersion"][3]
        @test abs(disp_value - 3.814915113516378e-5) < eps()
    end

    @test begin
        f1_value = data["SiO2"]["f1"][1]
        @test abs(f1_value - 29.448492570232553) < eps()
    end

    @test begin
        reSLD_value = data["SiO2"]["reSLD"][5]
        @test abs(reSLD_value - 1.8979332462750346e-5) < eps()
    end
end

@testset "SubRefrac Function Tests" begin
    si = SubRefrac("Si", [20.0], 2.33)

    @test begin
        disp_value = si["Dispersion"][1]
        @test abs(disp_value - 1.2096870811637957e-6) < eps()
    end

    @test begin
        f1_value = si["f1"][1]
        @test abs(f1_value - 14.048053002299712) < eps()
    end

    @test begin
        reSLD_value = si["reSLD"][1]
        @test abs(reSLD_value - 1.977791074150516e-5) < eps()
    end
end