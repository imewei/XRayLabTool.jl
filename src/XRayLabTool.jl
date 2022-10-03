module XRayLabTool
#=
XRayLabTool  Material property when interacted with x-rays.

Usage: result = Refrac([Chemical Formulas],[Energy],[MassDensity])
       result = SubRefrac(Chemical Formula, [Energy], MassDensity)

Input:
    1) Chemical formula. Can be either a vector (using Refrac for multiple formulas) or a single string (using SubRefrac for a single formula). The formula is case sensitive (e.g. CO for Carbon Monoxide vs Co for Cobalt);
    2) Energy (0.03KeV~30KeV). Have to be a vector, even one 1-element.
    3) List of mass densities (g/cm^3). Can be either a vector of numbers (using Refrac for multiple formulas) or a single number (using SubRefrac for a single formula).

Output:
    Can be an array of dictionaries or a single dictionary (if input is single string and using SubRefrac) with the following elements:
    1) Chemical formula;
    2) Molecular weight;
    3) Number of electrons per molecule;
    4) Mass density (g/cm^3);
    5) Electron density (1/A^3);
    6) Real part of the atomic scattering factor (f1)
    7) Imaginary part of the atomic scattering factor (f2);
    8) X-ray energy (KeV);
    9) Corresponding X-ray wavelength (A);
    10) Dispersion;
    11) Absorption;
    12) Critical angle (degree);
    13) Attenuation length (cm);
    14) Real part of scattering length density (SLD) (A^-2);
    15) Imaginary part of SLD (A^-2).

Example 1: >> result=Refrac(["H2O","Si3N4"],[8.04778],[1,3.44])
           Output: a Dict array
               result["Si3N4"] =
               Dict{String, Any}(
                    "Critical Angle" => [0.2707377852086883],
                    "reSLD" => [2.955439361638965e-5],
                    "MW" => 140.283,
                    "Electron Density" => 1.0337188334951493,
                    "Wavelength" => [1.5406009069367983],
                    "f2" => [1.0482169679625088],
                    "Formula" => "Si3N4",
                    "Number Of Electrons" => 70.0,
                    "f1" => [71.0208511434114],
                    "Attenuation Length" => [0.0074403349864984964],
                    "Energy" => [8.04778],
                    "Dispersion" => [1.116406825816022e-5],
                    "imSLD" => [4.362017121420607e-7],
                    "Density" => 3.44,
                    "Absorption" => [1.6477366282283455e-7])
                result["H2O"] =
                Dict{String, Any}(
                    "Critical Angle" => [0.1532479196629315],
                    "reSLD" => [9.469204261542127e-6],
                    "MW" => 18.015,
                    "Electron Density" => 0.3342848731612545,
                    "Wavelength" => [1.5406009069367983],
                    "f2" => [0.033709251629429504],
                    "Formula" => "H2O",
                    "Number Of Electrons" => 10.0,
                    "f1" => [10.052289064465802],
                    "Attenuation Length" => [0.10220737349289194],
                    "Energy" => [8.04778],
                    "Dispersion" => [3.576958610569931e-6],
                    "imSLD" => [3.175394053391685e-8],
                    "Density" => 1.0,
                    "Absorption" => [1.1994939371370336e-8])

Example 2: >>result=SubRefrac("SiO2",Vector(8:0.5:10),2.33)
          Output: Dict
           result =
           Dict{String, Any} with 15 entries:
              "Critical Angle"      => [0.224095, 0.210836, 0.199074, 0.188608, 0.179106]
              "reSLD"               => [2.00086e-5, 1.9994e-5, 1.99843e-5, 1.99867e-5, 1.99707e-5]
              "MW"                  => 60.083
              "Electron Density"    => 0.70061
              "Wavelength"          => [1.5498, 1.45864, 1.3776, 1.3051, 1.23984]
              "f2"                  => [0.396933, 0.352482, 0.314986, 0.283079, 0.255704]
              "Formula"             => "SiO2"
              "Number Of Electrons" => 30.0
              "f1"                  => [30.4039, 30.3817, 30.367, 30.3706, 30.3463]
              "Attenuation Length"  => [0.0123506, 0.0147774, 0.0175093, 0.0205652, 0.0239651]
              "Energy"              => [8.0, 8.5, 9.0, 9.5, 10.0]
              ⋮                     => ⋮

For more information about X-ray interactions with matter, go to
   http://www.cxro.lbl.gov
   http://www.nist.gov/

Atomic scattering factor table is taken from the above two websites.

This is translated from Zhang Jiang's MATLAB script
=#

using BenchmarkTools
using CSV
using DataFrames
using Interpolations
using LinearAlgebra
using PCHIPInterpolation
using PeriodicTable
using SpecialFunctions
using Unitful

export Refrac, SubRefrac

function Refrac(formulaList, energy, massDensityList)

    if !isa(formulaList, Vector{String})
        println("Need a styring for chemical formula input argument.")
        return
    end
    if !isa(energy, Vector) | isempty(energy)
        println("Invalid x-ray energy.")
        return
    end
    if count(!iszero, energy .< 0.03) | count(!iszero, energy .> 30) > 0
        println("Energy is out of range 0.03KeV~30KeV.")
        return
    end
    if !isa(massDensityList, Vector)
        println("Invalid mass density.")
        return
    end
    if length(formulaList) != length(massDensityList)
        println("Input arguements do not match.")
        return
    end

    # --- sort energy list and
    energy = sort!(energy)  # sort energy from min to max

    result = Dict()

    for (ii, formulaListItem) in enumerate(formulaList)
        result[formulaList[ii]] =
            SubRefrac(formulaList[ii], energy, massDensityList[ii])
    end
    return result
end


function SubRefrac(formulaStr, energy, massDensity)
    # Some physical constants
    thompson = 2.8179403227e-15 # m
    speedoflight = 299792458 # m/sec
    plank = 6.626068e-34 # m²kg/sec
    elementcharge = 1.60217646e-19 # Coulombs
    avogadro = 6.02214199e23 # mole⁻¹

    # convert energy to wavelength
    if !isa(energy, Vector)
        println("Invalid x-ray energy.")
    end
    wavelength = (speedoflight * plank / elementcharge) ./ (energy * 1000.0)

    # determine elements and number of atoms in the formula
    nElements = 0
    formulaElement = []
    nAtoms = []

    formulaChar = split(formulaStr, "")

    for (iFormula, formulaStrItem) in enumerate(formulaStr)
        if count(!iszero, formulaChar[iFormula] <= "Z") &
           count(!iszero, formulaChar[iFormula] >= "A") > 0
            nElements = nElements + 1
            push!(formulaElement, formulaChar[iFormula])
            push!(nAtoms, "0")
        elseif count(!iszero, formulaChar[iFormula] <= "z") &
               count(!iszero, formulaChar[iFormula] >= "a") &
               (
                   (
                       count(!iszero, string(formulaChar[iFormula-1]) <= "Z") &
                       count(!iszero, string(formulaChar[iFormula-1]) >= "A")
                   ) | (
                       count(!iszero, string(formulaChar[iFormula-1]) <= "z") &
                       count(!iszero, string(formulaChar[iFormula-1]) >= "a")
                   )
               ) > 0
            formulaElement[nElements] =
                join([formulaElement[nElements], formulaChar[iFormula]])
        elseif (
            count(!iszero, formulaChar[iFormula] <= "9") &
            count(!iszero, formulaChar[iFormula] >= "0")
        ) | count(!iszero, formulaChar[iFormula] == ".") > 0
            nAtoms[nElements] = join([nAtoms[nElements], formulaChar[iFormula]])
        else
            println("Invalid chemical formula.")
            return
        end
    end

    for iElements in 1:nElements
        nAtoms[iElements] = parse(Float64, nAtoms[iElements])
        if abs.(nAtoms[iElements]) .- 0.0 < eps()
            nAtoms[iElements] = 1.0
        end
    end

    # read f1 and f2 from tables
    f1f2Table = []
    for iElements in 1:nElements
        fname = join([lowercase(formulaElement[iElements]), ".nff"])
        file = normpath(joinpath(@__DIR__, "AtomicScatteringFactor", fname))
        # println(file)
        try
            table = CSV.File(file) |> DataFrame
            push!(f1f2Table, table)
        catch
            prtstr = formulaElement[iElements]
            println("Element $prtstr is NOT in the table list.")
        end
    end

    # determine the atomic number and atomic weight
    atomicNumber = []
    atomicWeight = []
    for iElements in 1:nElements
        AN = findall(
            x -> x == formulaElement[iElements],
            [elements[iAtomicnum].symbol for (iAtomicnum, elementsItem) in enumerate(elements)],
        )
        push!(atomicNumber, AN[1])
        push!(atomicWeight, elements[atomicNumber[iElements]].atomic_mass)
    end

    # determine molecular weight and number of electrons
    molecularWeight = 0
    numberOfElectrons = 0
    for iElements in 1:nElements
        molecularWeight =
            molecularWeight +
            nAtoms[iElements] * ustrip(atomicWeight[iElements])
        numberOfElectrons =
            numberOfElectrons + atomicNumber[iElements] * nAtoms[iElements]
    end

    # interpolate to get f1 and f2 for given energies
    interpf1 = []
    interpf2 = []
    for iElements in 1:nElements
        # use PCHIPInterpolation
        itp1 = Interpolator(f1f2Table[iElements].E, f1f2Table[iElements].f1)
        itp2 = Interpolator(f1f2Table[iElements].E, f1f2Table[iElements].f2)
        # use Interpolations
        # itp1 = interpolate((f1f2Table[iElements].E,), f1f2Table[iElements].f1, Gridded(Linear()))
        # itp2 = interpolate((f1f2Table[iElements].E,), f1f2Table[iElements].f2, Gridded(Linear()))
        push!(interpf1, [itp1(E) for E in energy * 1000])
        push!(interpf2, [itp2(E) for E in energy * 1000])
    end

    # calculate dispersion, absorption, critical angle and attenuation length
    Dispersion = zeros(length(energy))
    Absorption = zeros(length(energy))
    f1 = zeros(length(energy))
    f2 = zeros(length(energy))
    CriticalAngle = zeros(length(energy))
    ElectronDensity =
        1e6 * massDensity / molecularWeight * avogadro * numberOfElectrons /
        1e30
    for iElements in 1:nElements
        Dispersion =
            Dispersion +
            wavelength .^ 2 / (2 * pi) *
            thompson *
            avogadro *
            massDensity *
            1e6 / molecularWeight * nAtoms[iElements] .* interpf1[iElements]
        Absorption =
            Absorption +
            wavelength .^ 2 / (2 * pi) *
            thompson *
            avogadro *
            massDensity *
            1e6 / molecularWeight * nAtoms[iElements] .* interpf2[iElements]
        f1 = f1 + nAtoms[iElements] .* interpf1[iElements]
        f2 = f2 + nAtoms[iElements] .* interpf2[iElements]
    end
    CriticalAngle = sqrt.((2 * Dispersion)) * 180 / pi
    AttLength = wavelength ./ Absorption / (4 * pi) * 1e2
    reSLD = Dispersion * 2 * pi ./ wavelength .^ 2 / 1e20
    imSLD = Absorption * 2 * pi ./ wavelength .^ 2 / 1e20
    Xresult = Dict(
        "Formula" => formulaStr,
        "MW" => molecularWeight,
        "Number Of Electrons" => numberOfElectrons,
        "Density" => massDensity,
        "Electron Density" => ElectronDensity,
        "Energy" => energy,
        "Wavelength" => wavelength * 1e10,
        "Dispersion" => Dispersion,
        "Absorption" => Absorption,
        "f1" => f1,
        "f2" => f2,
        "Critical Angle" => CriticalAngle,
        "Attenuation Length" => AttLength,
        "reSLD" => reSLD,
        "imSLD" => imSLD,
    )
    return Xresult
end
end