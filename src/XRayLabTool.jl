module XRayLabTool

using CSV
using DataFrames
using PCHIPInterpolation
using PeriodicTable: elements
using Unitful

export Refrac, SubRefrac, XRayResult

# Define a custom struct to store X-ray results
struct XRayResult
    Formula::String
    MW::Float64
    Number_Of_Electrons::Float64
    Density::Float64
    Electron_Density::Float64
    Energy::Vector{Float64}
    Wavelength::Vector{Float64}
    Dispersion::Vector{Float64}
    Absorption::Vector{Float64}
    f1::Vector{Float64}
    f2::Vector{Float64}
    Critical_Angle::Vector{Float64}
    Attenuation_Length::Vector{Float64}
    reSLD::Vector{Float64}
    imSLD::Vector{Float64}
end

# Define constants as global variables
const thompson = 2.8179403227e-15 # m
const speedoflight = 299792458.0 # m/sec
const plank = 6.626068e-34 # m²kg/sec
const elementcharge = 1.60217646e-19 # Coulombs
const avogadro = 6.02214199e23 # mole⁻¹

# Helper function to parse a chemical formula
function parse_formula(formulaStr::String)
    elements = eachmatch(r"([A-Z][a-z]*)(\d*\.\d*|\d*)", formulaStr)
    element_symbols = [elem[1] for elem in elements]
    element_counts = parse.(Float64, [elem[2] == "" ? "1.0" : elem[2] for elem in elements])
    return element_symbols, element_counts
end

# Helper function to interpolate atomic scattering factors
function interpolate_f(E, f1, f2)
    itp1 = Interpolator(E, f1)
    itp2 = Interpolator(E, f2)
    return itp1, itp2
end

# Main function to calculate X-ray properties
function Refrac(formulaList::Vector{String}, energy::Vector{Float64}, massDensityList::Vector{Float64})
    if any(!isa(arg, Vector) for arg in [formulaList, energy, massDensityList])
        println("Invalid input: Arguments must be vectors.")
        return nothing
    end

    if any(isempty(arg) for arg in [formulaList, energy])
        println("Invalid input: Formula list and energy vector must not be empty.")
        return nothing
    end

    if any(energy .< 0.03) || any(energy .> 30)
        println("Energy is out of range 0.03KeV ~ 30KeV.")
        return nothing
    end

    if length(formulaList) != length(massDensityList)
        println("Input arguments do not match.")
        return nothing
    end

    energy = sort(energy)  # Sort energy from min to max

    results = Dict{String,XRayResult}()

    for (formula, massDensity) in zip(formulaList, massDensityList)
        results[formula] = SubRefrac(formula, energy, massDensity)
    end

    return results
end

# Subfunction to calculate X-ray properties for a single chemical formula
function SubRefrac(formulaStr::String, energy::Vector{Float64}, massDensity::Float64)
    formulaElement, element_counts = parse_formula(formulaStr)

    nElements = length(formulaElement)
    atomicNumber = []
    atomicWeight = []
    molecularWeight = 0.0
    numberOfElectrons = 0.0

    # Determine the atomic number and atomic weight
    for iElements in 1:nElements
        AN = findall(x -> x == formulaElement[iElements], [elements[iAtomicnum].symbol for iAtomicnum in eachindex(elements)],)
        push!(atomicNumber, AN[1])
        push!(atomicWeight, ustrip(elements[atomicNumber[iElements]].atomic_mass))
    end

    # Determine molecular weight and number of electrons
    for iElements in 1:nElements
        molecularWeight += element_counts[iElements] * ustrip(atomicWeight[iElements])
        numberOfElectrons += atomicNumber[iElements] * element_counts[iElements]
    end

    # Convert energy to wavelength
    wavelength = (speedoflight * plank / elementcharge) ./ (energy * 1000.0)

    # Initialize arrays
    Dispersion = zeros(length(energy))
    Absorption = zeros(length(energy))
    f1 = zeros(length(energy))
    f2 = zeros(length(energy))
    CriticalAngle = zeros(length(energy))
    ElectronDensity = 0.0
    f1f2Table = []

    # Read f1 and f2 from tables
    for iElements in 1:nElements
        fname = join([lowercase(formulaElement[iElements]), ".nff"])
        file = normpath(joinpath(@__DIR__, "AtomicScatteringFactor", fname))

        try
            table = CSV.File(file) |> DataFrame
            push!(f1f2Table, table)
        catch
            prtstr = formulaElement[iElements]
            println("Element $prtstr is NOT in the table list.")
        end
    end

    # Interpolate to get f1 and f2 for given energies
    interpf1 = []
    interpf2 = []
    for iElements in 1:nElements
        itp1, itp2 = interpolate_f(f1f2Table[iElements].E, f1f2Table[iElements].f1, f1f2Table[iElements].f2)
        push!(interpf1, [itp1(E) for E in energy * 1000])
        push!(interpf2, [itp2(E) for E in energy * 1000])
    end

    # Calculate contributions to dispersion and absorption
    for iElements in 1:nElements
        Dispersion = Dispersion .+ wavelength .^ 2 / (2 * π) * thompson * avogadro * massDensity *
                      1e6 / molecularWeight * element_counts[iElements] .* interpf1[iElements]
        Absorption = Absorption .+ wavelength .^ 2 / (2 * π) * thompson * avogadro * massDensity *
                      1e6 / molecularWeight * element_counts[iElements] .* interpf2[iElements]
        f1 = f1 .+  element_counts[iElements] .* interpf1[iElements]
        f2 = f2 .+  element_counts[iElements] .* interpf2[iElements]
    end

    ElectronDensity += 1e6 * massDensity / molecularWeight * avogadro * numberOfElectrons / 1e30

    CriticalAngle = sqrt.(2 * Dispersion) * 180 / π
    Attenuation_Length = wavelength ./ Absorption / (4 * π) * 1e2
    reSLD = Dispersion * 2 * π ./ wavelength .^ 2 / 1e20
    imSLD = Absorption * 2 * π ./ wavelength .^ 2 / 1e20

    # Create and return X-ray result struct
    result = XRayResult(formulaStr, molecularWeight, ElectronDensity, massDensity,
        ElectronDensity, energy, wavelength * 1e10, Dispersion, Absorption,
        f1, f2, CriticalAngle, Attenuation_Length, reSLD, imSLD)

    return result
end

end  # module