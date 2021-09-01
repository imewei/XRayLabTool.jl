# XRayLabTool: Material property when interacted with x-rays.

Usage: result = Refrac([Chemical Formulas],[Energy],[MassDensity])
       result = SubRefrac(Chemical Formula, [Energy], MassDensity)

Input:
    1. Chemical formula. Can be either a vector (using Refrac for multiple formulas) or a single string (using SubRefrac for a single formula). The formula is case sensitive (e.g. CO for Carbon Monoxide vs Co for Cobalt);
    2. Energy (0.03KeV~30KeV). Have to be a vector, even one 1-element.
    3. List of mass densities (g/cm^3). Can be either a vector of numbers (using Refrac for multiple formulas) or a single number (using SubRefrac for a single formula).

Output:
    Can be an array of dictionaries or a single dictionary (if input is single string and using SubRefrac) with the following elements:    
    1. Chemical formula;
    2. Molecular weight;
    3. Number of electrons per molecule;
    4. Mass density (g/cm^3);
    5. Electron density (1/A^3);
    6. Real part of the atomic scattering factor (f1)
    7. Imaginary part of the atomic scattering factor (f2);
    8. X-ray energy (KeV);
    9. Corresponding X-ray wavelength (A); 
    10. Dispersion;
    11. Absorption;
    12. Critical angle (degree);
    13. Attenuation length (cm);
    14. Real part of scattering length density (SLD) (A^-2);
    15. Imaginary part of SLD (A^-2).

Example 1: >> result=Refrac(["H2O","Si3N4"],[8.04778],[1,3.44])\n
           Output: a Dict array\n
               result["Si3N4"] =\n
               Dict{String, Any}(\n
                    "Critical Angle" => [0.2707377852086883],\n
                    "reSLD" => [2.955439361638965e-5],\n
                    "MW" => 140.283,\n
                    "Electron Density" => 1.0337188334951493,\n
                    "Wavelength" => [1.5406009069367983],\n
                    "f2" => [1.0482169679625088],\n
                    "Formula" => "Si3N4",\n
                    "Number Of Electrons" => 70.0,\n
                    "f1" => [71.0208511434114],\n
                    "Attenuation Length" => [0.0074403349864984964],\n
                    "Energy" => [8.04778],\n
                    "Dispersion" => [1.116406825816022e-5],\n
                    "imSLD" => [4.362017121420607e-7],\n
                    "Density" => 3.44,\n
                    "Absorption" => [1.6477366282283455e-7])\n
                result["H2O"] =\n
                Dict{String, Any}(\n
                    "Critical Angle" => [0.1532479196629315],\n
                    "reSLD" => [9.469204261542127e-6],\n
                    "MW" => 18.015,\n
                    "Electron Density" => 0.3342848731612545,\n
                    "Wavelength" => [1.5406009069367983],\n
                    "f2" => [0.033709251629429504],\n
                    "Formula" => "H2O",\n
                    "Number Of Electrons" => 10.0,\n
                    "f1" => [10.052289064465802],\n
                    "Attenuation Length" => [0.10220737349289194],\n
                    "Energy" => [8.04778],\n
                    "Dispersion" => [3.576958610569931e-6],\n
                    "imSLD" => [3.175394053391685e-8],\n
                    "Density" => 1.0,\n
                    "Absorption" => [1.1994939371370336e-8])\n

Example 2: >>result=SubRefrac("SiO2",Vector(8:0.5:10),2.33)\n
          Output: Dict\n
           result =\n
           Dict{String, Any} with 15 entries:\n
              "Critical Angle"      => [0.224095, 0.210836, 0.199074, 0.188608, 0.179106]\n
              "reSLD"               => [2.00086e-5, 1.9994e-5, 1.99843e-5, 1.99867e-5, 1.99707e-5]\n
              "MW"                  => 60.083\n
              "Electron Density"    => 0.70061\n
              "Wavelength"          => [1.5498, 1.45864, 1.3776, 1.3051, 1.23984]\n
              "f2"                  => [0.396933, 0.352482, 0.314986, 0.283079, 0.255704]\n
              "Formula"             => "SiO2"\n
              "Number Of Electrons" => 30.0\n
              "f1"                  => [30.4039, 30.3817, 30.367, 30.3706, 30.3463]\n
              "Attenuation Length"  => [0.0123506, 0.0147774, 0.0175093, 0.0205652, 0.0239651]\n
              "Energy"              => [8.0, 8.5, 9.0, 9.5, 10.0]\n
              ⋮                     => ⋮\n

For more information about X-ray interactions with matter, go to\n
   http://www.cxro.lbl.gov\n
   http://www.nist.gov/\n

Atomic scattering factor table is taken from the above two websites.\n

### This is translated from the MATLAB script written by Zhang Jiang at the Advanced Photon Source.\n
