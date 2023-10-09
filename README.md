# XRayLabTool

## Material property when interacting with X-rays

Usage:  
    result = Refrac([Chemical Formulas], [Energy], [MassDensity])  
    result = SubRefrac(Chemical Formula, [Energy], MassDensity)  

Input:  
    1. Chemical formula. Can be either a vector (using Refrac for multiple formulas) or a single string (using SubRefrac for a single formula). The formula is case-sensitive (e.g. CO for Carbon Monoxide vs. Co for Cobalt).  
    2. Energy (0.03KeV ~ 30KeV). It has to be a vector, even one 1-element.  
    3. List of mass densities (g/cm^3). Can be either a vector of numbers (using Refrac for multiple formulas) or a single number (using SubRefrac for a single formula).  

Output:  
    Can be an array of dictionaries or a single dictionary (if the input is a single string and using SubRefrac) with the following elements:  
        1. Chemical formula;  
        2. Molecular weight;  
        3. Number of electrons per molecule;  
        4. Mass density (g/cm^3);  
        5. Electron density (1/A^3);  
        6. Real part of the atomic scattering factor (f1);  
        7. Imaginary part of the atomic scattering factor (f2);  
        8. X-ray energy (KeV);  
        9. Corresponding X-ray wavelength (A);  
        10. Dispersion;  
        11. Absorption;  
        12. Critical angle (degree);  
        13. Attenuation length (cm);  
        14. Real part of scattering length density (SLD) (A^-2);  
        15. Imaginary part of SLD (A^-2).  

Example 1: >> result=Refrac(["H2O","Si3N4"], [8.04778], [1,3.44])  
                
              Output: Dict{String, XRayResult} with 2 entries:
                          "Si3N4" => XRayResult("Si3N4", 140.283, 1.03372, 3.44, 1.03372…
                          "H2O"   => XRayResult("H2O", 18.015, 0.334285, 1.0, 0.334285, …
                
              result["Si3N4"] = XRayResult(formulaStr, molecularWeight, ElectronDensity, massDensity,
                               ElectronDensity, energy, wavelength * 1e10, Dispersion, Absorption,
                               f1, f2, CriticalAngle, Attenuation_Length, reSLD, imSLD)
                
             result["H2O"] = XRayResult(formulaStr, molecularWeight, ElectronDensity, massDensity,
                               ElectronDensity, energy, wavelength * 1e10, Dispersion, Absorption,
                               f1, f2, CriticalAngle, Attenuation_Length, reSLD, imSLD)
                
             For example,    result["Si3N4"].reSLD[1] = 2.955439362415408e-5


Example 2: >>result=SubRefrac("SiO2",Vector(8:0.5:10),2.33)  
            
               Output: result = XRayResult(formulaStr, molecularWeight, ElectronDensity, massDensity,
                                ElectronDensity, energy, wavelength * 1e10, Dispersion, Absorption,
                                f1, f2, CriticalAngle, Attenuation_Length, reSLD, imSLD)
            
               For example,    result.f1[3] = 30.366953850108544


For more information about X-ray interactions with matter, go to  
        <http://www.cxro.lbl.gov>  
        <http://www.nist.gov/>  
The atomic scattering factor table is taken from the above two websites.

**This is translated from the MATLAB script written by Zhang Jiang at the Advanced Photon Source**
