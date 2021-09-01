# XRayLabTool
    ```
    ## Material property when interacted with x-rays.
    ```

Usage:<br/>
    ```
    result = Refrac([Chemical Formulas],[Energy],[MassDensity])<br/>
    result = SubRefrac(Chemical Formula, [Energy], MassDensity)<br/>
    ```

Input:<br/>
    ```
    1) Chemical formula. Can be either a vector (using Refrac for multiple formulas) or a single string (using SubRefrac for a single formula). The formula is case sensitive (e.g. CO for Carbon Monoxide vs Co for Cobalt).<br/>
    2) Energy (0.03KeV~30KeV). Have to be a vector, even one 1-element.<br/>
    3) List of mass densities (g/cm^3). Can be either a vector of numbers (using Refrac for multiple formulas) or a single number (using SubRefrac for a single formula).<br/>
    ```

Output:<br/>
    ```
    Can be an array of dictionaries or a single dictionary (if input is single string and using SubRefrac) with the following elements:<br/>   
    ```
        ```
        1) Chemical formula;<br/>
        2) Molecular weight;<br/>
        3) Number of electrons per molecule;<br/>
        4) Mass density (g/cm^3);<br/>
        5) Electron density (1/A^3);<br/>
        6) Real part of the atomic scattering factor (f1);<br/>
        7) Imaginary part of the atomic scattering factor (f2);<br/>
        8) X-ray energy (KeV);<br/>
        9) Corresponding X-ray wavelength (A); <br/>
        10) Dispersion;<br/>
        11) Absorption;<br/>
        12) Critical angle (degree);<br/>
        13) Attenuation length (cm);<br/>
        14) Real part of scattering length density (SLD) (A^-2);<br/>
        15) Imaginary part of SLD (A^-2).<br/>
        ```

Example 1: result=Refrac(["H2O","Si3N4"],[8.04778],[1,3.44])<br/>
            ```Output: a Dict array```<br/>
                ```
                result["Si3N4"] =<br/>
                Dict{String, Any}(<br/>
                ```
                    ```
                    "Critical Angle" => [0.2707377852086883],<br/>
                    reSLD" => [2.955439361638965e-5],<br/>
                    "MW" => 140.283,<br/>
                    "Electron Density" => 1.0337188334951493,<br/>
                    "Wavelength" => [1.5406009069367983],<br/>
                    "f2" => [1.0482169679625088],<br/>
                    "Formula" => "Si3N4",<br/>
                    "Number Of Electrons" => 70.0,<br/>
                    "f1" => [71.0208511434114],<br/>
                    "Attenuation Length" => [0.0074403349864984964],<br/>
                    "Energy" => [8.04778],<br/>
                    Dispersion" => [1.116406825816022e-5],<br/>
                    "imSLD" => [4.362017121420607e-7],<br/>
                    "Density" => 3.44,<br/>
                    "Absorption" => [1.6477366282283455e-7])<br/>
                    ```
                ```
                result["H2O"] =<br/>
                Dict{String, Any}(<br/>
                ```
                    ```
                    "Critical Angle" => [0.1532479196629315],<br/>
                    "reSLD" => [9.469204261542127e-6],<br/>
                    "MW" => 18.015,<br/>
                    "Electron Density" => 0.3342848731612545,<br/>
                    "Wavelength" => [1.5406009069367983],<br/>
                    "f2" => [0.033709251629429504],<br/>
                    "Formula" => "H2O",<br/>
                    "Number Of Electrons" => 10.0,<br/>
                    "f1" => [10.052289064465802],<br/>
                    "Attenuation Length" => [0.10220737349289194],<br/>
                    "Energy" => [8.04778],<br/>
                    "Dispersion" => [3.576958610569931e-6],<br/>
                    "imSLD" => [3.175394053391685e-8],<br/>
                    "Density" => 1.0,<br/>
                    "Absorption" => [1.1994939371370336e-8])<br/>
                    ```

Example 2: result=SubRefrac("SiO2",Vector(8:0.5:10),2.33)<br/>
            ```Output: Dict<br/>```
                ```
                result =<br/>
                Dict{String, Any} with 15 entries:<br/>
                ```
                    ```
                    "Critical Angle"      => [0.224095, 0.210836, 0.199074, 0.188608, 0.179106]<br/>
                    "reSLD"               => [2.00086e-5, 1.9994e-5, 1.99843e-5, 1.99867e-5, 1.99707e-5]<br/>
                    "MW"                  => 60.083<br/>
                    "Electron Density"    => 0.70061<br/>
                    "Wavelength"          => [1.5498, 1.45864, 1.3776, 1.3051, 1.23984]<br/>
                    "f2"                  => [0.396933, 0.352482, 0.314986, 0.283079, 0.255704]<br/>
                    "Formula"             => "SiO2"<br/>
                    "Number Of Electrons" => 30.0<br/>
                    "f1"                  => [30.4039, 30.3817, 30.367, 30.3706, 30.3463]<br/>
                    "Attenuation Length"  => [0.0123506, 0.0147774, 0.0175093, 0.0205652, 0.0239651]<br/>
                    "Energy"              => [8.0, 8.5, 9.0, 9.5, 10.0]<br/>
                    ⋮                     => ⋮<br/>
                    ```

For more information about X-ray interactions with matter, go to<br/>
    ```
    http://www.cxro.lbl.gov<br/>
    http://www.nist.gov/<br/>
    
    Atomic scattering factor table is taken from the above two websites.
    ```

        ```
        #### This is translated from the MATLAB script written by Zhang Jiang at the Advanced Photon Source.
        ```
